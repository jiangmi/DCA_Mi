// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Anderson/Pulay (DIIS) mixer for the cluster self-energy self-consistency loop.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SELF_ENERGY_MIXER_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SELF_ENERGY_MIXER_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <deque>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dca {
namespace phys {
namespace solver {
// dca::phys::solver::

// Anderson/Pulay (DIIS) mixer for the self-energy fixed-point iteration x -> F(x), with residual
// r = F(x) - x:
//   x_new = sum_i c_i (x_i + alpha * r_i),
// where the coefficients c_i minimize |sum_i c_i r_i|^2 subject to sum_i c_i = 1.
// For a single history entry this reduces to plain linear mixing, x + alpha * r.
//
// Parallel contract: this is a purely local, stateful object and performs no communication. It must
// be called exactly once per DCA iteration on every MPI rank and with identical inputs, so that the
// replicated histories stay in sync. Putting a call behind a rank guard, or calling it a
// rank-dependent number of times, silently makes the ranks diverge.
class SelfEnergyMixer {
public:
  SelfEnergyMixer() = default;

  void configure(const double alpha, const int max_history) {
    if (max_history < 1)
      throw std::logic_error("SelfEnergyMixer: max_history must be >= 1.");
    alpha_ = alpha;
    max_history_ = max_history;
  }

  void reset() {
    history_.clear();
    size_ = 0;
    last_coefficients_.clear();
    last_fallback_ = false;
  }

  std::size_t historySize() const {
    return history_.size();
  }

  // Coefficients used by the most recent mix() call: the DIIS weights c_j for a regular Anderson
  // step, or {1} for the single-entry / ill-conditioned linear-mixing fallback.
  const std::vector<double>& lastCoefficients() const {
    return last_coefficients_;
  }

  // True if the most recent mix() call fell back to linear mixing (history too short or the
  // DIIS system was ill-conditioned).
  bool lastFellBack() const {
    return last_fallback_;
  }

  // In-place safe: sigma_out may alias sigma_qmc and/or sigma_in.
  void mix(const std::complex<double>* const sigma_qmc,
           const std::complex<double>* const sigma_in, std::complex<double>* const sigma_out,
           const int size) {
    if (max_history_ < 1)
      throw std::logic_error("SelfEnergyMixer: configure() must be called before mix().");
    if (size <= 0)
      throw std::logic_error("SelfEnergyMixer: the size of the self-energy must be positive.");
    // Guard against a rank-dependent or iteration-dependent size, which would make the stored
    // history inconsistent with the vector that is being mixed.
    if (size_ == 0)
      size_ = size;
    else if (size_ != size)
      throw std::logic_error("SelfEnergyMixer: the size of the self-energy must not change.");

    const std::vector<std::complex<double>> y(sigma_qmc, sigma_qmc + size);
    const std::vector<std::complex<double>> x(sigma_in, sigma_in + size);

    std::vector<std::complex<double>> r(size);
    for (int i = 0; i < size; ++i)
      r[i] = y[i] - x[i];

    history_.push_back(Entry{x, r});
    while (static_cast<int>(history_.size()) > max_history_)
      history_.pop_front();

    std::vector<std::complex<double>> x_new(size);

    std::vector<double> c;
    if (history_.size() > 1)
      c = solveCoefficients(size);

    if (history_.size() == 1 || c.empty()) {
      // Single history entry (equivalent to linear mixing) or ill-conditioned fallback.
      for (int i = 0; i < size; ++i)
        x_new[i] = x[i] + alpha_ * r[i];
      last_coefficients_ = std::vector<double>{1.};
      last_fallback_ = true;
    }
    else {
      std::fill(x_new.begin(), x_new.end(), std::complex<double>(0., 0.));
      const int m = static_cast<int>(history_.size());
      for (int j = 0; j < m; ++j)
        for (int i = 0; i < size; ++i)
          x_new[i] += c[j] * (history_[j].x[i] + alpha_ * history_[j].r[i]);
      last_coefficients_ = c;
      last_fallback_ = false;
    }

    std::copy(x_new.begin(), x_new.end(), sigma_out);
  }

private:
  struct Entry {
    std::vector<std::complex<double>> x;
    std::vector<std::complex<double>> r;
  };

  // Solves the DIIS linear system for the coefficients c.
  // Returns an empty vector if the system is ill-conditioned.
  std::vector<double> solveCoefficients(const int size) const {
    const int m = static_cast<int>(history_.size());
    std::vector<double> a((m + 1) * (m + 1), 0.);
    double max_diag = 0.;

    // Gram matrix of the residual vectors, treated as real vectors of length 2 * size.
    for (int i = 0; i < m; ++i) {
      for (int j = i; j < m; ++j) {
        double dot = 0.;
        for (int k = 0; k < size; ++k) {
          const auto& ri = history_[i].r[k];
          const auto& rj = history_[j].r[k];
          dot += ri.real() * rj.real() + ri.imag() * rj.imag();
        }
        a[i * (m + 1) + j] = dot;
        a[j * (m + 1) + i] = dot;
      }
      max_diag = std::max(max_diag, std::abs(a[i * (m + 1) + i]));
    }

    // Tikhonov regularization.
    const double lambda = regularization_ * max_diag + 1.e-300;
    for (int i = 0; i < m; ++i)
      a[i * (m + 1) + i] += lambda;

    // Constraint: sum_i c_i = 1.
    for (int i = 0; i < m; ++i) {
      a[i * (m + 1) + m] = 1.;
      a[m * (m + 1) + i] = 1.;
    }
    a[m * (m + 1) + m] = 0.;

    std::vector<double> b(m + 1, 0.);
    b[m] = 1.;

    std::vector<double> solution(m + 1);
    if (!solveLinearSystem(a, b, solution, m + 1))
      return {};

    solution.resize(m);
    return solution;
  }

  // Dense Gaussian elimination with partial pivoting. Returns false if the system is singular.
  static bool solveLinearSystem(std::vector<double> a, std::vector<double> b,
                                std::vector<double>& x, const int n) {
    for (int col = 0; col < n; ++col) {
      int pivot = col;
      for (int row = col + 1; row < n; ++row)
        if (std::abs(a[row * n + col]) > std::abs(a[pivot * n + col]))
          pivot = row;
      if (std::abs(a[pivot * n + col]) < 1.e-14)
        return false;
      if (pivot != col) {
        for (int k = 0; k < n; ++k)
          std::swap(a[col * n + k], a[pivot * n + k]);
        std::swap(b[col], b[pivot]);
      }
      for (int row = col + 1; row < n; ++row) {
        const double factor = a[row * n + col] / a[col * n + col];
        for (int k = col; k < n; ++k)
          a[row * n + k] -= factor * a[col * n + k];
        b[row] -= factor * b[col];
      }
    }

    x.assign(n, 0.);
    for (int row = n - 1; row >= 0; --row) {
      double sum = b[row];
      for (int k = row + 1; k < n; ++k)
        sum -= a[row * n + k] * x[k];
      x[row] = sum / a[row * n + row];
      if (!std::isfinite(x[row]))
        return false;
    }
    return true;
  }

  double alpha_ = 1.;
  int max_history_ = 0;
  int size_ = 0;
  double regularization_ = 1.e-8;
  std::deque<Entry> history_;
  std::vector<double> last_coefficients_;
  bool last_fallback_ = false;
};

}  // solver
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_SELF_ENERGY_MIXER_HPP

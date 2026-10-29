// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// This file tests self_energy_mixer.hpp.

#include "dca/phys/dca_step/cluster_solver/self_energy_mixer.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "gtest/gtest.h"

namespace {

using dca::phys::solver::SelfEnergyMixer;
using Complex = std::complex<double>;

}  // namespace

TEST(SelfEnergyMixerTest, ThrowsOnInvalidConfiguration) {
  SelfEnergyMixer mixer;
  EXPECT_THROW(mixer.configure(1., 0), std::logic_error);

  std::vector<Complex> v{Complex(0., 0.)};
  EXPECT_THROW(mixer.mix(v.data(), v.data(), v.data(), 1), std::logic_error);
}

TEST(SelfEnergyMixerTest, SingleHistoryReducesToLinearMixing) {
  SelfEnergyMixer mixer;
  mixer.configure(0.5, 1);

  const std::vector<Complex> qmc{Complex(2., 1.), Complex(-1., 0.5), Complex(0.25, -3.)};
  const std::vector<Complex> in{Complex(0.5, 0.), Complex(0.25, 0.1), Complex(-1., 2.)};
  std::vector<Complex> out(in.size());

  mixer.mix(qmc.data(), in.data(), out.data(), static_cast<int>(out.size()));

  for (std::size_t i = 0; i < out.size(); ++i) {
    const Complex expected = in[i] + 0.5 * (qmc[i] - in[i]);
    EXPECT_NEAR(0., std::abs(out[i] - expected), 1.e-14);
  }
}

TEST(SelfEnergyMixerTest, IsInPlaceSafe) {
  SelfEnergyMixer mixer;
  mixer.configure(1., 1);

  std::vector<Complex> qmc{Complex(2., 1.), Complex(-1., 0.5)};
  const std::vector<Complex> in{Complex(0.5, 0.), Complex(0.25, 0.)};
  const std::vector<Complex> qmc_ref = qmc;

  // With alpha == 1 and a single history entry the output equals the QMC self-energy.
  mixer.mix(qmc.data(), in.data(), qmc.data(), static_cast<int>(qmc.size()));

  for (std::size_t i = 0; i < qmc.size(); ++i)
    EXPECT_NEAR(0., std::abs(qmc[i] - qmc_ref[i]), 1.e-14);
}

// On a linear contraction x -> x* + 0.5 (x - x*) the Anderson/Pulay mixer must converge to the
// fixed point in a couple of iterations (Aitken-like extrapolation), much faster than the plain
// linear mixing (which only halves the error each step).
TEST(SelfEnergyMixerTest, AndersonAcceleratesLinearContraction) {
  const int n = 4;
  const std::vector<Complex> fixed(n, Complex(1., -0.5));

  const auto compute_new = [&](const std::vector<Complex>& x) {
    std::vector<Complex> y(n);
    for (int i = 0; i < n; ++i)
      y[i] = fixed[i] + 0.5 * (x[i] - fixed[i]);
    return y;
  };

  SelfEnergyMixer mixer;
  mixer.configure(1., 4);

  std::vector<Complex> x(n, Complex(0., 0.));
  double residual = 0.;
  for (int it = 0; it < 4; ++it) {
    const std::vector<Complex> y = compute_new(x);
    mixer.mix(y.data(), x.data(), x.data(), n);

    residual = 0.;
    for (int i = 0; i < n; ++i)
      residual += std::norm(x[i] - fixed[i]);
    residual = std::sqrt(residual);
  }

  EXPECT_LT(residual, 1.e-6);
}

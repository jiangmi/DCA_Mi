// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class measures the equal time operator functions.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_TP_EQUAL_TIME_ACCUMULATOR_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_TP_EQUAL_TIME_ACCUMULATOR_HPP

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/math/interpolation/akima_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain_left_oriented.hpp"
#include "dca/phys/domains/time_and_frequency/vertex_time_domain.hpp"
#include "dca/util/plot.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <class parameters_type, class MOMS_type>
class TpEqualTimeAccumulator {
public:
  typedef double scalar_type;
  typedef vertex_singleton vertex_singleton_type;

  typedef typename parameters_type::profiler_type profiler_type;
  typedef typename parameters_type::concurrency_type concurrency_type;

  using t = func::dmn_0<domains::time_domain>;
  using t_VERTEX = func::dmn_0<domains::vertex_time_domain<domains::TP_TIME_DOMAIN_POSITIVE>>;
  using t_sp     = func::dmn_0<domains::vertex_time_domain<domains::SP_TIME_DOMAIN_POSITIVE>>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;
  using r_dmn_t = RClusterDmn;
  using k_dmn_t = KClusterDmn;

  typedef func::dmn_variadic<b, r_dmn_t, t_VERTEX> b_r_t_VERTEX_dmn_t;

  typedef func::dmn_0<domains::time_domain_left_oriented> shifted_t;
  typedef func::dmn_variadic<nu, nu, r_dmn_t, shifted_t> nu_nu_r_dmn_t_shifted_t;

  typedef func::dmn_0<func::dmn<4, int>> akima_dmn_t;
  typedef func::dmn_variadic<akima_dmn_t, nu, nu, r_dmn_t, shifted_t> akima_nu_nu_r_dmn_t_shifted_t;

public:
  TpEqualTimeAccumulator(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id);

  void resetAccumulation();

  void finalize();

  void sumTo(TpEqualTimeAccumulator<parameters_type, MOMS_type>& other) const;

  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& get_G_r_t() {
    return G_r_t;
  }
  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& get_G_r_t_stddev() {
    return G_r_t_stddev;
  }

  func::function<double, func::dmn_variadic<b, r_dmn_t>>& get_charge_cluster_moment() {
    return charge_cluster_moment;
  }
  func::function<double, func::dmn_variadic<b, r_dmn_t>>& get_magnetic_cluster_moment() {
    return magnetic_cluster_moment;
  }
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& get_spin_XX_chi() {
    return spin_XX_chi;
  }
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& get_spin_ZZ_chi() {
    return spin_ZZ_chi;
  }
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& get_charge_chi() {
    return charge_chi;
  }
  func::function<double, func::dmn_variadic<b, r_dmn_t>>& get_dwave_pp_correlator() {
    return dwave_pp_correlator;
  }

  template <class configuration_type, typename RealInp>
  void compute_G_r_t(const configuration_type& configuration_e_up,
                     const dca::linalg::Matrix<RealInp, dca::linalg::CPU>& M_up,
                     const configuration_type& configuration_e_dn,
                     const dca::linalg::Matrix<RealInp, dca::linalg::CPU>& M_dn);

  void accumulate_G_r_t(double sign);

  void accumulate_moments(double sign);

  void accumulate_chi(double sign);

  void accumulate_dwave_pp_correlator(double sign);

  // Accumulate all relevant quantities. This is equivalent to calling compute_G_r_t followed by all
  // the accumulation methods.
  template <class configuration_type, typename RealInp>
  void accumulateAll(const configuration_type& configuration_e_up,
                     const dca::linalg::Matrix<RealInp, dca::linalg::CPU>& M_up,
                     const configuration_type& configuration_e_dn,
                     const dca::linalg::Matrix<RealInp, dca::linalg::CPU>& M_dn, int sign);

  double get_GFLOP();

private:
  void initialize_my_configuration();
  void initialize_akima_coefficients();

  void initialize_G0_indices();
  void initialize_G0_original();
  void test_G0_original();

  void interpolate(func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& G_r_t,
                   func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& G_r_t_stddev,
                   func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& spin_XX_chi,
                   func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& spin_ZZ_chi,
                   func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& charge_chi);

  int find_first_non_interacting_spin(const std::vector<vertex_singleton_type>& configuration_e_spin);

  template <class configuration_type>
  void compute_G0_matrix(e_spin_states e_spin, const configuration_type& configuration,
                         dca::linalg::Matrix<float, dca::linalg::CPU>& G0_matrix);

  template <class configuration_type>
  void compute_G0_matrix_left(e_spin_states e_spin, const configuration_type& configuration,
                              dca::linalg::Matrix<float, dca::linalg::CPU>& G0_matrix);

  template <class configuration_type>
  void compute_G0_matrix_right(e_spin_states e_spin, const configuration_type& configuration,
                               dca::linalg::Matrix<float, dca::linalg::CPU>& G0_matrix);

  double interpolate_akima(int b_i, int s_i, int b_j, int s_j, int delta_r, double tau);

private:
  struct singleton_operator {
    int b_ind;
    int r_ind;
    int t_ind;

    double t_val;
  };

private:
  parameters_type& parameters;
  concurrency_type& concurrency;
  MOMS_type& MOMS;

  int thread_id;
  double GFLOP;

  b_r_t_VERTEX_dmn_t b_r_t_dmn;
  nu_nu_r_dmn_t_shifted_t nu_nu_r_dmn_t_t_shifted_dmn;

  func::function<double, akima_nu_nu_r_dmn_t_shifted_t> akima_coefficients;

  std::vector<singleton_operator> fixed_configuration;
  std::vector<singleton_operator> ctaux_configuration;

  func::function<float, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                           func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G0_sign_up;
  func::function<float, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                           func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G0_sign_dn;

  func::function<int, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                         func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G0_indices_up;
  func::function<int, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                         func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G0_indices_dn;

  func::function<float, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                           func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G0_integration_factor_up;
  func::function<float, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                           func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G0_integration_factor_dn;

  dca::linalg::Matrix<float, dca::linalg::CPU> G0_original_up;
  dca::linalg::Matrix<float, dca::linalg::CPU> G0_original_dn;

  dca::linalg::Matrix<float, dca::linalg::CPU> M_matrix_up;
  dca::linalg::Matrix<float, dca::linalg::CPU> M_matrix_dn;

  dca::linalg::Matrix<float, dca::linalg::CPU> G0_matrix_up;
  dca::linalg::Matrix<float, dca::linalg::CPU> G0_matrix_dn;

  dca::linalg::Matrix<float, dca::linalg::CPU> G0_matrix_up_left;
  dca::linalg::Matrix<float, dca::linalg::CPU> G0_matrix_dn_left;

  dca::linalg::Matrix<float, dca::linalg::CPU> G0_matrix_up_right;
  dca::linalg::Matrix<float, dca::linalg::CPU> G0_matrix_dn_right;

  dca::linalg::Matrix<float, dca::linalg::CPU> M_G0_matrix_up;
  dca::linalg::Matrix<float, dca::linalg::CPU> M_G0_matrix_dn;

  dca::linalg::Matrix<float, dca::linalg::CPU> G0_M_G0_matrix_up;
  dca::linalg::Matrix<float, dca::linalg::CPU> G0_M_G0_matrix_dn;

  func::function<float, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                           func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G_r_t_dn;
  func::function<float, func::dmn_variadic<func::dmn_variadic<b, r_dmn_t, t_VERTEX>,
                                           func::dmn_variadic<b, r_dmn_t, t_VERTEX>>>
      G_r_t_up;

  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>> G_r_t;
  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>> G_r_t_stddev;

  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t_VERTEX>> G_r_t_accumulated;
  func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t_VERTEX>> G_r_t_accumulated_squared;

  func::function<double, func::dmn_variadic<b, r_dmn_t>> charge_cluster_moment;
  func::function<double, func::dmn_variadic<b, r_dmn_t>> magnetic_cluster_moment;

  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>> spin_XX_chi;
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>> spin_ZZ_chi;
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>> charge_chi;
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_VERTEX>> spin_XX_chi_accumulated;
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_VERTEX>> spin_ZZ_chi_accumulated;
  func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_VERTEX>> charge_chi_accumulated;

  func::function<double, k_dmn_t> dwave_k_factor;
  func::function<double, r_dmn_t> dwave_r_factor;

  func::function<double, func::dmn_variadic<b, r_dmn_t>> dwave_pp_correlator;
};

template <class parameters_type, class MOMS_type>
TpEqualTimeAccumulator<parameters_type, MOMS_type>::TpEqualTimeAccumulator(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id)
    : parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      MOMS(MOMS_ref),

      thread_id(id),

      GFLOP(0),

      G_r_t_accumulated("G_r_t_accumulated"),
      G_r_t_accumulated_squared("G_r_t_accumulated_squared"),

      charge_cluster_moment("charge-cluster-moment"),
      magnetic_cluster_moment("magnetic-cluster-moment"),
      spin_XX_chi("spin-XX-susceptibility"),
      spin_ZZ_chi("spin-ZZ-susceptibility"),
      charge_chi("charge-susceptibility"),
      spin_XX_chi_accumulated("spin-XX-susceptibility_accumulated"),
      spin_ZZ_chi_accumulated("spin-ZZ-susceptibility_accumulated"),
      charge_chi_accumulated("charge-susceptibility_accumulated"),

      dwave_pp_correlator("dwave-pp-correlator") {
  for (int k_ind = 0; k_ind < k_dmn_t::dmn_size(); k_ind++)
    dwave_k_factor(k_ind) =
        cos(k_dmn_t::get_elements()[k_ind][0]) - cos(k_dmn_t::get_elements()[k_ind][1]);

  math::transform::FunctionTransform<k_dmn_t, r_dmn_t>::execute(dwave_k_factor, dwave_r_factor);

  initialize_my_configuration();

  initialize_akima_coefficients();

  initialize_G0_indices();

  initialize_G0_original();
}

template <class parameters_type, class MOMS_type>
double TpEqualTimeAccumulator<parameters_type, MOMS_type>::get_GFLOP() {
  double tmp = GFLOP;
  GFLOP = 0;
  return tmp;
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::resetAccumulation() {
  GFLOP = 0;

  G_r_t = 0;
  G_r_t_stddev = 0;

  G_r_t_accumulated = 0;
  G_r_t_accumulated_squared = 0;

  charge_cluster_moment = 0;
  magnetic_cluster_moment = 0;

  spin_XX_chi = 0;
  spin_ZZ_chi = 0;
  charge_chi = 0;
  spin_XX_chi_accumulated = 0;
  spin_ZZ_chi_accumulated = 0;
  charge_chi_accumulated = 0;

  dwave_pp_correlator = 0;
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::initialize_my_configuration() {
  fixed_configuration.resize(b::dmn_size() * r_dmn_t::dmn_size() * t_VERTEX::dmn_size());

  int index = 0;
  for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int t_ind = 0; t_ind < t_VERTEX::dmn_size(); t_ind++) {
        singleton_operator tmp;

        tmp.b_ind = b_ind;
        tmp.r_ind = r_ind;
        tmp.t_ind = t_ind;
        tmp.t_val = t_VERTEX::get_elements()[t_ind];

        fixed_configuration[index] = tmp;

        index += 1;
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::initialize_akima_coefficients() {
  int size = t::dmn_size() / 2;

  math::interpolation::akima_interpolation<double> ai_obj(size);

  double* x = new double[size];
  double* y = new double[size];

  for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
    x[t_ind] = t_ind;

  {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < nu::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < nu::dmn_size(); nu0_ind++) {
          for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
            y[t_ind] = MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);

          ai_obj.initialize(x, y);

          for (int t_ind = 0; t_ind < t::dmn_size() / 2 - 1; t_ind++)
            for (int l = 0; l < 4; l++)
              akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind) = ai_obj.get_alpha(l, t_ind);
        }
      }
    }
  }

  {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < nu::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < nu::dmn_size(); nu0_ind++) {
          for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size(); t_ind++)
            y[t_ind - t::dmn_size() / 2] =
                MOMS.G0_r_t_cluster_excluded(nu0_ind, nu1_ind, r_ind, t_ind);

          ai_obj.initialize(x, y);

          for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size() - 1; t_ind++)
            for (int l = 0; l < 4; l++)
              akima_coefficients(l, nu0_ind, nu1_ind, r_ind, t_ind - 1) =
                  ai_obj.get_alpha(l, t_ind - t::dmn_size() / 2);
        }
      }
    }
  }

  delete[] x;
  delete[] y;
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::initialize_G0_indices() {
  std::vector<double> multiplicities(t_VERTEX::dmn_size(), 0);
  {
    for (int i = 0; i < t_VERTEX::dmn_size(); i++) {
      for (int j = 0; j < t_VERTEX::dmn_size(); j++) {
        int t_ind = i - j;

        if (std::abs(t_VERTEX::get_elements().back() - parameters.get_beta()) < 1.e-6)
          t_ind = t_ind < 0 ? t_ind + t_VERTEX::dmn_size() - 1 : t_ind;
        else
          t_ind = t_ind < 0 ? t_ind + t_VERTEX::dmn_size() - 0 : t_ind;

        for (int l = 0; l < t_VERTEX::dmn_size(); l++)
          if (std::abs(t_VERTEX::get_elements()[l] - t_VERTEX::get_elements()[t_ind]) < 1.e-6)
            multiplicities[l] += 1;
      }
    }
  }

  int b_i, b_j;
  int r_i, r_j, delta_r;
  int t_i, t_j, delta_tau;

  G0_original_up.resizeNoCopy(
      std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), b_r_t_VERTEX_dmn_t::dmn_size()));
  G0_original_dn.resizeNoCopy(
      std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), b_r_t_VERTEX_dmn_t::dmn_size()));

  func::dmn_variadic<nu, nu, r_dmn_t, t_VERTEX> G_r_t_dmn;
  for (int j = 0; j < b_r_t_VERTEX_dmn_t::dmn_size(); j++) {
    b_j = fixed_configuration[j].b_ind;
    r_j = fixed_configuration[j].r_ind;
    t_j = fixed_configuration[j].t_ind;

    for (int i = 0; i < b_r_t_VERTEX_dmn_t::dmn_size(); i++) {
      b_i = fixed_configuration[i].b_ind;
      r_i = fixed_configuration[i].r_ind;
      t_i = fixed_configuration[i].t_ind;

      delta_r = RClusterDmn::parameter_type::subtract(r_j, r_i);

      delta_tau = t_i - t_j;

      G0_sign_dn(i, j) = delta_tau < 0 ? -1 : 1;
      G0_sign_up(i, j) = delta_tau < 0 ? -1 : 1;

      if (std::abs(t_VERTEX::get_elements().back() - parameters.get_beta()) < 1.e-6)
        delta_tau = delta_tau < 0 ? delta_tau + t_VERTEX::dmn_size() - 1 : delta_tau;
      else
        delta_tau = delta_tau < 0 ? delta_tau + t_VERTEX::dmn_size() - 0 : delta_tau;

      assert(multiplicities[delta_tau] > 0);

      G0_integration_factor_dn(i, j) = 1. / (r_dmn_t::dmn_size() * multiplicities[delta_tau]);
      G0_integration_factor_up(i, j) = 1. / (r_dmn_t::dmn_size() * multiplicities[delta_tau]);

      G0_indices_dn(i, j) = G_r_t_dmn(b_i, 0, b_j, 0, delta_r, delta_tau);
      G0_indices_up(i, j) = G_r_t_dmn(b_i, 1, b_j, 1, delta_r, delta_tau);
    }
  }
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::initialize_G0_original() {
  int r_ind, b_i, b_j, r_i, r_j;
  scalar_type t_i, t_j, delta_tau;  //, scaled_tau, f_tau;

  G0_original_dn.resizeNoCopy(
      std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), b_r_t_VERTEX_dmn_t::dmn_size()));
  G0_original_up.resizeNoCopy(
      std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), b_r_t_VERTEX_dmn_t::dmn_size()));

  G0_M_G0_matrix_dn.resizeNoCopy(
      std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), b_r_t_VERTEX_dmn_t::dmn_size()));
  G0_M_G0_matrix_up.resizeNoCopy(
      std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), b_r_t_VERTEX_dmn_t::dmn_size()));

  for (int j = 0; j < b_r_t_VERTEX_dmn_t::dmn_size(); j++) {
    b_j = fixed_configuration[j].b_ind;
    r_j = fixed_configuration[j].r_ind;
    t_j = fixed_configuration[j].t_val;

    for (int i = 0; i < b_r_t_VERTEX_dmn_t::dmn_size(); i++) {
      b_i = fixed_configuration[i].b_ind;
      r_i = fixed_configuration[i].r_ind;
      t_i = fixed_configuration[i].t_val;

      r_ind = RClusterDmn::parameter_type::subtract(r_j, r_i);

      delta_tau = t_i - t_j;

      G0_original_dn(i, j) = interpolate_akima(b_i, 0, b_j, 0, r_ind, delta_tau);
      G0_original_up(i, j) = interpolate_akima(b_i, 1, b_j, 1, r_ind, delta_tau);
    }
  }

  //       if(true)
  //        test_G0_original();
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::test_G0_original() {
  for (int i = 0; i < t_VERTEX::dmn_size(); i++) {
    for (int j = 0; j < t_VERTEX::dmn_size(); j++) {
      double t_val = t_VERTEX::get_elements()[i] - t_VERTEX::get_elements()[j];

      std::cout << "\t" << t_val;
    }
    std::cout << "\n";
  }
  std::cout << "\n";

  G_r_t_accumulated = 0;

  for (int j = 0; j < G0_original_dn.size().first; j++)
    for (int i = 0; i < G0_original_dn.size().first; i++)
      G_r_t_accumulated(G0_indices_dn(i, j)) +=
          G0_sign_dn(i, j) * G0_integration_factor_dn(i, j) * G0_original_dn(i, j);

  for (int j = 0; j < G0_original_up.size().first; j++)
    for (int i = 0; i < G0_original_up.size().first; i++)
      G_r_t_accumulated(G0_indices_up(i, j)) +=
          G0_sign_up(i, j) * G0_integration_factor_up(i, j) * G0_original_up(i, j);

  for (int i = 0; i < t_VERTEX::dmn_size(); i++)
    std::cout << "\t" << t_VERTEX::get_elements()[i] << "\t" << G_r_t_accumulated(0, 0, 0, i) << "\n";
  std::cout << std::endl;

  util::Plot::plotLinesPoints(MOMS.G0_r_t_cluster_excluded);

  util::Plot::plotLinesPoints(G_r_t_accumulated);

  interpolate(G_r_t, G_r_t_stddev, spin_XX_chi, spin_ZZ_chi, charge_chi);

  util::Plot::plotLinesPoints(G_r_t);

  G_r_t_accumulated = 0;

  throw std::logic_error(__FUNCTION__);
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::finalize() {
  // util::Plot::plotLinesPoints(G_r_t_accumulated);

  for (int l = 0; l < G_r_t_accumulated_squared.size(); l++)
    G_r_t_accumulated_squared(l) =
        std::sqrt(std::abs(G_r_t_accumulated_squared(l) - std::pow(G_r_t_accumulated(l), 2)));

  interpolate(G_r_t, G_r_t_stddev, spin_XX_chi, spin_ZZ_chi, charge_chi);

  // util::Plot::plotLinesPoints(G_r_t);
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::interpolate(
    func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& G_r_t,
    func::function<double, func::dmn_variadic<nu, nu, r_dmn_t, t>>& G_r_t_stddev,
    func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& spin_XX_chi,
    func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& spin_ZZ_chi,
    func::function<double, func::dmn_variadic<b, b, r_dmn_t, t_sp>>& charge_chi) {
  int size = t_VERTEX::dmn_size();

  math::interpolation::akima_interpolation<double> ai_obj(size);

  double* x = new double[size];
  double* y = new double[size];

  for (int t_ind = 0; t_ind < t_VERTEX::dmn_size(); t_ind++)
    x[t_ind] = t_VERTEX::get_elements()[t_ind];

  {
    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++) {
      for (int nu1_ind = 0; nu1_ind < nu::dmn_size(); nu1_ind++) {
        for (int nu0_ind = 0; nu0_ind < nu::dmn_size(); nu0_ind++) {
          {
            for (int t_ind = 0; t_ind < t_VERTEX::dmn_size(); t_ind++)
              y[t_ind] = G_r_t_accumulated(nu0_ind, nu1_ind, r_ind, t_ind);

            ai_obj.initialize(x, y);

            for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size(); t_ind++)
              G_r_t(nu0_ind, nu1_ind, r_ind, t_ind) = ai_obj.evaluate(t::get_elements()[t_ind]);
          }

          {
            for (int t_ind = 0; t_ind < t_VERTEX::dmn_size(); t_ind++)
              y[t_ind] = G_r_t_accumulated_squared(nu0_ind, nu1_ind, r_ind, t_ind);

            ai_obj.initialize(x, y);

            for (int t_ind = t::dmn_size() / 2; t_ind < t::dmn_size(); t_ind++)
              G_r_t_stddev(nu0_ind, nu1_ind, r_ind, t_ind) =
                  ai_obj.evaluate(t::get_elements()[t_ind]);
          }
        }
      }
    }

    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
      for (int nu1_ind = 0; nu1_ind < nu::dmn_size(); nu1_ind++)
        for (int nu0_ind = 0; nu0_ind < nu::dmn_size(); nu0_ind++)
          for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
            G_r_t(nu0_ind, nu1_ind, r_ind, t_ind) =
                -G_r_t(nu0_ind, nu1_ind, r_ind, t_ind + t::dmn_size() / 2);

    for (int r_ind = 0; r_ind < r_dmn_t::dmn_size(); r_ind++)
      for (int nu1_ind = 0; nu1_ind < nu::dmn_size(); nu1_ind++)
        for (int nu0_ind = 0; nu0_ind < nu::dmn_size(); nu0_ind++)
          for (int t_ind = 0; t_ind < t::dmn_size() / 2; t_ind++)
            G_r_t_stddev(nu0_ind, nu1_ind, r_ind, t_ind) =
                G_r_t_stddev(nu0_ind, nu1_ind, r_ind, t_ind + t::dmn_size() / 2);

    for(int b_j=0; b_j<b::dmn_size(); b_j++){
      for(int b_i=0; b_i<b::dmn_size(); b_i++){
        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
          for(int t_ind=0; t_ind<t_VERTEX::dmn_size(); t_ind++){

            y[t_ind] = spin_XX_chi_accumulated(b_i, b_j, r_ind, t_ind);

            ai_obj.initialize(x, y);

            for(int t_ind=0; t_ind<t_sp::dmn_size(); t_ind++)
              spin_XX_chi(b_i, b_j, r_ind, t_ind) = ai_obj.evaluate(t_sp::get_elements()[t_ind]);
          }
        }
      }
    }

    for(int b_j=0; b_j<b::dmn_size(); b_j++){
      for(int b_i=0; b_i<b::dmn_size(); b_i++){
        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
          for(int t_ind=0; t_ind<t_VERTEX::dmn_size(); t_ind++){

            y[t_ind] = spin_ZZ_chi_accumulated(b_i, b_j, r_ind, t_ind);

            ai_obj.initialize(x, y);

            for(int t_ind=0; t_ind<t_sp::dmn_size(); t_ind++)
              spin_ZZ_chi(b_i, b_j, r_ind, t_ind) = ai_obj.evaluate(t_sp::get_elements()[t_ind]);
          }
        }
      }
    }

    for(int b_j=0; b_j<b::dmn_size(); b_j++){
      for(int b_i=0; b_i<b::dmn_size(); b_i++){
        for(int r_ind=0; r_ind<r_dmn_t::dmn_size(); r_ind++){
          for(int t_ind=0; t_ind<t_VERTEX::dmn_size(); t_ind++){

            y[t_ind] = charge_chi_accumulated(b_i, b_j, r_ind, t_ind);

            ai_obj.initialize(x, y);

            for(int t_ind=0; t_ind<t_sp::dmn_size(); t_ind++)
              charge_chi(b_i, b_j, r_ind, t_ind) = ai_obj.evaluate(t_sp::get_elements()[t_ind]);
          }
        }
      }
    }
  }

  delete[] x;
  delete[] y;
}

template <class parameters_type, class MOMS_type>
template <class configuration_type, typename RealInp>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::compute_G_r_t(
    const configuration_type& configuration_e_up,
    const dca::linalg::Matrix<RealInp, linalg::CPU>& M_up,
    const configuration_type& configuration_e_dn,
    const dca::linalg::Matrix<RealInp, linalg::CPU>& M_dn) {
  {
    int configuration_size = find_first_non_interacting_spin(configuration_e_dn);

    M_matrix_dn.resizeNoCopy(std::pair<int, int>(configuration_size, configuration_size));

    G0_matrix_dn_left.resizeNoCopy(
        std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), configuration_size));
    G0_matrix_dn_right.resizeNoCopy(
        std::pair<int, int>(configuration_size, b_r_t_VERTEX_dmn_t::dmn_size()));

    M_G0_matrix_dn.resizeNoCopy(
        std::pair<int, int>(configuration_size, b_r_t_VERTEX_dmn_t::dmn_size()));

    for (int j = 0; j < configuration_size; j++)
      for (int i = 0; i < configuration_size; i++)
        M_matrix_dn(i, j) = M_dn(i, j);

    GFLOP += 2 * (1.e-9) * b_r_t_VERTEX_dmn_t::dmn_size() * std::pow(configuration_size, 2.);
    GFLOP += 2 * (1.e-9) * configuration_size * std::pow(b_r_t_VERTEX_dmn_t::dmn_size(), 2.);
  }

  {
    int configuration_size = find_first_non_interacting_spin(configuration_e_up);

    M_matrix_up.resizeNoCopy(std::pair<int, int>(configuration_size, configuration_size));

    G0_matrix_up_left.resizeNoCopy(
        std::pair<int, int>(b_r_t_VERTEX_dmn_t::dmn_size(), configuration_size));
    G0_matrix_up_right.resizeNoCopy(
        std::pair<int, int>(configuration_size, b_r_t_VERTEX_dmn_t::dmn_size()));

    M_G0_matrix_up.resizeNoCopy(
        std::pair<int, int>(configuration_size, b_r_t_VERTEX_dmn_t::dmn_size()));

    for (int j = 0; j < configuration_size; j++)
      for (int i = 0; i < configuration_size; i++)
        M_matrix_up(i, j) = M_up(i, j);

    GFLOP += 2. * (1.e-9) * b_r_t_VERTEX_dmn_t::dmn_size() * std::pow(configuration_size, 2.);
    GFLOP += 2. * (1.e-9) * configuration_size * std::pow(b_r_t_VERTEX_dmn_t::dmn_size(), 2.);
  }

  {
    compute_G0_matrix_left(e_DN, configuration_e_dn, G0_matrix_dn_left);
    compute_G0_matrix_left(e_UP, configuration_e_up, G0_matrix_up_left);

    compute_G0_matrix_right(e_DN, configuration_e_dn, G0_matrix_dn_right);
    compute_G0_matrix_right(e_UP, configuration_e_up, G0_matrix_up_right);
  }

  {
    dca::linalg::matrixop::gemm(M_matrix_dn, G0_matrix_dn_right, M_G0_matrix_dn);
    dca::linalg::matrixop::gemm(M_matrix_up, G0_matrix_up_right, M_G0_matrix_up);

    dca::linalg::matrixop::gemm(G0_matrix_dn_left, M_G0_matrix_dn, G0_M_G0_matrix_dn);
    dca::linalg::matrixop::gemm(G0_matrix_up_left, M_G0_matrix_up, G0_M_G0_matrix_up);
  }

  {
    for (int j = 0; j < G0_M_G0_matrix_dn.size().second; j++)
      for (int i = 0; i < G0_M_G0_matrix_dn.size().first; i++)
        G_r_t_dn(i, j) = G0_sign_dn(i, j) * (G0_original_dn(i, j) - G0_M_G0_matrix_dn(i, j));

    for (int j = 0; j < G0_M_G0_matrix_up.size().second; j++)
      for (int i = 0; i < G0_M_G0_matrix_up.size().first; i++)
        G_r_t_up(i, j) = G0_sign_up(i, j) * (G0_original_up(i, j) - G0_M_G0_matrix_up(i, j));
  }
}

template <class parameters_type, class MOMS_type>
//     template<class configuration_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::accumulate_G_r_t(double sign) {
  for (int j = 0; j < b_r_t_VERTEX_dmn_t::dmn_size(); j++) {
    for (int i = 0; i < b_r_t_VERTEX_dmn_t::dmn_size(); i++) {
      G_r_t_accumulated(G0_indices_dn(i, j)) +=
          sign * G0_integration_factor_dn(i, j) * G_r_t_dn(i, j);
      G_r_t_accumulated_squared(G0_indices_dn(i, j)) +=
          sign * G0_integration_factor_dn(i, j) * G_r_t_dn(i, j) * G_r_t_dn(i, j);

      G_r_t_accumulated(G0_indices_up(i, j)) +=
          sign * G0_integration_factor_up(i, j) * G_r_t_up(i, j);
      G_r_t_accumulated_squared(G0_indices_up(i, j)) +=
          sign * G0_integration_factor_up(i, j) * G_r_t_up(i, j) * G_r_t_up(i, j);
    }
  }
}

/*!
 *   <S_z> = (n_up-1/2)*(n_dn-1/2)
 */
template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::accumulate_moments(double sign) {
  for (int b_ind = 0; b_ind < b::dmn_size(); b_ind++) {
    for (int r_i = 0; r_i < r_dmn_t::dmn_size(); r_i++) {
      for (int t_ind = 0; t_ind < t_VERTEX::dmn_size(); t_ind++) {
        int i = b_r_t_dmn(b_ind, r_i, t_ind);
        int j = i;

        double charge_val = G_r_t_up(i, j) * G_r_t_dn(i, j);  // double occupancy = <n_d*n_u>
        double magnetic_val =
            1. - 2. * G_r_t_up(i, j) * G_r_t_dn(i, j);  // <m^2> = 1-2*<n_d*n_u> (T. Paiva, PRB 2001)

        charge_cluster_moment(b_ind, r_i) += sign * charge_val / t_VERTEX::dmn_size();
        magnetic_cluster_moment(b_ind, r_i) += sign * magnetic_val / t_VERTEX::dmn_size();
      }
    }
  }
}

/*
 *  spin_XX_Chi(r-r',tau-tau') = 1/2 * <S+(r,t)S-(r',t') + S-(r',t')S+(r,t)>  (this symmetrized form to average over the spin)
 *                             = - G_up(r-r',tau-tau') * G_down(r'-r,tau'-tau) -G_up(r'-r,tau-tau') * G_down(r-r',tau'-tau)
 *  spin_ZZ_Chi is almost the same as charge_Chi, except two sign change as below
 *  charge_Chi(r-r',tau-tau') = sum_s,s' <n(r,t,s)n(r',t',s')>

 *  Note there are (no) disconnected diagrams for charge(spin)_chi
 *  Chi(0,0) needs correction due to cc+ = 1-c+c
 *  Further FT to chi(q,iv) can be done by an external python script, although this way we do not have errorbar

 *  Note that below we use Chi(t) ~ Gdn(i,j)*Gup(j,i) instead of -Gdn(i,j)*Gup(j,i) for positive Chi values
 *  since G(t>0)<0 and G(t<0)>0
 *  charge(spin)_susceptibility_accumulated will be interpolated to sp time_intervals
 */
template<class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::accumulate_chi(double sign){

  int b_i, b_j, r_i, r_j, t_i, t_j, dr, dt;
  double upup, updn;
  double sfactor = 0.5*sign/((t_VERTEX::dmn_size()-1)*r_dmn_t::dmn_size());
  double cfactor = 2*sfactor;

  for(int j=0; j<b_r_t_VERTEX_dmn_t::dmn_size(); j++){
    b_j = fixed_configuration[j].b_ind;
    r_j = fixed_configuration[j].r_ind;
    t_j = fixed_configuration[j].t_ind;

    for(int i=0; i<b_r_t_VERTEX_dmn_t::dmn_size(); i++){
      b_i = fixed_configuration[i].b_ind;
      r_i = fixed_configuration[i].r_ind;
      t_i = fixed_configuration[i].t_ind;

      dr = RClusterDmn::parameter_type::subtract(r_j, r_i);
      dt = t_i-t_j;

      // chi(beta) considered seperately later
      if(t_i != t_VERTEX::dmn_size()-1 && t_j != t_VERTEX::dmn_size()-1)
      {
        dt = dt<0 ? dt+t_VERTEX::dmn_size()-1 : dt;

        // connected diagrams:
        // chi(0) ~ Gdn(0+)*Gup(0-) while Gup(0-) still has positive sign from G0_sign_up(i,j), so change sign here
        upup = G_r_t_up(i,j)*G_r_t_up(j,i) + G_r_t_dn(i,j)*G_r_t_dn(j,i);
        updn = G_r_t_dn(i,j)*G_r_t_up(j,i) + G_r_t_up(i,j)*G_r_t_dn(j,i);

        if(dt==0){
          spin_XX_chi_accumulated(b_i,b_j,dr,dt) -= sfactor* updn;
          spin_ZZ_chi_accumulated(b_i,b_j,dr,dt) -= sfactor* upup;
          charge_chi_accumulated (b_i,b_j,dr,dt) -= cfactor* upup;
        }
        else{
          spin_XX_chi_accumulated(b_i,b_j,dr,dt) += sfactor* updn;
          spin_ZZ_chi_accumulated(b_i,b_j,dr,dt) += sfactor* upup;
          charge_chi_accumulated (b_i,b_j,dr,dt) += cfactor* upup;
        }

        // disconnected diagrams:
        // note that (1-G_sigma)(1-G_sigma') switch to (1+G_sigma)(1+G_sigma') since G(dt=0)<0 needs changing sign
        upup = (1.0+G_r_t_up(i,i))*(1.0+G_r_t_up(j,j)) + (1.0+G_r_t_dn(i,i))*(1.0+G_r_t_dn(j,j));
        updn = (1.0+G_r_t_up(i,i))*(1.0+G_r_t_dn(j,j)) + (1.0+G_r_t_dn(i,i))*(1.0+G_r_t_up(j,j));
 
        spin_ZZ_chi_accumulated(b_i,b_j,dr,dt) += sfactor* (upup - updn);
        charge_chi_accumulated (b_i,b_j,dr,dt) += cfactor* (upup + updn);

        // correction due to cc+ = 1=c+c
        updn = G_r_t_up(j,j) + G_r_t_dn(j,j);
        if(b_i==b_j && dr==0 && dt==0){
          spin_XX_chi_accumulated(b_i,b_j,0,0) -= sfactor* updn;
          spin_ZZ_chi_accumulated(b_i,b_j,0,0) -= sfactor* updn;
          charge_chi_accumulated (b_i,b_j,0,0) -= cfactor* updn;
        }
      }
 
      // chi(beta) = chi(0)
      spin_XX_chi_accumulated(b_i,b_j,dr,t_VERTEX::dmn_size()-1) = spin_XX_chi_accumulated(b_i,b_j,dr,0);
      spin_ZZ_chi_accumulated(b_i,b_j,dr,t_VERTEX::dmn_size()-1) = spin_ZZ_chi_accumulated(b_i,b_j,dr,0);
      charge_chi_accumulated (b_i,b_j,dr,t_VERTEX::dmn_size()-1) = charge_chi_accumulated (b_i,b_j,dr,0);
    }
  }
}

/*
 * P_d
 */
template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::accumulate_dwave_pp_correlator(double sign) {
  double renorm = 1. / (t_VERTEX::dmn_size() * pow(r_dmn_t::dmn_size(), 2.));
  double factor = sign * renorm;

  for (int r_i = 0; r_i < r_dmn_t::dmn_size(); r_i++) {
    for (int r_j = 0; r_j < r_dmn_t::dmn_size(); r_j++) {
      for (int r_l = 0; r_l < r_dmn_t::dmn_size(); r_l++) {
        int l_minus_i = RClusterDmn::parameter_type::subtract(r_i, r_l);
        int l_minus_j = RClusterDmn::parameter_type::subtract(r_j, r_l);

        double struct_factor = dwave_r_factor(l_minus_i) * dwave_r_factor(l_minus_j);

        if (std::abs(struct_factor) > 1.e-6) {
          for (int b_i = 0; b_i < b::dmn_size(); b_i++) {
            for (int b_j = 0; b_j < b::dmn_size(); b_j++) {
              for (int b_l = 0; b_l < b::dmn_size(); b_l++) {
                double value = 0;

                for (int t_ind = 0; t_ind < t_VERTEX::dmn_size(); t_ind++) {
                  int i = b_r_t_dmn(b_i, r_i, t_ind);
                  int j = b_r_t_dmn(b_j, r_j, t_ind);
                  int l = b_r_t_dmn(b_l, r_l, t_ind);

                  double d_ij = i == j ? 1 : 0;
                  double d_il = i == l ? 1 : 0;
                  double d_lj = l == j ? 1 : 0;
                  double d_ll = 1;  // l==l? 1 : 0;

                  value += (d_ij - G_r_t_up(j, i)) * (d_ll - G_r_t_dn(l, l));
                  value += (d_ij - G_r_t_dn(j, i)) * (d_ll - G_r_t_up(l, l));

                  value += (d_il - G_r_t_up(l, i)) * (d_lj - G_r_t_dn(j, l));
                  value += (d_il - G_r_t_dn(l, i)) * (d_lj - G_r_t_up(j, l));
                }

                dwave_pp_correlator(b_l, r_l) += factor * struct_factor * value;
              }
            }
          }
        }
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
int TpEqualTimeAccumulator<parameters_type, MOMS_type>::find_first_non_interacting_spin(
    const std::vector<vertex_singleton_type>& configuration_e_spin) {
  int configuration_size = configuration_e_spin.size();

  int vertex_index = 0;
  while (vertex_index < configuration_size &&
         configuration_e_spin[vertex_index].get_HS_spin() != HS_ZERO)
    vertex_index++;

  assert(vertex_index == configuration_size ||
         configuration_e_spin[vertex_index].get_HS_spin() == HS_ZERO);

  return vertex_index;
}

template <class parameters_type, class MOMS_type>
template <class configuration_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::compute_G0_matrix(
    e_spin_states e_spin, const configuration_type& configuration,
    dca::linalg::Matrix<float, dca::linalg::CPU>& G0_matrix) {
  int spin_index = domains::electron_spin_domain::to_coordinate(e_spin);

  int r_ind, b_i, b_j, r_i, r_j;    //, s_i, s_j;
  scalar_type t_i, t_j, delta_tau;  //, scaled_tau, f_tau;

  int configuration_size = find_first_non_interacting_spin(configuration);
  for (int j = 0; j < configuration_size; j++) {
    const vertex_singleton_type& configuration_e_spin_j = configuration[j];

    b_j = configuration_e_spin_j.get_band();
    r_j = configuration_e_spin_j.get_r_site();
    t_j = configuration_e_spin_j.get_tau();

    for (int i = 0; i < b_r_t_VERTEX_dmn_t::dmn_size(); i++) {
      b_i = fixed_configuration[i].b_ind;
      r_i = fixed_configuration[i].r_ind;
      t_i = fixed_configuration[i].t_val;

      r_ind = RClusterDmn::parameter_type::subtract(r_j, r_i);

      delta_tau = t_i - t_j;

      G0_matrix(i, j) = interpolate_akima(b_i, spin_index, b_j, spin_index, r_ind, delta_tau);
    }
  }
}

template <class parameters_type, class MOMS_type>
template <class configuration_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::compute_G0_matrix_left(
    e_spin_states e_spin, const configuration_type& configuration,
    dca::linalg::Matrix<float, dca::linalg::CPU>& G0_matrix) {
  int spin_index = domains::electron_spin_domain::to_coordinate(e_spin);

  int r_ind, b_i, b_j, r_i, r_j;    //, s_i, s_j;
  scalar_type t_i, t_j, delta_tau;  //, scaled_tau, f_tau;

  int configuration_size = find_first_non_interacting_spin(configuration);
  for (int j = 0; j < configuration_size; j++) {
    const vertex_singleton_type& configuration_e_spin_j = configuration[j];

    b_j = configuration_e_spin_j.get_band();
    r_j = configuration_e_spin_j.get_r_site();
    t_j = configuration_e_spin_j.get_tau();

    for (int i = 0; i < b_r_t_VERTEX_dmn_t::dmn_size(); i++) {
      b_i = fixed_configuration[i].b_ind;
      r_i = fixed_configuration[i].r_ind;
      t_i = fixed_configuration[i].t_val;

      r_ind = RClusterDmn::parameter_type::subtract(r_j, r_i);

      delta_tau = t_i - t_j;

      G0_matrix(i, j) = interpolate_akima(b_i, spin_index, b_j, spin_index, r_ind, delta_tau);
    }
  }
}

template <class parameters_type, class MOMS_type>
template <class configuration_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::compute_G0_matrix_right(
    e_spin_states e_spin, const configuration_type& configuration,
    dca::linalg::Matrix<float, dca::linalg::CPU>& G0_matrix) {
  int spin_index = domains::electron_spin_domain::to_coordinate(e_spin);

  int r_ind, b_i, b_j, r_i, r_j;    //, s_i, s_j;
  scalar_type t_i, t_j, delta_tau;  //, scaled_tau, f_tau;

  int configuration_size = find_first_non_interacting_spin(configuration);

  for (int j = 0; j < b_r_t_VERTEX_dmn_t::dmn_size(); j++) {
    b_j = fixed_configuration[j].b_ind;
    r_j = fixed_configuration[j].r_ind;
    t_j = fixed_configuration[j].t_val;

    for (int i = 0; i < configuration_size; i++) {
      const vertex_singleton_type& configuration_e_spin_i = configuration[i];

      b_i = configuration_e_spin_i.get_band();
      r_i = configuration_e_spin_i.get_r_site();
      t_i = configuration_e_spin_i.get_tau();

      r_ind = RClusterDmn::parameter_type::subtract(r_j, r_i);

      delta_tau = t_i - t_j;

      G0_matrix(i, j) = interpolate_akima(b_i, spin_index, b_j, spin_index, r_ind, delta_tau);
    }
  }
}

template <class parameters_type, class MOMS_type>
inline double TpEqualTimeAccumulator<parameters_type, MOMS_type>::interpolate_akima(
    int b_i, int s_i, int b_j, int s_j, int delta_r, double tau) {
  const static double beta = parameters.get_beta();
  const static double N_div_beta = parameters.get_sp_time_intervals() / beta;

  int sign = 1;
  // Map tau to [0, beta).
  if (tau < 0) {
    tau += beta;
    sign = -1;
  }
  assert(0 <= tau && tau < beta);

  const double scaled_tau = tau * N_div_beta;
  // Find interpolation index of on the left of tau.
  const int t_ind = static_cast<int>(scaled_tau);

#ifndef NDEBUG
  const double* positive_times =
      shifted_t::get_elements().data() + shifted_t::get_elements().size() / 2;
  assert(positive_times[t_ind] <= tau && tau < positive_times[t_ind] + 1. / N_div_beta);
#endif  // NDEBUG

  const double delta_tau = scaled_tau - t_ind;
  assert(delta_tau >= 0 && delta_tau < 1);

  const int linind = 4 * nu_nu_r_dmn_t_t_shifted_dmn(b_i, s_i, b_j, s_j, delta_r, t_ind);

  const double* a_ptr = &akima_coefficients(linind);

  const double result =
      (a_ptr[0] + delta_tau * (a_ptr[1] + delta_tau * (a_ptr[2] + delta_tau * a_ptr[3])));

  return sign * result;
}

template <class parameters_type, class MOMS_type>
template <class configuration_type, typename RealInp>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::accumulateAll(
    const configuration_type& configuration_e_up,
    const dca::linalg::Matrix<RealInp, dca::linalg::CPU>& M_up,
    const configuration_type& configuration_e_dn,
    const dca::linalg::Matrix<RealInp, dca::linalg::CPU>& M_dn, int sign) {
  compute_G_r_t(configuration_e_up, M_up, configuration_e_dn, M_dn);

  accumulate_G_r_t(sign);

  accumulate_moments(sign);

  accumulate_chi(sign);

  accumulate_dwave_pp_correlator(sign);
}

template <class parameters_type, class MOMS_type>
void TpEqualTimeAccumulator<parameters_type, MOMS_type>::sumTo(
    dca::phys::solver::ctaux::TpEqualTimeAccumulator<parameters_type, MOMS_type>& other) const {
  other.G_r_t_accumulated += G_r_t_accumulated;
  other.G_r_t_accumulated_squared += G_r_t_accumulated_squared;
  other.charge_cluster_moment += charge_cluster_moment;
  other.magnetic_cluster_moment += magnetic_cluster_moment;
  other.spin_XX_chi_accumulated += spin_XX_chi_accumulated;
  other.spin_ZZ_chi_accumulated += spin_ZZ_chi_accumulated;
  other.charge_chi_accumulated += charge_chi_accumulated;
  other.dwave_pp_correlator += dwave_pp_correlator;
  other.GFLOP += GFLOP;
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_ACCUMULATOR_TP_TP_EQUAL_TIME_ACCUMULATOR_HPP

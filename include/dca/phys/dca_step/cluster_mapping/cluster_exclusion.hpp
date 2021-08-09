// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class performs the cluster exclusion step.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_CLUSTER_EXCLUSION_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_CLUSTER_EXCLUSION_HPP

#include <complex>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/linalg.hpp"
#include "dca/math/function_transform/function_transform.hpp"
#include "dca/phys/domains/cluster/cluster_domain.hpp"
#include "dca/phys/domains/quantum/electron_band_domain.hpp"
#include "dca/phys/domains/quantum/electron_spin_domain.hpp"
#include "dca/phys/domains/time_and_frequency/frequency_domain.hpp"
#include "dca/phys/domains/time_and_frequency/time_domain.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/util/plot.hpp"
#include "dca/phys/dca_step/symmetrization/symmetrize.hpp"
#include "dca/phys/dca_data/dca_data.hpp"

namespace dca {
namespace phys {
namespace clustermapping {
// dca::phys::clustermapping::

template <typename parameters_type, typename MOMS_type>
class cluster_exclusion {
public:
  using profiler_type = typename parameters_type::profiler_type;
  using concurrency_type = typename parameters_type::concurrency_type;

  using t = func::dmn_0<domains::time_domain>;
  using w = func::dmn_0<domains::frequency_domain>;

  using b = func::dmn_0<domains::electron_band_domain>;
  using s = func::dmn_0<domains::electron_spin_domain>;
  using nu = func::dmn_variadic<b, s>;  // orbital-spin index

  using CDA = ClusterDomainAliases<parameters_type::lattice_type::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;
  using KClusterDmn = typename CDA::KClusterDmn;
  using NuDmn = func::dmn_variadic<b, s>;  // orbital-spin index
  using NuNuDmn = func::dmn_variadic<NuDmn, NuDmn>;

public:
  cluster_exclusion(parameters_type& parameters_ref, MOMS_type& MOMS_ref);

  void execute();
  func::function<int, NuNuDmn> H_symmetry;
private:
  void compute_G0_K_w_cluster_excluded();
  void compute_G0_R_t_cluster_excluded();

  void plot_G0_R_t_cluster_excluded();

private:
  parameters_type& parameters;
  MOMS_type& MOMS;
};

template <typename parameters_type, typename MOMS_type>
cluster_exclusion<parameters_type, MOMS_type>::cluster_exclusion(parameters_type& parameters_ref,
                                                                 MOMS_type& MOMS_ref)
    : parameters(parameters_ref), MOMS(MOMS_ref) {}

template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::execute() {
  
  compute_G0_K_w_cluster_excluded();
  compute_G0_R_t_cluster_excluded();
}

/*
 *   G = G_0 + G_0*S*G = G_0 * (1 + S*G)
 *
 *   G_0 = G*(1 + S*G)^-1
 */
template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::compute_G0_K_w_cluster_excluded() {
  profiler_type profiler(__FUNCTION__, "cluster_exclusion", __LINE__);

  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G_matrix("G_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> S_matrix("S_matrix", nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> one_plus_S_G("one_plus_S_G",
                                                                           nu::dmn_size());
  dca::linalg::Matrix<std::complex<double>, dca::linalg::CPU> G0_matrix("G0_matrix", nu::dmn_size());

  // Allocate the work space for inverse only once.
  dca::linalg::Vector<int, dca::linalg::CPU> ipiv;
  dca::linalg::Vector<std::complex<double>, dca::linalg::CPU> work;

  for (int w_ind = 0; w_ind < w::dmn_size(); w_ind++) {
    for (int K_ind = 0; K_ind < KClusterDmn::dmn_size(); K_ind++) {
      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G_matrix(i, j) = MOMS.G_k_w(i, j, K_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          S_matrix(i, j) = MOMS.Sigma_cluster(i, j, K_ind, w_ind);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          one_plus_S_G(i, j) = 0;

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          for (int l = 0; l < nu::dmn_size(); l++)
            one_plus_S_G(i, j) += S_matrix(i, l) * G_matrix(l, j);

      for (int i = 0; i < nu::dmn_size(); i++)
        one_plus_S_G(i, i) += 1.;

      dca::linalg::matrixop::inverse(one_plus_S_G, ipiv, work);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          G0_matrix(i, j) = 0;

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          for (int l = 0; l < nu::dmn_size(); l++)
            G0_matrix(i, j) += G_matrix(i, l) * one_plus_S_G(l, j);

      for (int j = 0; j < nu::dmn_size(); j++)
        for (int i = 0; i < nu::dmn_size(); i++)
          MOMS.G0_k_w_cluster_excluded(i, j, K_ind, w_ind) = G0_matrix(i, j);
    }
  }
}

template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::compute_G0_R_t_cluster_excluded() {
  profiler_type profiler(__FUNCTION__, "cluster_exclusion", __LINE__);
  MOMS.G0_k_w_cluster_excluded -= MOMS.G0_k_w;
  H_symmetry = -1;
  {
    math::transform::FunctionTransform<w, t>::execute(MOMS.G0_k_w_cluster_excluded,
                                                      MOMS.G0_k_t_cluster_excluded);

    MOMS.G0_k_t_cluster_excluded += MOMS.G0_k_t;
    
    symmetrize::execute(MOMS.G0_k_t_cluster_excluded, H_symmetry, true);

    math::transform::FunctionTransform<KClusterDmn, RClusterDmn>::execute(MOMS.G0_k_t_cluster_excluded,
                                                              MOMS.G0_r_t_cluster_excluded);
  }

  MOMS.G0_k_w_cluster_excluded += MOMS.G0_k_w;
}

template <typename parameters_type, typename MOMS_type>
void cluster_exclusion<parameters_type, MOMS_type>::plot_G0_R_t_cluster_excluded() {
  {
    func::function<float, t> tmp("G0_k_t");

    util::Plot plot("lines");
    for (int R_ind = 0; R_ind < RClusterDmn::dmn_size(); R_ind++) {
      for (int t_ind = 0; t_ind < t::dmn_size(); t_ind++)
        tmp(t_ind) = MOMS.G0_k_t(0, 0, R_ind, t_ind);

      plot.plot(tmp);
    }
  }

  {
    func::function<float, t> tmp("G0_k_t_cluster_excluded");

    util::Plot plot("lines");
    for (int R_ind = 0; R_ind < RClusterDmn::dmn_size(); R_ind++) {
      for (int t_ind = 0; t_ind < t::dmn_size(); t_ind++)
        tmp(t_ind) = MOMS.G0_k_t_cluster_excluded(0, 0, R_ind, t_ind);

      plot.plot(tmp);
    }
  }

  {
    func::function<float, t> tmp("G0_r_t");

    util::Plot plot("lines");
    for (int R_ind = 0; R_ind < RClusterDmn::dmn_size(); R_ind++) {
      for (int t_ind = 0; t_ind < t::dmn_size(); t_ind++)
        tmp(t_ind) = MOMS.G0_r_t_cluster_excluded(0, 0, R_ind, t_ind);

      plot.plot(tmp);
    }
  }
}

}  // clustermapping
}  // phys
}  // dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_MAPPING_CLUSTER_EXCLUSION_HPP

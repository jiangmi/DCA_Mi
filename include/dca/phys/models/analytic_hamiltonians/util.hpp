// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//
// This file provides utility functions for tasks that are common between different lattices.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_UTIL_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_UTIL_HPP

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/cluster_operations.hpp"

namespace dca {
namespace phys {
namespace models {
namespace util {
// dca::phys::models::detail::

// Initializes the interaction part of the single-band real space Hubbard Hamiltonian with on-site
// and nearest-neighbor interaction.
// nn_vec contains the vectors that define the nearest neighbors.
// In: parameters
//     nn_vec
// Out: H_int
template <typename ParametersType, typename BandDmn, typename SpinDmn, typename RDmn>
void initializeSingleBandHint(
    const ParametersType& parameters,
    const std::vector<typename RDmn::parameter_type::element_type>& nn_vec,
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_int) {
  if (BandDmn::dmn_size() != 1)
    throw std::logic_error("Band domain size must be 1.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  // Get the index of the origin (0,0).
  const int origin = RDmn::parameter_type::origin_index();

  // Compute indices of nearest neighbors (nn) w.r.t. origin.
  std::vector<int> nn_index;

  const std::vector<typename RDmn::parameter_type::element_type>& super_basis =
      RDmn::parameter_type::get_super_basis_vectors();
  const std::vector<typename RDmn::parameter_type::element_type>& elements =
      RDmn::parameter_type::get_elements();

  for (const auto& vec : nn_vec) {
    std::vector<double> nn_vec_translated =
        domains::cluster_operations::translate_inside_cluster(vec, super_basis);
    int tmp = domains::cluster_operations::index(nn_vec_translated, elements, domains::BRILLOUIN_ZONE);
    nn_index.push_back(tmp);

    const int minus_r = RDmn::parameter_type::subtract(tmp, origin);
    nn_index.push_back(minus_r);

    // debug:
  //  for (auto index : nn_index) {
  //    std::cout << "nn_index = " << index << "\tvec = " << RDmn::get_elements()[index][0] << "\t" << RDmn::get_elements()[index][1] << "\n";
  //  }
  }

  // Set all elements to zero.
  H_int = 0.;

  // Nearest-neighbor opposite spin interaction
  const double V = parameters.get_V();
  for (auto index : nn_index) {
    H_int(0, 0, 0, 1, index) = V;
    H_int(0, 1, 0, 0, index) = V;
  }

  // Nearest-neighbor same spin interaction
  const double V_prime = parameters.get_V_prime();
  for (auto index : nn_index) {
    H_int(0, 0, 0, 0, index) = V_prime;
    H_int(0, 1, 0, 1, index) = V_prime;
  }

  // On-site interaction
  // This has to be set last since for small clusters a nearest neighbor might
  // be the same site and therefore V would overwrite U.
  const double U = parameters.get_U();
  H_int(0, 0, 0, 1, origin) = U;
  H_int(0, 1, 0, 0, origin) = U;
}


// Initializes the interaction part of the three-band real space Hubbard Hamiltonian with 
// nearest-neighbor interaction Vpp
// Note that here not only nn but also nnn (-1,1) separation is needed for Vpp
// nn_vec contains the vectors that define the nearest neighbors.
// In: parameters
//     nn_vec
// Out: H_int
template <typename ParametersType, typename BandDmn, typename SpinDmn, typename RDmn>
void initialize3BandHint(
    const ParametersType& parameters,
    const std::vector<typename RDmn::parameter_type::element_type>& nn_vec,
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_int) {
  if (BandDmn::dmn_size() != 3)
    throw std::logic_error("Band domain size must be 3.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  int BANDS = 3;
      
  // Get the index of the origin (0,0).
  const int origin = RDmn::parameter_type::origin_index();

  const std::vector<typename RDmn::parameter_type::element_type>& super_basis =
      RDmn::parameter_type::get_super_basis_vectors();
  const std::vector<typename RDmn::parameter_type::element_type>& elements =
      RDmn::parameter_type::get_elements();

  // debug:
  //for (int i=0; i<RDmn::dmn_size();i++) {
  //  std::cout << "rdmn " << i << "\t" << RDmn::get_elements()[i][0] << "\t" << RDmn::get_elements()[i][1] << "\n";
  //}
    
  // Set all elements to zero.
  H_int = 0.;
    
  // Nearest-neighbor opposite spin interaction between px and py
  const double V_pp = parameters.get_V_pp();
    
  // Nearest-neighbor same spin interaction between px and py
  const double V_pp_prime = parameters.get_V_pp_prime();
    
  // set V_pp within unit cell
  for (int i = 1; i < BANDS; i++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int j = 1; j < BANDS; j++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (i != j && s1 != s2)
            H_int(i, s1, j, s2, origin) = V_pp;

          if (i != j && s1 == s2)
            H_int(i, s1, j, s2, origin) = V_pp_prime;
        }
      }
    }
  }
    
  // set V_pp between unit cell with (1,0) separation between (px, py) orbitals
  std::vector<double> nn_vec_translated0 =
      domains::cluster_operations::translate_inside_cluster(nn_vec[0], super_basis);

  int i0 = domains::cluster_operations::index(nn_vec_translated0, elements, domains::BRILLOUIN_ZONE);
  int minus_i0 = RDmn::parameter_type::subtract(i0, origin);
    
  //std::cout << "i0 " << i0 << "\t" << RDmn::get_elements()[i0][0] << "\t" << RDmn::get_elements()[i0][1] << "\n";

  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      if (s1 != s2) {
        H_int(1, s1, 2, s2, i0) = V_pp;
        H_int(2, s2, 1, s1, minus_i0) = V_pp;
      }
      else if (s1 == s2) {
        H_int(1, s1, 2, s2, i0) = V_pp_prime;
        H_int(2, s2, 1, s1, minus_i0) = V_pp_prime;
      }
    }
  }
    
    
  // set V_pp between unit cell with (0,1) separation between (py, px) orbitals
  std::vector<double> nn_vec_translated1 =
      domains::cluster_operations::translate_inside_cluster(nn_vec[1], super_basis);
  int i1 = domains::cluster_operations::index(nn_vec_translated1, elements, domains::BRILLOUIN_ZONE);
  int minus_i1 = RDmn::parameter_type::subtract(i1, origin);

  //std::cout << "i1 " << i1 << "\t" << RDmn::get_elements()[i1][0] << "\t" << RDmn::get_elements()[i1][1] << "\n";

  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      if (s1 != s2) {
        H_int(2, s1, 1, s2, i1) = V_pp;
        H_int(1, s2, 2, s1, minus_i1) = V_pp;
      }
      else if (s1 == s2) {
        H_int(2, s1, 1, s2, i1) = V_pp_prime;
        H_int(1, s2, 2, s1, minus_i1) = V_pp_prime;
      }
    }
  }
    
    
  // set V_pp between unit cell with (-1,1) separation between (py, px) orbitals
  //                             and (1,-1) separation between (px, py) orbitals
  std::vector<double> nn_vec_translated2 =
      domains::cluster_operations::translate_inside_cluster(nn_vec[2], super_basis);
  int i2 = domains::cluster_operations::index(nn_vec_translated2, elements, domains::BRILLOUIN_ZONE);
  int minus_i2 = RDmn::parameter_type::subtract(i2, origin);

  //std::cout << "i2 " << i2 << "\t" << RDmn::get_elements()[i2][0] << "\t" << RDmn::get_elements()[i2][1] << "\n";

  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      if (s1 != s2) {
        H_int(2, s1, 1, s2, i2) = V_pp;
        H_int(1, s2, 2, s1, minus_i2) = V_pp;
      }
      else if (s1 == s2) {
        H_int(2, s1, 1, s2, i2) = V_pp_prime;
        H_int(1, s2, 2, s1, minus_i2) = V_pp_prime;
      }
    }
  }
    
    
  // On-site interaction
  // This has to be set last since for small clusters a nearest neighbor might
  // be the same site and therefore V would overwrite U.
  const double U_dd = parameters.get_U_dd();  // interaction in d band
  const double U_pp = parameters.get_U_pp();  // interaction in p bands
  
  for (int i = 0; i < BANDS; i++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int j = 0; j < BANDS; j++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (i == 0 && j == 0 && s1 != s2)
            H_int(i, s1, j, s2, origin) = U_dd;

          if (i == j && i != 0 && s1 != s2)
            H_int(i, s1, j, s2, origin) = U_pp;
        }
      }
    }
  }
}

}  // util
}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_UTIL_HPP

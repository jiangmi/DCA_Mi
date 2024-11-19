// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang (Nov. 18, 2024)
// 2x2x2 lattice for two-orbital model. Use 2D lattice to simulate realistic 3D lattice
// so that need the hopping between 1st and 4th layer to fulfil the PBC along z direction
//
// Bilayer lattice.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_EG_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_EG_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::

template <typename point_group_type>
class bilayer_eg {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
 //  typedef point_group_type DCA_point_group;

  // Aug.17, 2023 debug
  // for C4-C2 rotational sym breaking cases:
  typedef domains::no_symmetry<2> DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 4;

  static double* initialize_r_DCA_basis();
  static double* initialize_k_DCA_basis();

  static double* initialize_r_LDA_basis();
  static double* initialize_k_LDA_basis();

  static std::vector<int> get_flavors();
  static std::vector<std::vector<double>> get_a_vectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> get_orbital_permutations();

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initialize_H_interaction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class domain>
  static void initialize_H_symmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are two-dimensional (access through index 0 and 1).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initialize_H_0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);
};

template <typename point_group_type>
double* bilayer_eg<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
double* bilayer_eg<point_group_type>::initialize_k_DCA_basis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
double* bilayer_eg<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
double* bilayer_eg<point_group_type>::initialize_k_LDA_basis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> bilayer_eg<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> bilayer_eg<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> bilayer_eg<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void bilayer_eg<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer eg has four bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U1 = parameters.get_U1();              // Same band, opposite spin.
  const double U2 = parameters.get_U2();
  const double V  = parameters.get_V();              // Different band, opposite spin.
  const double V_prime = parameters.get_V_prime();  // Different band, same spin.

  H_interaction = 0.;

  H_interaction(0, 0, 0, 1, origin) = U1;
  H_interaction(0, 1, 0, 0, origin) = U1;
  H_interaction(1, 0, 1, 1, origin) = U2;
  H_interaction(1, 1, 1, 0, origin) = U2;
  H_interaction(2, 0, 2, 1, origin) = U1;
  H_interaction(2, 1, 2, 0, origin) = U1;
  H_interaction(3, 0, 3, 1, origin) = U2;
  H_interaction(3, 1, 3, 0, origin) = U2;

  // intralayer's inter-orbital Hund's J
  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      if (s1 != s2) {
        H_interaction(0, s1, 1, s2, origin) = V;
        H_interaction(1, s1, 0, s2, origin) = V;
        H_interaction(2, s1, 3, s2, origin) = V;
        H_interaction(3, s1, 2, s2, origin) = V;
      }
      if (s1 == s2) {
        H_interaction(0, s1, 1, s2, origin) = V_prime;
        H_interaction(1, s1, 0, s2, origin) = V_prime;
        H_interaction(2, s1, 3, s2, origin) = V_prime;
        H_interaction(3, s1, 2, s2, origin) = V_prime;
      }
    }
  }
}

template <typename point_group_type>
template <class domain>
void bilayer_eg<point_group_type>::initialize_H_symmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

/*  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;
*/
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void bilayer_eg<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto e1 = parameters.get_e1();
  const auto e2 = parameters.get_e2();
  const auto t1 = parameters.get_t1();
  const auto t2 = parameters.get_t2();
  const auto t1_prime = parameters.get_t1_prime();
  const auto t2_prime = parameters.get_t2_prime();
    
  // intralayer inter-orbital 
  const auto t_hyb = parameters.get_t_hyb();
  const auto t_hyb_prime = parameters.get_t_hyb_prime();
    
  // interlayer
  const auto t_perp1 = parameters.get_t_perp1();
  const auto t_perp2 = parameters.get_t_perp2();

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
      
    // the follwing do not use conventional -t (which is used for an electron-like band)
    // here directly use DFT parameters (with electron language), which are normally negative values
    const auto val1 =
        e1 +2. * t1 * (std::cos(k[0]) + std::cos(k[1])) + 4. * t1_prime * std::cos(k[0]) * std::cos(k[1]);
    const auto val2 =
        e2 +2. * t2 * (std::cos(k[0]) + std::cos(k[1])) + 4. * t2_prime * std::cos(k[0]) * std::cos(k[1]);
      
    // intralayer inter-orbital 
    const auto val3 = t_hyb + 2. * t_hyb_prime * (std::cos(k[0]) + std::cos(k[1]));
      
    // interlayer for two different orbitals
    const auto val41 = t_perp1;
    const auto val42 = t_perp2;

    H_0(0, 0, 0, 0, k_ind) = val1;
    H_0(0, 1, 0, 1, k_ind) = val1;
    H_0(1, 0, 1, 0, k_ind) = val2;
    H_0(1, 1, 1, 1, k_ind) = val2;
    H_0(2, 0, 2, 0, k_ind) = val1;
    H_0(2, 1, 2, 1, k_ind) = val1;
    H_0(3, 0, 3, 0, k_ind) = val2;
    H_0(3, 1, 3, 1, k_ind) = val2;

    // inter-orbital nn hybridization
    H_0(0, 0, 1, 0, k_ind) = val3;
    H_0(0, 1, 1, 1, k_ind) = val3;
    H_0(1, 0, 0, 0, k_ind) = val3;
    H_0(1, 1, 0, 1, k_ind) = val3;
    H_0(2, 0, 3, 0, k_ind) = val3;
    H_0(2, 1, 3, 1, k_ind) = val3;
    H_0(3, 0, 2, 0, k_ind) = val3;
    H_0(3, 1, 2, 1, k_ind) = val3;
      
    // inter-layer same orbital hybridization
    H_0(0, 0, 2, 0, k_ind) = val41;
    H_0(0, 1, 2, 1, k_ind) = val41;
    H_0(2, 0, 0, 0, k_ind) = val41;
    H_0(2, 1, 0, 1, k_ind) = val41;
    H_0(1, 0, 3, 0, k_ind) = val42;
    H_0(1, 1, 3, 1, k_ind) = val42;
    H_0(3, 0, 1, 0, k_ind) = val42;
    H_0(3, 1, 1, 1, k_ind) = val42;
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_EG_HPP

// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang
//
// Tetralayer lattice 

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TETRALAYER_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TETRALAYER_LATTICE_HPP

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

template <typename /*symmetry_group*/>
class tetralayer_lattice {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
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
double* tetralayer_lattice<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
double* tetralayer_lattice<point_group_type>::initialize_k_DCA_basis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
double* tetralayer_lattice<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
double* tetralayer_lattice<point_group_type>::initialize_k_LDA_basis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> tetralayer_lattice<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;
  flavors[2] = 2;
  flavors[3] = 3;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> tetralayer_lattice<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> tetralayer_lattice<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void tetralayer_lattice<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Tetralayer lattice has four bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U1 = parameters.get_U1();              // Same band, opposite spin.
  const double U2 = parameters.get_U2();
  const double U3 = parameters.get_U3();
  const double U4 = parameters.get_U4();
  const double V  = parameters.get_V();              // Different band, opposite spin.
  const double V_prime = parameters.get_V_prime();   // Different band, same spin.

  H_interaction = 0.;

  H_interaction(0, 0, 0, 1, origin) = U1;
  H_interaction(0, 1, 0, 0, origin) = U1;
  H_interaction(1, 0, 1, 1, origin) = U2;
  H_interaction(1, 1, 1, 0, origin) = U2;
  H_interaction(2, 0, 2, 1, origin) = U3;
  H_interaction(2, 1, 2, 0, origin) = U3;
  H_interaction(3, 0, 3, 1, origin) = U4;
  H_interaction(3, 1, 3, 0, origin) = U4;

  for (int b1 = 0; b1 < BANDS; b1++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int b2 = 0; b2 < BANDS; b2++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (b1 != b2 && s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = V;

          if (b1 != b2 && s1 == s2)
            H_interaction(b1, s1, b2, s2, origin) = V_prime;
        }
      }
    }
  }
}

template <typename point_group_type>
template <class domain>
void tetralayer_lattice<point_group_type>::initialize_H_symmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

/*  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;
*/
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void tetralayer_lattice<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Tetralayer lattice has four bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto e1 = parameters.get_e1();
  const auto e2 = parameters.get_e2();
  const auto e3 = parameters.get_e3();
  const auto e4 = parameters.get_e4();
  const auto t1 = parameters.get_t1();
  const auto t2 = parameters.get_t2();
  const auto t3 = parameters.get_t3();
  const auto t4 = parameters.get_t4();
  const auto t1_prime = parameters.get_t1_prime();
  const auto t2_prime = parameters.get_t2_prime();
  const auto t3_prime = parameters.get_t3_prime();
  const auto t4_prime = parameters.get_t4_prime();
  const auto t_perp12 = parameters.get_t_perp12();
  const auto t_perp_prime12 = parameters.get_t_perp_prime12();
  const auto t_perp23 = parameters.get_t_perp23();
  const auto t_perp_prime23 = parameters.get_t_perp_prime23();
  const auto t_perp34 = parameters.get_t_perp34();
  const auto t_perp_prime34 = parameters.get_t_perp_prime34();
  const auto coefp = parameters.get_coefp();
  const auto coefm = parameters.get_coefm();

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    const auto val1 =
        e1 -2. * t1 * (std::cos(k[0]) + std::cos(k[1])) - 4. * t1_prime * std::cos(k[0]) * std::cos(k[1]);
    const auto val2 =
        e2 -2. * t2 * (std::cos(k[0]) + std::cos(k[1])) - 4. * t2_prime * std::cos(k[0]) * std::cos(k[1]);
    const auto val3 =
        e3 -2. * t3 * (std::cos(k[0]) + std::cos(k[1])) - 4. * t3_prime * std::cos(k[0]) * std::cos(k[1]);
    const auto val4 =
        e4 -2. * t4 * (std::cos(k[0]) + std::cos(k[1])) - 4. * t4_prime * std::cos(k[0]) * std::cos(k[1]);


    const auto val12 = -t_perp12 -2. * t_perp_prime12 * coefp * (std::cos(k[0]) + std::cos(k[1])) -2. * t_perp_prime12 * coefm * (std::cos(k[0]) - std::cos(k[1]));
      
    const auto val23 = -t_perp23 -2. * t_perp_prime23 * coefp * (std::cos(k[0]) + std::cos(k[1])) -2. * t_perp_prime23 * coefm * (std::cos(k[0]) - std::cos(k[1]));
    const auto val34 = -t_perp34 -2. * t_perp_prime34 * coefp * (std::cos(k[0]) + std::cos(k[1])) -2. * t_perp_prime34 * coefm * (std::cos(k[0]) - std::cos(k[1]));

    H_0(0, 0, 0, 0, k_ind) = val1;
    H_0(0, 1, 0, 1, k_ind) = val1;
    H_0(1, 0, 1, 0, k_ind) = val2;
    H_0(1, 1, 1, 1, k_ind) = val2;
    H_0(2, 0, 2, 0, k_ind) = val3;
    H_0(2, 1, 2, 1, k_ind) = val3;
    H_0(3, 0, 3, 0, k_ind) = val4;
    H_0(3, 1, 3, 1, k_ind) = val4;

    H_0(0, 0, 1, 0, k_ind) = val12;
    H_0(0, 1, 1, 1, k_ind) = val12;
    H_0(1, 0, 0, 0, k_ind) = val12;
    H_0(1, 1, 0, 1, k_ind) = val12;

    H_0(1, 0, 2, 0, k_ind) = val23;
    H_0(1, 1, 2, 1, k_ind) = val23;
    H_0(2, 0, 1, 0, k_ind) = val23;
    H_0(2, 1, 1, 1, k_ind) = val23;

    H_0(2, 0, 3, 0, k_ind) = val34;
    H_0(2, 1, 3, 1, k_ind) = val34;
    H_0(3, 0, 2, 0, k_ind) = val34;
    H_0(3, 1, 2, 1, k_ind) = val34;
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_TETRALAYER_LATTICE_HPP

// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang (Aug. 29, 2025)
// Bilayer two-orbital model for La-327
//
// Bilayer lattice with totally 4 orbitals

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_LA3NI2O7_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_LA3NI2O7_HPP

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
class La3Ni2O7 {
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
double* La3Ni2O7<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
double* La3Ni2O7<point_group_type>::initialize_k_DCA_basis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
double* La3Ni2O7<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
double* La3Ni2O7<point_group_type>::initialize_k_LDA_basis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> La3Ni2O7<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> La3Ni2O7<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> La3Ni2O7<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void La3Ni2O7<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("La3Ni2O7 has four bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double Ux = parameters.get_Ux();             // Same band, opposite spin.
  const double Uz = parameters.get_Uz();
  const double V  = parameters.get_V();              // Different band, opposite spin.
  const double V_prime = parameters.get_V_prime();   // Different band, same spin.

  H_interaction = 0.;
    
  // 0,1,2,3 denote dx2,dz2,dx2,dz2 respectively (bilayer 2-orb)

  H_interaction(0, 0, 0, 1, origin) = Ux;
  H_interaction(0, 1, 0, 0, origin) = Ux;
  H_interaction(1, 0, 1, 1, origin) = Uz;
  H_interaction(1, 1, 1, 0, origin) = Uz;
  H_interaction(2, 0, 2, 1, origin) = Ux;
  H_interaction(2, 1, 2, 0, origin) = Ux;
  H_interaction(3, 0, 3, 1, origin) = Uz;
  H_interaction(3, 1, 3, 0, origin) = Uz;

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
void La3Ni2O7<point_group_type>::initialize_H_symmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

/*  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;
*/
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void La3Ni2O7<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Bilayer lattice has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  // 0,1,2,3 denote dx2,dz2,dx2,dz2 respectively (bilayer 2-orb)
    
  const auto ex = parameters.get_ex();
  const auto ez = parameters.get_ez();
    
  // below 1 and 2 denotes dx2 and dz2 orb
  // intralayer intra-orbital and inter-orbital 
  const auto t11x  = parameters.get_t11x();
  const auto t11xy = parameters.get_t11xy();
  const auto t11xx = parameters.get_t11xx();
    
  const auto t22x  = parameters.get_t22x();
  const auto t22xy = parameters.get_t22xy();
  const auto t22xx = parameters.get_t22xx();
    
  const auto t12x  = parameters.get_t12x();
  const auto t12xx = parameters.get_t12xx();
    
  // interlayer intra-orbital and inter-orbital 
  const auto s110  = parameters.get_s110();
  const auto s11x  = parameters.get_s11x();
  const auto s11xy = parameters.get_t11xy();
  const auto s11xx = parameters.get_t11xx();
    
  const auto s220  = parameters.get_s220();
  const auto s22x  = parameters.get_s22x();
  const auto s22xy = parameters.get_s22xy();
  const auto s22xx = parameters.get_s22xx();
    
  const auto s12x  = parameters.get_s12x();
  const auto s12xx = parameters.get_s12xx();
    

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
      
    // the follwing do not use conventional -t (which is used for an electron-like band)
    // here directly use DFT parameters (with electron language)
      
    // intralayer
    const auto Hx = ex + 2.* t11x * (std::cos(k[0]) + std::cos(k[1])) + 4. * t11xy * std::cos(k[0])*std::cos(k[1]) + 2. * t11xx * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
      
    const auto Hz = ez + 2.* t22x * (std::cos(k[0]) + std::cos(k[1])) + 4. * t22xy * std::cos(k[0])*std::cos(k[1]) + 2. * t22xx * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
      
    const auto V = 2.* t12x * (std::cos(k[0]) - std::cos(k[1])) + 2. * t12xx * (std::cos(2.*k[0]) - std::cos(2.*k[1]));
      
    // interlayer
    const auto Hxp = s110 + 2.* s11x * (std::cos(k[0]) + std::cos(k[1])) + 4. * s11xy * std::cos(k[0])*std::cos(k[1]) + 2. * s11xx * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
      
    const auto Hzp = s220 + 2.* s22x * (std::cos(k[0]) + std::cos(k[1])) + 4. * s22xy * std::cos(k[0])*std::cos(k[1]) + 2. * s22xx * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
      
    const auto Vp = 2.* s12x * (std::cos(k[0]) - std::cos(k[1])) + 2. * s12xx * (std::cos(2.*k[0]) - std::cos(2.*k[1]));
      
      
    for (int s = 0; s < 2; s++) {
        // intralayer
        H_0(0, s, 0, s, k_ind) = Hx;
        H_0(1, s, 1, s, k_ind) = Hz;
        H_0(0, s, 1, s, k_ind) = V;
        H_0(1, s, 0, s, k_ind) = V;
        
        H_0(2, s, 2, s, k_ind) = Hx;
        H_0(3, s, 3, s, k_ind) = Hz;
        H_0(2, s, 3, s, k_ind) = V;
        H_0(3, s, 2, s, k_ind) = V;
        
        // interlayer
        H_0(0, s, 2, s, k_ind) = Hxp;
        H_0(2, s, 0, s, k_ind) = Hxp;
        H_0(1, s, 3, s, k_ind) = Hzp;
        H_0(3, s, 1, s, k_ind) = Hzp;
        
        H_0(0, s, 3, s, k_ind) = Vp;
        H_0(3, s, 0, s, k_ind) = Vp;
        H_0(1, s, 2, s, k_ind) = Vp;
        H_0(2, s, 1, s, k_ind) = Vp;
    }
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_LA3NI2O7_HPP

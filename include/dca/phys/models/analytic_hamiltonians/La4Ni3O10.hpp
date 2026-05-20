// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang (Aug. 29, 2025)
// Trilayer two-orbital model for La-4310
//
// Trilayer lattice with totally 6 orbitals

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_LA4NI3O10_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_LA4NI3O10_HPP

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
class La4Ni3O10 {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
 //  typedef point_group_type DCA_point_group;

  // Aug.17, 2023 debug
  // for C4-C2 rotational sym breaking cases:
  typedef domains::no_symmetry<2> DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 6;

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
double* La4Ni3O10<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];

  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
double* La4Ni3O10<point_group_type>::initialize_k_DCA_basis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
double* La4Ni3O10<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
double* La4Ni3O10<point_group_type>::initialize_k_LDA_basis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> La4Ni3O10<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;
  flavors[2] = 2;
  flavors[3] = 3;
  flavors[4] = 4;
  flavors[5] = 5;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> La4Ni3O10<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs(BANDS, std::vector<double>(DIMENSION, 0.));

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> La4Ni3O10<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void La4Ni3O10<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("La4Ni3O10 has six bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double Ux = parameters.get_Ux();             // Same band, opposite spin.
  const double Uz = parameters.get_Uz();
  const double V  = parameters.get_V();              // Different band, opposite spin.
  const double V_prime = parameters.get_V_prime();   // Different band, same spin.

  H_interaction = 0.;
    
  // 0,1,2,3,4,5 denote dx2,dz2,dx2,dz2,dx2,dz2 respectively (trilayer 2-orb)

  H_interaction(0, 0, 0, 1, origin) = Ux;
  H_interaction(0, 1, 0, 0, origin) = Ux;
  H_interaction(1, 0, 1, 1, origin) = Uz;
  H_interaction(1, 1, 1, 0, origin) = Uz;
  H_interaction(2, 0, 2, 1, origin) = Ux;
  H_interaction(2, 1, 2, 0, origin) = Ux;
  H_interaction(3, 0, 3, 1, origin) = Uz;
  H_interaction(3, 1, 3, 0, origin) = Uz;
  H_interaction(4, 0, 4, 1, origin) = Ux;
  H_interaction(4, 1, 4, 0, origin) = Ux;
  H_interaction(5, 0, 5, 1, origin) = Uz;
  H_interaction(5, 1, 5, 0, origin) = Uz;

  // intralayer's inter-orbital Hund's J
  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      if (s1 != s2) {
        H_interaction(0, s1, 1, s2, origin) = V;
        H_interaction(1, s1, 0, s2, origin) = V;
        H_interaction(2, s1, 3, s2, origin) = V;
        H_interaction(3, s1, 2, s2, origin) = V;
        H_interaction(4, s1, 5, s2, origin) = V;
        H_interaction(5, s1, 4, s2, origin) = V;
      }
      if (s1 == s2) {
        H_interaction(0, s1, 1, s2, origin) = V_prime;
        H_interaction(1, s1, 0, s2, origin) = V_prime;
        H_interaction(2, s1, 3, s2, origin) = V_prime;
        H_interaction(3, s1, 2, s2, origin) = V_prime;
        H_interaction(4, s1, 5, s2, origin) = V_prime;
        H_interaction(5, s1, 4, s2, origin) = V_prime;
      }
    }
  }
}

template <typename point_group_type>
template <class domain>
void La4Ni3O10<point_group_type>::initialize_H_symmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

/*  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;
*/
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void La4Ni3O10<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("Trilayer lattice has six bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  // x,z denote dx2,dz2 respectively (trilayer 2-orb) in Hamiltonain
    
  const auto ex1 = parameters.get_ex1();
  const auto ez1 = parameters.get_ez1();
  const auto ex2 = parameters.get_ex2();
  const auto ez2 = parameters.get_ez2();
    
  // below 1, 2 denotes inner and outer layers
  // t, s denotes the intralaer, interlayer  hoppping, 
  // e.g. t1x1z2 denotes the nearest-neighbor hopping between inner dx2 and outer dz2

  // intralayer intra-orbital
  const auto t1x1x1 = parameters.get_t1x1x1();
  const auto t2x1x1 = parameters.get_t2x1x1();
  const auto t3x1x1 = parameters.get_t3x1x1();
  
  const auto t1x2x2 = parameters.get_t1x2x2();
  const auto t2x2x2 = parameters.get_t2x2x2();
  const auto t3x2x2 = parameters.get_t3x2x2();
 
  const auto t1z1z1 = parameters.get_t1z1z1();
  const auto t2z1z1 = parameters.get_t2z1z1();
  const auto t3z1z1 = parameters.get_t3z1z1();

  const auto t1z2z2 = parameters.get_t1z2z2();
  const auto t2z2z2 = parameters.get_t2z2z2();
  const auto t3z2z2 = parameters.get_t3z2z2();
    
  // intralayer inter-orbital
  const auto t1x1z1 = parameters.get_t1x1z1();
  const auto t2x1z1 = parameters.get_t2x1z1();

  const auto t1x2z2 = parameters.get_t1x2z2();
  const auto t2x2z2 = parameters.get_t2x2z2();
    
  // interlayer intra-orbital
  const auto s0x1x2 = parameters.get_s0x1x2();
  const auto s1x1x2 = parameters.get_s1x1x2();
  const auto s2x1x2 = parameters.get_s2x1x2();
  const auto s3x1x2 = parameters.get_s3x1x2();

  const auto s0x2x2 = parameters.get_s0x2x2();
  const auto s1x2x2 = parameters.get_s1x2x2();
  const auto s2x2x2 = parameters.get_s2x2x2();
  const auto s3x2x2 = parameters.get_s3x2x2();

  const auto s0z1z2 = parameters.get_s0z1z2();
  const auto s1z1z2 = parameters.get_s1z1z2();
  const auto s2z1z2 = parameters.get_s2z1z2();
  const auto s3z1z2 = parameters.get_s3z1z2();

  const auto s0z2z2 = parameters.get_s0z2z2();
  const auto s1z2z2 = parameters.get_s1z2z2();
  const auto s2z2z2 = parameters.get_s2z2z2();
  const auto s3z2z2 = parameters.get_s3z2z2();

  // interlayer inter-orbital
  const auto s1x1z2 = parameters.get_s1x1z2();
  const auto s2x1z2 = parameters.get_s2x1z2();
  
  const auto s1x2z2 = parameters.get_s1x2z2();
  const auto s2x2z2 = parameters.get_s2x2z2();
    

  H_0 = ScalarType(0);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
      
    // the follwing do not use conventional -t (which is used for an electron-like band)
    // here directly use DFT parameters (with electron language)
      
    // intralayer
    const auto Hx11 = ex1 + 2.* t1x1x1 * (std::cos(k[0]) + std::cos(k[1])) + 4. * t2x1x1 * std::cos(k[0])*std::cos(k[1]) + 2. * t3x1x1 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto Hz11 = ez1 + 2.* t1z1z1 * (std::cos(k[0]) + std::cos(k[1])) + 4. * t2z1z1 * std::cos(k[0])*std::cos(k[1]) + 2. * t3z1z1 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto Hx22 = ex2 + 2.* t1x2x2 * (std::cos(k[0]) + std::cos(k[1])) + 4. * t2x2x2 * std::cos(k[0])*std::cos(k[1]) + 2. * t3x2x2 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto Hz22 = ez2 + 2.* t1z2z2 * (std::cos(k[0]) + std::cos(k[1])) + 4. * t2z2z2 * std::cos(k[0])*std::cos(k[1]) + 2. * t3z2z2 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto V11 = 2.* t1x1z1 * (std::cos(k[0]) - std::cos(k[1])) + 2. * t2x1z1 * (std::cos(2.*k[0]) - std::cos(2.*k[1]));
    const auto V22 = 2.* t1x2z2 * (std::cos(k[0]) - std::cos(k[1])) + 2. * t2x2z2 * (std::cos(2.*k[0]) - std::cos(2.*k[1]));

    // interlayer
    const auto Hxp12 = s0x1x2 + 2.* s1x1x2 * (std::cos(k[0]) + std::cos(k[1])) + 4. * s2x1x2 * std::cos(k[0])*std::cos(k[1]) + 2. * s3x1x2 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto Hzp12 = s0z1z2 + 2.* s1z1z2 * (std::cos(k[0]) + std::cos(k[1])) + 4. * s2z1z2 * std::cos(k[0])*std::cos(k[1]) + 2. * s3z1z2 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto Hxp22 = s0x2x2 + 2.* s1x2x2 * (std::cos(k[0]) + std::cos(k[1])) + 4. * s2x2x2 * std::cos(k[0])*std::cos(k[1]) + 2. * s3x2x2 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto Hzp22 = s0z2z2 + 2.* s1z2z2 * (std::cos(k[0]) + std::cos(k[1])) + 4. * s2z2z2 * std::cos(k[0])*std::cos(k[1]) + 2. * s3z2z2 * (std::cos(2.*k[0]) + std::cos(2.*k[1]));
    const auto Vp12 = 2.* s1x1z2 * (std::cos(k[0]) - std::cos(k[1])) + 2. * s2x1z2 * (std::cos(2.*k[0]) - std::cos(2.*k[1]));
    const auto Vp22 = 2.* s1x2z2 * (std::cos(k[0]) - std::cos(k[1])) + 2. * s2x2z2 * (std::cos(2.*k[0]) - std::cos(2.*k[1]));
      
    // 0, 1 denotes x, z orbital in lower layer (outer)
    // 2, 3 denotes x, z orbital in middle layer (inner)
    // 4, 5 denotes x, z orbital in upper layer  (outer)

    for (int s = 0; s < 2; s++) {
        // intralayer
        H_0(0, s, 0, s, k_ind) = Hx22;
        H_0(1, s, 1, s, k_ind) = Hz22;
        H_0(0, s, 1, s, k_ind) = V22;
        H_0(1, s, 0, s, k_ind) = V22;
        
        H_0(2, s, 2, s, k_ind) = Hx11;
        H_0(3, s, 3, s, k_ind) = Hz11;
        H_0(2, s, 3, s, k_ind) = V11;
        H_0(3, s, 2, s, k_ind) = V11;
  
        H_0(4, s, 4, s, k_ind) = Hx22;
        H_0(5, s, 5, s, k_ind) = Hz22;
        H_0(4, s, 5, s, k_ind) = V22;
        H_0(5, s, 4, s, k_ind) = V22;

        // interlayer
        H_0(0, s, 2, s, k_ind) = Hxp12;
        H_0(2, s, 0, s, k_ind) = Hxp12;
        H_0(1, s, 3, s, k_ind) = Hzp12;
        H_0(3, s, 1, s, k_ind) = Hzp12;

        H_0(2, s, 4, s, k_ind) = Hxp12;
        H_0(4, s, 2, s, k_ind) = Hxp12;
        H_0(3, s, 5, s, k_ind) = Hzp12;
        H_0(5, s, 3, s, k_ind) = Hzp12;

        H_0(0, s, 4, s, k_ind) = Hxp22;
        H_0(4, s, 0, s, k_ind) = Hxp22;
        H_0(1, s, 5, s, k_ind) = Hzp22;
        H_0(5, s, 1, s, k_ind) = Hzp22;
        
        H_0(0, s, 3, s, k_ind) = Vp12;
        H_0(3, s, 0, s, k_ind) = Vp12;
        H_0(1, s, 2, s, k_ind) = Vp12;
        H_0(2, s, 1, s, k_ind) = Vp12;
        
        H_0(2, s, 5, s, k_ind) = Vp12;
        H_0(5, s, 2, s, k_ind) = Vp12;
        H_0(3, s, 4, s, k_ind) = Vp12;
        H_0(4, s, 3, s, k_ind) = Vp12;

        H_0(0, s, 5, s, k_ind) = Vp22;
        H_0(5, s, 0, s, k_ind) = Vp22;
        H_0(1, s, 4, s, k_ind) = Vp22;
        H_0(4, s, 1, s, k_ind) = Vp22;
    }
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_LA4NI3O10_HPP

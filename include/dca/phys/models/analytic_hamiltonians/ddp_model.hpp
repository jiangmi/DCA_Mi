// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Mi Jiang
//
// Three-band model added into dz2 orbital

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_DDP_MODEL_LATTICE_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_DDP_MODEL_LATTICE_HPP

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

// TODO: the symmetry of this system must be checked.
template <typename /*symmetry_group*/>
class ddp_model {
public:
  using LDA_point_group = domains::no_symmetry<2>;
  using DCA_point_group = domains::no_symmetry<2>;
  //typedef point_group_type DCA_point_group;

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
double* ddp_model<point_group_type>::initialize_r_DCA_basis() {
  static double* r_DCA = new double[4];
  r_DCA[0] = 1.0;
  r_DCA[1] = 0.0;
  r_DCA[2] = 0.0;
  r_DCA[3] = 1.0;

  return r_DCA;
}

template <typename point_group_type>
double* ddp_model<point_group_type>::initialize_k_DCA_basis() {
  static double* k_DCA = new double[4];

  k_DCA[0] = 2 * M_PI;
  k_DCA[1] = 0.;
  k_DCA[2] = 0.;
  k_DCA[3] = 2 * M_PI;

  return k_DCA;
}

template <typename point_group_type>
double* ddp_model<point_group_type>::initialize_r_LDA_basis() {
  static double* r_LDA = new double[4];

  r_LDA[0] = 1.;
  r_LDA[1] = 0.;
  r_LDA[2] = 0.;
  r_LDA[3] = 1.;

  return r_LDA;
}

template <typename point_group_type>
double* ddp_model<point_group_type>::initialize_k_LDA_basis() {
  static double* k_LDA = new double[4];

  k_LDA[0] = 2. * M_PI;
  k_LDA[1] = 0.;
  k_LDA[2] = 0.;
  k_LDA[3] = 2. * M_PI;

  return k_LDA;
}

template <typename point_group_type>
std::vector<int> ddp_model<point_group_type>::get_flavors() {
  static std::vector<int> flavors(BANDS);

  flavors[0] = 0;
  flavors[1] = 1;
  flavors[2] = 2;
  flavors[3] = 3;

  return flavors;
}

template <typename point_group_type>
std::vector<std::vector<double>> ddp_model<point_group_type>::get_a_vectors() {
  static std::vector<std::vector<double>> a_vecs{
      std::vector<double>{0, 0}, std::vector<double>{0.5, 0}, std::vector<double>{0, 0.5},
      std::vector<double>{0, 0}};

  return a_vecs;
}

template <typename point_group_type>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> ddp_model<
    point_group_type>::get_orbital_permutations() {
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> permutations(0);
  return permutations;
}

template <typename point_group_type>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void ddp_model<point_group_type>::initialize_H_interaction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("dp plus model has 4 bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U_dd  = parameters.get_U_dd();  // interaction in d band
  const double U_pp  = parameters.get_U_pp();  // interaction in p bands

  const double V  = parameters.get_V();             // Different band, opposite spin.
  const double V_prime = parameters.get_V_prime();  // Different band, same spin.

  H_interaction = 0.;

  for (int s1 = 0; s1 < 2; s1++) {
    for (int s2 = 0; s2 < 2; s2++) {
      if (s1 != s2) {
        H_interaction(0, s1, 0, s2, origin) = U_dd;
        H_interaction(1, s1, 1, s2, origin) = U_pp;
        H_interaction(2, s1, 2, s2, origin) = U_pp;
        H_interaction(3, s1, 3, s2, origin) = U_dd;
        H_interaction(0, s1, 3, s2, origin) = V;
        H_interaction(3, s1, 0, s2, origin) = V;
      }
      if (s1 == s2) {
        H_interaction(0, s1, 3, s2, origin) = V_prime;
        H_interaction(3, s1, 0, s2, origin) = V_prime;
      }
    }
  }
}

template <typename point_group_type>
template <class domain>
void ddp_model<point_group_type>::initialize_H_symmetry(
    func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  //  H_symmetry(i, s1, j, s2)
  //  H_symmetries(0, 0, 0, 0) = 0; // at b0, G of spin 0 or 1 has the same values.
  //  H_symmetries(0, 1, 0, 1) = 0;

  //  H_symmetries(1, 0, 1, 0) = 1;
  //  H_symmetries(1, 1, 1, 1) = 1; // at i, G of spin 0 or 1 has the same values.
}

template <typename point_group_type>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void ddp_model<point_group_type>::initialize_H_0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("dp plus model has 4 bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t_pd = parameters.get_t_pd();
  const auto t_pp = parameters.get_t_pp();
  const auto t_pdz2 = parameters.get_t_pdz2();

  const auto ep_d = parameters.get_ep_d();
  const auto ep_p = parameters.get_ep_p();
  const auto edz2 = parameters.get_edz2();

  H_0 = ScalarType(0);

  const ScalarType I(0, 1.);

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = k_vecs[k_ind];
    const auto valdpx =  2. * I * t_pd * std::sin(k[0] / 2.);
    const auto valdpy = -2. * I * t_pd * std::sin(k[1] / 2.);
    const auto valpxpy = 4. * t_pp * std::sin(k[0] / 2.) * std::sin(k[1] / 2.);
    const auto valdz2px =  2. * I * t_pdz2 * std::sin(k[0] / 2.);
    const auto valdz2py = -2. * I * t_pdz2 * std::sin(k[1] / 2.);

    for (int s = 0; s < 2; s++) {
      H_0(0, s, 0, s, k_ind) = ep_d;
      H_0(1, s, 1, s, k_ind) = ep_p;
      H_0(2, s, 2, s, k_ind) = ep_p;
      H_0(3, s, 3, s, k_ind) = edz2;

      H_0(0, s, 1, s, k_ind) =  valdpx;
      H_0(1, s, 0, s, k_ind) = -valdpx;

      H_0(0, s, 2, s, k_ind) =  valdpy;
      H_0(2, s, 0, s, k_ind) = -valdpy;

      H_0(1, s, 2, s, k_ind) = valpxpy;
      H_0(2, s, 1, s, k_ind) = valpxpy;

      // d3z2r2 has +,-,+ sign structure so that it is - in x-y plane
      // see hamiltonian.py of multiplet python code
      H_0(3, s, 1, s, k_ind) = -valdz2px;
      H_0(1, s, 3, s, k_ind) =  valdz2px;

      H_0(3, s, 2, s, k_ind) =  valdz2py;
      H_0(2, s, 3, s, k_ind) = -valdz2py;
    }
  }
}

}  // models
}  // phys
}  // dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_DDP_MODEL_LATTICE_HPP

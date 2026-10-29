// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// This file defines the schemes used to mix the cluster self-energy in the DCA self-consistency
// loop.

#ifndef DCA_PHYS_PARAMETERS_SELF_ENERGY_MIXING_SCHEME_HPP
#define DCA_PHYS_PARAMETERS_SELF_ENERGY_MIXING_SCHEME_HPP

#include <string>

namespace dca {
namespace phys {
// dca::phys::

// LINEAR:   Sigma^(n+1) = alpha * Sigma_QMC + (1 - alpha) * Sigma_cg
// ANDERSON: Anderson/Pulay (DIIS) acceleration of the same fixed-point iteration.
enum class SelfEnergyMixingScheme { LINEAR, ANDERSON };

SelfEnergyMixingScheme stringToSelfEnergyMixingScheme(const std::string& str);

std::string toString(SelfEnergyMixingScheme scheme);

}  // phys
}  // dca

#endif  // DCA_PHYS_PARAMETERS_SELF_ENERGY_MIXING_SCHEME_HPP

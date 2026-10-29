// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// This file implements the conversion between SelfEnergyMixingScheme and string.

#include "dca/phys/parameters/self_energy_mixing_scheme.hpp"

#include <stdexcept>

namespace dca {
namespace phys {
// dca::phys::

SelfEnergyMixingScheme stringToSelfEnergyMixingScheme(const std::string& str) {
  if (str == "LINEAR")
    return SelfEnergyMixingScheme::LINEAR;
  else if (str == "ANDERSON")
    return SelfEnergyMixingScheme::ANDERSON;
  else
    throw(std::logic_error("Invalid self-energy mixing scheme."));
}

std::string toString(const SelfEnergyMixingScheme scheme) {
  switch (scheme) {
    case SelfEnergyMixingScheme::LINEAR:
      return "LINEAR";
    case SelfEnergyMixingScheme::ANDERSON:
      return "ANDERSON";
    default:
      throw(std::logic_error("Invalid self-energy mixing scheme."));
  }
}

}  // phys
}  // dca

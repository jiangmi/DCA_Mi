// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// Tests the 2D NDFT performed by the cached_ndft class.

#include "dca/phys/dca_step/cluster_solver/shared_tools/accumulation/tp/ndft/cached_ndft_cpu.hpp"

#include <complex>
#include <test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/single_sector_accumulation_test.hpp>

#include "gtest/gtest.h"

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/function/util/difference.hpp"
#include "dca/profiling/events/time.hpp"
#include "test/unit/phys/dca_step/cluster_solver/shared_tools/accumulation/single_sector_accumulation_test.hpp"

constexpr int n_sites = 4;
constexpr int n_bands = 3;
constexpr int n_frqs = 16;

template <class Real>
using CachedNdftCpuTest = dca::testing::SingleSectorAccumulationTest<Real, n_bands, n_sites, n_frqs>;

template <typename Real>
double computeWithFastDNFT(const typename CachedNdftCpuTest<Real>::Configuration& config,
                           const typename CachedNdftCpuTest<Real>::Matrix& M,
                           typename CachedNdftCpuTest<Real>::F_w_w& f_w);

using TestTypes = ::testing::Types<float, double>;
TYPED_TEST_CASE(CachedNdftCpuTest, TestTypes);

// Compare the result provided by the CPU version of CachedNdft::execute with the definition of the
// DNFT f(w1, w2) = \sum_{t1, t2} f(t1, t2) exp(i * t1 * w1 - t2 w2) stored in f_baseline_.
TYPED_TEST(CachedNdftCpuTest, Execute) {
  constexpr int n_samples = 40;

  TestFixture::prepareConfiguration(TestFixture::configuration_, TestFixture::M_, n_samples);

  using Real = TypeParam;
  typename TestFixture::F_w_w f_w_fast("f_w_fast");
  const double time =
      computeWithFastDNFT<Real>(TestFixture::configuration_, TestFixture::M_, f_w_fast);

  auto f_baseline = TestFixture::compute2DFTBaseline();
  const auto err = dca::func::util::difference(f_baseline, f_w_fast);
  EXPECT_LT(err.l_inf, 100 * std::numeric_limits<Real>::epsilon());

  std::cout << "\nCached ndft time [sec]:\t " << time << "\n";
}

template <typename Real>
double computeWithFastDNFT(const typename CachedNdftCpuTest<Real>::Configuration& config,
                           const typename CachedNdftCpuTest<Real>::Matrix& M,
                           typename CachedNdftCpuTest<Real>::F_w_w& f_w) {
  using BDmn = typename CachedNdftCpuTest<Real>::BDmn;
  using RDmn = typename CachedNdftCpuTest<Real>::RDmn;
  using PosFreqDmn = typename CachedNdftCpuTest<Real>::PosFreqDmn;
  using FreqDmn = typename CachedNdftCpuTest<Real>::FreqDmn;
  dca::func::function<std::complex<double>,
                      dca::func::dmn_variadic<BDmn, BDmn, RDmn, RDmn, PosFreqDmn, FreqDmn>>
      f_b_b_r_r_w_w;
  dca::phys::solver::accumulator::CachedNdft<double, RDmn, FreqDmn, PosFreqDmn, dca::linalg::CPU> nft_obj;

  dca::profiling::WallTime start_time;
  nft_obj.execute(config, M, f_b_b_r_r_w_w);
  dca::profiling::WallTime end_time;

  // Rearrange output.
  const int n_w = PosFreqDmn::dmn_size();
  auto invert_w = [=](const int w) { return 2 * n_w - 1 - w; };
  for (int b2 = 0; b2 < BDmn::dmn_size(); ++b2)
    for (int b1 = 0; b1 < BDmn::dmn_size(); ++b1)
      for (int r2 = 0; r2 < RDmn::dmn_size(); ++r2)
        for (int r1 = 0; r1 < RDmn::dmn_size(); ++r1)
          for (int w2 = 0; w2 < FreqDmn::dmn_size(); ++w2)
            for (int w1 = 0; w1 < n_w; ++w1) {
              f_w(b1, b2, r1, r2, w1 + n_w, w2) = f_b_b_r_r_w_w(b1, b2, r1, r2, w1, w2);
              f_w(b1, b2, r1, r2, invert_w(w1 + n_w), invert_w(w2)) =
                  std::conj(f_b_b_r_r_w_w(b1, b2, r1, r2, w1, w2));
            }

  dca::profiling::Duration duration(end_time, start_time);
  return duration.sec + 1.e-6 * duration.usec;
}

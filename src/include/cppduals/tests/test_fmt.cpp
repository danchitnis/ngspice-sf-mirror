//===-- test_dual.cpp - test duals/dual -------------------------*- C++ -*-===//
//
// Part of the cppduals project.
// https://gitlab.com/tesch1/cppduals
//
// See https://gitlab.com/tesch1/cppduals/blob/master/LICENSE.txt for
// license information.
//
// (c)2019 Michael Tesch. tesch1@gmail.com
//
/**
 * \file test_dual Dual number libfmt formatting tests
 *
 * (c)2019 Michael Tesch. tesch1@gmail.com
 */
#include <fmt/format.h>
//#include <fmt/std.h>
#define CPPDUALS_LIBFMT
#define CPPDUALS_LIBFMT_COMPLEX
#include <duals/dual>
#include "gtest/gtest.h"

using duals::dualf;
using duals::duald;
using duals::dualld;
typedef std::complex<double> complexd;
typedef std::complex<float> complexf;
typedef std::complex<duald> cduald;
typedef std::complex<dualf> cdualf;
using namespace duals::literals;
using namespace std::complex_literals;

TEST(libfmt, float_) {
  std::string s = fmt::format("{:.1f}", 2.f);
  EXPECT_EQ(s, "2.0");
}
TEST(libfmt, double_) {
  std::string s = fmt::format("{:.1f}", 2.);
  EXPECT_EQ(s, "2.0");
}
TEST(libfmt, complex_) {
  std::string s = fmt::format("{:.1f}", 2.f + 3if);
  EXPECT_EQ(s, "(2.0+3.0if)");
  s = fmt::format("{:.2f}", 2.f + 3if);
  EXPECT_EQ(s, "(2.00+3.00if)");
}
TEST(libfmt, complex_el) {
  std::string s = fmt::format("{}", 2.l + 3il);
  EXPECT_EQ(s, "(2+3il)");
}
TEST(libfmt, complex_flags) {
  std::string s;

  s = fmt::format("{}", 2. + 3i);
  EXPECT_EQ(s, "(2+3i)");
}
TEST(libfmt, complex_flags_dollar) {
  std::string s;
  s = fmt::format("{:$}", 2. + 3i);
  EXPECT_EQ(s, "(2+3i)");
}
TEST(libfmt, complex_flags_star) {
  std::string s;
  s = fmt::format("{:*}", 2. + 3i);
  EXPECT_EQ(s, "(2+3*i)");
}
TEST(libfmt, complex_flags_star_f) {
  std::string s;
  s = fmt::format("{:*.1f}", 2. + 3i);
  EXPECT_EQ(s, "(2.0+3.0*i)");
}
TEST(libfmt, complex_flags_comma) {
  std::string s;
  s = fmt::format("{:,.1f}", 2. + 3i);
  EXPECT_EQ(s, "(2.0,3.0)");

  s = fmt::format("{:,}", 2. + 3i);
  EXPECT_EQ(s, "(2,3)");
}

TEST(libfmt, complex_flags_a) {
  std::string s;

  // + +
  s = fmt::format("{:$+}", 2. + 3i);
  EXPECT_EQ(s, "(+2+3i)");

  s = fmt::format("{:$+.1f}", 2. + 3i);
  EXPECT_EQ(s, "(+2.0+3.0i)");

  s = fmt::format("{:*+}", 2. + 3i);
  EXPECT_EQ(s, "(+2+3*i)");

  s = fmt::format("{:,+}", 2. + 3i);
  EXPECT_EQ(s, "(+2,+3)");
}

TEST(libfmt, complex_flags_b) {
  std::string s;

  // + -
  s = fmt::format("{:$+}", 2. - 3i);
  EXPECT_EQ(s, "(+2-3i)");

  s = fmt::format("{:*+}", 2. - 3i);
  EXPECT_EQ(s, "(+2-3*i)");

  s = fmt::format("{:,+}", 2. - 3i);
  EXPECT_EQ(s, "(+2,-3)");

  s = fmt::format("{:,+.1f}", 2. - 3i);
  EXPECT_EQ(s, "(+2.0,-3.0)");

}

TEST(libfmt, complex_all_real)  {
  std::string s;
  s = fmt::format("{}", 2. + 0i);
  EXPECT_EQ(s, "(2)");

  s = fmt::format("{:.1f}", 2. + 0i);
  EXPECT_EQ(s, "(2.0)");

  s = fmt::format("{:*}", 2. + 0i);
  EXPECT_EQ(s, "(2)");

  s = fmt::format("{:,}", 2. + 0i);
  EXPECT_EQ(s, "(2,0)");
  
  s = fmt::format("{:,.1f}", 2. + 0i);
  EXPECT_EQ(s, "(2.0,0.0)");
}

TEST(libfmt, complex_all_imag)  {
  std::string s;

  s = fmt::format("{}", 0. + 3i);
  EXPECT_EQ(s, "(3i)");

  s = fmt::format("{:.1f}", 0. + 3i);
  EXPECT_EQ(s, "(3.0i)");

  s = fmt::format("{:*}", 0. + 3i);
  EXPECT_EQ(s, "(3*i)");

  s = fmt::format("{:,}", 0. + 3i);
  EXPECT_EQ(s, "(0,3)");
}
TEST(libfmt, complex_plus)  {
  std::string s = fmt::format("{:+}", 1. + 3i);
  EXPECT_EQ(s, "(+1+3i)");

  s = fmt::format("{:+.1f}", 1. + 3i);
  EXPECT_EQ(s, "(+1.0+3.0i)");

  s = fmt::format("{:+}", 1. - 3i);
  EXPECT_EQ(s, "(+1-3i)");
}
TEST(libfmt, complex_g) {
  std::string s = fmt::format("{:g}", 2.f + 3if);
  EXPECT_EQ(s, "(2+3if)");
}
TEST(libfmt, complex_gs) {
  std::string s = fmt::format("{:*g}", 3i);
  EXPECT_EQ(s, "(3*i)");
}
TEST(libfmt, complex_gel) {
  std::string s = fmt::format("{:g}", 2.l + 3il);
  EXPECT_EQ(s, "(2+3il)");
}

TEST(libfmt, dual_) {
  std::string s = fmt::format("{}", 2 + 3_ef);
  EXPECT_EQ(s, "(2+3_ef)");

  s = fmt::format("{:.1f}", 2 + 3_ef);
  EXPECT_EQ(s, "(2.0+3.0_ef)");
}
TEST(libfmt, dual_el) {
  std::string s = fmt::format("{}", 2 + 3_el);
  EXPECT_EQ(s, "(2+3_el)");

  s = fmt::format("{:.1f}", 2 + 3_el);
  EXPECT_EQ(s, "(2.0+3.0_el)");
}
TEST(libfmt, dual_flags) {
  std::string s;

  s = fmt::format("{}", 2. + 3_e);
  EXPECT_EQ(s, "(2+3_e)");

  s = fmt::format("{:$}", 2. + 3_e);
  EXPECT_EQ(s, "(2+3_e)");

  s = fmt::format("{:$.1f}", 2. + 3_e);
  EXPECT_EQ(s, "(2.0+3.0_e)");

  s = fmt::format("{:*}", 2. + 3_e);
  EXPECT_EQ(s, "(2+3*e)");

  s = fmt::format("{:,}", 2. + 3_e);
  EXPECT_EQ(s, "(2,3)");

  s = fmt::format("{:,.1f}", 2. + 3_e);
  EXPECT_EQ(s, "(2.0,3.0)");

  // + +
  s = fmt::format("{:$+}", 2. + 3_e);
  EXPECT_EQ(s, "(+2+3_e)");

  s = fmt::format("{:*+}", 2. + 3_e);
  EXPECT_EQ(s, "(+2+3*e)");

  s = fmt::format("{:*+.1f}", 2. + 3_e);
  EXPECT_EQ(s, "(+2.0+3.0*e)");

  s = fmt::format("{:,+}", 2. + 3_e);
  EXPECT_EQ(s, "(+2,+3)");

  // + -
  s = fmt::format("{:$+}", 2. - 3_e);
  EXPECT_EQ(s, "(+2-3_e)");

  s = fmt::format("{:*+}", 2. - 3_e);
  EXPECT_EQ(s, "(+2-3*e)");

  s = fmt::format("{:*+.1f}", 2. - 3_e);
  EXPECT_EQ(s, "(+2.0-3.0*e)");

  s = fmt::format("{:,+}", 2. - 3_e);
  EXPECT_EQ(s, "(+2,-3)");
  s = fmt::format("{:,+}", 2. + 3_e);
  EXPECT_EQ(s, "(+2,+3)");
}
TEST(libfmt, dual_all_real)  {
  std::string s;

  s = fmt::format("{}", 2 + 0_e);
  EXPECT_EQ(s, "(2)");
  s = fmt::format("{:*}", 2 + 0_e);
  EXPECT_EQ(s, "(2)");
  s = fmt::format("{:,}", 2 + 0_e);
  EXPECT_EQ(s, "(2,0)");
}
TEST(libfmt, dual_all_dual)  {
  std::string s;

  s = fmt::format("a{}b", 0 + 3_e);
  EXPECT_EQ(s, "a(3_e)b");
  s = fmt::format("a{:.1f}b", 0 + 3_e);
  EXPECT_EQ(s, "a(3.0_e)b");
  s = fmt::format("a{:*}b", 0 + 3_e);
  EXPECT_EQ(s, "a(3*e)b");
  s = fmt::format("a{:,}b", 0 + 3_e);
  EXPECT_EQ(s, "a(0,3)b");
  s = fmt::format("a{:,.1f}b", 0 + 3_e);
  EXPECT_EQ(s, "a(0.0,3.0)b");
}

TEST(libfmt, dual_plus)  {
  std::string s = fmt::format("{:+}", 1. + 3_e);
  EXPECT_EQ(s, "(+1+3_e)");

  s = fmt::format("{:+.1f}", 1. + 3_e);
  EXPECT_EQ(s, "(+1.0+3.0_e)");

  s = fmt::format("{:+}", 1. - 3_e);
  EXPECT_EQ(s, "(+1-3_e)");
}
TEST(libfmt, dual_g) {
  std::string s = fmt::format("{:g}", 2 + 3_ef);
  EXPECT_EQ(s, "(2+3_ef)");
}
TEST(libfmt, dual_gs) {
  std::string s = fmt::format("a{:*g}b", 3_e);
  EXPECT_EQ(s, "a(3*e)b");

  s = fmt::format("{:*+g}", 3_e);
  EXPECT_EQ(s, "(+3*e)");
}
TEST(libfmt, dual_gel) {
  std::string s = fmt::format("{:g}", 2 + 3_el);
  EXPECT_EQ(s, "(2+3_el)");
}

TEST(libfmt, dual_cd) {
  std::string s = fmt::format("{}", cdualf(2 + 3_ef, 4 + 5_ef));
  EXPECT_EQ(s, "((2+3_ef)+(4+5_ef)i)");
}

TEST(libfmt, dual_cgt) {
  std::string s = fmt::format("{:g}", cdualf(2 + 3_ef, 4 + 5_ef));
  EXPECT_EQ(s, "((2+3_ef)+(4+5_ef)i)");
}

TEST(libfmt, dual_cgts) {
  std::string s;

  s = fmt::format("{:**g}", cdualf(2 + 3_ef, 4 + 5_ef));
  EXPECT_EQ(s, "((2+3*ef)+(4+5*ef)*i)"); // todo - should be *if

  s = fmt::format("{:,*g}", cdualf(2 + 3_ef, 4 + 5_ef));
  EXPECT_EQ(s, "((2+3*ef),(4+5*ef))");

  // nonsense but should still work
  s = fmt::format("{:*,g}", cdualf(2 + 3_ef, 4 + 5_ef));
  EXPECT_EQ(s, "((2,3)+(4,5)*i)");

  s = fmt::format("{:,,g}", cdualf(2 + 3_ef, 4 + 5_ef));
  EXPECT_EQ(s, "((2,3),(4,5))");
}

// Copyright (C) 2021-2026 The Simons Foundation, Inc. and the CSBQ authors (see NOTICE)
// Author: Dhairya Malhotra
// SPDX-License-Identifier: Apache-2.0

// Helpers shared by the unit tests in test/unit/.

#ifndef _CSBQ_TEST_UNIT_COMMON_HPP_
#define _CSBQ_TEST_UNIT_COMMON_HPP_

#include <csbq.hpp>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace sctl {
namespace unit_test {

  /**
   * Number of failed checks in this program.
   */
  inline Long& FailCount() {
    static Long fail_cnt = 0;
    return fail_cnt;
  }

  /**
   * Print the result of a check, and count it if it fails. The check passes
   * when err <= tol; a NaN error fails.
   */
  inline bool Check(const std::string& name, const double err, const double tol) {
    const bool pass = (err <= tol);
    if (!pass) FailCount()++;
    std::cout<<(pass ? "  pass  " : "  FAIL  ")<<std::left<<std::setw(60)<<name<<std::right;
    std::cout<<std::scientific<<std::setprecision(2)<<"  err = "<<std::setw(9)<<err<<"  (tol = "<<tol<<")\n";
    return pass;
  }

  /**
   * Print a summary of the checks, and return the exit status of the program
   * (0 when every check passed).
   */
  inline int Summary(const std::string& test_name) {
    if (FailCount()) std::cout<<test_name<<": "<<FailCount()<<" check(s) FAILED\n";
    else std::cout<<test_name<<": all checks passed\n";
    return (FailCount() ? 1 : 0);
  }

  /**
   * Max-norm of U - Uref divided by scal, or by the max-norm of Uref when scal
   * is zero. Returns NaN if any entry of U - Uref is NaN.
   */
  template <class Real> Real MaxErr(const Vector<Real>& U, const Vector<Real>& Uref, Real scal = 0) {
    SCTL_ASSERT(U.Dim() == Uref.Dim());
    Real max_err = 0, max_val = 0;
    bool nan = false;
    for (Long i = 0; i < U.Dim(); i++) {
      const Real err = fabs(U[i] - Uref[i]);
      if (err != err) nan = true;
      max_err = std::max<Real>(max_err, err);
      max_val = std::max<Real>(max_val, fabs(Uref[i]));
    }
    if (nan) return std::numeric_limits<Real>::quiet_NaN();
    return max_err / (scal > 0 ? scal : max_val);
  }

  /**
   * Max-norm of a vector.
   */
  template <class Real> Real MaxNorm(const Vector<Real>& U) {
    Real max_val = 0;
    for (const auto& u : U) max_val = std::max<Real>(max_val, fabs(u));
    return max_val;
  }

  /**
   * Exit with an error if a precomputed quadrature table for the kernel, the
   * element order or the Fourier order is missing from CSBQ_DATA_PATH: CSBQ
   * would build it on the fly in quadruple precision, which is very slow.
   */
  template <class Kernel> void RequireQuadTables(const Long FourierOrder, const Long ElemOrder = 10) {
    const std::string path = CSBQ_QUOTEME(CSBQ_DATA_PATH);
    const std::string fname_lst[3] = {
      path + "/log_quad",
      path + "/special_quad_q" + std::to_string(ElemOrder) + "_" + Kernel::Name(),
      path + "/toroidal_quad_rule_m" + std::to_string(FourierOrder/2+2) + "_" + Kernel::Name() // FourierModes + ModalUpsample = (FourierOrder/2+1) + 1
    };
    for (const auto& fname : fname_lst) {
      if (!std::ifstream(fname).good()) {
        std::cout<<"Missing quadrature table: "<<fname<<'\n';
        exit(1);
      }
    }
  }

  /**
   * Breakpoints 0, 1/N, ..., 1 of N equal panels.
   */
  template <class Real> Vector<Real> UniformPanels(const Long N) {
    Vector<Real> s_brk(N+1);
    for (Long i = 0; i <= N; i++) s_brk[i] = i/(Real)N;
    return s_brk;
  }

  /**
   * Build an element list from a parametrization geom_fn(x, y, z, r, s), s in
   * [0,1], of the centerline and the cross-sectional radius, with one element
   * of order ElemOrder and Fourier order FourierOrder on each panel
   * [s_brk[i], s_brk[i+1]].
   *
   * @param[out] s_nds (optional) parameter s at the centerline nodes.
   *
   * @param[out] Xc (optional) centerline coordinates at the nodes (AoS order).
   *
   * @param[out] radius (optional) cross-sectional radius at the nodes.
   */
  template <class Real, class GeomFn> void BuildElemList(SlenderElemList<Real>& elem_lst, const GeomFn& geom_fn, const Vector<Real>& s_brk, const Long ElemOrder, const Long FourierOrder, Vector<Real>* s_nds = nullptr, Vector<Real>* Xc = nullptr, Vector<Real>* radius = nullptr) {
    const Long Nelem = s_brk.Dim()-1;
    const auto& nds = SlenderElemList<Real>::CenterlineNodes(ElemOrder);

    Vector<Real> s_, X_, r_;
    for (Long i = 0; i < Nelem; i++) {
      for (Long j = 0; j < ElemOrder; j++) {
        const Real s = s_brk[i] + (s_brk[i+1]-s_brk[i]) * nds[j];
        Real x, y, z, r;
        geom_fn(x, y, z, r, s);
        s_.PushBack(s);
        X_.PushBack(x);
        X_.PushBack(y);
        X_.PushBack(z);
        r_.PushBack(r);
      }
    }

    Vector<Long> elem_order(Nelem), fourier_order(Nelem);
    elem_order = ElemOrder;
    fourier_order = FourierOrder;
    elem_lst.Init(elem_order, fourier_order, X_, r_);

    if (s_nds) (*s_nds) = s_;
    if (Xc) (*Xc) = X_;
    if (radius) (*radius) = r_;
  }

  /**
   * Points off the surface of an element, along the normal at the surface
   * points with parameters s (in [0,1]) and theta (the angle from the
   * orientation vector of the element), at distance dist*r from the surface,
   * where r is the cross-sectional radius at s (outside for dist > 0, inside
   * for dist < 0). The points are appended to Xt in AoS order, looping over s,
   * then theta, then dist.
   */
  template <class Real> void OffSurfacePoints(Vector<Real>& Xt, const SlenderElemList<Real>& elem_lst, const Long elem_idx, const Vector<Real>& s, const Vector<Real>& theta, const Vector<Real>& dist) {
    constexpr Integer COORD_DIM = 3;
    const Long Ns = s.Dim(), Nt = theta.Dim();
    Vector<Real> sin_theta(Nt), cos_theta(Nt);
    for (Long i = 0; i < Nt; i++) {
      sin_theta[i] = sin<Real>(theta[i]);
      cos_theta[i] = cos<Real>(theta[i]);
    }

    Vector<Real> X, Xn, Xc;
    elem_lst.GetGeom(&X, &Xn, nullptr, nullptr, nullptr, s, sin_theta, cos_theta, elem_idx);
    elem_lst.GetGeom(&Xc, nullptr, nullptr, nullptr, nullptr, s, Vector<Real>{(Real)0}, Vector<Real>{(Real)0}, elem_idx); // with sin(theta) = cos(theta) = 0, the point is on the centerline
    for (Long i = 0; i < Ns; i++) {
      for (Long j = 0; j < Nt; j++) {
        const Long idx = i*Nt+j;
        Real r2 = 0;
        for (Integer k = 0; k < COORD_DIM; k++) {
          const Real dx = X[idx*COORD_DIM+k] - Xc[i*COORD_DIM+k];
          r2 += dx*dx;
        }
        const Real r = sqrt<Real>(r2);
        for (const auto d : dist) {
          for (Integer k = 0; k < COORD_DIM; k++) {
            Xt.PushBack(X[idx*COORD_DIM+k] + Xn[idx*COORD_DIM+k] * (d*r));
          }
        }
      }
    }
  }

}
}

#endif // _CSBQ_TEST_UNIT_COMMON_HPP_

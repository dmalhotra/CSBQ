// Copyright (C) 2021-2026 The Simons Foundation, Inc. and the CSBQ authors (see NOTICE)
// Author: Dhairya Malhotra
// SPDX-License-Identifier: Apache-2.0

/**
 * Unit tests for the layer potentials of a closed loop (non-planar centerline,
 * varying cross-sectional radius) with the Laplace and the Stokes kernels,
 * against exact values:
 *
 * 1) the double-layer potential D[c] of a constant density c, which is -c
 * inside the body, -c/2 on the surface and 0 outside.
 *
 * 2) Green's identity S[du/dn] - D[u] for the field u of a point source (with
 * du/dn the traction for Stokes). With the source outside the body it is u
 * inside, u/2 on the surface and 0 outside; with the source inside the body it
 * is 0 inside, -u/2 on the surface and -u outside.
 *
 * 3) the accuracy of the self-interaction (targets on the surface) and of the
 * near-interaction (targets at 1e-6 to 4 radii from the surface, on both
 * sides), and how it follows the accuracy tolerance.
 *
 * The errors must be within ErrFactor times the accuracy tolerance.
 *
 * The program returns a non-zero exit status if a check fails. To compile and
 * run all the unit tests, start in the project root directory and run:
 * make check
 */

#include "common.hpp"

#include <sstream>

using namespace sctl;
using namespace sctl::unit_test;

constexpr Integer COORD_DIM = 3;
constexpr Long ElemOrder = 10;
constexpr Long FourierOrder = 16; // quadrature tables exist for element order 10 and Fourier orders with FourierOrder/2+2 in {4,6,8,9,10,12,...}
constexpr double ErrFactor = 20; // error threshold relative to the accuracy tolerance

/**
 * Centerline and cross-sectional radius of the loop, for s in [0,1].
 */
void LoopGeom(double& x, double& y, double& z, double& r, const double s) {
  const double t = 2*const_pi<double>()*s;
  x = cos(t);
  y = sin(t);
  z = 0.1*sin(2*t);
  r = 0.05*(1 + 0.25*cos(t+1));
}

/**
 * Panel breakpoints in [0,1] for the loop, refined geometrically (ratio 1.5)
 * toward s = 0, where the point source inside the body sits: from panels of
 * length 0.4 radius there to panels of length 1/16.
 */
Vector<double> PanelBreakpoints() {
  double x, y, z, r0;
  LoopGeom(x, y, z, r0, 0);
  const double h0 = 0.4*r0/(2*const_pi<double>()), h_max = 1/16.0; // |dx/ds| is about 2 pi
  Vector<double> h; // panel lengths in [0,1/2]
  double h_sum = 0;
  for (double h_ = h0; h_sum < 0.5; h_ *= 1.5) {
    h.PushBack(std::min(h_, h_max));
    h_sum += std::min(h_, h_max);
  }

  const Long N = h.Dim();
  Vector<double> s_brk(2*N+1);
  s_brk[0] = 0;
  for (Long i = 0; i < N; i++) s_brk[i+1] = s_brk[i] + h[i] * (0.5/h_sum);
  for (Long i = 0; i < N; i++) s_brk[2*N-i] = 1 - s_brk[i]; // mirror about s = 1/2
  return s_brk;
}

/**
 * Target points off the surface, in groups (for reporting errors), each point
 * marked inside or outside the body.
 */
struct Targets {
  Vector<double> X; // coordinates (AoS order)
  Vector<Long> inside; // 1 for a point inside the body
  Vector<std::string> grp_name;
  Vector<Long> grp_dsp{0}; // first point of each group, and the number of points at the end

  void AddGroup(const std::string& name) {
    grp_name.PushBack(name);
    grp_dsp.PushBack(grp_dsp[grp_dsp.Dim()-1]);
  }
  void Add(const Vector<double>& X_, const bool inside_) {
    for (const auto& x : X_) X.PushBack(x);
    for (Long i = 0; i < X_.Dim()/COORD_DIM; i++) inside.PushBack(inside_);
    grp_dsp[grp_dsp.Dim()-1] = inside.Dim();
  }
};

Targets BuildTargets(const SlenderElemList<double>& elem_lst) {
  Targets trg;
  const Long Nelem = elem_lst.Size();
  const Vector<Long> elem_idx{0, Nelem/3, (2*Nelem)/3, Nelem-1}; // elements 0 and Nelem-1 are next to the interior source
  const Vector<double> s{0.3, 0.8}, theta{0.5, 2.6, 4.4};
  const auto add_near = [&](const double dist, const bool both_sides) {
    for (const auto e : elem_idx) {
      Vector<double> X;
      OffSurfacePoints(X, elem_lst, e, s, theta, Vector<double>{dist});
      trg.Add(X, false);
      if (both_sides) {
        X.ReInit(0);
        OffSurfacePoints(X, elem_lst, e, s, theta, Vector<double>{-dist});
        trg.Add(X, true);
      }
    }
  };

  trg.AddGroup("near, d = 1e-1 r");
  add_near(1e-1, true);
  trg.AddGroup("near, d = 1e-3 r");
  add_near(1e-3, true);
  trg.AddGroup("near, d = 1e-6 r");
  add_near(1e-6, true);
  trg.AddGroup("d = r and 4 r, outside");
  add_near(1, false);
  add_near(4, false);

  trg.AddGroup("centerline, inside");
  for (const double s0 : {0.15, 0.4, 0.65, 0.9}) {
    Vector<double> X(COORD_DIM);
    double r;
    LoopGeom(X[0], X[1], X[2], r, s0);
    trg.Add(X, true);
  }

  trg.AddGroup("far, outside");
  trg.Add(Vector<double>{0,0,0, 0,0,0.6, 1.8,-0.4,0.5, -0.5,2.0,-1.0, 3,3,3}, false);
  return trg;
}

/**
 * Errors of the three identities (rows: D[c], Green's identity with the source
 * outside, with the source inside) for each group of targets (columns: on the
 * surface, then the groups of Targets), relative to the max-norm of c and of
 * u on the surface.
 */
template <class KerSL, class KerDL, class KerGrad> Matrix<double> IdentityErrors(const SlenderElemList<double>& elem_lst, const Targets& trg, const Vector<double>& X0_out, const Vector<double>& X0_in, const Vector<double>& F0, const double tol) {
  constexpr Integer KDIM = KerSL::SrcDim();
  RequireQuadTables<KerSL>(FourierOrder, ElemOrder);
  RequireQuadTables<KerDL>(FourierOrder, ElemOrder);

  const KerSL ker_sl;
  const KerDL ker_dl;
  const KerGrad ker_grad;
  BoundaryIntegralOp<double,KerSL> SL(ker_sl, false, Comm::Self());
  BoundaryIntegralOp<double,KerDL> DL(ker_dl, false, Comm::Self());
  SL.AddElemList(elem_lst);
  DL.AddElemList(elem_lst);
  SL.SetAccuracy(tol);
  DL.SetAccuracy(tol);

  Vector<double> X, Xn;
  elem_lst.GetNodeCoord(&X, &Xn, nullptr);
  const Long Nsurf = X.Dim()/COORD_DIM;
  const Long Ntrg = trg.X.Dim()/COORD_DIM;
  const Long Ngrp = trg.grp_name.Dim();

  // Densities: constant c = F0; for the point source at X0 with strength F0,
  // u (for the double-layer) and du/dn (for the single-layer) on the surface.
  Vector<double> c(Nsurf*KDIM);
  for (Long i = 0; i < Nsurf; i++) {
    for (Integer k = 0; k < KDIM; k++) c[i*KDIM+k] = F0[k];
  }
  const Vector<double> Xn0{0,0,0};
  const auto point_source = [&](Vector<double>& u, Vector<double>& du_dn, const Vector<double>& X0) {
    Vector<double> du;
    ker_sl  .Eval(u , X, X0, Xn0, F0);
    ker_grad.Eval(du, X, X0, Xn0, F0);
    du_dn.ReInit(Nsurf*KDIM);
    for (Long i = 0; i < Nsurf; i++) { // du_dn <-- du . Xn
      for (Integer j = 0; j < KDIM; j++) {
        double du_dot_n = 0;
        for (Integer k = 0; k < COORD_DIM; k++) du_dot_n += du[(i*KDIM+j)*COORD_DIM+k] * Xn[i*COORD_DIM+k];
        du_dn[i*KDIM+j] = du_dot_n;
      }
    }
  };
  Vector<double> u_out, du_dn_out, u_in, du_dn_in, u_out_trg, u_in_trg;
  point_source(u_out, du_dn_out, X0_out);
  point_source(u_in , du_dn_in , X0_in );
  ker_sl.Eval(u_out_trg, trg.X, X0_out, Xn0, F0);
  ker_sl.Eval(u_in_trg , trg.X, X0_in , Xn0, F0);

  // Layer potentials on the surface (Usurf) and at the targets (Utrg)
  Vector<double> Usl, Udl;
  Vector<Vector<double>> Usurf(3), Utrg(3);
  DL.ComputePotential(Usurf[0], c);
  SL.ComputePotential(Usl, du_dn_out); DL.ComputePotential(Udl, u_out); Usurf[1] = Usl - Udl;
  SL.ComputePotential(Usl, du_dn_in ); DL.ComputePotential(Udl, u_in ); Usurf[2] = Usl - Udl;
  SL.SetTargetCoord(trg.X);
  DL.SetTargetCoord(trg.X);
  DL.ComputePotential(Utrg[0], c);
  SL.ComputePotential(Usl, du_dn_out); DL.ComputePotential(Udl, u_out); Utrg[1] = Usl - Udl;
  SL.ComputePotential(Usl, du_dn_in ); DL.ComputePotential(Udl, u_in ); Utrg[2] = Usl - Udl;

  // Exact values
  Vector<Vector<double>> Usurf_ref(3), Utrg_ref(3);
  Usurf_ref[0] = c * (-0.5);
  Usurf_ref[1] = u_out * 0.5;
  Usurf_ref[2] = u_in * (-0.5);
  for (auto& U : Utrg_ref) {
    U.ReInit(Ntrg*KDIM);
    U = 0;
  }
  for (Long i = 0; i < Ntrg; i++) {
    for (Integer k = 0; k < KDIM; k++) {
      const Long idx = i*KDIM+k;
      if (trg.inside[i]) {
        Utrg_ref[0][idx] = -F0[k];
        Utrg_ref[1][idx] = u_out_trg[idx];
      } else {
        Utrg_ref[2][idx] = -u_in_trg[idx];
      }
    }
  }

  const double scal[3] = {MaxNorm(F0), MaxNorm(u_out), MaxNorm(u_in)};
  Matrix<double> err(3, Ngrp+1);
  for (Integer l = 0; l < 3; l++) {
    err[l][0] = MaxErr(Usurf[l], Usurf_ref[l], scal[l]);
    for (Long g = 0; g < Ngrp; g++) {
      const Long offset = trg.grp_dsp[g]*KDIM, len = (trg.grp_dsp[g+1]-trg.grp_dsp[g])*KDIM;
      const Vector<double> U(len, Utrg[l].begin()+offset, false);
      const Vector<double> U_ref(len, Utrg_ref[l].begin()+offset, false);
      err[l][g+1] = MaxErr(U, U_ref, scal[l]);
    }
  }
  return err;
}

/**
 * Check the errors of IdentityErrors for the accuracy tolerance tol.
 */
void CheckIdentities(const std::string& ker_name, const Matrix<double>& err, const Targets& trg, const double tol) {
  const std::string identity_name[3] = {"D[c]", "Green's identity, source outside", "Green's identity, source inside"};
  for (Integer l = 0; l < 3; l++) {
    for (Long g = 0; g < err.Dim(1); g++) {
      const std::string grp_name = (g ? trg.grp_name[g-1] : std::string("on surface"));
      Check(ker_name + ", " + identity_name[l] + ": " + grp_name, err[l][g], ErrFactor*tol);
    }
  }
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  int status;

  {
    SlenderElemList<double> elem_lst;
    BuildElemList(elem_lst, LoopGeom, PanelBreakpoints(), ElemOrder, FourierOrder);
    const Targets trg = BuildTargets(elem_lst);

    Vector<double> X0_in(COORD_DIM); // on the centerline at s = 0
    double r0;
    LoopGeom(X0_in[0], X0_in[1], X0_in[2], r0, 0);
    const Vector<double> X0_out{0.2, 1.8, 1.3}; // about 1.5 from the loop

    const double tol = 1e-10;
    const Matrix<double> err_laplace = IdentityErrors<Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(elem_lst, trg, X0_out, X0_in, Vector<double>{1.0}, tol);
    CheckIdentities("Laplace", err_laplace, trg, tol);

    const Matrix<double> err_stokes = IdentityErrors<Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(elem_lst, trg, X0_out, X0_in, Vector<double>{0.6, -0.8, 1.2}, tol);
    CheckIdentities("Stokes", err_stokes, trg, tol);

    for (const double tol_ : {1e-4, 1e-6, 1e-8}) { // Laplace: the error follows the accuracy tolerance
      const Matrix<double> err = IdentityErrors<Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(elem_lst, trg, X0_out, X0_in, Vector<double>{1.0}, tol_);
      double max_err = 0;
      for (const auto e : err) max_err = std::max(max_err, e);
      std::stringstream ss;
      ss<<"Laplace, tolerance "<<std::scientific<<std::setprecision(0)<<tol_<<": all targets";
      Check(ss.str(), max_err, ErrFactor*tol_);
    }

    status = Summary("layer-potentials");
  }

  Comm::MPI_Finalize();
  return status;
}

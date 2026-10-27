// Copyright (C) 2021-2026 The Simons Foundation, Inc. and the CSBQ authors (see NOTICE)
// Author: Dhairya Malhotra
// SPDX-License-Identifier: Apache-2.0

/**
 * Unit test for the invariance of the layer potentials under a change of the
 * length scale: the geometry and the target points scaled by a factor h, with
 * the same density, give single-layer potentials scaled by h and unchanged
 * double-layer potentials (Laplace and Stokes kernels).
 *
 * This is a regression test for the ring integral (toroidal_greens_fn_batched,
 * fixed in bf26f80): a target whose squared distance from the axis of a ring
 * was below machine epsilon, an absolute length, was given the reference angle
 * of the element. For cross-sectional radii below about 1e-8 that was every
 * target, so the self- and near-interactions were wrong at all other angles.
 * Here h = 1e-10 (radius 1e-11), and the targets are the surface nodes and
 * points near the surface at angles away from the reference direction (the
 * orientation vector, the default one or one passed to Init).
 *
 * The program returns a non-zero exit status if a check fails. To compile and
 * run all the unit tests, start in the project root directory and run:
 * make check
 */

#include "common.hpp"

#include <cmath>
#include <functional>

using namespace sctl;
using namespace sctl::unit_test;

constexpr Integer COORD_DIM = 3;
constexpr Long ElemOrder = 10;
constexpr Long FourierOrder = 12;

/**
 * A closed loop discretized at length scale 1, and target points near its
 * surface.
 */
struct Loop {
  std::string name;
  Vector<Long> elem_order, fourier_order;
  Vector<double> Xc, radius, orientation; // as passed to Init (orientation may be empty)
  Vector<double> X; // surface nodes
  Vector<double> Xt; // target points
};

/**
 * Discretize the loop geom_fn(x, y, z, r, s), s in [0,1], with Nelem elements,
 * with the orientation vectors orient_fn(o, s) (orthogonal to the tangent) if
 * given, else the default ones. The targets are at distances -0.5, -0.05, 0.05
 * and 0.5 times the radius from the surface, at angles theta from the
 * orientation vector between pi/2 and 3 pi/2.
 */
template <class GeomFn, class OrientFn> Loop BuildLoop(const std::string& name, const GeomFn& geom_fn, const Long Nelem, const OrientFn* orient_fn) {
  Loop loop;
  loop.name = name;
  Vector<double> s_nds;
  SlenderElemList<double> elem_lst;
  BuildElemList(elem_lst, geom_fn, UniformPanels<double>(Nelem), ElemOrder, FourierOrder, &s_nds, &loop.Xc, &loop.radius);
  loop.elem_order.ReInit(Nelem);
  loop.fourier_order.ReInit(Nelem);
  loop.elem_order = ElemOrder;
  loop.fourier_order = FourierOrder;

  if (orient_fn) {
    for (Long i = 0; i < s_nds.Dim(); i++) {
      StaticArray<double,COORD_DIM> o;
      (*orient_fn)(o, s_nds[i]);
      for (Integer k = 0; k < COORD_DIM; k++) loop.orientation.PushBack(o[k]);
    }
    elem_lst.Init(loop.elem_order, loop.fourier_order, loop.Xc, loop.radius, loop.orientation);
  }
  elem_lst.GetNodeCoord(&loop.X, nullptr, nullptr);

  const double pi = const_pi<double>();
  for (const Long e : {(Long)0, Nelem/3, (2*Nelem)/3}) {
    OffSurfacePoints(loop.Xt, elem_lst, e, Vector<double>{0.3, 0.75}, Vector<double>{pi/2, 2*pi/3, pi, 3*pi/2}, Vector<double>{-0.5, -0.05, 0.05, 0.5});
  }
  return loop;
}

/**
 * A smooth density with KDIM components at the points X.
 */
Vector<double> Density(const Vector<double>& X, const Integer KDIM) {
  const Long N = X.Dim()/COORD_DIM;
  Vector<double> F(N*KDIM);
  for (Long i = 0; i < N; i++) {
    const double x = X[i*COORD_DIM+0], y = X[i*COORD_DIM+1], z = X[i*COORD_DIM+2];
    for (Integer k = 0; k < KDIM; k++) {
      F[i*KDIM+k] = cos(2*x + k + 1) * sin(y - (k+1)*z) + 0.5;
    }
  }
  return F;
}

/**
 * Compare the potentials of the loop on its surface and at its targets at
 * length scales 1 and h: U_h / h^length_pow must equal U_1.
 */
template <class Kernel> void TestScaleInvariance(const std::string& ker_name, const Loop& loop, const double h, const Integer length_pow, const double err_tol) {
  RequireQuadTables<Kernel>(FourierOrder, ElemOrder);
  const Vector<double> F = Density(loop.X, Kernel::SrcDim());

  const auto compute_potential = [&loop,&F](Vector<double>& U_surf, Vector<double>& U_trg, const double scal) {
    const SlenderElemList<double> elem_lst(loop.elem_order, loop.fourier_order, loop.Xc*scal, loop.radius*scal, loop.orientation);
    const Kernel ker;
    BoundaryIntegralOp<double,Kernel> BIOp(ker, false, Comm::Self());
    BIOp.AddElemList(elem_lst);
    BIOp.SetAccuracy(1e-10);
    BIOp.ComputePotential(U_surf, F);
    BIOp.SetTargetCoord(loop.Xt*scal);
    BIOp.ComputePotential(U_trg, F);
  };
  Vector<double> U_surf1, U_trg1, U_surf, U_trg;
  compute_potential(U_surf1, U_trg1, 1);
  compute_potential(U_surf, U_trg, h);

  const double scal = std::pow(h, (double)length_pow);
  Check(loop.name + ", " + ker_name + ": on surface", MaxErr(U_surf*(1/scal), U_surf1), err_tol);
  Check(loop.name + ", " + ker_name + ": near surface", MaxErr(U_trg*(1/scal), U_trg1), err_tol);
}

template <class GeomFn, class OrientFn> void TestLoop(const std::string& name, const GeomFn& geom_fn, const Long Nelem, const OrientFn* orient_fn) {
  constexpr double h = 1e-10;
  constexpr double err_tol = 1e-11;
  const Loop loop = BuildLoop(name, geom_fn, Nelem, orient_fn);
  TestScaleInvariance<Laplace3D_FxU>("Laplace single-layer", loop, h, 1, err_tol);
  TestScaleInvariance<Laplace3D_DxU>("Laplace double-layer", loop, h, 0, err_tol);
  TestScaleInvariance<Stokes3D_FxU >("Stokes single-layer" , loop, h, 1, err_tol);
  TestScaleInvariance<Stokes3D_DxU >("Stokes double-layer" , loop, h, 0, err_tol);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  int status;

  {
    const double pi = const_pi<double>();
    using OrientFn = std::function<void(StaticArray<double,COORD_DIM>&, const double)>;

    { // Ring in the xy-plane with the default orientation vector (the z-axis)
      const auto ring = [pi](double& x, double& y, double& z, double& r, const double s) {
        x = cos(2*pi*s);
        y = sin(2*pi*s);
        z = 0;
        r = 0.1;
      };
      TestLoop("ring, default orientation", ring, 8, (OrientFn*)nullptr);
    }

    { // Non-planar loop with varying radius, with an orientation vector passed to Init: the z-axis without its component along the tangent
      const auto loop = [pi](double& x, double& y, double& z, double& r, const double s) {
        const double t = 2*pi*s;
        x = cos(t);
        y = sin(t);
        z = 0.1*sin(2*t);
        r = 0.1*(1 + 0.25*cos(t+1));
      };
      const OrientFn orient = [pi](StaticArray<double,COORD_DIM>& o, const double s) {
        const double t = 2*pi*s;
        const double dx[COORD_DIM] = {-sin(t), cos(t), 0.2*cos(2*t)};
        const double dx2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
        double o_norm = 0;
        for (Integer k = 0; k < COORD_DIM; k++) {
          o[k] = (k == 2 ? 1 : 0) - dx[k] * dx[2] / dx2;
          o_norm += o[k]*o[k];
        }
        for (Integer k = 0; k < COORD_DIM; k++) o[k] /= sqrt(o_norm);
      };
      TestLoop("loop, orientation given", loop, 8, &orient);
    }

    status = Summary("scale-invariance");
  }

  Comm::MPI_Finalize();
  return status;
}

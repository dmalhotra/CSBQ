// Copyright (C) 2021-2026 The Simons Foundation, Inc. and the CSBQ authors (see NOTICE)
// Author: Dhairya Malhotra
// SPDX-License-Identifier: Apache-2.0

/**
 * Unit tests for the geometry of SlenderElemList:
 *
 * 1) surface area, enclosed volume (int X.n dA / 3) and normals of a torus, a
 * tube around a trefoil knot, a torus with varying cross-section and a
 * prolate spheroid (closed ends), against exact values.
 *
 * 2) round trip of the geometry file through Write and Read.
 *
 * The program returns a non-zero exit status if a check fails. To compile and
 * run all the unit tests, start in the project root directory and run:
 * make check
 */

#include "common.hpp"

#include <cstdio>
#include <filesystem>
#include <functional>
#include <sstream>
#include <unistd.h>

using namespace sctl;
using namespace sctl::unit_test;

constexpr Integer COORD_DIM = 3;
using Vec3 = Tensor<double,true,COORD_DIM,1>;

double Dot(const Vec3& u, const Vec3& v) {
  return u(0,0)*v(0,0) + u(1,0)*v(1,0) + u(2,0)*v(2,0);
}

/**
 * A closed body given by its centerline x(s) and cross-sectional radius r(s)
 * for s in [0,1], with the derivatives in s.
 */
struct Shape {
  std::string name;
  std::function<void(Vec3& x, Vec3& dx, Vec3& d2x, double& r, double& dr, const double s)> geom;
  Long Nelem, FourierOrder;
  double area, volume; // exact values
  std::function<Vec3(const Vec3& X)> normal; // exact normal at a surface point (optional)
};

/**
 * Surface area, enclosed volume (divergence theorem: V = int X.n dA / 3) and
 * the max-norm of int n dA, by the far-field quadrature rule of the element
 * list (Gauss-Legendre in s, trapezoidal in theta, weights include the area
 * element).
 */
void SurfaceIntegrals(double& area, double& volume, double& normal_integral, const SlenderElemList<double>& elem_lst) {
  Vector<double> X, Xn, wts, dist_far;
  Vector<Long> node_cnt;
  elem_lst.GetFarFieldNodes(X, Xn, wts, dist_far, node_cnt, 1e-10);

  StaticArray<double,COORD_DIM> n_int{0,0,0};
  area = 0;
  volume = 0;
  for (Long i = 0; i < wts.Dim(); i++) {
    area += wts[i];
    for (Integer k = 0; k < COORD_DIM; k++) {
      volume += X[i*COORD_DIM+k] * Xn[i*COORD_DIM+k] * wts[i] / 3;
      n_int[k] += Xn[i*COORD_DIM+k] * wts[i];
    }
  }
  normal_integral = std::max(fabs(n_int[0]), std::max(fabs(n_int[1]), fabs(n_int[2])));
}

/**
 * Checks of the surface nodes and normals of an element list built by
 * BuildElemList from shape.geom, with the centerline parameter s_nds and the
 * coordinates Xc at the centerline nodes. A surface node X at centerline point
 * x(s) must be at distance r(s) from it, and its normal n must be a unit vector
 * pointing outward ((X-x).n > 0) and equal to the exact normal: for a tube
 * with tangent t = x'/|x'| and u = (X-x)/r,
 *
 *   n = (|x'| a u - r' t) / sqrt(|x'|^2 a^2 + r'^2),  a = 1 - r (u.x'') / |x'|^2,
 *
 * or shape.normal(X) when it is given.
 */
void CheckNodesNormals(const Shape& shape, const SlenderElemList<double>& elem_lst, const Vector<double>& s_nds, const Vector<double>& Xc, const Long ElemOrder) {
  Vector<double> X, Xn;
  Vector<Long> node_cnt;
  elem_lst.GetNodeCoord(&X, &Xn, &node_cnt);

  double max_r_err = 0, max_n_err = 0, max_unit_err = 0, min_outward = 1;
  Long node_idx = 0;
  for (Long e = 0; e < elem_lst.Size(); e++) {
    const Long Nt = node_cnt[e] / ElemOrder; // Fourier order of the element
    for (Long j = 0; j < ElemOrder; j++) {
      const Long c_idx = e*ElemOrder+j; // centerline node
      Vec3 x, dx, d2x;
      double r, dr;
      shape.geom(x, dx, d2x, r, dr, s_nds[c_idx]);
      for (Long k = 0; k < Nt; k++) {
        Vec3 y, n;
        for (Integer l = 0; l < COORD_DIM; l++) {
          y(l,0) = X[node_idx*COORD_DIM+l] - Xc[c_idx*COORD_DIM+l];
          n(l,0) = Xn[node_idx*COORD_DIM+l];
        }
        const double y_norm = sqrt(Dot(y,y));
        const Vec3 u = y * (1/y_norm);

        Vec3 n_exact;
        if (shape.normal) {
          Vec3 X_;
          for (Integer l = 0; l < COORD_DIM; l++) X_(l,0) = X[node_idx*COORD_DIM+l];
          n_exact = shape.normal(X_);
        } else {
          const double dx_norm = sqrt(Dot(dx,dx));
          const double a = 1 - r * Dot(u,d2x) / (dx_norm*dx_norm);
          n_exact = (u*(dx_norm*a) - dx*(dr/dx_norm)) * (1/sqrt(dx_norm*dx_norm*a*a + dr*dr));
        }

        const Vec3 n_err = n - n_exact;
        max_r_err = std::max(max_r_err, fabs(y_norm - r) / r);
        max_n_err = std::max(max_n_err, sqrt(Dot(n_err,n_err)));
        max_unit_err = std::max(max_unit_err, fabs(sqrt(Dot(n,n)) - 1));
        min_outward = std::min(min_outward, Dot(n,u));
        node_idx++;
      }
    }
  }
  SCTL_ASSERT(node_idx*COORD_DIM == X.Dim());

  Check(shape.name + ": distance of nodes from the centerline", max_r_err, 1e-12);
  Check(shape.name + ": normals are unit vectors", max_unit_err, 1e-14);
  Check(shape.name + ": normals point outward, (X-x).n > 0", (min_outward > 0 ? 0 : 1), 0);
  Check(shape.name + ": normals against exact normals", max_n_err, 1e-10);
}

/**
 * Area, volume and normals of shape, discretized with elements of order 10.
 */
void TestShape(const Shape& shape) {
  constexpr Long ElemOrder = 10;
  const auto geom_fn = [&shape](double& x, double& y, double& z, double& r, const double s) {
    Vec3 X, dX, d2X;
    double dr;
    shape.geom(X, dX, d2X, r, dr, s);
    x = X(0,0);
    y = X(1,0);
    z = X(2,0);
  };

  SlenderElemList<double> elem_lst;
  Vector<double> s_nds, Xc;
  BuildElemList(elem_lst, geom_fn, UniformPanels<double>(shape.Nelem), ElemOrder, shape.FourierOrder, &s_nds, &Xc);

  double area, volume, normal_integral;
  SurfaceIntegrals(area, volume, normal_integral, elem_lst);
  Check(shape.name + ": area", fabs(area - shape.area) / shape.area, 1e-12);
  Check(shape.name + ": volume", fabs(volume - shape.volume) / shape.volume, 1e-10); // depends on x'' through the curvature term, so less accurate than the area
  Check(shape.name + ": int n dA / area", normal_integral / shape.area, 1e-13);
  CheckNodesNormals(shape, elem_lst, s_nds, Xc, ElemOrder);
}

/**
 * Max-norm of the difference of the surface nodes and normals of two element
 * lists, or infinity if their elements differ in number or size.
 */
double GeomDiff(const SlenderElemList<double>& elem_lst0, const SlenderElemList<double>& elem_lst1) {
  Vector<double> X0, Xn0, X1, Xn1;
  Vector<Long> cnt0, cnt1;
  elem_lst0.GetNodeCoord(&X0, &Xn0, &cnt0);
  elem_lst1.GetNodeCoord(&X1, &Xn1, &cnt1);
  if (elem_lst0.Size() != elem_lst1.Size() || cnt0.Dim() != cnt1.Dim()) return std::numeric_limits<double>::infinity();
  for (Long i = 0; i < cnt0.Dim(); i++) {
    if (cnt0[i] != cnt1[i]) return std::numeric_limits<double>::infinity();
  }
  return std::max(MaxErr(X1, X0, 1.0), MaxErr(Xn1, Xn0, 1.0));
}

std::string ReadFileToString(const std::string& fname) {
  std::ifstream file(fname);
  std::stringstream ss;
  ss << file.rdbuf();
  return ss.str();
}

/**
 * Write an element list to a file and read it back: the geometry must be
 * reproduced exactly, and writing the list read back must give the same file.
 */
void TestWriteRead(const std::string& name, const SlenderElemList<double>& elem_lst) {
  static Long file_cnt = 0;
  const std::string fname0 = (std::filesystem::temp_directory_path() / ("csbq-unit-geometry-" + std::to_string(getpid()) + "-" + std::to_string(file_cnt++) + ".geom")).string();
  const std::string fname1 = fname0 + "1";

  elem_lst.Write(fname0);
  SlenderElemList<double> elem_lst_;
  elem_lst_.Read<double>(fname0);
  elem_lst_.Write(fname1);
  const std::string file0 = ReadFileToString(fname0);
  const std::string file1 = ReadFileToString(fname1);
  std::remove(fname0.c_str());
  std::remove(fname1.c_str());

  Check(name + ": Write/Read reproduces the geometry", (elem_lst.Size() ? GeomDiff(elem_lst, elem_lst_) : 1.0), 0);
  Check(name + ": writing it again gives the same file", (file0.size() && file0 == file1 ? 0 : 1), 0);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  int status;

  {
    const double pi = const_pi<double>();

    { // Torus: circle of radius R, constant radius r0
      const double R = 1.0, r0 = 0.1;
      Shape torus;
      torus.name = "torus";
      torus.geom = [R,r0,pi](Vec3& x, Vec3& dx, Vec3& d2x, double& r, double& dr, const double s) {
        const double t = 2*pi*s, dt = 2*pi;
        x   = Vec3{R*cos(t), R*sin(t), 0};
        dx  = Vec3{-R*sin(t), R*cos(t), 0} * dt;
        d2x = Vec3{-R*cos(t), -R*sin(t), 0} * (dt*dt);
        r = r0;
        dr = 0;
      };
      torus.Nelem = 16;
      torus.FourierOrder = 12;
      torus.area = 4*pi*pi*R*r0;
      torus.volume = 2*pi*pi*R*r0*r0;
      TestShape(torus);
    }

    { // Tube around a (2,3) torus knot (trefoil), constant radius r0: area = 2 pi r0 L, volume = pi r0^2 L
      const double R = 1.0, a = 0.5, r0 = 0.05;
      const auto knot = [R,a,pi](Vec3& x, Vec3& dx, Vec3& d2x, const double s) {
        const double t = 2*pi*s, dt = 2*pi;
        const double c2 = cos(2*t), s2 = sin(2*t), c3 = cos(3*t), s3 = sin(3*t);
        const double rho = R + a*c3, drho = -3*a*s3, d2rho = -9*a*c3;
        x   = Vec3{rho*c2, rho*s2, a*s3};
        dx  = Vec3{drho*c2 - 2*rho*s2, drho*s2 + 2*rho*c2, 3*a*c3} * dt;
        d2x = Vec3{d2rho*c2 - 4*drho*s2 - 4*rho*c2, d2rho*s2 + 4*drho*c2 - 4*rho*s2, -9*a*s3} * (dt*dt);
      };
      double L = 0; // length of the knot (trapezoidal rule, spectrally accurate for a periodic integrand)
      const Long N = 4096;
      for (Long i = 0; i < N; i++) {
        Vec3 x, dx, d2x;
        knot(x, dx, d2x, i/(double)N);
        L += sqrt(Dot(dx,dx)) / N;
      }

      Shape trefoil;
      trefoil.name = "trefoil tube";
      trefoil.geom = [knot,r0](Vec3& x, Vec3& dx, Vec3& d2x, double& r, double& dr, const double s) {
        knot(x, dx, d2x, s);
        r = r0;
        dr = 0;
      };
      trefoil.Nelem = 32;
      trefoil.FourierOrder = 12;
      trefoil.area = 2*pi*r0*L;
      trefoil.volume = pi*r0*r0*L;
      TestShape(trefoil);
    }

    { // Torus with varying cross-section: r(t) = r0 (1 + cos(t)/2) on a circle of radius R
      const double R = 1.0, r0 = 0.1;
      const auto radius = [r0](double& r, double& dr_dt, const double t) {
        r = r0 * (1 + 0.5*cos(t));
        dr_dt = -0.5 * r0 * sin(t);
      };
      double area = 0; // area = int int r sqrt(R^2 (1 + r cos(phi) / R)^2 + r'(t)^2) dphi dt (trapezoidal rule)
      const Long N = 256;
      for (Long i = 0; i < N; i++) {
        double r, dr;
        radius(r, dr, 2*pi*i/N);
        for (Long j = 0; j < N; j++) {
          const double a = 1 + r*cos(2*pi*j/N)/R;
          area += r * sqrt(R*R*a*a + dr*dr) * (2*pi/N) * (2*pi/N);
        }
      }

      Shape tapered;
      tapered.name = "torus, varying radius";
      tapered.geom = [R,radius,pi](Vec3& x, Vec3& dx, Vec3& d2x, double& r, double& dr, const double s) {
        const double t = 2*pi*s, dt = 2*pi;
        x   = Vec3{R*cos(t), R*sin(t), 0};
        dx  = Vec3{-R*sin(t), R*cos(t), 0} * dt;
        d2x = Vec3{-R*cos(t), -R*sin(t), 0} * (dt*dt);
        radius(r, dr, t);
        dr *= dt;
      };
      tapered.Nelem = 16;
      tapered.FourierOrder = 12;
      tapered.area = area;
      tapered.volume = 2.25*pi*pi*R*r0*r0; // int pi r(t)^2 R dt
      TestShape(tapered);
    }

    { // Prolate spheroid (closed ends) with semi-axes a along the axis d through c, and b
      const double a = 1.0, b = 0.25, e = sqrt(1 - b*b/(a*a));
      const Vec3 d = Vec3{1, 2, 2} * (1/3.0), c{0.1, -0.2, 0.3};

      Shape spheroid;
      spheroid.name = "prolate spheroid";
      spheroid.geom = [a,b,d,c,pi](Vec3& x, Vec3& dx, Vec3& d2x, double& r, double& dr, const double s) {
        x   = c - d * (a*cos(pi*s));
        dx  = d * (a*pi*sin(pi*s));
        d2x = d * (a*pi*pi*cos(pi*s));
        r   = b*sin(pi*s);
        dr  = b*pi*cos(pi*s);
      };
      spheroid.normal = [a,b,d,c](const Vec3& X) { // gradient of ((X-c).d)^2/a^2 + |(X-c) - ((X-c).d) d|^2/b^2
        const Vec3 Y = X - c;
        const double Yd = Dot(Y,d);
        const Vec3 n = d*(Yd/(a*a)) + (Y - d*Yd)*(1/(b*b));
        return n * (1/sqrt(Dot(n,n)));
      };
      spheroid.Nelem = 16;
      spheroid.FourierOrder = 12;
      spheroid.area = 2*pi*b*b * (1 + a/(b*e) * asin(e));
      spheroid.volume = 4*pi*a*b*b/3;
      TestShape(spheroid);
    }

    { // Write/Read: elements of different orders, with an orientation vector passed to Init
      const Vector<Long> elem_order{10, 12, 8, 14, 10, 6};
      const Vector<Long> fourier_order{12, 8, 16, 20, 4, 14};
      Vector<double> Xc, radius, orientation;
      double s0 = 0;
      for (Long e = 0; e < elem_order.Dim(); e++) {
        const auto& nds = SlenderElemList<double>::CenterlineNodes(elem_order[e]);
        for (Long j = 0; j < elem_order[e]; j++) {
          const double t = 2*const_pi<double>() * (s0 + nds[j]) / elem_order.Dim();
          const Vec3 x{cos(t), sin(t), 0.1*sin(2*t)}, dx{-sin(t), cos(t), 0.2*cos(2*t)};
          const Vec3 v{0.3, -0.5, 0.8}; // orientation: v without its component along the tangent
          Vec3 o = v - dx * (Dot(v,dx)/Dot(dx,dx));
          o = o * (1/sqrt(Dot(o,o)));
          for (Integer k = 0; k < COORD_DIM; k++) {
            Xc.PushBack(x(k,0));
            orientation.PushBack(o(k,0));
          }
          radius.PushBack(0.05*(2+sin(t+sqrt(2.0))));
        }
        s0 += 1;
      }
      SlenderElemList<double> elem_lst(elem_order, fourier_order, Xc, radius, orientation);
      TestWriteRead("mixed orders, orientation given", elem_lst);
    }

    { // Write/Read: default orientation
      const auto geom_fn = [](double& x, double& y, double& z, double& r, const double s) {
        const double t = 2*const_pi<double>()*s;
        x = cos(t);
        y = sin(t);
        z = 0.1*sin(2*t);
        r = 0.01*(2+sin(t+sqrt(2.0)));
      };
      SlenderElemList<double> elem_lst;
      BuildElemList(elem_lst, geom_fn, UniformPanels<double>(16), 10, 8);
      TestWriteRead("default orientation", elem_lst);
    }

    { // Write/Read: geometry file in the data directory
      SlenderElemList<double> elem_lst;
      elem_lst.Read<double>(std::string(CSBQ_QUOTEME(CSBQ_DATA_PATH)) + "/loop.geom");
      Check("data/loop.geom: read", (elem_lst.Size() > 0 ? 0 : 1), 0);
      TestWriteRead("data/loop.geom", elem_lst);
    }

    status = Summary("geometry");
  }

  Comm::MPI_Finalize();
  return status;
}

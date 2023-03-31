#include "utils.hpp"
#include "sctl.hpp"
using namespace sctl;

const int ChebOrder = 10; // do not change this
const int FourierOrder = 16; // must be multiple of 4

/**
 * Solve Stokes Dirichlet BVP on a half-turn helix and print force-per-length
 * along the centerline.
 *
 * \param[in] r0 cross-sectional radius of the fiber.
 */
template <class Real> void test(const Comm& comm, const Real r0, const bool ellipsoid) {
  static constexpr Integer COORD_DIM = 3;
  using Vec3 = Tensor<Real,true,COORD_DIM,1>;

  const Real quad_tol = 1e-14;
  const Real gmres_tol = 1e-13;
  const Long gmres_max_iter = 400;
  const Real pi = const_pi<Real>();
  const Real scale = 1 / (sqrt<Real>(2.0)*pi);
  Real s0 = 0, s1 = 1;

  Vector<Real> panel_len; // define length of panels (in parameter s)
  SlenderElemList<Real> elem_lst0; // Initialize geometry from geom function
  if (ellipsoid) { // ellipsoidally tapered half-turn helix
    const Long Npanels = 32;
    panel_len.ReInit(Npanels);
    panel_len = 1 / (Real)Npanels;

    // Geometry function for half-turn helix (ellipsoidally tapered)
    auto geom = [&r0,&scale,&pi](Real& x, Real& y, Real& z, Real& r, const Real s) {
      // input: s \in [0,1] is the parameterization along the centerline.
      // output: x, y, z, r : the position of the centerline and the cross-sectional radius.

      Real s_ = (1 - cos<Real>(s*pi)) / 2;
      r = r0 * sin<Real>(s*pi);

      x = (cos<Real>(s_*pi)) * scale;
      y = (sin<Real>(s_*pi)) * scale;
      z = 0;//          (s_*pi)  * scale;
    };

    GenericGeom(elem_lst0, geom, panel_len, Vector<Long>(), comm, ChebOrder, FourierOrder);
  } else {
    const Real end_len = r0/scale/sqrt<Real>(8.0); // length of end-caps (in parameter space)
    s0 = 0 - end_len; // range of s
    s1 = 1 + end_len;

    // Geometry function for half-turn helix (with hemispherical end-caps)
    auto geom = [&r0,&s0,&s1,&scale,&end_len,&pi](Real& x, Real& y, Real& z, Real& r, const Real s) {
      // input: s \in [0,1] is the parameterization along the centerline.
      // output: x, y, z, r : the position of the centerline and the cross-sectional radius.

      Real s_ = s0 + s;
      if (s_ < s0 + end_len) { // round the end at s=0
        r = sin<Real>((s_-s0)/end_len*pi/2)*r0;
        s_ = s0 + end_len - cos<Real>((s_-s0)/end_len*pi/2) * end_len*2/pi;
      } else if (s_ > s1 - end_len) { // round the end at s=1
        r = sin<Real>((s1-s_)/end_len*pi/2)*r0;
        s_ = s1 - end_len + cos<Real>((s1-s_)/end_len*pi/2) * end_len*2/pi;
      } else {
        r = r0;
      }

      x = (cos<Real>(s_*pi)) * scale;
      y = (sin<Real>(s_*pi)) * scale;
      z =           (s_*pi)  * scale;
    };

    { // dyadically refine panels at the ends
      const Real L = (s1-s0)-2*end_len;
      const Long min_depth = 4;
      const Long max_depth = std::max<Long>(min_depth, (Long)(log(end_len/L)/log(0.5))+3);

      for (Long i = 0; i < 4; i++) panel_len.PushBack(end_len/4);
      panel_len.PushBack(pow<Real>(0.5,max_depth)*L);
      for (Long d = max_depth; d >= min_depth+1; d--) panel_len.PushBack(pow<Real>(0.5,d)*L);
      for (Long i = 0; i < ((1<<min_depth)-2); i++) panel_len.PushBack(pow<Real>(0.5,min_depth)*L);
      for (Long d = min_depth+1; d <= max_depth; d++) panel_len.PushBack(pow<Real>(0.5,d)*L);
      panel_len.PushBack(pow<Real>(0.5,max_depth)*L);
      for (Long i = 0; i < 4; i++) panel_len.PushBack(end_len/4);
    }
    GenericGeom(elem_lst0, geom, panel_len, Vector<Long>(), comm, ChebOrder, FourierOrder);
  }
  const Long Npanels = panel_len.Dim();

  const auto cross_prod = [](const Vec3& u, const Vec3& v) -> Vec3 {
    Vec3 uxv;
    uxv(0,0) = u(1,0) * v(2,0) - u(2,0) * v(1,0);
    uxv(1,0) = u(2,0) * v(0,0) - u(0,0) * v(2,0);
    uxv(2,0) = u(0,0) * v(1,0) - u(1,0) * v(0,0);
    return uxv;
  };
  const auto dot_prod = [](const Vec3& u, const Vec3& v) {
    Real u_dot_v = 0;
    u_dot_v += u(0,0) * v(0,0);
    u_dot_v += u(1,0) * v(1,0);
    u_dot_v += u(2,0) * v(2,0);
    return u_dot_v;
  };
  const auto vec3_mag = [&dot_prod](const Vec3& u) {
    return sqrt<Real>(dot_prod(u,u));
  };

  Vector<Real> Ws(Npanels * ChebOrder); // weights to integrate along the centerline
  Vector<Real> Xc(Npanels * ChebOrder * COORD_DIM); // coordinates of centerline
  Vector<Real> dXc(Npanels * ChebOrder * COORD_DIM); // tangent vector to the centerline
  Vector<Real> U(Npanels * ChebOrder * FourierOrder * COORD_DIM); // fluid velocity on the surface (the boundary condition)
  Vector<Real> r(Npanels * ChebOrder * FourierOrder * COORD_DIM); // the vector X-Xc
  Vector<Real> Wa(Npanels * ChebOrder * FourierOrder); // weights to average force on surface to the centerline
  Vector<Real> PanelNodes = SlenderElemList<Real>::CenterlineNodes(ChebOrder); // Chebyshev discretization nodes
  Vector<Real> ChebQuadWts = ChebQuadRule<Real>::ComputeWts(ChebOrder); // Clenshaw-Curtis quadrature weights on [0,1]
  { // Set fluid velocity on the boundary U and the weights Wa
    Vector<Real> sin_theta(FourierOrder), cos_theta(FourierOrder);
    for (Long i = 0; i < FourierOrder; i++) {
      Real theta = 2 * const_pi<Real>() * i / FourierOrder;
      sin_theta[i] = sin<Real>(theta);
      cos_theta[i] = cos<Real>(theta);
    }

    for (Long panel_idx = 0; panel_idx < Npanels; panel_idx++) { // loop over panels
      Vector<Real> X, dX_ds, dX_dt;
      elem_lst0.SlenderElemList<Real>::GetGeom(&X, nullptr, nullptr, &dX_ds, &dX_dt, PanelNodes, sin_theta, cos_theta, panel_idx);

      Vector<Real> U_(ChebOrder * FourierOrder * COORD_DIM, U.begin() + panel_idx * ChebOrder * FourierOrder * COORD_DIM, false); // sub-array of U
      U_ = dX_dt / (r0*r0); // the tangential derivative of the surface in angular direction scaled by (1/r0^2)

      { // Set Xc, dXc
        Vector<Real> X0(COORD_DIM*ChebOrder), X0s(COORD_DIM*ChebOrder); X0 = 0; // Centerline coordinates and derivatives
        for (Long i = 0; i < ChebOrder; i++) {
          for (Long j = 0; j < FourierOrder; j++) {
            for (Long k = 0; k < COORD_DIM; k++) {
              X0[k*ChebOrder+i] += X[(i*FourierOrder+j)*COORD_DIM+k] / FourierOrder;
            }
          }
        }
        LagrangeInterp<Real>::Derivative(X0s, X0, PanelNodes);
        for (Long i = 0; i < ChebOrder; i++) {
          for (Long k = 0; k < COORD_DIM; k++) {
            Xc[(panel_idx*ChebOrder+i)*COORD_DIM+k] = X0[k*ChebOrder+i];
            dXc[(panel_idx*ChebOrder+i)*COORD_DIM+k] = X0s[k*ChebOrder+i];
          }
        }
      }

      for (Long i = 0; i < ChebOrder; i++) {
        const Vec3 dXc_(dXc.begin() + (panel_idx*ChebOrder+i)*COORD_DIM);
        Ws[panel_idx*ChebOrder+i] = vec3_mag(dXc_) * ChebQuadWts[i];
        for (Long j = 0; j < FourierOrder; j++) {
          const Long node_idx = i * FourierOrder + j;
          const auto& Xs = Vec3(dX_ds.begin() + node_idx * COORD_DIM);
          const auto& Xt = Vec3(dX_dt.begin() + node_idx * COORD_DIM);
          Wa[panel_idx * ChebOrder*FourierOrder + node_idx] = vec3_mag(cross_prod(Xs, Xt)) / vec3_mag(dXc_);

          Vec3 r_ = cross_prod(dXc_, Xt) / vec3_mag(dXc_);
          for (Long k = 0; k < COORD_DIM; k++) {
            r[(panel_idx * ChebOrder*FourierOrder + node_idx) * COORD_DIM + k] = r_(k,0);
          }
        }
      }
    }
  }

  // write visualization to VTK file
  elem_lst0.WriteVTK("vis/U-helix", U, comm);

  const Stokes3D_FxU ker_FxU; // single-layer kernel
  BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesFxU(ker_FxU, false, comm); // boundary integral operator
  BIOp_StokesFxU.SetFMMKer(ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU);
  BIOp_StokesFxU.AddElemList(elem_lst0);
  BIOp_StokesFxU.SetAccuracy(quad_tol);

  if (0) { // print quadrature error
    const Stokes3D_DxU ker_DxU;
    const Stokes3D_FSxU ker_FSxU;
    BoundaryIntegralOp<Real,Stokes3D_DxU> BIOp_StokesDxU(ker_DxU, false, comm); // boundary integral operator
    BIOp_StokesDxU.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
    BIOp_StokesDxU.AddElemList(elem_lst0);
    BIOp_StokesDxU.SetAccuracy(quad_tol);

    long N = BIOp_StokesDxU.Dim(0);
    Vector<Real> U, F(N); F = 1;
    BIOp_StokesDxU.ComputePotential(U, F);
    Real err = 0;
    for (long i = 0; i < N; i++) {
      U[i] += 0.5;
      err = std::max<Real>(err, fabs(U[i]));
    }
    elem_lst0.WriteVTK("vis/E-helix", U, comm);
    std::cout<<"Double-layer quadrature error = "<<err<<'\n';
  }

  Vector<Real> F;
  ParallelSolver<Real> solver(comm); // linear solver
  auto BIOp = [&BIOp_StokesFxU](Vector<Real>* Ax, const Vector<Real>& x) {
    Ax->SetZero();
    BIOp_StokesFxU.ComputePotential(*Ax, x);
  };
  solver(&F, BIOp, U, gmres_tol, gmres_max_iter); // solve for surface force F

  // write visualization to VTK file
  elem_lst0.WriteVTK("vis/F-helix", F, comm);

  Real s_offset = s0;
  Vector<Real> s_vec(Npanels * ChebOrder); // the parameterization s in [0,1]
  Vector<Real> F0(Npanels * ChebOrder * COORD_DIM); // force per-unit-length
  Vector<Real> T0(Npanels * ChebOrder); // parallel torque per-unit-length
  for (Long panel_idx = 0; panel_idx < Npanels; panel_idx++) { // loop over panels
    for (Long i = 0; i < ChebOrder; i++) {
      const Long idx = panel_idx*ChebOrder+i;
      s_vec[idx] = s_offset + PanelNodes[i]*panel_len[panel_idx];

      Vec3 F0_, T0_;
      F0_ = (Real)0;
      T0_ = (Real)0;
      for (Long j = 0; j < FourierOrder; j++) {
        const Long idx_ = idx * FourierOrder + j;
        const auto F_ = Vec3(F.begin() + idx_*COORD_DIM);
        const auto r_ = Vec3(r.begin() + idx_*COORD_DIM);
        T0_ = T0_ + cross_prod(r_, F_) * Wa[idx_];
        F0_ = F0_ + F_ * Wa[idx_];
      }
      for (Long k = 0; k < COORD_DIM; k++) {
        F0[idx * COORD_DIM + k] = F0_(k,0);
      }
      const auto dXc_ = Vec3(dXc.begin() + idx*COORD_DIM);
      T0[idx] = dot_prod(T0_, dXc_) / vec3_mag(dXc_);
    }
    s_offset += panel_len[panel_idx];
  }

  // print result
  printf("                   s                    X                    Y                    Z                   Fx                   Fy                   Fz      parallel-torque                   Ws\n");
  for (Long i = 0; i < s_vec.Dim(); i++) {
    printf("%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", s_vec[i], Xc[i*COORD_DIM+0], Xc[i*COORD_DIM+1], Xc[i*COORD_DIM+2], F0[i*COORD_DIM+0], F0[i*COORD_DIM+1], F0[i*COORD_DIM+2], T0[i], Ws[i]);
  }
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    Comm comm = Comm::Self();
    commandline_option_start(argc, argv, nullptr, comm);
    double eps = strtod(commandline_option(argc, argv, "-eps", "1e-2" , true, "the cross-sectional redius of the fiber", comm), nullptr);
    Long geom = strtol(commandline_option(argc, argv, "-geom", "1" , false, "geometry: 1) ellipsoidally tapered 2) hemispherical end-caps", comm), nullptr, 10);
    commandline_option_end(argc, argv);

    test<double>(comm, eps, (geom==1));
  }

  Comm::MPI_Finalize();
  return 0;
}


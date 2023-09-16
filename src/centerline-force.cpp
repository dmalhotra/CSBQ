#include "utils.hpp"
#include "sctl.hpp"
using namespace sctl;

const bool DL_fix = true;
const bool sqrt_scaling = true;
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
  Real s0 = 0;//, s1 = 1;

  Vector<Real> panel_len; // define length of panels (in parameter s)
  SlenderElemList<Real> elem_lst0; // Initialize geometry from geom function
  if (ellipsoid) { // ellipsoidally tapered half-turn helix
    using ValueType = long double;
    const ValueType pi = const_pi<ValueType>();
    const ValueType scale = 1 / (sqrt<ValueType>(2.0)*pi);

    Vector<ValueType> panel_len_;
    if (1) { // dyadic panel refinement
      const Integer MinDepth = 6;
      const Integer MaxDepth = 9;
      ValueType ds = pow<ValueType>(0.5,MaxDepth), ds_max = pow<ValueType>(0.5,MinDepth), sum = 0;
      panel_len_.PushBack(ds); sum += ds;
      while (sum < 1/(ValueType)2) {
        panel_len_.PushBack(ds); sum += ds;
        if (ds<ds_max) ds*=(ValueType)2;
      }
      SCTL_ASSERT(sum == 1/(ValueType)2);

      const Long N0 = panel_len_.Dim();
      for (Long i = 0; i < N0; i++) {
        panel_len_.PushBack(panel_len_[N0-1-i]);
      }
    }

    // Geometry function for half-turn helix (ellipsoidally tapered)
    auto geom = [&r0,&scale,&pi](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType s) {
      // input: s \in [0,1] is the parameterization along the centerline.
      // output: x, y, z, r : the position of the centerline and the cross-sectional radius.

      ValueType s_ = (1 - cos<ValueType>(s*pi)) / 2;
      r = r0 * sin<ValueType>(s*pi);

      x = (cos<ValueType>(s_*pi)) * scale;
      y = (sin<ValueType>(s_*pi)) * scale;
      z =                (s_*pi)  * scale;
    };

    panel_len.ReInit(0);
    for (const auto x : panel_len_) panel_len.PushBack((Real)x);
    GenericGeom(elem_lst0, geom, panel_len_, Vector<Long>(), comm, ChebOrder, FourierOrder);
  } else {
    using ValueType = long double;
    const ValueType pi = const_pi<ValueType>();
    const ValueType scale = 1 / (sqrt<ValueType>(2.0)*pi);

    const ValueType end_len = r0/scale/sqrt<ValueType>(8.0); // length of end-caps (in parameter space)
    s0 = 0 - (Real)end_len;
    //s1 = 1 + (Real)end_len;
    const ValueType s0 = (ValueType)0 - end_len; // range of s
    const ValueType s1 = (ValueType)1 + end_len;

    // Geometry function for half-turn helix (with hemispherical end-caps)
    auto geom = [&r0,&s0,&s1,&scale,&end_len,&pi](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType s) {
      // input: s \in [0,1] is the parameterization along the centerline.
      // output: x, y, z, r : the position of the centerline and the cross-sectional radius.

      ValueType s_ = s0 + s;
      if (s_ < s0 + end_len) { // round the end at s=0
        r = sin<ValueType>((s_-s0)/end_len*pi/2)*r0;
        s_ = s0 + end_len - cos<ValueType>((s_-s0)/end_len*pi/2) * end_len*2/pi;
      } else if (s_ > s1 - end_len) { // round the end at s=1
        r = sin<ValueType>((s1-s_)/end_len*pi/2)*r0;
        s_ = s1 - end_len + cos<ValueType>((s1-s_)/end_len*pi/2) * end_len*2/pi;
      } else {
        r = r0;
      }

      x = (cos<ValueType>(s_*pi)) * scale;
      y = (sin<ValueType>(s_*pi)) * scale;
      z =                (s_*pi)  * scale;
    };

    Vector<ValueType> panel_len_;
    { // dyadically refine panels at the ends
      const ValueType L = (s1-s0)-2*end_len;
      const Long min_depth = 4;
      const Long max_depth = std::max<Long>(min_depth, (Long)(log(end_len/L)/log(0.5))+3);

      for (Long i = 0; i < 4; i++) panel_len_.PushBack(end_len/4);
      panel_len_.PushBack(pow<ValueType>(0.5,max_depth)*L);
      for (Long d = max_depth; d >= min_depth+1; d--) panel_len_.PushBack(pow<ValueType>(0.5,d)*L);
      for (Long i = 0; i < ((1<<min_depth)-2); i++) panel_len_.PushBack(pow<ValueType>(0.5,min_depth)*L);
      for (Long d = min_depth+1; d <= max_depth; d++) panel_len_.PushBack(pow<ValueType>(0.5,d)*L);
      panel_len_.PushBack(pow<ValueType>(0.5,max_depth)*L);
      for (Long i = 0; i < 4; i++) panel_len_.PushBack(end_len/4);
    }

    panel_len.ReInit(0);
    for (const auto x : panel_len_) panel_len.PushBack((Real)x);
    GenericGeom<ValueType>(elem_lst0, geom, panel_len_, Vector<Long>(), comm, ChebOrder, FourierOrder);
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

  Vector<Real> Wa_surf(Npanels * ChebOrder * FourierOrder); // surface area element
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
      Vector<Real> X, dX_ds, dX_dt, Wa_surf_;
      elem_lst0.SlenderElemList<Real>::GetGeom(&X, nullptr, &Wa_surf_, &dX_ds, &dX_dt, PanelNodes, sin_theta, cos_theta, panel_idx);

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
          const Long idx = panel_idx * ChebOrder*FourierOrder + node_idx;

          //const auto Xs = Vec3(dX_ds.begin() + node_idx * COORD_DIM);
          const auto Xt = Vec3(dX_dt.begin() + node_idx * COORD_DIM);
          const Real PTR_weight = (2*pi/FourierOrder);
          Wa_surf[idx] = Wa_surf_[node_idx];
          Wa[idx] = Wa_surf_[node_idx] / vec3_mag(dXc_) * PTR_weight;
          //Wa[idx] = vec3_mag(cross_prod(Xs, Xt)) / vec3_mag(dXc_) * PTR_weight;

          Vec3 r_ = cross_prod(dXc_, Xt) / vec3_mag(dXc_);
          for (Long k = 0; k < COORD_DIM; k++) {
            r[idx * COORD_DIM + k] = r_(k,0);
          }
        }
      }
    }
  }

  // write visualization to VTK file
  elem_lst0.WriteVTK("vis/U-helix", U, comm);

  const Stokes3D_FxU ker_FxU; // single-layer kernel
  const Stokes3D_FSxU ker_FSxU; // stokeslet + source kernel
  BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesFxU(ker_FxU, false, comm); // boundary integral operator
  BIOp_StokesFxU.SetFMMKer(ker_FxU, ker_FxU, ker_FxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
  BIOp_StokesFxU.AddElemList(elem_lst0);
  BIOp_StokesFxU.SetAccuracy(quad_tol);

  const Stokes3D_DxU ker_DxU; // double-layer kernel
  BoundaryIntegralOp<Real,Stokes3D_DxU> BIOp_StokesDxU(ker_DxU, false, comm); // boundary integral operator (double-layer)
  BIOp_StokesDxU.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
  BIOp_StokesDxU.AddElemList(elem_lst0);
  BIOp_StokesDxU.SetAccuracy(quad_tol);

  if (1) { // print quadrature error
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

  Vector<Real> sqrt_Wa, invsqrt_Wa; // weights for sqrt scaling
  if (sqrt_scaling) {
    for (const auto& x : Wa_surf) sqrt_Wa.PushBack(sqrt<Real>(x));
    for (const auto& x : Wa_surf) invsqrt_Wa.PushBack(1/sqrt<Real>(x));
  }

  Vector<Real> F, DU;
  ParallelSolver<Real> solver(comm); // linear solver
  auto scaling = [](const Vector<Real>& V_, const Vector<Real>& scal) {
    Vector<Real> V = V_;
    const Long N = scal.Dim();
    const Long dof = V.Dim()/std::max<Long>(N,1);
    for (Long i = 0; i < N; i++) {
      for (Long k = 0; k < dof; k++) {
        V[i*dof+k] *= scal[i];
      }
    }
    return V;
  };
  auto BIOp = [&BIOp_StokesFxU, &scaling, &sqrt_Wa, &invsqrt_Wa](Vector<Real>* Ax, const Vector<Real>& x) {
    Ax->SetZero();
    BIOp_StokesFxU.ComputePotential(*Ax, scaling(x, invsqrt_Wa));
    *Ax = scaling(*Ax, sqrt_Wa);
  };
  BIOp_StokesDxU.ComputePotential(DU, U); DU += 0.5* U;
  solver(&F, BIOp, scaling(U-DU*(DL_fix?1:0), sqrt_Wa), gmres_tol, gmres_max_iter); // solve for surface force F
  F = scaling(F, invsqrt_Wa);

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
    printf("%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n", (double)s_vec[i], (double)Xc[i*COORD_DIM+0], (double)Xc[i*COORD_DIM+1], (double)Xc[i*COORD_DIM+2], (double)F0[i*COORD_DIM+0], (double)F0[i*COORD_DIM+1], (double)F0[i*COORD_DIM+2], (double)T0[i], (double)Ws[i]);
  }

  auto L2_norm = [&scaling](const Vector<Real>& V, const Vector<Real>& W) {
    const Long N = W.Dim();
    const Long dof = V.Dim()/W.Dim();
    Real sum = 0;
    for (Long i = 0; i < N; i++) {
      for (Long k = 0; k < dof; k++) {
        sum += V[i*dof+k] * V[i*dof+k] * W[i];
      }
    }
    return sqrt<Real>(sum);
  };
  auto Linf_norm = [&scaling](const Vector<Real>& V) {
    Real max_val = 0;
    for (const auto& x : V) max_val = std::max<Real>(max_val, fabs(x));
    return max_val;
  };
  if (1) { // Check Green's identity on the surface
    BIOp_StokesFxU.SetTargetCoord(Vector<Real>());
    BIOp_StokesDxU.SetTargetCoord(Vector<Real>());

    Vector<Real> DU, SF;
    BIOp_StokesFxU.ComputePotential(SF, F);
    BIOp_StokesDxU.ComputePotential(DU, U); DU -= 0.5*U;
    Vector<Real> E = SF + DU;

    std::cout<<"Green's identity L2-error (on-surface) = "<<L2_norm(E,Wa_surf)/L2_norm(U,Wa_surf)<<" (same as GMRES residual)\n";
    std::cout<<"Green's identity Linf-error (on-surface) = "<<Linf_norm(E)/Linf_norm(U)<<'\n';
    elem_lst0.WriteVTK("vis/EE-helix", E/Linf_norm(U), comm);
  }
  if (1) { // Check Green's identity on the centerline
    BIOp_StokesFxU.SetTargetCoord(Xc);
    BIOp_StokesDxU.SetTargetCoord(Xc);

    Vector<Real> DU, SF;
    BIOp_StokesFxU.ComputePotential(SF, F);
    BIOp_StokesDxU.ComputePotential(DU, U);
    Vector<Real> E = SF + DU;

    std::cout<<"Green's identity L2-error (off-surface) = "<<L2_norm(E,Ws)/L2_norm(U,Ws)<<'\n';
    std::cout<<"Green's identity Linf-error (off-surface) = "<<Linf_norm(E)/Linf_norm(U)<<'\n';
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


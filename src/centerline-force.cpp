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
template <class Real> void test(const Comm& comm, Real r0) {
  static constexpr Integer COORD_DIM = 3;
  const Real quad_tol = 1e-14;
  const Real gmres_tol = 1e-13;
  const Long gmres_max_iter = 400;
  const Real pi = const_pi<Real>();

  const Real scale = 1 / (sqrt<Real>(2.0)*pi);
  const Real end_len = r0/scale/sqrt<Real>(8.0); // length of end-caps (in parameter space)
  const Real s0 = 0 - end_len; // range of s
  const Real s1 = 1 + end_len;

  // Geometry function for half-turn helix
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

  Vector<Real> panel_len; // define length of panels (in parameter s)
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
  const Integer Npanels = panel_len.Dim();

  SlenderElemList<Real> elem_lst0; // Initialize geometry from geom function
  GenericGeom(elem_lst0, geom, panel_len, Vector<Long>(), comm, ChebOrder, FourierOrder);

  Vector<Real> U(Npanels * ChebOrder * FourierOrder * COORD_DIM);
  Vector<Real> Wa(Npanels * ChebOrder * FourierOrder); // weights to average force on surface to the centerline
  Vector<Real> PanelNodes = SlenderElemList<Real>::CenterlineNodes(ChebOrder); // Chebyshev discretization nodes
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
        for (Long j = 0; j < FourierOrder; j++) {
          const Long node_idx = i * FourierOrder + j;
          const auto& Xs = dX_ds.begin() + node_idx * COORD_DIM;
          const auto& Xt = dX_dt.begin() + node_idx * COORD_DIM;
          Wa[panel_idx * ChebOrder*FourierOrder + node_idx] = sqrt<Real>(Xs[0]*Xs[0]+Xs[1]*Xs[1]+Xs[2]*Xs[2]) * sqrt<Real>(Xt[0]*Xt[0]+Xt[1]*Xt[1]+Xt[2]*Xt[2]) / sqrt<Real>(X0s[0*ChebOrder+i]*X0s[0*ChebOrder+i] + X0s[1*ChebOrder+i]*X0s[1*ChebOrder+i] + X0s[2*ChebOrder+i]*X0s[2*ChebOrder+i]);
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
  Vector<Real> F0(Npanels * ChebOrder * COORD_DIM), s_vec(Npanels * ChebOrder);
  for (Long panel_idx = 0; panel_idx < Npanels; panel_idx++) { // loop over panels
    for (Long i = 0; i < ChebOrder; i++) {
      s_vec[panel_idx*ChebOrder+i] = s_offset + PanelNodes[i]*panel_len[panel_idx];

      for (Long k = 0; k < COORD_DIM; k++) {
        Real F0_ = 0;
        for (Long j = 0; j < FourierOrder; j++) {
          F0_ += F[((panel_idx * ChebOrder + i) * FourierOrder + j) * COORD_DIM + k] * Wa[(panel_idx * ChebOrder + i) * FourierOrder + j];
        }
        F0[(panel_idx * ChebOrder + i) * COORD_DIM + k] = F0_;
      }
    }
    s_offset += panel_len[panel_idx];
  }

  // print result
  printf("                   s                   Fx                   Fy                   Fz\n");
  for (Long i = 0; i < s_vec.Dim(); i++) {
    printf("%20.12e %20.12e %20.12e %20.12e\n", s_vec[i], F0[i*COORD_DIM+0], F0[i*COORD_DIM+1], F0[i*COORD_DIM+2]);
  }
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    Profile::Enable(true);
    Comm comm = Comm::Self();
    commandline_option_start(argc, argv, nullptr, comm);
    double eps = strtod(commandline_option(argc, argv, "-eps", "1e-2" , true, "the cross-sectional redius of the fiber", comm), nullptr);
    commandline_option_end(argc, argv);

    test<double>(comm, eps);
  }

  Comm::MPI_Finalize();
  return 0;
}


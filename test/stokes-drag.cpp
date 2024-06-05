#include <csbq.hpp>
using namespace sctl;

template <class Real> void test(const Comm& comm, const std::string fname, const Real Ux, const Real Uy, const Real Uz, const Real gmres_tol, const Real quad_tol, const Real Rmaj, const Real Rmin, const Real eps, const Real SL_scal) {
  static constexpr Integer COORD_DIM = 3;

  SlenderElemList<Real> elem_lst0;
  if (!fname.empty()) {
    elem_lst0.template Read<QuadReal>(fname, comm);
  }
  if (!elem_lst0.Size()) { // use default geom
    const Long Npanel = 40;
    Vector<Real> panel_len(Npanel); panel_len = 1/(Real)Npanel;
    Vector<Long> FourierOrder(Npanel); FourierOrder = 24;
    GeomEllipse(elem_lst0, panel_len, FourierOrder, comm, Rmaj, Rmin, (Real)eps*2, 10);
  }

  Long N = 0;
  { // Set N
    Vector<Real> X;
    elem_lst0.GetNodeCoord(&X, nullptr, nullptr);
    N = X.Dim()/COORD_DIM;
  }

  Vector<Real> U(N * COORD_DIM);
  for (Long i = 0; i < N; i++) {
    U[i*COORD_DIM+0] = Ux;
    U[i*COORD_DIM+1] = Uy;
    U[i*COORD_DIM+2] = Uz;
  }

  const Stokes3D_FxU ker_FxU;
  const Stokes3D_DxU ker_DxU;
  const Stokes3D_FSxU ker_FSxU;
  BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesFxU(ker_FxU, false, comm);
  BoundaryIntegralOp<Real,Stokes3D_DxU> BIOp_StokesDxU(ker_DxU, false, comm);
  BIOp_StokesFxU.SetFMMKer(ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU);
  BIOp_StokesDxU.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
  BIOp_StokesFxU.AddElemList(elem_lst0);
  BIOp_StokesDxU.AddElemList(elem_lst0);
  BIOp_StokesFxU.SetAccuracy(quad_tol);
  BIOp_StokesDxU.SetAccuracy(quad_tol);

  auto BIOp = [&BIOp_StokesFxU,&BIOp_StokesDxU,&SL_scal](Vector<Real>* Ax, const Vector<Real>& x) {
    Vector<Real> Udl, Usl;
    BIOp_StokesFxU.ComputePotential(Usl, x);
    BIOp_StokesDxU.ComputePotential(Udl, x);
    (*Ax) = (x*0.5 + Udl) + Usl*SL_scal;
    //(*Ax) = Usl*SL_scal;
  };

  Vector<Real> q;
  GMRES<Real> solver(comm);
  solver(&q, BIOp, U, gmres_tol, 100);

  auto SurfIntegral = [&comm](const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst) {
    Vector<Long> cnt, dsp;
    Vector<Real> X, Xn, wts, dist_far, F;
    elem_lst.GetFarFieldNodes(X, Xn, wts, dist_far, cnt, (Real)1);
    elem_lst.GetFarFieldDensity(F, F_cheb);
    comm.PartitionN(F, (comm.Rank()==0?1:0));
    comm.PartitionN(wts, (comm.Rank()==0?1:0));
    comm.PartitionN(cnt, (comm.Rank()==0?1:0));
    if (!wts.Dim()) return Vector<Real>();

    dsp.ReInit(cnt.Dim()); dsp = 0;
    omp_par::scan(cnt.begin(), dsp.begin(), cnt.Dim());

    const Integer dof = F.Dim() / wts.Dim();
    Vector<Real> IntegralF(dof); IntegralF = 0;
    for (Long elem = 0; elem < cnt.Dim(); elem++) {
      for (Long i = 0; i < cnt[elem]; i++) {
        for(Long k = 0; k < dof; k++) {
          IntegralF[k] += F[(dsp[elem]+i)*dof+k] * wts[dsp[elem]+i];
        }
      }
    }
    return IntegralF;
  };
  Vector<Real> F = SurfIntegral(q*SL_scal, elem_lst0);
  if (!comm.Rank()) {
    std::cout<<std::setprecision(14);
    std::cout<<"Fx = "<<F[0]<<'\n';
    std::cout<<"Fy = "<<F[1]<<'\n';
    std::cout<<"Fz = "<<F[2]<<'\n';
  }
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    #ifdef SCTL_HAVE_PVFMM
    pvfmm::Profile::Enable(true);
    #endif
    using Real = double;

    Profile::Enable(true);
    Comm comm = Comm::World();
    commandline_option_start(argc, argv, "\
Solve the exterior Stokes flow for a rigid ring with a given\n\
translation velocity (Ux,Uy,Uz) and return the total drag on the ring.\n", comm);
    std::string fname = commandline_option(argc, argv, "-fname", "", false, nullptr, comm);
    Real gmres_tol = strtod(commandline_option(argc, argv, "-gmres_tol", "1e-20", false, nullptr, comm), nullptr);
    Real quad_tol = strtod(commandline_option(argc, argv, "-quad_tol", "1e-20", false, nullptr, comm), nullptr);
    Real Ux = strtod(commandline_option(argc, argv, "-Ux", "0", false, nullptr, comm), nullptr);
    Real Uy = strtod(commandline_option(argc, argv, "-Uy", "0", false, nullptr, comm), nullptr);
    Real Uz = strtod(commandline_option(argc, argv, "-Uz", "1", false, nullptr, comm), nullptr);
    Real Rmaj = strtod(commandline_option(argc, argv, "-Rmaj", "1.0", false, nullptr, comm), nullptr);
    Real Rmin = strtod(commandline_option(argc, argv, "-Rmin", "1.0", false, nullptr, comm), nullptr);
    Real eps = strtod(commandline_option(argc, argv, "-eps", "1e-3", false, nullptr, comm), nullptr);
    Long omp_p   = strtol(commandline_option(argc, argv, "-omp"   , "1"    , false, nullptr, comm), nullptr, 10);
    commandline_option_end(argc, argv);
    omp_set_num_threads(omp_p);

    if (!comm.Rank()) {
      std::cout<<"Command:\n";
      for (Integer i = 0; i < argc; i++) {
        std::cout<<argv[i]<<' ';
      }
      std::cout<<'\n';
    }

    test<Real>(comm, fname, Ux, Uy, Uz, gmres_tol, quad_tol, Rmaj, Rmin, eps, 1.0/eps);
  }

  Comm::MPI_Finalize();
  return 0;
}


#include <csbq.hpp>
using namespace sctl;

enum TestCase {
  AxisAligned = 1,
  Twisted = 2
};

template <class Real> void test(const Long test_case, const Real eps, const Real delta) {
  constexpr Integer COORD_DIM = 3;
  const Comm comm = Comm::Self();
  const Real geom_tol = 1e-8;
  const Real quad_tol = 1e-12;
  const Real gmres_tol = 1e-10;

  const Stokes3D_FxU ker_FxU;
  const Stokes3D_DxU ker_DxU;
  const Stokes3D_FxUP ker_FxUP;
  const Stokes3D_FSxU ker_FSxU;
  BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesFxU(ker_FxU, false, comm);
  BoundaryIntegralOp<Real,Stokes3D_DxU> BIOp_StokesDxU(ker_DxU, false, comm);
  BIOp_StokesDxU.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
  BIOp_StokesFxU.SetFMMKer(ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU);
  BIOp_StokesFxU.SetAccuracy(quad_tol);
  BIOp_StokesDxU.SetAccuracy(quad_tol);

  const Long Nobj = 2;
  const Real SL_scal = 1/eps;
  const Long DefaultChebOrder = 10;
  const Long DefaultFourierOrder = 48;
  const Long DefaultElemCount = 128;

  RigidBodyList<Real> geom(comm);
  { // Init geom
    Vector<Long> obj_elem_cnt(Nobj), ChebOrder(Nobj*DefaultElemCount), FourierOrder(Nobj*DefaultElemCount);
    obj_elem_cnt = DefaultElemCount;
    ChebOrder = DefaultChebOrder;
    FourierOrder = DefaultFourierOrder;

    Vector<Real> Mr_lst, panel_len(Nobj * DefaultElemCount);
    RigidBodyList<Real>::RotationMatrix(Mr_lst, Vector<Real>(Nobj*COORD_DIM), 0);
    panel_len = Nobj / (Real)panel_len.Dim();

    Vector<Real> R(Nobj * DefaultElemCount * DefaultChebOrder); R = eps;
    Vector<Real> X(Nobj * DefaultElemCount * DefaultChebOrder * COORD_DIM);
    Vector<Real> OrientVec(Nobj * DefaultElemCount * DefaultChebOrder * COORD_DIM);
    const Vector<Real>& CenterlineNds = SlenderElemList<Real>::CenterlineNodes(DefaultChebOrder);
    for (Long i = 0; i < Nobj; i++) { // Set R, X, OrientVec
      for (Long j = 0; j < DefaultElemCount; j++) {
        for (Long k = 0; k < DefaultChebOrder; k++) {
          const Long node_idx = (i * DefaultElemCount + j) * DefaultChebOrder + k;
          const Real theta =  2 * const_pi<Real>() * (j + CenterlineNds[k]) / DefaultElemCount;

          const Real x_dsp = i * (2*(1+eps) + delta);
          const Real omega = i * (test_case == TestCase::Twisted ? const_pi<Real>()/3 : 0);

          X[node_idx*COORD_DIM+0] = cos<Real>(theta) + x_dsp;
          X[node_idx*COORD_DIM+1] = cos<Real>(omega) * sin<Real>(theta);
          X[node_idx*COORD_DIM+2] = sin<Real>(omega) * sin<Real>(theta);

          OrientVec[node_idx*COORD_DIM+0] = 0;
          OrientVec[node_idx*COORD_DIM+1] = -sin<Real>(omega);
          OrientVec[node_idx*COORD_DIM+2] =  cos<Real>(omega);
        }
      }
    }

    geom.Init(obj_elem_cnt, Mr_lst, ChebOrder, FourierOrder, panel_len, X, R, OrientVec);
  }

  Vector<Real> F;
  Long Nunknown = 0, Nunknown_ = 0;
  for (Long iter = 0; iter < 20; iter++) { // Solve Dirichlet BVP
    if (iter) { // Adaptive refinement
      Vector<Vector<Real>> F_(1); F_[0] = F;
      geom.RefineAdaptive(F_, geom_tol);
    }

    const SlenderElemList<Real>& elem_lst0 = geom.GetElemList();
    BIOp_StokesFxU.AddElemList(elem_lst0);
    BIOp_StokesDxU.AddElemList(elem_lst0);
    const Long N = BIOp_StokesFxU.Dim(1);

    Real DL_scal = 1.0;
    auto BIOp = [&BIOp_StokesFxU,&BIOp_StokesDxU,&DL_scal,&SL_scal](Vector<Real>* Ax, const Vector<Real>& x) {
      Vector<Real> Udl, Usl;
      BIOp_StokesFxU.ComputePotential(Usl, x);
      BIOp_StokesDxU.ComputePotential(Udl, x);
      (*Ax) = (x*0.5 + Udl)*DL_scal + Usl*SL_scal;
    };

    Vector<Real> U0(N);
    { // Set U0
      const auto& obj_elem_cnt = geom.obj_elem_cnt;
      const auto& obj_elem_dsp = geom.obj_elem_dsp;
      const auto& FourierOrder = geom.FourierOrder;
      const auto& ChebOrder = geom.ChebOrder;

      U0 = 0;
      Long node_dsp = 0;
      for (Long i = 0; i < Nobj; i++) {
        for (Long j = 0; j < obj_elem_cnt[i]; j++) {
          Long elem_idx = obj_elem_dsp[i] + j;
          Long Nnds = ChebOrder[elem_idx] * FourierOrder[elem_idx];
          for (Long k = 0; k < Nnds; k++) {
            if (test_case == TestCase::AxisAligned) {
              U0[(node_dsp + k)*COORD_DIM + 0] = 0.8*i;
              U0[(node_dsp + k)*COORD_DIM + 1] =-0.5*i;
              U0[(node_dsp + k)*COORD_DIM + 2] = 1.0*i;
            } else if (test_case == TestCase::Twisted) {
              U0[(node_dsp + k)*COORD_DIM + 0] = 0.8*i;
              U0[(node_dsp + k)*COORD_DIM + 1] =-0.5*i;
              U0[(node_dsp + k)*COORD_DIM + 2] = 1.0*i;
            } else {
              SCTL_ASSERT_MSG(false, "Test case not implemented");
            }
          }
          node_dsp += Nnds;
        }
      }
    }
    geom.GetElemList().WriteVTK("vis/U", U0, comm);

    F.ReInit(0);
    Long Ngmres;
    GMRES<Real> solver(comm);
    solver(&F, BIOp, U0, gmres_tol, 200, false, &Ngmres);
    geom.GetElemList().WriteVTK("vis/F", F, comm);

    BIOp_StokesFxU.template DeleteElemList<SlenderElemList<Real>>();
    BIOp_StokesDxU.template DeleteElemList<SlenderElemList<Real>>();

    std::cout<<"Nunknown = "<<N<<'\n';
    if (Nunknown == N || Nunknown_ == N) break;
    Nunknown_ = Nunknown;
    Nunknown = N;
  }

  Vector<Real> Fc;
  Vector<QuadReal> Fc_;
  RigidBodyList<Real>::ObjIntegral(Fc, F, geom.GetElemList(), geom.obj_elem_cnt, geom.obj_elem_dsp, comm);
  for (const auto& a : Fc) Fc_.PushBack(a*SL_scal);
  std::cout<<"Force: "<<eps<<' '<<delta<<' '<<Vector<QuadReal>(3,Fc_.begin());
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    #ifdef SCTL_HAVE_PVFMM
    pvfmm::Profile::Enable(true);
    #endif
    Profile::Enable(true);
    Comm comm = Comm::World();
    commandline_option_start(argc, argv, "\
Solve the Stokes BVP for two rigid rings each with cross-sectional radius eps\n\
and deparated by a distance delta. Reports the total drag on the first ring\n\
when the second ring is translated with a velocity (0.8,-0.5,1.0)\n", comm);
    double eps   = strtod(commandline_option(argc, argv, "-eps"   , "1e-2" , false, nullptr, comm), nullptr);
    double delta = strtod(commandline_option(argc, argv, "-delta" , "1e-2" , false, nullptr, comm), nullptr);
    Long omp_p   = strtol(commandline_option(argc, argv, "-omp"   , "1"    , false, nullptr, comm), nullptr, 10);
    Long t_case  = strtol(commandline_option(argc, argv, "-test"  , "1"    , false, "1) Aligned tori\n\
                                   2) pi/3 twisted tori", comm), nullptr, 10);
    commandline_option_end(argc, argv);
    omp_set_num_threads(omp_p);

    if (!comm.Rank()) {
      std::cout<<"Command:\n";
      for (Integer i = 0; i < argc; i++) {
        std::cout<<argv[i]<<' ';
      }
      std::cout<<'\n';
    }

    test<double>(t_case, eps, delta);

    Profile::print(&comm);
  }

  Comm::MPI_Finalize();
  return 0;
}


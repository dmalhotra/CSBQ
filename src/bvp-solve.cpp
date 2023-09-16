#include <sctl.hpp>
#include "utils.hpp"
#include "mobility.hpp"

using namespace sctl;

using RefValueType = long double;

const std::string ref_path = "ref/";

constexpr bool UNIF_TANGLE = false;
const Long CubeResolution = 100;

template <Long ScalExp> struct Laplace3D_FDxU_Scal_ {
  static const std::string& Name() {
    static const std::string name = "Laplace3D-FDxU";
    //static const std::string name = ("Laplace3D-"+std::to_string(ScalExp) + "FDxU").c_str();
    return name;
  }
  static constexpr Integer FLOPS() {
    return 16;
  }
  template <class Real> static constexpr Real uKerScaleFactor() {
    return 1 / (4 * const_pi<Real>());
  }
  template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[1][1], const VecType (&r)[3], const VecType (&n)[3], const void* ctx_ptr) {
    constexpr auto SL_scal = pow<ScalExp>((typename VecType::ScalarType)2);
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv3 = rinv * rinv * rinv;
    VecType rdotn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
    u[0][0] = rinv * VecType(SL_scal);
    u[0][0] += rdotn * rinv3;
  }
};

template <Long ScalExp> struct Stokes3D_FDxU_Scal_ {
  static const std::string& Name() {
    static const std::string name = "Stokes3D-FDxU";
    //static const std::string name = ("Stokes3D-"+std::to_string(ScalExp) + "FDxU").c_str();
    return name;
  }
  static constexpr Integer FLOPS() {
    return 50;
  }
  template <class Real> static constexpr Real uKerScaleFactor() {
    return 1 / (8 * const_pi<Real>());
  }
  template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[3][3], const VecType (&r)[3], const VecType (&n)[3], const void* ctx_ptr) {
    constexpr auto SL_scal = pow<ScalExp>((typename VecType::ScalarType)2);
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv2 = rinv*rinv;
    VecType rinv3 = rinv2*rinv;
    VecType rinv5 = rinv3*rinv2;
    VecType rdotn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
    const VecType rdotn_rinv5_6 = VecType((typename VecType::ScalarType)6) * rdotn*rinv5;
    for (Integer i = 0; i < 3; i++) {
      for (Integer j = 0; j < 3; j++) {
        const VecType ri_rj = r[i]*r[j];
        u[i][j] = ((i==j ? rinv : VecType((typename VecType::ScalarType)0)) + ri_rj*rinv3) * VecType(SL_scal);
        u[i][j] += ri_rj * rdotn_rinv5_6;
      }
    }
  }
};

template <Long ScalExp> struct Laplace3D_FDxU_Scal : public GenericKernel<Laplace3D_FDxU_Scal_<ScalExp>> {};
template <Long ScalExp> struct Stokes3D_FDxU_Scal : public GenericKernel<Stokes3D_FDxU_Scal_<ScalExp>> {};

template <class Real, class KerSL, class KerDL, class KerSLDL, class KerM2M=KerSL, class KerM2L=KerSL, class KerM2T=KerSL, class KerL2L=KerSL, class KerL2T=KerSL> Real bvp_solve_conv(const std::string& fname, const std::string& ker_name, const Real tol, const Real gmres_tol, const Real SLScaling, const Comm& comm = Comm::World()) {
  SlenderElemList<RefValueType> elem_lst0;
  SCTL_ASSERT(!fname.empty());
  elem_lst0.template Read<RefValueType>(fname, comm);
  if (UNIF_TANGLE) { // Initialize elem_lst0 (uniform discretization of tangle)
    const Long Npanel = 8192;
    const Long ChebOrder_ = 10;
    const Long FourierOrder_ = 44; // 12

    Vector<Long> ChebOrder(Npanel); ChebOrder = ChebOrder_;
    Vector<Long> FourierOrder(Npanel); FourierOrder = FourierOrder_;
    Vector<RefValueType> panel_len(Npanel); panel_len = 1/(RefValueType)Npanel;
    GeomTangle<RefValueType>(elem_lst0, panel_len, FourierOrder, comm, ChebOrder_);
  }

  Vector<Real> Xt;
  CubeVolumeVis<Real> cube_vol(CubeResolution, 3.0);
  { // Set Xt
    if (!comm.Rank()) Xt = cube_vol.GetCoord();
    const Long N = cube_vol.GetCoord().Dim()/3;
    const Long i0 = N*(comm.Rank()+0)/comm.Size();
    const Long i1 = N*(comm.Rank()+1)/comm.Size();
    comm.PartitionN(Xt, i1-i0);
  }

  Vector<Real> Uref, Iref;
  { // Load reference solution from file
    std::string ref_fname = ref_path + std::string("Uref-") + ker_name;

    StaticArray<Long,2> Ntrg{Xt.Dim()/3*KerSL::TrgDim(), 0};
    comm.Allreduce(Ntrg+0, Ntrg+1, 1, Comm::CommOp::SUM);

    Vector<double> Uref0, Iref0;
    Uref0.Read(ref_fname.c_str());
    Iref0.Read((ref_fname+"id").c_str());
    bool compute_Uref = (Uref0.Dim()!=Ntrg[1]) || (Iref0.Dim()!=Ntrg[1]);

    if (compute_Uref) {
      Vector<RefValueType> Xtrg;
      for (const auto& a : Xt) Xtrg.PushBack((RefValueType)a);
      Vector<RefValueType> U0 = bvp_solve<RefValueType, KerSL, KerDL, KerM2M, KerM2L, KerM2T, KerL2L, KerL2T>(elem_lst0, (RefValueType)1e-17, (RefValueType)1e-15, (RefValueType)SLScaling, Vector<RefValueType>(), Xtrg, comm);
      comm.PartitionN(U0, (comm.Rank()?0:1));

      Uref0.ReInit(0);
      for (const auto& a : U0) Uref0.PushBack((double)a);
      if (!comm.Rank()) Uref0.Write(ref_fname.c_str());

      { // Compute indicator fuction using double-layer kernel
        KerDL ker;
        BoundaryIntegralOp<RefValueType,KerDL> BIOp(ker, false, comm);
        BIOp.AddElemList(elem_lst0);
        BIOp.SetAccuracy(1e-17);

        Vector<RefValueType> I0, F(BIOp.Dim(0)); F = (RefValueType)1;
        BIOp.ComputePotential(I0,F);
        comm.PartitionN(I0, (comm.Rank()?0:1));
        if (!comm.Rank()) { // Print error
          RefValueType err = 0;
          for (const auto& x : I0) err = std::max<RefValueType>(err, fabs(fabs(x)-(RefValueType)0.5)-(RefValueType)0.5);
          std::cout<<"Ref-quad error (on-surface) = "<<err<<'\n';
        }

        BIOp.SetTargetCoord(Xtrg);
        BIOp.ComputePotential(I0,F);
        comm.PartitionN(I0, (comm.Rank()?0:1));

        if (!comm.Rank()) { // Print error
          RefValueType err = 0;
          for (const auto& x : I0) err = std::max<RefValueType>(err, fabs(fabs(x)-(RefValueType)0.5)-(RefValueType)0.5);
          std::cout<<"Ref-quad error (off-surface) = "<<err<<'\n';
        }

        Iref0.ReInit(0);
        for (const auto& a : I0) Iref0.PushBack((double)a);
        if (!comm.Rank()) Iref0.Write((ref_fname+"id").c_str());
      }
    }
    if (!comm.Rank()) {
      for (const auto& a : Uref0) Uref.PushBack((Real)a);
      for (const auto& a : Iref0) Iref.PushBack((Real)a);
    }
  }

  SlenderElemList<Real> elem_lst;
  elem_lst0.Copy(elem_lst);

  #ifdef SCTL_HAVE_PVFMM
  Vector<Real> U = bvp_solve<Real, KerSL, KerDL, KerM2M, KerM2L, KerM2T, KerL2L, KerL2T>(elem_lst, tol, gmres_tol, SLScaling, Vector<Real>(), Xt, comm);
  #else
  Vector<Real> U = bvp_solve_combined<Real, KerSLDL>(elem_lst, tol, gmres_tol, Vector<Real>(), Xt, comm);
  #endif
  comm.PartitionN(U, (comm.Rank()?0:1));

  auto eval_dbl_layer = [&comm,elem_lst,tol](const Vector<Real>& Xt) {
    KerDL ker;
    KerM2M ker_m2m;
    KerM2L ker_m2l;
    KerM2T ker_m2t;
    KerL2L ker_l2l;
    KerL2T ker_l2t;
    BoundaryIntegralOp<Real,KerDL> BIOp(ker, false, comm);
    BIOp.SetFMMKer(ker, ker, ker, ker_m2m, ker_m2l, ker_m2t, ker_l2l, ker_l2t);
    BIOp.AddElemList(elem_lst);
    BIOp.SetTargetCoord(Xt);
    BIOp.SetAccuracy(tol);

    Vector<Real> I, F(BIOp.Dim(0)); F = (Real)1;
    BIOp.ComputePotential(I,F);
    comm.PartitionN(I, (comm.Rank()?0:1));
    return I;
  };
  Vector<Real> I0 = eval_dbl_layer(Vector<Real>());
  Vector<Real> I = eval_dbl_layer(Xt);

  Profile::print(&comm);
  Real max_err = 0;
  if (!comm.Rank() && U.Dim()) { // Print error
    Vector<Real> Uerr = U - Uref;
    for (Long i = 0; i < Uerr.Dim(); i++) {
      //if (fabs(Uerr[i]) > std::max(tol, gmres_tol)*10) {
      //  std::cout<<i<<' '<<Uerr[i]<<' '<<Iref[i]<<'\n';
      //}
      if (fabs(Iref[i]) > 0.5 || fabs(I[i]) > 0.5) { // for interior points
        Uerr[i] = 0;
        I[i] = 0;
      }
    }
    for (auto& x : I0) x += 0.5;
    if (1) { // print top 20 largest off-surface quadrature errors
      Vector<Real> II = I;
      for (auto& x : II) x = fabs(x);
      std::sort(II.begin(), II.end());
      for (Long j = 0; j < 20; j++) std::cout<<II[II.Dim() - j-1]<<'\n';
      for (Long i = 0; i < I.Dim(); i++) {
        if (fabs(I[i]) == II[II.Dim()-1]) {
          std::cout<<i<<'\n';
        }
      }
    }

    cube_vol.WriteVTK("vis/I-cube", I, Comm::Self());
    cube_vol.WriteVTK("vis/U-cube", U, Comm::Self());
    cube_vol.WriteVTK("vis/err-cube", Uerr, Comm::Self());
    { // Write quad-err-on-surf
      SlenderElemList<Real> elem_lst0;
      elem_lst0.template Read<RefValueType>(fname, Comm::Self());
      if (UNIF_TANGLE) { // Initialize elem_lst0 (uniform discretization of tangle)
        const Long Npanel = 8192;
        const Long ChebOrder_ = 10;
        const Long FourierOrder_ = 44; // 12

        Vector<Long> ChebOrder(Npanel); ChebOrder = ChebOrder_;
        Vector<Long> FourierOrder(Npanel); FourierOrder = FourierOrder_;
        Vector<RefValueType> panel_len(Npanel); panel_len = 1/(RefValueType)Npanel;
        GeomTangle<RefValueType>(elem_lst0, panel_len, FourierOrder, Comm::Self(), ChebOrder_);
      }
      elem_lst0.WriteVTK("vis/quad-err-on-surf", I0, Comm::Self());
    }

    for (auto x : Uerr) max_err = std::max<Real>(max_err, fabs(x));
    std::cout<<"Max-Error = "<<max_err<<'\n';

    Real l2_error = 0;
    for (auto x : Uerr) l2_error += x*x;
    l2_error = sqrt(l2_error/Uerr.Dim());
    std::cout<<"L2-Error = "<<l2_error<<'\n';

    Real quad_off_surf_err = 0, quad_on_surf_err = 0;
    for (auto x : I) quad_off_surf_err = std::max<Real>(quad_off_surf_err, fabs(x));
    for (auto x : I0) quad_on_surf_err = std::max<Real>(quad_on_surf_err, fabs(x));
    std::cout<<"Quad-Error (off-surf) = "<<quad_off_surf_err<<'\n';
    std::cout<<"Quad-Error (on-surf) = "<<quad_on_surf_err<<'\n';
    std::cout<<"Nunknown = "<<I0.Dim()<<'\n';
  }
  { // all-gather max_err
    const Long Np = comm.Size();
    Vector<Long> rcnt(Np), rdsp(Np);
    rcnt = 0; rcnt[0] = 1;
    rdsp = 1; rdsp[0] = 0;
    StaticArray<Real,2> err{max_err,0};
    comm.Allgatherv(err+0, (comm.Rank()==0?1:0), err+1, rcnt.begin(), rdsp.begin());
    max_err = err[1];
  }
  return max_err;
}

// Deprecated: use src/tangle-adap-geom.cpp
template <class Real> void tangle_adaptive_discretize(const Long Npanel, Real geom_tol) {
  Profile::Enable(true);
  const Comm comm = Comm::World();

  const std::string out_path("vis/");
  const Long ts_order = 1, start_idx = 0;
  const Real T = 1e-27, ts_tol = 0, gmres_tol = 1e-14, quad_tol = 1e-18;
  Real dt = 1e-28;

  RigidBodyList<Real> geom(comm);
  { // Set geom
    constexpr Integer COORD_DIM = 3;
    const Long ChebOrder_ = 10;
    const Long FourierOrder_ = 12;

    Vector<Long> ChebOrder(Npanel); ChebOrder = ChebOrder_;
    Vector<Long> FourierOrder(Npanel); FourierOrder = FourierOrder_;
    Vector<Real> panel_len(Npanel); panel_len = 1/(Real)Npanel;

    Vector<Real> X, R, OrientVec;
    { // Set X, R, OrientVec
      Vector<Real> X_surf, Xn;
      SlenderElemList<Real> elem_lst;
      GeomTangle(elem_lst, panel_len, FourierOrder, comm, ChebOrder_);
      elem_lst.GetNodeCoord(&X_surf, &Xn, nullptr);

      R.ReInit(Npanel*ChebOrder_);
      X.ReInit(Npanel*ChebOrder_*COORD_DIM);
      OrientVec.ReInit(Npanel*ChebOrder_*COORD_DIM);
      for (Long i = 0; i < Npanel*ChebOrder_; i++) {
        Real R2 = 0;
        for (Integer k = 0; k < COORD_DIM; k++) {
          Real sum = 0;
          for (Integer j = 0; j < FourierOrder_; j++) {
            sum += X_surf[(i*FourierOrder_+j)*COORD_DIM+k];
          }
          X[i*COORD_DIM+k] = sum/FourierOrder_;
          OrientVec[i*COORD_DIM+k] = X_surf[i*FourierOrder_*COORD_DIM+k] - X[i*COORD_DIM+k];
          R2 += OrientVec[i*COORD_DIM+k]*OrientVec[i*COORD_DIM+k];
        }
        R[i] = sctl::sqrt<Real>(R2);
      }
    }

    Vector<Real> Xc(COORD_DIM); Xc = 0;
    Vector<Real> Mr_lst(COORD_DIM*COORD_DIM); Mr_lst = 0;
    for (Integer i = 0; i < COORD_DIM; i++) Mr_lst[i*COORD_DIM+i] = 1;
    Vector<Long> obj_elem_cnt(1); obj_elem_cnt = Npanel;

    geom.Init(obj_elem_cnt, Xc, Mr_lst, ChebOrder, FourierOrder, panel_len, X, R, OrientVec);
  }

  const BgFlow<Real> bg_flow(comm);
  const Mobility<Real> stokes_mobility(comm, 1.0);
  stokes_mobility.AdaptiveTimeStep(geom, dt, T, bg_flow, ts_order, ts_tol, quad_tol, gmres_tol, geom_tol, out_path, start_idx);

  std::ostringstream geom_tol_str;
  geom_tol_str << Npanel << "_" << geom_tol;
  geom.GetElemList().Write(std::string("./tangle/tangle_")+geom_tol_str.str(), comm);

  Profile::print(&comm);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  using Real = double;

  #ifdef SCTL_HAVE_PVFMM
  pvfmm::Profile::Enable(true);
  #endif

  { // Solve BVP
    const Comm comm = Comm::World();
    commandline_option_start(argc, argv, "\
      Solve the exterior Laplace or Stokes Dirichlet BVP:\n\
      (1/2 + D + S * SLScaling) sigma = V0, where V0 = 1. \n\
      Default geometry is a circular torus. The solution on a cube is written to a VTU file.\n", comm);
    std::string fname = commandline_option(argc, argv, "-geom", "", false, "Input geometry filename", comm);
    std::string ker   = commandline_option(argc, argv, "-ker" , "Stokes", false, "Stokes/Laplace"         , comm);
    bool conv = to_bool(commandline_option(argc, argv, "-conv" , "false", false, "Run convergence tests", comm));

    double thickness = strtod(commandline_option(argc, argv, "-r"        , "1e-03", false, "Thickness of default geometry" , comm), nullptr);
    double tol       = strtod(commandline_option(argc, argv, "-tol"      , "1e-10", false, ""                              , comm), nullptr);
    double gmres_tol = strtod(commandline_option(argc, argv, "-gmres_tol", "1e-08", false, ""                              , comm), nullptr);
    double SLScaling = strtod(commandline_option(argc, argv, "-scale", "1e3", false, "Single-layer operator scaling factor", comm), nullptr);
    Long omp_p     = strtol(commandline_option(argc, argv, "-omp"      , "1"   , false, nullptr, comm), nullptr, 10);
    omp_set_num_threads(omp_p);
    commandline_option_end(argc, argv);

    if (conv == false) {
      SlenderElemList<Real> elem_lst0;
      CubeVolumeVis<Real> cube_vol(50, 1.5);
      if (fname.empty()) { // Initialize elem_lst0
        GeomEllipse(elem_lst0, Vector<Real>(), Vector<Long>(), comm, (Real)1, (Real)1, (Real)thickness, 10);
      } else {
        elem_lst0.template Read<RefValueType>(fname, comm);
      }
      if ((!strcmp(ker.c_str(), "Stokes"))) {
        if (!comm.Rank()) std::cout<<"\n\n## Solving : (1/2 + D + S * SLscaling) sigma = V0    (Stokes exterior Dirichlet BVP)\n";
        Vector<Real> U = bvp_solve<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>(elem_lst0, tol, gmres_tol, SLScaling, Vector<Real>(), cube_vol.GetCoord(), comm);
        cube_vol.WriteVTK("vis/BVP-cube", U, comm);
      } else if ((!strcmp(ker.c_str(), "Laplace"))) {
        if (!comm.Rank()) std::cout<<"\n\n## Solving : (1/2 + D + S * scaling) sigma = V0    (Laplace exterior Dirichlet BVP)\n";
        Vector<Real> U = bvp_solve<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxU, Laplace3D_FxU, Laplace3D_FxU, Laplace3D_FxU, Laplace3D_FxU>(elem_lst0, tol, gmres_tol, SLScaling, Vector<Real>(), cube_vol.GetCoord(), comm);
        cube_vol.WriteVTK("vis/BVP-cube", U, comm);
      } else {
        if (!comm.Rank()) std::cout<<"Unknown kernel "<<ker<<'\n';
      }

    } else { // convergence tests
      if (0) { // Tangle - Laplace
        constexpr Long ScalExp = 5;
        Profile::Enable(true);

        auto param_search = [&comm](const std::string& fname, const std::string& label) {
          Real err0 = bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>(fname, label, 1e-15, 1e-15, 100, comm);
          Integer qtol[2] = {0,15};
          Integer gtol[2] = {0,15};
          while (qtol[0] < qtol[1]) {
            if (!comm.Rank()) std::cout<<"qtol = "<<qtol[0]<<' '<<qtol[1]<<'\n';
            Integer q0 = (qtol[0]+qtol[1])/2;
            Real err = bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>(fname, label, pow<Real>(0.1,q0), 1e-15, 100, comm);
            if (err > 2*err0) {
              qtol[0] = q0+1;
            } else {
              qtol[1] = q0;
            }
          }
          while (gtol[0] < gtol[1]) {
            if (!comm.Rank()) std::cout<<"gtol = "<<gtol[0]<<' '<<gtol[1]<<'\n';
            Integer g0 = (gtol[0]+gtol[1])/2;
            Real err = bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>(fname, label, 1e-15, pow<Real>(0.1,g0), 100, comm);
            if (err > 2*err0) {
              gtol[0] = g0+1;
            } else {
              gtol[1] = g0;
            }
          }
          Real err = bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>(fname, label, pow<Real>(0.1,qtol[0]), pow<Real>(0.1,gtol[0]), 100, comm);
          if (!comm.Rank()) std::cout<<"Param = "<<err0<<' '<<err<<' '<<qtol[0]<<' '<<gtol[0]<<'\n';
        };

        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<-5>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<-5>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<-3>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<-3>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<-1>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<-1>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal< 1>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow< 1>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal< 3>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow< 3>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal< 4>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow< 4>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal< 5>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow< 5>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal< 6>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow< 6>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal< 7>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow< 7>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal< 9>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow< 9>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<11>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<11>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<13>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<13>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<15>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<15>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<17>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<17>((Real)2), comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<19>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, pow<19>((Real)2), comm);

        //const Real SL_scal = pow<ScalExp>((Real)2);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-12", "Laplace", 1e+01, 1e+01, SL_scal, comm); // compute reference solution
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-01", "Laplace", 1e-02, 1e-02, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-01", "Laplace", 1e-02, 1e-02, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-02", "Laplace", 1e-03, 1e-02, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-02", "Laplace", 1e-03, 1e-02, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-03", "Laplace", 1e-03, 1e-02, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-03", "Laplace", 1e-04, 1e-04, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-04", "Laplace", 1e-06, 1e-05, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-04", "Laplace", 1e-06, 1e-05, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-05", "Laplace", 1e-07, 1e-06, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-05", "Laplace", 1e-07, 1e-07, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-06", "Laplace", 1e-08, 1e-07, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-06", "Laplace", 1e-09, 1e-08, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-07", "Laplace", 1e-09, 1e-08, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-07", "Laplace", 1e-09, 1e-09, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-08", "Laplace", 1e-10, 1e-10, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-08", "Laplace", 1e-10, 1e-10, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-09", "Laplace", 1e-10, 1e-10, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-09", "Laplace", 1e-11, 1e-11, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-10", "Laplace", 1e-13, 1e-11, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-10", "Laplace", 1e-13, 1e-11, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-11", "Laplace", 1e-14, 1e-12, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-11", "Laplace", 1e-14, 1e-12, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_3e-12", "Laplace", 1e-14, 1e-13, SL_scal, comm);
        //bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU_Scal<ScalExp>>("tangle/tangle_1e-12", "Laplace", 1e-14, 1e-13, SL_scal, comm);

        //         geom   gmres_tol      tol       N     Nelem Max-FourierOrder        alpha   iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        // tangle_3e-07       1e-08    1e-09   23200       204               24      3.1e-02    200     5.6e-06        4.8e-08     1.4057   1.6504e+04   104.5540     0.1193     5.9556
        // tangle_3e-07       1e-08    1e-09   23200       204               24      1.3e-01    115     1.4e-06        1.3e-08     1.4104   1.6449e+04    60.0337     0.1194     3.3694
        // tangle_3e-07       1e-08    1e-09   23200       204               24      5.0e-01     63     3.5e-07        2.9e-09     1.4101   1.6453e+04    32.7515     0.1188     1.8219
        // tangle_3e-07       1e-08    1e-09   23200       204               24      2.0e+00     35     7.7e-08        1.4e-09     1.4084   1.6473e+04    18.1157     0.1188     1.0149
        // tangle_3e-07       1e-08    1e-09   23200       204               24      8.0e+00     21     2.6e-08        6.8e-10     1.4090   1.6466e+04    10.8710     0.1189     0.6074
        // tangle_3e-07       1e-08    1e-09   23200       204               24      1.6e+01     21     2.0e-08        6.1e-10     1.4068   1.6491e+04    10.8614     0.1190     0.6072
        // tangle_3e-07       1e-08    1e-09   23200       204               24      3.2e+01     22     1.2e-08        5.8e-10     1.4041   1.6523e+04    11.4197     0.1180     0.6310
        // tangle_3e-07       1e-08    1e-09   23200       204               24      6.4e+01     23     1.3e-08        8.6e-10     1.4120   1.6431e+04    11.9160     0.1187     0.6661
        // tangle_3e-07       1e-08    1e-09   23200       204               24      1.3e+02     26     1.5e-08        8.1e-10     1.4059   1.6502e+04    13.4956     0.1183     0.7486
        // tangle_3e-07       1e-08    1e-09   23200       204               24      5.1e+02     36     1.9e-08        6.0e-10     1.4112   1.6440e+04    18.6854     0.1183     1.0429
        // tangle_3e-07       1e-08    1e-09   23200       204               24      2.0e+03     46     1.3e-08        4.5e-10     1.4103   1.6450e+04    23.9419     0.1187     1.3306
        // tangle_3e-07       1e-08    1e-09   23200       204               24      8.2e+03     53     1.1e-08        4.3e-10     1.4066   1.6494e+04    27.4460     0.1184     1.5423
        // tangle_3e-07       1e-08    1e-09   23200       204               24      3.3e+04     55     1.4e-08        4.3e-10     1.4091   1.6464e+04    28.5808     0.1187     1.6057
        // tangle_3e-07       1e-08    1e-09   23200       204               24      1.3e+05     55     1.5e-08        4.6e-10     1.4074   1.6484e+04    28.5196     0.1184     1.8459
        // tangle_3e-07       1e-08    1e-09   23200       204               24      5.2e+05     56     1.0e-08        3.6e-10     1.4207   1.6330e+04    29.0958     0.1183     1.6181

        //         geom   gmres_tol      tol       N     Nelem Max-FourierOrder        alpha   iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        // tangle_3e-01       1e-02    1e-02     640        16                4      3.2e+01      4     8.3e-02        1.3e-02     0.0167   3.8323e+04     0.0014     0.0062     0.0059
        // tangle_1e-01       1e-02    1e-02     640        16                4      3.2e+01      4     8.3e-02        1.3e-02     0.0167   3.8323e+04     0.0013     0.0054     0.0054
        // tangle_3e-02       1e-02    1e-03     800        17                8      3.2e+01      4     2.9e-02        1.2e-03     0.0338   2.3669e+04     0.0021     0.0088     0.0057
        // tangle_1e-02       1e-02    1e-03    1240        26                8      3.2e+01      4     2.6e-02        1.8e-03     0.0435   2.8506e+04     0.0043     0.0123     0.0059
        // tangle_3e-03       1e-02    1e-03    3480        49                8      3.2e+01      4     1.9e-02        1.2e-03     0.0934   3.7259e+04     0.0291     0.0248     0.0087
        // tangle_1e-03       1e-04    1e-04    4920        68                8      3.2e+01      9     9.1e-04        2.6e-05     0.1423   3.4575e+04     0.1758     0.0208     0.0210
        // tangle_3e-04       1e-05    1e-06    6920        88               12      3.2e+01     13     4.7e-05        7.4e-07     0.2921   2.3691e+04     0.4957     0.0375     0.0399
        // tangle_1e-04       1e-05    1e-06    8640       103               12      3.2e+01     13     1.9e-05        7.8e-07     0.3332   2.5930e+04     0.7681     0.0291     0.0564
        // tangle_3e-05       1e-06    1e-07   10280       117               16      3.2e+01     16     1.6e-06        9.3e-08     0.4497   2.2860e+04     1.3400     0.0356     0.0892
        // tangle_1e-05       1e-07    1e-07   12440       136               16      3.2e+01     19     5.8e-07        1.2e-08     0.5248   2.3704e+04     2.3282     0.0448     0.1360
        // tangle_3e-06       1e-07    1e-08   15560       157               20      3.2e+01     19     1.5e-07        1.1e-08     0.8058   1.9310e+04     4.4228     0.0714     0.2748
        // tangle_1e-06       1e-08    1e-09   18480       177               20      3.2e+01     22     2.3e-08        5.8e-10     1.1678   1.5825e+04     7.2966     0.0934     0.3987
        // tangle_3e-07       1e-08    1e-09   23200       204               24      3.2e+01     22     1.2e-08        5.8e-10     1.4101   1.6453e+04    11.4274     0.1181     0.6327
        // tangle_1e-07       1e-09    1e-09   27520       227               24      3.2e+01     24     5.3e-09        1.5e-10     1.6038   1.7159e+04    17.4778     0.1287     0.8143
        // tangle_3e-08       1e-10    1e-10   32600       249               28      3.2e+01     27     7.4e-10        2.0e-11     2.2462   1.4513e+04    27.5027     0.1961     1.5310
        // tangle_1e-08       1e-10    1e-10   40040       289               32      3.2e+01     27     7.4e-10        1.8e-11     2.8276   1.4160e+04    41.4370     0.3246     2.4062
        // tangle_3e-09       1e-10    1e-10   48600       331               32      3.2e+01     27     5.5e-10        1.5e-11     3.4950   1.3906e+04    61.0337     0.3727     3.4072
        // tangle_1e-09       1e-11    1e-11   57680       372               36      3.2e+01     30     1.3e-10        1.2e-12     5.4166   1.0649e+04    95.9387     0.6075     5.3853
        // tangle_3e-10       1e-11    1e-13   74600       457               40      3.2e+01     30     4.2e-11        1.0e-12    11.1902   6.6670e+03   160.1744     1.7415    10.0256
        // tangle_1e-10       1e-11    1e-13   94120       557               44      3.2e+01     30     4.5e-11        1.1e-12    15.8894   5.9230e+03   254.8467     2.5093    15.7502
        // tangle_3e-11       1e-12    1e-14  122920       693               48      3.2e+01     33     2.5e-12        8.6e-14    26.6056   4.6200e+03   479.1154     4.2832    29.5689
        // tangle_1e-11       1e-12    1e-14  160400       893               48      3.2e+01     33     1.9e-12        8.2e-14    35.7700   4.4840e+03   813.7500     5.3888    47.4621
        // tangle_3e-12       1e-13    1e-14  219000      1215               52      3.2e+01     35     9.1e-13        3.2e-14    51.2236   4.2750e+03  1604.6252     5.4192    92.0567
        // tangle_1e-12       1e-13    1e-14  294080      1604               56      3.2e+01     35     6.9e-13        1.8e-14    72.4466   4.0590e+03  2890.0566     9.4950   163.4421
        // With FMM acceleration
        // tangle_3e-01       1e-02    1e-02     640        16                4      3.2e+01      4     5.6e-02        7.7e-03     0.0323   1.9814e+04     0.0459     0.0113     0.1123
        // tangle_1e-01       1e-02    1e-02     640        16                4      3.2e+01      4     5.6e-02        7.7e-03     0.0324   1.9753e+04     0.0479     0.0124     0.1175
        // tangle_3e-02       1e-02    1e-03     800        17                8      3.2e+01      4     2.9e-02        1.1e-03     0.0653   1.2251e+04     0.0505     0.0177     0.1944
        // tangle_1e-02       1e-02    1e-03    1240        26                8      3.2e+01      4     2.5e-02        1.5e-03     0.0842   1.4727e+04     0.0581     0.0244     0.1911
        // tangle_3e-03       1e-02    1e-03    3480        49                8      3.2e+01      4     1.9e-02        1.1e-03     0.1806   1.9269e+04     0.1404     0.0505     0.1940
        // tangle_1e-03       1e-04    1e-04    4920        68                8      3.2e+01      9     4.1e-04        3.0e-05     0.2833   1.7367e+04     0.4564     0.0424     0.2680
        // tangle_3e-04       1e-05    1e-06    6920        88               12      3.2e+01     13     4.7e-05        8.7e-07     0.5888   1.1753e+04     1.2523     0.0779     0.7780
        // tangle_1e-04       1e-05    1e-06    8640       103               12      3.2e+01     13     1.9e-05        9.1e-07     0.6728   1.2842e+04     1.8573     0.0599     0.8025
        // tangle_3e-05       1e-06    1e-07   10280       117               16      3.2e+01     16     1.6e-06        9.4e-08     0.9174   1.1206e+04     3.2892     0.0759     1.9701
        // tangle_1e-05       1e-07    1e-07   12440       136               16      3.2e+01     19     6.4e-07        1.1e-08     1.0624   1.1709e+04     5.6154     0.0944     2.5839
        // tangle_3e-06       1e-07    1e-08   15560       157               20      3.2e+01     19     1.4e-07        1.1e-08     1.6141   9.6400e+03     8.6487     0.1478     2.1987
        // tangle_1e-06       1e-08    1e-09   18480       177               20      3.2e+01     22     2.3e-08        5.6e-10     2.3625   7.8220e+03    14.3233     0.1941     4.1998
        // tangle_3e-07       1e-08    1e-09   23200       204               24      3.2e+01     22     1.3e-08        5.6e-10     2.8582   8.1170e+03    16.4175     0.2423     6.7847
        // tangle_1e-07       1e-09    1e-09   27520       227               24      3.2e+01     24     7.7e-09        1.5e-10     3.2441   8.4830e+03    16.5445     0.2622     9.0959
        // tangle_3e-08       1e-10    1e-10   32600       249               28      3.2e+01     27     1.0e-09        1.1e-11     4.4729   7.2880e+03    28.3380     0.3974     9.0879
        // tangle_1e-08       1e-10    1e-10   40040       289               32      3.2e+01     27     9.9e-10        1.0e-11     5.6220   7.1220e+03    38.6874     0.6469    11.0412
        // tangle_3e-09       1e-10    1e-10   48600       331               32      3.2e+01     27     9.8e-10        1.1e-11     6.9676   6.9750e+03    52.0068     0.7443    10.9798
        // tangle_1e-09       1e-11    1e-11   57680       372               36      3.2e+01     30     9.9e-11        1.6e-12    10.7911   5.3450e+03    74.1201     1.1990    21.7926
        // tangle_3e-10       1e-11    1e-13   74600       457               40      3.2e+01     30     3.4e-11        1.1e-12    22.1130   3.3740e+03   120.7445     3.4448    31.4061
        // tangle_1e-10       1e-11    1e-13   94120       557               44      3.2e+01     30     5.6e-11        1.1e-12    31.2978   3.0070e+03   148.0406     4.9899    40.3828
        // tangle_3e-11       1e-12    1e-14  122920       693               48      3.2e+01     33     3.5e-11        1.6e-13    52.9124   2.3230e+03   265.2149     8.5526    67.2871
        // tangle_1e-11       1e-12    1e-14  160400       893               48      3.2e+01     33     6.8e-11        1.9e-13    70.9052   2.2620e+03   353.6618    10.5377    96.8143
        // tangle_3e-12       1e-13    1e-14  219000      1215               52      3.2e+01     35     1.2e-10        2.8e-13   102.4772   2.1370e+03   487.8219    10.8642    97.1319
        // tangle_1e-12       1e-13    1e-14  294080      1604               56      3.2e+01     35     1.6e-10        4.2e-13   144.7232   2.0320e+03   555.6083    18.8379    95.6947
      }
      if (0) { // Tangle - Stokes
        constexpr Long ScalExp = 6;
        Profile::Enable(true);

        auto param_search = [&comm](const std::string& fname, const std::string& label) {
          Real err0 = bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>>(fname, label, 1e-15, 1e-15, 100, comm);
          Integer qtol[2] = {0,15};
          Integer gtol[2] = {0,15};
          while (qtol[0] < qtol[1]) {
            if (!comm.Rank()) std::cout<<"qtol = "<<qtol[0]<<' '<<qtol[1]<<'\n';
            Integer q0 = (qtol[0]+qtol[1])/2;
            Real err = bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>>(fname, label, pow<Real>(0.1,q0), 1e-15, 100, comm);
            if (err > 2*err0) {
              qtol[0] = q0+1;
            } else {
              qtol[1] = q0;
            }
          }
          while (gtol[0] < gtol[1]) {
            if (!comm.Rank()) std::cout<<"gtol = "<<gtol[0]<<' '<<gtol[1]<<'\n';
            Integer g0 = (gtol[0]+gtol[1])/2;
            Real err = bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>>(fname, label, 1e-15, pow<Real>(0.1,g0), 100, comm);
            if (err > 2*err0) {
              gtol[0] = g0+1;
            } else {
              gtol[1] = g0;
            }
          }
          Real err = bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>>(fname, label, pow<Real>(0.1,qtol[0]), pow<Real>(0.1,gtol[0]), 100, comm);
          if (!comm.Rank()) std::cout<<"Param = "<<err0<<' '<<err<<' '<<qtol[0]<<' '<<gtol[0]<<'\n';
        };

        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<-5>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)-5), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<-3>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)-3), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<-1>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)-1), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 1>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real) 1), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 3>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real) 3), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 4>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real) 4), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 5>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real) 5), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 6>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real) 6), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 7>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real) 7), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 9>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real) 9), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<11>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)11), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<13>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)13), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<15>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)15), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<17>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)17), comm);
        //bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<19>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, pow<ScalExp>((Real)19), comm);

        const Real SL_scal = pow<ScalExp>((Real)2);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-12", "Stokes", 1e+01, 1e+01, SL_scal, comm); // compute reference solution
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-01", "Stokes", 1e-03, 1e-02, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-01", "Stokes", 1e-03, 1e-02, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-02", "Stokes", 1e-03, 1e-02, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-02", "Stokes", 1e-03, 1e-02, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-03", "Stokes", 1e-03, 1e-02, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-03", "Stokes", 1e-04, 1e-03, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-04", "Stokes", 1e-05, 1e-04, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-04", "Stokes", 1e-05, 1e-05, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-05", "Stokes", 1e-06, 1e-05, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-05", "Stokes", 1e-07, 1e-06, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-06", "Stokes", 1e-07, 1e-07, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-06", "Stokes", 1e-07, 1e-07, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-07", "Stokes", 1e-08, 1e-07, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-07", "Stokes", 1e-09, 1e-08, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-08", "Stokes", 1e-09, 1e-08, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-08", "Stokes", 1e-10, 1e-09, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-09", "Stokes", 1e-10, 1e-10, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-09", "Stokes", 1e-11, 1e-10, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-10", "Stokes", 1e-11, 1e-10, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-10", "Stokes", 1e-12, 1e-10, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-11", "Stokes", 1e-13, 1e-11, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-11", "Stokes", 1e-13, 1e-11, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_3e-12", "Stokes", 1e-14, 1e-12, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("tangle/tangle_1e-12", "Stokes", 1e-14, 1e-12, SL_scal, comm);

        //         geom   gmres_tol      tol       N     Nelem Max-FourierOrder        alpha   iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        // tangle_1e-07       1e-08    1e-09   82560       227               24      3.1e-02    200     3.6e-04        2.6e-05     3.6332   2.2724e+04   419.0546     0.3170    20.1481
        // tangle_1e-07       1e-08    1e-09   82560       227               24      1.3e-01    200     1.6e-05        5.6e-07     3.6173   2.2824e+04   418.8950     0.3197    20.2051
        // tangle_1e-07       1e-08    1e-09   82560       227               24      5.0e-01    189     8.6e-07        1.0e-08     3.6260   2.2769e+04   395.6605     0.3155    19.0964
        // tangle_1e-07       1e-08    1e-09   82560       227               24      2.0e+00     95     2.2e-07        2.7e-09     3.6034   2.2912e+04   198.0879     0.3522     9.8020
        // tangle_1e-07       1e-08    1e-09   82560       227               24      8.0e+00     49     5.9e-08        1.1e-09     3.6001   2.2933e+04   101.7530     0.3165     4.7676
        // tangle_1e-07       1e-08    1e-09   82560       227               24      1.6e+01     38     4.5e-08        1.1e-09     3.6062   2.2894e+04    78.8655     0.3179     3.7653
        // tangle_1e-07       1e-08    1e-09   82560       227               24      3.2e+01     38     4.6e-08        1.2e-09     3.6175   2.2822e+04    78.8716     0.3174     3.7403
        // tangle_1e-07       1e-08    1e-09   82560       227               24      6.4e+01     38     4.5e-08        1.2e-09     3.6056   2.2898e+04    78.8662     0.3187     3.7590
        // tangle_1e-07       1e-08    1e-09   82560       227               24      1.3e+02     38     4.4e-08        1.3e-09     3.6042   2.2907e+04    78.8637     0.3174     3.6804
        // tangle_1e-07       1e-08    1e-09   82560       227               24      5.1e+02     46     4.3e-08        9.3e-10     3.6064   2.2893e+04    95.5167     0.3169     4.5353
        // tangle_1e-07       1e-08    1e-09   82560       227               24      2.0e+03     63     4.0e-08        6.6e-10     3.6116   2.2860e+04   130.9461     0.3174     6.1158
        // tangle_1e-07       1e-08    1e-09   82560       227               24      8.2e+03     94     4.3e-08        4.9e-10     3.6005   2.2930e+04   195.8292     0.3181     9.2162
        // tangle_1e-07       1e-08    1e-09   82560       227               24      3.3e+04    143     7.7e-08        4.3e-10     3.6220   2.2794e+04   298.9685     0.3165    14.5431
        // tangle_1e-07       1e-08    1e-09   82560       227               24      1.3e+05    200     1.8e-07        5.0e-10     3.6171   2.2825e+04   418.8483     0.3166    20.1264
        // tangle_1e-07       1e-08    1e-09   82560       227               24      5.2e+05    200     6.5e-07        2.0e-09     3.6103   2.2868e+04   418.9345     0.3182    20.1452

        //         geom   gmres_tol      tol       N     Nelem Max-FourierOrder        alpha   iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        // tangle_3e-01       1e-02    1e-03    1920        16                4      6.4e+01      6     7.6e-02        2.3e-03     0.0432   4.4444e+04     0.0070     0.0101     0.0097
        // tangle_1e-01       1e-02    1e-03    1920        16                4      6.4e+01      6     7.6e-02        2.3e-03     0.0430   4.4651e+04     0.0070     0.0096     0.0087
        // tangle_3e-02       1e-02    1e-03    2400        17                8      6.4e+01      6     4.5e-02        2.3e-03     0.0568   4.2254e+04     0.0107     0.0133     0.0096
        // tangle_1e-02       1e-02    1e-03    3720        26                8      6.4e+01      6     4.8e-02        2.4e-03     0.0755   4.9272e+04     0.0221     0.0192     0.0113
        // tangle_3e-03       1e-02    1e-03   10440        49                8      6.4e+01      5     3.5e-02        2.8e-03     0.1918   5.4432e+04     0.1302     0.0420     0.0168
        // tangle_1e-03       1e-03    1e-04   14760        68                8      6.4e+01     10     2.5e-03        1.0e-04     0.2976   4.9597e+04     0.6107     0.0330     0.0420
        // tangle_3e-04       1e-04    1e-05   20760        88               12      6.4e+01     16     4.1e-04        1.0e-05     0.4761   4.3604e+04     1.9071     0.0492     0.1050
        // tangle_1e-04       1e-05    1e-05   25920       103               12      6.4e+01     22     5.5e-05        1.1e-06     0.5720   4.5315e+04     4.0387     0.0466     0.2145
        // tangle_3e-05       1e-05    1e-06   30840       117               16      6.4e+01     22     1.9e-05        1.1e-06     0.8015   3.8478e+04     5.7045     0.0677     0.3249
        // tangle_1e-05       1e-06    1e-07   37320       136               16      6.4e+01     27     4.3e-06        1.4e-07     1.1256   3.3156e+04    10.2226     0.0925     0.5435
        // tangle_3e-06       1e-07    1e-07   46680       157               20      6.4e+01     33     6.6e-07        1.1e-08     1.4157   3.2973e+04    19.5177     0.1341     1.1621
        // tangle_1e-06       1e-07    1e-07   55440       177               20      6.4e+01     33     6.0e-07        1.1e-08     1.6850   3.2902e+04    27.3703     0.1485     1.4526
        // tangle_3e-07       1e-07    1e-08   69600       204               24      6.4e+01     33     2.0e-07        9.1e-09     2.5730   2.7050e+04    49.0301     0.2500     2.6560
        // tangle_1e-07       1e-08    1e-09   82560       227               24      6.4e+01     38     4.5e-08        1.2e-09     3.6229   2.2788e+04    78.9074     0.3237     3.6885
        // tangle_3e-08       1e-08    1e-09   97800       249               28      6.4e+01     38     3.3e-08        1.2e-09     4.4114   2.2170e+04   110.3651     0.4483     6.1473
        // tangle_1e-08       1e-09    1e-10  120120       289               32      6.4e+01     43     4.3e-09        1.5e-10     6.5818   1.8250e+04   189.4604     0.9223    11.3647
        // tangle_3e-09       1e-10    1e-10  145800       331               32      6.4e+01     49     1.4e-09        1.2e-11     8.2596   1.7652e+04   315.5162     1.1152    18.1918
        // tangle_1e-09       1e-10    1e-11  173040       372               36      6.4e+01     49     7.0e-10        9.8e-12    13.1195   1.3190e+04   447.3592     1.8374    25.7692
        // tangle_3e-10       1e-10    1e-11  223800       457               40      6.4e+01     49     2.9e-10        9.8e-12    21.9485   1.0197e+04   746.9663     4.4575    48.4942
        // tangle_1e-10       1e-10    1e-12  282360       557               44      6.4e+01     49     2.4e-10        9.4e-12    36.1811   7.8040e+03  1190.4534     7.1906    76.7934
        // tangle_3e-11       1e-11    1e-13  368760       693               48      6.4e+01     54     2.4e-11        1.0e-12    64.1282   5.7500e+03  2238.2410    12.7075   143.5433
        // tangle_1e-11       1e-11    1e-13  481200       893               48      6.4e+01     54     2.4e-11        9.7e-13    84.3626   5.7040e+03  3788.9482    15.1770   227.7466
        // tangle_3e-12       1e-12    1e-14  657000      1215               52      6.4e+01     59     2.0e-12        8.5e-14   136.6449   4.8080e+03  7694.4373    17.1662   447.8937
        // tangle_1e-12       1e-12    1e-14  882240      1604               56      6.4e+01     59     2.1e-12        8.7e-14   190.4861   4.6320e+03 13821.4466    30.3955   793.1211
        // With FMM acceleration
        // tangle_3e-01       1e-02    1e-03    1920        16                4      6.4e+01      6     7.7e-02        2.3e-03     0.0863   2.2248e+04     0.0643     0.0191     0.7956
        // tangle_1e-01       1e-02    1e-03    1920        16                4      6.4e+01      6     7.7e-02        2.3e-03     0.0869   2.2094e+04     0.0643     0.0189     0.8058
        // tangle_3e-02       1e-02    1e-03    2400        17                8      6.4e+01      6     4.6e-02        2.3e-03     0.1137   2.1108e+04     0.0737     0.0277     0.7744
        // tangle_1e-02       1e-02    1e-03    3720        26                8      6.4e+01      6     4.9e-02        2.4e-03     0.1510   2.4636e+04     0.0994     0.0383     0.8749
        // tangle_3e-03       1e-02    1e-03   10440        49                8      6.4e+01      5     3.5e-02        2.8e-03     0.3734   2.7959e+04     0.3354     0.0896     0.6026
        // tangle_1e-03       1e-03    1e-04   14760        68                8      6.4e+01     10     2.6e-03        2.3e-04     0.5729   2.5764e+04     1.1055     0.0690     1.1046
        // tangle_3e-04       1e-04    1e-05   20760        88               12      6.4e+01     16     4.1e-04        1.1e-05     0.9671   2.1466e+04     3.4146     0.1042     9.5836
        // tangle_1e-04       1e-05    1e-05   25920       103               12      6.4e+01     22     5.6e-05        4.5e-06     1.1861   2.1853e+04     7.0591     0.0964     8.4118
        // tangle_3e-05       1e-05    1e-06   30840       117               16      6.4e+01     22     4.3e-05        4.4e-06     1.5912   1.9382e+04    10.7050     0.1361    11.9869
        // tangle_1e-05       1e-06    1e-07   37320       136               16      6.4e+01     27     4.3e-06        1.9e-07     2.2478   1.6603e+04    20.2690     0.1858    27.7722
        // tangle_3e-06       1e-07    1e-07   46680       157               20      6.4e+01     33     4.0e-06        1.3e-07     2.8330   1.6477e+04    31.6939     0.2873    22.1203
        // tangle_1e-06       1e-07    1e-07   55440       177               20      6.4e+01     33     4.0e-06        1.3e-07     3.4003   1.6304e+04    35.8504     0.3061    50.9908
        // tangle_3e-07       1e-07    1e-08   69600       204               24      6.4e+01     33     4.0e-06        1.3e-07     5.1884   1.3415e+04    45.9656     0.4970    38.6054
        // tangle_1e-07       1e-08    1e-09   82560       227               24      6.4e+01     38     7.2e-08        4.3e-09     7.1397   1.1564e+04    88.1763     0.6484    57.4392
        // tangle_3e-08       1e-08    1e-09   97800       249               28      6.4e+01     38     7.2e-08        4.3e-09     9.0328   1.0827e+04   104.4971     0.9131    60.5884
        // tangle_1e-08       1e-09    1e-10  120120       289               32      6.4e+01     43     4.4e-09        1.6e-10    13.0224   9.2240e+03   197.8166     1.8392   102.1248
        // tangle_3e-09       1e-10    1e-10  145800       331               32      6.4e+01     49     1.6e-09        5.6e-11    16.9161   8.6190e+03   270.5505     2.1645   143.2406
        // tangle_1e-09       1e-10    1e-11  173040       372               36      6.4e+01     49     1.2e-09        5.5e-11    26.2770   6.5850e+03   311.3909     3.6172   117.3024
        // tangle_3e-10       1e-10    1e-11  223800       457               40      6.4e+01     49     1.2e-09        5.5e-11    44.1771   5.0660e+03   488.4402     8.7472   216.9929
        // tangle_1e-10       1e-10    1e-12  282360       557               44      6.4e+01     49     2.4e-10        9.5e-12    71.9050   3.9270e+03   737.1814    14.2891   335.9796
        // tangle_3e-11       1e-11    1e-13  368760       693               48      6.4e+01     54     1.7e-10        1.5e-12   126.9316   2.9050e+03  1137.6443    25.2768   441.2389
        // tangle_1e-11       1e-11    1e-13  481200       893               48      6.4e+01     54     1.4e-10        1.5e-12   167.5552   2.8720e+03  1818.0844    30.2350   387.2485
        // tangle_3e-12       1e-12    1e-14  657000      1215               52      6.4e+01     59     1.4e-10        2.7e-13   268.6389   2.4457e+03  3500.9450    39.6764   826.6070
        // tangle_1e-12       1e-12    1e-14  882240      1604               56      6.4e+01     59     2.2e-10        3.2e-13   379.5377   2.3245e+03  4455.5450    75.8798  1034.5002
      }
      if (1) { // Close-to-touching - Stokes
        constexpr Long ScalExp = 1;
        Profile::Enable(true);

        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<-3>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real)-3), comm); // GMREs-iter = 88, Err = 6.05006e-07
        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<-2>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real)-2), comm); // GMREs-iter = 65, Err = 3.11569e-07
        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<-1>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real)-1), comm); // GMREs-iter = 52, Err = 1.53957e-07
        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 0>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real) 0), comm); // GMREs-iter = 40, Err = 8.02818e-08
        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 1>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real) 1), comm); // GMREs-iter = 35, Err = 4.27814e-08
        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 2>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real) 2), comm); // GMREs-iter = 37, Err = 3.08961e-08
        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 3>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real) 3), comm); // GMREs-iter = 46, Err = 2.20968e-08
        //// bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal< 4>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, pow<ScalExp>((Real) 4), comm); // GMREs-iter = 58, Err = 3.60305e-08

        const Real SL_scal = pow<ScalExp>((Real)2);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-02, 1e-01, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-03, 1e-02, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-04, 1e-03, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-04, 1e-04, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-06, 1e-05, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-07, 1e-06, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-08, 1e-07, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-11, 1e-09, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-12, 1e-10, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-13, 1e-11, SL_scal, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU_Scal<ScalExp>, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FSxU, Stokes3D_FxU, Stokes3D_FxU>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-14, 1e-12, SL_scal, comm);

        //      N  gmres_tol      tol      iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        //  64560      1e-01    1e-02         2     1.3e-01        3.4e-02     5.8751        10989     2.0165     0.9595     0.2651
        //  64560      1e-02    1e-03         4     2.0e-02        3.1e-03     6.7928         9504     4.1026     1.0957     0.5388
        //  64560      1e-03    1e-04         8     1.4e-03        1.8e-04     8.2092         7864    10.0741     1.2504     1.3027
        //  64560      1e-04    1e-04        13     4.6e-04        1.8e-05     8.2341         7841    16.3690     1.2523     2.1177
        //  64560      1e-05    1e-06        21     1.5e-05        1.5e-06    12.2541         5268    27.9678     1.8013     3.6426
        //  64560      1e-06    1e-07        28     1.8e-06        1.4e-07    14.3414         4502    38.3228     2.0359     5.0246
        //  64560      1e-07    1e-08        31     3.1e-07        2.2e-08    16.4855         3916    48.4967     2.2385     6.2109
        //  64560      1e-08    1e-09        35     3.5e-08        2.1e-09    20.0756         3216    56.4711     2.5640     7.1991
        //  64560      1e-09    1e-11        40     1.8e-09        1.1e-10    27.5760         2341    69.2363     3.0647     8.6031
        //  64560      1e-10    1e-12        45     3.5e-10        1.7e-11    32.3967         1993    82.2548     3.3874     9.9720
        //  64560      1e-11    1e-13        49     2.8e-11        1.6e-12    39.5530         1632    95.5994     3.7386    11.5925
        //  64560      1e-12    1e-14        52     6.7e-12        1.9e-13    47.5180         1359   108.0743     4.0715    12.2181
      }
    }
  }

  Comm::MPI_Finalize();
  return 0;
}


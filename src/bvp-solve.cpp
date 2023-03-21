#include <sctl.hpp>
#include "utils.hpp"
#include "mobility.hpp"

using namespace sctl;

using RefValueType = long double;

struct Laplace3D_FDxU100_ {
  static const std::string& Name() {
    static const std::string name = "Laplace3D-FxU";
    return name;
  }
  static constexpr Integer FLOPS() {
    return 0;
  }
  template <class Real> static constexpr Real uKerScaleFactor() {
    return 1 / (4 * const_pi<Real>());
  }
  template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[1][1], const VecType (&r)[3], const VecType (&n)[3], const void* ctx_ptr) {
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv3 = rinv * rinv * rinv;
    VecType rdotn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
    u[0][0] = rinv * VecType((typename VecType::ScalarType)100);
    u[0][0] += rdotn * rinv3;
  }
};

struct Stokes3D_FDxU100_ {
  static const std::string& Name() {
    static const std::string name = "Stokes3D-FxU";
    return name;
  }
  static constexpr Integer FLOPS() {
    return 0;
  }
  template <class Real> static constexpr Real uKerScaleFactor() {
    return 1 / (8 * const_pi<Real>());
  }
  template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[3][3], const VecType (&r)[3], const VecType (&n)[3], const void* ctx_ptr) {
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv2 = rinv*rinv;
    VecType rinv3 = rinv2*rinv;
    VecType rinv5 = rinv3*rinv2;
    VecType rdotn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
    for (Integer i = 0; i < 3; i++) {
      for (Integer j = 0; j < 3; j++) {
        u[i][j] = ((i==j ? rinv : VecType((typename VecType::ScalarType)0)) + r[i]*r[j]*rinv3) * VecType((typename VecType::ScalarType)100);
        u[i][j] += VecType((typename VecType::ScalarType)6) * r[i]*r[j]*rdotn*rinv5;
      }
    }
  }
};

struct Stokes3D_FDxU10_ {
  static const std::string& Name() {
    static const std::string name = "Stokes3D-FxU";
    return name;
  }
  static constexpr Integer FLOPS() {
    return 0;
  }
  template <class Real> static constexpr Real uKerScaleFactor() {
    return 1 / (8 * const_pi<Real>());
  }
  template <Integer digits, class VecType> static void uKerMatrix(VecType (&u)[3][3], const VecType (&r)[3], const VecType (&n)[3], const void* ctx_ptr) {
    VecType r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    VecType rinv = approx_rsqrt<digits>(r2, r2 > VecType::Zero());
    VecType rinv2 = rinv*rinv;
    VecType rinv3 = rinv2*rinv;
    VecType rinv5 = rinv3*rinv2;
    VecType rdotn = r[0]*n[0] + r[1]*n[1] + r[2]*n[2];
    for (Integer i = 0; i < 3; i++) {
      for (Integer j = 0; j < 3; j++) {
        u[i][j] = ((i==j ? rinv : VecType((typename VecType::ScalarType)0)) + r[i]*r[j]*rinv3) * VecType((typename VecType::ScalarType)10);
        u[i][j] += VecType((typename VecType::ScalarType)6) * r[i]*r[j]*rdotn*rinv5;
      }
    }
  }
};

struct Laplace3D_FDxU100 : public GenericKernel<Laplace3D_FDxU100_> {};
struct Stokes3D_FDxU100 : public GenericKernel<Stokes3D_FDxU100_> {};
struct Stokes3D_FDxU10 : public GenericKernel<Stokes3D_FDxU10_> {};

template <class Real, class KerSL, class KerDL, class KerSLDL> void bvp_solve_conv(const std::string& fname, const std::string& ker_name, const Real tol, const Real gmres_tol, const Real SLScaling, const Comm& comm = Comm::World()) {
  SlenderElemList<Real> elem_lst0;
  SCTL_ASSERT(!fname.empty());
  elem_lst0.Read(fname, comm);

  Vector<Real> Xt;
  CubeVolumeVis<Real> cube_vol(100, 3.0);
  { // Set Xt
    if (!comm.Rank()) Xt = cube_vol.GetCoord();
    const Long N = cube_vol.GetCoord().Dim()/3;
    const Long i0 = N*(comm.Rank()+0)/comm.Size();
    const Long i1 = N*(comm.Rank()+1)/comm.Size();
    comm.PartitionN(Xt, i1-i0);
  }

  Vector<Real> Uref;
  { // Load reference solution from file
    std::string ref_fname = std::string("vis/Uref-") + ker_name;

    StaticArray<Long,2> Ntrg{Xt.Dim()/3*KerSL::TrgDim(), 0};
    comm.Allreduce(Ntrg+0, Ntrg+1, 1, Comm::CommOp::SUM);

    Vector<double> Uref0;
    Uref0.Read(ref_fname.c_str());
    bool compute_Uref = (Uref0.Dim()!=Ntrg[1]);

    if (compute_Uref) {
      SlenderElemList<RefValueType> elem_lst;
      elem_lst0.Copy(elem_lst);

      Vector<RefValueType> Xtrg;
      for (const auto& a : Xt) Xtrg.PushBack((RefValueType)a);
      Vector<RefValueType> U0 = bvp_solve<RefValueType, KerSL, KerDL>(elem_lst, (RefValueType)1e-17, (RefValueType)1e-15, (RefValueType)SLScaling, Vector<RefValueType>(), Xtrg, comm);
      comm.PartitionN(U0, (comm.Rank()?0:1));

      Uref0.ReInit(0);
      for (const auto& a : U0) Uref0.PushBack((double)a);
      if (!comm.Rank()) Uref0.Write(ref_fname.c_str());
    }
    if (!comm.Rank()) {
      for (const auto& a : Uref0) Uref.PushBack((Real)a);
    }
  }

  //Vector<Real> U = bvp_solve<Real, KerSL, KerDL>(elem_lst0, tol, gmres_tol, SLScaling, Vector<Real>(), Xt, comm);
  Vector<Real> U = bvp_solve_combined<Real, KerSLDL>(elem_lst0, tol, gmres_tol, Vector<Real>(), Xt, comm);
  comm.PartitionN(U, (comm.Rank()?0:1));
  Profile::print(&comm);

  if (!comm.Rank() && U.Dim()) { // Print error
    Vector<Real> Uerr = U - Uref;

    cube_vol.WriteVTK("vis/U-cube", U, Comm::Self());
    cube_vol.WriteVTK("vis/err-cube", Uerr, Comm::Self());

    Real max_err = 0;
    for (auto x : Uerr) max_err = std::max<Real>(max_err, fabs(x));
    std::cout<<"Max-Error = "<<max_err<<'\n';

    Real l2_error = 0;
    for (auto x : Uerr) l2_error += x*x;
    l2_error = sqrt(l2_error/Uerr.Dim());
    std::cout<<"L2-Error = "<<l2_error<<'\n';
  }
}

template <class Real> void tangle_adaptive_discretize(const Long Npanel, Real geom_tol) {
  #ifdef SCTL_HAVE_PVFMM
  pvfmm::Profile::Enable(true);
  #endif
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
    commandline_option_end(argc, argv);

    if (conv == false) {
      SlenderElemList<Real> elem_lst0;
      CubeVolumeVis<Real> cube_vol(50, 1.5);
      if (fname.empty()) { // Initialize elem_lst0
        GeomEllipse(elem_lst0, Vector<Real>(), Vector<Long>(), comm, (Real)1, (Real)1, thickness, 10);
      } else {
        elem_lst0.Read(fname, comm);
      }
      if ((!strcmp(ker.c_str(), "Stokes"))) {
        if (!comm.Rank()) std::cout<<"\n\n## Solving : (1/2 + D + S * SLscaling) sigma = V0    (Stokes exterior Dirichlet BVP)\n";
        Vector<Real> U = bvp_solve<Real, Stokes3D_FxU, Stokes3D_DxU>(elem_lst0, tol, gmres_tol, SLScaling, Vector<Real>(), cube_vol.GetCoord(), comm);
        cube_vol.WriteVTK("vis/BVP-cube", U, comm);
      } else if ((!strcmp(ker.c_str(), "Laplace"))) {
        if (!comm.Rank()) std::cout<<"\n\n## Solving : (1/2 + D + S * scaling) sigma = V0    (Laplace exterior Dirichlet BVP)\n";
        Vector<Real> U = bvp_solve<Real, Laplace3D_FxU, Laplace3D_DxU>(elem_lst0, tol, gmres_tol, SLScaling, Vector<Real>(), cube_vol.GetCoord(), comm);
        cube_vol.WriteVTK("vis/BVP-cube", U, comm);
      } else {
        if (!comm.Rank()) std::cout<<"Unknown kernel "<<ker<<'\n';
      }

    } else { // convergence tests
      //tangle_adaptive_discretize<long double>(32, 1e-4);
      //tangle_adaptive_discretize<long double>(32, 1e-6);
      //tangle_adaptive_discretize<long double>(64, 1e-8);
      //tangle_adaptive_discretize<long double>(64, 1e-10);
      //tangle_adaptive_discretize<long double>(64, 1e-12);
      //tangle_adaptive_discretize<long double>(64, 1e-14);
      //tangle_adaptive_discretize<long double>(128, 1e-14);

      if (1) { // Tangle - Laplace
        Profile::Enable(true);
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle254.geom", "Laplace", 1e-12, 1e-11, 100, comm); // Compute reference solution
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle50.geom" , "Laplace", 1e-03, 1e-02, 100, comm);
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle100.geom", "Laplace", 1e-04, 1e-03, 100, comm);
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle150.geom", "Laplace", 1e-05, 1e-04, 100, comm);
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle230.geom", "Laplace", 1e-06, 1e-05, 100, comm);
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle240.geom", "Laplace", 1e-08, 1e-07, 100, comm);
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle250.geom", "Laplace", 1e-10, 1e-09, 100, comm);
        bvp_solve_conv<Real, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FDxU100>("./tangle/tangle254.geom", "Laplace", 1e-12, 1e-11, 100, comm);

        // Tangle BVP - Laplace
        //      geom   gmres_tol      tol       N     Nelem     FourierOrder         iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        // tangle50        1e-02    1e-03    2800        70                4            4      4.2e-2         9.2e-4     0.1302        21505     0.0314     0.0200     0.0131
        // tangle100       1e-03    1e-04    4880       122                4            7      4.9e-3         6.9e-5     0.2338        20873     0.1617     0.0195     0.0272
        // tangle150       1e-04    1e-05   13760       172                8           10      1.0e-3         8.5e-6     0.7216        19069     1.8098     0.0514     0.0940
        // tangle230       1e-05    1e-06   30240       252               12           14      3.1e-5         8.1e-7     1.8162        16650    12.2452     0.0905     2.5270
        // tangle240       1e-07    1e-08   31440       262               12           20      2.4e-7         8.2e-9     2.4693        12732    18.9716     0.2125     4.2385
        // tangle250       1e-09    1e-10   65280       272               24           28      1.1e-9        4.5e-11     7.7427         8431   114.0527     0.3250     7.1356
        // tangle254       1e-11    1e-12   77280       276               28           35     6.6e-11        5.5e-13    11.7547         6574   200.0480     0.5391    10.6896
      }
      if (1) { // Tangle - Stokes
        Profile::Enable(true);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle254.geom", "Stokes", 1e-12, 1e-11, 100, comm); // Compute reference solution
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle50.geom" , "Stokes", 1e-03, 1e-02, 100, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle100.geom", "Stokes", 1e-04, 1e-03, 100, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle150.geom", "Stokes", 1e-05, 1e-04, 100, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle230.geom", "Stokes", 1e-06, 1e-05, 100, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle240.geom", "Stokes", 1e-08, 1e-07, 100, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle250.geom", "Stokes", 1e-10, 1e-09, 100, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU100>("./tangle/tangle254.geom", "Stokes", 1e-12, 1e-11, 100, comm);

        // Tangle BVP - Stokes
        //      geom   gmres_tol      tol       N     Nelem     FourierOrder         iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        // tangle50         1e-2     1e-3    2800        70                4            6     2.1e-01        2.2e-03     0.1856        15086     0.1589     0.0248     0.0234
        // tangle100        1e-3     1e-4    4880       122                4           10     1.9e-02        1.0e-04     0.3313        14730     0.7745     0.0243     0.0565
        // tangle150        1e-4     1e-5   13760       172                8           16     1.7e-02        1.6e-05     1.2295        11192     9.8059     0.0770     1.8448
        // tangle230        1e-5     1e-6   30240       252               12           21     1.7e-04        9.0e-07     3.3138         9125    61.2092     0.1975     5.2584
        // tangle240        1e-7     1e-8   31440       262               12           33     4.1e-06        7.6e-09     4.4355         7088   104.3853     0.2241     7.6990
        // tangle250        1e-9    1e-10   65280       272               24           43     1.4e-08        1.1e-10    17.7085         3686   586.0695     0.7960    22.9405
        // tangle254       1e-11    1e-12   77280       276               28           54     4.1e-09        6.9e-12    27.6771         2792  1034.2305     1.2298    38.8589
      }
      if (1) { // Close-to-touching - Stokes
        Profile::Enable(true);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-02, 1e-01, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-03, 1e-02, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-04, 1e-03, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-05, 1e-04, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-06, 1e-05, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-07, 1e-06, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-08, 1e-07, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-09, 1e-08, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-10, 1e-09, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-11, 1e-10, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-12, 1e-11, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-13, 1e-12, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-14, 1e-13, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-15, 1e-14, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-16, 1e-15, 10, comm);
        bvp_solve_conv<Real, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FDxU10>("./data/close-to-touching.geom", "Stokes-close-touching", 1e-17, 1e-16, 10, comm);

        //      N  gmres_tol      tol      iter    MaxError       L2-error    T_setup   setup-rate    T_solve    T_setup    T_solve
        //  21520      1e-01     1e-2         2     1.3e-01        3.2e-02     5.6700                  3.1944     0.8531     0.4806
        //  21520      1e-02     1e-3         4     2.1e-02        2.5e-03     8.1061                  6.5360     1.2818     1.3614
        //  21520      1e-03     1e-4         7     1.6e-02        3.1e-04    10.8099                 11.8118     1.7274     2.2869
        //  21520      1e-04     1e-5        13     9.3e-03        2.4e-05    13.6997                 22.5707     2.1291     4.8351
        //  21520      1e-05     1e-6        24     2.4e-03        3.7e-06    16.8026                 42.8992     2.5001     7.6538
        //  21520      1e-06     1e-7        34     3.4e-05        2.2e-07    19.9488                 62.5492     2.8044    10.8931
        //  21520      1e-07     1e-8        43     2.8e-06        1.4e-08    23.5213                 81.6355     3.3077    12.7662
        //  21520      1e-08     1e-9        49     2.6e-07        1.9e-09    27.4751                 96.2095     3.7236    14.7706
        //  21520      1e-09    1e-10        54     9.3e-08        5.5e-10    31.4113                109.2922     3.9118    16.2876
        //  21520      1e-10    1e-11        59     5.4e-08        2.3e-10    35.6971                122.8530     4.0588    19.2035
        //  21520      1e-11    1e-12        64     5.0e-09        2.2e-11    40.5914                137.0600     4.5563    20.2282
        //  21520      1e-12    1e-13        69     5.0e-10        2.5e-12    45.6508                152.2238     4.9972    22.3425
        //  21520      1e-13    1e-14        72     1.3e-10        1.5e-12    49.9494                162.6172     5.2653    23.2362
      }
    }
  }

  Comm::MPI_Finalize();
  return 0;
}


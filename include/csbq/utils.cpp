#include <random>

namespace sctl {

template <class Real, class Kernel> void double_layer_test(const SlenderElemList<Real>& elem_lst0, const Comm& comm, Real tol) { // Double-layer identity test
  const Long pid = comm.Rank();

  const Kernel ker_fn;
  BoundaryIntegralOp<Real,Kernel> BIOp(ker_fn, false, comm);
  BIOp.AddElemList(elem_lst0);
  BIOp.SetAccuracy(tol);

  // Warm-up run
  Vector<Real> F(BIOp.Dim(0)), U; F = 1;
  BIOp.ComputePotential(U,F);
  BIOp.ClearSetup();
  U = 0;

  Profile::Enable(true);
  Profile::Tic("Setup+Eval", &comm, true);
  BIOp.ComputePotential(U,F);
  Profile::Toc();

  Vector<Real> Uerr = U + 0.5;
  elem_lst0.WriteVTK("vis/Uerr_", Uerr, comm); // Write VTK
  { // Print error
    StaticArray<Real,2> max_err{0,0};
    for (auto x : Uerr) max_err[0] = std::max<Real>(max_err[0], fabs(x));
    comm.Allreduce(max_err+0, max_err+1, 1, Comm::CommOp::MAX);
    if (!pid) std::cout<<"Error = "<<max_err[1]<<'\n';
  }
  Profile::Enable(false);
  Profile::print(&comm);
}

template <class Real, class KerSL, class KerDL, class KerGrad> void test_greens_identity(const SlenderElemList<Real>& elem_lst0, const Comm& comm, const Real tol, const Vector<Real> X0) {
  static constexpr Integer COORD_DIM = 3;
  const Long pid = comm.Rank();

  const KerSL kernel_sl;
  const KerDL kernel_dl;
  const KerGrad kernel_grad;
  BoundaryIntegralOp<Real,KerSL> BIOpSL(kernel_sl, false, comm);
  BoundaryIntegralOp<Real,KerDL> BIOpDL(kernel_dl, false, comm);
  BIOpSL.AddElemList(elem_lst0);
  BIOpDL.AddElemList(elem_lst0);
  BIOpSL.SetAccuracy(tol);
  BIOpDL.SetAccuracy(tol);

  Vector<Real> X, Xn, Fs, Fd, Uref, Us, Ud;
  { // Set Fs, Fd, Uref
    elem_lst0.GetNodeCoord(&X, &Xn, nullptr);

    Vector<Real> Xn0{0,0,0}, F0(KerSL::SrcDim()), dU;
    for (auto& x : F0) x = drand48()-0.5;
    kernel_sl  .Eval(Uref, X, X0, Xn0, F0);
    kernel_grad.Eval(dU  , X, X0, Xn0, F0);

    Fd = Uref;
    { // Set Fs <-- -dot_prod(dU, Xn)
      constexpr Integer KDIM0 = KerSL::SrcDim();
      const Long N = X.Dim()/COORD_DIM;
      Fs.ReInit(N * KDIM0);
      for (Long i = 0; i < N; i++) {
        for (Integer j = 0; j < KDIM0; j++) {
          Real dU_dot_Xn = 0;
          for (Long k = 0; k < COORD_DIM; k++) {
            dU_dot_Xn += dU[(i*KDIM0+j)*COORD_DIM+k] * Xn[i*COORD_DIM+k];
          }
          Fs[i*KDIM0+j] = dU_dot_Xn;
        }
      }
    }
  }

  // Warm-up run
  BIOpSL.ComputePotential(Us,Fs);
  BIOpDL.ComputePotential(Ud,Fd);
  BIOpSL.ClearSetup();
  BIOpDL.ClearSetup();
  Us = 0; Ud = 0;

  comm.Barrier();
  sctl::Profile::Enable(true);
  Profile::Tic("Setup+Eval", &comm);
  BIOpSL.ComputePotential(Us,Fs);
  BIOpDL.ComputePotential(Ud,Fd);
  Profile::Toc();

  Vector<Real> Uerr = (Fd*0.5 + Us - Ud) - Uref;
  elem_lst0.WriteVTK("vis/Uerr0", Uerr, comm);
  { // Print error
    StaticArray<Real,2> max_err{0,0};
    StaticArray<Real,2> max_val{0,0};
    for (auto x : Uerr) max_err[0] = std::max<Real>(max_err[0], fabs(x));
    for (auto x : Uref) max_val[0] = std::max<Real>(max_val[0], fabs(x));
    comm.Allreduce(max_err+0, max_err+1, 1, Comm::CommOp::MAX);
    comm.Allreduce(max_val+0, max_val+1, 1, Comm::CommOp::MAX);
    if (!pid) std::cout<<"Error = "<<max_err[1]/max_val[1]<<'\n';
  }

  sctl::Profile::print(&comm);
  sctl::Profile::Enable(false);
}

template <class Real, class KerSL, class KerDL, class KerM2M, class KerM2L, class KerM2T, class KerL2L, class KerL2T> Vector<Real> bvp_solve(const SlenderElemList<Real>& elem_lst0, const Real tol, const Real gmres_tol, const Real SLScaling, Vector<Real> V0, const Vector<Real> Xt, const Comm& comm, Vector<Real>* sigma_) {
  const KerSL ker_sl;
  const KerDL ker_dl;
  const KerM2M ker_m2m;
  const KerM2L ker_m2l;
  const KerM2T ker_m2t;
  const KerL2L ker_l2l;
  const KerL2T ker_l2t;
  BoundaryIntegralOp<Real,KerSL> SLOp(ker_sl, false, comm);
  BoundaryIntegralOp<Real,KerDL> DLOp(ker_dl, false, comm);
  #ifdef SCTL_HAVE_PVFMM
  SLOp.SetFMMKer(ker_sl, ker_sl, ker_sl, ker_m2m, ker_m2l, ker_m2t, ker_l2l, ker_l2t);
  DLOp.SetFMMKer(ker_dl, ker_dl, ker_dl, ker_m2m, ker_m2l, ker_m2t, ker_l2l, ker_l2t);
  #endif
  SLOp.AddElemList(elem_lst0);
  DLOp.AddElemList(elem_lst0);
  SLOp.SetAccuracy(tol);
  DLOp.SetAccuracy(tol);

  { // Warmup
    bool prof_enable = Profile::Enable(false);
    Vector<Real> F(SLOp.Dim(0)), U; F = 1;
    SLOp.ComputePotential(U,F);
    DLOp.ComputePotential(U,F);
    SLOp.ClearSetup();
    DLOp.ClearSetup();
    Profile::Enable(prof_enable);
    comm.Barrier();
  }

  auto BIOp = [&SLOp,&DLOp,&SLScaling](Vector<Real>* Ax, const Vector<Real>& x){
    Vector<Real> Usl, Udl;
    SLOp.ComputePotential(Usl,x);
    DLOp.ComputePotential(Udl,x);
    (*Ax) = x*0.5 + Udl + Usl*(SLScaling);
  };

  const Long N = SLOp.Dim(0);
  Vector<Real> sigma;
  if (V0.Dim() != N) {
    V0.ReInit(N);
    V0 = 1;
  }
  Profile::Tic("Solve", &comm, true);
  ParallelSolver<Real> solver(comm);
  solver(&sigma, BIOp, V0, gmres_tol, 200);
  if (0) { // Write solution to file
    auto sigma_ = sigma;
    comm.PartitionN(sigma_, (comm.Rank()?0:1));
    if (!comm.Rank()) sigma_.Write((std::string("ref/sigma")+KerSL::Name()).c_str());
  }
  Profile::Toc();

  if (sigma_) (*sigma_) = sigma;
  elem_lst0.WriteVTK("vis/sigma", sigma, comm);
  Vector<Real> Utrg;
  { // Evaluate at Xt
    SLOp.SetTargetCoord(Xt);
    DLOp.SetTargetCoord(Xt);

    Vector<Real> Usl, Udl;
    SLOp.ComputePotential(Usl,sigma);
    DLOp.ComputePotential(Udl,sigma);
    Utrg = Udl + Usl*(SLScaling);
  }
  return Utrg;
}

template <class Real, class Ker> Vector<Real> bvp_solve_combined(const SlenderElemList<Real>& elem_lst0, const Real tol, const Real gmres_tol, Vector<Real> V0, const Vector<Real> Xt, const Comm& comm, Vector<Real>* sigma_) {
  const Ker ker;
  BoundaryIntegralOp<Real,Ker> KerOp(ker, false, comm);
  KerOp.AddElemList(elem_lst0);
  KerOp.SetAccuracy(tol);

  { // Warmup
    bool prof_enable = Profile::Enable(false);
    Vector<Real> F(KerOp.Dim(0)), U; F = 1;
    KerOp.ComputePotential(U,F);
    KerOp.ClearSetup();
    Profile::Enable(prof_enable);
  }

  auto BIOp = [&KerOp](Vector<Real>* Ax, const Vector<Real>& x){
    Vector<Real> U;
    KerOp.ComputePotential(U,x);
    (*Ax) = x*0.5 + U;
  };

  const Long N = KerOp.Dim(0);
  Vector<Real> sigma;
  if (V0.Dim() != N) {
    V0.ReInit(N);
    V0 = 1;
  }
  //comm.Barrier();
  //KerOp.Setup();
  //KerOp.ClearSetup();
  //comm.Barrier();
  //KerOp.Setup();
  //KerOp.ClearSetup();
  //comm.Barrier();
  //KerOp.Setup();
  Profile::Tic("Solve", &comm, true);
  ParallelSolver<Real> solver(comm);
  solver(&sigma, BIOp, V0, gmres_tol, 200);
  Profile::Toc();

  if (sigma_) (*sigma_) = sigma;
  elem_lst0.WriteVTK("vis/sigma", sigma, comm);
  Vector<Real> Utrg;
  if (Xt.Dim()) { // Evaluate at Xt
    KerOp.SetTargetCoord(Xt);

    Vector<Real> U;
    KerOp.ComputePotential(U,sigma);
    Utrg = U;
  }
  return Utrg;
}

template <class ValueType, class Real, class GeomFn> void GenericGeom(SlenderElemList<Real>& elem_lst0, const GeomFn& geom_fn, Vector<ValueType> panel_len, Vector<Long> FourierOrder, const Comm& comm, const Long ChebOrder, const Long DefaultFourierOrder, const Long DefaultNpanel, const ValueType s_len) {
  int Npanel = panel_len.Dim();
  if (!Npanel) { // use DefaultNpanel of equal length
    panel_len.ReInit(DefaultNpanel);
    panel_len = s_len/(ValueType)panel_len.Dim();
    GenericGeom<ValueType,Real,GeomFn>(elem_lst0, geom_fn, panel_len, FourierOrder, comm, ChebOrder, DefaultFourierOrder, DefaultNpanel, s_len);
    return;
  }
  if (FourierOrder.Dim() != Npanel) { // use DefaultFourierOrder
    FourierOrder.ReInit(Npanel);
    FourierOrder = DefaultFourierOrder;
  }

  const Long pid = comm.Rank();
  const Long Np = comm.Size();
  const Long k0 = Npanel*(pid+0)/Np;
  const Long k1 = Npanel*(pid+1)/Np;

  { // Set elem_lst0
    Vector<ValueType> panel_dsp(Npanel); panel_dsp = 0;
    omp_par::scan(panel_len.begin(), panel_dsp.begin(), Npanel);

    Vector<ValueType> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    for (Long k = k0; k < k1; k++) { // Init elem_lst0
      const auto& nds = SlenderElemList<ValueType>::CenterlineNodes(ChebOrder);
      for (Long i = 0; i < nds.Dim(); i++) {
        ValueType x, y, z, r;
        geom_fn(x, y, z, r, panel_dsp[k]+nds[i]*panel_len[k]);
        radius.PushBack(r);
        coord.PushBack(x);
        coord.PushBack(y);
        coord.PushBack(z);
      }
      cheb_order.PushBack(ChebOrder);
      fourier_order.PushBack(FourierOrder[k]);
    }
    elem_lst0.template Init<ValueType>(cheb_order, fourier_order, coord, radius);
  }
}

template <class ValueType, class Real, class GeomFn> void GenericGeomAdap(SlenderElemList<Real>& elem_lst0, const GeomFn& geom_fn, ValueType tol, const Comm& comm, const Long ChebOrder, const Long DefaultFourierOrder, const ValueType s_len) {
  const auto& nds1 = SlenderElemList<ValueType>::CenterlineNodes(ChebOrder);
  const auto& nds2 = SlenderElemList<ValueType>::CenterlineNodes(2*ChebOrder);

  Matrix<ValueType> Minterp(nds1.Dim(), nds2.Dim());
  Vector<ValueType> Minterp_(nds1.Dim()*nds2.Dim(), Minterp.begin(), false);
  LagrangeInterp<ValueType>::Interpolate(Minterp_, nds1, nds2);

  ValueType err = tol + 1;
  Vector<ValueType> max_err(4), max_val(4);
  Matrix<ValueType> X1(4,ChebOrder), X2(4,2*ChebOrder);
  Vector<ValueType> panel_len(1); panel_len = s_len;
  while(err > tol) {
    err = 0;
    ValueType panel_dsp = 0;
    Vector<ValueType> panel_len_new;
    for (Long k = 0; k < panel_len.Dim(); k++) {
      max_err = 0;
      max_val = 0;
      for (Long i = 0; i < nds1.Dim(); i++) {
        geom_fn(X1[0][i], X1[1][i], X1[2][i], X1[3][i], panel_dsp+nds1[i]*panel_len[k]);
      }
      Matrix<ValueType>::GEMM(X2, X1, Minterp);
      for (Long i = 0; i < nds2.Dim(); i++) {
        ValueType X_[4];
        geom_fn(X_[0], X_[1], X_[2], X_[3], panel_dsp+nds2[i]*panel_len[k]);
        for (Long j = 0; j < 4; j++) {
          max_err[j] = std::max<ValueType>(max_err[j], fabs(X_[j]-X2[j][i]));
          max_val[j] = std::max<ValueType>(max_val[j], fabs(X_[j]));
        }
      }

      ValueType rel_err_ = 0;
      for (Long j = 0; j < 4; j++) {
        rel_err_ = std::max<ValueType>(rel_err_, max_err[j]/max_val[j]);
      }
      if (rel_err_ > tol) {
        panel_len_new.PushBack(panel_len[k]/2);
        panel_len_new.PushBack(panel_len[k]/2);
      } else {
        panel_len_new.PushBack(panel_len[k]);
      }

      panel_dsp += panel_len[k];
      err = std::max<ValueType>(err, rel_err_);
    }
    std::cout<<"Npanels = "<<panel_len.Dim()<<"    Error = "<<err<<'\n';
    panel_len.Swap(panel_len_new);
  }

  GenericGeom<ValueType>(elem_lst0, geom_fn, panel_len, Vector<Long>(), comm, ChebOrder, DefaultFourierOrder, 0, s_len);
}

template <class ValueType, class Real> void GeomEllipse(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len, Vector<Long> FourierOrder, const Comm& comm, const ValueType Rmaj, const ValueType Rmin, const ValueType thickness, const Long ChebOrder) {
  auto loop_geom = [Rmaj,Rmin,thickness](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType s){
    ValueType theta = 2*const_pi<ValueType>()*s;
    x = Rmaj * cos<ValueType>(theta);
    y = Rmin * sin<ValueType>(theta);
    z = 0;
    r = thickness/2;
  };
  GenericGeom<ValueType>(elem_lst0, loop_geom, panel_len, FourierOrder, comm, ChebOrder, 14, 16, (ValueType)1);
}

template <class ValueType, class Real> void GeomTouchingTori(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len, Vector<Long> FourierOrder, const Comm& comm, const ValueType separation, const Long ChebOrder) {
  auto geom_fn = [separation](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType s){
    auto loop_geom1 = [](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType theta, ValueType x_shift){
      x = 2*cos<ValueType>(theta)+x_shift;
      y = 2*sin<ValueType>(theta);
      z = 0;
      r = 0.125;
    };
    auto loop_geom2 = [](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType theta, ValueType x_shift){
      x = 2*cos<ValueType>(theta)+x_shift;
      y = 0;
      z = 2*sin<ValueType>(theta);
      r = 0.125;
    };
    if (s<1) loop_geom1(x, y, z, r, 2*const_pi<ValueType>()*s, -(1.875-separation/2));
    else     loop_geom2(x, y, z, r, 2*const_pi<ValueType>()*s,  (1.875-separation/2));
  };
  GenericGeom<ValueType>(elem_lst0, geom_fn, panel_len, FourierOrder, comm, ChebOrder, 14, 32, (ValueType)2);
}

template <class ValueType, class Real> void GeomTangle(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len, Vector<Long> FourierOrder, const Comm& comm, const Long ChebOrder, ValueType tol) {
  std::default_random_engine g;
  g.seed(0);         // fix random seed
  std::normal_distribution<double> randn(0.0,1.0);

  // build a tangle complicated curve as Fourier cos & sin series...
  int tangle_freq = 10;   // max Fourier freq, 20 in linequad paper w/ Ludvig
  Vector<ValueType> co(6*tangle_freq);       // choose coeffs
  for (int j=1; j<=tangle_freq; ++j) {
    // pick curve vector j'th Fourier coeffs (real: 3 for cos, 3 for sin)
    ValueType ampl = 1.0 / (1.0 +(ValueType)j/3.0);  // 1/(5+|j|) in linequad
    for (int i=0;i<6;++i)
      co[i + 6*(j-1)] = ampl*randn(g);
  }

  auto geom_fn = [&tangle_freq,&co](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType s) { // Initialize elem_lst0
    x = 0; y = 0; z = 0;
    const ValueType theta = 2*const_pi<ValueType>()*s;
    for (int j=1; j<=tangle_freq; ++j) {
      // add in j'th Fourier mode w/ these coeffs, for this panel...
      ValueType cmode = cos<ValueType>(j*theta), smode = sin<ValueType>(j*theta);
      int o = 6*(j-1);     // offset index to j'th mode coeffs
      x += cmode*co[0+o] + smode*co[3+o];
      y += cmode*co[1+o] + smode*co[4+o];
      z += cmode*co[2+o] + smode*co[5+o];
    }
    r = ((ValueType)0.005) * (((ValueType)2)+sin<ValueType>(theta+sqrt<ValueType>(2)));
  };
  if (tol <= 0) {
    GenericGeom<ValueType>(elem_lst0, geom_fn, panel_len, FourierOrder, comm, ChebOrder, 14, 40, (ValueType)1);
  } else {
    GenericGeomAdap<ValueType>(elem_lst0, geom_fn, tol, comm, ChebOrder, 14, (ValueType)1);
  }
}

template <class ValueType, class Real> void GeomSphere(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len, Vector<Long> FourierOrder, const Comm& comm, const ValueType R, const Long ChebOrder) {
  auto geom_fn = [R](ValueType& x, ValueType& y, ValueType& z, ValueType& r, const ValueType theta){
    x = R * cos<ValueType>(theta/2);
    y = 0;
    z = 0;
    r = R * sin<ValueType>(theta/2);
  };
  GenericGeom<ValueType>(elem_lst0, geom_fn, panel_len, FourierOrder, comm, ChebOrder, 14, 24, (ValueType)1);
}


template <class Real> CubeVolumeVis<Real>::CubeVolumeVis(const Long N_, Real L, const Comm& comm) : N(N_) {
  const Long pid = comm.Rank();
  const Long Np = comm.Size();

  const Long NN = pow<COORD_DIM-1,Long>(N);
  const Long a = (N-1)*(pid+0)/Np;
  const Long b = (N-1)*(pid+1)/Np;
  N0 = b-a+1;
  if (N0<2) return;

  coord.ReInit(N0 * NN * COORD_DIM);
  for (Long i = 0; i < N0; i++) {
    for (Long j = 0; j < NN; j++) {
      for (Long k = 0; k < COORD_DIM; k++) {
        Long idx = ((i+a)*NN+j);
        coord[(i*NN+j)*COORD_DIM+k] = (((idx/pow<Long>(N,k)) % N)/(Real)(N-1)*2 - 1) * L;
      }
    }
  }
}

template <class Real> const Vector<Real>& CubeVolumeVis<Real>::GetCoord() const {
  return coord;
}

template <class Real> void CubeVolumeVis<Real>::GetVTUData(VTUData& vtu_data, const Vector<Real>& F) const {
  for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
  for (const auto& x :     F) vtu_data.value.PushBack((float)x);
  for (Long i = 0; i < N0-1; i++) {
    for (Long j = 0; j < N-1; j++) {
      for (Long k = 0; k < N-1; k++) {
        auto idx = [this](Long i, Long j, Long k) {
          return (i*N+j)*N+k;
        };
        vtu_data.connect.PushBack(idx(i+0,j+0,k+0));
        vtu_data.connect.PushBack(idx(i+0,j+0,k+1));
        vtu_data.connect.PushBack(idx(i+0,j+1,k+1));
        vtu_data.connect.PushBack(idx(i+0,j+1,k+0));
        vtu_data.connect.PushBack(idx(i+1,j+0,k+0));
        vtu_data.connect.PushBack(idx(i+1,j+0,k+1));
        vtu_data.connect.PushBack(idx(i+1,j+1,k+1));
        vtu_data.connect.PushBack(idx(i+1,j+1,k+0));
        vtu_data.offset.PushBack(vtu_data.connect.Dim());;
        vtu_data.types.PushBack(12);
      }
    }
  }
}
template <class Real> void CubeVolumeVis<Real>::WriteVTK(const std::string& fname, const Vector<Real>& F, const Comm& comm) const {
  VTUData vtu_data;
  GetVTUData(vtu_data, F);
  vtu_data.WriteVTK(fname, comm);
}



void commandline_option_start(int argc, char** argv, const char* help_text, const Comm& comm) {
  if (comm.Rank()) return;
  char help[] = "--help";
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], help)) {
      if (help_text != NULL) std::cout<<help_text<<'\n';
      std::cout<<"Usage:\n";
      std::cout<<"  "<<argv[0]<<" [options]\n\n";
    }
  }
}

const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, const char* err_msg, const Comm& comm){
  char help[] = "--help";
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], help)) {
      if (!comm.Rank()) std::cout<<"        "<<std::setw(10)<<opt<<" = ["<<std::setw(10)<<def_val<<"]  "<<(err_msg?err_msg:"")<<'\n';
      return def_val;
    }
  }

  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], opt)) {
      return argv[(i+1)%argc];
    }
  }
  if (required) {
    if (!comm.Rank()) std::cout<<"Missing: required option\n"<<"    "<<opt<<" : "<<(err_msg?err_msg:"")<<"\n\n";
    if (!comm.Rank()) std::cout<<"To see usage options\n"<<"    "<<argv[0]<<" --help\n\n";
    exit(0);
  }
  return def_val;
}

void commandline_option_end(int argc, char** argv) {
  char help[] = "--help";
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], help)) {
      exit(0);
    }
  }
}

bool to_bool(std::string str) {
  std::transform(str.begin(), str.end(), str.begin(), ::tolower);
  std::istringstream is(str);
  bool b;
  is >> std::boolalpha >> b;
  return b;
}

// Unused functions
#if 0
template <class Real> class SlenderVolumeVis {
    static constexpr Integer COORD_DIM = 3;
    static constexpr Integer Nlevels = 15;
  public:

    SlenderVolumeVis() = default;

    SlenderVolumeVis(const Vector<Long>& cheb_order_, const Vector<Long>& fourier_order_, const Vector<Real>& coord_, const Vector<Real>& radius_, const Vector<Real>& orientation_ = Vector<Real>(), Integer cheb_upsample = 1, Integer fourier_upsample = 1) {
      const auto concat_vecs = [](Vector<Real>& v, const Vector<Vector<Real>>& vec_lst) {
        const Long N = vec_lst.Dim();
        Vector<Long> dsp(N+1); dsp[0] = 0;
        for (Long i = 0; i < N; i++) {
          dsp[i+1] = dsp[i] + vec_lst[i].Dim();
        }
        if (v.Dim() != dsp[N]) v.ReInit(dsp[N]);
        for (Long i = 0; i < N; i++) {
          Vector<Real> v_(vec_lst[i].Dim(), v.begin()+dsp[i], false);
          v_ = vec_lst[i];
        }
      };
      fourier_order = fourier_order_ * fourier_upsample;
      cheb_order = cheb_order_ * cheb_upsample;

      Vector<Vector<Real>> X;
      SlenderElemList<Real> elem_lst;
      Vector<Real> radius = radius_*0.99;
      for (Long i = 0; i < Nlevels; i++) {
        Vector<Real> X_, s_param, sin_theta, cos_theta;
        elem_lst.Init(cheb_order, fourier_order, coord_, radius, orientation_);
        for (Long j = 0; j < cheb_order.Dim(); j++) {
          s_param.ReInit(cheb_order[j]);
          sin_theta.ReInit(fourier_order[j]);
          cos_theta.ReInit(fourier_order[j]);
          for (Long k = 0; k < cheb_order[j]; k++) {
            s_param[k] = 0.5 - cos<Real>(const_pi<Real>()*k/(cheb_order[j]-1))*0.5;
          }
          for (Long k = 0; k < fourier_order[j]; k++) {
            sin_theta[k] = sin<Real>(2*const_pi<Real>()*k/fourier_order[j]);
            cos_theta[k] = cos<Real>(2*const_pi<Real>()*k/fourier_order[j]);
          }
          elem_lst.GetGeom(&X_, nullptr, nullptr, nullptr, nullptr, s_param, sin_theta, cos_theta, j);
          X.PushBack(X_);
        }
        radius *= 0.95;
      }
      concat_vecs(coord, X);
    }

    const Vector<Real>& GetCoord() const {
      return coord;
    }

    void GetVTUData(VTUData& vtu_data, const Vector<Real>& F) const {
      Long level_offset = 0;
      Vector<Long> elem_offset(cheb_order.Dim());
      for (Long i = 0; i < cheb_order.Dim(); i++) {
        elem_offset[i] = level_offset;
        level_offset += cheb_order[i] * fourier_order[i];
      }

      for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
      for (const auto& x :     F) vtu_data.value.PushBack((float)x);
      for (Long i = 0; i < cheb_order.Dim(); i++) {
        for (Long l = 0; l < Nlevels-1; l++) {
          for (Long j = 0; j < cheb_order[i]-1; j++) {
            for (Long k = 0; k < fourier_order[i]; k++) {
              auto idx = [&level_offset,&elem_offset,this](Integer elem_idx, Integer level, Integer cheb_idx, Integer fourier_idx) {
                return elem_offset[elem_idx] + level_offset*level + cheb_idx*fourier_order[elem_idx] + fourier_idx%fourier_order[elem_idx];
              };
              vtu_data.connect.PushBack(idx(i, l+0,j+0,k+0));
              vtu_data.connect.PushBack(idx(i, l+0,j+0,k+1));
              vtu_data.connect.PushBack(idx(i, l+0,j+1,k+1));
              vtu_data.connect.PushBack(idx(i, l+0,j+1,k+0));
              vtu_data.connect.PushBack(idx(i, l+1,j+0,k+0));
              vtu_data.connect.PushBack(idx(i, l+1,j+0,k+1));
              vtu_data.connect.PushBack(idx(i, l+1,j+1,k+1));
              vtu_data.connect.PushBack(idx(i, l+1,j+1,k+0));
              vtu_data.offset.PushBack(vtu_data.connect.Dim());;
              vtu_data.types.PushBack(12);
            }
          }
        }
      }
    }
    void WriteVTK(const std::string& fname, const Vector<Real>& F, const Comm& comm) const {
      VTUData vtu_data;
      GetVTUData(vtu_data, F);
      vtu_data.WriteVTK(fname, comm);
    }

  private:

    Vector<Long> cheb_order;
    Vector<Long> fourier_order;
    Vector<Real> coord;
};
#endif

}


#ifndef SLENDERBODY_UTILS
#define SLENDERBODY_UTILS

#include <sctl.hpp>
#include <random>

namespace sctl {

template <class Real, class Kernel> void double_layer_test(const SlenderElemList<Real>& elem_lst0, const Comm& comm, Real tol) { // Double-layer identity test
  const Long pid = comm.Rank();

  Kernel ker_fn;
  BoundaryIntegralOp<Real,Kernel> BIOp(ker_fn, comm);
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

template <class Real, class KerSL, class KerDL, class KerGrad> void test_greens_identity(const SlenderElemList<Real>& elem_lst0, const Comm& comm, Real tol) {
  static constexpr Integer COORD_DIM = 3;
  const Long pid = comm.Rank();

  KerSL kernel_sl;
  KerDL kernel_dl;
  KerGrad kernel_grad;
  BoundaryIntegralOp<Real,KerSL> BIOpSL(kernel_sl, comm);
  BoundaryIntegralOp<Real,KerDL> BIOpDL(kernel_dl, comm);
  BIOpSL.AddElemList(elem_lst0);
  BIOpDL.AddElemList(elem_lst0);
  BIOpSL.SetAccuracy(tol);
  BIOpDL.SetAccuracy(tol);

  Vector<Real> X, Xn, Fs, Fd, Uref, Us, Ud;
  { // Set Fs, Fd, Uref
    elem_lst0.GetNodeCoord(&X, &Xn, nullptr);

    Vector<Real> X0{0.3,0.6,0.2}, Xn0{0,0,0}, F0(KerSL::SrcDim()), dU;
    for (auto& x : F0) x = drand48();
    kernel_sl  .Eval(Uref, X, X0, Xn0, F0*kernel_sl  .template ScaleFactor<Real>());
    kernel_grad.Eval(dU  , X, X0, Xn0, F0*kernel_grad.template ScaleFactor<Real>());

    Fd = -Uref;
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

  sctl::Profile::Enable(true);
  Profile::Tic("Setup+Eval", &comm);
  BIOpSL.ComputePotential(Us,Fs);
  BIOpDL.ComputePotential(Ud,Fd);
  Profile::Toc();

  Vector<Real> Uerr = (-Fd*0.5 + Us + Ud) - Uref;
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

template <class Real, class KerSL, class KerDL> Vector<Real> bvp_solve(const SlenderElemList<Real>& elem_lst0, const Real tol, const Real gmres_tol, const Real SLScaling = 1.0, Vector<Real> V0 = Vector<Real>(), const Vector<Real> Xt = Vector<Real>(), const Comm& comm = Comm::World()) { // Solve exterior Dirichlet BVP: (1/2 + D + S * SLScaling) sigma = V0
  KerSL ker_sl;
  KerDL ker_dl;
  BoundaryIntegralOp<Real,KerSL> SLOp(ker_sl, comm);
  BoundaryIntegralOp<Real,KerDL> DLOp(ker_dl, comm);
  SLOp.AddElemList(elem_lst0);
  DLOp.AddElemList(elem_lst0);
  SLOp.SetAccuracy(tol);
  DLOp.SetAccuracy(tol);
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
  ParallelSolver<Real> solver(comm);
  solver(&sigma, BIOp, V0, gmres_tol);
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


template <class Real> void GeomTangle(SlenderElemList<Real>& elem_lst0, Vector<Real> panel_len = Vector<Real>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const Long ChebOrder = 10) {
  int npan = panel_len.Dim();
  if (!npan) {
    panel_len.ReInit(40);
    panel_len = 1/(Real)panel_len.Dim();
    return GeomTangle(elem_lst0, panel_len, FourierOrder, comm, ChebOrder);
  }
  if (FourierOrder.Dim() != npan) {
    FourierOrder.ReInit(npan);
    FourierOrder = 14;
  }

  const Long pid = comm.Rank();
  const Long Np = comm.Size();
  const Long k0 = npan*(pid+0)/Np;
  const Long k1 = npan*(pid+1)/Np;

  std::default_random_engine g;
  g.seed(0);         // fix random seed
  std::normal_distribution<double> randn(0.0,1.0);

  // build a tangle complicated curve as Fourier cos & sin series...
  int tangle_freq = 10;   // max Fourier freq, 20 in linequad paper w/ Ludvig
  Vector<Real> co(6*tangle_freq);       // choose coeffs
  for (int j=1; j<=tangle_freq; ++j) {
    // pick curve vector j'th Fourier coeffs (real: 3 for cos, 3 for sin)
    Real ampl = 1.0 / (1.0 +(Real)j/3.0);  // 1/(5+|j|) in linequad
    for (int i=0;i<6;++i)
      co[i + 6*(j-1)] = ampl*randn(g);
  }

  { // Initialize elem_lst0
    Vector<Real> panel_dsp(npan); panel_dsp = 0;
    omp_par::scan(panel_len.begin(), panel_dsp.begin(), npan);

    Vector<Real> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    for (sctl::Long k = k0; k < k1; k++) {   // loop over panels...
      cheb_order.PushBack(ChebOrder);
      fourier_order.PushBack(FourierOrder[k]);

      const auto& nds = SlenderElemList<Real>::CenterlineNodes(ChebOrder);
      Vector<Real> elem_coord(nds.Dim()*3), elem_radius(nds.Dim());   // just for this pan
      elem_coord = 0.0;              // init ctrline nodes for this pan
      for (int j=1; j<=tangle_freq; ++j) {
        // add in j'th Fourier mode w/ these coeffs, for this panel...
        for (Long i = 0; i < nds.Dim(); i++) {     // loop over ctrline nodes
          // param (theta) domain for the i'th node in panel...
          Real theta = 2*const_pi<Real>()*(panel_dsp[k]+nds[i]*panel_len[k]);
          Real cmode = cos<Real>(j*theta), smode = sin<Real>(j*theta);
          int o = 6*(j-1);     // offset index to j'th mode coeffs
          elem_coord[i*3+0] += cmode*co[0+o] + smode*co[3+o];
          elem_coord[i*3+1] += cmode*co[1+o] + smode*co[4+o];
          elem_coord[i*3+2] += cmode*co[2+o] + smode*co[5+o];
        }
      }
      for (Long i = 0; i < nds.Dim(); i++) {         // fill radii, this pan
        Real theta = 2*const_pi<Real>()*(panel_dsp[k]+nds[i]*panel_len[k]);
        elem_radius[i] = ((Real)0.005) * (((Real)2)+sin<Real>(theta+sqrt<Real>(2)));
      }

      // append to coord and radius
      for (const auto x: elem_coord) coord.PushBack(x);
      for (const auto r: elem_radius) radius.PushBack(r);
    }
    elem_lst0.Init(cheb_order, fourier_order, coord, radius);     // send coords into this pan
  }
}

template <class Real> void GeomTorus(SlenderElemList<Real>& elem_lst0, Vector<Real> panel_len = Vector<Real>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const Real thickness = 0.01, const Long ChebOrder = 10) {
  int npan = panel_len.Dim();
  if (!npan) {
    panel_len.ReInit(16);
    panel_len = 1/(Real)panel_len.Dim();
    return GeomTorus(elem_lst0, panel_len, FourierOrder, comm, thickness, ChebOrder);
  }
  if (FourierOrder.Dim() != npan) {
    FourierOrder.ReInit(npan);
    FourierOrder = 14;
  }

  const Long pid = comm.Rank();
  const Long Np = comm.Size();
  const Long k0 = npan*(pid+0)/Np;
  const Long k1 = npan*(pid+1)/Np;

  { // Set elem_lst0
    auto loop_geom = [](Real& x, Real& y, Real& z, Real& r, const Real theta){
      x = cos<Real>(theta);
      y = sin<Real>(theta);
      z = 0.1*sin<Real>(theta-sqrt<Real>(2));
      r = 0.01*(2+sin<Real>(theta+sqrt<Real>(2)));
    };
    Vector<Real> panel_dsp(npan); panel_dsp = 0;
    omp_par::scan(panel_len.begin(), panel_dsp.begin(), npan);

    Vector<Real> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    for (Long k = k0; k < k1; k++) { // Init elem_lst0
      const auto& nds = SlenderElemList<Real>::CenterlineNodes(ChebOrder);
      for (Long i = 0; i < nds.Dim(); i++) {
        Real x, y, z, r;
        loop_geom(x, y, z, r, 2*const_pi<Real>()*(panel_dsp[k]+nds[i]*panel_len[k]));
        coord.PushBack(x);
        coord.PushBack(y);
        coord.PushBack(z);
        radius.PushBack(r);
      }
      cheb_order.PushBack(ChebOrder);
      fourier_order.PushBack(FourierOrder[k]);
    }
    elem_lst0.Init(cheb_order, fourier_order, coord, radius);
  }
}

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

template <class Real> class CubeVolumeVis {
    static constexpr Integer COORD_DIM = 3;
  public:

    CubeVolumeVis() = default;

    CubeVolumeVis(const Long N_, Real L) : N(N_) {
      coord.ReInit(pow<COORD_DIM,Long>(N) * COORD_DIM);
      for (Long i = 0; i < coord.Dim()/COORD_DIM; i++) {
        for (Long k = 0; k < COORD_DIM; k++) {
          coord[i*COORD_DIM+k] = (((i/pow<Long>(N,k)) % N)/(Real)(N-1)*2 - 1) * L;
        }
      }
    }

    const Vector<Real>& GetCoord() const {
      return coord;
    }

    void GetVTUData(VTUData& vtu_data, const Vector<Real>& F) const {
      for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
      for (const auto& x :     F) vtu_data.value.PushBack((float)x);
      for (Long i = 0; i < N-1; i++) {
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
    void WriteVTK(const std::string& fname, const Vector<Real>& F, const Comm& comm) const {
      VTUData vtu_data;
      GetVTUData(vtu_data, F);
      vtu_data.WriteVTK(fname, comm);
    }

  private:

    Long N;
    Vector<Real> coord;
};



void commandline_option_start(int argc, char** argv, const char* help_text = nullptr, const Comm& comm = Comm::Self()) {
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

const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, const char* err_msg, const Comm& comm = Comm::Self()){
  char help[] = "--help";
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], help)) {
      if (!comm.Rank()) std::cout<<"        "<<std::setw(10)<<opt<<" = ["<<std::setw(10)<<def_val<<"]  "<<err_msg<<'\n';
      return def_val;
    }
  }

  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], opt)) {
      return argv[(i+1)%argc];
    }
  }
  if (required) {
    if (!comm.Rank()) std::cout<<"Missing: required option\n"<<"    "<<err_msg<<"\n\n";
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

}

#endif

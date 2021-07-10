#include "sctl.hpp"

using namespace sctl;

template <class Real, class Kernel> void double_layer_test(const Comm& comm, Real tol) { // Double-layer identity test
  const Long pid = comm.Rank();
  const Long Np = comm.Size();

  SlenderElemList<Real> elem_lst0;
  //elem_lst0.Read("data/geom.data"); // Read geometry from file
  if (1) { // Initialize elem_lst0 in code
    const Long Nelem = 16;
    const Long ChebOrder = 10;
    const Long FourierOrder = 8;

    Vector<Real> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    const Long k0 = (Nelem*(pid+0))/Np;
    const Long k1 = (Nelem*(pid+1))/Np;
    for (Long k = k0; k < k1; k++) {
      cheb_order.PushBack(ChebOrder);
      fourier_order.PushBack(FourierOrder);
      const auto& nds = SlenderElemList<Real>::CenterlineNodes(ChebOrder);
      for (Long i = 0; i < nds.Dim(); i++) {
        Real theta = 2*const_pi<Real>()*(k+nds[i])/Nelem;
        coord.PushBack(cos<Real>(theta));
        coord.PushBack(sin<Real>(theta));
        coord.PushBack(0.1*sin<Real>(2*theta));
        radius.PushBack(0.01*(2+sin<Real>(theta+sqrt<Real>(2))));
      }
    }
    elem_lst0.Init(cheb_order, fourier_order, coord, radius);
  }

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
  elem_lst0.WriteVTK("Uerr_", Uerr, comm); // Write VTK
  { // Print error
    StaticArray<Real,2> max_err{0,0};
    for (auto x : Uerr) max_err[0] = std::max<Real>(max_err[0], fabs(x));
    comm.Allreduce(max_err+0, max_err+1, 1, Comm::CommOp::MAX);
    if (!pid) std::cout<<"Error = "<<max_err[1]<<'\n';
  }
  Profile::Enable(false);
  Profile::print(&comm);
}

template <class Real, class KerSL, class KerDL, class KerGrad> void test_greens_identity(const Comm& comm, Real tol) {
  auto loop_geom = [](Real& x, Real& y, Real& z, Real& r, const Real theta){
    x = cos<Real>(theta);
    y = sin<Real>(theta);
    z = 0.1*sin<Real>(theta-sqrt<Real>(2));
    r = 0.01*(2+sin<Real>(theta+sqrt<Real>(2)));
  };
  static constexpr Integer COORD_DIM = 3;
  const Long pid = comm.Rank();
  const Long Np = comm.Size();

  SlenderElemList<Real> elem_lst0;
  { // Set elem_lst0
    const Integer ChebOrder = 10;
    const Integer FourierOrder = 8;

    const Long Nelem = 32;
    const Long idx0 = Nelem*(pid+0)/Np;
    const Long idx1 = Nelem*(pid+1)/Np;

    Vector<Real> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    for (Long k = idx0; k < idx1; k++) { // Init elem_lst0
      const auto& nds = SlenderElemList<Real>::CenterlineNodes(ChebOrder);
      for (Long i = 0; i < nds.Dim(); i++) {
        Real x, y, z, r;
        loop_geom(x, y, z, r, 2*const_pi<Real>()*(k+nds[i])/Nelem);
        coord.PushBack(x);
        coord.PushBack(y);
        coord.PushBack(z);
        radius.PushBack(r);
      }
      cheb_order.PushBack(ChebOrder);
      fourier_order.PushBack(FourierOrder);
    }
    elem_lst0.Init(cheb_order, fourier_order, coord, radius);
  }

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
  elem_lst0.WriteVTK("Uerr0", Uerr, comm);
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

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    Comm comm = Comm::World();

    if (!comm.Rank()) std::cout<<"\n\n###########  DOUBLE-LAYER IDENTITY TEST WITH LAPLACE KERNEL  ###########\n";
    double_layer_test<double, Laplace3D_DxU>(comm, 1e-10);

    if (!comm.Rank()) std::cout<<"\n\n#############  GREEN'S IDENTITY TEST WITH LAPLACE KERNEL  ##############\n";
    test_greens_identity<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-10);

    if (!comm.Rank()) std::cout<<"\n\n###########  DOUBLE-LAYER IDENTITY TEST WITH STOKES KERNEL  ############\n";
    double_layer_test<double, Stokes3D_DxU>(comm, 1e-10);

    if (!comm.Rank()) std::cout<<"\n\n#############  GREEN'S IDENTITY TEST WITH Stokes KERNEL  ###############\n";
    test_greens_identity<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-10);
  }

  Comm::MPI_Finalize();
  return 0;
}


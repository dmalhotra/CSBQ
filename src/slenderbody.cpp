#include "sctl.hpp"

using namespace sctl;

template <class Real, class Kernel> void double_layer_test() { // Double-layer identity test
  SlenderElemList<Real> elem_lst0;
  elem_lst0.Read("data/geom.data"); // Read geometry from file
  if (0) { // Initialize elem_lst0 in code
    const Long Nelem = 16;
    const Long ChebOrder = 10;
    const Long FourierOrder = 14;

    Vector<Real> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    for (Long k = 0; k < Nelem; k++) {
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
  BoundaryIntegralOp<Real,Kernel> BIOp(ker_fn);
  BIOp.AddElemList(elem_lst0);
  BIOp.SetAccuracy(1e-10);

  // Warm-up run
  Vector<Real> F(BIOp.Dim(0)), U; F = 1;
  BIOp.ComputePotential(U,F);
  BIOp.ClearSetup();

  Profile::Enable(true);
  Profile::Tic("Setup+Eval");
  BIOp.ComputePotential(U,F);
  Profile::Toc();

  U = 0;
  BIOp.ComputePotential(U,F);

  Vector<Real> Uerr = U + 0.5;
  elem_lst0.WriteVTK("Uerr_", Uerr); // Write VTK
  { // Print error
    Real max_err = 0;
    for (auto x : Uerr) max_err = std::max<Real>(max_err, fabs(x));
    std::cout<<"Error = "<<max_err<<'\n';
  }
  Profile::Enable(false);
  Profile::print();
}

int main(int argc, char** argv) {
  std::cout<<"\n\n###########  DOUBLE-LAYER IDENTITY TEST WITH STOKES KERNEL  ############\n";
  double_layer_test<double, Stokes3D_DxU>();

  std::cout<<"\n\n###########  DOUBLE-LAYER IDENTITY TEST WITH LAPLACE KERNEL  ###########\n";
  double_layer_test<double, Laplace3D_DxU>();

  // SlenderElemList<double>::test_greens_identity();

  return 0;
}


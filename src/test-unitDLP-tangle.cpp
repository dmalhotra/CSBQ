// test unit DLP on tangled complicated loop.
// edit of slenderbody.cpp.  Barnett 4/30/21

#include "sctl.hpp"
#include <random>

using namespace sctl;

template <class Real> void double_layer_test() { // Double-layer identity test

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
  SlenderElemList<Real> elem_lst0;
  { // Initialize elem_lst0
    int npan = 40;
    int elem_cheb_order = 10;
    int elem_fourier_order = 8;

    Vector<Real> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    for (sctl::Long k = 0; k < npan; k++) {   // loop over panels...
      cheb_order.PushBack(elem_cheb_order);
      fourier_order.PushBack(elem_fourier_order);

      const auto& nds = SlenderElemList<Real>::CenterlineNodes(elem_cheb_order);
      Vector<Real> elem_coord(nds.Dim()*3), elem_radius(nds.Dim());   // just for this pan
      elem_coord = 0.0;              // init ctrline nodes for this pan
      for (int j=1; j<=tangle_freq; ++j) {
        // add in j'th Fourier mode w/ these coeffs, for this panel...
        for (Long i = 0; i < nds.Dim(); i++) {     // loop over ctrline nodes
          // param (theta) domain for the i'th node in panel...
          Real theta = 2*const_pi<Real>()*(k+nds[i])/npan;
          Real cmode = cos<Real>(j*theta), smode = sin<Real>(j*theta);
          int o = 6*(j-1);     // offset index to j'th mode coeffs
          elem_coord[i*3+0] += cmode*co[0+o] + smode*co[3+o];
          elem_coord[i*3+1] += cmode*co[1+o] + smode*co[4+o];
          elem_coord[i*3+2] += cmode*co[2+o] + smode*co[5+o];
        }
      }
      for (Long i = 0; i < nds.Dim(); i++) {         // fill radii, this pan
        Real theta = 2*const_pi<Real>()*(k+nds[i])/npan;
        elem_radius[i] = 0.005 * (2+sin<Real>(theta+sqrt<Real>(2)));
      }

      // append to coord and radius
      for (const auto x: elem_coord) coord.PushBack(x);
      for (const auto r: elem_radius) radius.PushBack(r);
    }
    elem_lst0.Init(cheb_order, fourier_order, coord, radius);     // send coords into this pan
  }

  Laplace3D_DxU laplace_dl;
  BoundaryIntegralOp<Real,Laplace3D_DxU> LapDL(laplace_dl);
  LapDL.AddElemList(elem_lst0);
  LapDL.SetAccuracy(1e-10);

  // Warm-up run
  Vector<Real> F(LapDL.Dim(0)), U; F = 1;
  LapDL.ComputePotential(U,F);

  sctl::Profile::Enable(true);
  Profile::Tic("Setup+Eval");
  LapDL.ClearSetup();
  LapDL.ComputePotential(U,F);
  Profile::Toc();

  U = 0;
  Profile::Tic("Eval");
  LapDL.ComputePotential(U,F);
  Profile::Toc();

  Vector<Real> Uerr = U + 0.5;
  elem_lst0.WriteVTK("tangleUerr", Uerr); // Write VTK
  { // Print error
    Real max_err = 0;
    for (auto x : Uerr) max_err = std::max<Real>(max_err, fabs(x));
    std::cout<<"Error = "<<max_err<<'\n';
  }
  sctl::Profile::print();
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  double_layer_test<double>();

  Comm::MPI_Finalize();
  return 0;
}


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
   int npan = 40;
   Vector<SlenderElement<Real>> elem_lst0(npan);
   for (sctl::Long k = 0; k < elem_lst0.Dim(); k++) {   // loop over panels...
     const auto& nds = SlenderElement<Real>::CenterlineNodes();
     Vector<Real> coord(nds.Dim()*3), radius(nds.Dim());   // just for this pan
     coord = 0.0;              // init ctrline nodes for this pan
     for (int j=1; j<=tangle_freq; ++j) {
       // add in j'th Fourier mode w/ these coeffs, for this panel...
       for (Long i = 0; i < nds.Dim(); i++) {     // loop over ctrline nodes
         // param (theta) domain for the i'th node in panel...
         Real theta = 2*const_pi<Real>()*(k+nds[i])/(Real)elem_lst0.Dim();
         Real cmode = cos<Real>(j*theta), smode = sin<Real>(j*theta);
         int o = 6*(j-1);     // offset index to j'th mode coeffs
         coord[i*3+0] += cmode*co[0+o] + smode*co[3+o];
         coord[i*3+1] += cmode*co[1+o] + smode*co[4+o];
         coord[i*3+2] += cmode*co[2+o] + smode*co[5+o];
       }
     }
     for (Long i = 0; i < nds.Dim(); i++) {         // fill radii, this pan
       Real theta = 2*const_pi<Real>()*(k+nds[i])/(Real)elem_lst0.Dim();
       radius[i] = 0.005 * (2+sin<Real>(theta+sqrt<Real>(2)));
     }
     elem_lst0[k].Init(coord,radius);     // send coords into this pan
   }
   
    GenericKernel<sctl::Laplace3D_DxU> laplace_dl;
    BoundaryIntegralOp<Real,GenericKernel<sctl::Laplace3D_DxU>> LapDL(laplace_dl);
    LapDL.AddElemList(elem_lst0);

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

    Vector<Real> Uerr = U*(1/(4*const_pi<Real>())) - 0.5;
    SlenderElement<Real>::WriteVTK("tangleUerr", elem_lst0, Uerr); // Write VTK
    { // Print error
      Real max_err = 0;
      for (auto x : Uerr) max_err = std::max<Real>(max_err, fabs(x));
      std::cout<<"Error = "<<max_err<<'\n';
    }
    sctl::Profile::print();
  }

int main(int argc, char** argv) {
  double_layer_test<double>();
  return 0;
}


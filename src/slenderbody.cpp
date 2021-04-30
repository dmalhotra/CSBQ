#include "sctl.hpp"

using namespace sctl;

template <class Real> void double_layer_test() { // Double-layer identity test
    Vector<SlenderElement<Real>> elem_lst0(8);
    for (sctl::Long k = 0; k < elem_lst0.Dim(); k++) { // Init elem_lst0
      const auto& nds = SlenderElement<Real>::CenterlineNodes();
      Vector<Real> coord(nds.Dim()*3), radius(nds.Dim());
      for (Long i = 0; i < nds.Dim(); i++) {
        Real theta = 2*const_pi<Real>()*(k+nds[i])/(Real)elem_lst0.Dim();
        coord[i*3+0] = cos<Real>(theta);
        coord[i*3+1] = sin<Real>(theta);
        coord[i*3+2] = 0.1*sin<Real>(2*theta);
        radius[i] = 0.01*(2+sin<Real>(theta+sqrt<Real>(2)));
      }
      elem_lst0[k].Init(coord,radius);
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
    SlenderElement<Real>::WriteVTK("Uerr", elem_lst0, Uerr); // Write VTK
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


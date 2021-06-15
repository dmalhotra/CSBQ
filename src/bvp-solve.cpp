#include "sctl.hpp"

using namespace sctl;

template <class Real, bool EnableSLScaling> void bvp_solve(const Real thickness) { // Solve: (1/2 + D + S/thickness) sigma = V0, where V0 = 1
  SlenderElemList<Real> elem_lst0;
  { // Initialize elem_lst0 in code
    const Long Nelem = 12;
    const Long ChebOrder = 10;
    const Long FourierOrder = 8;

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
        radius.PushBack(thickness*(2+sin<Real>(theta+sqrt<Real>(2))));
      }
    }
    elem_lst0.Init(cheb_order, fourier_order, coord, radius);
  }

  Stokes3D_FxU ker_sl;
  Stokes3D_DxU ker_dl;
  BoundaryIntegralOp<Real,Stokes3D_FxU> SLOp(ker_sl);
  BoundaryIntegralOp<Real,Stokes3D_DxU> DLOp(ker_dl);
  SLOp.AddElemList(elem_lst0);
  DLOp.AddElemList(elem_lst0);
  auto BIOp = [&SLOp,&DLOp,&thickness](Vector<Real>* Ax, const Vector<Real>& x){
    Vector<Real> Usl, Udl;
    SLOp.ComputePotential(Usl,x);
    DLOp.ComputePotential(Udl,x);
    if (EnableSLScaling) {
      (*Ax) = x*0.5 + Udl + Usl*(1/thickness);
    } else {
      (*Ax) = x*0.5 + Udl + Usl;
    }
  };

  const Long N = SLOp.Dim(0);
  Vector<Real> sigma, V0(N); V0 = 1;
  ParallelSolver<Real> solver;
  solver(&sigma, BIOp, V0, 1e-6);

  // TODO: generate visualization of solution
}

int main(int argc, char** argv) {
  std::cout<<"\n\n## Solving : (1/2 + D + S/thickness) sigma = V0    (Stokes exterior Dirichlet BVP with SL-scaling)\n";
  bvp_solve<double, true>(0.001);

  std::cout<<"\n\n## Solving : (1/2 + D + S) sigma = V0              (Stokes exterior Dirichlet BVP without SL-scaling)\n";
  bvp_solve<double, false>(0.001);

  return 0;
}


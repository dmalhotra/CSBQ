#include "utils.hpp"
#include "sctl.hpp"

using namespace sctl;

template <class Real, class Kernel> void double_layer_test_(const Comm& comm, Real tol) { // Double-layer identity test
  SlenderElemList<Real> elem_lst0;
  GeomTorus(elem_lst0, Vector<Real>(), Vector<Long>(), comm, 0.01, 10);
  //elem_lst0.Read("data/geom.data");

  double_layer_test<Real, Kernel>(elem_lst0, comm, tol);
}

template <class Real, class KerSL, class KerDL, class KerGrad> void test_greens_identity_(const Comm& comm, Real tol) {
  SlenderElemList<Real> elem_lst0;
  GeomTorus(elem_lst0, Vector<Real>(), Vector<Long>(), comm, 0.01, 10);
  //elem_lst0.Read("data/geom.data");

  test_greens_identity<Real, KerSL, KerDL, KerGrad>(elem_lst0, comm, tol);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    Comm comm = Comm::World();

    if (!comm.Rank()) std::cout<<"\n\n###########  DOUBLE-LAYER IDENTITY TEST WITH LAPLACE KERNEL  ###########\n";
    double_layer_test_<double, Laplace3D_DxU>(comm, 1e-10);

    if (!comm.Rank()) std::cout<<"\n\n#############  GREEN'S IDENTITY TEST WITH LAPLACE KERNEL  ##############\n";
    test_greens_identity_<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-10);

    if (!comm.Rank()) std::cout<<"\n\n###########  DOUBLE-LAYER IDENTITY TEST WITH STOKES KERNEL  ############\n";
    double_layer_test_<double, Stokes3D_DxU>(comm, 1e-10);

    if (!comm.Rank()) std::cout<<"\n\n#############  GREEN'S IDENTITY TEST WITH Stokes KERNEL  ###############\n";
    test_greens_identity_<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-10);
  }

  Comm::MPI_Finalize();
  return 0;
}


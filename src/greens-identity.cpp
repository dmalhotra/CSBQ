#include <sctl.hpp>
#include "utils.hpp"

using namespace sctl;

template <class Real, class KerSL, class KerDL, class KerGrad> void test_greens_identity_torus(const Comm& comm, Real tol, Long Nelem, Long FourierOrder) {
  SlenderElemList<Real> elem_lst0;
  Vector<Real> panel_len(Nelem); panel_len = 1/(Real)Nelem;
  Vector<Long> FourierOrder_(Nelem); FourierOrder_ = FourierOrder;
  GeomEllipse(elem_lst0, panel_len, FourierOrder_, comm, (Real)1, (Real)1, (Real)0.5, 10);
  //GeomSphere(elem_lst0, panel_len, FourierOrder_, comm, (Real)1.0, 10); // Sphere
  test_greens_identity<Real, KerSL, KerDL, KerGrad>(elem_lst0, comm, tol, Vector<Real>{1,2,3});
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  if (1) { // Comparison with BIEST
    Comm comm = Comm::World();

    // Laplace                                                                                             // T_setup       T_eval          error
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-01,  2,  4); //  0.0009       0.0000    2.40383e-01
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-02,  2,  4); //  0.0017       0.0000    1.93179e-02
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-03,  4,  8); //  0.0079       0.0000    4.98477e-05
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-04,  4,  8); //  0.0094       0.0001    1.13478e-05
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-05,  6, 12); //  0.0260       0.0003    2.13988e-07
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-06,  6, 12); //  0.0297       0.0003    2.32068e-07
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-07,  6, 12); //  0.0363       0.0003    5.20561e-09
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-08,  8, 16); //  0.0733       0.0010    9.30810e-10
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-09, 10, 20); //  0.1381       0.0023    2.03763e-10
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-10, 10, 20); //  0.1564       0.0024    1.07048e-11
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-11, 12, 24); //  0.2669       0.0049    2.81498e-12
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-12, 14, 28); //  0.4255       0.0088    2.02177e-12

    // Stokes                                                                                          // T_setup       T_eval          error
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-01,  2,  4); //  0.0015       0.0000    2.56128e-01
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-02,  2,  4); //  0.0025       0.0000    2.44248e-02
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-03,  4,  8); //  0.0141       0.0003    8.99652e-05
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-04,  4,  8); //  0.0164       0.0003    1.94202e-05
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-05,  6, 12); //  0.0536       0.0014    1.39755e-06
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-06,  6, 12); //  0.0600       0.0014    2.36850e-08
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-07,  6, 12); //  0.0725       0.0016    1.82834e-08
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-08,  8, 16); //  0.1634       0.0044    1.58734e-09
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-09, 10, 20); //  0.3173       0.0100    1.23922e-10
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-10, 10, 20); //  0.3694       0.0106    6.99807e-12
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-11, 12, 40); //  1.4705       0.0493    4.94875e-12
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-12, 14, 28); //  1.1100       0.0357    1.71970e-12
  }

  Comm::MPI_Finalize();
  return 0;
}


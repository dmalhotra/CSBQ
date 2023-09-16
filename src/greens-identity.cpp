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

    // Laplace                                                                                                   N   T_setup  Setup-rate   T_eval         Error
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-01,  2,  4); //   80    0.0035       22857   0.0006   2.18095e-01
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-02,  2,  4); //   80    0.0042       19048   0.0006   3.23968e-02
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-03,  4,  8); //  320    0.0166       19277   0.0007   3.25571e-03
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-04,  4,  8); //  320    0.0205       15610   0.0007   2.91289e-05
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-05,  6, 12); //  720    0.0559       12880   0.0015   7.57798e-07
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-06,  6, 12); //  720    0.0666       10811   0.0015   2.35023e-07
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-07,  6, 12); //  720    0.0808        8911   0.0015   2.16675e-07
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-08,  8, 16); // 1280    0.1660        7711   0.0037   1.05686e-09
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-09, 10, 20); // 2000    0.3090        6472   0.0082   1.79366e-09
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-10, 10, 20); // 2000    0.3572        5599   0.0083   2.07719e-10
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-11, 12, 24); // 2880    0.6047        4763   0.0168   1.16382e-11
    test_greens_identity_torus<double, Laplace3D_FxU, Laplace3D_DxU, Laplace3D_FxdU>(comm, 1e-12, 14, 28); // 3920    0.9395        4172   0.0292   1.96540e-12

    // Stokes                                                                                                    N   T_setup  Setup-rate   T_eval         Error
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-01,  2,  4); //      240    0.0049       48980   0.0006   2.23282e-01
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-02,  2,  4); //      240    0.0066       36364   0.0006   2.77515e-02
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-03,  4,  8); //      960    0.0323       29721   0.0012   5.55396e-03
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-04,  4,  8); //      960    0.0398       24121   0.0013   1.72239e-04
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-05,  6, 12); //     2160    0.1266       17062   0.0041   2.22350e-07
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-06,  6, 12); //     2160    0.1493       14468   0.0043   1.31611e-06
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-07,  6, 12); //     2160    0.1795       12033   0.0050   2.17364e-08
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-08,  8, 16); //     3840    0.4019        9555   0.0131   7.60145e-09
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-09, 10, 20); //     6000    0.7554        7943   0.0283   1.46205e-10
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-10, 10, 20); //     6000    0.8640        6944   0.0304   3.56598e-11
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-11, 12, 24); //     8640    1.6651        5189   0.0577   2.09852e-12
    test_greens_identity_torus<double, Stokes3D_FxU, Stokes3D_DxU, Stokes3D_FxT>(comm, 1e-12, 14, 28); //    11760    2.6551        4429   0.0992   1.72616e-12
  }

  Comm::MPI_Finalize();
  return 0;
}


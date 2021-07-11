// test unit DLP on tangled complicated loop.
// edit of slenderbody.cpp.  Barnett 4/30/21

#include "utils.hpp"
#include "sctl.hpp"

using namespace sctl;

template <class Real> void double_layer_test_(const Comm& comm = Comm::World()) { // Double-layer identity test
  SlenderElemList<Real> elem_lst0;
  GeomTangle(elem_lst0, Vector<Real>(), Vector<Long>(), comm, 10);
  //elem_lst0.Read("data/tangle-adap.geom", comm);

  double_layer_test<Real,Laplace3D_DxU>(elem_lst0, comm, 1e-10);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  double_layer_test_<double>();

  Comm::MPI_Finalize();
  return 0;
}


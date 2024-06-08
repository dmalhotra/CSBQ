/**
 * This demo code shows how to use the class sctl::BoundaryIntegralOp. The
 * class is declared in SCTL/include/sctl/boundary_integral.hpp. In particular,
 * we show how to:
 *
 * 1) build a boundary integral operator.
 *
 * 2) evaluate potentials at on- and off-surface target points.
 *
 * To compile and run the code, start in the project root directory and run:
 * make bin/demo2-bio && ./bin/demo2-bio
 */

#include <csbq.hpp>
using namespace sctl;

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    const Comm comm = Comm::World(); // comm must have local scope since it needs to be destroyed before Comm::MPI_Finalize()

    Stokes3D_DxU ker; // Stokes double-layer kernel (@see SCTL/include/sctl/kernel_function.hpp)
    BoundaryIntegralOp<double,Stokes3D_DxU> BIOp(ker, false, comm); // construct the boundary integral operator
    BIOp.SetAccuracy(1e-10); // set quadrature accuracy

    SlenderElemList<double> elem_lst;
    elem_lst.Read<double>("data/ring.geom", comm); // read geometry data from file
    BIOp.AddElemList(elem_lst); // add geometry to the boundary integral operator

    // The target points can be specified as follows. If not set or Xt is empty,
    // then the default target points are the surface discretization nodes.
    //BIOp.SetTargetCoord(Xt);

    const Long Ninput = BIOp.Dim(0); // (local) input dimension of the operator
    const Long Noutput = BIOp.Dim(1); // (local) output dimension of the operator

    // Set the density function sigma at each surface discretization node.
    // The nodes can be queried using elem_lst.GetNodeCoord() (@see ex1_geometry.cpp).
    Vector<double> sigma(Ninput);
    sigma = 1;

    Vector<double> U;
    BIOp.ComputePotential(U, sigma); // compute the potential U
    SCTL_ASSERT(U.Dim() == Noutput); // U will be of size Noutput

    // visualizing the geometry and the potential
    elem_lst.WriteVTK("vis/U", U, comm);
  }

  Comm::MPI_Finalize();
  return 0;
}

/**
 * This demo code shows how to use the class sctl::SlenderElementList. The
 * class is declared in SCTL/include/sctl/slender_element.hpp. In particular,
 * we show how to:
 *
 * 1) construct a slender-body geometry from the coordinates of the centerline
 * and cross-sectional radius at each point on the centerline.
 *
 * 2) read/write the geometry from/to files (human readable).
 *
 * 3) get the surface discretization nodes and normals.
 *
 * 4) write VTK visualization of the geometry.
 *
 * To compile and run the code, start in the project root directory and run:
 * make bin/demo1-geometry && ./bin/demo1-geometry
 */

#include <csbq.hpp>
using namespace sctl;

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    const Comm comm = Comm::World();
    SCTL_ASSERT_MSG(comm.Size()==1, "\
        This demo is sequential. In a distributed memory implementation, each process\n\
        would build only it's local section of the geometry.");

    // Quadrature files have been precomputed for 10th order elements (in 's')
    // and Fourier order (in 'theta') in multiples of 4 up to 100.
    const Long ElemOrder = 10;
    const Long FourierOrder = 12;

    const Long Nelem = 8; // number of elements
    Vector<double> Xc, eps; // centerline coordinates, and cross-sectional radius
    Vector<Long> ElemOrderVec(Nelem), FourierOrderVec(Nelem); // order in 's' and 'theta' for each element
    for (Long i = 0; i < Nelem; i++) { // loop over panels (build the centerline for a circle)
      ElemOrderVec[i] = ElemOrder;
      FourierOrderVec[i] = FourierOrder;

      // discretization nodes within a element in [0,1] interval
      const Vector<double>& nodes = SlenderElemList<double>::CenterlineNodes(ElemOrderVec[i]);

      for (Long j = 0; j < ElemOrderVec[i]; j++) { // loop over panel nodes
        const double phi = 2 * M_PI * (i + nodes[j]) / Nelem; // circle parameterization phi in [0,2pi]
        Xc.PushBack(cos(phi)); // x-coordinate
        Xc.PushBack(sin(phi)); // y-coordinate
        Xc.PushBack(0       ); // z-coordinate
        eps.PushBack((cos(2 * phi) + 2) * 0.1); // cross-sectional radius
      }
    }

    // Initialize the element list
    SlenderElemList<double> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps);

    elem_lst.Write("vis/ring.geom", comm); // write geometry data to file
    elem_lst.Read<double>("vis/ring.geom", comm); // read geometry data from file

    Vector<double> X, Xn; // the surface discretization nodes and normals
    Vector<Long> element_wise_node_cnt; // number of nodes per element
    elem_lst.GetNodeCoord(&X, &Xn, &element_wise_node_cnt); // get the surface geometry

    // visualizing the geometry and surface normals
    elem_lst.WriteVTK("vis/ring", Xn, comm);
  }

  Comm::MPI_Finalize();
  return 0;
}


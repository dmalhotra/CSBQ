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
  const Long COORD_DIM = 3; // 3D
  using Real = double;

  {
    const Comm comm = Comm::World();
    SCTL_ASSERT_MSG(comm.Size()==1, "\
        This demo is sequential. In a distributed memory implementation, each process\n\
        would build only it's local section of the geometry.");

    // Quadrature files have been precomputed for 10th order elements (in 's')
    // and Fourier order (in 'theta') in multiples of 4 up to 100.
    const Long ElemOrder = 10;
    const Long FourierOrder = 12;

    // discretization nodes within a element in [0,1] interval
    const Vector<Real>& elem_nodes = SlenderElemList<Real>::CenterlineNodes(ElemOrder);

    const Long Nelem = 8; // number of elements
    Vector<Real> eps(Nelem*ElemOrder); // cross-sectional radius
    Vector<Real> Xc(Nelem*ElemOrder*COORD_DIM); // centerline coordinates
    Vector<Long> ElemOrderVec(Nelem), FourierOrderVec(Nelem); // order in 's' and 'theta' for each element
    for (Long i = 0; i < Nelem; i++) { // loop of elements (build the centerline for a circle)
      ElemOrderVec[i] = ElemOrder;
      FourierOrderVec[i] = FourierOrder;
      for (Long j = 0; j < ElemOrder; j++) { // loop over panel nodes
        const Real phi = 2*const_pi<Real>() * (i + elem_nodes[j]) / Nelem; // circle parameterization phi in [0,2pi]
        const Real x = cos(phi);
        const Real y = sin(phi);
        const Real z = 0;

        const Long node_idx = i * ElemOrder + j;
        Xc[node_idx*COORD_DIM+0] = x;
        Xc[node_idx*COORD_DIM+1] = y;
        Xc[node_idx*COORD_DIM+2] = z;
        eps[node_idx] = 0.1;
      }
    }

    // Initialize the element list
    SlenderElemList<Real> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps);

    elem_lst.Write("vis/ring.geom", comm); // write geometry data to file
    elem_lst.Read<Real>("vis/ring.geom", comm); // read geometry data from file

    Vector<Real> X, Xn; // the surface discretization nodes and normals
    Vector<Long> element_wise_node_cnt; // number of nodes per element
    elem_lst.GetNodeCoord(&X, &Xn, &element_wise_node_cnt); // get the surface geometry

    // visualizing the geometry and surface normals
    elem_lst.WriteVTK("vis/ring", Xn, comm);
  }

  Comm::MPI_Finalize();
  return 0;
}


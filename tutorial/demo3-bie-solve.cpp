/**
 * This demo code shows how to solve BIEs (boundary integral equations). In
 * particular, we show how to:
 *
 * 1) Define a custom kernel type.
 *
 * 2) Build a BIO (boundary integral operator) with a custom kernel.
 *
 * 3) Solve a boundary integral equation using the BIO and GMRES.
 *
 * 4) Visualize the computed potentials.
 *
 * To compile and run the code, start in the project root directory and run:
 * make bin/demo3-bie-solve && ./bin/demo3-bie-solve
 */

#include <csbq.hpp>
using namespace sctl;

// Define a custom kernel for Stokes combined field (double-layer + scaled single-layer potential)
template <Long SL_scal> struct Stokes3D_CF_ {
  // Function to return the name of the kernel
  static const std::string& Name() {
    // Name determines what quadrature tables to use.
    // Single-layer quadrature tables also works for combined fields.
    // Not compatible with PVFMM, since it is not scale invariant
    // (single- and double-layer kernels have different scales)
    static const std::string name = "Stokes3D-FxU";
    return name;
  }

  // Function to return the number of floating point operations
  static constexpr Integer FLOPS() {
    return 50;
  }

  // Function to return the scale factor for the kernel
  template <class Real> static constexpr Real uKerScaleFactor() {
    return 1 / (8 * const_pi<Real>());
  }

  // Function to compute the kernel matrix
  template <Integer digits, class VecType>
  static void uKerMatrix(VecType (&u)[3][3], const VecType (&r)[3], const VecType (&n)[3], const void* ctx_ptr) {
    using Real = typename VecType::ScalarType;
    const auto SL_scal_ = VecType((Real)SL_scal);
    const auto r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]; // Compute squared distance
    const auto rinv = approx_rsqrt<digits>(r2, r2 > VecType::Zero()); // Compute inverse square root
    const auto rinv2 = rinv * rinv;
    const auto rinv3 = rinv2 * rinv;
    const auto rinv5 = rinv3 * rinv2;
    const auto rdotn = r[0] * n[0] + r[1] * n[1] + r[2] * n[2];
    const auto rdotn_rinv5_6 = VecType((Real)6) * rdotn * rinv5;
    for (Integer i = 0; i < 3; i++) {
      for (Integer j = 0; j < 3; j++) {
        const auto ri_rj = r[i] * r[j];
        const auto ker_dl = ri_rj * rdotn_rinv5_6; // Double-layer kernel
        const auto ker_sl = (i == j ? rinv + ri_rj * rinv3 : ri_rj * rinv3); // Single-layer kernel
        u[i][j] = ker_dl + ker_sl * SL_scal_; // Combine kernels
      }
    }
  }
};

// Define a type alias for the custom kernel
using Stokes3D_CF = GenericKernel<Stokes3D_CF_<4>>;

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    const Comm comm = Comm::World();
    double tol = 1e-8; // Tolerance for accuracy

    SlenderElemList<double> elem_lst;
    elem_lst.Read<double>("data/ring.geom", comm); // read geometry data from file
    SCTL_ASSERT_MSG(elem_lst.Size() > 0, "Could not read geometry file.");

    Stokes3D_CF ker;
    BoundaryIntegralOp<double,Stokes3D_CF> LayerPotenOp(ker, false, comm); // Create the layer-potential
    LayerPotenOp.SetAccuracy(tol);                                         // Set the accuracy tolerance
    LayerPotenOp.AddElemList(elem_lst);                                    // Add the element list to the operator

    // Define the boundary-integral operator: (I/2 + D + S * scale_factor)
    const auto BIO = [&LayerPotenOp](Vector<double>* U, const Vector<double>& sigma) {
      LayerPotenOp.ComputePotential(*U, sigma);
      (*U) += sigma * 0.5;
    };

    Vector<double> sigma, U0(LayerPotenOp.Dim(0)); // Vectors for sigma (density), and U0 (boundary condition)
    U0 = 1;                                        // Set boundary condition to (Ux, Uy, Uz) = (1, 1, 1)
    GMRES<double> solver(comm);                    // Define GMRES solver
    solver(&sigma, BIO, U0, tol);                  // Solve the linear system
    elem_lst.WriteVTK("vis/sigma", sigma, comm);   // Write sigma to VTK file

    // Visualize the fluid velocity in the volume
    CubeVolumeVis<double> cube(50, 2.0, comm);   // Define a cubic volume with 50x50x50 points and side length 2
    LayerPotenOp.SetTargetCoord(cube.GetCoord());// Set target coordinates for the operator
    Vector<double> U;                            // Vector for storing the velocity field
    LayerPotenOp.ComputePotential(U, sigma);     // Evaluate the velocity field
    cube.WriteVTK("vis/bie-solution", U - 1);    // Write to a VTK file
  }

  Comm::MPI_Finalize();
  return 0;
}

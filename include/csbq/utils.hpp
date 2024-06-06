#ifndef CSBQ_UTILS
#define CSBQ_UTILS

#include <sctl.hpp>
#include "csbq/slender_element.hpp"

namespace sctl {

/**
 * @brief Double-layer identity test.
 *
 * This function performs a double-layer identity test.
 *
 * @tparam Real The data type for real numbers, typically double or float.
 * @tparam Kernel The kernel function used in the test.
 * @param elem_lst0 The list of slender elements.
 * @param comm The communication object.
 * @param tol The tolerance level for the test.
 */
template <class Real, class Kernel> void double_layer_test(const SlenderElemList<Real>& elem_lst0, const Comm& comm, Real tol); // Double-layer identity test

/**
 * @brief Test Green's identity for Laplace and Stokes kernels.
 *
 * This function tests Green's identity for Laplace and Stokes kernels.
 *
 * @tparam Real The data type for real numbers, typically double or float.
 * @tparam KerSL The single-layer kernel.
 * @tparam KerDL The double-layer kernel.
 * @tparam KerGrad The gradient kernel.
 * @param elem_lst0 The list of slender elements.
 * @param comm The communication object.
 * @param tol The tolerance level for the test.
 * @param X0 The vector for testing coordinates. Default is {0.3, 0.6, 0.2}.
 */
template <class Real, class KerSL, class KerDL, class KerGrad> void test_greens_identity(const SlenderElemList<Real>& elem_lst0, const Comm& comm, const Real tol, const Vector<Real> X0 = Vector<Real>{0.3,0.6,0.2});

/**
 * @brief Solve exterior Dirichlet boundary value problem.
 *
 * This function solves the exterior Dirichlet boundary value problem: (1/2 + D + S * SLScaling) sigma = V0.
 * The solution at off-surface points Xt is returned.
 *
 * @tparam Real The data type for real numbers, typically double or float.
 * @tparam KerSL The single-layer kernel.
 * @tparam KerDL The double-layer kernel.
 * @tparam KerM2M The multipole-to-multipole kernel.
 * @tparam KerM2L The multipole-to-local kernel.
 * @tparam KerM2T The multipole-to-target kernel.
 * @tparam KerL2L The local-to-local kernel.
 * @tparam KerL2T The local-to-target kernel.
 * @param elem_lst0 The list of slender elements.
 * @param tol The tolerance level for the problem.
 * @param gmres_tol The GMRES solver tolerance.
 * @param SLScaling The scaling factor for the single-layer kernel. Default is 1.0.
 * @param V0 The right-hand side vector. Default is an empty vector.
 * @param Xt The vector of off-surface points. Default is an empty vector.
 * @param comm The communication object. Default is Comm::World().
 * @param sigma_ (Optional) Pointer to store the solution sigma.
 * @return The solution at off-surface points.
 */
template <class Real, class KerSL, class KerDL, class KerM2M, class KerM2L, class KerM2T, class KerL2L, class KerL2T> Vector<Real> bvp_solve(const SlenderElemList<Real>& elem_lst0, const Real tol, const Real gmres_tol, const Real SLScaling = 1.0, Vector<Real> V0 = Vector<Real>(), const Vector<Real> Xt = Vector<Real>(), const Comm& comm = Comm::World(), Vector<Real>* sigma_ = nullptr);

/**
 * @brief Solve exterior Dirichlet boundary value problem with combined kernels.
 *
 * This function solves the exterior Dirichlet boundary value problem: (1/2 + K) sigma = V0.
 * The solution at off-surface points Xt is returned.
 *
 * @tparam Real The data type for real numbers, typically double or float.
 * @tparam Ker The combined kernel function.
 * @param elem_lst0 The list of slender elements.
 * @param tol The tolerance level for the problem.
 * @param gmres_tol The GMRES solver tolerance.
 * @param V0 The right-hand side vector. Default is an empty vector.
 * @param Xt The vector of off-surface points. Default is an empty vector.
 * @param comm The communication object. Default is Comm::World().
 * @param sigma_ (Optional) Pointer to store the solution sigma.
 * @return The solution at off-surface points.
 */
template <class Real, class Ker> Vector<Real> bvp_solve_combined(const SlenderElemList<Real>& elem_lst0, const Real tol, const Real gmres_tol, Vector<Real> V0 = Vector<Real>(), const Vector<Real> Xt = Vector<Real>(), const Comm& comm = Comm::World(), Vector<Real>* sigma_ = nullptr);

/**
 * @brief Build geometry from a given geometry function.
 *
 * This function builds the geometry from a given geometry function: geom_fn(Real& x, Real& y, Real& z, Real& r, const Real s).
 *
 * @tparam ValueType The type of the value.
 * @tparam Real The data type for real numbers, typically double or float.
 * @tparam GeomFn The geometry function.
 * @param elem_lst0 The list of slender elements.
 * @param geom_fn The geometry function.
 * @param panel_len (Optional) Vector of panel lengths. Default is an empty vector.
 * @param FourierOrder (Optional) Vector of Fourier orders. Default is an empty vector.
 * @param comm The communication object. Default is Comm::Self().
 * @param ElemOrder The Chebyshev order. Default is 10.
 * @param DefaultFourierOrder The default Fourier order. Default is 14.
 * @param DefaultNpanel The default number of panels. Default is 16.
 * @param s_len The length in parameter space. Default is 1.
 */
template <class ValueType, class Real, class GeomFn> void GenericGeom(SlenderElemList<Real>& elem_lst0, const GeomFn& geom_fn, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const Long ElemOrder = 10, const Long DefaultFourierOrder = 14, const Long DefaultNpanel = 16, const ValueType s_len = (ValueType)1);

/**
 * @brief Adaptively build geometry from a given geometry function.
 *
 * This function adaptively builds the geometry from a given geometry function.
 *
 * @tparam ValueType The type of the value.
 * @tparam Real The data type for real numbers, typically double or float.
 * @tparam GeomFn The geometry function.
 * @param elem_lst0 The list of slender elements.
 * @param geom_fn The geometry function.
 * @param tol The tolerance level for adaptation.
 * @param comm The communication object. Default is Comm::Self().
 * @param ElemOrder The Chebyshev order. Default is 10.
 * @param DefaultFourierOrder The default Fourier order. Default is 14.
 * @param s_len The length in parameter space. Default is 1.
 */
template <class ValueType, class Real, class GeomFn> void GenericGeomAdap(SlenderElemList<Real>& elem_lst0, const GeomFn& geom_fn, ValueType tol, const Comm& comm = Comm::Self(), const Long ElemOrder = 10, const Long DefaultFourierOrder = 14, const ValueType s_len = (ValueType)1);

/**
 * @brief Build a torus geometry.
 *
 * This function builds a torus geometry of radius 1 and given thickness.
 *
 * @tparam ValueType The type of the value.
 * @tparam Real The data type for real numbers, typically double or float.
 * @param elem_lst0 The list of slender elements.
 * @param panel_len (Optional) Vector of panel lengths. Default is an empty vector.
 * @param FourierOrder (Optional) Vector of Fourier orders. Default is an empty vector.
 * @param comm The communication object. Default is Comm::Self().
 * @param Rmaj The major radius of the torus. Default is 2.
 * @param Rmin The minor radius of the torus. Default is 0.5.
 * @param thickness The thickness of the torus. Default is 0.001.
 * @param ElemOrder The Chebyshev order. Default is 10.
 */
template <class ValueType, class Real> void GeomEllipse(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const ValueType Rmaj = 2, const ValueType Rmin = 0.5, const ValueType thickness = 0.001, const Long ElemOrder = 10);

/**
 * @brief Build geometry with two nearly touching tori.
 *
 * This function builds geometry with two nearly touching tori with specified separation.
 *
 * @tparam ValueType The type of the value.
 * @tparam Real The data type for real numbers, typically double or float.
 * @param elem_lst0 The list of slender elements.
 * @param panel_len (Optional) Vector of panel lengths. Default is an empty vector.
 * @param FourierOrder (Optional) Vector of Fourier orders. Default is
 * @param comm MPI communicator.
 * @param separation Separation between the tori.
 * @param ElemOrder Element order.
 */
template <class ValueType, class Real> void GeomTouchingTori(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const ValueType separation = 0.01, const Long ElemOrder = 10);

/**
 * @brief Generate AB tangle geometry.
 *
 * @tparam ValueType Data type for geometric values.
 * @tparam Real Data type for real numbers.
 * @param elem_lst0 Slender element list.
 * @param panel_len Panel lengths.
 * @param FourierOrder Fourier orders.
 * @param comm MPI communicator.
 * @param ElemOrder Element order.
 * @param tol Tolerance for the geometry generation.
 */
template <class ValueType, class Real> void GeomTangle(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const Long ElemOrder = 10, ValueType tol = -1);

/**
 * @brief Generate geometry of a sphere.
 *
 * @tparam ValueType Data type for geometric values.
 * @tparam Real Data type for real numbers.
 * @param elem_lst0 Slender element list.
 * @param panel_len Panel lengths.
 * @param FourierOrder Fourier orders.
 * @param comm MPI communicator.
 * @param R Radius of the sphere.
 * @param ElemOrder Element order.
 */
template <class ValueType, class Real> void GeomSphere(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const ValueType R = 1, const Long ElemOrder = 10);

/**
 * @brief Represents a uniformly discretized cube volume.
 *
 * @tparam Real Data type for real numbers.
 */
template <class Real> class CubeVolumeVis {
    static constexpr Integer COORD_DIM = 3;
  public:

    CubeVolumeVis() = default;

    /**
     * @brief Construct a new CubeVolumeVis object.
     *
     * @param N_ Number of discretization points along one edge of the cube.
     * @param L Length of one edge of the cube.
     * @param comm MPI communicator.
     */
    CubeVolumeVis(const Long N_, Real L, const Comm& comm = Comm::Self());

    /**
     * @brief Get the coordinates of the discretization points.
     *
     * @return const Vector<Real>& Vector containing the coordinates.
     */
    const Vector<Real>& GetCoord() const;

    /**
     * @brief Write the cube volume to a VTK file.
     *
     * @param fname File name.
     * @param F Data associated with the discretization points.
     */
    void WriteVTK(const std::string& fname, const Vector<Real>& F) const;

    /**
     * @brief Get VTU data.
     *
     * @param vtu_data VTU data object.
     * @param F Data associated with the discretization points.
     */
    void GetVTUData(VTUData& vtu_data, const Vector<Real>& F) const;

  private:

    Long N, N0;
    Comm comm;
    Vector<Real> coord;
};

/**
 * @brief Start parsing command-line options.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments array.
 * @param help_text Help text to be displayed (optional).
 * @param comm MPI communicator.
 */
void commandline_option_start(int argc, char** argv, const char* help_text = nullptr, const Comm& comm = Comm::Self());

/**
 * @brief Parse a command-line option.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments array.
 * @param opt Option string.
 * @param def_val Default value.
 * @param required Whether the option is required.
 * @param err_msg Error message to display if the option is not found (required).
 * @param comm MPI communicator.
 * @return const char* Value of the option.
 */
const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, const char* err_msg, const Comm& comm = Comm::Self());

/**
 * @brief End parsing command-line options.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments array.
 */
void commandline_option_end(int argc, char** argv);

bool to_bool(std::string str);

}

#include <csbq/utils.cpp>

#endif


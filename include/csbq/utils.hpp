#ifndef SLENDERBODY_UTILS
#define SLENDERBODY_UTILS

#include <sctl.hpp>

namespace sctl {

template <class Real, class Kernel> void double_layer_test(const SlenderElemList<Real>& elem_lst0, const Comm& comm, Real tol); // Double-layer identity test

/**
 * Test Green's identity for Laplace and Stokes kernels.
 */
template <class Real, class KerSL, class KerDL, class KerGrad> void test_greens_identity(const SlenderElemList<Real>& elem_lst0, const Comm& comm, const Real tol, const Vector<Real> X0 = Vector<Real>{0.3,0.6,0.2});

/**
 * Solve exterior Dirichlet BVP: (1/2 + D + S * SLScaling) sigma = V0
 * Solution at off-surface points Xt is returned.
 */
template <class Real, class KerSL, class KerDL, class KerM2M, class KerM2L, class KerM2T, class KerL2L, class KerL2T> Vector<Real> bvp_solve(const SlenderElemList<Real>& elem_lst0, const Real tol, const Real gmres_tol, const Real SLScaling = 1.0, Vector<Real> V0 = Vector<Real>(), const Vector<Real> Xt = Vector<Real>(), const Comm& comm = Comm::World(), Vector<Real>* sigma_ = nullptr);

/**
 * Solve exterior Dirichlet BVP: (1/2 + K) sigma = V0
 * Solution at off-surface points Xt is returned.
 */
template <class Real, class Ker> Vector<Real> bvp_solve_combined(const SlenderElemList<Real>& elem_lst0, const Real tol, const Real gmres_tol, Vector<Real> V0 = Vector<Real>(), const Vector<Real> Xt = Vector<Real>(), const Comm& comm = Comm::World(), Vector<Real>* sigma_ = nullptr);

/**
 * Build geometry from a given geometry function:
 * geom_fn(Real& x, Real& y, Real& z, Real& r, const Real s)
 */
template <class ValueType, class Real, class GeomFn> void GenericGeom(SlenderElemList<Real>& elem_lst0, const GeomFn& geom_fn, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const Long ChebOrder = 10, const Long DefaultFourierOrder = 14, const Long DefaultNpanel = 16, const ValueType s_len = (ValueType)1);

template <class ValueType, class Real, class GeomFn> void GenericGeomAdap(SlenderElemList<Real>& elem_lst0, const GeomFn& geom_fn, ValueType tol, const Comm& comm = Comm::Self(), const Long ChebOrder = 10, const Long DefaultFourierOrder = 14, const ValueType s_len = (ValueType)1);

/**
 * Build a torus geometry of radius 1 and given thickness.
 */
template <class ValueType, class Real> void GeomEllipse(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const ValueType Rmaj = 2, const ValueType Rmin = 0.5, const ValueType thickness = 0.001, const Long ChebOrder = 10);

/**
 * Build geometry with two nearly touching tori with specified separation.
 */
template <class ValueType, class Real> void GeomTouchingTori(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const ValueType separation = 0.01, const Long ChebOrder = 10);

/**
 * AB tangle geometry.
 */
template <class ValueType, class Real> void GeomTangle(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const Long ChebOrder = 10, ValueType tol = -1);

template <class ValueType, class Real> void GeomSphere(SlenderElemList<Real>& elem_lst0, Vector<ValueType> panel_len = Vector<ValueType>(), Vector<Long> FourierOrder = Vector<Long>(), const Comm& comm = Comm::Self(), const ValueType R = 1, const Long ChebOrder = 10);


/**
 * Uniform discretization of a cube, along with visualization routines.
 */
template <class Real> class CubeVolumeVis {
    static constexpr Integer COORD_DIM = 3;
  public:

    CubeVolumeVis() = default;

    CubeVolumeVis(const Long N_, Real L, const Comm& comm = Comm::Self());

    const Vector<Real>& GetCoord() const;

    void WriteVTK(const std::string& fname, const Vector<Real>& F) const;

    void GetVTUData(VTUData& vtu_data, const Vector<Real>& F) const;

  private:

    Long N, N0;
    Comm comm;
    Vector<Real> coord;
};



void commandline_option_start(int argc, char** argv, const char* help_text = nullptr, const Comm& comm = Comm::Self());

const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, const char* err_msg, const Comm& comm = Comm::Self());

void commandline_option_end(int argc, char** argv);

bool to_bool(std::string str);

}

#include <csbq/utils.cpp>

#endif


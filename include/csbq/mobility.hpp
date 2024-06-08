#ifndef CSBQ_MOBILITY_HPP
#define CSBQ_MOBILITY_HPP

#include <sctl.hpp>
#include "csbq/slender_element.hpp"

#define WRITE_VTK 1

namespace sctl {

  /**
   * @class RigidBodyList
   * @brief A class for representing a set of rigid slender bodies.
   *
   * Each object is partitioned into elements along its centerline.
   * Each element is described by the positions of a set of discretization
   * nodes along its centerline and the cross-sectional radius of the fiber at
   * each node.
   *
   * Note: For all member functions, the input/output arrays store data in AoS (Array of Structures) order.
   * For example, the array of coordinates is stored as {x1,y1,z1, x2,y2,z2, ..., xn,yn,zn}.
   *
   * @tparam Real The data type for real numbers, typically double or float.
   */
  template <class Real> class RigidBodyList {
    static constexpr Integer COORD_DIM = 3;  ///< Coordinate dimension, always 3 for 3D space.

    public:
    /**
     * @enum Geom
     * @brief Predefined geometries for the rigid bodies.
     */
    enum class Geom {
      Loop,     ///< Loop geometry.
      Bacteria  ///< Bacteria geometry.
    };

    /**
     * @brief Constructor. Initialize all objects to either loop or bacteria geometry.
     *
     * @param comm The communication object. Default is Comm::Self().
     * @param Nobj The number of objects to initialize. Default is 1.
     * @param loop_rad The radius for loop geometry. Default is 0.45.
     * @param geom_type The type of geometry to initialize. Default is Geom::Loop.
     */
    RigidBodyList(const Comm& comm = Comm::Self(), const Long Nobj = 1, const Real loop_rad = 0.45, const Geom& geom_type = Geom::Loop);

    /**
     * @brief Initialize a set of N rigid slender bodies.
     *
     * @tparam ValueType The type of values used for initialization.
     * @param obj_elem_cnt Vector of length containing the number of elements in each object.
     * @param Mr_lst (Optional) A 3x3 rotation matrix for each object. Default is the identity matrix for each object.
     * @param ElemOrder Vector containing the polynomial order of each element.
     * @param FourierOrder Vector containing the Fourier order of each element.
     * @param panel_len Length of each element in parameter space. Total length of each object in parameter space is 1.
     * @param X Position of each node on the centerline.
     * @param R Cross-sectional radius at each node on the centerline.
     * @param OrientVec (Optional) An orientation vector perpendicular to the centerline.
     */
    template <class ValueType>
    void Init(const Vector<Long>& obj_elem_cnt, const Vector<ValueType>& Mr_lst, const Vector<Long>& ElemOrder, const Vector<Long>& FourierOrder, const Vector<ValueType>& panel_len, const Vector<ValueType>& X, const Vector<ValueType>& R, const Vector<ValueType>& OrientVec);

    /**
     * @brief Write the rigid object to a file.
     *
     * @param fname The name of the file to write to.
     */
    void Write(const std::string& fname) const;

    /**
     * @brief Read rigid object from a file.
     *
     * @param fname The name of the file to read from.
     * @param Nobj (Optional) Load only the first Nobj objects. Default is -1 for all objects.
     */
    void Read(const std::string& fname, const Long Nobj = -1);

    /**
     * @brief Return a const reference to the underlying SlenderElemList, to use for layer-potential evaluation or querying geometry.
     *
     * @return const reference to the SlenderElemList.
     */
    const SlenderElemList<Real>& GetElemList() const;

    /**
     * @brief Get the object centers and orientations (given by the 3x3 rotation matrix).
     *
     * @param Xc Pointer to store the object centers. Can be nullptr.
     * @param Mr Pointer to store the object orientations. Can be nullptr.
     */
    void GetObjPosition(Vector<Real>* Xc = nullptr, Vector<Real>* Mr = nullptr) const;

    /**
     * @brief Given a velocity field (at each surface discretization node), return its projection to the space of rigid body motions.
     *
     * @param u_proj The projected velocity field.
     * @param u The input velocity field.
     */
    void RigidBodyMotionProj(Vector<Real>& u_proj, const Vector<Real>& u) const;

    /**
     * @brief Get the object center velocity and angular velocity from the velocity at the surface discretization nodes.
     *
     * @param Uc_loc The object center velocity.
     * @param Omega_c_loc The angular velocity.
     * @param U_loc The velocity at the surface discretization nodes.
     */
    void GetRigidBodyMotion(Vector<Real>& Uc_loc, Vector<Real>& Omega_c_loc, const Vector<Real>& U_loc) const;

    /**
     * @brief Update the position of the objects given the translation in the object centers and a rotation matrix.
     *
     * @param dXc_loc The translation of the object centers.
     * @param Mr_loc The rotation matrix.
     */
    void RigidBodyUpdate(const Vector<Real>& dXc_loc, const Vector<Real>& Mr_loc);

    /**
     * @brief Adaptively refine/coarsen the geometry.
     *
     * @param v The data vector to be resolved on the refined/coarsend mesh.
     * @param tol The tolerance for refinement.
     * @param relative (Optional) Whether the tolerance is relative. Default is true.
     */
    void RefineAdaptive(const Vector<Vector<Real>>& v, const Real tol, const bool relative = true);

    /**
     * @brief [[deprecated]] Static function to create a rotation matrix.
     *
     * @param Mr_lst The list of rotation matrices.
     * @param Omega The angular velocity.
     * @param dt The time step.
     */
    static void RotationMatrix(Vector<Real>& Mr_lst, const Vector<Real>& Omega, const Real dt);

    /**
     * @brief Apply a preconditioner to each object. Each object is rotated to the orientation and refinement of the reference geometry.
     * The operator matrix is applied, and then each object is rotated back and refined to the original orientation and refinement.
     *
     * @param Ax0 The result of applying matrix A to each object in x0.
     * @param x0 The input vector containing the values at each surface discretization node.
     * @param A The operator to apply for one object in the reference geometry.
     * @param ref_geom The reference geometry.
     */
    void ApplyPrecond(Vector<Real>& Ax0, const Vector<Real> x0, const Matrix<Real>& A, const RigidBodyList<Real>& ref_geom) const;

    /**
     * @brief Get the minimum separation between the rigid objects.
     *
     * @return The minimum separation distance.
     */
    Real GetMinDist() const;

    //private:

    template <Integer dof> static void Resample(Vector<Real>& x_trg0, const Vector<Real>& panel_len_trg, const Vector<Long>& ElemOrder_trg, const Vector<Long>& FourierOrder_trg, const Vector<Real>& x_src, const Vector<Real>& panel_len_src, const Vector<Long>& ElemOrder_src, const Vector<Long>& FourierOrder_src);

    static const Matrix<Real>& FourierResample(const Integer N0, const Integer N1);

    static const Matrix<Real>& FourierErrMatrix(const Integer FourierOrder);

    static const Matrix<Real>& ChebErrMatrix(const Integer ElemOrder);

    [[deprecated]]
    void UpdatePointwise(const Vector<Real>& Y_loc);

    static void InitGeom(Vector<Real>& X, Vector<Real>& R, Vector<Real>& OrientVec, Vector<Long>& ElemOrder, Vector<Long>& FourierOrder, Vector<Real>& panel_len, Vector<Real>& Mr_lst, Vector<Long>& cnt, Vector<Long>& dsp, const Long Nobj, const Real loop_rad, const Geom& geom_type);

    static void SlenderElemIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst);

    static void ObjIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm comm);

    template <class ValueType> static void InitElemList(Long& loc_elem_cnt, Long& loc_elem_dsp, SlenderElemList<Real>& elem_lst, const Vector<Long>& ElemOrder, const Vector<Long>& FourierOrder, const Vector<ValueType>& X, const Vector<ValueType>& R, const Vector<ValueType>& OrientVec, const Comm comm);

    /**
     * Returns the center of mass of each object boundary.
     */
    static void GetXc(Vector<Real>& Xc, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm& comm);

    static void GetNodeObjIdx(Vector<Long>& node_obj_idx, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm comm);

    /**
     * Returns the rigid body rotation velocities corresponding to given
     * translational and rotational velocities.
     */
    static void RigidBodyVelocity(Vector<Real>& U, const Vector<Real> Uobj, const Vector<Real> Omega_obj, const Vector<Real> Xc, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm comm);

    Comm comm_;
    Vector<Long> obj_elem_cnt, obj_elem_dsp; // panel count for each object (global)
    Vector<Real> Xc, Mr_lst; // Center of mass and rotation matrix of each object (global)

    Vector<Long> ElemOrder, FourierOrder; // Chebyshev and Fourier order for each panel (global)
    Vector<Real> panel_len; // s-parameter span for each panel (global)

    Vector<Real> X, R, OrientVec; // panel position, radius and orientation at centerline nodes (global)

    Long loc_elem_cnt, loc_elem_dsp; // local element count and displacement
    SlenderElemList<Real> elem_lst;
  };

  /**
   * @class Mobility
   * @brief Solve Stokes mobility problem on a RigidBodyList.
   *
   * @tparam Real The data type for real numbers, typically double or float.
   */
  template <class Real> class Mobility {
    static constexpr Integer COORD_DIM = 3;  ///< Coordinate dimension, always 3 for 3D space.

    public:
    /**
     * @brief Constructor to initialize the Mobility class.
     *
     * @param comm The communicator object.
     * @param gravity The value of the gravitational force. Default is 0.
     */
    Mobility(const Comm& comm, Real gravity = 0);

    /**
     * @brief Set scaling factor for the single-layer potential in the combined field formulation.
     *
     * The integral operator is u = (I/2 + D + scal*S).
     *
     * @param SL_scal The scaling factor for the single-layer potential.
     */
    void SetSLScaling(const Real SL_scal) const;

    /**
     * @brief Construct a preconditioner for a given geometry and store it in a file.
     *
     * @param geom The geometry of the rigid bodies.
     * @param precond_fname The filename to store the preconditioner.
     * @param quad_tol The quadrature accuracy tolerance. Default is 1e-14.
     */
    void BuildPrecond(const RigidBodyList<Real>& geom, const std::string& precond_fname, const Real quad_tol = 1e-14) const;

    /**
     * @brief Solve the mobility problem and get the solution velocity.
     *
     * @param geom The geometry of the rigid bodies.
     * @param Ubg The background velocity field at each surface discretization node.
     * @param tol The accuracy tolerance for the GMRES solver. Default is 1e-8.
     * @param quad_tol The quadrature accuracy tolerance. Default is 1e-14.
     * @param q (Optional) The solution vector for the integral equation. Default is nullptr.
     * @param Ngmres (Optional) The number of GMRES iterations required. Default is nullptr.
     * @return The velocity of each surface discretization node.
     */
    Vector<Real> ComputeVelocity(const RigidBodyList<Real>& geom, const Vector<Real>& Ubg, const Real tol = 1e-8, const Real quad_tol = 1e-14, Vector<Real>* q = nullptr, Long* Ngmres = nullptr) const;

    /**
     * @brief Compute one timestep.
     *
     * @tparam BgFlow The type of the background velocity field functor.
     * @param geom0 The geometry (will be overwritten by the updated geometry).
     * @param bg_flow The background velocity field functor.
     * @param dt The time-step size.
     * @param ode_solver The instance of spectral deferred correction solver.
     * @param time_step_tol The accuracy of time-stepping.
     * @param quad_tol The accuracy of quadrature.
     * @param gmres_tol The GMRES accuracy tolerance.
     * @return The time-step error estimate.
     */
    template <class BgFlow> Real TimeStep(RigidBodyList<Real>& geom0, const BgFlow& bg_flow, const Real dt, const SDC<Real>& ode_solver, const Real time_step_tol, const Real quad_tol, const Real gmres_tol) const;

    /**
     * @brief Solve the mobility problem with adaptive time-stepping in the time interval [0, T].
     *
     * @tparam BgFlow The type of the background velocity field functor.
     * @param geom The geometry (will be overwritten by the updated geometry).
     * @param dt The starting time-step size.
     * @param T The total time horizon.
     * @param bg_flow The background velocity field functor.
     * @param time_step_order The order of the time-stepper.
     * @param time_step_tol The accuracy of time-stepping.
     * @param quad_tol The accuracy of quadrature.
     * @param gmres_tol The GMRES accuracy tolerance.
     * @param geom_tol The tolerance for geometry updates.
     * @param out_path The output path for visualization and checkpointing.
     * @param idx The starting index for the frames. Default is 0.
     * @return The time-step error estimate.
     */
    template <class BgFlow> Real AdaptiveTimeStep(RigidBodyList<Real>& geom, Real& dt, const Real T, const BgFlow& bg_flow, const Integer time_step_order, const Real time_step_tol, const Real quad_tol, const Real gmres_tol, const Real geom_tol, const std::string& out_path, Long idx = 0) const;

    static void test(const Comm& comm, const std::string& geom_file, const typename RigidBodyList<Real>::Geom& geom_type, const Long Nobj, const Real loop_rad, const Long start_idx = 0, const Long ts_order = 5, Real dt = 0.1, const Real T = 1000, const Real time_step_tol = 1e-7, const Real gmres_tol = 1e-10, const Real quad_tol = 1e-10, const Real geom_tol = 1e-8, const std::string& precond = "", const std::string& out_path = "vis/");

    static void test_(const Comm& comm, Real gmres_tol = 1e-12, Real quad_tol = 1e-12, Real geom_tol = 1e-8);

    private:

    Vector<Real> ApplyPrecond(const Vector<Real>& x, const RigidBodyList<Real>& geom) const;

    const Comm comm_;
    const Stokes3D_FxU ker_FxU;
    const Stokes3D_FxT ker_FxT;
    const Stokes3D_DxU ker_DxU;
    const Stokes3D_FxUP ker_FxUP;
    const Stokes3D_FSxU ker_FSxU;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesFxU;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxT> BIOp_StokesFxT;
    mutable BoundaryIntegralOp<Real,Stokes3D_DxU> BIOp_StokesDxU;

    Real gravity_;
    mutable Real SL_scal, DL_scal;
    mutable Matrix<Real> Mprecond;
    mutable RigidBodyList<Real> geom_precond;
  };

}

#include <csbq/mobility.cpp>

#endif

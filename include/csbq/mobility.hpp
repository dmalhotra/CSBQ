#ifndef MOBILITY_HPP
#define MOBILITY_HPP

#include <sctl.hpp>

#define WRITE_VTK 1

namespace sctl {


  template <class Real> class BgFlow {
    static constexpr Integer COORD_DIM = 3;

    public:
      BgFlow(const sctl::Comm& comm) : comm_(comm) {
      }

      Vector<Real> Velocity(const Vector<Real>& Xt) const {
        Vector<Real> U(Xt.Dim()); U = 0;
        return U;
      }

    private:
      Comm comm_;
  };

  template <class Real> class RigidBodyList {

    static constexpr Integer COORD_DIM = 3;

    public:

    enum class Geom {
      Loop,
      Bacteria
    };

    RigidBodyList(const Comm& comm = Comm::Self(), const Long Nobj = 1, const Real loop_rad = 0.45, const Geom& geom_type = Geom::Loop);

    template <class ValueType> void Init(const Vector<Long>& obj_elem_cnt_, const Vector<ValueType>& Mr_lst_, const Vector<Long>& ChebOrder_, const Vector<Long>& FourierOrder_, const Vector<ValueType>& panel_len_, const Vector<ValueType>& X_, const Vector<ValueType>& R_, const Vector<ValueType>& OrientVec_);

    void Write(const std::string& fname) const;

    void Read(const std::string& fname, const Long Nobj_ = -1);

    const SlenderElemList<Real>& GetElemList() const;

    void GetObjPosition(Vector<Real>* Xc_ = nullptr, Vector<Real>* Mr_ = nullptr) const;

    void RigidBodyMotionProj(Vector<Real>& u_proj, const Vector<Real>& u) const;

    void GetRigidBodyMotion(Vector<Real>& Uc_loc, Vector<Real>& Omega_c_loc, const Vector<Real>& U_loc) const;

    void RigidBodyUpdate(const Vector<Real>& dXc_loc, const Vector<Real>& Mr_loc);

    void RefineAdaptive(const Vector<Vector<Real>>& v_, const Real tol, const bool relative = true);

    //[[deprecated]]
    static void RotationMatrix(Vector<Real>& Mr_lst, const Vector<Real>& Omega, const Real dt);

    void ApplyPrecond(Vector<Real>& Ax0, const Vector<Real> x0, const Matrix<Real>& A, const RigidBodyList<Real>& ref_geom) const;

    Real GetMinDist() const;

    //private:

    template <Integer dof> static void Resample(Vector<Real>& x_trg0, const Vector<Real>& panel_len_trg, const Vector<Long>& ChebOrder_trg, const Vector<Long>& FourierOrder_trg, const Vector<Real>& x_src, const Vector<Real>& panel_len_src, const Vector<Long>& ChebOrder_src, const Vector<Long>& FourierOrder_src);

    static const Matrix<Real>& FourierResample(const Integer N0, const Integer N1);

    static const Matrix<Real>& FourierErrMatrix(const Integer FourierOrder);

    static const Matrix<Real>& ChebErrMatrix(const Integer ChebOrder);

    [[deprecated]]
    void UpdatePointwise(const Vector<Real>& Y_loc);

    static void InitGeom(Vector<Real>& X, Vector<Real>& R, Vector<Real>& OrientVec, Vector<Long>& ChebOrder, Vector<Long>& FourierOrder, Vector<Real>& panel_len, Vector<Real>& Mr_lst, Vector<Long>& cnt, Vector<Long>& dsp, const Long Nobj, const Real loop_rad, const Geom& geom_type);

    static void SlenderElemIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst);

    static void ObjIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm comm);

    template <class ValueType> static void InitElemList(Long& loc_elem_cnt, Long& loc_elem_dsp, SlenderElemList<Real>& elem_lst, const Vector<Long>& ChebOrder, const Vector<Long>& FourierOrder, const Vector<ValueType>& X, const Vector<ValueType>& R, const Vector<ValueType>& OrientVec, const Comm comm);

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

    Vector<Long> ChebOrder, FourierOrder; // Chebyshev and Fourier order for each panel (global)
    Vector<Real> panel_len; // s-parameter span for each panel (global)

    Vector<Real> X, R, OrientVec; // panel position, radius and orientation at centerline nodes (global)

    Long loc_elem_cnt, loc_elem_dsp; // local element count and displacement
    SlenderElemList<Real> elem_lst;
  };

  template <class Real> class Mobility {

    static constexpr Integer COORD_DIM = 3;

    public:

    Mobility(const Comm& comm, Real gravity = 0);

    void SetSLScaling(const Real SL_scal_) const;

    void BuildPrecond(const RigidBodyList<Real>& geom, const std::string& precond_fname, const Real quad_tol = 1e-14) const;

    Vector<Real> ComputeVelocity(const RigidBodyList<Real>& geom, const Vector<Real>& Ubg, const Real tol = 1e-8, const Real quad_tol = 1e-14, Vector<Real>* q_ = nullptr, Long* Ngmres_ = nullptr) const;

    template <class BgFlow> Real TimeStep(RigidBodyList<Real>& geom0, const BgFlow& bg_flow, const Real dt, const SDC<Real>& ode_solver, const Real time_step_tol, const Real quad_tol, const Real gmres_tol) const;

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

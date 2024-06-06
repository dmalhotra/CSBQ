.. _mobility_hpp:

mobility.hpp
============

This header file provides classes and methods for representing rigid slender bodies and solving the Stokes mobility problem for these bodies. The main classes are `RigidBodyList` and `Mobility`.

Classes and Types
-----------------

.. doxygenclass:: sctl::RigidBodyList
..   :members:
..

    **Constructor**:

    - ``RigidBodyList(comm = Comm::Self(), Nobj = 1, loop_rad = 0.45, geom_type = Geom::Loop)``: Initialize all objects to either loop or bacteria geometry.

    **Methods**:

    - ``Init(obj_elem_cnt, Mr_lst, ElemOrder, FourierOrder, panel_len, X, R, OrientVec)``: Initialize a set of N rigid slender bodies.

    - ``Write(fname) const``: Write the rigid object to a file.

    - ``Read(fname, Nobj = -1)``: Read rigid object from a file.

    - ``GetElemList() const``: Return a const reference to the underlying SlenderElemList, to use for layer-potential evaluation or querying geometry.

    - ``GetObjPosition(Xc = nullptr, Mr = nullptr) const``: Get the object centers and orientations (given by the 3x3 rotation matrix).

    - ``RigidBodyMotionProj(u_proj, u) const``: Given a velocity field (at each surface discretization node), return its projection to the space of rigid body motions.

    - ``GetRigidBodyMotion(Uc_loc, Omega_c_loc, U_loc) const``: Get the object center velocity and angular velocity from the velocity at the surface discretization nodes.

    - ``RigidBodyUpdate(dXc_loc, Mr_loc)``: Update the position of the objects given the translation in the object centers and a rotation matrix.

    - ``RefineAdaptive(v, tol, relative = true)``: Adaptively refine/coarsen the geometry.

    - ``RotationMatrix(Mr_lst, Omega, dt)``: Static function to create a rotation matrix.

    - ``ApplyPrecond(Ax0, x0, A, ref_geom) const``: Apply a preconditioner to each object.

    - ``GetMinDist() const``: Get the minimum separation between the rigid objects.

|

.. doxygenclass:: sctl::Mobility
..   :members:
..

    **Constructor**:

    - ``Mobility(comm, gravity = 0)``: Initialize the Mobility class.

    **Methods**:

    - ``SetSLScaling(SL_scal) const``: Set scaling factor for the single-layer potential in the combined field formulation.

    - ``BuildPrecond(geom, precond_fname, quad_tol = 1e-14) const``: Construct a preconditioner for a given geometry and store it in a file.

    - ``ComputeVelocity(geom, Ubg, tol = 1e-8, quad_tol = 1e-14, q = nullptr, Ngmres = nullptr) const``: Solve the mobility problem and get the solution velocity.

    - ``TimeStep(geom0, bg_flow, dt, ode_solver, time_step_tol, quad_tol, gmres_tol) const``: Compute one timestep.

    - ``AdaptiveTimeStep(geom, dt, T, bg_flow, time_step_order, time_step_tol, quad_tol, gmres_tol, geom_tol, out_path, idx = 0) const``: Solve the mobility problem with adaptive time-stepping in the time interval [0, T].

|

.. raw:: html

   <div style="border-top: 1px solid"></div>
   <br>

.. literalinclude:: ../../include/csbq/mobility.hpp
   :language: c++

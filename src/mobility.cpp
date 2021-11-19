#include "utils.hpp"
#include "sctl.hpp"

using namespace sctl;

template <class Real> class BgFlow {
  static constexpr Integer COORD_DIM = 3;

  public:
    BgFlow(const sctl::Comm& comm) : comm_(comm), fmm(comm){
      if (!comm.Rank()) {
        Xs.PushBack(-1000);
        Xs.PushBack(-1000);
        Xs.PushBack(0);

        Xs.PushBack(1000);
        Xs.PushBack(1000);
        Xs.PushBack(0);

        Fs.PushBack(1000000);
        Fs.PushBack(0);
        Fs.PushBack(0);

        Fs.PushBack(-1000000);
        Fs.PushBack(0);
        Fs.PushBack(0);
      }

      fmm.SetKernels(ker_FxUP, ker_FxUP, ker_FxUP);
      fmm.SetAccuracy(14);

      fmm.AddSrc("StokesSL", ker_FxUP, ker_FxUP);
      fmm.SetSrcCoord("StokesSL", Xs);
      fmm.SetSrcDensity("StokesSL", Fs);

      fmm.AddTrg("StokesVel", ker_FxU, ker_FxU);
      fmm.SetKernelS2T("StokesSL", "StokesVel", ker_FxU);

      fmm.AddTrg("StokesStress", ker_FxT, ker_FxT);
      fmm.SetKernelS2T("StokesSL", "StokesStress", ker_FxT);
    }

    Vector<Real> Velocity(const Vector<Real>& Xt) const {
      Vector<Real> U;
      fmm.SetTrgCoord("StokesVel", Xt);
      fmm.EvalDirect(U, "StokesVel");
      return U;
    }

    Vector<Real> Traction(const Vector<Real>& Xt) const {
      Vector<Real> U;
      fmm.SetTrgCoord("StokesStress", Xt);
      fmm.EvalDirect(U, "StokesStress");
      return U;
    }

  private:
    Comm comm_;
    Vector<Real> Xs, Fs;
    const Stokes3D_FxU ker_FxU;
    const Stokes3D_FxT ker_FxT;
    const Stokes3D_FxUP ker_FxUP;
    mutable ParticleFMM<Real,COORD_DIM> fmm;
};

template <class Real> class RigidGeom {

    static constexpr Integer COORD_DIM = 3;

  public:

    RigidGeom(const Comm& comm) : comm_(comm) {
      const Long Nobj = 2;
      const Long Npanel = 32;

      obj_elem_cnt.ReInit(Nobj);
      obj_elem_cnt = Npanel;

      Vector<Real> X, R, panel_len; // panel position, radius and s-parameter span for each panel
      panel_len.ReInit(Nobj * Npanel);
      FourierOrder.ReInit(Nobj * Npanel);
      ChebOrder.ReInit(Nobj * Npanel);
      panel_len = 1/(Real)Npanel;
      FourierOrder = 14;
      ChebOrder = 10;
      Real separation  = 2.0;
      InitGeom(X, R, obj_elem_cnt, ChebOrder, panel_len, separation);

      InitElemList(LocChebOrder, LocFourierOrder, elem_lst, Xc, ChebOrder, FourierOrder, X, R, obj_elem_cnt, comm_);
    }

    const SlenderElemList<Real>& GetElemList() const {
      return elem_lst;
    }

    Vector<Long> GetChebOrder() const {
      return LocChebOrder;
    }

    Vector<Long> GetFourierOrder() const {
      return LocFourierOrder;
    }

    [[deprecated]]
    void RigidBodyUpdate_(const Vector<Real>& U_loc, const Real dt) { // deprecated
      auto cross_prod = [](StaticArray<Real,COORD_DIM>& AxB, const StaticArray<Real,COORD_DIM>& A, const StaticArray<Real,COORD_DIM>& B) {
        AxB[0] = A[1]*B[2] - A[2]*B[1];
        AxB[1] = A[2]*B[0] - A[0]*B[2];
        AxB[2] = A[0]*B[1] - A[1]*B[0];
      };

      Vector<Long> node_obj_idx;
      Vector<Long> obj_elem_dsp(obj_elem_cnt.Dim()); obj_elem_dsp = 0;
      GetNodeObjIdx(node_obj_idx, elem_lst, obj_elem_cnt, comm_);
      omp_par::scan(obj_elem_cnt.begin(), obj_elem_dsp.begin(), obj_elem_cnt.Dim());
      const Long Nloc = node_obj_idx.Dim();
      const Long Nobj = obj_elem_cnt.Dim();

      Vector<Long> node_cnt, node_dsp;
      Vector<Real> U, Uc, X0, X_loc, Xc, area;
      {
        Vector<Long> node_cnt_loc;
        RigidBodyMotionProj(U, U_loc);
        elem_lst.GetNodeCoord(&X_loc, nullptr, &node_cnt_loc);
        ObjIntegral(Xc, X_loc, elem_lst, obj_elem_cnt, comm_);
        ObjIntegral(Uc, U_loc, elem_lst, obj_elem_cnt, comm_);
        ObjIntegral(area, X_loc*0+1, elem_lst, obj_elem_cnt, comm_);
        X0 = Allgather(X_loc, comm_);
        U = Allgather(U, comm_);
        Uc /= area;
        Xc /= area;

        node_cnt = Allgather(node_cnt_loc, comm_);
        node_dsp.ReInit(node_cnt.Dim()); node_dsp = 0;
        omp_par::scan(node_cnt.begin(), node_dsp.begin(), node_cnt.Dim());
      }


      Vector<Real> Ur_loc(Nloc * COORD_DIM);
      for (Long i = 0; i < Nloc; i++) {
        const Long obj = node_obj_idx[i];
        for (Long k = 0; k < COORD_DIM; k++) {
          Ur_loc[i*COORD_DIM+k] = U_loc[i*COORD_DIM+k] - Uc[obj*COORD_DIM+k];
        }
      }
      Vector<Real> Ur = Allgather(Ur_loc, comm_);
      Vector<Real> X  = Allgather( X_loc, comm_);

      Vector<Real> Omega_c(Nobj*COORD_DIM);
      for (Long obj = 0; obj < Nobj; obj++) { // Compute Omega_c
        StaticArray<Real, COORD_DIM> Omega_orient;
        { // Compute Omega_orient
          Matrix<Real> M(3,3), U, S, Vt; M = 0;
          for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
            const Long elem = obj_elem_dsp[obj]+i;
            for (Long j = 0; j < node_cnt[elem]; j++) {
              const Long node_idx = node_dsp[elem] + j;
              for (Long k0 = 0; k0 < COORD_DIM; k0++) {
                for (Long k1 = 0; k1 < COORD_DIM; k1++) {
                  M[k0][k1] += Ur[node_idx*COORD_DIM+k0] * Ur[node_idx*COORD_DIM+k1];
                }
              }
            }
          }
          M.SVD(U, S, Vt);

          Long min_idx = 0;
          Real min_val = fabs(S[0][0]);
          for (Long k = 0; k < COORD_DIM; k++) {
            if (fabs(S[k][k]) < min_val) {
              min_val = fabs(S[k][k]);
              min_idx = k;
            }
          }
          for (Long k = 0; k < COORD_DIM; k++) {
            Omega_orient[k] = U[k][min_idx];
          }
        }

        Real V_dot_U = 0, V_dot_V = 0;
        for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
          const Long elem = obj_elem_dsp[obj]+i;
          for (Long j = 0; j < node_cnt[elem]; j++) {
            const Long node_idx = node_dsp[elem] + j;

            StaticArray<Real,COORD_DIM> R, V;
            for (Long k = 0; k < COORD_DIM; k++) {
              R[k] = X[node_idx*COORD_DIM+k] - Xc[obj*COORD_DIM+k];
            }
            cross_prod(V, Omega_orient, R);

            for (Long k = 0; k < COORD_DIM; k++) {
              V_dot_U += Ur[node_idx*COORD_DIM+k] * V[k];
              V_dot_V += V[k] * V[k];
            }
          }
        }
        for (Long k = 0; k < COORD_DIM; k++) {
          Omega_c[obj*COORD_DIM+k] = Omega_orient[k] * V_dot_U/V_dot_V;
        }
      }


      Vector<Real> Y(U.Dim());
      for (Long obj = 0; obj < obj_elem_cnt.Dim(); obj++) {
        StaticArray<Real,COORD_DIM> omega;
        for (Long k = 0; k < COORD_DIM; k++) omega[k] = Omega_c[obj*COORD_DIM+k];

        Matrix<Real> Mr;
        { // Build rotation matrix Mr
          Matrix<Real> P1(3,3), P2(3,3), R(3,3), W0(1,3, omega), W1, W2;
          W0 *= dt;

          Real theta = atan2(W0[0][1], W0[0][0]);
          P1[0][0] = cos<Real>(theta); P1[0][1] =-sin<Real>(theta); P1[0][2] = 0;
          P1[1][0] = sin<Real>(theta); P1[1][1] = cos<Real>(theta); P1[1][2] = 0;
          P1[2][0] =                0; P1[2][1] =                0; P1[2][2] = 1;
          W1 = W0 * P1;

          Real phi = atan2(W1[0][2], W1[0][0]);
          P2[0][0] = cos<Real>(phi); P2[0][1] = 0; P2[0][2] =-sin<Real>(phi);
          P2[1][0] =              0; P2[1][1] = 1; P2[1][2] =              0;
          P2[2][0] = sin<Real>(phi); P2[2][1] = 0; P2[2][2] = cos<Real>(phi);
          W2 = W1 * P2;

          Real w = W2[0][0];
          R[0][0] = 1; R[0][1] =            0; R[0][2] =            0;
          R[1][0] = 0; R[1][1] = cos<Real>(w); R[1][2] = sin<Real>(w);
          R[2][0] = 0; R[2][1] =-sin<Real>(w); R[2][2] = cos<Real>(w);

          Mr = P1 * P2 * R * P2.Transpose() * P1.Transpose();
        }

        for (Long i = 0; i < obj_elem_cnt[obj]; i++) { // Update Y
          const Long elem = obj_elem_dsp[obj]+i;
          for (Long j = 0; j < node_cnt[elem]; j++) {
            const Long node_idx = node_dsp[elem] + j;

            StaticArray<Real,COORD_DIM> dR;
            StaticArray<Real,COORD_DIM> Mr_dR;
            for (Long k = 0; k < COORD_DIM; k++) {
              dR[k] = X0[node_idx*COORD_DIM+k] - Xc[obj*COORD_DIM+k];
            }
            Matrix<Real> dR_(1,COORD_DIM, dR, false);
            Matrix<Real> Mr_dR_(1,COORD_DIM, Mr_dR, false);
            Matrix<Real>::GEMM(Mr_dR_, dR_, Mr);

            for (Long k = 0; k < COORD_DIM; k++) {
              Y[node_idx*COORD_DIM+k] = Xc[obj*COORD_DIM+k] + Uc[obj*COORD_DIM+k]*dt + Mr_dR[k];
            }
          }
        }
      }

      Vector<Real> Uerr;
      { // Set Uerr
        Vector<Real> U0;
        RigidBodyVelocity(U0, Uc, Omega_c, elem_lst, obj_elem_cnt, comm_);
        Uerr = Allgather(U0 - U_loc, comm_);// - U;

        //Vector<Real> U1;
        //RigidBodyMotionProj(U1, U_loc);
        //Uerr = Allgather(U1, comm_) - U;

        //Uerr = U0-U1;
      }
      if (!comm_.Rank()) { // Print error
        Real max_err = 0;
        Real max_val = 0;
        for (const auto& a : Uerr) max_err = std::max<Real>(max_err, fabs(a));
        for (const auto& a : U  ) max_val = std::max<Real>(max_val, fabs(a));
        std::cout<<"Stokes err : "<<max_err/max_val<<'\n';
      } else {
        Uerr.ReInit(0);
      }
      comm_.PartitionN(Uerr, U_loc.Dim());
      elem_lst.WriteVTK("vis/Uerr", Uerr, comm_); // Write VTK

      if (comm_.Rank()) Y.ReInit(0);
      comm_.PartitionN(Y, U_loc.Dim());
      UpdatePointwise(Y);
    }
    void RigidBodyUpdate(const Vector<Real>& U_loc_, const Real dt) {
      //TODO: remove all Allgather and parallelize correctly
      auto cross_prod = [](StaticArray<Real,COORD_DIM>& AxB, const StaticArray<Real,COORD_DIM>& A, const StaticArray<Real,COORD_DIM>& B) {
        AxB[0] = A[1]*B[2] - A[2]*B[1];
        AxB[1] = A[2]*B[0] - A[0]*B[2];
        AxB[2] = A[0]*B[1] - A[1]*B[0];
      };

      Vector<Real> U_loc;
      RigidBodyMotionProj(U_loc, U_loc_);

      Vector<Real> X_loc, Xc, Uc, Omega_c;
      { // Set X_loc, Xc, Uc, Omega_c
        Vector<Long> obj_elem_dsp(obj_elem_cnt.Dim()); obj_elem_dsp = 0;
        omp_par::scan(obj_elem_cnt.begin(), obj_elem_dsp.begin(), obj_elem_cnt.Dim());
        const Long Nobj = obj_elem_cnt.Dim();

        Vector<Real> area;
        Vector<Long> node_cnt_loc;
        elem_lst.GetNodeCoord(&X_loc, nullptr, &node_cnt_loc);
        ObjIntegral(Xc, X_loc, elem_lst, obj_elem_cnt, comm_);
        ObjIntegral(Uc, U_loc, elem_lst, obj_elem_cnt, comm_);
        ObjIntegral(area, X_loc*0+1, elem_lst, obj_elem_cnt, comm_);
        Xc /= area;
        Uc /= area;

        const Vector<Real> X0 = Allgather(X_loc, comm_);
        const Vector<Real> U  = Allgather(U_loc, comm_);

        Vector<Long> node_cnt, node_dsp;
        node_cnt = Allgather(node_cnt_loc, comm_);
        node_dsp.ReInit(node_cnt.Dim()); node_dsp = 0;
        omp_par::scan(node_cnt.begin(), node_dsp.begin(), node_cnt.Dim());

        Omega_c.ReInit(Nobj*COORD_DIM);
        for (Long obj = 0; obj < Nobj; obj++) { // Compute Omega_c
          StaticArray<Real, COORD_DIM> Omega_orient;
          { // Compute Omega_orient
            Matrix<Real> M(3,3); M = 0;
            for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
              const Long elem = obj_elem_dsp[obj]+i;
              for (Long j = 0; j < node_cnt[elem]; j++) {
                const Long node_idx = node_dsp[elem] + j;
                for (Long k0 = 0; k0 < COORD_DIM; k0++) {
                  for (Long k1 = 0; k1 < COORD_DIM; k1++) {
                    const Real Ur_k0 = U[node_idx*COORD_DIM+k0] - Uc[obj*COORD_DIM+k0];
                    const Real Ur_k1 = U[node_idx*COORD_DIM+k1] - Uc[obj*COORD_DIM+k1];
                    M[k0][k1] += Ur_k0 * Ur_k1;
                  }
                }
              }
            }
            Matrix<Real> U, S, Vt;
            M.SVD(U, S, Vt);

            Long min_idx = 0;
            Real min_val = fabs(S[0][0]);
            for (Long k = 0; k < COORD_DIM; k++) {
              if (fabs(S[k][k]) < min_val) {
                min_val = fabs(S[k][k]);
                min_idx = k;
              }
            }
            for (Long k = 0; k < COORD_DIM; k++) {
              Omega_orient[k] = U[k][min_idx];
            }
          }

          Real V_dot_U = 0, V_dot_V = 0;
          for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
            const Long elem = obj_elem_dsp[obj]+i;
            for (Long j = 0; j < node_cnt[elem]; j++) {
              const Long node_idx = node_dsp[elem] + j;

              StaticArray<Real,COORD_DIM> R, V;
              for (Long k = 0; k < COORD_DIM; k++) {
                R[k] = X0[node_idx*COORD_DIM+k] - Xc[obj*COORD_DIM+k];
              }
              cross_prod(V, Omega_orient, R);

              for (Long k = 0; k < COORD_DIM; k++) {
                const Real Ur_k = U[node_idx*COORD_DIM+k] - Uc[obj*COORD_DIM+k];
                V_dot_U += Ur_k * V[k];
                V_dot_V += V[k] * V[k];
              }
            }
          }
          for (Long k = 0; k < COORD_DIM; k++) {
            Omega_c[obj*COORD_DIM+k] = Omega_orient[k] * V_dot_U/V_dot_V;
          }
        }
      }

      Vector<Long> node_obj_idx;
      GetNodeObjIdx(node_obj_idx, elem_lst, obj_elem_cnt, comm_);
      const Long Nloc = node_obj_idx.Dim();

      Vector<Real> Y_loc(X_loc.Dim());
      Vector<Matrix<Real>> Mr_lst(Xc.Dim()/COORD_DIM);
      for (Long node_idx = 0; node_idx < Nloc; node_idx++) {
        const Long obj = node_obj_idx[node_idx];

        Matrix<Real>& Mr = Mr_lst[obj];
        if (Mr.Dim(0) * Mr.Dim(1) == 0) { // Build rotation matrix Mr
          StaticArray<Real,COORD_DIM> omega;
          for (Long k = 0; k < COORD_DIM; k++) omega[k] = Omega_c[obj*COORD_DIM+k];
          Matrix<Real> P1(3,3), P2(3,3), R(3,3), W0(1,3, omega), W1, W2;
          W0 *= dt;

          Real theta = atan2(W0[0][1], W0[0][0]);
          P1[0][0] = cos<Real>(theta); P1[0][1] =-sin<Real>(theta); P1[0][2] = 0;
          P1[1][0] = sin<Real>(theta); P1[1][1] = cos<Real>(theta); P1[1][2] = 0;
          P1[2][0] =                0; P1[2][1] =                0; P1[2][2] = 1;
          W1 = W0 * P1;

          Real phi = atan2(W1[0][2], W1[0][0]);
          P2[0][0] = cos<Real>(phi); P2[0][1] = 0; P2[0][2] =-sin<Real>(phi);
          P2[1][0] =              0; P2[1][1] = 1; P2[1][2] =              0;
          P2[2][0] = sin<Real>(phi); P2[2][1] = 0; P2[2][2] = cos<Real>(phi);
          W2 = W1 * P2;

          Real w = W2[0][0];
          R[0][0] = 1; R[0][1] =            0; R[0][2] =            0;
          R[1][0] = 0; R[1][1] = cos<Real>(w); R[1][2] = sin<Real>(w);
          R[2][0] = 0; R[2][1] =-sin<Real>(w); R[2][2] = cos<Real>(w);

          Mr = P1 * P2 * R * P2.Transpose() * P1.Transpose();
        }

        StaticArray<Real,COORD_DIM> dR;
        StaticArray<Real,COORD_DIM> Mr_dR;
        for (Long k = 0; k < COORD_DIM; k++) {
          dR[k] = X_loc[node_idx*COORD_DIM+k] - Xc[obj*COORD_DIM+k];
        }
        Matrix<Real> dR_(1,COORD_DIM, dR, false);
        Matrix<Real> Mr_dR_(1,COORD_DIM, Mr_dR, false);
        Matrix<Real>::GEMM(Mr_dR_, dR_, Mr);

        for (Long k = 0; k < COORD_DIM; k++) {
          Y_loc[node_idx*COORD_DIM+k] = Xc[obj*COORD_DIM+k] + Uc[obj*COORD_DIM+k]*dt + Mr_dR[k];
        }
      }

      { // Print error
        Vector<Real> U0;
        RigidBodyVelocity(U0, Uc, Omega_c, elem_lst, obj_elem_cnt, comm_);
        Vector<Real> Uerr = Allgather(U_loc_-U0, comm_);
        Vector<Real> U    = Allgather(U_loc_, comm_);

        Real max_err = 0, max_val = 0;
        for (const auto& a : U   ) max_val = std::max<Real>(max_val, fabs(a));
        for (const auto& a : Uerr) max_err = std::max<Real>(max_err, fabs(a));
        if (!comm_.Rank()) std::cout<<"Stokes err : "<<max_err/max_val<<'\n';
        elem_lst.WriteVTK("vis/Uerr0", U_loc_-U0, comm_); // Write VTK
      }
      UpdatePointwise(Y_loc);
    }
    void UpdatePointwise(const Vector<Real>& Y_loc) {
      const Long Nelem = ChebOrder.Dim();
      SCTL_ASSERT(FourierOrder.Dim() == Nelem);
      if (!Nelem) return;

      static constexpr Integer COORD_DIM = 3;
      Vector<Long> cnt0(Nelem), dsp0(Nelem); dsp0[0] = 0;
      Vector<Long> cnt1(Nelem), dsp1(Nelem); dsp1[0] = 0;
      for (Long i = 0; i < Nelem; i++) cnt0[i] = ChebOrder[i];
      for (Long i = 0; i < Nelem; i++) cnt1[i] = ChebOrder[i] * FourierOrder[i];
      omp_par::scan(cnt0.begin(), dsp0.begin(), Nelem);
      omp_par::scan(cnt1.begin(), dsp1.begin(), Nelem);

      const Vector<Real> Y = Allgather(Y_loc, comm_);
      SCTL_ASSERT(Y.Dim() == (cnt1[Nelem-1]+dsp1[Nelem-1])*COORD_DIM);

      Vector<Real> X, R;
      for (Long elem = 0; elem < Nelem; elem++) {
        for (Long i = 0; i < ChebOrder[elem]; i++) {
          StaticArray<Real,COORD_DIM> x;
          for (Long k = 0; k < COORD_DIM; k++) x[k] = 0;
          for (Long j = 0; j < FourierOrder[elem]; j++) {
            for (Long k = 0; k < COORD_DIM; k++) {
              x[k] += Y[(dsp1[elem] + i*FourierOrder[elem] + j) * COORD_DIM+k];
            }
          }
          for (Long k = 0; k < COORD_DIM; k++) x[k] /= FourierOrder[elem];
          for (Long k = 0; k < COORD_DIM; k++) X.PushBack(x[k]);

          Real r = 0;
          for (Long j = 0; j < FourierOrder[elem]; j++) {
            Real r2 = 0;
            for (Long k = 0; k < COORD_DIM; k++) {
              Real y_k = Y[(dsp1[elem] + i*FourierOrder[elem] + j) * COORD_DIM+k];
              r2 += (x[k] - y_k) * (x[k] - y_k);
            }
            r += sqrt<Real>(r2);
          }
          R.PushBack(r/FourierOrder[elem]);
        }
      }
      InitElemList(LocChebOrder, LocFourierOrder, elem_lst, Xc, ChebOrder, FourierOrder, X, R, obj_elem_cnt, comm_);
    }

    void ComputeObjIntegral(Vector<Real>& Lf_cheb, const Vector<Real>& f_cheb) const {
      ObjIntegral(Lf_cheb, f_cheb, elem_lst, obj_elem_cnt, comm_);
    }

    [[deprecated]]
    void MobilityNullSpaceCorrection(Vector<Real>& Lf_cheb, const Vector<Real>& f_cheb) const { // deprecated TODO: replace by RigidBodyMOtionProj?
      Vector<Real> Fobj, Tobj;
      ComputeObjForceTorque(Fobj, Tobj, f_cheb, elem_lst, obj_elem_cnt, comm_);
      RigidBodyVelocity(Lf_cheb, Fobj, Tobj, elem_lst, obj_elem_cnt, comm_);
    }

    void PrintObjForceTorque(const Vector<Real>& f_cheb) const {
      Vector<Real> Fobj, Tobj;
      ComputeObjForceTorque(Fobj, Tobj, f_cheb, elem_lst, obj_elem_cnt, comm_);
      if (!comm_.Rank()) {
        std::cout<<"F = ";
        for (const auto& x : Fobj) std::cout<<' '<<x;
        std::cout<<"\n";

        std::cout<<"T = ";
        for (const auto& x : Tobj) std::cout<<' '<<x;
        std::cout<<"\n";
      }
    }

    void RigidBodyMotionProj(Vector<Real>& u_proj, const Vector<Real>& u) {
      auto inner_prod = [this](const Vector<Real>& A, const Vector<Real>& B) {
        SCTL_ASSERT(A.Dim() == B.Dim());
        Vector<Real> AB(A.Dim()), IntegralAB;
        for (Long i = 0; i < A.Dim(); i++) AB[i] = A[i]*B[i];
        ObjIntegral(IntegralAB, AB, elem_lst, obj_elem_cnt, comm_);

        const Long Nobj = obj_elem_cnt.Dim();
        const Long dof = IntegralAB.Dim() / Nobj;
        Vector<Real> AdotB(Nobj); AdotB = 0;
        for (Long i = 0; i < Nobj; i++) {
          for (Long j = 0; j < dof; j++) {
            AdotB[i] += IntegralAB[i*dof+j];
          }
        }
        return AdotB;
      };

      Vector<Long> elem_node_cnt, elem_node_dsp, node_obj_idx;
      { // Set elem_node_cnt, elem_node_dsp, node_obj_idx
        Vector<Real> X;
        elem_lst.GetNodeCoord(&X, nullptr, &elem_node_cnt);
        elem_node_dsp.ReInit(elem_node_cnt.Dim()); elem_node_dsp = 0;
        omp_par::scan(elem_node_cnt.begin(), elem_node_dsp.begin(), elem_node_cnt.Dim());
        GetNodeObjIdx(node_obj_idx, elem_lst, obj_elem_cnt, comm_);
      }
      const Long Nloc = elem_node_cnt[elem_node_cnt.Dim()-1] + elem_node_dsp[elem_node_dsp.Dim()-1];
      const Long Nobj = obj_elem_cnt.Dim();

      Matrix<Real> Q(2*COORD_DIM, u.Dim());
      for (Integer k = 0; k < 2*COORD_DIM; k++) { // Set Q
        Vector<Real> Uc(Nobj*COORD_DIM), Tc(Nobj*COORD_DIM);
        Uc.SetZero();
        Tc.SetZero();
        if (k < COORD_DIM) {
          for (Long i = 0; i < Nobj; i++) Uc[i*COORD_DIM+k] = 1;
        } else {
          for (Long i = 0; i < Nobj; i++) Tc[i*COORD_DIM+k-COORD_DIM] = 1;
        }

        Vector<Real> q_k(Q.Dim(1), Q[k], false);
        RigidBodyVelocity(q_k, Uc, Tc, elem_lst, obj_elem_cnt, comm_);

        for (Integer j = 0; j < k; j++) { // q_k <-- orthonormalize(q_k, Q)
          Vector<Real> q_j(Q.Dim(1), Q[j], false);
          Vector<Real> qk_dot_qj = inner_prod(q_k, q_j);
          for (Long node_idx = 0; node_idx < Nloc; node_idx++) {
            const Long obj = node_obj_idx[node_idx];
            for(Long l = 0; l < COORD_DIM; l++) {
              q_k[node_idx*COORD_DIM+l] -= q_j[node_idx*COORD_DIM+l] * qk_dot_qj[obj];
            }
          }
        }
        Vector<Real> qk_dot_qk = inner_prod(q_k, q_k);
        for (Long node_idx = 0; node_idx < Nloc; node_idx++) { // Q[k] <-- normalize(q_k)
          const Long obj = node_obj_idx[node_idx];
          for(Long l = 0; l < COORD_DIM; l++) {
            q_k[node_idx*COORD_DIM+l] /= sqrt<Real>(qk_dot_qk[obj]);
          }
        }
      }

      if (u_proj.Dim() != u.Dim()) u_proj.ReInit(u.Dim());
      u_proj.SetZero();
      for (Integer k = 0; k < 2*COORD_DIM; k++) { // Compute u_proj
        const Vector<Real> q_k(Q.Dim(1), Q[k], false);
        Vector<Real> qk_dot_u = inner_prod(u, q_k);
        for (Long node_idx = 0; node_idx < Nloc; node_idx++) { // u_proj <-- Q qk_dot_u
          const Long obj = node_obj_idx[node_idx];
          for(Long l = 0; l < COORD_DIM; l++) {
            u_proj[node_idx*COORD_DIM+l] += q_k[node_idx*COORD_DIM+l] * qk_dot_u[obj];
          }
        }
      }
    }

  private:

    template <class ValueType> static Vector<ValueType> Allgather(const Vector<ValueType>& Y_loc, const Comm& comm) {
      const Long Np = comm.Size();
      const StaticArray<Long,1> cnt_loc{Y_loc.Dim()};

      Vector<Long> cnt(Np), dsp(Np); dsp = 0;
      comm.Allgather((ConstIterator<Long>)cnt_loc, 1, cnt.begin(), 1);
      omp_par::scan(cnt.begin(), dsp.begin(), cnt.Dim());

      Vector<ValueType> Y(dsp[Np-1] + cnt[Np-1]);
      comm.Allgatherv(Y_loc.begin(), Y_loc.Dim(), Y.begin(), cnt.begin(), dsp.begin());
      return Y;
    }

    static void InitGeom(Vector<Real>& X, Vector<Real>& R, const Vector<Long>& cnt, const Vector<Long>& ChebOrder, const Vector<Real>& panel_len, const Real separation) {
      const Long Nobj = cnt.Dim();

      Vector<Long> dsp;
      dsp.ReInit(Nobj); dsp = 0;
      omp_par::scan(cnt.begin(), dsp.begin(), cnt.Dim());

      auto loop_geom1 = [](Real& x, Real& y, Real& z, Real& r, const Real theta, Real x_shift){
        x = 2*cos<Real>(theta)+x_shift;
        y = 2*sin<Real>(theta);
        z = 0;
        r = 0.125;
      };
      auto loop_geom2 = [](Real& x, Real& y, Real& z, Real& r, const Real theta, Real x_shift){
        x = 2*cos<Real>(theta)+x_shift;
        y = 0;
        z = 2*sin<Real>(theta);
        r = 0.125;
      };

      X.ReInit(0);
      R.ReInit(0);
      for (Long i = 0; i < Nobj; i++) {
        Real s_dsp = 0;
        for (Long j = 0; j < cnt[i]; j++) {
          const Long ChebOrder_ = ChebOrder[dsp[i]+j];
          const auto& nds = SlenderElemList<Real>::CenterlineNodes(ChebOrder_);
          for (Long k = 0; k < ChebOrder_; k++) {
            Real x, y, z, r;
            Real s = s_dsp + nds[k]*panel_len[dsp[i]+j];
            if (i == 0) loop_geom1(x, y, z, r, 2*const_pi<Real>()*s, -(1.875-separation/2));
            else        loop_geom2(x, y, z, r, 2*const_pi<Real>()*s,  (1.875-separation/2));
            X.PushBack(x);
            X.PushBack(y);
            X.PushBack(z);
            R.PushBack(r);
          }
          s_dsp += panel_len[dsp[i]+j];
        }
      }
    }

    static void SlenderElemIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst) {
      Vector<Long> cnt, dsp;
      Vector<Real> X, Xn, wts, dist_far, F;
      elem_lst.GetFarFieldNodes(X, Xn, wts, dist_far, cnt, (Real)1);
      dsp.ReInit(cnt.Dim()); dsp = 0;
      omp_par::scan(cnt.begin(), dsp.begin(), cnt.Dim());
      elem_lst.GetFarFieldDensity(F, F_cheb);
      if (!wts.Dim()) {
        IntegralF.ReInit(0);
        return;
      }

      const Integer dof = F.Dim() / wts.Dim();
      IntegralF.ReInit(cnt.Dim() * dof); IntegralF = 0;
      for (Long elem = 0; elem < cnt.Dim(); elem++) {
        for (Long i = 0; i < cnt[elem]; i++) {
          for(Long k = 0; k < dof; k++) {
            IntegralF[elem*dof+k] += F[(dsp[elem]+i)*dof+k] * wts[dsp[elem]+i];
          }
        }
      }
    }

    static void ObjIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Comm comm) {
      Vector<Real> Fe, Fe_local;
      SlenderElemIntegral(Fe_local, F_cheb, elem_lst);
      Fe = Allgather(Fe_local, comm);

      const Long Nobj = obj_elem_cnt.Dim();
      Vector<Long> obj_elem_dsp(Nobj); obj_elem_dsp = 0;
      omp_par::scan(obj_elem_cnt.begin(), obj_elem_dsp.begin(), Nobj);
      const Long Nelem = obj_elem_dsp[Nobj-1] + obj_elem_cnt[Nobj-1];
      const Integer dof = Fe.Dim() / Nelem;

      IntegralF.ReInit(Nobj * dof); IntegralF = 0;
      for (Long obj = 0; obj < Nobj; obj++) {
        for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
          for (Long k = 0; k < dof; k++) {
            IntegralF[obj*dof+k] += Fe[(obj_elem_dsp[obj]+i)*dof+k];
          }
        }
      }
    }

    static void InitElemList(Vector<Long>& LocChebOrder, Vector<Long>& LocFourierOrder, SlenderElemList<Real>& elem_lst, Vector<Real>& Xc, const Vector<Long>& ChebOrder, const Vector<Long>& FourierOrder, const Vector<Real>& X, const Vector<Real>& R, const Vector<Long>& obj_elem_cnt, const Comm comm) {
      Vector<Real> X_, R_;
      if (!comm.Rank()) {
        LocChebOrder = ChebOrder;
        LocFourierOrder = FourierOrder;
        X_ = X;
        R_ = R;
      } else {
        LocChebOrder.ReInit(0);
        LocFourierOrder.ReInit(0);
      }

      Long Nnodes = 0;
      comm.PartitionW(LocChebOrder);
      comm.PartitionN(LocFourierOrder, LocChebOrder.Dim());
      for (const auto& a : LocChebOrder) Nnodes += a;
      comm.PartitionN(X_, Nnodes);
      comm.PartitionN(R_, Nnodes);

      elem_lst.Init(LocChebOrder, LocFourierOrder, X_, R_);
      //elem_lst.Read("low-res-debug/vis/geom972", comm); ///////////////////////////////////////////

      { // Set Xc
        Vector<Real> X_cheb_node_loc;
        elem_lst.GetNodeCoord(&X_cheb_node_loc, nullptr, nullptr);
        ObjIntegral(Xc, X_cheb_node_loc, elem_lst, obj_elem_cnt, comm);

        Vector<Real> Aobj, A_cheb_node_loc(X_cheb_node_loc.Dim()/COORD_DIM); A_cheb_node_loc = 1;
        ObjIntegral(Aobj, A_cheb_node_loc, elem_lst, obj_elem_cnt, comm);
        for (Long obj = 0; obj < obj_elem_cnt.Dim(); obj++) {
          for (Long k = 0; k < COORD_DIM; k++) {
            Xc[obj*COORD_DIM+k] /= Aobj[obj];
          }
        }

        //Vector<Real> Xe, Xe_local;
        //Vector<Real> Ae, Ae_local, AA(XX.Dim()/COORD_DIM); AA = 1;
        //SlenderElemIntegral(Xe_local, XX, elem_lst);
        //SlenderElemIntegral(Ae_local, AA, elem_lst);
        //Xe = Allgather(Xe_local, comm);
        //Ae = Allgather(Ae_local, comm);

        //Vector<Long> obj_elem_dsp(obj_elem_cnt.Dim()); obj_elem_dsp = 0;
        //omp_par::scan(obj_elem_cnt.begin(), obj_elem_dsp.begin(), obj_elem_cnt.Dim());

        //const Long Nobj = obj_elem_cnt.Dim();
        //Xc.ReInit(obj_elem_cnt.Dim() * COORD_DIM); Xc = 0;
        //Vector<Real> A(Nobj); A = 0;
        //for (Long obj = 0; obj < Nobj; obj++) {
        //  for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
        //    for (Long k = 0; k < COORD_DIM; k++) {
        //      Xc[obj*COORD_DIM+k] += Xe[(obj_elem_dsp[obj]+i)*COORD_DIM+k];
        //    }
        //    A[obj] += Ae[obj_elem_dsp[obj]+i];
        //  }
        //}
        //for (Long obj = 0; obj < Nobj; obj++) {
        //  for (Long k = 0; k < COORD_DIM; k++) {
        //    Xc[obj*COORD_DIM+k] /= A[obj];
        //  }
        //}
      }
      //if (!comm.Rank()) std::cout<<Xc<<'\n';
    }

    static void GetNodeObjIdx(Vector<Long>& node_obj_idx, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Comm comm) {
      auto sum = [](const Vector<Long>& V) {
        Long sum = 0;
        for (const auto& x : V) sum+= x;
        return sum;
      };
      Vector<Long> elem_node_dsp, elem_node_cnt, elem_node_cnt_loc;
      elem_lst.GetNodeCoord(nullptr, nullptr, &elem_node_cnt_loc);
      const Long Nloc = sum(elem_node_cnt_loc);

      elem_node_cnt = Allgather(elem_node_cnt_loc, comm);
      const Long Nelem = elem_node_cnt.Dim();
      const Long Nobj = obj_elem_cnt.Dim();

      elem_node_dsp.ReInit(Nelem); elem_node_dsp = 0;
      Vector<Long> obj_elem_dsp(Nobj); obj_elem_dsp = 0;
      omp_par::scan(elem_node_cnt.begin(), elem_node_dsp.begin(), Nelem);
      omp_par::scan(obj_elem_cnt.begin(), obj_elem_dsp.begin(), Nobj);

      node_obj_idx.ReInit(elem_node_dsp[Nelem-1] + elem_node_cnt[Nelem-1]);
      for (Long obj = 0; obj < Nobj; obj++) {
        for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
          const Long elem = obj_elem_dsp[obj] + i;
          for (Long j = 0; j < elem_node_cnt[elem]; j++) {
            const Long node_idx = elem_node_dsp[elem] + j;
            node_obj_idx[node_idx] = obj;
          }
        }
      }

      if (comm.Rank()) node_obj_idx.ReInit(0);
      comm.PartitionN(node_obj_idx, Nloc);
    }

    static void ComputeObjForceTorque(Vector<Real>& Fobj, Vector<Real>& Tobj, const Vector<Real>& f_cheb, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Comm comm) {
      Vector<Long> node_obj_idx;
      GetNodeObjIdx(node_obj_idx, elem_lst, obj_elem_cnt, comm);
      const Long Nloc = node_obj_idx.Dim();
      SCTL_ASSERT(f_cheb.Dim() == Nloc * COORD_DIM);

      Vector<Real> X_cheb, T_cheb, Xc;
      elem_lst.GetNodeCoord(&X_cheb, nullptr, nullptr);
      ObjIntegral(Xc, X_cheb, elem_lst, obj_elem_cnt, comm);
      ObjIntegral(Fobj, f_cheb, elem_lst, obj_elem_cnt, comm);
      { // Set T_cheb (torque)
        T_cheb.ReInit(X_cheb.Dim());
        for (Long i = 0; i < Nloc; i++) {
          const Long obj = node_obj_idx[i];

          StaticArray<Real,COORD_DIM> R;
          R[0] = X_cheb[i*COORD_DIM+0] - Xc[obj*COORD_DIM+0];
          R[1] = X_cheb[i*COORD_DIM+1] - Xc[obj*COORD_DIM+1];
          R[2] = X_cheb[i*COORD_DIM+2] - Xc[obj*COORD_DIM+2];

          StaticArray<Real,COORD_DIM> f;
          f[0] = f_cheb[i*COORD_DIM+0];
          f[1] = f_cheb[i*COORD_DIM+1];
          f[2] = f_cheb[i*COORD_DIM+2];

          StaticArray<Real,COORD_DIM> T;
          T[0] = R[1] * f[2] - R[2] * f[1];
          T[1] = R[2] * f[0] - R[0] * f[2];
          T[2] = R[0] * f[1] - R[1] * f[0];

          T_cheb[i*COORD_DIM+0] = T[0];
          T_cheb[i*COORD_DIM+1] = T[1];
          T_cheb[i*COORD_DIM+2] = T[2];
        }
      }
      ObjIntegral(Tobj, T_cheb, elem_lst, obj_elem_cnt, comm);
    }

    static void RigidBodyVelocity(Vector<Real>& U, const Vector<Real> Uobj, const Vector<Real> Omega_obj, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Comm comm) {
      Vector<Long> node_obj_idx;
      GetNodeObjIdx(node_obj_idx, elem_lst, obj_elem_cnt, comm);
      const Long Nloc = node_obj_idx.Dim();

      Vector<Real> X_cheb, area, Xc;
      elem_lst.GetNodeCoord(&X_cheb, nullptr, nullptr);
      ObjIntegral(Xc, X_cheb, elem_lst, obj_elem_cnt, comm);
      ObjIntegral(area, X_cheb*0+1, elem_lst, obj_elem_cnt, comm);
      Xc /= area;

      if (U.Dim() != Nloc*COORD_DIM) U.ReInit(Nloc*COORD_DIM);
      for (Long i = 0; i < Nloc; i++) {
        const Long obj = node_obj_idx[i];

        StaticArray<Real,COORD_DIM> R;
        R[0] = X_cheb[i*COORD_DIM+0] - Xc[obj*COORD_DIM+0];
        R[1] = X_cheb[i*COORD_DIM+1] - Xc[obj*COORD_DIM+1];
        R[2] = X_cheb[i*COORD_DIM+2] - Xc[obj*COORD_DIM+2];

        StaticArray<Real,COORD_DIM> TxR;
        TxR[0] = Omega_obj[obj*COORD_DIM+1] * R[2] - Omega_obj[obj*COORD_DIM+2] * R[1];
        TxR[1] = Omega_obj[obj*COORD_DIM+2] * R[0] - Omega_obj[obj*COORD_DIM+0] * R[2];
        TxR[2] = Omega_obj[obj*COORD_DIM+0] * R[1] - Omega_obj[obj*COORD_DIM+1] * R[0];

        for (Long k = 0; k < COORD_DIM; k++) {
          U[i*COORD_DIM+k] = Uobj[obj*COORD_DIM+k] + TxR[k];
        }
      }
    }

    const Comm comm_;
    Vector<Long> obj_elem_cnt; // panel count for each object
    Vector<Long> ChebOrder, FourierOrder; // Chebyshev and Fourier order for each panel
    Vector<Real> Xc; // Center of mass of each object

    Vector<Long> LocChebOrder, LocFourierOrder;
    SlenderElemList<Real> elem_lst;
};

template <class Real> class Mobility_ {
  public:

    Mobility_(const Comm& comm) : comm_(comm), BIOp_StokesFxU(ker_FxU, false, comm), BIOp_StokesFxT(ker_FxT, true, comm), BIOp_StokesDxU(ker_DxU, false, comm) {
      BIOp_StokesFxT.SetFMMKer(ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxT);
      BIOp_StokesDxU.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
    }

    Vector<Real> ComputeVelocity(const RigidGeom<Real>& geom, const Vector<Real>& Ubg, const Vector<Real>& Tbg, const Real tol = 1e-8, const Real quad_tol = 1e-14) const {
      Vector<Real> Xn;
      const SlenderElemList<Real>& elem_lst0 = geom.GetElemList();
      elem_lst0.GetNodeCoord(nullptr, &Xn, nullptr);

      BIOp_StokesFxU.AddElemList(elem_lst0);
      BIOp_StokesFxT.AddElemList(elem_lst0);
      BIOp_StokesFxU.SetAccuracy(quad_tol);
      BIOp_StokesFxT.SetAccuracy(quad_tol);

      if (1) { ////////////////////////////////////////////////////
        auto dot_prod = [](const Vector<Real>& Y, const Vector<Real>& Xn) {
          const Long N = Xn.Dim()/3;
          if (!N) return Vector<Real>();
          const Long dof = Y.Dim()/N/3;
          SCTL_ASSERT(Y.Dim() == N*dof*3);

          Vector<Real> Y_dot_Xn;
          for (Long i = 0; i < N; i++) {
            for (Long k = 0; k < dof; k++) {
              Real y_dot_xn = 0;
              y_dot_xn += Y[(i*dof+k)*3+0] * Xn[i*3+0];
              y_dot_xn += Y[(i*dof+k)*3+1] * Xn[i*3+1];
              y_dot_xn += Y[(i*dof+k)*3+2] * Xn[i*3+2];
              Y_dot_Xn.PushBack(y_dot_xn);
            }
          }
          return Y_dot_Xn;
        };

        BoundaryIntegralOp<Real,Stokes3D_FxT> BIOp_StokesFxT_(ker_FxT, false, comm_);
        BIOp_StokesFxT_.SetFMMKer(ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxT);
        BIOp_StokesFxT_.AddElemList(elem_lst0);
        BIOp_StokesFxT_.SetAccuracy(quad_tol*0.1);

        Vector<Real> TdotXn, TdotXn_;
        BIOp_StokesFxT.ComputePotential(TdotXn, Tbg);
        BIOp_StokesFxT_.ComputePotential(TdotXn_, Tbg);
        //elem_lst0.WriteVTK("T", TdotXn_, comm_); // Write VTK
        TdotXn_ = dot_prod(TdotXn_, Xn);

        elem_lst0.WriteVTK("TdotXn1", TdotXn, comm_); // Write VTK
        elem_lst0.WriteVTK("TdotXn2", TdotXn_, comm_); // Write VTK

        elem_lst0.WriteVTK("Terr", TdotXn_-TdotXn, comm_); // Write VTK

        Real err = 0;
        for (const auto a : (TdotXn_-TdotXn)) err = std::max(err, fabs(a));
        std::cout<<err<<'\n';
        exit(0);
      }

      Vector<Real> F;
      auto TractionOp = [this,&geom](Vector<Real>* Ax, const Vector<Real>& x){
        Vector<Real> TdotXn, Lf;
        BIOp_StokesFxT.ComputePotential(TdotXn, x);
        geom.MobilityNullSpaceCorrection(Lf, x);
        (*Ax) = x*0.5 + TdotXn + Lf;
      };

      ParallelSolver<Real> solver(comm_);
      solver(&F, TractionOp, Tbg*(Real)-1, tol, 500);
      { ////////////////////////////////////////////////////////////////////////////////
        Vector<Real> res;
        TractionOp(&res, F); res += Tbg;
        elem_lst0.WriteVTK("vis/Residual", res, comm_); // Write VTK
        elem_lst0.WriteVTK("vis/F", F, comm_); // Write VTK
        geom.PrintObjForceTorque(F);
      }
      if (0) { // build operator matrix ////////////////////////////////////////////////
        auto build_matrix = [this,&tol](std::function<void(Vector<Real>*, const Vector<Real>&)> Op, Vector<Real> b) {
          Long N = b.Dim(), Nglb;
          { // Set Nglb
            StaticArray<Long,1>  Nloc_{N}, Nglb_;
            comm_.Allreduce<Long>(Nloc_, Nglb_, 1, Comm::CommOp::SUM);
            Nglb = Nglb_[0];
          }

          Matrix<double> M;
          Vector<Real> x, Ax;
          if (!comm_.Rank()) M.ReInit(Nglb,Nglb);
          for (Long i = 0; i < Nglb; i++) {
            x.ReInit(comm_.Rank() ? 0 : Nglb); x = 0;
            if (!comm_.Rank()) x[i] = 1;
            comm_.PartitionN(x, N);
            Ax.ReInit(N); Ax = 0;
            Op(&Ax, x);

            comm_.PartitionN(Ax, (comm_.Rank()?0:1));
            if (!comm_.Rank()) {
              std::cout<<i<<'\n';
              for (Long j = 0; j < Nglb; j++) {
                M[i][j] = (double)Ax[j];
              }
            }
          }
          if (!comm_.Rank()) M.Write("vis/Mobility_mat");

          Vector<Real> F;
          ParallelSolver<Real> solver(comm_);
          solver(&F, Op, b, tol, 500);
          comm_.PartitionN(b, (comm_.Rank()?0:1));
          comm_.PartitionN(F, (comm_.Rank()?0:1));
          if (!comm_.Rank()) b.Write("vis/b_vec");
          if (!comm_.Rank()) F.Write("vis/x_vec");
          comm_.Barrier();
          exit(0);
        };
        build_matrix(TractionOp, Tbg*(Real)-1);
      }

      Vector<Real> U;
      BIOp_StokesFxU.ComputePotential(U, F);

      BIOp_StokesFxU.template DeleteElemList<SlenderElemList<Real>>();
      BIOp_StokesFxT.template DeleteElemList<SlenderElemList<Real>>();
      return U + Ubg;
    }

  private:
    const Comm comm_;
    const Stokes3D_FxU ker_FxU;
    const Stokes3D_FxT ker_FxT;
    const Stokes3D_FxU ker_DxU;
    const Stokes3D_FxUP ker_FxUP;
    const Stokes3D_FSxU ker_FSxU;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesFxU;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxT> BIOp_StokesFxT;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesDxU;
};

template <class Real> class Mobility {
  public:

    Mobility(const Comm& comm) : comm_(comm), BIOp_StokesFxU(ker_FxU, false, comm), BIOp_StokesFxT(ker_FxT, true, comm), BIOp_StokesDxU(ker_DxU, false, comm) {
      BIOp_StokesFxT.SetFMMKer(ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxT);
      BIOp_StokesDxU.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
    }

    Vector<Real> ComputeVelocity(const RigidGeom<Real>& geom, const Vector<Real>& Ubg, const Vector<Real>& Tbg, const Real tol = 1e-8, const Real quad_tol = 1e-14) const {
      Vector<Real> Xn;
      const SlenderElemList<Real>& elem_lst0 = geom.GetElemList();
      elem_lst0.GetNodeCoord(nullptr, &Xn, nullptr);

      BIOp_StokesFxU.AddElemList(elem_lst0);
      BIOp_StokesFxT.AddElemList(elem_lst0);
      BIOp_StokesDxU.AddElemList(elem_lst0);
      BIOp_StokesFxU.SetAccuracy(quad_tol);
      BIOp_StokesFxT.SetAccuracy(quad_tol);
      BIOp_StokesDxU.SetAccuracy(quad_tol);

      Vector<Real> F;
      auto TractionOp = [this,&geom](Vector<Real>* Ax, const Vector<Real>& x){
        Vector<Real> TdotXn, Lf;
        BIOp_StokesFxT.ComputePotential(TdotXn, x);
        geom.MobilityNullSpaceCorrection(Lf, x);
        (*Ax) = x*0.5 + TdotXn + Lf;
      };

      ParallelSolver<Real> solver(comm_);
      solver(&F, TractionOp, Tbg*(Real)-1, tol, 500);
      { ////////////////////////////////////////////////////////////////////////////////
        Vector<Real> res;
        TractionOp(&res, F); res += Tbg;
        elem_lst0.WriteVTK("vis/Residual", res, comm_); // Write VTK
        elem_lst0.WriteVTK("vis/F", F, comm_); // Write VTK
        geom.PrintObjForceTorque(F);
      }
      if (0) { // build operator matrix ////////////////////////////////////////////////
        auto build_matrix = [this,&tol](std::function<void(Vector<Real>*, const Vector<Real>&)> Op, Vector<Real> b) {
          Long N = b.Dim(), Nglb;
          { // Set Nglb
            StaticArray<Long,1>  Nloc_{N}, Nglb_;
            comm_.Allreduce<Long>(Nloc_, Nglb_, 1, Comm::CommOp::SUM);
            Nglb = Nglb_[0];
          }

          Matrix<double> M;
          Vector<Real> x, Ax;
          if (!comm_.Rank()) M.ReInit(Nglb,Nglb);
          for (Long i = 0; i < Nglb; i++) {
            x.ReInit(comm_.Rank() ? 0 : Nglb); x = 0;
            if (!comm_.Rank()) x[i] = 1;
            comm_.PartitionN(x, N);
            Ax.ReInit(N); Ax = 0;
            Op(&Ax, x);

            comm_.PartitionN(Ax, (comm_.Rank()?0:1));
            if (!comm_.Rank()) {
              std::cout<<i<<'\n';
              for (Long j = 0; j < Nglb; j++) {
                M[i][j] = (double)Ax[j];
              }
            }
          }
          if (!comm_.Rank()) M.Write("vis/Mobility_mat");

          Vector<Real> F;
          ParallelSolver<Real> solver(comm_);
          solver(&F, Op, b, tol, 500);
          comm_.PartitionN(b, (comm_.Rank()?0:1));
          comm_.PartitionN(F, (comm_.Rank()?0:1));
          if (!comm_.Rank()) b.Write("vis/b_vec");
          if (!comm_.Rank()) F.Write("vis/x_vec");
          comm_.Barrier();
          exit(0);
        };
        build_matrix(TractionOp, Tbg*(Real)-1);
      }

      Vector<Real> U;
      BIOp_StokesFxU.ComputePotential(U, F);

      BIOp_StokesFxU.template DeleteElemList<SlenderElemList<Real>>();
      BIOp_StokesFxT.template DeleteElemList<SlenderElemList<Real>>();
      return U + Ubg;
    }

  private:
    const Comm comm_;
    const Stokes3D_FxU ker_FxU;
    const Stokes3D_FxT ker_FxT;
    const Stokes3D_FxU ker_DxU;
    const Stokes3D_FxUP ker_FxUP;
    const Stokes3D_FSxU ker_FSxU;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesFxU;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxT> BIOp_StokesFxT;
    mutable BoundaryIntegralOp<Real,Stokes3D_FxU> BIOp_StokesDxU;
};

template <class Real> void test(const Comm& comm, Real tol, Real quad_eps) {
  auto dot_prod = [](const Vector<Real>& Y, const Vector<Real>& Xn) {
    const Long N = Xn.Dim()/3;
    if (!N) return Vector<Real>();
    const Long dof = Y.Dim()/N/3;
    SCTL_ASSERT(Y.Dim() == N*dof*3);

    Vector<Real> Y_dot_Xn;
    for (Long i = 0; i < N; i++) {
      for (Long k = 0; k < dof; k++) {
        Real y_dot_xn = 0;
        y_dot_xn += Y[(i*dof+k)*3+0] * Xn[i*3+0];
        y_dot_xn += Y[(i*dof+k)*3+1] * Xn[i*3+1];
        y_dot_xn += Y[(i*dof+k)*3+2] * Xn[i*3+2];
        Y_dot_Xn.PushBack(y_dot_xn);
      }
    }
    return Y_dot_Xn;
  };

  Real dt = 0.1;
  RigidGeom<Real> geom(comm);

  const BgFlow<Real> bg_flow(comm);
  const Mobility<Real> stokes_mobility(comm);
  for (Long i = 0; i < 100000; i++) {
    const SlenderElemList<Real> elem_lst0 = geom.GetElemList();

    Vector<Real> X, Xn;
    elem_lst0.GetNodeCoord(&X, &Xn, nullptr);
    const Vector<Real> Ubg = bg_flow.Velocity(X);
    const Vector<Real> Tbg = dot_prod(bg_flow.Traction(X), Xn);
    elem_lst0.WriteVTK(std::string("vis/Ubg") + std::to_string(i), Ubg, comm); // Write VTK
    //elem_lst0.WriteVTK(std::string("vis/Tbg") + std::to_string(i), Ubg, comm); // Write VTK
    elem_lst0.Write(std::string("vis/geom") + std::to_string(i), comm);

    Vector<Real> U = stokes_mobility.ComputeVelocity(geom, Ubg, Tbg, tol, quad_eps);
    elem_lst0.WriteVTK(std::string("vis/U") + std::to_string(i), U, comm); // Write VTK
    geom.RigidBodyUpdate(U, dt);
    break;
    //exit(0);
  }
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    Profile::Enable(true);
    Comm comm = Comm::World();
    commandline_option_start(argc, argv, nullptr, comm);
    double quad_eps = strtod(commandline_option(argc, argv, "-qeps", "1e-10", false, nullptr, comm), nullptr);

    //ParticleFMM<double,3>::test(comm);
    test<double>(comm, 1e-13, quad_eps);

    Profile::print(&comm);
    commandline_option_end(argc, argv);
  }

  Comm::MPI_Finalize();
  return 0;
}


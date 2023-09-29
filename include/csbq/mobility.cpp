
namespace sctl {

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



  template <class Real> RigidBodyList<Real>::RigidBodyList(const Comm& comm, const Long Nobj, const Real loop_rad, const Geom& geom_type) : comm_(comm) {
    InitGeom(X, R, OrientVec, ChebOrder, FourierOrder, panel_len, Mr_lst, obj_elem_cnt, obj_elem_dsp, Nobj, loop_rad, geom_type);
    InitElemList(loc_elem_cnt, loc_elem_dsp, elem_lst, ChebOrder, FourierOrder, X, R, OrientVec, comm_);
    GetXc(Xc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);
  }

  template <class Real> template <class ValueType> void RigidBodyList<Real>::Init(const Vector<Long>& obj_elem_cnt_, const Vector<ValueType>& Mr_lst_, const Vector<Long>& ChebOrder_, const Vector<Long>& FourierOrder_, const Vector<ValueType>& panel_len_, const Vector<ValueType>& X_, const Vector<ValueType>& R_, const Vector<ValueType>& OrientVec_) {
    const Long Nobj = obj_elem_cnt_.Dim();
    obj_elem_cnt = obj_elem_cnt_;
    obj_elem_dsp.ReInit(Nobj); obj_elem_dsp = 0;
    omp_par::scan(obj_elem_cnt.begin(), obj_elem_dsp.begin(), Nobj);

    if (Mr_lst.Dim() != Mr_lst_.Dim()) Mr_lst.ReInit(Mr_lst_.Dim());
    for (Long i = 0; i < Mr_lst.Dim(); i++) Mr_lst[i] = (Real)Mr_lst_[i];
    SCTL_ASSERT(Mr_lst.Dim() == Nobj*COORD_DIM*COORD_DIM);

    ChebOrder = ChebOrder_;
    FourierOrder = FourierOrder_;

    if (panel_len.Dim() != panel_len_.Dim()) panel_len.ReInit(panel_len_.Dim());
    if (        X.Dim() !=         X_.Dim())         X.ReInit(        X_.Dim());
    if (        R.Dim() !=         R_.Dim())         R.ReInit(        R_.Dim());
    if (OrientVec.Dim() != OrientVec_.Dim()) OrientVec.ReInit(OrientVec_.Dim());
    for (Long i = 0; i < panel_len.Dim(); i++) panel_len[i] = (Real)panel_len_[i];
    for (Long i = 0; i <         X.Dim(); i++)         X[i] = (Real)        X_[i];
    for (Long i = 0; i <         R.Dim(); i++)         R[i] = (Real)        R_[i];
    for (Long i = 0; i < OrientVec.Dim(); i++) OrientVec[i] = (Real)OrientVec_[i];

    InitElemList<ValueType>(loc_elem_cnt, loc_elem_dsp, elem_lst, ChebOrder, FourierOrder, X_, R_, OrientVec_, comm_);
    GetXc(Xc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);
  }

  template <class Real> void RigidBodyList<Real>::Write(const std::string& fname) const {
    if (comm_.Rank()) return;
    FILE* f1 = fopen(fname.c_str(), "wb+");
    if (f1 == nullptr) {
      std::cout << "Unable to open file for writing:" << fname << '\n';
      return;
    }
    const uint64_t Nobj = obj_elem_cnt.Dim();
    const uint64_t Nelem = ChebOrder.Dim();
    const uint64_t Nnds = R.Dim();
    fwrite(&Nobj, sizeof(uint64_t), 1, f1);
    fwrite(&Nelem, sizeof(uint64_t), 1, f1);
    fwrite(&Nnds, sizeof(uint64_t), 1, f1);
    if (Nelem) {
      SCTL_ASSERT(obj_elem_cnt.Dim() == (Long)Nobj);
      SCTL_ASSERT(obj_elem_dsp.Dim() == (Long)Nobj);
      SCTL_ASSERT(Xc.Dim() == (Long)Nobj*COORD_DIM);
      SCTL_ASSERT(Mr_lst.Dim() == (Long)Nobj*COORD_DIM*COORD_DIM);

      SCTL_ASSERT(ChebOrder.Dim() == (Long)Nelem);
      SCTL_ASSERT(FourierOrder.Dim() == (Long)Nelem);
      SCTL_ASSERT(panel_len.Dim() == (Long)Nelem);

      SCTL_ASSERT(OrientVec.Dim() == (Long)Nnds*COORD_DIM);
      SCTL_ASSERT(X.Dim() == (Long)Nnds*COORD_DIM);
      SCTL_ASSERT(R.Dim() == (Long)Nnds);

      fwrite(&obj_elem_cnt[0], sizeof(Long), Nobj, f1);
      fwrite(&obj_elem_dsp[0], sizeof(Long), Nobj, f1);
      fwrite(&Xc[0], sizeof(Real), Nobj*COORD_DIM, f1);
      fwrite(&Mr_lst[0], sizeof(Real), Nobj*COORD_DIM*COORD_DIM, f1);

      fwrite(&ChebOrder[0], sizeof(Long), Nelem, f1);
      fwrite(&FourierOrder[0], sizeof(Long), Nelem, f1);
      fwrite(&panel_len[0], sizeof(Real), Nelem, f1);

      fwrite(&OrientVec[0], sizeof(Real), Nnds*COORD_DIM, f1);
      fwrite(&X[0], sizeof(Real), Nnds*COORD_DIM, f1);
      fwrite(&R[0], sizeof(Real), Nnds, f1);
    }
    fclose(f1);
  }

  template <class Real> void RigidBodyList<Real>::Read(const std::string& fname, const Long Nobj_) {
    FILE* f1 = fopen(fname.c_str(), "r");
    if (f1 == nullptr) {
      std::cout << "Unable to open file for reading:" << fname << '\n';
      return;
    }

    uint64_t Nobj, Nelem, Nnds, readlen;
    readlen = fread(&Nobj, sizeof(uint64_t), 1, f1); assert(readlen == 1);
    readlen = fread(&Nelem, sizeof(uint64_t), 1, f1); assert(readlen == 1);
    readlen = fread(&Nnds, sizeof(uint64_t), 1, f1); assert(readlen == 1);

    obj_elem_cnt.ReInit(Nobj);
    obj_elem_dsp.ReInit(Nobj);
    Xc.ReInit(Nobj*COORD_DIM);
    Mr_lst.ReInit(Nobj*COORD_DIM*COORD_DIM);

    ChebOrder.ReInit(Nelem);
    FourierOrder.ReInit(Nelem);
    panel_len.ReInit(Nelem);

    OrientVec.ReInit(Nnds*COORD_DIM);
    X.ReInit(Nnds*COORD_DIM);
    R.ReInit(Nnds);

    readlen = fread(&obj_elem_cnt[0], sizeof(Long), Nobj, f1); SCTL_ASSERT(readlen == Nobj);
    readlen = fread(&obj_elem_dsp[0], sizeof(Long), Nobj, f1); SCTL_ASSERT(readlen == Nobj);
    readlen = fread(&Xc[0], sizeof(Real), Nobj*COORD_DIM, f1); SCTL_ASSERT(readlen == Nobj*COORD_DIM);
    readlen = fread(&Mr_lst[0], sizeof(Real), Nobj*COORD_DIM*COORD_DIM, f1); SCTL_ASSERT(readlen == Nobj*COORD_DIM*COORD_DIM);

    readlen = fread(&ChebOrder[0], sizeof(Long), Nelem, f1); SCTL_ASSERT(readlen == Nelem);
    readlen = fread(&FourierOrder[0], sizeof(Long), Nelem, f1); SCTL_ASSERT(readlen == Nelem);
    readlen = fread(&panel_len[0], sizeof(Real), Nelem, f1); SCTL_ASSERT(readlen == Nelem);

    readlen = fread(&OrientVec[0], sizeof(Real), Nnds*COORD_DIM, f1); SCTL_ASSERT(readlen == Nnds*COORD_DIM);
    readlen = fread(&X[0], sizeof(Real), Nnds*COORD_DIM, f1); SCTL_ASSERT(readlen == Nnds*COORD_DIM);
    readlen = fread(&R[0], sizeof(Real), Nnds, f1); SCTL_ASSERT(readlen == Nnds);

    SCTL_UNUSED(readlen);
    fclose(f1);

    if (Nobj_ >= 0) { // Select the first Nobj_ objects
      Vector<Long> obj_elem_cnt_, obj_elem_dsp_;
      Vector<Real> Xc_, Mr_lst_;
      Vector<Long> ChebOrder_, FourierOrder_;
      Vector<Real> panel_len_;
      Vector<Real> X_, R_, OrientVec_;

      Long Nelem_ = 0, Nnds_ = 0;
      for (Long i = 0; i < Nobj_; i++) {
        Nelem_ += obj_elem_cnt[i];
        obj_elem_cnt_.PushBack(obj_elem_cnt[i]);
        obj_elem_dsp_.PushBack(obj_elem_dsp[i]);
        for (Long j = 0; j < COORD_DIM; j++) Xc_.PushBack(Xc[i*COORD_DIM+j]);
        for (Long j = 0; j < COORD_DIM*COORD_DIM; j++) Mr_lst_.PushBack(Mr_lst[i*COORD_DIM*COORD_DIM+j]);
      }
      for (Long i = 0; i < Nelem_; i++) {
        Nnds_ += ChebOrder[i];
        ChebOrder_.PushBack(ChebOrder[i]);
        FourierOrder_.PushBack(FourierOrder[i]);
        panel_len_.PushBack(panel_len[i]);
      }
      for (Long i = 0; i < Nnds_; i++) R_.PushBack(R[i]);
      for (Long i = 0; i < Nnds_*COORD_DIM; i++) X_.PushBack(X[i]);
      for (Long i = 0; i < Nnds_*COORD_DIM; i++) OrientVec_.PushBack(OrientVec[i]);

      obj_elem_cnt = obj_elem_cnt_;
      obj_elem_dsp = obj_elem_dsp_;
      Xc = Xc_;
      Mr_lst = Mr_lst_;
      ChebOrder = ChebOrder_;
      FourierOrder = FourierOrder_;
      panel_len = panel_len_;
      R = R_;
      X = X_;
      OrientVec = OrientVec_;
    }

    InitElemList(loc_elem_cnt, loc_elem_dsp, elem_lst, ChebOrder, FourierOrder, X, R, OrientVec, comm_);
  }

  template <class Real> const SlenderElemList<Real>& RigidBodyList<Real>::GetElemList() const {
    return elem_lst;
  }

  template <class Real> void RigidBodyList<Real>::GetObjPosition(Vector<Real>* Xc_, Vector<Real>* Mr_) const {
    const Long Nobj = obj_elem_cnt.Dim();
    const Long a = (comm_.Rank()+0)*Nobj / comm_.Size();
    const Long b = (comm_.Rank()+1)*Nobj / comm_.Size();

    if (Xc_) Xc_[0] = Vector<Real>((b-a)*COORD_DIM, (Iterator<Real>)Xc.begin() + a*COORD_DIM, false);
    if (Mr_) Mr_[0] = Vector<Real>((b-a)*COORD_DIM*COORD_DIM, (Iterator<Real>)Mr_lst.begin() + a*COORD_DIM*COORD_DIM, false);
  }

  template <class Real> void RigidBodyList<Real>::RigidBodyMotionProj(Vector<Real>& u_proj, const Vector<Real>& u) const { // TODO: optimize
    Profile::Tic("RigidBodyProj", &comm_);
    auto inner_prod = [this](const Vector<Real>& A, const Vector<Real>& B) {
      SCTL_ASSERT(A.Dim() == B.Dim());
      Vector<Real> AB(A.Dim()), IntegralAB;
      for (Long i = 0; i < A.Dim(); i++) AB[i] = A[i]*B[i];
      ObjIntegral(IntegralAB, AB, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);

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
      GetNodeObjIdx(node_obj_idx, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);
    }
    const Long Nloc = (elem_node_cnt.Dim() ? elem_node_cnt[elem_node_cnt.Dim()-1] + elem_node_dsp[elem_node_dsp.Dim()-1] : 0);
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
      RigidBodyVelocity(q_k, Uc, Tc, Xc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);

      for (Integer j = 0; j < k; j++) { // q_k <-- orthogonalize(q_k, Q)
        Vector<Real> q_j(Q.Dim(1), Q[j], false);
        Vector<Real> qk_dot_qj = inner_prod(q_k, q_j);
        for (Long node_idx = 0; node_idx < Nloc; node_idx++) {
          const Long obj = node_obj_idx[node_idx];
          for(Long l = 0; l < COORD_DIM; l++) {
            q_k[node_idx*COORD_DIM+l] -= q_j[node_idx*COORD_DIM+l] * qk_dot_qj[obj];
          }
        }
      }

      Vector<Real> rsqrt_qk_dot_qk = inner_prod(q_k, q_k);
      for (auto& x : rsqrt_qk_dot_qk) x = 1/sqrt<Real>(x);

      for (Long node_idx = 0; node_idx < Nloc; node_idx++) { // Q[k] <-- normalize(q_k)
        const Long obj = node_obj_idx[node_idx];
        for(Long l = 0; l < COORD_DIM; l++) {
          q_k[node_idx*COORD_DIM+l] *= rsqrt_qk_dot_qk[obj];
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
    Profile::Toc();
  }

  template <class Real> void RigidBodyList<Real>::GetRigidBodyMotion(Vector<Real>& Uc_loc, Vector<Real>& Omega_c_loc, const Vector<Real>& U_loc) const {
    auto cross_prod = [](StaticArray<Real,COORD_DIM>& AxB, const StaticArray<Real,COORD_DIM>& A, const StaticArray<Real,COORD_DIM>& B) {
      AxB[0] = A[1]*B[2] - A[2]*B[1];
      AxB[1] = A[2]*B[0] - A[0]*B[2];
      AxB[2] = A[0]*B[1] - A[1]*B[0];
    };

    Vector<Real> X_loc;
    Vector<Long> node_cnt_loc, node_dsp_loc;
    { // Set X_loc, node_cnt_loc, node_dsp_loc
      elem_lst.GetNodeCoord(&X_loc, nullptr, &node_cnt_loc);
      node_dsp_loc.ReInit(node_cnt_loc.Dim());
      if (node_cnt_loc.Dim()) node_dsp_loc[0] = 0;
      omp_par::scan(node_cnt_loc.begin(), node_dsp_loc.begin(), node_cnt_loc.Dim());
    }

    Vector<Real> Uc, Omega_c;
    { // Set Omega_c, Uc
      Vector<Real> area;
      ObjIntegral(Uc, U_loc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);
      ObjIntegral(area, X_loc*0+1, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);
      Uc /= area;
    }
    { // Set Omega_c
      const Long Nobj = obj_elem_cnt.Dim();
      const Long N = U_loc.Dim()/COORD_DIM;
      Vector<Real> UUt, UUt_loc(N*COORD_DIM*COORD_DIM);
      for (Long i = 0; i < loc_elem_cnt; i++) {
        const Long elem_idx = loc_elem_dsp + i;
        const Long obj_idx = std::upper_bound(obj_elem_dsp.begin(), obj_elem_dsp.end(), elem_idx) - obj_elem_dsp.begin() - 1;
        for (Long j = 0; j < node_cnt_loc[i]; j++) {
          const Long node_idx = node_dsp_loc[i]+j;
          for (Integer k0 = 0; k0 < COORD_DIM; k0++) {
            for (Integer k1 = 0; k1 < COORD_DIM; k1++) {
              const Real Ur_k0 = U_loc[node_idx*COORD_DIM+k0] - Uc[obj_idx*COORD_DIM+k0];
              const Real Ur_k1 = U_loc[node_idx*COORD_DIM+k1] - Uc[obj_idx*COORD_DIM+k1];
              UUt_loc[(node_idx*COORD_DIM+k0)*COORD_DIM+k1] = Ur_k0 * Ur_k1;
            }
          }
        }
      }
      ObjIntegral(UUt, UUt_loc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);

      Omega_c.ReInit(Nobj*COORD_DIM);
      for (Long obj = 0; obj < Nobj; obj++) { // Compute Omega_c direction
        Matrix<Real> M(3,3, UUt.begin()+obj*COORD_DIM*COORD_DIM);
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
          Omega_c[obj*COORD_DIM+k] = U[k][min_idx];
        }
      }
      { // Compute Omega_c magnitude
        Vector<Real> V_dot_U_loc(N), V_dot_V_loc(N);
        V_dot_U_loc = 0;
        V_dot_V_loc = 0;
        for (Long i = 0; i < loc_elem_cnt; i++) {
          const Long elem_idx = loc_elem_dsp + i;
          const Long obj_idx = std::upper_bound(obj_elem_dsp.begin(), obj_elem_dsp.end(), elem_idx) - obj_elem_dsp.begin() - 1;
          for (Long j = 0; j < node_cnt_loc[i]; j++) {
            const Long node_idx = node_dsp_loc[i]+j;

            StaticArray<Real,COORD_DIM> Omega_orient, R, V;
            for (Long k = 0; k < COORD_DIM; k++) {
              Omega_orient[k] = Omega_c[obj_idx*COORD_DIM+k];
              R[k] = X_loc[node_idx*COORD_DIM+k] - Xc[obj_idx*COORD_DIM+k];
            }
            cross_prod(V, Omega_orient, R);

            for (Long k = 0; k < COORD_DIM; k++) {
              const Real Ur_k = U_loc[node_idx*COORD_DIM+k] - Uc[obj_idx*COORD_DIM+k];
              V_dot_U_loc[node_idx] += Ur_k * V[k];
              V_dot_V_loc[node_idx] += V[k] * V[k];
            }
          }
        }

        Vector<Real> V_dot_U, V_dot_V;
        ObjIntegral(V_dot_V, V_dot_V_loc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);
        ObjIntegral(V_dot_U, V_dot_U_loc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm_);
        for (Long obj_idx = 0; obj_idx < Nobj; obj_idx++) {
          for (Long k = 0; k < COORD_DIM; k++) {
            Omega_c[obj_idx*COORD_DIM+k] *= V_dot_U[obj_idx]/V_dot_V[obj_idx];
          }
        }
      }
    }

    { // Set Uc_loc, Omega_c_loc
      const Long Nobj = obj_elem_cnt.Dim();
      const Long a = (comm_.Rank()+0)*Nobj / comm_.Size();
      const Long b = (comm_.Rank()+1)*Nobj / comm_.Size();
      Uc_loc = Vector<Real>((b-a)*COORD_DIM, Uc.begin() + a*COORD_DIM, false);
      Omega_c_loc = Vector<Real>((b-a)*COORD_DIM, Omega_c.begin() + a*COORD_DIM, false);
    }
  }

  template <class Real> void RigidBodyList<Real>::RigidBodyUpdate(const Vector<Real>& dXc_loc, const Vector<Real>& Mr_loc) {
    const Vector<Real> dXc = Allgather(dXc_loc, comm_);
    const Vector<Real> Mr_lst0 = Allgather(Mr_loc, comm_);

    const Long Nobj = obj_elem_cnt.Dim();
    Vector<Long> cnt(Nobj), dsp(Nobj); cnt = 0; dsp = 0;
    for (Long obj = 0; obj < Nobj; obj++) {
      for (Long i = 0; i < obj_elem_cnt[obj]; i++) {
        cnt[obj] += ChebOrder[obj_elem_dsp[obj] + i];
      }
    }
    omp_par::scan(cnt.begin(), dsp.begin(), Nobj);

    for (Long obj = 0; obj < Nobj; obj++) { // Updatee X, OrientVec
      const Matrix<Real> Mr0(COORD_DIM, COORD_DIM, (Iterator<Real>)Mr_lst0.begin() + obj*COORD_DIM*COORD_DIM, false);
      Matrix<Real> X_(cnt[obj], COORD_DIM, X.begin() + dsp[obj]*COORD_DIM, false);
      Matrix<Real> OrientVec_(cnt[obj], COORD_DIM, OrientVec.begin() + dsp[obj]*COORD_DIM, false);

      for (Long j = 0; j < X_.Dim(0); j++) {
        for (Long k = 0; k < COORD_DIM; k++) {
          X_[j][k] -= Xc[obj*COORD_DIM+k];
        }
      }
      X_ = X_ * Mr0;;
      OrientVec_ = OrientVec_ * Mr0;;
      for (Long j = 0; j < X_.Dim(0); j++) {
        for (Long k = 0; k < COORD_DIM; k++) {
          X_[j][k] += Xc[obj*COORD_DIM+k] + dXc[obj*COORD_DIM+k];
        }
      }
    }
    InitElemList(loc_elem_cnt, loc_elem_dsp, elem_lst, ChebOrder, FourierOrder, X, R, OrientVec, comm_);

    for (Long obj = 0; obj < Nobj; obj++) { // Update Mr
      const Tensor<Real,false,COORD_DIM,COORD_DIM> Mr0((Iterator<Real>)Mr_lst0.begin() + obj*COORD_DIM*COORD_DIM);
      Tensor<Real,false,COORD_DIM,COORD_DIM> Mr(Mr_lst.begin() + obj*COORD_DIM*COORD_DIM);
      Mr = Mr * Mr0;
    }
    Xc += dXc;
  }

  template <class Real> void RigidBodyList<Real>::RefineAdaptive(const Vector<Vector<Real>>& v_, const Real tol, const bool relative) {
    // TODO: implements element-wise relative error. Would global relative error be better?
    static constexpr Integer FOURIER_ORDER_MIN = 4;
    static constexpr Integer FOURIER_ORDER_STEP = 4;
    static constexpr Integer FOURIER_ORDER_MAX = 100;

    Vector<Real> X_loc;
    Vector<Long> node_cnt_loc, node_dsp_loc;
    { // Set X_loc, node_cnt_loc, node_dsp_loc
      elem_lst.GetNodeCoord(&X_loc, nullptr, &node_cnt_loc);
      node_dsp_loc.ReInit(node_cnt_loc.Dim());
      if (node_cnt_loc.Dim()) node_dsp_loc[0] = 0;
      omp_par::scan(node_cnt_loc.begin(), node_dsp_loc.begin(), node_cnt_loc.Dim());
    }
    Vector<Vector<Real>> v = v_;
    v.PushBack(X_loc);

    Vector<Long> FourierOrderLoc(loc_elem_cnt);
    Vector<Real> elem_fourier_err_loc(loc_elem_cnt), elem_cheb_err_loc(loc_elem_cnt);
    for (Long i = 0; i < loc_elem_cnt; i++) {
      const Long elem_idx = loc_elem_dsp + i;
      const Long ChebOrder_ = ChebOrder[elem_idx];
      const Long FourierOrder_ = FourierOrder[elem_idx];
      SCTL_ASSERT(node_cnt_loc[i] == ChebOrder_ * FourierOrder_);

      Real max_cheb_err = 0, max_fourier_err = 0;
      for (Long j = 0; j < v.Dim(); j++) {
        const Long dof = v[j].Dim()/(X_loc.Dim()/COORD_DIM);

        Real max_val_inv = 0;
        if (relative) {
          Real max_val = 0;
          const Vector<Real> vv(ChebOrder_*FourierOrder_*dof, v[j].begin() + node_dsp_loc[i]*dof, false);
          for (const auto& x : vv) max_val = std::max<Real>(max_val, fabs(x));
          max_val_inv = 1/max_val;
        } else {
          max_val_inv = 1;
        }

        const Matrix<Real> v_(ChebOrder_, FourierOrder_*dof, v[j].begin() + node_dsp_loc[i]*dof, false);
        const Matrix<Real> cheb_err = v_.Transpose() * ChebErrMatrix(ChebOrder_);
        for (const auto& x : cheb_err) max_cheb_err = std::max<Real>(max_cheb_err, fabs(x)*max_val_inv);

        const Matrix<Real> v_tmp = Matrix<Real>(ChebOrder_*FourierOrder_, dof, (Iterator<Real>)v_.begin(), false).Transpose();
        const Matrix<Real> v__(dof*ChebOrder_, FourierOrder_, (Iterator<Real>)v_tmp.begin(), false);
        const Matrix<Real> fourier_err = v__ * FourierErrMatrix(FourierOrder_);
        for (const auto& x : fourier_err) max_fourier_err = std::max<Real>(max_fourier_err, fabs(x)*max_val_inv);
      }
      elem_fourier_err_loc[i] = pow<Real>(max_fourier_err, (FourierOrder_+1)/(Real)(FourierOrder_-1));
      elem_cheb_err_loc[i] = pow<Real>(max_cheb_err, (ChebOrder_+1)/(Real)(ChebOrder_-1));

      FourierOrderLoc[i] = std::min<Long>(std::max<Long>(FOURIER_ORDER_MIN, (Long)(FourierOrder_*log<Real>(tol)/log<Real>(max_fourier_err)+FOURIER_ORDER_STEP)), FOURIER_ORDER_MAX);
      FourierOrderLoc[i] -= FourierOrderLoc[i] % FOURIER_ORDER_STEP;
    }
    FourierOrder = Allgather(FourierOrderLoc, comm_);

    Vector<Long> ChebOrder_new, FourierOrder_new;
    Vector<Long> obj_elem_cnt_new, obj_elem_dsp_new;
    Vector<Real> X_new, R_new, OrientVec_new, panel_len_new;
    { // refine/coarsen panels
      Vector<Long> node_cnt = ChebOrder, node_dsp = ChebOrder*0;
      omp_par::scan(node_cnt.begin(), node_dsp.begin(), node_cnt.Dim());

      Vector<Real> cheb_err = Allgather(elem_cheb_err_loc, comm_);
      for (Long obj = 0; obj < obj_elem_cnt.Dim(); obj++) {
        const Long Nelem = obj_elem_cnt[obj];
        Vector<Real> panel_length(Nelem, panel_len.begin() + obj_elem_dsp[obj], false);
        Vector<Real> panel_offset(Nelem); panel_offset = 0;
        omp_par::scan(panel_length.begin(), panel_offset.begin(), Nelem);

        Long obj_elem_cnt_new_ = 0;
        for (Long i = 0; i < Nelem; i++) {
          const Long elem = obj_elem_dsp[obj] + i;
          const Long FourierOrder_ = FourierOrder[elem];
          const Long ChebOrder_ = ChebOrder[elem];

          const auto RefinementMatrix = [](const Integer ChebOrder_) {
            Vector<Real> trg_nds(2*ChebOrder_);
            const Vector<Real>& src_nds = ChebQuadRule<Real>::nds(ChebOrder_);
            for (Long j = 0; j < ChebOrder_; j++) {
              trg_nds[j] = src_nds[j]*0.5;
              trg_nds[ChebOrder_+j] = 0.5 + src_nds[j]*0.5;
            }

            Vector<Real> M_(ChebOrder_ * 2*ChebOrder_);
            LagrangeInterp<Real>::Interpolate(M_, src_nds, trg_nds);
            return Matrix<Real>(ChebOrder_, 2*ChebOrder_, M_.begin(), false).Transpose();
          };
          const auto CoarseningMatrix = [](const Integer ChebOrder_) {
            const Vector<Real>& src_nds = ChebQuadRule<Real>::nds(ChebOrder_);
            Vector<Real> trg_nds0(ChebOrder_/2), trg_nds1(ChebOrder_-ChebOrder_/2);
            for (Long j = 0; j < ChebOrder_/2; j++) trg_nds0[j] = src_nds[j]*2;
            for (Long j = ChebOrder_/2; j < ChebOrder_; j++) trg_nds1[j-ChebOrder_/2] = src_nds[j]*2-1.0;

            Matrix<Real> Mcoarsen(ChebOrder_, 2*ChebOrder_); Mcoarsen = 0;
            Vector<Real> M0(ChebOrder_ * trg_nds0.Dim()), M1(ChebOrder_ * trg_nds1.Dim());
            LagrangeInterp<Real>::Interpolate(M0, src_nds, trg_nds0);
            LagrangeInterp<Real>::Interpolate(M1, src_nds, trg_nds1);

            for (Long k0 = 0; k0 < ChebOrder_; k0++) {
              for (Long k1 = 0; k1 < trg_nds0.Dim(); k1++) {
                Mcoarsen[k1][k0] = M0[k0*trg_nds0.Dim()+k1];
              }
            }
            for (Long k0 = 0; k0 < ChebOrder_; k0++) {
              for (Long k1 = 0; k1 < trg_nds1.Dim(); k1++) {
                Mcoarsen[trg_nds0.Dim()+k1][ChebOrder_+k0] = M1[k0*trg_nds1.Dim()+k1];
              }
            }
            return Mcoarsen;
          };

          bool if_coarsen = false;
          { // Set if_coarsen
            const Real tol_coarsen = 0.1 * tol * pow<Real>(0.5, ChebOrder_);
            if (i<Nelem-1 && panel_length[i]==panel_length[i+1] && ChebOrder[elem]==ChebOrder[elem+1]) {
              if (((Integer)(panel_offset[i]/panel_length[i]+0.5)) % 2 == 0) { // panels are siblings in binary tree
                if (cheb_err[elem]<tol_coarsen && cheb_err[elem+1]<tol_coarsen) {
                  if_coarsen = true;
                }
              }
            }
          }

          if (if_coarsen) {
            const Matrix<Real> Mcoarsen = CoarseningMatrix(ChebOrder_);
            panel_len_new.PushBack(panel_len[elem]*2);
            FourierOrder_new.PushBack(std::max<Long>(FourierOrder[elem],FourierOrder[elem+1]));
            ChebOrder_new.PushBack(ChebOrder_);
            { // Set X_new, R_new, OrientVec_new
              Matrix<Real> R_ = Mcoarsen * Matrix<Real>(2*ChebOrder_, 1, R.begin() + node_dsp[elem], false);
              Matrix<Real> X_ = Mcoarsen * Matrix<Real>(2*ChebOrder_, COORD_DIM, X.begin() + node_dsp[elem]*COORD_DIM, false);
              Matrix<Real> OrientVec_ = Mcoarsen * Matrix<Real>(2*ChebOrder_, COORD_DIM, OrientVec.begin() + node_dsp[elem]*COORD_DIM, false);
              for (Long j = 0; j < ChebOrder_; j++) {
                R_new.PushBack(R_[j][0]);
                for (Integer k = 0; k < COORD_DIM; k++) {
                  X_new.PushBack(X_[j][k]);
                  OrientVec_new.PushBack(OrientVec_[j][k]);
                }
              }
            }
            obj_elem_cnt_new_++;
            i++;
          } else if (cheb_err[elem] > tol) {
            const Matrix<Real> M = RefinementMatrix(ChebOrder_);
            panel_len_new.PushBack(panel_len[elem]/2);
            panel_len_new.PushBack(panel_len[elem]/2);
            FourierOrder_new.PushBack(FourierOrder_);
            FourierOrder_new.PushBack(FourierOrder_);
            ChebOrder_new.PushBack(ChebOrder_);
            ChebOrder_new.PushBack(ChebOrder_);
            { // Set X_new, R_new, OrientVec_new
              Matrix<Real> R_ = M * Matrix<Real>(ChebOrder_, 1, R.begin() + node_dsp[elem], false);
              Matrix<Real> X_ = M * Matrix<Real>(ChebOrder_, COORD_DIM, X.begin() + node_dsp[elem]*COORD_DIM, false);
              Matrix<Real> OrientVec_ = M * Matrix<Real>(ChebOrder_, COORD_DIM, OrientVec.begin() + node_dsp[elem]*COORD_DIM, false);
              for (Long j = 0; j < 2*ChebOrder_; j++) {
                R_new.PushBack(R_[j][0]);
                for (Integer k = 0; k < COORD_DIM; k++) {
                  X_new.PushBack(X_[j][k]);
                  OrientVec_new.PushBack(OrientVec_[j][k]);
                }
              }
            }
            obj_elem_cnt_new_ += 2;
          } else {
            panel_len_new.PushBack(panel_len[elem]);
            FourierOrder_new.PushBack(FourierOrder_);
            ChebOrder_new.PushBack(ChebOrder_);
            for (Long j = 0; j < ChebOrder_; j++) {
              const Long node_idx = node_dsp[elem] + j;
              R_new.PushBack(R[node_idx]);
              for (Integer k = 0; k < COORD_DIM; k++) {
                X_new.PushBack(X[node_idx*COORD_DIM+k]);
                OrientVec_new.PushBack(OrientVec[node_idx*COORD_DIM+k]);
              }
            }
            obj_elem_cnt_new_++;
          }
        }
        obj_elem_cnt_new.PushBack(obj_elem_cnt_new_);
      }
    }

    X = X_new;
    R = R_new;
    OrientVec = OrientVec_new;

    panel_len = panel_len_new;
    ChebOrder = ChebOrder_new;
    FourierOrder = FourierOrder_new;

    obj_elem_cnt = obj_elem_cnt_new;
    obj_elem_dsp.ReInit(obj_elem_cnt.Dim()); obj_elem_dsp = 0;
    omp_par::scan(obj_elem_cnt.begin(), obj_elem_dsp.begin(), obj_elem_cnt.Dim());

    InitElemList(loc_elem_cnt, loc_elem_dsp, elem_lst, ChebOrder, FourierOrder, X, R, OrientVec, comm_);

    //if (!comm_.Rank()) std::cout<<FourierOrderNew<<'\n';
    //std::cout<<elem_fourier_err<<'\n';
    //std::cout<<elem_cheb_err<<'\n';
  }

  template <class Real> void RigidBodyList<Real>::RotationMatrix(Vector<Real>& Mr_lst, const Vector<Real>& Omega, const Real dt) {
    const Long Nobj = Omega.Dim() / COORD_DIM;
    if (Mr_lst.Dim() != Nobj*COORD_DIM*COORD_DIM) Mr_lst.ReInit(Nobj*COORD_DIM*COORD_DIM);
    for (Long obj = 0; obj < Nobj; obj++) {
      Tensor<Real,true,1,COORD_DIM> W0, W1, W2;
      Tensor<Real,true,COORD_DIM,COORD_DIM> P1, P2, R;
      for (Long k = 0; k < COORD_DIM; k++) W0(0,k) = Omega[obj*COORD_DIM+k]*dt;

      Real theta = atan2(W0(0,1), W0(0,0));
      P1(0,0) = cos<Real>(theta); P1(0,1) =-sin<Real>(theta); P1(0,2) = 0;
      P1(1,0) = sin<Real>(theta); P1(1,1) = cos<Real>(theta); P1(1,2) = 0;
      P1(2,0) =                0; P1(2,1) =                0; P1(2,2) = 1;
      W1 = W0 * P1;

      Real phi = atan2(W1(0,2), W1(0,0));
      P2(0,0) = cos<Real>(phi); P2(0,1) = 0; P2(0,2) =-sin<Real>(phi);
      P2(1,0) =              0; P2(1,1) = 1; P2(1,2) =              0;
      P2(2,0) = sin<Real>(phi); P2(2,1) = 0; P2(2,2) = cos<Real>(phi);
      W2 = W1 * P2;

      Real w = W2(0,0);
      R(0,0) = 1; R(0,1) =            0; R(0,2) =            0;
      R(1,0) = 0; R(1,1) = cos<Real>(w); R(1,2) = sin<Real>(w);
      R(2,0) = 0; R(2,1) =-sin<Real>(w); R(2,2) = cos<Real>(w);

      Tensor<Real,false,COORD_DIM,COORD_DIM> Mr(Mr_lst.begin() + obj*COORD_DIM*COORD_DIM);
      Mr = P1 * P2 * R * P2.RotateRight() * P1.RotateRight();
    }
  }

  template <class Real> void RigidBodyList<Real>::ApplyPrecond(Vector<Real>& Ax0, const Vector<Real> x0, const Matrix<Real>& A, const RigidBodyList<Real>& ref_geom) const {
    constexpr Long dof = COORD_DIM;
    Long Nelem, Nobj;

    // for each object
    Vector<Real> Mr_lst_;
    Vector<Long> obj_elem_cnt_, obj_elem_dsp_;

    // for each element
    Vector<Real> panel_len_;
    Vector<Long> FourierOrder_, ChebOrder_;
    Vector<Long> node_cnt_, node_dsp_;
    { // Set
      const Long a = (comm_.Rank()+0)*obj_elem_cnt.Dim() / comm_.Size();
      const Long b = (comm_.Rank()+1)*obj_elem_cnt.Dim() / comm_.Size();
      Nelem = (b>0 ? obj_elem_cnt[b-1] + obj_elem_dsp[b-1] - obj_elem_dsp[a] : 0);
      Nobj = b-a;

      obj_elem_dsp_.ReInit(Nobj); obj_elem_dsp_ = 0;
      obj_elem_cnt_.ReInit(Nobj, (Iterator<Long>)obj_elem_cnt.begin() + a);
      omp_par::scan(obj_elem_cnt_.begin(), obj_elem_dsp_.begin(), obj_elem_cnt_.Dim());

      Mr_lst_.ReInit(Nobj*COORD_DIM*COORD_DIM, (Iterator<Real>)Mr_lst.begin() + a*COORD_DIM*COORD_DIM);
      FourierOrder_.ReInit(Nelem, (Iterator<Long>)FourierOrder.begin() + obj_elem_dsp[a]);
      ChebOrder_.ReInit(Nelem, (Iterator<Long>)ChebOrder.begin() + obj_elem_dsp[a]);
      panel_len_.ReInit(Nelem, (Iterator<Real>)panel_len.begin() + obj_elem_dsp[a]);

      node_cnt_.ReInit(Nelem);
      node_dsp_.ReInit(Nelem); node_dsp_ = 0;
      for (Long i = 0; i < Nelem; i++) node_cnt_[i] = ChebOrder_[i] * FourierOrder_[i];
      omp_par::scan(node_cnt_.begin(), node_dsp_.begin(), node_cnt_.Dim());
    }

    Vector<Real> x1 = x0;
    { // Parition x1 across processors
      const Long Nnds = (Nelem ? node_cnt_[Nelem-1] + node_dsp_[Nelem-1] : 0);
      comm_.PartitionN(x1, Nnds);
      SCTL_ASSERT(Nnds==0 || x1.Dim()==Nnds*dof);
    }

    Vector<Real> x2;
    for (Long obj = 0; obj < Nobj; obj++) { // x2 <-- resample(x1)
      const Long elem_offset = obj_elem_dsp_[obj];
      const Long elem_count = obj_elem_cnt_[obj];
      SCTL_ASSERT_MSG(elem_count, "Each object must have at least one element");

      const Long nds_offset = node_dsp_[elem_offset];
      const Long nds_count = node_dsp_[elem_offset+elem_count-1] + node_cnt_[elem_offset+elem_count-1] - nds_offset;

      const Vector<Long> ChebOrder_src(elem_count, ChebOrder_.begin() + elem_offset, false);
      const Vector<Long> FourierOrder_src(elem_count, FourierOrder_.begin() + elem_offset, false);
      const Vector<Real> panel_len_src(elem_count, panel_len_.begin() + elem_offset, false);
      const Vector<Real> x_src(nds_count*dof, x1.begin() + nds_offset*dof, false);

      Vector<Real> x_trg;
      Resample<dof>(x_trg, ref_geom.panel_len, ref_geom.ChebOrder, ref_geom.FourierOrder, x_src, panel_len_src, ChebOrder_src, FourierOrder_src);
      for (Long i = 0; i < x_trg.Dim(); i++) x2.PushBack(x_trg[i]);
    }

    for (Long obj = 0; obj < Nobj; obj++) { // x2 <-- rotate(x2)
      const Long N = x2.Dim()/Nobj/dof;
      const Matrix<Real> Mr = Matrix<Real>(COORD_DIM, COORD_DIM, Mr_lst_.begin()+obj*COORD_DIM*COORD_DIM).pinv(sctl::machine_eps<Real>());
      Matrix<Real> U(N, dof, x2.begin() + obj*N*dof, false);
      U = U * Mr;
    }

    Vector<Real> Ax1(Nobj * A.Dim(1));
    if (A.Dim(0) * A.Dim(1) == 0) { //  Set Ax1 <-- x2 * A
      Ax1 = x2;
    } else {
      const Matrix<Real> x_(Nobj, A.Dim(0), x2.begin(), false);
      Matrix<Real> Ax_(Nobj, A.Dim(1), Ax1.begin(), false);
      Ax_ = x_ * A;
    }

    for (Long obj = 0; obj < Nobj; obj++) { // Ax1 <-- rotate(Ax1)
      const Long N = Ax1.Dim()/Nobj/dof;
      const Matrix<Real> Mr = Matrix<Real>(COORD_DIM, COORD_DIM, Mr_lst_.begin()+obj*COORD_DIM*COORD_DIM);
      Matrix<Real> U(N, dof, Ax1.begin() + obj*N*dof, false);
      U = U * Mr;
    }

    Vector<Real> Ax2;
    for (Long obj = 0; obj < Nobj; obj++) { // Ax2 <-- resample(Ax1)
      const Long N = Ax1.Dim()/Nobj;
      const Vector<Real> x_src(N, Ax1.begin() + obj*N, false);

      const Long elem_offset = obj_elem_dsp_[obj];
      const Long elem_count = obj_elem_cnt_[obj];
      SCTL_ASSERT_MSG(elem_count, "Each object must have at least one element");

      //const Long nds_offset = node_dsp_[elem_offset];
      //const Long nds_count = node_dsp_[elem_offset+elem_count-1] + node_cnt_[elem_offset+elem_count-1];

      const Vector<Long> ChebOrder_trg(elem_count, ChebOrder_.begin() + elem_offset, false);
      const Vector<Long> FourierOrder_trg(elem_count, FourierOrder_.begin() + elem_offset, false);
      const Vector<Real> panel_len_trg(elem_count, panel_len_.begin() + elem_offset, false);

      Vector<Real> x_trg;
      Resample<dof>(x_trg, panel_len_trg, ChebOrder_trg, FourierOrder_trg, x_src, ref_geom.panel_len, ref_geom.ChebOrder, ref_geom.FourierOrder);
      for (Long i = 0; i < x_trg.Dim(); i++) Ax2.PushBack(x_trg[i]);
    }

    comm_.PartitionN(Ax2, x0.Dim()); // Parition Ax2 across processors
    Ax0 = Ax2;
  }

  template <class Real> Real RigidBodyList<Real>::GetMinDist() const {
    const Long Nobj = obj_elem_cnt.Dim();
    Vector<Long> obj_node_cnt(Nobj), obj_node_dsp(Nobj);
    obj_node_cnt = 0; obj_node_dsp = 0;
    for (Long i = 0; i < Nobj; i++) {
      for (Long j = 0; j < obj_elem_cnt[i]; j++) {
        obj_node_cnt[i] += ChebOrder[obj_elem_dsp[i]+j];
      }
    }
    omp_par::scan(obj_node_cnt.begin(), obj_node_dsp.begin(), obj_node_cnt.Dim());

    Real r_min = -1;
    for (Long i = 0; i < Nobj; i++) {
      for (Long j = 0; j < Nobj; j++) {
        if (i == j) continue;
        for (Long ii = 0; ii < obj_node_cnt[i]; ii++) {
          for (Long jj = 0; jj < obj_node_cnt[j]; jj++) {
            const Long idx_i = obj_node_dsp[i]+ii;
            const Long idx_j = obj_node_dsp[j]+jj;
            const StaticArray<Real,COORD_DIM> dx{X[idx_i*COORD_DIM+0]-X[idx_j*COORD_DIM+0],
                                                 X[idx_i*COORD_DIM+1]-X[idx_j*COORD_DIM+1],
                                                 X[idx_i*COORD_DIM+2]-X[idx_j*COORD_DIM+2]};
            const Real r = sqrt<Real>(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) - R[idx_i] - R[idx_j];
            r_min = (r_min < 0 ? r : std::min<Real>(r_min, r));
          }
        }
      }
    }
    return r_min;
  }

  template <class Real> template <Integer dof> void RigidBodyList<Real>::Resample(Vector<Real>& x_trg0, const Vector<Real>& panel_len_trg, const Vector<Long>& ChebOrder_trg, const Vector<Long>& FourierOrder_trg, const Vector<Real>& x_src, const Vector<Real>& panel_len_src, const Vector<Long>& ChebOrder_src, const Vector<Long>& FourierOrder_src) {
    const Long Ns = panel_len_src.Dim();
    const Long Nt = panel_len_trg.Dim();
    if (!Ns || !Nt) return;

    if (Ns > 1 && Nt > 1) {
      Vector<Real> x_trg;
      Long s = 0, t = 0;
      Long src_nds_dsp = 0, trg_nds_dsp = 0;
      while (s < Ns && t < Nt) {
        Real src_len = 0, trg_len = 0;
        Long src_cnt = 0, trg_cnt = 0;
        Long src_nds_cnt = 0, trg_nds_cnt = 0;
        if (panel_len_src[s] > panel_len_trg[t]) {
          src_cnt = 1;
          src_len = panel_len_src[s];
          src_nds_cnt = ChebOrder_src[s] * FourierOrder_src[s];
          while (trg_len < panel_len_src[s] && t+trg_cnt<Nt) {
            trg_nds_cnt += ChebOrder_trg[t+trg_cnt] * FourierOrder_trg[t+trg_cnt];
            trg_len += panel_len_trg[t+trg_cnt];
            trg_cnt++;
          }
        } else {
          trg_cnt = 1;
          trg_len = panel_len_trg[t];
          trg_nds_cnt = ChebOrder_trg[t] * FourierOrder_trg[t];
          while (src_len < panel_len_trg[t] && s+src_cnt<Ns) {
            src_nds_cnt += ChebOrder_src[s+src_cnt] * FourierOrder_src[s+src_cnt];
            src_len += panel_len_src[s+src_cnt];
            src_cnt++;
          }
        }
        SCTL_ASSERT(src_len == trg_len);

        Vector<Real> x_trg_;
        const Vector<Real> x_src_(src_nds_cnt*dof, (Iterator<Real>)x_src.begin() + src_nds_dsp*dof, false);
        Resample<dof>(x_trg_, Vector<Real>(trg_cnt,(Iterator<Real>)panel_len_trg.begin()+t,false), Vector<Long>(trg_cnt,(Iterator<Long>)ChebOrder_trg.begin()+t,false), Vector<Long>(trg_cnt,(Iterator<Long>)FourierOrder_trg.begin()+t,false),
                      x_src_, Vector<Real>(src_cnt,(Iterator<Real>)panel_len_src.begin()+s,false), Vector<Long>(src_cnt,(Iterator<Long>)ChebOrder_src.begin()+s,false), Vector<Long>(src_cnt,(Iterator<Long>)FourierOrder_src.begin()+s,false));
        for (Long i = 0; i < x_trg_.Dim(); i++) x_trg.PushBack(x_trg_[i]);

        src_nds_dsp += src_nds_cnt;
        trg_nds_dsp += trg_nds_cnt;
        s+=src_cnt;
        t+=trg_cnt;
      }
      x_trg0 = x_trg;
      return;
    }

    if (Ns == 1 && Nt == 1 && ChebOrder_src[0]==ChebOrder_trg[0] && FourierOrder_src[0]==FourierOrder_trg[0]) {
      SCTL_ASSERT(panel_len_src[0] == panel_len_trg[0]);
      x_trg0 = x_src;
      return;
    }

    if (Ns == 1) {
      Matrix<Real> Mcheb;
      { // Set Mcheb
        Real trg_offset = 0;
        Vector<Real> trg_nds, src_nds = ChebQuadRule<Real>::nds(ChebOrder_src[0]) * panel_len_src[0];
        for (Long t = 0; t < Nt; t++) {
          const auto& nds = ChebQuadRule<Real>::nds(ChebOrder_trg[t]);
          for (Long i = 0; i < nds.Dim(); i++) trg_nds.PushBack(trg_offset + nds[i] * panel_len_trg[t]);
          trg_offset += panel_len_trg[t];
        }

        Mcheb.ReInit(src_nds.Dim(), trg_nds.Dim());
        Vector<Real> Vcheb(src_nds.Dim() * trg_nds.Dim(), Mcheb.begin(), false);
        LagrangeInterp<Real>::Interpolate(Vcheb, src_nds, trg_nds);
      }
      Matrix<Real> x0 = Mcheb.Transpose()*Matrix<Real>(ChebOrder_src[0], FourierOrder_src[0]*dof, (Iterator<Real>)x_src.begin(), false);

      Vector<Real> x_;
      Long x0_offset = 0;
      for (Long i = 0; i < Nt; i++) {
        { // Rearrange x0
          Matrix<Real> y0(ChebOrder_trg[i] * FourierOrder_src[0], dof, x0.begin() + x0_offset*dof, false);
          Matrix<Real> y1(dof, ChebOrder_trg[i] * FourierOrder_src[0], x0.begin() + x0_offset*dof, false);
          y1 = y0.Transpose();
        }

        Vector<Real> x1(dof * ChebOrder_trg[i] * FourierOrder_trg[i]);
        { // x1 <-- Fourier-resample(x0)
          const auto& Mfourier = FourierResample(FourierOrder_src[0], FourierOrder_trg[i]);
          Matrix<Real> y0(dof * ChebOrder_trg[i], FourierOrder_src[0], x0.begin() + x0_offset*dof, false);
          Matrix<Real> y1(dof * ChebOrder_trg[i], FourierOrder_trg[i], x1.begin(), false);
          y1 = y0 * Mfourier;
        }

        { // Rearrange x1
          Matrix<Real> y0(dof, ChebOrder_trg[i] * FourierOrder_trg[i], x1.begin(), false);
          Matrix<Real> y1(ChebOrder_trg[i] * FourierOrder_trg[i], dof, x1.begin(), false);
          y1 = y0.Transpose();
        }
        for (Long j = 0; j < x1.Dim(); j++) x_.PushBack(x1[j]);

        x0_offset += ChebOrder_trg[i] * FourierOrder_src[0];
      }
      x_trg0 = x_;
      return;
    } else { // Nt == 1
      Vector<Real> src_nds, trg_nds = ChebQuadRule<Real>::nds(ChebOrder_trg[0]) * panel_len_trg[0];

      Vector<Real> x_;
      Real s_offset = 0;
      Long x_src_offset = 0;
      for (Long i = 0; i < Ns; i++) {
        Vector<Real> x0(ChebOrder_src[i] * FourierOrder_src[i] * dof, (Iterator<Real>)x_src.begin() + x_src_offset*dof);
        { // Rearrange x0
          Matrix<Real> y0(ChebOrder_src[i] * FourierOrder_src[i], dof, x0.begin(), false);
          Matrix<Real> y1(dof, ChebOrder_src[i] * FourierOrder_src[i], x0.begin(), false);
          y1 = y0.Transpose();
        }

        Vector<Real> x1(dof * ChebOrder_src[i] * FourierOrder_trg[0]);
        { // x1 <-- Fourier-resample(x0)
          const auto& Mfourier = FourierResample(FourierOrder_src[i], FourierOrder_trg[0]);
          Matrix<Real> y0(dof * ChebOrder_src[i], FourierOrder_src[i], x0.begin(), false);
          Matrix<Real> y1(dof * ChebOrder_src[i], FourierOrder_trg[0], x1.begin(), false);
          y1 = y0 * Mfourier;
        }
        { // Rearrange x1
          Matrix<Real> y0(dof, ChebOrder_src[i] * FourierOrder_trg[0], x1.begin(), false);
          Matrix<Real> y1(ChebOrder_src[i] * FourierOrder_trg[0], dof, x1.begin(), false);
          y1 = y0.Transpose();
        }

        Vector<Real> x2;
        { // Chebyshev-resample(x1)
          const Long a = std::lower_bound(trg_nds.begin(), trg_nds.end(), s_offset) - trg_nds.begin();
          const Long b = std::lower_bound(trg_nds.begin(), trg_nds.end(), s_offset + panel_len_src[i]) - trg_nds.begin();
          Vector<Real> src_nds_ = ChebQuadRule<Real>::nds(ChebOrder_src[i]) * panel_len_src[i] + s_offset;
          Vector<Real> trg_nds_(b-a, trg_nds.begin() + a, false);

          Matrix<Real> Mcheb(src_nds_.Dim(), trg_nds_.Dim());
          Vector<Real> Vcheb(src_nds_.Dim()* trg_nds_.Dim(), Mcheb.begin(), false);
          LagrangeInterp<Real>::Interpolate(Vcheb, src_nds_, trg_nds_);

          x2.ReInit(trg_nds_.Dim()*FourierOrder_trg[0]*dof);
          Matrix<Real> Mx1(src_nds_.Dim(), FourierOrder_trg[0]*dof, x1.begin(), false);
          Matrix<Real> Mx2(trg_nds_.Dim(), FourierOrder_trg[0]*dof, x2.begin(), false);
          Mx2 = Mcheb.Transpose() * Mx1;
        }

        for (Long j = 0; j < x2.Dim(); j++) x_.PushBack(x2[j]);
        x_src_offset += ChebOrder_src[i] * FourierOrder_src[i];
        s_offset += panel_len_src[i];
      }
      x_trg0 = x_;
      return;
    }
  }

  template <class Real> const Matrix<Real>& RigidBodyList<Real>::FourierResample(const Integer N0, const Integer N1) {
    constexpr Integer MaxOrder = 100;
    static Matrix<Matrix<Real>> MM(MaxOrder+1, MaxOrder+1);
    SCTL_ASSERT(N0 <= MaxOrder);
    SCTL_ASSERT(N1 <= MaxOrder);

    Matrix<Real>& M = MM[N0][N1];
    #pragma omp critical(SCTL_FOURIER_RESAMPLE_MAT)
    if (M.Dim(0) * M.Dim(1) == 0) {
      const Integer N0_ = N0/2+1;
      const Integer N1_ = N1/2+1;
      const Integer NN = std::min(N0_, N1_);
      Matrix<Real> M0(NN*2, N0), M1(NN*2, N1);
      for (Long k = 0; k < NN; k++) {
        for (Long i = 0; i < N0; i++) {
          M0[k*2+0][i] = cos<Real>(2*const_pi<Real>()*k*i/N0);
          M0[k*2+1][i] = sin<Real>(2*const_pi<Real>()*k*i/N0);
        }
        for (Long i = 0; i < N1; i++) {
          M1[k*2+0][i] = cos<Real>(2*const_pi<Real>()*k*i/N1);
          M1[k*2+1][i] = sin<Real>(2*const_pi<Real>()*k*i/N1);
        }
      }
      M = M0.pinv(sctl::machine_eps<Real>()*32) * M1;
    }
    return M;
  }

  template <class Real> const Matrix<Real>& RigidBodyList<Real>::FourierErrMatrix(const Integer FourierOrder) {
    constexpr Integer MaxOrder = 100;
    static Vector<Matrix<Real>> M_lst(MaxOrder+1);
    SCTL_ASSERT(FourierOrder <= MaxOrder);

    Matrix<Real>& M = M_lst[FourierOrder];
    #pragma omp critical(SCTL_FOURIER_ERR_MAT)
    if (M.Dim(0) * M.Dim(1) == 0) {
      M = FourierResample(FourierOrder, FourierOrder-2) * FourierResample(FourierOrder-2, FourierOrder);
      for (Long k = 0; k < FourierOrder; k++) M[k][k] -= 1;
    }
    return M;
  }

  template <class Real> const Matrix<Real>& RigidBodyList<Real>::ChebErrMatrix(const Integer ChebOrder) {
    constexpr Integer MaxOrder = 100;
    static Vector<Matrix<Real>> M_lst(MaxOrder+1);
    SCTL_ASSERT(ChebOrder <= MaxOrder);

    Matrix<Real>& M = M_lst[ChebOrder];
    #pragma omp critical(SCTL_CHEB_ERR_MAT)
    if (M.Dim(0) * M.Dim(1) == 0) {
      Vector<Real> Mcheb0, Mcheb1;
      LagrangeInterp<Real>::Interpolate(Mcheb0, ChebQuadRule<Real>::nds(ChebOrder), ChebQuadRule<Real>::nds(ChebOrder-2));
      LagrangeInterp<Real>::Interpolate(Mcheb1, ChebQuadRule<Real>::nds(ChebOrder-2), ChebQuadRule<Real>::nds(ChebOrder));
      M = Matrix<Real>(ChebOrder, ChebOrder-2, Mcheb0.begin(), false) * Matrix<Real>(ChebOrder-2, ChebOrder, Mcheb1.begin(), false);
      for (Long k = 0; k < ChebOrder; k++) M[k][k] -= 1;
    }
    return M;
  }

  template <class Real> void RigidBodyList<Real>::UpdatePointwise(const Vector<Real>& Y_loc) {
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

    const Vector<Real> Y = Allgather(Y_loc, comm_); // parallelize
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
    SCTL_ASSERT(false); // update OrientVec
    InitElemList(loc_elem_cnt, loc_elem_dsp, elem_lst, Xc, ChebOrder, FourierOrder, X, R, OrientVec, obj_elem_cnt, comm_);
  }

  template <class Real> void RigidBodyList<Real>::InitGeom(Vector<Real>& X, Vector<Real>& R, Vector<Real>& OrientVec, Vector<Long>& ChebOrder, Vector<Long>& FourierOrder, Vector<Real>& panel_len, Vector<Real>& Mr_lst, Vector<Long>& cnt, Vector<Long>& dsp, const Long Nobj, const Real loop_rad, const Geom& geom_type) {
    srand48(2);
    const Long ChebOrder0 = 10;
    if (1) { // Set ChebOrder, FourierOrder, panel_len, Mr_lst, cnt, dsp (Uniform discretization)
      const Long Npanel = 16;
      const Long FourierOrder0 = 16;

      //const Long Npanel = 32;
      //const Long FourierOrder0 = 32;

      cnt.ReInit(Nobj); cnt = Npanel;
      dsp.ReInit(Nobj); dsp = 0;
      omp_par::scan(cnt.begin(), dsp.begin(), Nobj);

      ChebOrder.ReInit(Nobj * Npanel);
      FourierOrder.ReInit(Nobj * Npanel);
      panel_len.ReInit(Nobj * Npanel);

      ChebOrder = ChebOrder0;
      FourierOrder = FourierOrder0;
      panel_len = 1/(Real)Npanel;
    } else { // Adaptive discretization
      Vector<Real> panel_len0;
      Vector<Long> FourierOrder0;
      if (1) { // Set panel_len0, FourierOrder0
        panel_len0.PushBack(1.0/(Real)128);
        panel_len0.PushBack(1.0/(Real)128);
        panel_len0.PushBack(1.0/(Real)64 );
        panel_len0.PushBack(1.0/(Real)32 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)32 );
        panel_len0.PushBack(1.0/(Real)64 );
        panel_len0.PushBack(1.0/(Real)128);
        panel_len0.PushBack(1.0/(Real)128);
        panel_len0.PushBack(1.0/(Real)128);
        panel_len0.PushBack(1.0/(Real)128);
        panel_len0.PushBack(1.0/(Real)64 );
        panel_len0.PushBack(1.0/(Real)32 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)16 );
        panel_len0.PushBack(1.0/(Real)32 );
        panel_len0.PushBack(1.0/(Real)64 );
        panel_len0.PushBack(1.0/(Real)128);
        panel_len0.PushBack(1.0/(Real)128);

        FourierOrder0.PushBack(100);
        FourierOrder0.PushBack(100);
        FourierOrder0.PushBack( 72);
        FourierOrder0.PushBack( 40);
        FourierOrder0.PushBack( 24);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 24);
        FourierOrder0.PushBack( 40);
        FourierOrder0.PushBack( 72);
        FourierOrder0.PushBack(100);
        FourierOrder0.PushBack(100);
        FourierOrder0.PushBack(100);
        FourierOrder0.PushBack(100);
        FourierOrder0.PushBack( 72);
        FourierOrder0.PushBack( 40);
        FourierOrder0.PushBack( 24);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 16);
        FourierOrder0.PushBack( 24);
        FourierOrder0.PushBack( 40);
        FourierOrder0.PushBack( 72);
        FourierOrder0.PushBack(100);
        FourierOrder0.PushBack(100);
      } else {
        panel_len0.PushBack(1.0/(Real)64);
        panel_len0.PushBack(1.0/(Real)64);
        panel_len0.PushBack(1.0/(Real)32);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)32);
        panel_len0.PushBack(1.0/(Real)64);
        panel_len0.PushBack(1.0/(Real)64);
        panel_len0.PushBack(1.0/(Real)64);
        panel_len0.PushBack(1.0/(Real)64);
        panel_len0.PushBack(1.0/(Real)32);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)16);
        panel_len0.PushBack(1.0/(Real)32);
        panel_len0.PushBack(1.0/(Real)64);
        panel_len0.PushBack(1.0/(Real)64);

        FourierOrder0.PushBack(72);
        FourierOrder0.PushBack(40);
        FourierOrder0.PushBack(24);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(24);
        FourierOrder0.PushBack(40);
        FourierOrder0.PushBack(72);
        FourierOrder0.PushBack(72);
        FourierOrder0.PushBack(40);
        FourierOrder0.PushBack(24);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(16);
        FourierOrder0.PushBack(24);
        FourierOrder0.PushBack(40);
        FourierOrder0.PushBack(72);
      }

      cnt.ReInit(Nobj); cnt = panel_len0.Dim();
      dsp.ReInit(Nobj); dsp = 0;
      omp_par::scan(cnt.begin(), dsp.begin(), Nobj);

      ChebOrder.ReInit(0);
      FourierOrder.ReInit(0);
      panel_len.ReInit(0);
      for (Long i = 0; i < Nobj; i++) {
        for (Long j = 0; j < cnt[i]; j++) {
          ChebOrder.PushBack(ChebOrder0);
          FourierOrder.PushBack(FourierOrder0[j]);
          panel_len.PushBack(panel_len0[j]);
        }
      }
    }

    auto loop_geom = [&loop_rad](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta){
      x = loop_rad * cos<Real>(theta);
      y = loop_rad * sin<Real>(theta);
      z = 0;
      ex = 0;
      ey = 0;
      ez = 1;
      r = 0.025;
    };
    auto bacteria_geom = [&loop_rad](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta){
      Real t = theta/const_pi<Real>()-1; // -1:1
      Real aspect = const_pi<Real>()*3/2+1;

      Real L = aspect-1+const_pi<Real>()/2;
      Real scal = loop_rad/(1+L-const_pi<Real>()/2) * 0.7;
      if (L*(1+t) < const_pi<Real>()/2) z = scal * (-cos<Real>(L*(1+t)) - L+const_pi<Real>()/2);
      else if (L*(1-t) < const_pi<Real>()/2) z = scal * (cos<Real>(L*(1-t)) + L-const_pi<Real>()/2);
      else z = scal * L * t;

      y = 0;
      x = 0;
      ex = 1;
      ey = 0;
      ez = 0;

      if (L*(1+t) < const_pi<Real>()/2) r = scal * sin<Real>(L*(1+t));
      else if (L*(1-t) < const_pi<Real>()/2) r = scal * sin<Real>(L*(1-t));
      else r = scal;
    };

    X.ReInit(0);
    R.ReInit(0);
    Mr_lst.ReInit(0);
    for (Long i = 0; i < Nobj; i++) {
      Real X0, Y0, Z0;
      { // Set offsets X0, Y0, Z0
        const Long N = (Long)ceil((double)pow<Real>((Real)Nobj,1/(Real)3));
        X0 = (i/pow<0>(N))%N;
        Y0 = (i/pow<1>(N))%N;
        Z0 = (i/pow<2>(N))%N;

        if (geom_type == Geom::Loop) {
          if (Nobj>2) Z0 += drand48()*0.5;
        } else if (geom_type == Geom::Bacteria) {
          if (Nobj>2) {
            X0 = X0*0.8 + drand48()*0.4;
            Y0 = Y0*0.8 + drand48()*0.4;
            Z0 = Z0*1.6 + drand48()*0.4;
          } else {
            X0 = X0*0.2;
          }
        } else {
          SCTL_ASSERT(false); // not implemented
        }
      }

      Real s_dsp = 0;
      for (Long j = 0; j < cnt[i]; j++) { // Set X, OrientVec, R
        const Long ChebOrder_ = ChebOrder[dsp[i]+j];
        const auto& nds = SlenderElemList<Real>::CenterlineNodes(ChebOrder_);
        for (Long k = 0; k < ChebOrder_; k++) {
          Real x, y, z, ex, ey, ez, r;
          Real s = s_dsp + nds[k]*panel_len[dsp[i]+j];
          if (geom_type == Geom::Loop) {
            loop_geom(x, y, z, ex, ey, ez, r, 2*const_pi<Real>()*s);
          } else if (geom_type == Geom::Bacteria) {
            bacteria_geom(x, y, z, ex, ey, ez, r, 2*const_pi<Real>()*s);
          } else {
            SCTL_ASSERT(false); // not implemented
          }
          X.PushBack(x+X0);
          X.PushBack(y+Y0);
          X.PushBack(z+Z0);
          OrientVec.PushBack(ex);
          OrientVec.PushBack(ey);
          OrientVec.PushBack(ez);
          R.PushBack(r);
        }
        s_dsp += panel_len[dsp[i]+j];
      }

      for (Long j = 0; j < COORD_DIM; j++) { // Set Mr
        for (Long k = 0; k < COORD_DIM; k++) {
          Mr_lst.PushBack(j==k?1:0);
        }
      }
    }
  }

  template <class Real> void RigidBodyList<Real>::SlenderElemIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst) {
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

  template <class Real> void RigidBodyList<Real>::ObjIntegral(Vector<Real>& IntegralF, const Vector<Real>& F_cheb, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm comm) {
    Vector<Real> Fe, Fe_local;
    SlenderElemIntegral(Fe_local, F_cheb, elem_lst);
    Fe = Allgather(Fe_local, comm);

    const Long Nobj = obj_elem_cnt.Dim();
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

  template <class Real> template <class ValueType> void RigidBodyList<Real>::InitElemList(Long& loc_elem_cnt, Long& loc_elem_dsp, SlenderElemList<Real>& elem_lst, const Vector<Long>& ChebOrder, const Vector<Long>& FourierOrder, const Vector<ValueType>& X, const Vector<ValueType>& R, const Vector<ValueType>& OrientVec, const Comm comm) {
    // TODO: parallelize

    const Long Nelem = ChebOrder.Dim();
    if (Nelem) { // Set loc_elem_cnt, loc_elem_dsp // TODO: use better cost estimate
      Vector<Long> node_cnt(Nelem), node_dsp(Nelem); node_dsp = 0;
      for (Long i = 0; i < Nelem; i++) node_cnt[i] = ChebOrder[i] * FourierOrder[i] * FourierOrder[i];
      omp_par::scan(node_cnt.begin(), node_dsp.begin(), Nelem);
      const Long Nnodes = node_cnt[Nelem-1] + node_dsp[Nelem-1];

      const Long Np = comm.Size();
      const Long rank = comm.Rank();
      Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+0)/Np) - node_dsp.begin();
      Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+1)/Np) - node_dsp.begin();
      if (rank == Np - 1) b = Nelem;
      if (rank == 0) a = 0;
      loc_elem_cnt = b-a;
      loc_elem_dsp = a;

      if (0 && !comm.Rank()) { // Print partitioning
        std::cout<<"Partitioning: ";
        for (Long i = 0; i < comm.Size(); i++) {
          Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+0)/Np) - node_dsp.begin();
          Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+1)/Np) - node_dsp.begin();
          if (i == Np - 1) b = Nelem;
          if (i == 0) a = 0;
          std::cout<<b-a<<' ';
        }
        std::cout<<'\n';

        std::cout<<"Weight: ";
        for (Long i = 0; i < comm.Size(); i++) {
          Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+0)/Np) - node_dsp.begin();
          Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+1)/Np) - node_dsp.begin();
          if (i == Np - 1) b = Nelem;
          if (i == 0) a = 0;
          std::cout<<node_dsp[b-1]+node_cnt[b-1]-node_dsp[a]<<' ';
        }
        std::cout<<'\n';
      }
    } else {
      loc_elem_cnt = 0;
      loc_elem_dsp = 0;
    }

    { // Set elem_lst
      const Vector<Long> LocChebOrder(loc_elem_cnt, (Iterator<Long>)ChebOrder.begin() + loc_elem_dsp, false);
      const Vector<Long> LocFourierOrder(loc_elem_cnt, (Iterator<Long>)FourierOrder.begin() + loc_elem_dsp, false);

      Long dsp = 0, cnt = 0;
      for (Long i = 0; i < loc_elem_dsp; i++) dsp += ChebOrder[i];
      for (Long i = 0; i < loc_elem_cnt; i++) cnt += ChebOrder[loc_elem_dsp+i];
      const Vector<ValueType> X_(cnt*COORD_DIM, (Iterator<ValueType>)X.begin() + dsp*COORD_DIM, false);
      const Vector<ValueType> R_(cnt, (Iterator<ValueType>)R.begin() + dsp, false);
      const Vector<ValueType> OrientVec_(cnt*COORD_DIM, (Iterator<ValueType>)OrientVec.begin() + dsp*COORD_DIM, false);
      elem_lst.template Init<ValueType>(LocChebOrder, LocFourierOrder, X_, R_, OrientVec_);
    }
  }

  template <class Real> void RigidBodyList<Real>::GetXc(Vector<Real>& Xc, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm& comm) {
    Vector<Real> X_cheb_node_loc;
    elem_lst.GetNodeCoord(&X_cheb_node_loc, nullptr, nullptr);
    ObjIntegral(Xc, X_cheb_node_loc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm);

    Vector<Real> Aobj, A_cheb_node_loc(X_cheb_node_loc.Dim()/COORD_DIM); A_cheb_node_loc = 1;
    ObjIntegral(Aobj, A_cheb_node_loc, elem_lst, obj_elem_cnt, obj_elem_dsp, comm);
    for (Long obj = 0; obj < obj_elem_cnt.Dim(); obj++) {
      for (Long k = 0; k < COORD_DIM; k++) {
        Xc[obj*COORD_DIM+k] /= Aobj[obj];
      }
    }
  }

  template <class Real> void RigidBodyList<Real>::GetNodeObjIdx(Vector<Long>& node_obj_idx, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm comm) {
    // TODO: parallelize
    auto sum = [](const Vector<Long>& V) {
      Long sum = 0;
      for (const auto& x : V) sum+= x;
      return sum;
    };
    Vector<Long> elem_node_dsp, elem_node_cnt, elem_node_cnt_loc;
    elem_lst.GetNodeCoord(nullptr, nullptr, &elem_node_cnt_loc);
    const Long Nloc = sum(elem_node_cnt_loc);

    elem_node_cnt = Allgather(elem_node_cnt_loc, comm); // parallelize
    const Long Nelem = elem_node_cnt.Dim();
    const Long Nobj = obj_elem_cnt.Dim();

    elem_node_dsp.ReInit(Nelem); elem_node_dsp = 0;
    omp_par::scan(elem_node_cnt.begin(), elem_node_dsp.begin(), Nelem);

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

  template <class Real> void RigidBodyList<Real>::RigidBodyVelocity(Vector<Real>& U, const Vector<Real> Uobj, const Vector<Real> Omega_obj, const Vector<Real> Xc, const SlenderElemList<Real>& elem_lst, const Vector<Long>& obj_elem_cnt, const Vector<Long>& obj_elem_dsp, const Comm comm) {
    Vector<Long> node_obj_idx;
    GetNodeObjIdx(node_obj_idx, elem_lst, obj_elem_cnt, obj_elem_dsp, comm);
    const Long Nloc = node_obj_idx.Dim();

    Vector<Real> X_cheb;
    elem_lst.GetNodeCoord(&X_cheb, nullptr, nullptr);

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




  template <class Real> Mobility<Real>::Mobility(const Comm& comm, Real gravity) : comm_(comm), BIOp_StokesFxU(ker_FxU, false, comm), BIOp_StokesFxT(ker_FxT, true, comm), BIOp_StokesDxU(ker_DxU, false, comm), gravity_(gravity), SL_scal((Real)1), DL_scal((Real)1) {
    BIOp_StokesFxT.SetFMMKer(ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxUP, ker_FxT, ker_FxUP, ker_FxT);
    BIOp_StokesDxU.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
    BIOp_StokesFxU.SetFMMKer(ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU);
  }

  template <class Real> void Mobility<Real>::SetSLScaling(const Real SL_scal_) const {
    SL_scal = SL_scal_;
  }

  template <class Real> void Mobility<Real>::BuildPrecond(const RigidBodyList<Real>& geom, const std::string& precond_fname, const Real quad_tol) const {
    geom_precond = geom;
    Mprecond.ReInit(0,0);

    Long precond_size = 0;
    { // Try to read from file
      Vector<Real> x0, x1;
      geom_precond.GetElemList().GetNodeCoord(&x0, nullptr, nullptr);
      x0 = Allgather(x0, comm_);
      precond_size = x0.Dim();

      x1.Read((precond_fname+std::string("-X-")+std::to_string(precond_size)).c_str());
      if (x1.Dim() == x0.Dim()) {
        Real max_err = 0, max_val = 0;
        for (Long i = 0; i < x0.Dim(); i++) {
          max_err = std::max<Real>(max_err, fabs(x0[i]-x1[i]));
          max_val = std::max<Real>(max_val, fabs(x0[i]));
        }
        if (max_err/max_val < quad_tol) {
          Mprecond.Read((precond_fname+std::string("-M-")+std::to_string(precond_size)).c_str());
          if (Mprecond.Dim(0) == x0.Dim() && Mprecond.Dim(1) == x0.Dim()) {
            return;
          } else {
            Mprecond.ReInit(0,0);
          }
        } else {
          SCTL_ASSERT_MSG(false, "Preconditioner file already exists but does not match the current geometry");
        }
      }
      if (Mprecond.Dim(0)*Mprecond.Dim(1) == 0) {
        if (!comm_.Rank()) x0.Write((precond_fname+std::string("-X-")+std::to_string(precond_size)).c_str());
        if (!comm_.Rank()) geom_precond.Write(precond_fname+std::string("-geom-")+std::to_string(precond_size));
      }
    }

    const SlenderElemList<Real>& elem_lst0 = geom_precond.GetElemList();
    BIOp_StokesFxU.AddElemList(elem_lst0);
    BIOp_StokesDxU.AddElemList(elem_lst0);
    BIOp_StokesFxU.SetAccuracy(quad_tol);
    BIOp_StokesDxU.SetAccuracy(quad_tol);
    const Long N = BIOp_StokesFxU.Dim(1);

    auto MobilityOp = [this](Vector<Real>* Ax, const Vector<Real>& x) {
      Vector<Real> Udl, Usl, Kx;
      geom_precond.RigidBodyMotionProj(Kx, x);

      BIOp_StokesFxU.ComputePotential(Usl, (x-Kx));
      BIOp_StokesDxU.ComputePotential(Udl, (x-Kx));
      (*Ax) = ((x-Kx)*0.5 + Udl)*DL_scal + Usl*SL_scal + Kx;
    };

    auto build_matrix = [this](Matrix<Real>& M, std::function<void(Vector<Real>*, const Vector<Real>&)> Op, const Long N) {
      Long Nglb;
      { // Set Nglb
        StaticArray<Long,2>  N_{N,0};
        comm_.Allreduce<Long>(N_+0, N_+1, 1, Comm::CommOp::SUM);
        Nglb = N_[1];
      }
      if (M.Dim(0) != Nglb || M.Dim(1) != Nglb) M.ReInit(Nglb,Nglb);

      Vector<Real> x, Ax;
      for (Long i = 0; i < Nglb; i++) {
        x.ReInit(comm_.Rank() ? 0 : Nglb); x = 0;
        if (!comm_.Rank()) x[i] = 1;
        comm_.PartitionN(x, N);
        Ax.ReInit(N); Ax = 0;
        Op(&Ax, x);

        comm_.PartitionN(Ax, (comm_.Rank()?0:1));
        if (!comm_.Rank()) {
          std::cout<<i<<"/"<<Nglb<<'\n';
          for (Long j = 0; j < Nglb; j++) {
            M[i][j] = Ax[j];
          }
        }
      }

      Vector<Real> Vglb(Nglb*Nglb, M.begin(), false);
      Vglb = Allgather((comm_.Rank()?Vector<Real>():Vglb), comm_);
      SCTL_ASSERT(Vglb.Dim() == Nglb * Nglb);
    };
    auto print_cond_num = [](Matrix<Real> M) {
      static const Real eps_sqrt = sqrt<Real>(machine_eps<Real>());
      Matrix<Real> U, S, Vt;
      M.SVD(U, S, Vt);

      Real S_max=fabs(S[0][0]), S_min=fabs(S[0][0]);
      for (Long i = 0; i < S.Dim(0); i++) {
        if (S[i][i] > eps_sqrt) {
          S_max = std::max<Real>(S_max, fabs(S[i][i]));
          S_min = std::min<Real>(S_min, fabs(S[i][i]));
        }
      }
      std::cout<<"condition-number = "<<S_max/S_min<<'\n';
    };
    SCTL_UNUSED(print_cond_num);

    Matrix<Real> M;
    build_matrix(M, MobilityOp, N);
    //if (!comm_.Rank()) print_cond_num(M);

    Mprecond = M.pinv(sctl::machine_eps<Real>()*32);
    if (!comm_.Rank()) Mprecond.Write((precond_fname+std::string("-M-")+std::to_string(precond_size)).c_str());

    BIOp_StokesFxU.template DeleteElemList<SlenderElemList<Real>>();
    BIOp_StokesDxU.template DeleteElemList<SlenderElemList<Real>>();
  }

  template <class Real> Vector<Real> Mobility<Real>::ComputeVelocity(const RigidBodyList<Real>& geom, const Vector<Real>& Ubg, const Real tol, const Real quad_tol, Vector<Real>* q_, Long* Ngmres_) const {
    const SlenderElemList<Real>& elem_lst0 = geom.GetElemList();
    BIOp_StokesFxU.AddElemList(elem_lst0);
    BIOp_StokesDxU.AddElemList(elem_lst0);
    BIOp_StokesFxU.SetAccuracy(quad_tol);
    BIOp_StokesDxU.SetAccuracy(quad_tol);
    const Long N = BIOp_StokesFxU.Dim(1);

    auto MobilityOp = [this,&geom](Vector<Real>* Ax, const Vector<Real>& x_) {
      Vector<Real> Udl, Usl, Kx, x = x_;
      BIOp_StokesFxU.InvSqrtScaling(x);
      geom.RigidBodyMotionProj(Kx, x);

      BIOp_StokesFxU.ComputePotential(Usl, (x-Kx));
      BIOp_StokesDxU.ComputePotential(Udl, (x-Kx));
      (*Ax) = ApplyPrecond( ((x-Kx)*0.5 + Udl)*DL_scal + Usl*SL_scal + Kx , geom);
      BIOp_StokesFxU.SqrtScaling(*Ax);
    };

    Vector<Real> rhs(N); rhs = 0;
    if (gravity_) {
      Vector<Real> F(N); F = 0;
      for (Long i = 0; i < N/COORD_DIM; i++) {
        F[i*COORD_DIM+2] = gravity_;
      }
      BIOp_StokesFxU.ComputePotential(rhs, F);
    }
    if (Ubg.Dim()) rhs += Ubg;
    rhs = ApplyPrecond(rhs, geom);
    BIOp_StokesFxU.SqrtScaling(rhs);

    Long Ngmres;
    Vector<Real> q;
    ParallelSolver<Real> solver(comm_);
    solver(&q, MobilityOp, rhs, tol, 300, false, &Ngmres);
    BIOp_StokesFxU.InvSqrtScaling(q);

    //Vector<Real> err, Mq;
    //MobilityOp(&Mq, q);
    //err = rhs - Mq;

    BIOp_StokesFxU.template DeleteElemList<SlenderElemList<Real>>();
    BIOp_StokesDxU.template DeleteElemList<SlenderElemList<Real>>();

    Vector<Real> U;
    geom.RigidBodyMotionProj(U, q);
    if (Ngmres_) (*Ngmres_) = Ngmres;
    if (q_) (*q_) = q;
    return U;
  }

  template <class Real> template <class BgFlow> Real Mobility<Real>::TimeStep(RigidBodyList<Real>& geom0, const BgFlow& bg_flow, const Real dt, const SDC<Real>& ode_solver, const Real time_step_tol, const Real quad_tol, const Real gmres_tol) const {
    auto append_vecs = [](const Vector<Real>& A, const Vector<Real>& B) {
      Vector<Real> C(A.Dim() + B.Dim());
      for (Long i = 0; i < A.Dim(); i++) C[i] = A[i];
      for (Long i = 0; i < B.Dim(); i++) C[A.Dim() + i] = B[i];
      return C;
    };
    auto update_geom = [](RigidBodyList<Real>& geom, const Vector<Real>& dY) { // Set geom += dY
      const Long N = dY.Dim() / (COORD_DIM+COORD_DIM*COORD_DIM);
      const Vector<Real> Uc(N*COORD_DIM, (Iterator<Real>)dY.begin(), false);
      Vector<Real> Mr(N*COORD_DIM*COORD_DIM, (Iterator<Real>)dY.begin() + N*COORD_DIM);
      for (Long i = 0; i < N; i++) { // project Mr to orthonormal matrix
        Matrix<Real> M(COORD_DIM,COORD_DIM, Mr.begin() + i*COORD_DIM*COORD_DIM, false);
        Matrix<Real> U, S, Vt, M_ = M;
        M_.SVD(U,S,Vt);
        M = U*Vt;
      }
      geom.RigidBodyUpdate(Uc, Mr);
    };
    auto fn = [this,&append_vecs,&update_geom, &geom0,&bg_flow, &gmres_tol,&quad_tol](Vector<Real>* dYdt, const Vector<Real>& dY) {
      RigidBodyList<Real> geom = geom0;
      update_geom(geom, dY);

      Vector<Real> Uc, Omega;
      { // [Uc, Omega] <-- Mobility(geom)
        Profile::Tic("MobilSolve", &comm_, true);
        bool prof_state = Profile::Enable(false);

        Vector<Real> X, sigma;
        geom.GetElemList().GetNodeCoord(&X, nullptr, nullptr);
        const Vector<Real> Ubg = bg_flow.Velocity(X);
        const Vector<Real> U = ComputeVelocity(geom, Ubg, gmres_tol, quad_tol, &sigma);
        geom.GetRigidBodyMotion(Uc, Omega, U);

        Profile::Enable(prof_state);
        Profile::Toc();
      }

      { // Set dYdt
        auto cross_prod = [](const Tensor<Real,true,3,1>& u, const Tensor<Real,true,3,1>& v) {
          Tensor<Real,true,3,1> uxv;
          uxv(0,0) = u(1,0) * v(2,0) - u(2,0) * v(1,0);
          uxv(1,0) = u(2,0) * v(0,0) - u(0,0) * v(2,0);
          uxv(2,0) = u(0,0) * v(1,0) - u(1,0) * v(0,0);
          return uxv;
        };

        const Long Nobj = Uc.Dim() / COORD_DIM;
        Vector<Real> Mr(Nobj*COORD_DIM*COORD_DIM);
        for (Long obj = 0; obj < Nobj; obj++) {
          Tensor<Real,true,COORD_DIM,1> X, U, W(Omega.begin() + obj*COORD_DIM);
          for (Integer j = 0; j < COORD_DIM; j++) {
            X = 0; X(j,0) = 1;
            U = cross_prod(W, X);
            for (Integer k = 0; k < COORD_DIM; k++) {
              Mr[(obj*COORD_DIM+j)*COORD_DIM+k] = U(k,0);
            }
          }
        }
        dYdt[0] = append_vecs(Uc, Mr);
      }
    };

    Vector<Real> Y0, dY;
    { // Set Y0 // TODO remove
      Vector<Real> X, Uc, Omega, Mr;
      const SlenderElemList<Real> elem_lst0 = geom0.GetElemList();
      elem_lst0.GetNodeCoord(&X, nullptr, nullptr);
      geom0.GetRigidBodyMotion(Uc, Omega, X * 0);

      //geom0.RotationMatrix(Mr, Omega, 0);
      const Long Nobj = Omega.Dim()/COORD_DIM;
      Mr.ReInit(Nobj*COORD_DIM*COORD_DIM);
      for (Long i = 0; i < Nobj; i++) {
        for (Integer k0 = 0; k0 < COORD_DIM; k0++) {
          for (Integer k1 = 0; k1 < COORD_DIM; k1++) {
            Mr[(i*COORD_DIM+k0)*COORD_DIM+k1] = (k0==k1?1:0);
          }
        }
      }

      Y0 = append_vecs(Uc, Mr);
    }

    Real error_interp, error_picard;
    ode_solver(&dY, dt, Y0, fn, -1, time_step_tol, &error_interp, &error_picard);
    update_geom(geom0, dY);

    if (0) { // Print minimum distance in two ring setup
      Vector<Real> x;
      geom0.GetElemList().GetNodeCoord(&x, nullptr, nullptr);
      x = Allgather(x, comm_);

      Vector<Real> min_x(COORD_DIM); min_x = 1e10;
      for (Long i = 0; i < x.Dim()/COORD_DIM; i++) {
        for (Long k = 0; k < COORD_DIM; k++) {
          min_x[k] = std::min<Real>(min_x[k], fabs(x[i*COORD_DIM+k]-0.5));
        }
      }
      if (!comm_.Rank()) {
        std::cout<<"Min-dist = "<<2*min_x[0]<<'\n';
        std::cout<<geom0.panel_len;
        std::cout<<geom0.FourierOrder;
      }
    }

    if (!comm_.Rank()) std::cout<<"error = "<<error_interp/dt<<"  "<<error_picard/dt<<'\n';
    return std::max<Real>(error_interp, error_picard);
  }

  template <class Real> template <class BgFlow> Real Mobility<Real>::AdaptiveTimeStep(RigidBodyList<Real>& geom, Real& dt, const Real T, const BgFlow& bg_flow, const Integer time_step_order, const Real time_step_tol, const Real quad_tol, const Real gmres_tol, const Real geom_tol, const std::string& out_path, Long idx) const {
    const SDC<Real> ode_solver(std::max<Integer>(2,time_step_order), comm_);
    const Real eps = machine_eps<Real>();

    { // Dry run
      bool prof_state = Profile::Enable(false);
      Vector<Real> X, sigma;
      geom.GetElemList().GetNodeCoord(&X, nullptr, nullptr);
      const Vector<Real> Ubg = bg_flow.Velocity(X);
      ComputeVelocity(geom, Ubg, gmres_tol, quad_tol);
      Profile::Enable(prof_state);
    }

    Real t = 0;
    while (t < T && dt > eps*T) {
      Profile::Tic("TimeStep", &comm_, true);
      RigidBodyList<Real> geom_ = geom;

      Vector<Real> U, sigma;
      Long Ngmres = 0, Nunknown = 0;
      { // Compute U, sigma, Ngmres, Nunknown
        Vector<Real> X;
        geom_.GetElemList().GetNodeCoord(&X, nullptr, nullptr);
        const Vector<Real> Ubg = bg_flow.Velocity(X);
        U = ComputeVelocity(geom_, Ubg, gmres_tol, quad_tol, &sigma, &Ngmres);

        Nunknown = [this,&sigma]() {
          StaticArray<Long,2> len{sigma.Dim(),0};
          comm_.Allreduce(len+0, len+1, 1, Comm::CommOp::SUM);
          return len[1];
        }();
      }
      if (WRITE_VTK) [this,&geom_,&U,&sigma,&idx,&out_path]() { // Write output
        RigidBodyList<Real> geom0 = geom_;
        geom0.GetElemList().WriteVTK(out_path + "./S" + std::to_string(idx), sigma, comm_);
        geom0.GetElemList().WriteVTK(out_path + "./U" + std::to_string(idx),     U, comm_);
        geom0.Write(out_path + "./geom" + std::to_string(idx));

        if (1) { // Shift geom0
          Vector<Real> Xc, X0(COORD_DIM); X0 = 0;
          geom0.GetObjPosition(&Xc);
          const Vector<Real> Xc_ = Allgather(Xc, comm_);
          const Long Nobj = Xc.Dim()/COORD_DIM;
          const Long Nobj_ = Xc_.Dim()/COORD_DIM;
          for (Long i = 0; i < Nobj_; i++) {
            for (Long k = 0; k < COORD_DIM; k++) {
              X0[k] += Xc_[i*COORD_DIM+k]/Nobj_;
            }
          }

          Vector<Real> dXc(Nobj*COORD_DIM), Mr(Nobj*COORD_DIM*COORD_DIM);
          for (Long i = 0; i < Nobj; i++) {
            for (Long k0 = 0; k0 < COORD_DIM; k0++) {
              dXc[i*COORD_DIM+k0] = -X0[k0];
              for (Long k1 = 0; k1 < COORD_DIM; k1++) {
                Mr[(i*COORD_DIM+k0)*COORD_DIM+k1] = (k0==k1?1:0);
              }
            }
          }

          geom0.RigidBodyUpdate(dXc*0.95, Mr);
          geom0.GetElemList().WriteVTK(out_path + "./SS" + std::to_string(idx), sigma, comm_);
          geom0.GetElemList().WriteVTK(out_path + "./UU" + std::to_string(idx),     U, comm_);

          geom0 = geom_;
          geom0.RigidBodyUpdate(dXc, Mr);
          geom0.GetElemList().WriteVTK(out_path + "./SSS" + std::to_string(idx), sigma, comm_);
          geom0.GetElemList().WriteVTK(out_path + "./UUU" + std::to_string(idx),     U, comm_);

          //geom.RigidBodyUpdate(dXc, Mr); // TODO: add this shift correctly
        }
      }();

      Real max_err = 0;
      if (time_step_order > 1) {
        max_err = TimeStep(geom_, bg_flow, dt, ode_solver, time_step_tol*dt*0.1, quad_tol, gmres_tol);
      } else { // First order explicit, no adaptivity in time
        Vector<Real> Uc, Omega_c, Mr_lst;
        geom_.GetRigidBodyMotion(Uc, Omega_c, U);
        geom_.RotationMatrix(Mr_lst, Omega_c, dt);
        geom_.RigidBodyUpdate(Uc*dt, Mr_lst);
      }

      bool accept_solution = (max_err <= time_step_tol*dt);
      const Long MaxFourierOrder = [geom_]() {
        Long max_order = 0;
        for (const auto& Nf : geom_.FourierOrder) max_order = std::max<Long>(max_order, Nf);
        return max_order;
      }();
      if (!comm_.Rank()) {
        std::cout<<(accept_solution?"Accepted":"Rejected")<<": idx = "<<idx<<"     t0 = "<<t<<"     dt = "<<dt<<"     err = "<<max_err/dt<<"     Ngmres = "<<Ngmres<<"     Nunknown = "<<Nunknown<<"     MaxFourierOrder = "<<MaxFourierOrder;
        //std::cout<<"     min_dist = "<<geom_.GetMinDist();
        std::cout<<'\n';
      }
      if (accept_solution) { // Accept solution
        geom = geom_;
        t = t + dt;
        idx++;
      }

      // Adjust time-step size (Quaife, Biros - JCP 2016)
      if (time_step_order > 1) dt = std::min<Real>(T-t, 0.7*dt*pow<Real>((time_step_tol*dt)/max_err, 1/(Real)(time_step_order)));
      { // Adaptive refinement
        Vector<Vector<Real>> density_function(1);
        density_function[0].ReInit(sigma.Dim(), sigma.begin(), false);
        geom.RefineAdaptive(density_function, geom_tol);
      }
      Profile::Toc();
      Profile::print(&comm_);
    }
    return t;
  }

  template <class Real> void Mobility<Real>::test(const Comm& comm, const std::string& geom_file, const typename RigidBodyList<Real>::Geom& geom_type, const Long Nobj, const Real loop_rad, const Long start_idx, const Long ts_order, Real dt, const Real T, const Real time_step_tol, const Real gmres_tol, const Real quad_tol, const Real geom_tol, const std::string& precond, const std::string& out_path) {
    RigidBodyList<Real> geom(comm, Nobj, loop_rad, geom_type);
    if (!geom_file.empty()) geom.Read(geom_file);

    RigidBodyList<Real> precond_geom(comm, 1, loop_rad, geom_type);
    if (WRITE_VTK && !precond.empty()) {
      precond_geom.Read(precond, 1); // include loop_rad in filename
      precond_geom.GetElemList().WriteVTK(out_path + "./Sprecond", Vector<Real>(), comm);
    }

    const BgFlow<Real> bg_flow(comm);
    const Mobility<Real> stokes_mobility(comm, 1.0);
    if (geom_type == RigidBodyList<Real>::Geom::Loop) {
      stokes_mobility.SetSLScaling(30.0);
    } else {
      stokes_mobility.SetSLScaling(1.0);
    }
    if (!precond.empty()) stokes_mobility.BuildPrecond(precond_geom, precond);
    stokes_mobility.AdaptiveTimeStep(geom, dt, T, bg_flow, ts_order, time_step_tol, quad_tol, gmres_tol, geom_tol, out_path, start_idx);
  }

  template <class Real> void Mobility<Real>::test_(const Comm& comm, Real gmres_tol, Real quad_tol, Real geom_tol) { // fixed size, explicit first-order time-stepping
    Real dt = 1;
    RigidBodyList<Real> geom(comm, 2);

    const Mobility<Real> stokes_mobility(comm, 1.0);
    for (Long i = 0; i < 20; i++) {
      const SlenderElemList<Real> elem_lst0 = geom.GetElemList();

      Vector<Real> U, sigma;
      { // Compute U, sigma
        Profile::Tic("MobilSolve", &comm, true);
        bool prof_state = Profile::Enable(false);

        U = stokes_mobility.ComputeVelocity(geom, Vector<Real>(), gmres_tol, quad_tol, &sigma);
        if (0) { // print error
          Vector<Real> X;
          elem_lst0.GetNodeCoord(&X, nullptr, nullptr);
          auto max_norm = [](const Vector<Real>& X, const Comm& comm) {
            StaticArray<Real,2> max_val{0,0};
            for (const auto& x : X) max_val[0] = std::max<Real>(max_val[0], fabs(x));
            comm.Allreduce(max_val+0, max_val+1, 1, Comm::CommOp::MAX);
            return max_val[1];
          };
          Vector<Real> Uref = stokes_mobility.ComputeVelocity(geom, Vector<Real>(), 1e-14, 1e-14);
          Vector<Real> Uerr = Uref - U;
          Real max_err = max_norm(Uerr, comm)/max_norm(Uref, comm);
          if (!comm.Rank()) std::cout<<max_err<<'\n';
        }

        Profile::Enable(prof_state);
        Profile::Toc();
      }
      //elem_lst0.Write(std::string("vis/geom") + std::to_string(i), comm);
      //elem_lst0.WriteVTK(std::string("vis/U") + std::to_string(i), U, comm); // Write VTK
      elem_lst0.WriteVTK(std::string("vis/sigma") + std::to_string(i), sigma, comm); // Write VTK

      Vector<Real> Uc, Omega_c, Mr_lst;
      geom.GetRigidBodyMotion(Uc, Omega_c, U);
      geom.RotationMatrix(Mr_lst, Omega_c, dt);
      geom.RigidBodyUpdate(Uc*dt, Mr_lst);

      if (0) { // adaptive refinement
        Vector<Vector<Real>> v(1); v[0] = sigma;
        geom.RefineAdaptive(v, geom_tol);
        geom.RefineAdaptive(Vector<Vector<Real>>(), geom_tol);
      }
    }
  }

  template <class Real> Vector<Real> Mobility<Real>::ApplyPrecond(const Vector<Real>& x, const RigidBodyList<Real>& geom) const {
    if (Mprecond.Dim(0) * Mprecond.Dim(1) == 0) return x;

    Vector<Real> Ax0, x0;
    geom.ApplyPrecond(Ax0, x, Mprecond, geom_precond);
    geom.ApplyPrecond(x0, x, Matrix<Real>(), geom_precond);
    return (x-x0) + Ax0;
  }

}

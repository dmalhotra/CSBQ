#include <csbq.hpp>
using namespace sctl;

const Integer max_depth = 12;

template <class Real> void test(Real geom_tol) {
  #ifdef SCTL_HAVE_PVFMM
  pvfmm::Profile::Enable(true);
  #endif
  Profile::Enable(true);
  const Comm comm = Comm::Self();

  const std::string out_path("vis/");
  const Real min_panel_len = 1/(Real)(1<<max_depth);

  const Long ts_order = 1, start_idx = 0;
  const Real T = 1e-28, ts_tol = 0, gmres_tol = 1e-14, quad_tol = 1e-16;
  Real dt = 1e-28;

  auto build_geom = [&comm, geom_tol](Vector<Real> panel_len_, Vector<Long> FourierOrder) { // Set geom
    constexpr Integer COORD_DIM = 3;
    const Long ElemOrder_ = 10;
    const Long FourierOrder_ = 14;

    Vector<QuadReal> panel_len, X, R, OrientVec;
    { // Set X, R, OrientVec
      panel_len.ReInit(panel_len_.Dim());
      for (Long i = 0; i < panel_len.Dim(); i++) {
        panel_len[i] = (QuadReal)panel_len_[i];
      }

      SlenderElemList<QuadReal> elem_lst;
      if (panel_len.Dim()) {
        if (FourierOrder.Dim() == 0) {
          FourierOrder.ReInit(panel_len.Dim());
          FourierOrder = FourierOrder_;
        }
        GeomTangle(elem_lst, panel_len, FourierOrder, comm, ElemOrder_);
      } else {
        GeomTangle<QuadReal>(elem_lst, Vector<QuadReal>(), Vector<Long>(), comm, ElemOrder_, geom_tol*1e-2);
      }

      Vector<QuadReal> X_surf, Xn;
      Vector<Long> elem_wise_node_cnt;
      elem_lst.GetNodeCoord(&X_surf, &Xn, &elem_wise_node_cnt);
      const Long Npanel = elem_wise_node_cnt.Dim();

      R.ReInit(Npanel*ElemOrder_);
      X.ReInit(Npanel*ElemOrder_*COORD_DIM);
      OrientVec.ReInit(Npanel*ElemOrder_*COORD_DIM);
      Long node_offset = 0;
      for (Long i0 = 0; i0 < Npanel; i0++) {
        const Long FourierOrder_ = elem_wise_node_cnt[i0] / ElemOrder_;
        for (Long i1 = 0; i1 < ElemOrder_; i1++) {
          Long i = i0 * ElemOrder_ + i1;
          QuadReal R2 = 0;
          for (Integer k = 0; k < COORD_DIM; k++) {
            QuadReal sum = 0;
            for (Integer j = 0; j < FourierOrder_; j++) {
              sum += X_surf[(node_offset + i1*FourierOrder_ + j) * COORD_DIM + k];
            }
            X[i*COORD_DIM+k] = sum/FourierOrder_;
            OrientVec[i*COORD_DIM+k] = X_surf[(node_offset + i1*FourierOrder_ + 0) * COORD_DIM + k] - X[i*COORD_DIM+k];
            R2 += OrientVec[i*COORD_DIM+k]*OrientVec[i*COORD_DIM+k];
          }
          R[i] = sctl::sqrt<QuadReal>(R2);
        }
        node_offset += ElemOrder_ * FourierOrder_;
      }
    }
    const Long Npanel = R.Dim()/ElemOrder_;

    Vector<QuadReal> Mr_lst(COORD_DIM*COORD_DIM); Mr_lst = 0;
    for (Integer i = 0; i < COORD_DIM; i++) Mr_lst[i*COORD_DIM+i] = 1;
    Vector<Long> obj_elem_cnt(1); obj_elem_cnt = Npanel;

    Vector<Long> ElemOrder(Npanel); ElemOrder = ElemOrder_;
    if (FourierOrder.Dim() != Npanel) {
      FourierOrder.ReInit(Npanel);
      FourierOrder = FourierOrder_;
    }
    if (panel_len.Dim() != Npanel) {
      panel_len.ReInit(Npanel);
      panel_len = 1;
    }
    RigidBodyList<Real> geom(comm);
    geom.template Init<QuadReal>(obj_elem_cnt, Mr_lst, ElemOrder, FourierOrder, panel_len, X, R, OrientVec);
    return geom;
  };
  RigidBodyList<Real> geom;
  { // Init geom
    Vector<Real> panel_len(16); panel_len = 1/(Real)16;
    geom = build_geom(panel_len, Vector<Long>());
  }

  const auto bg_flow = [](const Vector<Real>& Xt){
    Vector<Real> U(Xt.Dim()); U = 0;
    return U;
  };
  const Mobility<Real> stokes_mobility(comm, 1.0);
  stokes_mobility.SetSLScaling(25);

  Long Npanel = 0;
  for (Long i = 0; i < 20; i++) {
    stokes_mobility.AdaptiveTimeStep(geom, dt, T, bg_flow, ts_order, ts_tol, quad_tol, gmres_tol, geom_tol, out_path, start_idx);
    const auto panel_len = geom.panel_len;
    const auto FourierOrder = geom.FourierOrder;
    { // re-init geom
      Vector<Real> panel_len_;
      Vector<Long> FourierOrder_;
      for (Long i = 0; i < panel_len.Dim(); i++) {
        const Long N = panel_len_.Dim();
        if (N == 0 || panel_len_[N-1] >= min_panel_len) {
          panel_len_.PushBack(panel_len[i]);
          FourierOrder_.PushBack(FourierOrder[i]);
        } else {
          panel_len_[N-1] += panel_len[i];
          FourierOrder_[N-1] = std::max(FourierOrder_[N-1], FourierOrder[i]);
        }
      }
      std::cout<<panel_len_<<'\n';
      std::cout<<FourierOrder_<<'\n';
      std::cout<<panel_len.Dim()<<'\n';
      geom = build_geom(panel_len_, FourierOrder_);
    }
    if (panel_len.Dim() == Npanel) {
      break;
    } else {
      Npanel = panel_len.Dim();
    }
  }

  std::ostringstream geom_tol_str;
  geom_tol_str << geom_tol;
  SlenderElemList<QuadReal> elem_lst;
  { // Set elem_lst
    const Long ElemOrder_ = 10;
    const auto FourierOrder = geom.FourierOrder;
    const auto panel_len_ = geom.panel_len;
    Vector<QuadReal> panel_len(panel_len_.Dim());
    for (Long i = 0; i < panel_len.Dim(); i++) {
      panel_len[i] = (QuadReal)panel_len_[i];
    }
    GeomTangle(elem_lst, panel_len, FourierOrder, comm, ElemOrder_);
  }
  elem_lst.Write(std::string("data/tangle/tangle_")+geom_tol_str.str(), comm);

  Profile::print(&comm);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);
  using Real = long double;

  {
    Comm comm = Comm::Self();
    SCTL_ASSERT_MSG(Comm::World().Size() == 1, "This is a sequential code. Please run with 1 MPI process.");
    commandline_option_start(argc, argv, "\
Generate adaptively deiscretized tangle geometry. The code works by running the\n\
mobility code with the required gemoetry tolereance and dt=0 for a few\n\
iterations unitl adaptive refinement stabilizes. The result is writted to\n\
data/tangle/tangle_${tol}\n", comm);
    double tol = strtod(commandline_option(argc, argv, "-tol", "-1" , false, nullptr, comm), nullptr);
    commandline_option_end(argc, argv);

    if (!comm.Rank()) {
      if (tol > 0) {
        test<Real>(tol);
      } else {
        test<Real>(3e-1);
        test<Real>(1e-1);
        test<Real>(3e-2);
        test<Real>(1e-2);
        test<Real>(3e-3);
        test<Real>(1e-3);
        test<Real>(3e-4);
        test<Real>(1e-4);
        test<Real>(3e-5);
        test<Real>(1e-5);
        test<Real>(3e-6);
        test<Real>(1e-6);
        test<Real>(3e-7);
        test<Real>(1e-7);
        test<Real>(3e-8);
        test<Real>(1e-8);
        test<Real>(3e-9);
        test<Real>(1e-9);
        test<Real>(3e-10);
        test<Real>(1e-10);
        test<Real>(3e-11);
        test<Real>(1e-11);
        test<Real>(3e-12);
        test<Real>(1e-12);
      }
    }
  }

  Comm::MPI_Finalize();
  return 0;
}


#include "utils.hpp"
#include "sctl.hpp"

using namespace sctl;

static constexpr bool cylinder_case = true;
static constexpr Integer FOURIER_ORDER = 20;
static constexpr Integer PERIODIC_REPEAT_COUNT = 3000;

template <class Real> void GeomCylinder(SlenderElemList<Real>& elem_lst0, const Comm& comm, const Real thickness, const int npan) {
  Vector<Long> FourierOrder(npan); FourierOrder = 2*FOURIER_ORDER+4;
  const Long ChebOrder = 10;

  const Long pid = comm.Rank();
  const Long Np = comm.Size();
  const Long k0 = npan*(pid+0)/Np;
  const Long k1 = npan*(pid+1)/Np;

  { // Set elem_lst0
    auto geom = [&thickness](Real& x, Real& y, Real& z, Real& r, const Real theta){
      x = 0;
      y = 0;
      z = theta;
      r = thickness/2;
    };

    Vector<Real> coord, radius;
    Vector<Long> cheb_order, fourier_order;
    for (Long k = k0; k < k1; k++) { // Init elem_lst0
      const auto& nds = SlenderElemList<Real>::CenterlineNodes(ChebOrder);
      for (Long i = 0; i < nds.Dim(); i++) {
        Real x, y, z, r;
        geom(x, y, z, r, 2*const_pi<Real>()*((k-npan/2)+nds[i])*(1/(Real)(2*FOURIER_ORDER)));
        coord.PushBack(x);
        coord.PushBack(y);
        coord.PushBack(z);
        radius.PushBack(r);
      }
      cheb_order.PushBack(ChebOrder);
      fourier_order.PushBack(FourierOrder[k]);
    }
    elem_lst0.Init(cheb_order, fourier_order, coord, radius);
  }
}

template <class Real> void test(const Comm& comm, Real quad_eps, Real radius) {
  constexpr Integer COORD_DIM = 3;
  const Laplace3D_FxU ker_FxU;
  const Laplace3D_DxU ker_DxU;
  BoundaryIntegralOp<Real,Laplace3D_FxU> BIOp_LaplaceFxU(ker_FxU, false, comm);
  BoundaryIntegralOp<Real,Laplace3D_DxU> BIOp_LaplaceDxU(ker_DxU, false, comm);

  Vector<Real> Xtrg;
  SlenderElemList<Real> elem_lst_src, elem_lst_trg;
  if (cylinder_case) { // cylinder
    GeomCylinder<Real>(elem_lst_src, comm, radius*2, 2*FOURIER_ORDER*PERIODIC_REPEAT_COUNT);
    GeomCylinder<Real>(elem_lst_trg, comm, radius*2, 2*FOURIER_ORDER);
  } else { // torus
    Vector<Real> panel_len(4*FOURIER_ORDER); panel_len = 1/(Real)(4*FOURIER_ORDER);
    Vector<Long> FourierOrder(4*FOURIER_ORDER); FourierOrder = 2*FOURIER_ORDER+4;
    GeomEllipse<Real>(elem_lst_src, panel_len, FourierOrder, comm, (Real)1, (Real)1, radius*2);
    GeomEllipse<Real>(elem_lst_trg, panel_len, FourierOrder, comm, (Real)1, (Real)1, radius*2);
  }
  elem_lst_trg.GetNodeCoord(&Xtrg, nullptr, nullptr);
  elem_lst_src.WriteVTK("vis/Xsrc", Vector<Real>(), comm);
  elem_lst_trg.WriteVTK("vis/Xtrg", Vector<Real>(), comm);

  BIOp_LaplaceFxU.AddElemList(elem_lst_src);
  BIOp_LaplaceDxU.AddElemList(elem_lst_src);
  BIOp_LaplaceFxU.SetAccuracy(quad_eps);
  BIOp_LaplaceDxU.SetAccuracy(quad_eps);
  BIOp_LaplaceFxU.SetTargetCoord(Xtrg);
  BIOp_LaplaceDxU.SetTargetCoord(Xtrg);

  if (1) { // build operator matrix in Fourier basis
    auto build_matrix = [&comm, &elem_lst_src, &elem_lst_trg, &BIOp_LaplaceFxU, &BIOp_LaplaceDxU](const std::string& fname) {
      Vector<Real> Xsrc, Xtrg, F;
      elem_lst_src.GetNodeCoord(&Xsrc, nullptr, nullptr);
      elem_lst_trg.GetNodeCoord(&Xtrg, nullptr, nullptr);
      const Long Nsrc = Xsrc.Dim() / COORD_DIM;
      const Long Ntrg = Xtrg.Dim() / COORD_DIM;

      Matrix<Real>  Mfourier_basis, Mfourier_transform;
      { // Set Mfourier
        Vector<Real> wts;
        { // Set wts
          Vector<Real> X, Xn, dist_far;
          Vector<Long> element_wise_node_cnt;
          elem_lst_trg.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        }
        const Long NN = wts.Dim();
        Matrix<Real> Mwts(NN, NN); Mwts = 0;
        for (Long i = 0; i < NN; i++) Mwts[i][i] = wts[i];

        Mfourier_basis.ReInit(pow<2>(FOURIER_ORDER), NN);
        for (Long k0 = 0; k0 < FOURIER_ORDER; k0++) {
          for (Long k1 = 0; k1 < FOURIER_ORDER; k1++) {
            Vector<Real> F(Ntrg);
            for (Long i = 0; i < Ntrg; i++) {
              if (cylinder_case) {
                Real s = Xtrg[i*COORD_DIM+2];
                Real cos_theta = Xtrg[i*COORD_DIM+0];
                Real sin_theta = Xtrg[i*COORD_DIM+1];
                Real theta = atan2(sin_theta, cos_theta);
                F[i] = cos<Real>(theta*k0) * cos<Real>(s*k1);
              } else {
                Real s = atan2(Xtrg[i*COORD_DIM+1],Xtrg[i*COORD_DIM+0]);
                Real R = sqrt<Real>(Xtrg[i*COORD_DIM+0]*Xtrg[i*COORD_DIM+0] + Xtrg[i*COORD_DIM+1]*Xtrg[i*COORD_DIM+1]);
                Real theta = atan2(Xtrg[i*COORD_DIM+2], R-1);
                F[i] = cos<Real>(s*k0) * cos<Real>(theta*k1);
              }
            }
            Vector<Real> F_(NN, Mfourier_basis[k0*FOURIER_ORDER+k1], false);
            elem_lst_trg.GetFarFieldDensity(F_, F);
          }
        }
        Mfourier_transform = Mwts*Mfourier_basis.Transpose();

        Matrix<Real> Mtmp(pow<2>(FOURIER_ORDER), pow<2>(FOURIER_ORDER)); Mtmp = 0;
        if (NN) Mtmp = Mfourier_basis * Mfourier_transform;
        Matrix<Real> Mtmp_glb(pow<2>(FOURIER_ORDER), pow<2>(FOURIER_ORDER));
        comm.Allreduce(Mtmp.begin(), Mtmp_glb.begin(), pow<4>(FOURIER_ORDER), Comm::CommOp::SUM);
        for (Long i = 0; i < pow<2>(FOURIER_ORDER); i++) {
          Real scal = 1/Mtmp_glb[i][i];
          for (Long j = 0; j < NN; j++) {
            Mfourier_transform[j][i] *= scal;
          }
        }
      }
      const Long NN = Mfourier_transform.Dim(0);

      F.ReInit(Nsrc);
      Matrix<Real> Ms(pow<2>(FOURIER_ORDER), NN); Ms = 0;
      Matrix<Real> Md(pow<2>(FOURIER_ORDER), NN); Md = 0;
      for (Long k0 = 0; k0 < FOURIER_ORDER; k0++) {
        for (Long k1 = 0; k1 < FOURIER_ORDER; k1++) {
          for (Long i = 0; i < Nsrc; i++) {
            if (cylinder_case) {
              Real s = Xsrc[i*COORD_DIM+2];
              Real cos_theta = Xsrc[i*COORD_DIM+0];
              Real sin_theta = Xsrc[i*COORD_DIM+1];
              Real theta = atan2(sin_theta, cos_theta);
              F[i] = cos<Real>(theta*k0) * cos<Real>(s*k1);
            } else {
              Real s = atan2(Xsrc[i*COORD_DIM+1],Xsrc[i*COORD_DIM+0]);
              Real R = sqrt<Real>(Xsrc[i*COORD_DIM+0]*Xsrc[i*COORD_DIM+0] + Xsrc[i*COORD_DIM+1]*Xsrc[i*COORD_DIM+1]);
              Real theta = atan2(Xsrc[i*COORD_DIM+2], R-1);
              F[i] = cos<Real>(s*k0) * cos<Real>(theta*k1);
            }
          }
          { // single-layer eval
            Vector<Real> V, V_(Ntrg, Ms[k0*FOURIER_ORDER+k1], false);
            BIOp_LaplaceFxU.ComputePotential(V, F);
            elem_lst_trg.GetFarFieldDensity(V_, V);
          }
          { // double-layer eval
            Vector<Real> V, V_(Ntrg, Md[k0*FOURIER_ORDER+k1], false);
            BIOp_LaplaceDxU.ComputePotential(V, F);
            elem_lst_trg.GetFarFieldDensity(V_, V);
          }
        }
      }
      Md += Mfourier_basis*0.5;

      Matrix<Real> Ms_(pow<2>(FOURIER_ORDER), pow<2>(FOURIER_ORDER)); Ms_ = 0;
      Matrix<Real> Md_(pow<2>(FOURIER_ORDER), pow<2>(FOURIER_ORDER)); Md_ = 0;
      { // Set Ms_
        Matrix<Real> Mtmp(pow<2>(FOURIER_ORDER), pow<2>(FOURIER_ORDER)); Mtmp = 0;
        if (NN) Mtmp = Ms*Mfourier_transform;
        comm.Allreduce(Mtmp.begin(), Ms_.begin(), pow<4>(FOURIER_ORDER), Comm::CommOp::SUM);
      }
      { // Set Md_
        Matrix<Real> Mtmp(pow<2>(FOURIER_ORDER), pow<2>(FOURIER_ORDER)); Mtmp = 0;

        Mtmp.SetZero();
        if (NN) Mtmp = Md*Mfourier_transform;
        comm.Allreduce(Mtmp.begin(), Md_.begin(), pow<4>(FOURIER_ORDER), Comm::CommOp::SUM);
      }
      if (!comm.Rank()) Ms_.Write((fname+"s").c_str());
      if (!comm.Rank()) Md_.Write((fname+"d").c_str());

      if (!comm.Rank()) {
        std::cout<<Ms_<<'\n';
        std::cout<<Md_<<'\n';
      }
    };
    build_matrix("OpMat/M");
  }

  if (1) { // build operator matrix in nodal basis [deprecated, not maintained]
    //Real SL_scal = 1.0;
    //Real DL_scal = 1.0;
    //auto MobilityOp = [&BIOp_LaplaceFxU,&BIOp_LaplaceDxU,&SL_scal,&DL_scal](Vector<Real>* Ax, const Vector<Real>& x) {
    //  Vector<Real> Ud, Us;
    //  BIOp_LaplaceFxU.ComputePotential(Us, x);
    //  BIOp_LaplaceDxU.ComputePotential(Ud, x);
    //  (*Ax) = (x*0.5 + Ud)*DL_scal + Us*SL_scal;
    //};

    //auto build_matrix = [&comm, &elem_lst_src, &elem_lst_trg, &BIOp_LaplaceFxU, &BIOp_LaplaceDxU](const std::string& fname) {
    //  Matrix<double> M;
    //  Vector<Real> x, Ax;
    //  if (!comm.Rank()) M.ReInit(Nglb,Nglb);
    //  for (Long i = 0; i < Nglb; i++) {
    //    x.ReInit(comm.Rank() ? 0 : Nglb); x = 0;
    //    if (!comm.Rank()) x[i] = 1;
    //    comm.PartitionN(x, N);
    //    Ax.ReInit(N); Ax = 0;

    //    comm.PartitionN(Ax, (comm.Rank()?0:1));
    //    if (!comm.Rank()) {
    //      std::cout<<i<<"/"<<Nglb<<'\n';
    //      for (Long j = 0; j < Nglb; j++) {
    //        M[i][j] = (double)Ax[j];
    //      }
    //    }
    //  }
    //  if (!comm.Rank()) M.Write(fname.c_str());
    //  return M;
    //};

    //SL_scal = 1; DL_scal = 0;
    //auto Ms = build_matrix("OpMat/Ms", MobilityOp, BIOp_LaplaceFxU.Dim(0));
    //SL_scal = 0; DL_scal = 1;
    //auto Md = build_matrix("OpMat/Md", MobilityOp, BIOp_LaplaceDxU.Dim(0));

    //auto print_cond_num = [](Matrix<Real> M) {
    //  Matrix<Real> U, S, Vt;
    //  M.SVD(U, S, Vt);
    //  Real S_max=S[0][0], S_min=S[0][0];
    //  for (Long i = 0; i < S.Dim(0); i++) {
    //    if (S[i][i] > 1e-8) {
    //      S_max = std::max<Real>(S_max, fabs(S[i][i]));
    //      S_min = std::min<Real>(S_min, fabs(S[i][i]));
    //    }
    //  }
    //  std::cout<<S_max/S_min<<'\n';
    //};
    //if (!comm.Rank()) {
    //  //omp_set_num_threads(16);
    //  print_cond_num(Md);
    //  print_cond_num(Ms*.125 + Md);
    //  print_cond_num(Ms*.250 + Md);
    //  print_cond_num(Ms*.500 + Md);
    //  print_cond_num(Ms*1.00 + Md);
    //  print_cond_num(Ms*2.00 + Md);
    //  print_cond_num(Ms*4.00 + Md);
    //  print_cond_num(Ms*8.00 + Md);
    //  print_cond_num(Ms*16.0 + Md);
    //  print_cond_num(Ms*32.0 + Md);
    //  print_cond_num(Ms*64.0 + Md);
    //  print_cond_num(Ms);
    //}
  }

  BIOp_LaplaceFxU.template DeleteElemList<SlenderElemList<Real>>();
  BIOp_LaplaceDxU.template DeleteElemList<SlenderElemList<Real>>();
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    #ifdef SCTL_HAVE_PVFMM
    pvfmm::Profile::Enable(true);
    #endif
    Profile::Enable(true);
    Comm comm = Comm::World();
    commandline_option_start(argc, argv, nullptr, comm);
    double quad_eps = strtod(commandline_option(argc, argv, "-quad_eps", "1e-14", false, nullptr, comm), nullptr);
    double radius = strtod(commandline_option(argc, argv, "-eps", "0.1", false, nullptr, comm), nullptr);
    commandline_option_end(argc, argv);

    test<double>(comm, quad_eps, radius);
    //Profile::print(&comm);
  }

  Comm::MPI_Finalize();
  return 0;
}


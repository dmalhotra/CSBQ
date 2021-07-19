#include "utils.hpp"
#include "sctl.hpp"

using namespace sctl;

template <class Real, class KerSL, class KerDL> void bvp_solve_(const std::string& fname, const Real thickness, const Real tol, const Real gmres_tol, const Real SLScaling, const Comm& comm = Comm::World()) {
  SlenderElemList<Real> elem_lst0;
  if (fname.empty()) { // Initialize elem_lst0 in code
    GeomTorus(elem_lst0, Vector<Real>(), Vector<Long>(), comm, thickness, 10);
  } else {
    elem_lst0.Read(fname, comm);
  }

  CubeVolumeVis<Real> cube_vol(50, 1.5);
  Vector<Real> U = bvp_solve<Real, KerSL, KerDL>(elem_lst0, tol, gmres_tol, SLScaling, Vector<Real>(), cube_vol.GetCoord(), comm);
  cube_vol.WriteVTK("vis/BVP-cube", U, comm);
}

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    const Comm comm = Comm::World();
    commandline_option_start(argc, argv, nullptr, comm);
    std::string fname = commandline_option(argc, argv, "-geom", ""      , false, "Input geometry filename", comm);
    std::string ker   = commandline_option(argc, argv, "-ker" , "Stokes", false, "Stokes/Laplace"         , comm);

    double thickness = strtod(commandline_option(argc, argv, "-r"        , "1e-03", false, "Thickness of default geometry" , comm), nullptr);
    double tol       = strtod(commandline_option(argc, argv, "-tol"      , "1e-10", false, ""                              , comm), nullptr);
    double gmres_tol = strtod(commandline_option(argc, argv, "-gmres_tol", "1e-08", false, ""                              , comm), nullptr);
    double scaling   = strtod(commandline_option(argc, argv, "-scale", "1e3", false, "Single-layer operator scaling factor", comm), nullptr);
    commandline_option_end(argc, argv);

    if ((!strcmp(ker.c_str(), "Stokes"))) {
      if (!comm.Rank()) std::cout<<"\n\n## Solving : (1/2 + D + S * scaling) sigma = V0    (Stokes exterior Dirichlet BVP)\n";
      bvp_solve_<double, Stokes3D_FxU, Stokes3D_DxU>(fname, thickness, tol, gmres_tol, scaling, comm);
    } else if ((!strcmp(ker.c_str(), "Laplace"))) {
      if (!comm.Rank()) std::cout<<"\n\n## Solving : (1/2 + D + S * scaling) sigma = V0    (Laplace exterior Dirichlet BVP)\n";
      bvp_solve_<double, Laplace3D_FxU, Laplace3D_DxU>(fname, thickness, tol, gmres_tol, scaling, comm);
    } else {
      if (!comm.Rank()) std::cout<<"Unknown kernel "<<ker<<'\n';
    }
  }

  Comm::MPI_Finalize();
  return 0;
}


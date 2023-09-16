#include "utils.hpp"
#include "mobility.hpp"
#include "sctl.hpp"

using namespace sctl;

int main(int argc, char** argv) {
  Comm::MPI_Init(&argc, &argv);

  {
    #ifdef SCTL_HAVE_PVFMM
    pvfmm::Profile::Enable(true);
    #endif
    Profile::Enable(true);
    Comm comm = Comm::World();
    commandline_option_start(argc, argv, nullptr, comm);
    Long Nobj       = strtol(commandline_option(argc, argv, "-Nobj"     , "2"    , false, nullptr, comm), nullptr, 10);
    double ts_tol   = strtod(commandline_option(argc, argv, "-ts_tol"   , "1e-7" , false, nullptr, comm), nullptr);
    double gmres_tol= strtod(commandline_option(argc, argv, "-gmres_tol", "1e-10", false, nullptr, comm), nullptr);
    double quad_tol = strtod(commandline_option(argc, argv, "-quad_tol" , "1e-10", false, nullptr, comm), nullptr);
    double geom_tol = strtod(commandline_option(argc, argv, "-geom_tol" , "1e-8" , false, nullptr, comm), nullptr);
    double loop_rad = strtod(commandline_option(argc, argv, "-loop_rad" , "0.45" , false, nullptr, comm), nullptr);
    std::string out_path =   commandline_option(argc, argv, "-out_path" , "vis/" , false, nullptr, comm);
    std::string precond  =   commandline_option(argc, argv, "-precond"  , ""     , false, nullptr, comm);

    std::string geom =      commandline_option(argc, argv, "-geom"     , ""    , false, nullptr, comm);
    Long ts_order  = strtol(commandline_option(argc, argv, "-ts_order" , "5"   , false, nullptr, comm), nullptr, 10);
    Long start_idx = strtol(commandline_option(argc, argv, "-start_idx", "0"   , false, nullptr, comm), nullptr, 10);
    double dt      = strtod(commandline_option(argc, argv, "-dt"       , "0.1" , false, nullptr, comm), nullptr);
    double T       = strtod(commandline_option(argc, argv, "-T"        , "1000", false, nullptr, comm), nullptr);
    Long omp_p     = strtol(commandline_option(argc, argv, "-omp"      , "1"   , false, nullptr, comm), nullptr, 10);
    commandline_option_end(argc, argv);
    omp_set_num_threads(omp_p);

    if (!comm.Rank()) {
      for (Integer i = 0; i < argc; i++) {
        std::cout<<argv[i]<<'\n';
      }
    }

    comm.Barrier();
    //Mobility<double>::test(comm, geom, RigidBodyList<double>::Geom::Bacteria, Nobj, loop_rad, start_idx, ts_order, dt, T, ts_tol, gmres_tol, quad_tol, geom_tol, precond, out_path);
    Mobility<double>::test(comm, geom, RigidBodyList<double>::Geom::Loop, Nobj, loop_rad, start_idx, ts_order, dt, T, ts_tol, gmres_tol, quad_tol, geom_tol, precond, out_path);

    //Mobility<double>::test_(comm);

    Profile::print(&comm);
  }

  Comm::MPI_Finalize();
  return 0;
}


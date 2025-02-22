// Author: Clara Bayley
// File: main.cpp
/* This file runs the entire superdrop model (SDM)
coupled with a CVODE ode solver for the thermodynamics
(p, temp, qv and qc) over time */

// after make/compiling, execute for example via:
// ./src/adiabatic0D "../src/config/condconfig.txt" "../libs/claras_SDconstants.hpp"

/* standard library packages */
#include <vector>
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>

/* constants and equations */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* sdm gridboxes setup */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/runsdmstep.hpp"
#include "sdmgridboxes/sdmotion.hpp"
#include "sdmgridboxes/detectors.hpp"
#include "sdmgridboxes/detectors_ptr.hpp"
#include "sdmgridboxes/logbooks.hpp"

/* sdm observers setup */
#include "observers/observers.hpp"
#include "observers/observegbxs.hpp"
#include "observers/observelbks.hpp"
#include "observers/gridboxes_intostore.hpp"

#include "singlevarstorage.hpp"
#include "zarrstorage/sdattributes_intostore.hpp"
#include "zarrstorage/contigraggedsdstorage.hpp"
#include "zarrstorage/thermostatestorage.hpp"
#include "zarrstorage/zarrstores.hpp"

/* sdm superdroplets setup */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/condensation.hpp"

/* thermodynamics solver and coupled model setup */
#include "cvodecoupld/run_cvodecoupld.hpp"

namespace dlc = dimless_constants;

template <SuperdropIntoStoreViaBuffer S>
struct SomeZarrStores
{
  ThermoStateStorage thermozarr;
  ContiguousRaggedSDStorage<S> sdzarr;
  CoordinateStorage<double> timezarr;

  SomeZarrStores(FSStore &store, const int maxchunk,
                 const unsigned int ngbxs, S sdattrs)
      : thermozarr(store, maxchunk, ngbxs),
        sdzarr(store, sdattrs, maxchunk),
        timezarr(make_timezarr(store, maxchunk)) {}
};

SuperdropIntoStoreViaBuffer auto sdattrs_to_observe()
/* choose which methods are used to write attributes
of a superdroplet into zarr storage */
{
  SuperdropIntoStoreViaBuffer auto id = IdIntoStore();
  SuperdropIntoStoreViaBuffer auto eps = EpsIntoStore();
  SuperdropIntoStoreViaBuffer auto radius = RadiusIntoStore();
  SuperdropIntoStoreViaBuffer auto m_sol = M_solIntoStore();

  SuperdropIntoStoreViaBuffer auto attrs = id >> eps >> radius >> m_sol;

  return attrs;
}

template <SuperdropIntoStoreViaBuffer S>
Observer auto create_observer(SomeZarrStores<S> &stores, const int obsstep)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const ObserveGBxs auto obs3 = ObserveThermoState(stores.thermozarr);
  const ObserveGBxs auto obs2 = ObserveSDsAttributes(stores.sdzarr);
  const ObserveGBxs auto obs1 = ObserveTime(stores.timezarr);
  const auto obsgbxs = obs3 >> obs2 >> obs1;

  const Observer auto observer = PrintObserver(obsstep) >>
                                 ConstIntervalGBxsObserver(obsstep, obsgbxs);

  return observer;
}

int main(int argc, char *argv[])
{
  Kokkos::Timer kokkostimer;

  if (argc < 3)
  {
    throw std::invalid_argument("config and/or constants files not specified");
  }

  /* object containing input parameters from configuration file */
  const std::string configfilepath = argv[1];    // path to configuration (.txt file)
  const std::string constantsfilepath = argv[2]; // path to constants (.hpp file)
  const Config config(configfilepath, constantsfilepath);

  /* object for time-stepping parameters of coupled model */
  const SDMTimesteps mdlsteps(config.CONDTSTEP, config.COLLTSTEP,
                              config.MOTIONTSTEP, config.COUPLTSTEP,
                              config.OBSTSTEP, config.T_END);

  /* create map from gridbox index to its coordinate boundaries */
  const Maps4GridBoxes gbxmaps(config.SDnspace, config.grid_filename);

  /* create superdroplet model (SDM) process from combination of chosen SDM processes */
  const double cond_subtstep = realtime2dimless(config.cond_SUBTSTEP);
  const auto sdmprocess(CondensationProcess(mdlsteps.condsubstep, &step2dimlesstime,
                                            config.doAlterThermo, config.cond_iters,
                                            cond_subtstep, config.cond_rtol,
                                            config.cond_atol));
  const MoveSuperdropsInDomain sdmmotion(NullMotion{});

  /* create observer from combination of chosen observers */
  FSStore fsstore(config.zarrbasedir);
  SomeZarrStores zarrstores(fsstore,
                            config.maxchunk,
                            gbxmaps.ngridboxes,
                            sdattrs_to_observe());
  const auto observer = create_observer(zarrstores, mdlsteps.obsstep);

  const NullDetectorsPtr dtrs{};

  const RunSDMStep sdm(gbxmaps, sdmmotion, sdmprocess, observer);

  Kokkos::initialize(argc, argv);
  {
    /* RUN SDM MODEL COUPLED TO CVODE ODE SOLVER */
    run_cvodecoupld(config, sdm, dtrs, mdlsteps.t_end, mdlsteps.couplstep);
  }
  Kokkos::finalize();
  std::cout << "  ------ Total Duration: " << kokkostimer.seconds() << "s ----- \n";

  return 0;
}
// Author: Clara Bayley
// File: main.cpp
/* This file runs the entire superdrop model (SDM)
coupled with a CVODE ode solver for the thermodynamics
(p, temp, qv and qc) over time */

// after make/compiling, execute for example via:
// ./src/colls0D "../src/config/config.txt" "../libs/claras_SDconstants.hpp"

/* standard library packages */
#include <vector>
#include <stdexcept>
#include <string>

/* constants and equations */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* Coupled model setup */
#include "cvodesdm/run_cvodesdm.hpp"

/* sdm gridboxes setup */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/run_sdmstep.hpp"
#include "observers/observers.hpp"
#include "observers/intostore_observers.hpp"
#include "observers/sdattributes_intostore.hpp"
#include "observers/contigraggedsdstorage.hpp"
#include "observers/thermostatestorage.hpp"
#include "observers/zarrstores.hpp"

/* Superdroplet Model (SDM) */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/sdmotion.hpp"
#include "superdrop_solver/coalescencekernel.hpp"
#include "superdrop_solver/collisionsmethod.hpp"

namespace dlc = dimless_constants;

template <SuperdropIntoStoreViaBuffer S>
struct SomeZarrStores
{
  ThermoStateStorage thermozarr;
  ContiguousRaggedSDStorage<S> sdzarr;
  CoordinateStorage<double> timezarr;

SomeZarrStores(FSStore &fsstore, const int maxcsize,
              const unsigned int ngridboxes, S sdattrs)
      : thermozarr(fsstore, maxcsize, ngridboxes),
        sdzarr(fsstore, sdattrs, maxcsize),
        timezarr(fsstore, maxcsize, "time",
                 "<f8", "s", dlc::TIME0) {}
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
Observer auto create_observer(SomeZarrStores<S> &stores)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const Observer auto obs3 = ThermoStateObserver(stores.thermozarr);
  const Observer auto obs2 = SDsAttributeObserver(stores.sdzarr);
  const Observer auto obs1 = TimeObserver(stores.timezarr);
  
  const auto observer = obs3 >> obs2 >> obs1 >> PrintObserver{};

  return observer;
}

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    throw std::invalid_argument("config and/or constants files not specified");
    return -1;
  }

  /* object containing input parameters from configuration file */
  const std::string configfilepath = argv[1];    // path to configuration (.txt file)
  const std::string constantsfilepath = argv[2]; // path to constants (.hpp file)
  const Config config(configfilepath, constantsfilepath);

  /* object for time-stepping parameters of coupled model */
  const SDMTimesteps mdlsteps(config.CONDTSTEP, config.COLLTSTEP,
                                config.MOTIONTSTEP, config.COUPLTSTEP,
                                config.T_END);

  /* create map from gridbox index to its coordinate boundaries */
  const Maps4GridBoxes gbxmaps(config.SDnspace, config.grid_filename);

  /* create superdroplet model (SDM) process from combination of chosen SDM processes */
  const auto sdmprocess(CollisionsProcess(mdlsteps.collsubstep,
                                       &step2realtime,
                                       GolovinProb(dlc::R0)));
    
  const auto sdmotion(NullMotion{});

  /* create observer from combination of chosen observers */
  FSStore fsstore(config.zarrbasedir);
  SomeZarrStores zarrstores(fsstore, config.maxcsize,
                            gbxmaps.gbxidxs.size(),
                            sdattrs_to_observe());
  const auto observer = create_observer(zarrstores);

  const RunSDMStep sdm(gbxmaps, sdmotion, sdmprocess, observer);

  /* RUN SDM MODEL COUPLED TO CVODE ODE SOLVER */
  run_cvodesdm(config, sdm, mdlsteps.t_end, mdlsteps.couplstep);

  return 0;
}
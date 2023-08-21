// Author: Clara Bayley
// File: main_breakup.cpp
/* This file runs the entire superdrop model (SDM)
thermodynamics read from binaries files
(p, temp, qv and qc) over time for SU18
collision-coalescence-breakup-rebound kernel */

// after make/compiling, execute for example via:
// ./src/runbreakup "../src/config/config.txt" "../libs/claras_SDconstants.hpp"

/* standard library packages */
#include <stdexcept>
#include <string>

#include <Kokkos_Core.hpp>

/* constants and configuration */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* sdm gridboxes setup */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/runsdmstep.hpp"
#include "sdmgridboxes/sdmotion.hpp"
#include "sdmgridboxes/detectors.hpp"
#include "sdmgridboxes/detectors_ptr.hpp"

/* sdm observers setup */
#include "observers/observers.hpp"
#include "observers/observegbxs.hpp"
#include "observers/gridboxes_intostore.hpp"

#include "singlevarstorage.hpp"
#include "massmomentsstorage.hpp"
#include "zarrstorage/sdattributes_intostore.hpp"
#include "zarrstorage/contigraggedsdstorage.hpp"
#include "zarrstorage/zarrstores.hpp"

/* sdm superdroplets setup */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/coal_breakup_rebound.hpp"
#include "superdrop_solver/terminalvelocity.hpp"

/* thermodynamics solver and coupled model setup */
#include "thermofromfile/run_thermofromfile.hpp"
#include "thermofromfile/prescribedmotion.hpp"

namespace dlc = dimless_constants;

template <SuperdropIntoStoreViaBuffer S>
struct SomeZarrStores
{
  ContiguousRaggedSDStorage<S> sdzarr;
  MomentsStorages massmomszarr;
  CoordinateStorage<double> timezarr;
  TwoDStorage<size_t> nsuperszarr;

  SomeZarrStores(FSStore &store, const int maxchunk,
                 const unsigned int ngbxs, S sdattrs)
      : sdzarr(store, sdattrs, maxchunk),
        massmomszarr(store, maxchunk, ngbxs),
        timezarr(make_timezarr(store, maxchunk)),
        nsuperszarr(make_nsuperszarr(store, maxchunk, ngbxs)) {}
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

SdmProcess auto create_sdmprocess(const Config &config,
                                  const SDMTimesteps &mdlsteps)
/* return an SdmProcess type from SU18 
collision-coalescence-breakup-rebound method */
{
  // const auto probs_coal(CollCoalProb_Long());
  // const auto sdmprocess(CollisionCoalescenceProcess(mdlsteps.collsubstep,
  //                                             &step2realtime,
  //                                             probs_coal));

  const auto probs_coll(CollCoalProb_Long());
  const auto terminalv(SimmelTerminalVelocity{});
  const auto sdmprocess(CollisionAllProcess(mdlsteps.collsubstep,
                                         &step2realtime,
                                         probs_coll,
                                         terminalv,
                                         config.nfrags));
  return sdmprocess;
}

ObserveGBxs auto
create_observegbx_massmoments(MomentsStorages &mms,
                              const size_t ngbxs)
{
  const auto mom0 = ObserveNthMassMoment(mms.mom0zarr, 0, ngbxs);
  const auto mom1 = ObserveNthMassMoment(mms.mom1zarr, 1, ngbxs);
  const auto mom2 = ObserveNthMassMoment(mms.mom2zarr, 2, ngbxs);

  return mom2 >> mom1 >> mom0; 
}

template <SuperdropIntoStoreViaBuffer S>
Observer auto create_observer(SomeZarrStores<S> &stores, const int obsstep,
                              const size_t ngbxs)
/* return an Observer type from an amalgamation of other observer types.
For example return an observer that observes both the thermostate and the
superdroplets from combination of those two seperate observers */
{
  const ObserveGBxs auto og1 = ObserveTime(stores.timezarr);
  const ObserveGBxs auto og2 = ObserveSDsAttributes(stores.sdzarr);
  const ObserveGBxs auto og3 = ObserveNsupersPerGridBox(stores.nsuperszarr,
                                                        ngbxs);
  const ObserveGBxs auto og4 = create_observegbx_massmoments(stores.massmomszarr,
                                                             ngbxs);
  const auto obsgbxs = og4 >> og3 >> og2 >> og1;

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
  
  /* create process for SDM */
  const auto sdmprocess(create_sdmprocess(config, mdlsteps));

  /* create method for SD motion in SDM */
  const MoveSuperdropsInDomain sdmmotion(NullMotion{});

  /* create observer from combination of chosen observers */
  FSStore fsstore(config.zarrbasedir);
  SomeZarrStores zarrstores(fsstore, config.maxchunk,
                            gbxmaps.ngridboxes, sdattrs_to_observe());
  const auto observer = create_observer(zarrstores,
                                        mdlsteps.obsstep,
                                        gbxmaps.ngridboxes);

  const NullDetectorsPtr dtrs{};

  const RunSDMStep sdm(gbxmaps, sdmmotion, sdmprocess, observer);

  Kokkos::initialize(argc, argv);
  {
    /* RUN SDM MODEL WITH THERMODYNAMICS FROM FILE */
    run_thermofromfile(config, sdm, dtrs,
                       mdlsteps.t_end, mdlsteps.couplstep);
  }
  Kokkos::finalize();
  std::cout << "  ------ Total Duration: " << kokkostimer.seconds() << "s ----- \n";

  return 0;
}
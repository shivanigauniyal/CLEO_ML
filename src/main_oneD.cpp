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
#include <iostream>

#include <Kokkos_Core.hpp>

/* constants and configuration */
#include "claras_SDconstants.hpp"
#include "initialisation/config.hpp"

/* sdm gridboxes setup */
#include "sdmgridboxes/maps4gridboxes.hpp"
#include "sdmgridboxes/sdmtimesteps.hpp"
#include "sdmgridboxes/runsdmstep.hpp"
#include "sdmgridboxes/sdmotion.hpp"
#include "sdmgridboxes/detectors_ptr.hpp"

/* sdm observers setup */
#include "observers/observers.hpp"
#include "observers/observegbxs.hpp"
#include "observers/gridboxes_intostore.hpp"

#include "singlevarstorage.hpp"
#include "massmomentsstorage.hpp"
#include "zarrstorage/sdattributes_intostore.hpp"
#include "zarrstorage/contigraggedsdstorage.hpp"
#include "zarrstorage/thermostatestorage.hpp"
#include "zarrstorage/zarrstores.hpp"

/* sdm superdroplets setup */
#include "superdrop_solver/thermodynamic_equations.hpp"
#include "superdrop_solver/sdmprocess.hpp"
#include "superdrop_solver/coalescence.hpp"
#include "superdrop_solver/condensation.hpp"
#include "superdrop_solver/terminalvelocity.hpp"

/* thermodynamics solver and coupled model setup */
#include "thermofromfile/run_thermofromfile.hpp"

namespace dlc = dimless_constants;

template <SuperdropIntoStoreViaBuffer S>
struct SomeZarrStores
/* structures for outputing data (edit if you are advanced user only) */
{
  ThermoStateStorage thermoz;
  ContiguousRaggedSDStorage<S> sdz;
  ContiguousRaggedSDStorage<SdgbxIntoStore> sdgbxz;
  MomentsStorages mmomsz;
  RainMomentsStorages rainmmomsz;
  CoordinateStorage<double> timez;
  CoordinateStorage<unsigned int> gbxz;
  TwoDStorage<size_t> nsdsz;
  TwoDStorage<size_t> nrainsdsz;

  SomeZarrStores(FSStore &store, const int maxchunk,
                 const unsigned int ngbxs, S sdattrs);
};

SdmProcess auto
create_sdmprocess(const Config &config, const SDMTimesteps &mdlsteps)
/* return an SdmProcess type for condensation and
collision-coalescence in Superdroplet Model */
{
  /* create process for condensation in SDM */
  const double cond_subtstep(realtime2dimless(config.cond_SUBTSTEP));
  const auto cond(CondensationProcess(mdlsteps.condsubstep,
                                      &step2dimlesstime,
                                      config.doAlterThermo,
                                      config.cond_iters,
                                      cond_subtstep,
                                      config.cond_rtol,
                                      config.cond_atol));

  /* create process for collision-coalescence in SDM */
  const auto probs_coal(CollCoalProb_Long());
  const auto coal(CollisionCoalescenceProcess(mdlsteps.collsubstep,
                                              &step2realtime,
                                              probs_coal));
  
  /* combine condensation and collision-coalescence
  processes a single process called 'sdmprocess' */
  const auto sdmprocess = cond >> coal;
  return sdmprocess;
}

SdMotion auto create_sdmotion(const int motionstep)
/* return an SdMotion type for movement of
superdroplets in Superdroplet Model */
{
  /* superdroplet velocity = wind velocity + terminal 
  velocity of droplet from Simmel et al. 2002 formula */
  const auto terminalv(SimmelTerminalVelocity{});
  const SdMotion auto
      movewithsedi = MoveWithSedimentation(motionstep,
                                           &step2dimlesstime,
                                           terminalv);
  return movewithsedi;
}

ObserveGBxs auto
observe_totalmassmoms(MomentsStorages &mms, const size_t ngbxs)
/* choose which moments of the entire droplet
distribution to observe */
{
  const auto mom0(ObserveNthMassMoment(mms.mom0zarr, 0, ngbxs));
  const auto mom1(ObserveNthMassMoment(mms.mom1zarr, 1, ngbxs));

  const auto moms = mom1 >> mom0;
  return moms;
}

ObserveGBxs auto
observe_rainmassmoms(RainMomentsStorages &rmms, const size_t ngbxs)
/* choose which moments of the raindrop distribution
to observe. (Note: raindrop = droplet >40 microns) */
{
  const auto rain0(ObserveNthRainMassMoment(rmms.mom0zarr, 0, ngbxs));
  const auto rain1(ObserveNthRainMassMoment(rmms.mom1zarr, 1, ngbxs));

  const auto rainmoms = rain1 >> rain0;
  return rainmoms;
}

SuperdropIntoStoreViaBuffer auto
sdattrs_to_observe()
/* choose which attributes of a superdroplets to observe */
{
  SuperdropIntoStoreViaBuffer auto id = IdIntoStore();
  SuperdropIntoStoreViaBuffer auto eps = EpsIntoStore();
  SuperdropIntoStoreViaBuffer auto radius = RadiusIntoStore();
  SuperdropIntoStoreViaBuffer auto m_sol = M_solIntoStore();
  SuperdropIntoStoreViaBuffer auto coord3 = Coord3IntoStore();

  return id >> eps >> radius >> m_sol >> coord3;
}

template <SuperdropIntoStoreViaBuffer S>
Observer auto
create_observer(SomeZarrStores<S> &stores,
                const int obsstep,
                const size_t ngbxs)
/* return an Observer type from an amalgamation of observers.
For example, obs = obs1 >> obs2 returns a single observer
that does combination of the observers 'obs1' and 'obs2' */
{
  const auto og1 = ObserveTime(stores.timez);
  const auto og2 = ObserveGridBoxIndex(stores.gbxz);
  const auto og3 = ObserveThermoState(stores.thermoz);

  const auto og4 = observe_totalmassmoms(stores.mmomsz, ngbxs);
  const auto og5 = observe_rainmassmoms(stores.rainmmomsz, ngbxs);
  const auto og6 = ObserveNsupersPerGridBox(stores.nsdsz, ngbxs);
  const auto og7 = ObserveNRainsupersPerGridBox(stores.nrainsdsz, ngbxs);

  const auto og8a = ObserveSDsAttributes(stores.sdz);
  const auto og8b = ObserveSDsGbxindex(stores.sdgbxz);

  const ObserveGBxs auto ogs = og8b >> og8a >>
                                   og7 >> og6 >> og5 >> og4 >>
                                   og3 >> og2 >> og1;
  const Observer auto obs1 = ConstIntervalGBxsObserver(obsstep, ogs);
  const Observer auto obs2 = PrintObserver(obsstep);

  const Observer auto observer = obs2 >> obs1;
  // const Observer auto observer = obs1;
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
  const MoveSuperdropsInDomain
      sdmmotion(create_sdmotion(mdlsteps.motionstep));

  /* create observer from combination of chosen observers */
  const NullDetectorsPtr dtrs{};
  FSStore fsstore(config.zarrbasedir);
  SomeZarrStores zarrstores(fsstore,
                            config.maxchunk,
                            gbxmaps.ngridboxes,
                            sdattrs_to_observe());
  const auto observer = create_observer(zarrstores,
                                        mdlsteps.obsstep,
                                        gbxmaps.ngridboxes);

  /* run SDM with thermodynamics from file according to
  the sdm motion, process and observer constructed above */
  const RunSDMStep sdm(gbxmaps, sdmmotion, sdmprocess, observer);
  Kokkos::initialize(argc, argv);
  {
    run_thermofromfile(config, sdm, dtrs,
                       mdlsteps.t_end, mdlsteps.couplstep);
  }
  Kokkos::finalize();
  std::cout << "  ------ Total Duration: " << kokkostimer.seconds() << "s ----- \n";

  return 0;
}

template <SuperdropIntoStoreViaBuffer S>
SomeZarrStores<S>::SomeZarrStores(FSStore &store, const int maxchunk,
                 const unsigned int ngbxs, S sdattrs)
/* constructor for outputing data (edit if you are advanced user only) */
      : thermoz(store, maxchunk, ngbxs),
        sdz(store, sdattrs, maxchunk),
        sdgbxz(store, SdgbxIntoStore(), maxchunk),
        mmomsz(store, maxchunk, ngbxs),
        rainmmomsz(store, maxchunk, ngbxs),
        timez(make_timezarr(store, maxchunk)),
        gbxz(make_gbxzarr(store, maxchunk)),
        nsdsz(make_nsuperszarr(store, maxchunk, ngbxs)),
        nrainsdsz(make_nrainsuperszarr(store, maxchunk, ngbxs)) {}
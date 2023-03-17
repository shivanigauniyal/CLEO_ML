// Author: Clara Bayley
// File: run_coupledmodel.hpp
/* Header file for functions that run SDM coupled to
CVODE ODE thermodynamics solver. Requires also it's supplement:
'run_coupledmodel_implement.hpp' which is just 2nd half
of this header file that's been moved into a new file for clarity */

#ifndef RUN_COUPLEDMODEL_HPP
#define RUN_COUPLEDMODEL_HPP

#include "run_coupledmodel_implement.hpp"

void run_cvodeSDM_coupledmodel(const Config &config,
                               const Timesteps &mdlsteps,
                               const Maps4GridBoxes mdlmaps,
                               const SdmProcess auto &sdmprocess,
                               const Observer auto &observer)
/* create CVODE thermodynamics solver, superdroplets and gridboxes and
then run superdroplet model (SDM) coupled to the thermodynamics solver */
{
  /* CVODE thermodynamics solver */
  const unsigned int ngridboxes = mdlmaps.idx2bounds_z.size();
  CvodeThermoSolver cvode(config, init_thermodynamics(ngridboxes, config));

  /* vector containing all superdroplets within a struct that also holds their
  associated gridbox index. (all superdroplets have same solute properties) */
  const auto solute(std::make_shared<const SoluteProperties>());
  std::vector<SuperdropWithGridbox>
      SDsInGBxs = superdrops_from_initSDsfile(config.initSDs_filename,
                                              config.NSUPERS,
                                              config.SDnspace, solute,
                                              mdlmaps.idx2bounds_z);

  /* vector containing all gridboxes that makeup the SDM domain */
  std::vector<GridBox> gridboxes = create_gridboxes(ngridboxes, mdlmaps.idx2vol, 
                                                    SDsInGBxs);

  /* prepare, launch, and end coupled model */
  auto gen = prepare_coupledmodel(mdlsteps, cvode, gridboxes);

  timestep_coupledmodel(mdlsteps, mdlmaps, sdmprocess, observer, config.doCouple,
                        cvode, gen, gridboxes, SDsInGBxs);

  printfinish_coupledmodel();
}

void timestep_coupledmodel(const Timesteps &mdlsteps,
                           const Maps4GridBoxes mdlmaps,
                           const SdmProcess auto &sdmprocess,
                           const Observer auto &observer,
                           const bool doCouple,
                           CvodeThermoSolver &cvode,
                           std::mt19937 &gen,
                           std::vector<GridBox> &gridboxes,
                           std::vector<SuperdropWithGridbox> &SDsInGBxs)
/* timestep coupled model from t=0 to t=tend. Each coupled step is
length 'outstep' and decomposed into 4 parts: 1) start of step (coupled)
2) run SDM step (independent) 3) run CVODE step (independent)
4) proceed to next step (coupled) */
{
  int t_out = 0; // time that is incremented by length 'outstep' between coupling communication
  while (t_out <= mdlsteps.tend)
  {
    /* begin coupled step */
    const std::vector<ThermoState> previousstates = start_coupledstep(observer,
                                                                      gridboxes,
                                                                      cvode);
    /* advance SDM by outstep (parallel to CVODE section) */
    run_sdmstep(t_out, mdlsteps.outstep, mdlsteps.xchangestep,
                sdmprocess, mdlmaps, gen, gridboxes, SDsInGBxs);

    /* advance CVODE solver by outstep (parallel to SDM) */
    cvode.run_cvodestep(timestep2dimlesstime(t_out + mdlsteps.outstep));

    /* prepare for next coupled step */
    t_out = proceeed_tonext_coupledstep(t_out, mdlsteps.outstep, doCouple,
                                        previousstates, gridboxes, cvode);
  }
}

#endif // RUN_COUPLEDMODEL_HPP
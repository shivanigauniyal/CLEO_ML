// Author: Clara Bayley
// File: coal_breakup.hpp
/* Header file for class that enacts
collision events in which either coalescence,
or breakup occurs with fixed efficiency
Ecoal (and 1 - Ecoal). ConstCoalBreakup
struct satisfies SDPairEnactX concept used in
CollisionX struct */

#ifndef COAL_BREAKUP_HPP
#define COAL_BREAKUP_HPP

#include <functional>
#include <concepts>
#include <stdexcept>

#include "../claras_SDconstants.hpp"
#include "./superdrop.hpp"
#include "./terminalvelocity.hpp"
#include "./collisionxkernels.hpp"
#include "./collisionx.hpp"
#include "./coalescence.hpp"
#include "./breakup.hpp"

class CoalBreakupConstEff
/* class is method for coalescence / breakup between
two superdroplets with constant efficiency of
coalescence coaleff. (Can be used in collisionsx
struct to enact collision-coalescence or
collision-breakup events in the superdroplet model) */
{
private:
  Coalescence coal;
  Breakup breakup;
  double coaleff;
  double bueff;

public:
  CoalBreakupConstEff(const double infrags, const double coaleff)
      : coal(Coalescence{}), breakup(Breakup(infrags)),
        coaleff(coaleff), bueff(1.0 - coaleff)
  {
    if ((coaleff > 1.0) || (coaleff < 0.0))
    {
      throw std::invalid_argument("Invalid coalescence efficiency, coaleff");
    }
  }

  void operator()(Superdrop &drop1, Superdrop &drop2,
                  const double probcoll, const double phi) const
  /* this operator is used as an "adaptor" for using
  CoalBreakupConstEff as a function in CollisionsX
  that satistfies the SDPairEnactX concept.
  *note* operator uses probcoll, probability of collision,
  NOT probability of collision-coalescence! */
  {
    /* 1. calculate gamma factor for collision-coalescence  */
    const double probcoal(probcoll * coaleff);
    const unsigned long long gamma_coal(coal.coalescence_gamma(drop1.eps,
                                                               drop2.eps,
                                                               probcoal,
                                                               phi));
    /* 2. enact collision-coalescence between pair
      of superdroplets if gamma is not zero */
    if (gamma_coal != 0)
    {
      coal.coalesce_superdroplet_pair(drop1, drop2, gamma_coal);
    }

    else // if not coalescence, check for breakup
    {
      /* 3. calculate gamma factor for collision-breakup  */
      const double probbu(probcoll * bueff);
      const unsigned long long gamma_bu(breakup.breakup_gamma(drop1.eps,
                                                              drop2.eps,
                                                              probbu,
                                                              phi));
      /* 4. enact collision-breakup between pair
        of superdroplets if gamma is not zero */
      if (gamma_bu != 0)
      {
        breakup.breakup_superdroplet_pair(drop1, drop2);
      }
    }
  }
};

SdmProcess auto
CollisionDeJongValid(const int interval,
                     const std::function<double(int)> int2time,
                     const double nfrags,
                     const double coaleff)
/* SDM process for collisions of superdroplets
followed by coalescence or breakup with constant
coaleff (similar to de Jong et al. 2023 sect. 3) */
{
  const double realtstep = int2time(interval);
  const auto terminalv(SimmelTerminalVelocity{});
  const auto collprob(HydrodynamicProb(LongKernelEff{1.0}, terminalv));

  CollisionX<HydrodynamicProb<LongKernelEff, SimmelTerminalVelocity>,
             CoalBreakupRebound<SimmelTerminalVelocity>>
      collall(realtstep, collprob, CoalBreakupConstEff(nfrags, coaleff));

  return ConstTstepProcess{interval, collall};
};

#endif // COAL_BREAKUP_HPP
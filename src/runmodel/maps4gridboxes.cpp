// Author: Clara Bayley
// File: "maps4gridboxes.cpp"
/* functionality for creating and using
the map between a grid box index
and its coordinate boundaries */

#include "maps4gridboxes.hpp"

Maps4GridBoxes::Maps4GridBoxes(const unsigned int SDnspace,
                               std::string_view gridfile)
/* initilaises idx2bounds_[i] maps (for i = x, y or z) which map
from every gridbox index to its boundaries in domain coordinates.
Also initialises idx2vol map whose values are the volume of a gridbox
given the gridbox's index as key. The keys of idx2bounds_[i] map's
are also gridbox indexes. The corresponding value is that gridbox's
{upper boundary, lower boundary}. In a non-3D case, coordinates of the
gridbox boundaries for unused dimensions are the min/max possible
doubles of computer (numerical limits), however the volume remains
finite. E.g. In the 0-D case, the idx2bounds maps have 1 {key, value}
for gridbox 0 which are the upper and lower numerical limits,
whilst the volume is determind by reading the gridfile */
{
  if (SDnspace == 0)
  {
    const double domainvol = domainvol_from_gridfile(gridfile);
    set_0Dmodel_maps(domainvol);
  }

  else if (SDnspace == 1)
  {
    const GridBoxBoundaries gbxbounds(read_gbxboundaries(gridfile));
    set_1Dmodel_maps(gbxbounds);
  }

  else
  {
    const std::string errmsg("SDnspace > 1, no method exists "
                             " for constructing Maps4GridBoxes object");
    throw std::invalid_argument(errmsg);
  }
}

void Maps4GridBoxes::set_0Dmodel_maps(const double domainvol)
/* set idx2bounds_[i] maps to numeical limits. Set volume
 map using coords read from gridfile */
{
  idx2bounds_z[0] = {std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max()};

  idx2bounds_x[0] = {std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max()};

  idx2bounds_y[0] = {std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max()};

  idx2vol[0] = domainvol; // dimensionless volume of 0D model
}

void Maps4GridBoxes::set_1Dmodel_maps(const GridBoxBoundaries &gbxbounds)
/* Set idx2bounds x and y maps to numerical limits. Set z and volume maps using coords
from gridfile. It is assumed zhalf coords read from gbxbounds are 
monotonically decreasing. */
{
  idx2bounds_x[0] = {std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max()};

  idx2bounds_y[0] = {std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max()};

  unsigned int i = 0;
  for (auto it = gbxbounds.zhalf.begin() + 1; it != gbxbounds.zhalf.end(); ++it)
  {
    const double zup = *(it - 1);
    const double zlow = *it;
    idx2bounds_z[i] = {zup, zlow};

    const double vol = (zup - zlow) * gbxbounds.domainarea();
    idx2vol[i] = vol;
    ++i;
  }
}
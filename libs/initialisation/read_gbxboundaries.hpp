// Author: Clara Bayley
// File: read_gbxboundaries.hpp
/* Header for initialisation of Maps4GridBoxes
struct from binary file */

#ifndef READ_GBXBOUNDARIES_HPP
#define READ_GBXBOUNDARIES_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string_view>

#include "./readbinary.hpp"

struct GridBoxBoundaries
/* holds vectors containing z, x and y half coords
(ie. gridbox boundaries) which are read from
gridfile and used in construction of Maps4GridBoxes */
{
  std::vector<double> zhalf;
  std::vector<double> xhalf;
  std::vector<double> yhalf;

  double domainarea() const
  /* returns horizontal area of entire domain */
  {
    const double xdelta = std::abs(xhalf.front() - xhalf.back());
    const double ydelta = std::abs(yhalf.front() - yhalf.back());

    return xdelta * ydelta;
  }

  double domainvol() const
  /* returns volume of entire domain */
  {
    const double zdelta = std::abs(zhalf.front() - zhalf.back());

    return zdelta * domainarea();
  }
};

GridBoxBoundaries read_gbxboundaries(std::string_view gridfile);
/* read metadata and data in binary file called 'gridfile', then
return GridBoxBoundaries instance created from that data */

inline double domainvol_from_gridfile(std::string_view gridfile)
/* return the volume of the domain determined by
reading data from 'gridfile' binary */
{
  const GridBoxBoundaries gbxbounds(read_gbxboundaries(gridfile));

  return gbxbounds.domainvol();
}

#endif // READ_GBXBOUNDARIES_HPP 
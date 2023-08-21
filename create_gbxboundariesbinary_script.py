import sys
import numpy as np
from pathlib import Path

from pySD.gbxboundariesbinary_src.create_gbxboundaries import *
from pySD.gbxboundariesbinary_src.read_gbxboundaries import *

### ----------------------- INPUT PARAMETERS ----------------------- ###
### Essential Absolute or Relative Paths and Filenames ###
path2CLEO     = sys.argv[1]
path2build    = sys.argv[2]
configfile    = sys.argv[3]
constsfile    = path2CLEO+"libs/claras_SDconstants.hpp"
binariespath  = path2build+"/share/"
savefigpath   = path2build+"/bin/"
gridfile      = binariespath+"/dimlessGBxboundaries.dat" # Note: this should match config_oneD.txt

### Booleans for Plotting Figures of Initial Conditions ###
isfigures      = [True, True]                               # [making+showing, saving]

### input parameters for coordinates of gridbox boundaries
zmax          = 5000     # maximum z coord [m]
zmin          = 0        # minimum z coord [m]
zdelta        = 20     # uniform gridbox spacing [m]
zgrid         = [zmin, zmax, zdelta] 

xgrid         = np.asarray([0, 20]) # [m]
ygrid         = np.asarray([0, 20]) # [m]
### ---------------------------------------------------------------- ###

### -------------------- BINARY FILE GENERATION--------------------- ###
### ensure build, share and bin directories exist ###
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True) 
  Path(binariespath).mkdir(exist_ok=True) 

### write gridbox boundaries binary
write_gridboxboundaries_binary(gridfile, zgrid, 
                               xgrid, ygrid, constsfile)
print_domain_info(constsfile, gridfile)

### plot gridbox boundaries binary
if isfigures[0]:
  if isfigures[1]:
    Path(savefigpath).mkdir(exist_ok=True) 
  plot_gridboxboundaries(constsfile, gridfile,
                         savefigpath, isfigures[1])
### ---------------------------------------------------------------- ###
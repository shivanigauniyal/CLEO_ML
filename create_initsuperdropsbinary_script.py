import sys
import numpy as np
from pathlib import Path

from pySD.initsuperdropsbinary_src import initattributes as iSDs
from pySD.initsuperdropsbinary_src import radiiprobdistribs as rprobs
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers 
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers 

### ----------------------- INPUT PARAMETERS ----------------------- ###
### Essential Absolute or Relative Paths and Filenames ###
path2CLEO     = sys.argv[1]
path2build    = sys.argv[2]
configfile    = sys.argv[3]
constsfile     = path2CLEO+"libs/claras_SDconstants.hpp"
binariespath   = path2build+"/share/"
savefigpath    = path2build+"/bin/"
gridfile       = binariespath+"/dimlessGBxboundaries.dat"   # Note: this should match config_oneD.txt
initSDsfile    = binariespath+"/dimlessSDsinit.dat"         # Note: this should match config_oneD.txt

### Booleans for Plotting Figures of Initial Conditions ###
isfigures      = [True, True]                               # [making+showing, saving]
gbxs2plt       = [0]                                        # GBx index of SDs to plot (nb. "all" can be very slow)

### Number of Superdroplets (an int or dict of ints) ###
npergbx = 8                                         # number of SDs per Gridbox (for gridboxes with SDs)
zlim    = 4000                                      # only have SDs in gridboxes with lower boundary > zlim [m]
nsupers = iSDs.nsupers_at_domain_top(gridfile, constsfile, npergbx, zlim)
# nsupers = 1024                                    # all GBxs have this many SDs initially

### Choice of Superdroplet Radii Generator ###
rspan                = [1e-8, 9e-5]                 # min and max range of radii to sample [m]
randomr              = True                         # sample radii range randomly or not
radiigen = iSDs.SampleDryradiiGen(rspan, randomr)   # radii are sampled from rspan [m]

### Choice of Droplet Radius Probability Distribution ###      
geomeans             = [0.02e-6, 0.2e-6, 3.5e-6]    # mean radius of each mode [m]        
geosigs              = [1.55, 2.3, 2]               # logarithmic sigma of each mode
scalefacs            = [1, 0.3, 0.025]              # relative weights of each mode
numconc              = np.sum(scalefacs)*1e8        # [number of droplets / m^3] (in gridboxes with superdroplets initially)
radiiprobdist = rprobs.LnNormal(geomeans, geosigs, scalefacs)

### Choice of Superdroplet Spatial Coordinates Generator ###
coord3gen            = iSDs.SampleCoordGen(True)   # sample coord3 range randomly or not
coord1gen            = None                        # do not generate superdroplet coord1s
coord2gen            = None                        # do not generate superdroplet coord2s
### ---------------------------------------------------------------- ###


### -------------------- BINARY FILE GENERATION--------------------- ###
### ensure build, share and bin directories exist
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True) 
  Path(binariespath).mkdir(exist_ok=True) 

### write initial superdrops binary
initattrsgen = iSDs.InitManyAttrsGen(radiigen, radiiprobdist,
                                      coord3gen, coord1gen, coord2gen)
csupers.write_initsuperdrops_binary(initSDsfile, initattrsgen, 
                                    configfile, constsfile,
                                    gridfile, nsupers, numconc)

### plot initial superdrops binary
if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True) 
    rsupers.plot_initGBxsdistribs(configfile, constsfile, initSDsfile,
                                   gridfile, savefigpath, isfigures[1],
                                   gbxs2plt)
### ---------------------------------------------------------------- ###
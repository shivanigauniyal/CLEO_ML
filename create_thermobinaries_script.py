import sys
import numpy as np
from pathlib import Path

from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

### ----------------------- INPUT PARAMETERS ----------------------- ###
### Essential Absolute or Relative Paths and Filenames ###
path2CLEO     = sys.argv[1]
path2build    = sys.argv[2]
configfile    = sys.argv[3]
constsfile    = path2CLEO+"libs/claras_SDconstants.hpp"
binariespath  = path2build+"/share/"
savefigpath   = path2build+"/bin/"
gridfile      = binariespath+"/dimlessGBxboundaries.dat"    # Note: this should match config_oneD.txt
thermofiles   =  binariespath+"/dimlessthermo.dat"          # note this should be consitent with config_oneD.txt

### Booleans for Plotting Figures of Initial Conditions ###
isfigures     = [True, True]                                # [making+showing, saving]

### Choose Initial Thermodynamic Conditions for Gridboxes ###
PRESS0        = 101500                                      # [Pa]
THETA         = 289                                         # [K]
qvap          = "sratio"                                    # key for calculating initial vapour mass mixing ratio
Zbase         = 2000                                        # [m]
relhratios    = [0.95, 1.0001]                              # relative humidity [below, above] Zbase
moistlayer    = False                                       # no moist layer at Zbase
qcond         = 0                                           # [Kg/Kg]
WMAX          = None                                        # [m/s]
VVEL          = None                                        # [m/s]
Zlength       = 0                                           # [m]
Xlength       = 0                                           # [m]
thermodyngen = thermogen.ConstHydrostaticAdiabat(configfile, constsfile,
                                                 PRESS0, THETA, qvap,
                                                 relhratios, Zbase,
                                                 qcond, WMAX, Zlength,
                                                 Xlength, VVEL, moistlayer)
### ---------------------------------------------------------------- ###

### -------------------- BINARY FILE GENERATION--------------------- ###
### ensure build, share and bin directories exist ###
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True) 
  Path(binariespath).mkdir(exist_ok=True) 
  
cthermo.write_thermodynamics_binary(thermofiles, thermodyngen, configfile,
                                    constsfile, gridfile)

if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True) 
    rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
                                          thermofiles, savefigpath,
                                          isfigures[1])
### ---------------------------------------------------------------- ###
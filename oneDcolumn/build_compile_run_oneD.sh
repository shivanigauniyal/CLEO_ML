#!/bin/bash
#SBATCH --job-name=runoneD
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=/work/mh1126/m300950/breakup/build/tmp/runoneD_out.%j.out
#SBATCH --error=/work/mh1126/m300950/breakup/build/tmp/runoneD_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/cleoenv 
path2CLEO="/home/m/m300950/CLEO_ML/"
path2build="/home/m/m300950/CLEO_ML/oneDcolumn/build/"
configfile="/home/m/m300950/CLEO_ML/oneDcolumn/config_oneD.txt"
### ---------------------------------------------------- ###

### ------------------- compile_run.sh ----------------- ###
### build CLEO using cmake (with openMP thread parallelism through Kokkos)
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON"  # openMP parallelism enabled
CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags}

### ensure these directories exist
mkdir ${path2build}bin
mkdir ${path2build}share
mkdir ${path2build}tmp

### generate config and input files within build directory
tmpconfig="${path2build}/tmp/config.txt"
echo cp ${configfile} ${tmpconfig}
cp ${configfile} ${tmpconfig}
${python} ${path2CLEO}/oneDcolumn/create_oneD_inputbinaries.py ${path2CLEO} ${path2build} ${tmpconfig}

### compile CLEO in build directory
cd ${path2build} && pwd 
make runoneD -j 16

### run CLEO
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
constsfile="${path2CLEO}libs/claras_SDconstants.hpp"
runcmd="${path2build}/src/runoneD ${tmpconfig} ${constsfile}"
echo ${runcmd}
${runcmd}
### ---------------------------------------------------- ###
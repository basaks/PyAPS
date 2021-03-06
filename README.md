# PyAPS 
This repo starts from the PyAPS code I could download from google code archive from here: https://code.google.com/archive/p/pyaps/downloads

More recent code may exist here: http://earthdef.caltech.edu/

## Creating a conda `virtualenv` for `PyAPS`

####  Install anaconda:

On 64bit linux:
    
    bash Anaconda2-2.4.1-Linux-x86_64.sh

Please accept all defaults. This will install anaconda in the `~/anaconda2` directory.
 
You will need to download and install the appropriate version for your OS. The following instructions are for ubuntu 14.04 and are exptect to work for most versions of ubuntu.

#### Install `conda-env`

    source ~/anaconda2/bin/activate anaconda2/
    conda install -c conda conda-env

#### Create the `pyaps` environment
    
    conda env create -f environment.yml
    ~/anaconda2/bin/conda env create -f path/to/PyAPS/environment.yml

#### Activate the `pyaps` environment

    source ~/anaconda2/bin/activate pyaps
    
#### Back to main anaconda if you need to:
    
    source activate /home/user/anaconda2
    
#### Deactivate
Once you are done using `pyaps` you could deactivate from the conda env using: 

    source deactivate

## Run eraptr2rdr.py

This is how I can run the `eraptr2rdr.py`: 

    python bin/eraptr2rdr.py examples/20061016.dat examples/dem_16rlks.hgt output 5.6 23
    
This creates the `output` file.

#### Explanation
Copied from the PyAPS.pdf I download with the original code:


    The eraptr2rdr.py is the exact translation of the original codes from Jolivet
    et al 2011 into python. However, the formatting of inputs is slightly different:
        
        > eraptr2rdr.py inputfile dem.hgt outname wvl inc
    
    The inputs to the program are:
    
    1. inputfile: The input file is a text file containing all the data for all the needed
    grid points. Special formatting is required. An example is shown in
    the example directory (20061016.dat).
    
    2. dem.hgt: Radar simulation in radar geometry as produced by ROI PAC.
    An example is shown in the example directory (`dem_16rlks.hgt`). The
    associated `rsc` file is also needed.
    
    3. outname: Name of the output file. The output is a binary file (float32) that
    has the same width and length as the `dem.hgt` file.
    
    4. wvl: Radar wavelength in cm.
    
    5. inc: Radar signal incidence angle, with respect to the vertical (ex: Envisat∼23 degrees).


## Run era2rdr.py

This is how I can run the `era2rdr.py`:

    python bin/era2rdr.py examples/ei.oper.an.pl.regn128sc.2006101606 examples/dem_16rlks.hgt output 5.6 23
    
This creates the `output` file.

Please note that the file `ei.oper.an.pl.regn128sc.2006101606` is not included in this repo as it is rather large (>100MB).

#### Explanation
Copied from the PyAPS.pdf I download with the original code:


    The era2rdr.py routine computes a delay map, in radar geometry, directly from the GRIB file given by ERA-Interim. Such input files can be found at http://dss.ucar.edu/datasets/ds627.0/.
        > era2rdr.py inputfile dem.hgt outname wvl inc
    
    The inputs to the program are:
    
    1. inputfile The input file is a solution GRIB file provided by ERA-I. Usual name format is ei.oper.an.pl.regn128sc.yyyymmddhh.
    
    2. dem.hgt: Radar simulation in radar geometry as produced by ROI PAC.
    An example is shown in the example directory (`dem_16rlks.hgt`). The
    associated `rsc` file is also needed.
    
    3. outname: Name of the output file. The output is a binary file (float32) that
    has the same width and length as the `dem.hgt` file.
    
    4. wvl: Radar wavelength in cm.
    
    5. inc: Radar signal incidence angle, with respect to the vertical (ex: Envisat∼23 degrees).




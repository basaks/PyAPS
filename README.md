# PyAPS 
This repo starts from the PyAPS code I could download from google code archive from here: https://code.google.com/archive/p/pyaps/downloads

More recent code may exist here: http://earthdef.caltech.edu/


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



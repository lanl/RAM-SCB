#!/bin/tcsh

if (! -f config.txt) then
   echo 'Please create a config.txt in this directory'
   exit 1
endif

if (! -d vtk_files) then
	mkdir vtk_files
        python makeCustomSource.py
endif

if (! -d images) then
   mkdir images
endif

module load paraview
python convertRAMrestart.py nc_files #populates vtk_files
pvpython visualizeRAM.py #populates images
exit 0

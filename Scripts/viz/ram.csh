#!/bin/tcsh

if (! -f config.txt) then
   echo 'Please create a config.txt in this directory'
   exit 1
endif

if (-d vts_files) then
	rm vts_files/restart_*
else
	mkdir vts_files
        python makeCustomSource.py
endif

if (! -d images) then
   mkdir images
endif

module load paraview
python convertRAMrestart.py nc_files #populates vts_files
pvpython visualizeRAM.py #populates images
exit 0

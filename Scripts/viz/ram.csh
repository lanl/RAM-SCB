#!/bin/tcsh

if (! -f config.txt) then
   echo 'Please create a config.txt in this directory'
   exit 1
endif

if (-d vts_files) then
	rm vts_files/*
else
	mkdir vts_files
endif

if (! -d images) then
   mkdir images
endif

module load paraview
python ram_automate1.py nc_files #populates vts_files
pvpython ram_automate2.py #populates images
exit 0

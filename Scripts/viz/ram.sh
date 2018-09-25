#!/bin/bash

if [[ ! -f 'config.txt' ]]
then
   echo 'Please create a config.txt in this directory'
   exit 1
fi

if [[ ! -d 'vtk_files' ]]
then
	mkdir vtk_files
        python makeCustomSource.py
fi

if [[ ! -d 'images' ]]
then
   mkdir images
fi

module load paraview
python convertRAMrestart.py nc_files #populates vtk_files
pvpython visualizeRAM.py #populates images
exit 0

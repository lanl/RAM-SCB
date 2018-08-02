#!/bin/bash

if [[ ! -f 'config.txt' ]]
then
   echo 'Please create a config.txt in this directory'
   exit 1
fi

if [[ -d 'vts_files' ]]
then
	rm vts_files/*
else
	mkdir vts_files
fi

if [[ ! -d 'images' ]]
then
   mkdir images
fi

module load paraview
python ram_automate1.py nc_files #populates vts_files
pvpython ram_automate2.py #populates images
rm -r vts_files
exit 0

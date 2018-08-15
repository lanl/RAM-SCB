'''Creates the .vts files in vts_files directory for generating RAM-SCB visualizations
   1) Generates fn_field.vts and fn_pressure.vts
	2) Generates the sphere source to get streamlines, if streamline source is sphere
	Usage: python ram_automate1.py <directory_with_NETCDF_files>
	Input: directory containing .nc files'''

import os, sys, itertools
import numpy as np
import spacepy.datamodel as dm

#=================================================================================================
def gen_sphere():
	'''Generating the sphere source for streamlines as sphere.vts in vts_files directory'''
	n = 1 #how much denser points on sphere are, compared to the magnetic field
	#----------------------------generate the points describing the sphere--------------------------
	r = 1.1
	points = []

	Ls = np.arange(2,7.0,0.25)
        Lats = np.arccos(np.sqrt(1./Ls)) #N.Hem latitudes in radians
        Lats = np.concatenate([-1*Lats[::-1], Lats]) #add S.Hem
        Colats = (np.pi/2.)-Lats #Colatitudes in radians
        Lons = np.deg2rad(np.arange(0,360,20))
        npts = len(Lats)*len(Lons)
        xyz = np.zeros([npts,3])
	for idx, LatLon in enumerate(itertools.product(Lons,Colats)):
        	xyz[idx,0] = r*np.cos(LatLon[0])*np.sin(LatLon[1])
                xyz[idx,1] = r*np.sin(LatLon[0])*np.sin(LatLon[1])
                xyz[idx,2] = r*np.cos(LatLon[1])

	#----------------------------generate the sphere.vts file--------------------------
	to_write = ''

	to_write += '<?xml version="1.0"?>\n'  #1
	to_write += '<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian">\n' #2
	to_write += '\t<PolyData>\n' #3
	to_write += '\t\t<Piece NumberOfPoints="{}" NumberOfVerts="{}" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">\n'.format(npts, 1) #4

	to_write += '\t\t\t<PointData>\n' #5
	to_write += '\t\t\t</PointData>\n' #9
	to_write += '\t\t\t<CellData>\n'
	to_write += '\t\t\t</CellData>\n'

	to_write += '\t\t\t<Points>\n'
	to_write += '\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'
	for p in xyz:#points:
		to_write += '\t\t\t\t\t' + str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2]) + '\n'
	to_write += '\t\t\t\t</DataArray>\n'
	to_write += '\t\t\t</Points>\n'

	to_write += '<Verts>\n'
        to_write += '<DataArray type="Int64" Name="connectivity" RangeMin="0" RangeMax="{}">\n'.format(npts-1)
        for i in range(npts): to_write += '{} '.format(i)
        to_write += '\n</DataArray>\n'
        to_write += '<DataArray type="Int64" Name="offsets" RangeMin="{}" RangeMax="{}">\n'.format(npts, npts)
        to_write += '{}\n'.format(npts)
        to_write += '</DataArray>\n'
	to_write += '</Verts>\n'
	to_write += '<Lines>\n'
        to_write += '<DataArray type="Int64" Name="connectivity" RangeMin="0" RangeMax="{}">\n'.format(npts)
        to_write += '</DataArray>\n'
        to_write += '<DataArray type="Int64" Name="offsets" RangeMin="{}" RangeMax="{}">\n'.format(npts, npts)
        to_write += '</DataArray>\n'
	to_write += '</Lines>\n'
	to_write += '<Strips>\n'
        to_write += '<DataArray type="Int64" Name="connectivity" RangeMin="0" RangeMax="{}">\n'.format(npts)
        to_write += '</DataArray>\n'
        to_write += '<DataArray type="Int64" Name="offsets" RangeMin="{}" RangeMax="{}">\n'.format(npts, npts)
        to_write += '</DataArray>\n'
	to_write += '</Strips>\n'
	to_write += '<Polys>\n'
        to_write += '<DataArray type="Int64" Name="connectivity" RangeMin="0" RangeMax="{}">\n'.format(npts)
        to_write += '</DataArray>\n'
        to_write += '<DataArray type="Int64" Name="offsets" RangeMin="{}" RangeMax="{}">\n'.format(npts, npts)
        to_write += '</DataArray>\n'
	to_write += '</Polys>\n'

	to_write += '\t\t</Piece>\n'
	to_write += '\t</PolyData>\n'
	to_write += '</VTKFile>\n'

	with open('vts_files/sphere.vtp', 'w') as fh:
		fh.write(to_write)
	#return xyz, Lats, Lons

#=================================================================================================
if __name__ == '__main__':
	gen_sphere() #generates sphere.vts in vts_files directory

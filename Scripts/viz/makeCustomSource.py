#!/usr/bin/env python3
'''Creates the .vts files in vts_files directory for generating RAM-SCB visualizations
   1) Generates fn_field.vts and fn_pressure.vts
    2) Generates the sphere source to get streamlines, if streamline source is sphere
    Usage: python ram_automate1.py <directory_with_NETCDF_files>
    Input: directory containing .nc files'''

import os, sys, itertools, io
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
    Lons = np.deg2rad(np.arange(0,360,15))
    npts = len(Lats)*len(Lons)
    xyz = np.zeros([npts,3])
    for idx, LatLon in enumerate(itertools.product(Lons,Colats)):
        xyz[idx,0] = r*np.cos(LatLon[0])*np.sin(LatLon[1])
        xyz[idx,1] = r*np.sin(LatLon[0])*np.sin(LatLon[1])
        xyz[idx,2] = r*np.cos(LatLon[1])

    #----------------------------generate the sphere.vtx file--------------------------
    to_write = xmlPolyGen(xyz)
    with open('vts_files/sphere.vtp', 'wb') as fh:
        fh.write(to_write)


def xmlPolyGen(xyz):
    import lxml.etree
    npts = len(xyz)

    #Set up element tree for XML
    fulltree = lxml.etree.ElementTree(lxml.etree.Element('VTKFile', type='PolyData', version='1.0', byte_order='LittleEndian'))
    parent = fulltree.getroot()
    #Everything goes into PolyData as it's a PolyData filetype
    polydata = lxml.etree.Element('PolyData')
    parent.append(polydata)
    #All points are in a single "Piece"
    piece = lxml.etree.Element('Piece',  NumberOfPoints='{}'.format(npts), NumberOfVerts='1', NumberOfLines='0', NumberOfStrips='0', NumberOfPolys='0')
    polydata.append(piece)
    #No PointData or CellData as we just want point locations (with no associated values)
    pointdata = lxml.etree.Element('PointData')
    celldata = lxml.etree.Element('CellData')
    piece.append(pointdata)
    piece.append(celldata)
    #Now we add the actual point locations
    points = lxml.etree.Element('Points') #points needs a "DataArray" element with the point locations
    xyzlocations = lxml.etree.Element('DataArray', type='Float64', NumberOfComponents='3', format='ascii')
    xyzlocations.text = '\n'+''.join(['\t\t\t{} {} {}\n'.format(a,b,c) for (a,b,c) in xyz.tolist()])
    points.append(xyzlocations)
    piece.append(points)
    #PolyData requires other elements: Verts, Lines, Strips, Polys. These can mostly be empty.
    verts = lxml.etree.Element('Verts')
    verts_c = lxml.etree.Element('DataArray', type='Int64', Name='connectivity', RangeMax='{}'.format(npts-1), RangeMin='0')
    verts_o = lxml.etree.Element('DataArray', type='Int64', Name='offsets', RangeMin='{}'.format(npts), RangeMax='{}'.format(npts))
    verts_c.text = '\n' + ' '.join(['{}'.format(i) for i in range(npts)]) + '\n'
    verts_o.text = '\n{}'.format(npts) +'\n'
    verts.append(verts_c)
    verts.append(verts_o)
    piece.append(verts)
    for partname in ['Lines', 'Strips', 'Polys']:
        dac = lxml.etree.Element('DataArray', type='Int64', Name='connectivity', RangeMax='{}'.format(npts), RangeMin='0')
        dao = lxml.etree.Element('DataArray', type='Int64', Name='offsets', RangeMin='{}'.format(npts), RangeMax='{}'.format(npts))
        part = lxml.etree.Element(partname)
        part.append(dac)
        part.append(dao)
        piece.append(part)

    #Print to a string, include XML declaration and "pretty print" to get newlines and indentation.
    out = lxml.etree.tostring(fulltree, xml_declaration=True, pretty_print=True)
    return out


#=================================================================================================
if __name__ == '__main__':
    gen_sphere() #generates sphere.vtp in vts_files directory

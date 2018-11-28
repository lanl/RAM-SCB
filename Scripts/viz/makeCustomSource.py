#!/usr/bin/env python3
'''Creates the .vtp files with seed locations in vtk_files directory 
   for generating RAM-SCB visualizations

   Usage: python makeCustomSource.py
'''

import os, sys, itertools, io
import lxml.etree
import numpy as np
import spacepy.datamodel as dm

#=================================================================================================
def gen_sphere():
    '''Generating a spherical source for streamlines as sphere.vtp in vtk_files directory'''
    #----------------------------generate the points describing the sphere--------------------------
    xyz, nLons, nLats = pointsGenSphere(r=1.5, minL=3.0, maxL=7.0, dL=0.25)

    #----------------------------generate the sphere.vtp file--------------------------
    to_write = xmlPolyGen(xyz, dim1=nLons, dim2=nLats)
    outdir = 'vtk_files'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(os.path.join(outdir, 'sphere.vtp'), 'wb') as fh:
        fh.write(to_write)


def gen_disc():
    '''Generating a disc source for streamlines as disc.vtp in vtk_files directory'''
    #----------------------------generate the points describing the disc--------------------------
    xyz, nLons, nL = pointsGenDisc(minL=2.5, maxL=8.5, dL=0.5)

    #----------------------------generate the disc.vtp file--------------------------
    to_write = xmlPolyGen(xyz, dim1=nLons, dim2=nL)
    outdir = 'vtk_files'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    with open(os.path.join(outdir, 'disc.vtp'), 'wb') as fh:
        fh.write(to_write)


def pointsGenSphere(r=1.1, minL=2.0, maxL=7.0, dL=0.25):
    '''generates Cartesian coordinates for points on a spherical surface given extent/resolution'''
    points = []
    Ls = np.arange(minL, maxL, dL)
    Lats = np.arccos(np.sqrt(1./Ls)) #N.Hem latitudes in radians
    #Lats = np.concatenate([-1*Lats[::-1], Lats]) #add S.Hem
    Colats = (np.pi/2.)-Lats #Colatitudes in radians
    Lons = np.deg2rad(np.arange(0,360,15))
    npts = len(Lats)*len(Lons)
    xyz = np.zeros([npts,3])
    for idx, LonLat in enumerate(itertools.product(Lons,Colats)):
        xyz[idx,0] = r*np.cos(LonLat[0])*np.sin(LonLat[1])
        xyz[idx,1] = r*np.sin(LonLat[0])*np.sin(LonLat[1])
        xyz[idx,2] = r*np.cos(LonLat[1])
    return xyz, len(Lons), len(Lats)


def pointsGenDisc(minL=2.0, maxL=7.0, dL=0.25):
    '''generates Cartesian coordinates for points on an annular disc given extent/resolution'''
    points = []
    Ls = np.arange(minL, maxL, dL)
    Lats = np.zeros(len(Ls))
    Colats = (np.pi/2.)-Lats #Colatitudes in radians
    Lons = np.deg2rad(np.arange(0,360,10))
    npts = len(Lats)*len(Lons)
    xyz = np.zeros([npts,3])
    for idx, LonL in enumerate(itertools.product(Lons,Ls)):
        xyz[idx,0] = LonL[1]*np.cos(LonL[0])
        xyz[idx,1] = LonL[1]*np.sin(LonL[0])
        xyz[idx,2] = 0.0
    return xyz, len(Lons), len(Ls)


def xmlPolyGen(xyz, dim1, dim2):
    '''writes Cartesian coordinates of 2D surface to VTK PolyData XML, including connectivity'''
    npts = len(xyz)

    #convenience handles
    indent3 = '\n\t\t\t'
    indent2 = '\n\t\t'

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
    xyzlocations.text = '\n'+''.join(['\t\t\t{} {} {}\n'.format(a,b,c) for (a,b,c) in xyz.tolist()]) + '\t\t'
    points.append(xyzlocations)
    piece.append(points)
    #PolyData requires other elements: Verts, Lines, Strips, Polys. These can mostly be empty.
    verts = lxml.etree.Element('Verts')
    verts_c = lxml.etree.Element('DataArray', type='Int64', Name='connectivity', RangeMax='{}'.format(npts-1), RangeMin='0')
    verts_o = lxml.etree.Element('DataArray', type='Int64', Name='offsets', RangeMin='{}'.format(npts), RangeMax='{}'.format(npts))
    vertnums = ['{}'.format(i) for i in range(npts)]
    for i in range(npts//8)[::-1]:
        vertnums.insert((i+1)*8, indent3)
    verts_c.text = indent3 + ' ' + ' '.join(vertnums) + indent2
    verts_o.text = indent3+ '{}'.format(npts) + indent2
    verts.append(verts_c)
    verts.append(verts_o)
    piece.append(verts)
    for partname in ['Lines', 'Strips', 'Polys']:
        dac = lxml.etree.Element('DataArray', type='Int64', Name='connectivity', RangeMax='{}'.format(npts), RangeMin='0')
        dao = lxml.etree.Element('DataArray', type='Int64', Name='offsets', RangeMin='{}'.format(0), RangeMax='{}'.format(npts))
        if partname is 'Polys':
            #write the polygons that connect the points together, defining the grid
            polyconn = getPolyVertOrder(len(xyz), dim1, dim2)
            nconn = len(polyconn)
            npolys = dim1*(dim2-1)
            for i in range(nconn//4)[::-1]:
                polyconn.insert((i+1)*4, indent3)
            connstr = ' ' + ' '.join(polyconn)
            #offsets are 4, 8, 12, ...
            polyoffset = ['{}'.format(4*(n+1)) for n in range(nconn//4)]
            for i in range(nconn//4)[::-1]:
                polyoffset.insert((i+1)*4, indent3)
            offstr = ' ' + ' '.join(polyoffset)
            dac.text = indent3 + connstr.rstrip() + indent2
            dao.text = indent3 + offstr.rstrip() + indent2
            dao.set('RangeMax', '{}'.format(nconn))
            piece.set('NumberOfPolys', '{}'.format(npolys))
        part = lxml.etree.Element(partname)
        part.append(dac)
        part.append(dao)
        piece.append(part)

    #Print to a string, EXclude XML declaration and "pretty print" to get newlines and indentation.
    out = lxml.etree.tostring(fulltree, xml_declaration=False, pretty_print=True)
    return out

def getPolyVertOrder(ncoords, dim1, dim2):
    # x, y, z are calculated using itertools.product
    # given (for disc) Lon, L, the loop order is all L for each Lon, iterating over Lon
    # i.e. (Lon1, L1), (Lon1, L2), (Lon1, L3), ..., (Lon2, L1), (Lon2, L2), ...
    cells = []
    for idx1, idx2 in itertools.product(range(dim1-1), range(dim2-1)):
        vert1 = idx1*dim2+idx2
        vert2 = vert1+1
        vert3 = (idx1+1)*dim2+idx2+1
        vert4 = (idx1+1)*dim2+idx2
        cells.append([vert1, vert2, vert3, vert4])
    #add cells that wrap in Longitude
    for idx in range(dim2-1):
        vert1 = idx
        vert2 = idx+1
        vert3 = ncoords-dim2+idx+1
        vert4 = ncoords-dim2+idx
        cells.append([vert1, vert2, vert3, vert4])
    flat_cells = [str(item) for sublist in cells for item in sublist]
    return flat_cells

def getCellVertOrder(ncoords, dim1, dim2, dim3, wrap=False, verbose=False):
    #order is (Lon1, L1, lat1), (Lon1, L1, Lat2), ...
    cells = []
    if wrap:
        ivert1, fvert1 = np.zeros((dim2,dim3)).astype(int), np.zeros((dim2,dim3)).astype(int)
        ivert2, fvert2 = np.zeros((dim2,dim3)).astype(int), np.zeros((dim2,dim3)).astype(int)
        ivert3, fvert3 = np.zeros((dim2,dim3)).astype(int), np.zeros((dim2,dim3)).astype(int)
        ivert4, fvert4 = np.zeros((dim2,dim3)).astype(int), np.zeros((dim2,dim3)).astype(int)
    for idx1, idx2, idx3 in itertools.product(range(dim1-1), range(dim2-1), range(dim3-1)):
        vert1 = idx1*dim2*dim3 +idx2*dim3 + idx3 #(npts to skip to reach LT) + (nFLs to skip to reach FL) + (npts along FL)
        vert2 = vert1+1 #pts 1 and 2 are adjacent in 3rd dimension
        vert3 = idx1*dim2*dim3 +(idx2+1)*dim3 + idx3 +1 #same LT, next FL out in radius
        vert4 = vert3-1 #this is 1 less than vert 3 to get correct cell ordering for VTK hexahedra
        vert5 = (idx1+1)*dim2*dim3 +idx2*dim3 + idx3 #as above, but 1 LT further around
        vert6 = vert5+1
        vert7 = (idx1+1)*dim2*dim3 +(idx2+1)*dim3 + idx3+1
        vert8 = vert7-1
        if wrap and idx1==0:
            ivert1[idx2,idx3] = vert1
            ivert2[idx2,idx3] = vert2
            ivert3[idx2,idx3] = vert3
            ivert4[idx2,idx3] = vert4
        if wrap and idx1==dim1-2:
            fvert1[idx2,idx3] = vert5
            fvert2[idx2,idx3] = vert6
            fvert3[idx2,idx3] = vert7
            fvert4[idx2,idx3] = vert8
        cells.append([vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8])
    #add cells that wrap in Longitude
    if wrap:
        wrap_counter = 0
        for idx2,idx3 in itertools.product(range(dim2-1), range(dim3-1)):
            vert1 = fvert1[idx2,idx3]
            vert2 = fvert2[idx2,idx3]
            vert3 = fvert3[idx2,idx3]
            vert4 = fvert4[idx2,idx3]
            vert5 = ivert1[idx2,idx3]
            vert6 = ivert2[idx2,idx3]
            vert7 = ivert3[idx2,idx3]
            vert8 = ivert4[idx2,idx3]
            cells.append([vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8])
            wrap_counter += 1
        if verbose: print('Adding {} cells at the wraparound'.format(wrap_counter))
    flat_cells = [str(item) for sublist in cells for item in sublist]
    return flat_cells

#=================================================================================================
if __name__ == '__main__':
    gen_sphere() #generates sphere.vtp in vtk_files directory
    gen_disc()

'''Creates the .vtx files in vtk_files directory for generating RAM-SCB visualizations
   1) Generates fn_field.vtu and fn_pressure.vtp
   Usage: python convertRAMrestart.py <directory_with_NETCDF_files>
   Input: directory containing .nc files'''

import os, sys, itertools, glob
import lxml.etree
import numpy as np
import spacepy.datamodel as dm
import netCDF4
from makeCustomSource import getPolyVertOrder, getCellVertOrder

#========================================================================================================
def fromNC4(fname):
    '''Uses netCDF4 Python library to load .nc restart files into SpaceData/dmarray'''
    outdata = dm.SpaceData()
    with netCDF4.Dataset(fname, "r", format="NETCDF4") as ncf:
        #global attrs
        for attname in ncf.ncattrs():
            outdata.attrs[attname] = ncf.getncattr(attname)
        #get copies of all variables
        for varname in ncf.variables:
            outdata[varname] = dm.dmarray(ncf.variables[varname][...].filled())
            #RAM-SCB variables don't have metadata, so skip
            #RAM-SCB also doesn't have nested groups, so don't recurse
        #and add netCDF4 "dimensions" as variables
        for dimname in ncf.dimensions:
            outdata[dimname] = dm.dmarray(np.atleast_1d(ncf.dimensions[dimname].size))
    return outdata

#========================================================================================================
def read_config():
    '''Reads configurations from config.txt'''
    global properties
    property_labels = ['Streamline source type', 'Plasma pressure display', 'Camera Position',  
                       'Camera Focal Point', 'Camera View Up', 'Scale', 'Movie', 'SaveState']
    with open('config.txt', 'r') as to_read:
        lines = to_read.readlines()
    if len(lines) - 1 != len(property_labels):
        raise ValueError('Number of properties in the config file should be ' + len(properties))
        sys.exit(1)
    lines = [x[x.find(':')+1:].strip() for x in lines]
    properties = dict()
    for i in range(len(property_labels)):
        properties[property_labels[i]] = lines[i+1]
                
#========================================================================================================
def gen_vtx(fileName, pressure=True, field=True, verbose=False):
    '''Generates fn_pressure.vtp and fn_field.vtu in the vtk_files directory, given the filename fn'''

    data = fromNC4(fileName)
    outfn = os.path.split(os.path.splitext(fileName)[0])[1]
    indent3 = '\n\t\t\t'
    indent2 = '\n\t\t'

    if pressure:
        nT = data['nT'][0]-1 #RAM uses a ghost point in local time
        nR = data['nR'][0]
        nS = data['nS'][0]
        rnT = range(nT)
        rnR = range(nR)
        npts = nT*nR
        par_data = data['PParT']; per_data = data['PPerT']

        #-----------------------Generating fn_pressure.vtp----------------------------------
        #Set up XML tree structure
        fulltree = lxml.etree.ElementTree(lxml.etree.Element('VTKFile', type='PolyData', version='1.0', byte_order='LittleEndian'))
        parent = fulltree.getroot()
        RAMdata = lxml.etree.Element('PolyData')
        parent.append(RAMdata)
        piece = lxml.etree.Element('Piece', NumberOfPoints='{}'.format(npts), NumberOfVerts='1', NumberOfLines='0', NumberOfStrips='0', NumberOfPolys='0')
        RAMdata.append(piece)
        pointdata = lxml.etree.Element('PointData', Scalars='proton pressure')
        piece.append(pointdata)
        epress  = lxml.etree.Element('DataArray', type='Float64', Name='electron pressure', NumberOfComponents='1', format='ascii')
        hpress  = lxml.etree.Element('DataArray', type='Float64', Name='proton pressure', NumberOfComponents='1', format='ascii')
        hepress = lxml.etree.Element('DataArray', type='Float64', Name='heliumion pressure', NumberOfComponents='1', format='ascii')
        opress  = lxml.etree.Element('DataArray', type='Float64', Name='oxygenion pressure', NumberOfComponents='1', format='ascii')
        anisot  = lxml.etree.Element('DataArray', type='Float64', Name='pressure anisotropy', NumberOfComponents='1', format='ascii')
        #write data to XML tree
        to_write = '\n'
        for i, j in itertools.product(rnT,rnR):
            to_write +=  '\t\t\t\t\t' + str(par_data[i,j,0] + 2*per_data[i,j,0]) + '\n'
        epress.text = to_write
        to_write = '\n'
        for i, j in itertools.product(rnT,rnR):
            to_write +=  '\t\t\t\t\t' + str(par_data[i,j,1] + 2*per_data[i,j,1]) + '\n'
        hpress.text = to_write
        to_write = '\n'
        for i, j in itertools.product(rnT,rnR):
            to_write +=  '\t\t\t\t\t' + str(par_data[i,j,2] + 2*per_data[i,j,2]) + '\n'
        hepress.text = to_write
        to_write = '\n'
        for i, j in itertools.product(rnT,rnR):
            to_write +=  '\t\t\t\t\t' + str(par_data[i,j,3] + 2*per_data[i,j,3]) + '\n'
        opress.text = to_write
        for el in [epress, hpress, hepress, opress]:
            pointdata.append(el)
        to_write = '\n'
        for i, j in itertools.product(rnT,rnR):
            para = par_data[i,j,0] + par_data[i,j,1] + par_data[i,j,2] + par_data[i,j,3]
            perp = per_data[i,j,0] + per_data[i,j,1] + per_data[i,j,2] + per_data[i,j,3]
            if para==0 and perp==0:
                to_write +=  '\t\t\t\t\t-99.0\n'
            else:
                to_write +=  '\t\t\t\t\t' + str(para/perp) + '\n'
        anisot.text = to_write
        for el in [epress, hpress, hepress, opress, anisot]:
            pointdata.append(el)
        celldata = lxml.etree.Element('CellData')
        piece.append(celldata)

        #Now we add the actual point locations
        points = lxml.etree.Element('Points') #points needs a "DataArray" element with the point locations
        xyzlocations = lxml.etree.Element('DataArray', type='Float64', NumberOfComponents='3', format='ascii')
        #theta_list = np.arange(np.pi, (3*np.pi) + (2*np.pi/(nT-1)), 2*np.pi/(nT-1))
        theta_list = [np.pi+np.deg2rad(15*idx*24.0/nT) for idx in range(nT)]
        step = (6.75 - 1.75)/nR
        r_list = np.arange(1.75 + step, 6.75 + step, step)
        to_write = '\n'
        for i in range(npts):
            r = r_list[i%nR]; theta = theta_list[i//nR]
            to_write += '\t\t\t\t\t' + str(r*np.cos(theta)) + '  ' + str(r*np.sin(theta)) + '  0.00\n'
        xyzlocations.text = to_write
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
                polyconn = getPolyVertOrder(npts, nT, nR)
                nconn = len(polyconn)
                npolys = nT*(nR-1)
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

        out = lxml.etree.tostring(fulltree, xml_declaration=False, pretty_print=True).decode()
        outdir = 'vtk_files'
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        with open(os.path.join(outdir, outfn + '_pressure.vtp'), 'w') as fh:
            fh.write(out)

    if field:
        #-----------------------Generating fn_field.vtu----------------------------------
        nZeta = data['nZeta'][0]-1 #SCB uses a ghost point in local timea
        #nZeta = nZeta//2 #only use half...
        nPsi = data['nPsi'][0]
        nTheta = data['nTheta'][0]
        rnZeta = range(nZeta)
        rnPsi = range(nPsi)
        rnTheta = range(nTheta)

        Bx = data['Bx'].tolist()
        By = data['By'].tolist()
        Bz = data['Bz'].tolist()
        x = data['x'].tolist()
        y = data['y'].tolist()
        z = data['z'].tolist()

        npts = nZeta*nPsi*nTheta
        ncells = (nZeta-1)*(nPsi-1)*(nTheta-1)+((nPsi-1)*(nTheta-1))
        if verbose: print('nPts = {}\nnCells = {}'.format(npts, ncells))

        #Set up XML tree structure
        if verbose: print('Setting up field XML structure')
        fulltree = lxml.etree.ElementTree(lxml.etree.Element('VTKFile', type='UnstructuredGrid', version='1.0', byte_order='LittleEndian'))
        parent = fulltree.getroot()
        SCBdata = lxml.etree.Element('UnstructuredGrid')
        parent.append(SCBdata)
        piece = lxml.etree.Element('Piece', NumberOfPoints='{}'.format(npts), NumberOfCells='{}'.format(ncells))
        SCBdata.append(piece)
        # Point data; any data sets (scalar or vector) at the points need to be given here
        pointdata = lxml.etree.Element('PointData', Vectors='B')
        piece.append(pointdata)
        # Cell data - leave blank as data is all at points
        celldata = lxml.etree.Element('CellData')
        piece.append(celldata)
        #Populate the point data
        Bvec = lxml.etree.Element('DataArray', type='Float64', Name='B', NumberOfComponents='3', format='ascii')
        if verbose: print('Writing B vector to XML structure')
        to_write = '\n'
        for i, j, k in itertools.product(rnZeta, rnPsi, rnTheta):
            to_write += '\t\t\t\t\t{} {} {}\n'.format(Bx[i][j][k], By[i][j][k], Bz[i][j][k])
        Bvec.text = to_write
        pointdata.append(Bvec)
        #Now we add the actual point locations
        points = lxml.etree.Element('Points') #points needs a "DataArray" element with the point locations
        xyzlocations = lxml.etree.Element('DataArray', type='Float64', NumberOfComponents='3', format='ascii')
        to_write = '\n'
        if verbose: print('Writing grid points to XML structure')
        for i, j, k in itertools.product(rnZeta, rnPsi, rnTheta):
            to_write += '\t\t\t\t\t{} {} {} \n'.format(x[i][j][k], y[i][j][k], z[i][j][k])
        xyzlocations.text = to_write
        points.append(xyzlocations)
        piece.append(points)

        #Now define the cells...
        cells = lxml.etree.Element('Cells') # needs a "DataArray" element with the point locations
        cellconn = lxml.etree.Element('DataArray', Name='connectivity', type='Int32')
        celloffs = lxml.etree.Element('DataArray', Name='offsets', type='Int32')
        celltype = lxml.etree.Element('DataArray', Name='types', type='Int32')
        #write the hexahedrons that connect the points together, defining the grid
        if verbose: print('Calculating cell connectivity')
        connarr = getCellVertOrder(npts, nZeta, nPsi, nTheta, wrap=True)
        nconn = len(connarr)
        if verbose: print('connarr has length {}'.format(len(connarr)))
        if verbose: print('nCells (nConn//8) = {}'.format(nconn//8))
        #for i in range(nconn)[::-1]:
        #    connarr.insert((i+1)*8, indent3) ###TODO: This is absurdly slow--- needs fixing. removal just makes the files illegible
        connstr = ' ' + ' '.join(connarr)
        if verbose: print('Writing cell connectivity to XML structure')
        #offsets are 4, 8, 12, ...
        celloffset = ['{}'.format(8*(n+1)) for n in range(nconn//8)]
        #for i in range(nconn//8)[::-1]:
        #    celloffset.insert((i+1)*8, indent3)
        offstr = ' ' + ' '.join(celloffset)
        cellconn.text = indent3 + connstr.rstrip() + indent2
        celloffs.text = indent3 + offstr.rstrip() + indent2
        celltype.text = '\n' + ' '.join(['12']*(nconn//8)) + '\n'

        cells.append(cellconn)
        cells.append(celloffs)
        cells.append(celltype)
        piece.append(cells)
        piece.set('NumberOfCells', '{}'.format(nconn//8))

        out = lxml.etree.tostring(fulltree, xml_declaration=False, pretty_print=True).decode()
        outdir = 'vtk_files'
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        with open(os.path.join(outdir, outfn + '_field.vtu'), 'w') as fh:
            fh.write(out)

#=================================================================================================
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python convertRAMrestart.py <directory_with_NETCDF_files>')
        sys.exit(1)
    read_config() #read configurations from config.txt
    #Make vtk files for all netcdf files in the given directory:
    files = glob.glob(os.path.join(sys.argv[1], '*.nc'))
 
    if properties['Movie'] == 'no':
        for item in files:
            print(item)
            gen_vtx(item)
    else:
        grp = properties['Movie']
        for item in files:
            if item[:len(grp)] == grp:
                gen_vtx(item)

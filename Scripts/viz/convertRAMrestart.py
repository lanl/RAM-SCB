'''Creates the .vts files in vts_files directory for generating RAM-SCB visualizations
   1) Generates fn_field.vts and fn_pressure.vts
        2) Generates the sphere source to get streamlines, if streamline source is sphere
        Usage: python ram_automate1.py <directory_with_NETCDF_files>
        Input: directory containing .nc files'''

import os, sys, itertools, glob
import lxml.etree
import numpy as np
import spacepy.datamodel as dm
from makeCustomSource import gen_sphere, getPolyVertOrder

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
    lines = map(lambda x: x[x.find(':')+1:].strip(), lines)
    properties = dict()
    for i in range(len(property_labels)):
        properties[property_labels[i]] = lines[i+1]
                
#========================================================================================================
def gen_vts(fileName, pressure=True, field=True):
    '''Generates fn_pressure.vtp and fn_field.vts in the vts_files directory, given the filename'''
    #TODO: RAM pressure is now a VTP, SCB field should be converted to VTU

    data = dm.fromHDF5(fileName)
    outfn = os.path.split(os.path.splitext(fileName)[0])[1]
    indent3 = '\n\t\t\t'
    indent2 = '\n\t\t'

    if pressure:
        nT = len(data['nT'])-1 #RAM uses a ghost point in local time
        nR = len(data['nR'])
        nS = len(data['nS'])
        rnT = range(nT)
        rnR = range(nR)
        npts = nT*nR
        par_data = data['PParT']; per_data = data['PPerT']

        #-----------------------Generating fn_pressure.vts----------------------------------
        #Set up XML tree structure
        fulltree = lxml.etree.ElementTree(lxml.etree.Element('VTKFile', type='PolyData', version='1.0', byte_order='LittleEndian'))
        parent = fulltree.getroot()
        RAMdata = lxml.etree.Element('PolyData')
        parent.append(RAMdata)
        piece = lxml.etree.Element('Piece', NumberOfPoints='{}'.format(npts), NumberOfVerts='1', NumberOfLines='0', NumberOfStrips='0', NumberOfPolys='0')
        RAMdata.append(piece)
        pointdata = lxml.etree.Element('PointData', Scalars='electron pressure')
        piece.append(pointdata)
        epress  = lxml.etree.Element('DataArray', type='Float64', Name='electron pressure', NumberOfComponents='1', format='ascii')
        hpress  = lxml.etree.Element('DataArray', type='Float64', Name='proton pressure', NumberOfComponents='1', format='ascii')
        hepress = lxml.etree.Element('DataArray', type='Float64', Name='heliumion pressure', NumberOfComponents='1', format='ascii')
        opress  = lxml.etree.Element('DataArray', type='Float64', Name='oxygenion pressure', NumberOfComponents='1', format='ascii')
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
        celldata = lxml.etree.Element('CellData')
        piece.append(celldata)

        #Now we add the actual point locations
        points = lxml.etree.Element('Points') #points needs a "DataArray" element with the point locations
        xyzlocations = lxml.etree.Element('DataArray', type='Float64', NumberOfComponents='3', format='ascii')
        theta_list = np.arange(0, (2*np.pi) + (2*np.pi/(nT-1)), 2*np.pi/(nT-1))
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

        out = lxml.etree.tostring(fulltree, xml_declaration=False, pretty_print=True)

        with open('vts_files/' + outfn + '_pressure.vtp', 'w') as fh:
            fh.write(out)

    if field:
        #-----------------------Generating fn_field.py----------------------------------
        nZeta = len(data['nZeta'])-1 #SCB uses a ghost point in local time
        nPsi = len(data['nPsi'])
        nTheta = len(data['nTheta'])
        rnZeta = range(nZeta)
        rnPsi = range(nPsi)
        rnTheta = range(nTheta)

        Bx = data['Bx'].tolist()
        By = data['By'].tolist()
        Bz = data['Bz'].tolist()
        x = data['x'].tolist()
        y = data['y'].tolist()
        z = data['z'].tolist()

        to_write = '<?xml version="1.0"?>\n'
        to_write += '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n'
        to_write += '\t<StructuredGrid WholeExtent="0 {} 0 {} 0 {}">\n'.format(nZeta-1, nPsi-1, nTheta-1)
        to_write += '\t\t<Piece Extent="0 {} 0 {} 0 {}">\n'.format(nZeta-1, nPsi-1, nTheta-1)
        to_write += '\t\t\t<PointData Vectors="B">\n'
        to_write += '\t\t\t\t<DataArray type="Float32" Name="B" NumberOfComponents="3" format="ascii">\n'

        for i, j, k in itertools.product(rnZeta, rnPsi, rnTheta):
            to_write += '\t\t\t\t\t{} {} {}\n'.format(Bx[i][j][k], By[i][j][k], Bz[i][j][k]) + '\n'
        to_write += '\t\t\t\t</DataArray>\n'

        to_write += '\t\t\t</PointData>\n'
        to_write += '\t\t\t<CellData>\n'
        to_write += '\t\t\t</CellData>\n'
        to_write += '\t\t\t<Points>\n'
        to_write += '\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'

        for i, j, k in itertools.product(rnZeta, rnPsi, rnTheta):
            to_write += '\t\t\t\t\t{} {} {} \n'.format(x[i][j][k], y[i][j][k], z[i][j][k])
        to_write += '\t\t\t\t</DataArray>\n'

        to_write += '\t\t\t</Points>\n'
        to_write += '\t\t</Piece>\n'
        to_write += '\t</StructuredGrid>\n'
        to_write += '\t</VTKFile>\n'

        with open('vts_files/' + outfn + '_field.vts', 'w') as fh:
            fh.write(to_write)

#=================================================================================================
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python ram_automate1.py <directory_with_NETCDF_files>')
        sys.exit(1)
    read_config() #read configurations from config.txt
    #Get vts files for all netcdf files in the given directory:
    files = glob.glob(os.path.join(sys.argv[1], '*.nc'))
 
    if properties['Movie'] == 'no':
        for item in files:
            gen_vts(item)
    else:
        grp = properties['Movie']
        for item in files:
            if item[:len(grp)] == grp:
                gen_vts(item)

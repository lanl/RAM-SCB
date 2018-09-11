'''Creates the .vts files in vts_files directory for generating RAM-SCB visualizations
   1) Generates fn_field.vts and fn_pressure.vts
        2) Generates the sphere source to get streamlines, if streamline source is sphere
        Usage: python ram_automate1.py <directory_with_NETCDF_files>
        Input: directory containing .nc files'''

import os, sys, itertools
import lxml.etree
import numpy as np
import spacepy.datamodel as dm
from makeCustomSource import gen_sphere

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
def gen_vts(fileName):
    '''Generates fn_pressure.vts and fn_field.vts in the vts_files directory, given the filename'''

    data = dm.fromHDF5(sys.argv[1] + '/' + fileName)
    nT = len(data['nT'])
    nR = len(data['nR'])
    nS = len(data['nS'])
    rnT = range(nT)
    rnR = range(nR)
    par_data = data['PParT']; per_data = data['PPerT']

    #-----------------------Generating fn_pressure.vts----------------------------------
    #Set up XML tree structure
    fulltree = lxml.etree.ElementTree(lxml.etree.Element('VTKFile', type='StructuredGrid', version='1.0', byte_order='LittleEndian'))
    parent = fulltree.getroot()
    RAMdata = lxml.etree.Element('StructuredGrid', WholeExtent='0 {} 0 {} 0 0'.format(nT-1, nR-1))
    parent.append(RAMdata)
    piece = lxml.etree.Element('Piece', Extent='0 {} 0 {} 0 0'.format(nT-1, nR-1))
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
        to_write +=  '\t\t\t\t\t' + str(par_data[i][j][0] + 2*per_data[i][j][0]) + '\n'
    epress.text = to_write
    to_write = '\n'
    for i, j in itertools.product(rnT,rnR):
        to_write +=  '\t\t\t\t\t' + str(par_data[i][j][1] + 2*per_data[i][j][1]) + '\n'
    hpress.text = to_write
    to_write = '\n'
    for i, j in itertools.product(rnT,rnR):
        to_write +=  '\t\t\t\t\t' + str(par_data[i][j][2] + 2*per_data[i][j][2]) + '\n'
    hepress.text = to_write
    to_write = '\n'
    for i, j in itertools.product(rnT,rnR):
        to_write +=  '\t\t\t\t\t' + str(par_data[i][j][3] + 2*per_data[i][j][3]) + '\n'
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
    for i in range(nT*nR):
        r = r_list[i%nR]; theta = theta_list[i//nR]
        to_write += '\t\t\t\t\t' + str(r*np.cos(theta)) + '  ' + str(r*np.sin(theta)) + '  0.00\n'
    xyzlocations.text = to_write
    points.append(xyzlocations)
    piece.append(points)

    out = lxml.etree.tostring(fulltree, xml_declaration=False, pretty_print=True)

    with open('vts_files/' + fileName[:-3] + '_pressure.vts', 'w') as fh:
        fh.write(out)

    #-----------------------Generating fn_field.py----------------------------------
    nZeta = len(data['nZeta'].tolist())
    nPsi = len(data['nPsi'].tolist())
    nTheta = len(data['nTheta'].tolist())

    Bx_list = data['Bx'].tolist()
    By_list = data['By'].tolist()
    Bz_list = data['Bz'].tolist()
    x_list = data['x'].tolist()
    y_list = data['y'].tolist()
    z_list = data['z'].tolist()

    to_write = '<?xml version="1.0"?>\n'
    to_write += '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n'
    to_write += '\t<StructuredGrid WholeExtent="0 {} 0 {} 0 {}">\n'.format(nZeta-1, nPsi-1, nTheta-1)
    to_write += '\t\t<Piece Extent="0 {} 0 {} 0 {}">\n'.format(nZeta-1, nPsi-1, nTheta-1)
    to_write += '\t\t\t<PointData Vectors="B">\n'
    to_write += '\t\t\t\t<DataArray type="Float32" Name="B" NumberOfComponents="3" format="ascii">\n'

    for i in range(nZeta):
        for j in range(nPsi):
            for k in range(nTheta):
                to_write += '\t\t\t\t\t' + str(Bx_list[i][j][k]) + ' ' + str(By_list[i][j][k]) + ' ' + str(Bz_list[i][j][k]) + '\n'
    to_write += '\t\t\t\t</DataArray>\n'

    to_write += '\t\t\t</PointData>\n'
    to_write += '\t\t\t<CellData>\n'
    to_write += '\t\t\t</CellData>\n'
    to_write += '\t\t\t<Points>\n'
    to_write += '\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'

    for i in range(nZeta):
        for j in range(nPsi):
            for k in range(nTheta):
                to_write += '\t\t\t\t\t' + str(x_list[i][j][k]) + ' ' + str(y_list[i][j][k]) + ' ' + str(z_list[i][j][k]) + '\n'
    to_write += '\t\t\t\t</DataArray>\n'

    to_write += '\t\t\t</Points>\n'
    to_write += '\t\t</Piece>\n'
    to_write += '\t</StructuredGrid>\n'
    to_write += '\t</VTKFile>\n'

    with open('vts_files/' + fileName[:-3] + '_field.vts', 'w') as fh:
        fh.write(to_write)

#=================================================================================================
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python ram_automate1.py <directory_with_NETCDF_files>')
        sys.exit(1)
    read_config() #read configurations from config.txt
    #Get vts files for all netcdf files in the given directory:
    files = os.listdir(sys.argv[1])
 
    if properties['Movie'] == 'no':
        for item in files:
            if item[-3:] == '.nc':
                gen_vts(item)
    else:
        grp = properties['Movie']
        for item in files:
            if item[-3:] == '.nc' and item[:len(grp)] == grp:
                gen_vts(item)

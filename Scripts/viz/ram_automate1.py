'''Creates the .vts files in vts_files directory for generating RAM-SCB visualizations
   1) Generates fn_field.vts and fn_pressure.vts
	2) Generates the sphere source to get streamlines, if streamline source is sphere
	Usage: python ram_automate1.py <directory_with_NETCDF_files>
	Input: directory containing .nc files'''

import os, sys
import numpy as np
import spacepy.datamodel as dm
from makeCustomSource import gen_sphere

#========================================================================================================
def read_config():
	'''Reads configurations from config.txt'''
	global properties
	property_labels = ['Streamline source type', 'Plasma pressure display', 'Camera Position',  'Camera Focal Point', 'Camera View Up', 'Scale', 'Movie', 'SaveState']
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
	'''Generates fn_pressure.vts and fn_field.vts in the vts_files directory, given the fileName'''
	data = dm.fromHDF5(sys.argv[1] + '/' + fileName)
	nT = len(data['nT'].tolist())
	nR = len(data['nR'].tolist())
	nS = len(data['nS'].tolist())
	#-----------------------Generating fn_pressure.vts----------------------------------
	to_write = '<?xml version="1.0"?>\n'
	to_write += '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n'
	to_write += '\t<StructuredGrid WholeExtent="0 {} 0 {} 0 0">\n'.format(nT-1, nR-1)
	to_write += '\t\t<Piece Extent="0 {} 0 {} 0 0">\n'.format(nT-1, nR-1)
	to_write += '\t\t\t<PointData Scalars="electron pressure">\n'
	par_data = data['PParT'].tolist(); per_data = data['PPerT'].tolist()

	to_write += '\t\t\t\t<DataArray type="Float32" Name="electron pressure" NumberOfComponents="1" format="ascii">\n'

	for i in range(nT):
		for j in range(nR):
			to_write +=  '\t\t\t\t\t' + str(par_data[i][j][0] + 2*per_data[i][j][0]) + '\n'
	to_write += '\t\t\t\t</DataArray>\n'

	to_write += '\t\t\t\t<DataArray type="Float32" Name="proton pressure" NumberOfComponents="1" format="ascii">\n'

	for i in range(nT):
		for j in range(nR):
			to_write +=  '\t\t\t\t\t' + str(par_data[i][j][1] + 2*per_data[i][j][1]) + '\n'
	to_write += '\t\t\t\t</DataArray>\n'

	to_write += '\t\t\t\t<DataArray type="Float32" Name="heliumion pressure" NumberOfComponents="1" format="ascii">\n'

	for i in range(nT):
		for j in range(nR):
			to_write +=  '\t\t\t\t\t' + str(par_data[i][j][2] + 2*per_data[i][j][2]) + '\n'
	to_write += '\t\t\t\t</DataArray>\n'

	to_write += '\t\t\t\t<DataArray type="Float32" Name="oxygenion pressure" NumberOfComponents="1" format="ascii">\n'

	for i in range(nT):
		for j in range(nR):
			to_write +=  '\t\t\t\t\t' + str(par_data[i][j][3] + 2*per_data[i][j][3]) + '\n'
	to_write += '\t\t\t\t</DataArray>\n'

	to_write += '\t\t\t</PointData>\n'
	to_write += '\t\t\t<CellData>\n'
	to_write += '\t\t\t</CellData>\n'
	to_write += '\t\t\t<Points>\n'
	to_write += '\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'

	theta_list = np.arange(0, (2*np.pi) + (2*np.pi/(nT-1)), 2*np.pi/(nT-1))
	step = (6.75 - 1.75)/nR
	r_list = np.arange(1.75 + step, 6.75 + step, step)

	for i in range(nT*nR):
		r = r_list[i%nR]; theta = theta_list[i//nR]
		to_write += '\t\t\t\t\t' + str(r*np.cos(theta)) + '  ' + str(r*np.sin(theta)) + '  0.00\n'
	to_write += '\t\t\t\t</DataArray>\n'

	to_write += '\t\t\t</Points>\n'
	to_write += '\t\t</Piece>\n'
	to_write += '\t</StructuredGrid>\n'
	to_write += '</VTKFile>\n'

	with open('vts_files/' + fileName[:-3] + '_pressure.vts', 'w') as fh:
		fh.write(to_write)
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
#def gen_sphere():
#	'''Generating the sphere source for streamlines as sphere.vts in vts_files directory'''
#	n = 1 #how much denser points on sphere are, compared to the magnetic field
#	#----------------------------generate the points describing the sphere--------------------------
#	theta = np.linspace(0, 2*np.pi, 873*n).tolist() #latitudes
#	dtor = np.pi/180.0
#	delta = 90.0/(253.0*n)
#
#	#generate 1st hemisphere
#	L_lo = 1.0/np.sin(delta*dtor) #90/(253*n)  {0}
#	L_hi = 1.0/np.sin(90*dtor)  #{90}
#	Lvalues = np.linspace(L_lo,L_hi,253.0*n*n).tolist()
#	phi = map(lambda x: np.arcsin(1.0/x),Lvalues)
#	#generate 2nd hemisphere
#	L_lo = 1.0/np.cos((90 + delta)*dtor)  #{90}
#	L_hi = 1.0/np.cos((-180 + delta)*dtor) #90/(253*n) - 180  {180}
#	Lvalues2 = np.linspace(L_lo,L_hi,252*n*n).tolist()
#	phi2 = map(lambda x: np.arccos(1.0/x),Lvalues2)
#	phi2 = map(lambda x: x+(np.pi/2), phi2)
#
#	#combine the two hemispheres
#	phi.extend(phi2)
#
#	r = 1
#	points = []
#	for t in theta:
#		for p in phi:
#			points.append((r*np.cos(t)*np.sin(p), r*np.sin(t)*np.sin(p), r*np.cos(p)))
#	#----------------------------generate the sphere.vts file--------------------------
#	to_write = ''
#
#	to_write += '<?xml version="1.0"?>\n'  #1
#	to_write += '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n' #2
#	to_write += '\t<StructuredGrid WholeExtent="0 {} 0 {} 0 {}">\n'.format(97*n-1, 45*n-1, 101*n-1) #3
#	to_write += '\t\t<Piece Extent="0 {} 0 {} 0 {}">\n'.format(97*n-1, 45*n-1, 101*n-1) #4
#
#	to_write += '\t\t\t<PointData Scalars="dummy">\n' #5
#	to_write += '\t\t\t\t<DataArray type="Int32" Name="dummy" NumberOfComponents="1" format="ascii">\n' #6
#	for i in range(97*45*101*(n**3)):#440865
#		to_write += '1 ' #7
#	to_write += '\n\t\t\t\t</DataArray>\n' #8
#	to_write += '\t\t\t</PointData>\n' #9
#	to_write += '\t\t\t<CellData>\n'
#	to_write += '\t\t\t</CellData>\n'
#	to_write += '\t\t\t<Points>\n'
#	to_write += '\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'
#	for p in points:
#		to_write += '\t\t\t\t\t' + str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2]) + '\n'
#	to_write += '\t\t\t\t</DataArray>\n'
#	to_write += '\t\t\t</Points>\n'
#
#	to_write += '\t\t</Piece>\n'
#	to_write += '\t</StructuredGrid>\n'
#	to_write += '</VTKFile>\n'
#
#	with open('vts_files/sphere.vts', 'w') as fh:
#		fh.write(to_write)
#
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
		
	if properties['Streamline source type'] == 'sphere':
		gen_sphere() #generates sphere.vts in vts_files directory

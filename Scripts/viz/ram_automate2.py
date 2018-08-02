'''Creates the .png visualizations in the images directory'''
#### import the simple module from the paraview
from paraview.simple import *
import os, sys

#========================================================================================================
def read_config():
	'''Reads configurations from config.txt'''
	global properties
	property_labels = ['Streamline source type', 'Plasma pressure display', 'Camera angle', 'Scale']
	with open('config.txt', 'r') as to_read:
		lines = to_read.readlines()
	lines = map(lambda x: x[x.find(':')+1:].strip(), lines)
	properties = dict()
	for i in range(len(property_labels)):
		properties[property_labels[i]] = lines[i+1]
#========================================================================================================
def gen_viz(fileName):
	#### disable automatic camera reset on 'Show'
	paraview.simple._DisableFirstRenderCameraReset()

	# create a new 'XML Structured Grid Reader'
	field_20130317_T04D_RSCE_GEO_t02000vts = XMLStructuredGridReader(FileName=['vts_files/' + fileName + '_field.vts'])
	field_20130317_T04D_RSCE_GEO_t02000vts.PointArrayStatus = ['B']

	# get active view
	renderView1 = GetActiveViewOrCreate('RenderView')
	renderView1.Background = [1.0, 1.0, 1.0]
	# uncomment following to set a specific view size
	# renderView1.ViewSize = [989, 703]

	# show data in view
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay = Show(field_20130317_T04D_RSCE_GEO_t02000vts, renderView1)
	# trace defaults for the display properties.
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.Representation = 'Outline'
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.AmbientColor = [0.0, 0.0, 0.0]
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.ColorArrayName = [None, '']
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.OSPRayScaleArray = 'B'
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.SelectOrientationVectors = 'B'
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.ScaleFactor = 1.461107873916626
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.SelectScaleArray = 'B'
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.GlyphType = 'Arrow'
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.ScalarOpacityUnitDistance = 0.28164403662270693

	# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
	field_20130317_T04D_RSCE_GEO_t02000vtsDisplay.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# reset view to fit data
	renderView1.ResetCamera()

	#changing interaction mode based on data extents
	renderView1.InteractionMode = '3D'

	#if slice:
	if properties['Streamline source type'] == 'slice':
		# create a new 'Slice'
		slice1 = Slice(Input=field_20130317_T04D_RSCE_GEO_t02000vts)
		slice1.SliceType = 'Plane'
		slice1.SliceOffsetValues = [0.0]

		# init the 'Plane' selected for 'SliceType'
		slice1.SliceType.Origin = [-0.6482846736907959, 0.25454020500183105, -0.00968027114868164]

		# toggle 3D widget visibility (only when running from the GUI)
		Show3DWidgets(proxy=slice1.SliceType)

		# Properties modified on slice1.SliceType
		#slice1.SliceType.Center = [0.0, 0.0, 0.0]
		#slice1.SliceType.Radius = 1.01
		slice1.add_attribute('Center', [0.0, 0.0, 0.0])
		slice1.add_attribute('Radius', 1.01)

		# Properties modified on slice1
		slice1.SliceType = 'Sphere'

		# Properties modified on slice1.SliceType
		slice1.SliceType.Center = [0.0, 0.0, 0.0]
		slice1.SliceType.Radius = 1.01

		# show data in view
		slice1Display = Show(slice1, renderView1)
		# trace defaults for the display properties.
		slice1Display.AmbientColor = [0.0, 0.0, 0.0]
		slice1Display.ColorArrayName = [None, '']
		slice1Display.OSPRayScaleArray = 'B'
		slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
		slice1Display.SelectOrientationVectors = 'B'
		slice1Display.ScaleFactor = 0.18505290746688843
		slice1Display.SelectScaleArray = 'B'
		slice1Display.GlyphType = 'Arrow'
		slice1Display.GaussianRadius = 0.09252645373344422
		slice1Display.SetScaleArray = [None, '']
		slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
		slice1Display.OpacityArray = [None, '']
		slice1Display.OpacityTransferFunction = 'PiecewiseFunction'

		# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
		slice1Display.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
		slice1Display.ScaleTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
		slice1Display.OpacityTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# set active source
		SetActiveSource(field_20130317_T04D_RSCE_GEO_t02000vts)

		# hide data in view
		Hide(slice1, renderView1)

		# create a new 'Stream Tracer With Custom Source'
		streamTracerWithCustomSource1 = StreamTracerWithCustomSource(Input=field_20130317_T04D_RSCE_GEO_t02000vts,
			 SeedSource=slice1)
		streamTracerWithCustomSource1.Vectors = ['POINTS', 'B']
		streamTracerWithCustomSource1.MaximumStreamlineLength = 14.61107873916626

		# Properties modified on streamTracerWithCustomSource1
		streamTracerWithCustomSource1.IntegratorType = 'Runge-Kutta 4'
		streamTracerWithCustomSource1.SurfaceStreamlines = 1


		# show data in view
		streamTracerWithCustomSource1Display = Show(streamTracerWithCustomSource1, renderView1)
		# trace defaults for the display properties.
		streamTracerWithCustomSource1Display.AmbientColor = [0.0, 0.0, 0.0]
		streamTracerWithCustomSource1Display.ColorArrayName = [None, '']
		streamTracerWithCustomSource1Display.OSPRayScaleArray = 'AngularVelocity'
		streamTracerWithCustomSource1Display.OSPRayScaleFunction = 'PiecewiseFunction'
		streamTracerWithCustomSource1Display.SelectOrientationVectors = 'Normals'
		streamTracerWithCustomSource1Display.ScaleFactor = 0.9971208572387695
		streamTracerWithCustomSource1Display.SelectScaleArray = 'AngularVelocity'
		streamTracerWithCustomSource1Display.GlyphType = 'Arrow'
		streamTracerWithCustomSource1Display.GaussianRadius = 0.49856042861938477
		streamTracerWithCustomSource1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
		streamTracerWithCustomSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
		streamTracerWithCustomSource1Display.OpacityArray = ['POINTS', 'AngularVelocity']
		streamTracerWithCustomSource1Display.OpacityTransferFunction = 'PiecewiseFunction'

		# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
		streamTracerWithCustomSource1Display.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
		streamTracerWithCustomSource1Display.ScaleTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
		streamTracerWithCustomSource1Display.OpacityTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# hide data in view
		Hide(slice1, renderView1)

	elif properties['Streamline source type'] == 'sphere':
		# create a new 'XML Structured Grid Reader'
		spherevts = XMLStructuredGridReader(FileName=['vts_files/sphere.vts'])
		spherevts.PointArrayStatus = ['dummy']

		# show data in view
		spherevtsDisplay = Show(spherevts, renderView1)
		# trace defaults for the display properties.
		spherevtsDisplay.Representation = 'Outline'
		spherevtsDisplay.ColorArrayName = ['POINTS', '']
		spherevtsDisplay.OSPRayScaleArray = 'dummy'
		spherevtsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
		spherevtsDisplay.SelectOrientationVectors = 'None'
		spherevtsDisplay.ScaleFactor = 0.20199999809265137
		spherevtsDisplay.SelectScaleArray = 'dummy'
		spherevtsDisplay.GlyphType = 'Arrow'
		spherevtsDisplay.ScalarOpacityUnitDistance = 0.04663026724290165

		# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
		spherevtsDisplay.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# set active source
		SetActiveSource(field_20130317_T04D_RSCE_GEO_t02000vts)

		# create a new 'Stream Tracer With Custom Source'
		streamTracerWithCustomSource1 = StreamTracerWithCustomSource(Input=field_20130317_T04D_RSCE_GEO_t02000vts, SeedSource=spherevts)
		streamTracerWithCustomSource1.Vectors = ['POINTS', 'B']
		streamTracerWithCustomSource1.MaximumStreamlineLength = 14.61107873916626

		# Properties modified on streamTracerWithCustomSource1
		#streamTracerWithCustomSource1.SurfaceStreamlines = 1
		streamTracerWithCustomSource1.IntegratorType = 'Runge-Kutta 4'

		# show data in view
		streamTracerWithCustomSource1Display = Show(streamTracerWithCustomSource1, renderView1)
		# trace defaults for the display properties.
		streamTracerWithCustomSource1Display.ColorArrayName = [None, '']
		streamTracerWithCustomSource1Display.OSPRayScaleArray = 'AngularVelocity'
		streamTracerWithCustomSource1Display.OSPRayScaleFunction = 'PiecewiseFunction'
		streamTracerWithCustomSource1Display.SelectOrientationVectors = 'Normals'
		streamTracerWithCustomSource1Display.ScaleFactor = 0.8783905982971192
		streamTracerWithCustomSource1Display.SelectScaleArray = 'AngularVelocity'
		streamTracerWithCustomSource1Display.GlyphType = 'Arrow'
		streamTracerWithCustomSource1Display.GaussianRadius = 0.4391952991485596
		streamTracerWithCustomSource1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
		streamTracerWithCustomSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
		streamTracerWithCustomSource1Display.OpacityArray = ['POINTS', 'AngularVelocity']
		streamTracerWithCustomSource1Display.OpacityTransferFunction = 'PiecewiseFunction'

		# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
		streamTracerWithCustomSource1Display.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
		streamTracerWithCustomSource1Display.ScaleTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
		streamTracerWithCustomSource1Display.OpacityTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

		# hide data in view
		Hide(field_20130317_T04D_RSCE_GEO_t02000vts, renderView1)

		# hide data in view
		Hide(spherevts, renderView1)
	else:
		raise ValueError('Streamline source should be either sphere or slice')
		sys.exit(1)
	#------------------------------------------------------------------------------------------------------

	# set scalar coloring
	ColorBy(streamTracerWithCustomSource1Display, ('POINTS', 'B'))

	# rescale color and/or opacity maps used to include current data range
	streamTracerWithCustomSource1Display.RescaleTransferFunctionToDataRange(True, False)

	# show color bar/color legend
	streamTracerWithCustomSource1Display.SetScalarBarVisibility(renderView1, True)

	# get color transfer function/color map for 'B'
	bLUT = GetColorTransferFunction('B')

	# get opacity transfer function/opacity map for 'B'
	bPWF = GetOpacityTransferFunction('B')

	# hide data in view
	Hide(field_20130317_T04D_RSCE_GEO_t02000vts, renderView1)

	# create a new 'Tube'
	tube1 = Tube(Input=streamTracerWithCustomSource1)
	tube1.Scalars = ['POINTS', 'AngularVelocity']
	tube1.Vectors = ['POINTS', 'Normals']
	tube1.Radius = 0.09971208572387695

	# Properties modified on tube1
	tube1.Radius = 0.016951054573059083

	# show data in view
	tube1Display = Show(tube1, renderView1)
	# trace defaults for the display properties.
	tube1Display.AmbientColor = [0.0, 0.0, 0.0]
	tube1Display.ColorArrayName = ['POINTS', 'B']
	tube1Display.LookupTable = bLUT
	tube1Display.OSPRayScaleArray = 'AngularVelocity'
	tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
	tube1Display.SelectOrientationVectors = 'Normals'
	tube1Display.ScaleFactor = 1.0000883102416993
	tube1Display.SelectScaleArray = 'AngularVelocity'
	tube1Display.GlyphType = 'Arrow'
	tube1Display.GaussianRadius = 0.5000441551208497
	tube1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
	tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
	tube1Display.OpacityArray = ['POINTS', 'AngularVelocity']
	tube1Display.OpacityTransferFunction = 'PiecewiseFunction'

	# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
	tube1Display.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
	tube1Display.ScaleTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
	tube1Display.OpacityTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# hide data in view
	Hide(streamTracerWithCustomSource1, renderView1)

	# show color bar/color legend
	tube1Display.SetScalarBarVisibility(renderView1, True)

	# create a new 'XML Structured Grid Reader'
	pressure_20130317_T04D_RSCE_GEO_t02000vts = XMLStructuredGridReader(FileName=['vts_files/' + fileName + '_pressure.vts'])
	pressure_20130317_T04D_RSCE_GEO_t02000vts.PointArrayStatus = ['electron pressure', 'proton pressure', 'helium ion pressure', 'oxygen ion pressure']

	# get color transfer function/color map for 'electronpressure'
	#pressureLUT = GetColorTransferFunction('electronpressure')
	pressureLUT = GetColorTransferFunction(properties['Plasma pressure display'] + 'pressure')
	pressureLUT.RescaleOnVisibilityChange = 1

	# show data in view
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay = Show(pressure_20130317_T04D_RSCE_GEO_t02000vts, renderView1)
	# trace defaults for the display properties.
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.AmbientColor = [0.0, 0.0, 0.0]
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.ColorArrayName = ['POINTS', properties['Plasma pressure display'] + 'pressure']
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.LookupTable = pressureLUT
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.OSPRayScaleArray = properties['Plasma pressure display'] + 'pressure'
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.SelectOrientationVectors = 'None'
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.ScaleFactor = 1.35
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.SelectScaleArray = properties['Plasma pressure display'] + 'pressure'
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.GlyphType = 'Arrow'
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.ScalarOpacityUnitDistance = 2.480431009474945

	# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# show color bar/color legend
	pressure_20130317_T04D_RSCE_GEO_t02000vtsDisplay.SetScalarBarVisibility(renderView1, True)

	# Rescale transfer function
	bLUT.RescaleTransferFunction(1.97797401027, 591.487212231)
	pressureLUT.MapControlPointsToLogSpace()

	# Rescale transfer function
	bPWF.RescaleTransferFunction(1.97797401027, 591.487212231)

	# get opacity transfer function/opacity map for 'electronpressure'
	pressurePWF = GetOpacityTransferFunction(properties['Plasma pressure display'] + 'pressure')

	# create a new 'Point Volume Interpolator'
	pointVolumeInterpolator1 = PointVolumeInterpolator(Input=pressure_20130317_T04D_RSCE_GEO_t02000vts,
		 Source='Bounded Volume')
	pointVolumeInterpolator1.Kernel = 'VoronoiKernel'
	pointVolumeInterpolator1.Locator = 'Static Point Locator'

	# init the 'Bounded Volume' selected for 'Source'
	pointVolumeInterpolator1.Source.Origin = [-6.75, -6.75, 0.0]
	pointVolumeInterpolator1.Source.Scale = [13.5, 13.5, 0.0]

	# show data in view
	pointVolumeInterpolator1Display = Show(pointVolumeInterpolator1, renderView1)
	# trace defaults for the display properties.
	pointVolumeInterpolator1Display.Representation = 'Outline'
	pointVolumeInterpolator1Display.AmbientColor = [0.0, 0.0, 0.0]
	pointVolumeInterpolator1Display.ColorArrayName = ['POINTS', properties['Plasma pressure display'] + 'pressure']
	pointVolumeInterpolator1Display.LookupTable = pressureLUT
	pointVolumeInterpolator1Display.OSPRayScaleArray = properties['Plasma pressure display'] + 'pressure'
	pointVolumeInterpolator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
	pointVolumeInterpolator1Display.SelectOrientationVectors = 'None'
	pointVolumeInterpolator1Display.ScaleFactor = 1.35
	pointVolumeInterpolator1Display.SelectScaleArray = properties['Plasma pressure display'] + 'pressure'
	pointVolumeInterpolator1Display.GlyphType = 'Arrow'
	pointVolumeInterpolator1Display.ScalarOpacityUnitDistance = 0.1909188309203679
	pointVolumeInterpolator1Display.Slice = 50

	# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
	pointVolumeInterpolator1Display.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# hide data in view
	Hide(pressure_20130317_T04D_RSCE_GEO_t02000vts, renderView1)

	# show color bar/color legend
	pointVolumeInterpolator1Display.SetScalarBarVisibility(renderView1, True)

	# create a new 'Clip'
	clip1 = Clip(Input=pointVolumeInterpolator1)
	clip1.ClipType = 'Plane'
	clip1.Scalars = ['POINTS', properties['Plasma pressure display'] + 'pressure']
	clip1.Value = 5.546324253082275

	# Rescale transfer function
	bLUT.RescaleTransferFunction(1.97797401027, 591.487212231)

	# toggle 3D widget visibility (only when running from the GUI)
	Show3DWidgets(proxy=clip1.ClipType)

	# Properties modified on clip1
	clip1.ClipType = 'Sphere'
	clip1.ClipType.Radius = 6.75

	# show data in view
	clip1Display = Show(clip1, renderView1)
	# trace defaults for the display properties.
	clip1Display.AmbientColor = [0.0, 0.0, 0.0]
	clip1Display.ColorArrayName = ['POINTS', properties['Plasma pressure display'] + 'pressure']
	clip1Display.LookupTable = pressureLUT
	clip1Display.OSPRayScaleArray = properties['Plasma pressure display'] + 'pressure'
	clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
	clip1Display.SelectOrientationVectors = 'None'
	clip1Display.ScaleFactor = 1.35
	clip1Display.SelectScaleArray = properties['Plasma pressure display'] + 'pressure'
	clip1Display.GlyphType = 'Arrow'
	clip1Display.ScalarOpacityUnitDistance = 0.3035512141561865
	clip1Display.GaussianRadius = 0.675
	clip1Display.SetScaleArray = ['POINTS', properties['Plasma pressure display'] + 'pressure']
	clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
	clip1Display.OpacityArray = ['POINTS', properties['Plasma pressure display'] + 'pressure']
	clip1Display.OpacityTransferFunction = 'PiecewiseFunction'

	if properties['Plasma pressure display'][-3:] == 'ion':
		species = properties['Plasma pressure display'][:-4] + ' ion pressure'  
	else:
		species = properties['Plasma pressure display'] + ' pressure'
	ColorBy(clip1Display, ('POINTS', species))

	# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
	clip1Display.OSPRayScaleFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
	clip1Display.ScaleTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
	clip1Display.OpacityTransferFunction.Points = [0.00253301385078493, 0.0, 0.5, 0.0, 553.421725972081, 1.0, 0.5, 0.0]

	# show color bar/color legend
	clip1Display.SetScalarBarVisibility(renderView1, True)

	# hide data in view
	Hide(pointVolumeInterpolator1, renderView1)

	# Rescale transfer function
	bLUT.RescaleTransferFunction(1.97797401027, 591.487212231)

	# convert to log space
	bLUT.MapControlPointsToLogSpace()

	# Properties modified on bLUT
	bLUT.UseLogScale = 1
	bLUT.ApplyPreset('Blues', True)

	# Rescale transfer function
	pressureLUT.RescaleTransferFunction(0.178525596857, 1.18874013424)

	# Rescale transfer function
	pressurePWF.RescaleTransferFunction(0.178525596857, 1.18874013424)

	# Properties modified on clip1
	clip1.InsideOut = 1

	# Rescale transfer function
	pressureLUT.RescaleTransferFunction(0.0, 11.0926485062)

	# Rescale transfer function
	pressurePWF.RescaleTransferFunction(0.0, 11.0926485062)

	# toggle 3D widget visibility (only when running from the GUI)
	Hide3DWidgets(proxy=clip1.ClipType)

	# convert to log space
	pressureLUT.MapControlPointsToLogSpace()

	# Properties modified on pressureLUT
	pressureLUT.UseLogScale = 1
	pressureLUT.ApplyPreset('Inferno (matplotlib)', True)

	#### saving camera placements for all active views

	# current camera placement for renderView1
	renderView1.CameraPosition = [-15.328799375388362, -23.58278282367084, 28.387079424203996]
	renderView1.CameraFocalPoint = [-0.43224358558654785, 0.200559139251709, -0.006810903549194211]
	renderView1.CameraViewUp = [-0.019392088078452663, 0.7714426265022851, 0.6360033183366368]
	renderView1.CameraParallelScale = 5.332579277749032

	#Adjusting color legend properties
	bLUTColorBar = GetScalarBar(bLUT, renderView1)
	bLUTColorBar.TitleFontSize = 7
	bLUTColorBar.LabelFontSize = 7
	bLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
	bLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
	#bLUTColorBar.AspectRatio = 25
	bLUTColorBar.add_attribute('AspectRatio', 130)
	pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)
	pressureLUTColorBar.TitleFontSize = 7
	pressureLUTColorBar.LabelFontSize = 7
	pressureLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
	pressureLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
	#pressureLUTColorBar.AspectRatio = 25
	pressureLUTColorBar.add_attribute('AspectRatio', 130)

	bLUT.add_attribute('Position2', [200,900])

	renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
	renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]

	if properties['Scale'] == 'on':
		renderView1.AxesGrid.Visibility = 1
		renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
		renderView1.AxesGrid.XTitleFontSize = 9
		renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
		renderView1.AxesGrid.YTitleFontSize = 9
		renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
		renderView1.AxesGrid.ZTitleFontSize = 9
		renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
		renderView1.AxesGrid.XLabelFontSize = 9
		renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
		renderView1.AxesGrid.YLabelFontSize = 9
		renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
		renderView1.AxesGrid.ZLabelFontSize = 9

	renderView1.ViewSize = [600, 600]
	#renderView1.add_attribute('ViewSize', [400, 400])
	SaveScreenshot('images/' +  fileName + '_viz.png', magnification=1.75, quality=600, view=renderView1)
#========================================================================================================
if __name__ == '__main__':
	read_config()	

	files = os.listdir('vts_files')
	for item in files:
		if item[-10:] == '_field.vts':
			gen_viz(item[:-10])

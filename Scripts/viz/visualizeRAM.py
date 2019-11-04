'''Creates the .png (image)/ .avi (video) visualizations in the images directory'''

#### import the simple module from the paraview
from paraview.simple import *
import os, sys

#========================================================================================================
def read_config():
    '''Reads configurations from config.txt'''
    global properties, state_written
    state_written = False
    property_labels = ['Streamline source type', 'Plasma pressure display', 'Camera Position', 'Camera Focal Point', 'Camera View Up', 'Scale', 'Movie', 'SaveState']

    with open('config.txt', 'r') as to_read:
            lines = to_read.readlines()
    lines = [x[x.find(':')+1:].strip() for x in lines]
    properties = dict()
    for i in range(len(property_labels)):
                properties[property_labels[i]] = lines[i+1]

#========================================================================================================
def generateVisualization(pressurefile, fieldfile, pointsfile, opacity=True):
    '''Reads the VTK files, applies filters and saves the visualizations'''
    global state_written

    paraview.simple._DisableFirstRenderCameraReset()

    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    
    # Create a new view
    renderView1 = CreateView('RenderView')
    renderView1.ViewSize = [1100, 600]
    renderView1.AxesGrid = 'GridAxes3DActor'
    renderView1.OrientationAxesLabelColor = [0.3254901960784314, 0.3254901960784314, 0.3254901960784314]
    renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.0]
    renderView1.HiddenLineRemoval = 1
    renderView1.StereoType = 0
    # 
    renderView1.CameraPosition = [-13.526998574661805, -27.266201179785636, 5.937741253914423]
    renderView1.CameraFocalPoint = [0.8335685946434308, 0.36600320041523715, 0.5361380602473481]
    renderView1.CameraViewUp = [0.08620298865162888, 0.14779944232659809, 0.98525345449558]
    renderView1.CameraParallelScale = 8.56237586950119
    renderView1.Background = [1.0, 1.0, 1.0]
    renderView1.OSPRayMaterialLibrary = materialLibrary1
    
    # init the 'GridAxes3DActor' selected for 'AxesGrid'
    renderView1.AxesGrid.Visibility = 1
    renderView1.AxesGrid.XTitle = 'X SM'
    renderView1.AxesGrid.YTitle = 'Y SM'
    renderView1.AxesGrid.ZTitle = 'Z SM'
    renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.XTitleFontFile = ''
    renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.YTitleFontFile = ''
    renderView1.AxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ZTitleFontFile = ''
    renderView1.AxesGrid.GridColor = [0.5451, 0.5451, 0.5451]
    renderView1.AxesGrid.ShowGrid = 1
    renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.XLabelFontFile = ''
    renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.YLabelFontFile = ''
    renderView1.AxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
    renderView1.AxesGrid.ZLabelFontFile = ''
    renderView1.AxesGrid.XAxisNotation = 'Fixed'
    renderView1.AxesGrid.YAxisNotation = 'Fixed'
    renderView1.AxesGrid.ZAxisNotation = 'Fixed'
    renderView1.AxesGrid.DataScale = [6.5, 1.0, 1.0]
    renderView1.AxesGrid.DataBoundsScaleFactor = 1.0
    
    # ----------------------------------------------------------------
    # restore active view
    SetActiveView(renderView1)
    # ----------------------------------------------------------------
    
    # ----------------------------------------------------------------
    # setup the data processing pipelines
    # ----------------------------------------------------------------
    
    # create a new 'XML PolyData Reader'
    # RAM domain data
    pressurevtp = XMLPolyDataReader(FileName=[pressurefile])
    pressurevtp.PointArrayStatus = ['electron pressure', 'proton pressure', 'heliumion pressure', 'oxygenion pressure', 'pressure anisotropy']
    
    # create a new 'XML PolyData Reader'
    # Seed points for streamline tracing
    discvtp = XMLPolyDataReader(FileName=[pointsfile])
    
    # create a new 'XML Unstructured Grid Reader'
    # SCB domain data
    fieldvtu = XMLUnstructuredGridReader(FileName=[fieldfile])
    fieldvtu.PointArrayStatus = ['B']
    
    # create a new 'Stream Tracer With Custom Source'
    streamTracerWithCustomSource1 = StreamTracerWithCustomSource(Input=fieldvtu,
        SeedSource=discvtp)
    streamTracerWithCustomSource1.Vectors = ['POINTS', 'B']
    streamTracerWithCustomSource1.IntegrationStepUnit = 'Length'
    streamTracerWithCustomSource1.InitialStepLength = 0.025
    streamTracerWithCustomSource1.MinimumStepLength = 0.005
    streamTracerWithCustomSource1.MaximumStepLength = 0.075
    streamTracerWithCustomSource1.MaximumSteps = 25000
    streamTracerWithCustomSource1.MaximumStreamlineLength = 40.0
    streamTracerWithCustomSource1.MaximumError = 1e-06
    streamTracerWithCustomSource1.ComputeVorticity = 0
    
    # create a new 'Calculator'
    calculator1 = Calculator(Input=pressurevtp)
    calculator1.ResultArrayName = 'Total pressure'
    calculator1.Function = '(electron pressure+heliumion pressure+oxygenion pressure+proton pressure)*0.16'
    
    # ----------------------------------------------------------------
    # set up the visualization in view 'renderView1'
    # ----------------------------------------------------------------
    
    # show data from streamTracerWithCustomSource1
    streamTracerWithCustomSource1Display = Show(streamTracerWithCustomSource1, renderView1)
    
    # get color transfer function/color map for 'B'
    Bmin = 5.0#nT
    Bmax = 4500.0 #nT
    bLUT = GetColorTransferFunction('B')
    bLUT.RescaleTransferFunction(Bmin, Bmax)
    bLUT.MapControlPointsToLogSpace()
    bLUT.UseLogScale = 1
    bLUT.ApplyPreset('Linear YGB 1211g', True)
    bLUT.ColorSpace = 'Lab'
    #bLUT.NanColor = [0.25, 0.0, 0.0]
    bLUT.ScalarRangeInitialized = 1.0
    
    # trace defaults for the display properties.
    streamTracerWithCustomSource1Display.Representation = 'Surface'
    streamTracerWithCustomSource1Display.ColorArrayName = ['POINTS', 'B']
    streamTracerWithCustomSource1Display.LookupTable = bLUT
    streamTracerWithCustomSource1Display.LineWidth = 2.5
    streamTracerWithCustomSource1Display.RenderLinesAsTubes = 1
    streamTracerWithCustomSource1Display.Specular = 0.41
    streamTracerWithCustomSource1Display.SpecularPower = 71.0
    streamTracerWithCustomSource1Display.Luminosity = 100.0
    streamTracerWithCustomSource1Display.Ambient = 0.3
    streamTracerWithCustomSource1Display.Diffuse = 0.74
    streamTracerWithCustomSource1Display.OSPRayScaleArray = 'B'
    streamTracerWithCustomSource1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    streamTracerWithCustomSource1Display.SelectOrientationVectors = 'B'
    streamTracerWithCustomSource1Display.ScaleFactor = 1.3390784740448
    streamTracerWithCustomSource1Display.SelectScaleArray = 'B'
    streamTracerWithCustomSource1Display.GlyphType = 'Arrow'
    streamTracerWithCustomSource1Display.GlyphTableIndexArray = 'B'
    streamTracerWithCustomSource1Display.GaussianRadius = 0.06695392370224
    streamTracerWithCustomSource1Display.SetScaleArray = ['POINTS', 'B']
    streamTracerWithCustomSource1Display.ScaleTransferFunction = 'PiecewiseFunction'
    streamTracerWithCustomSource1Display.OpacityArray = ['POINTS', 'B']
    streamTracerWithCustomSource1Display.OpacityTransferFunction = 'PiecewiseFunction'
    streamTracerWithCustomSource1Display.DataAxesGrid = 'GridAxesRepresentation'
    streamTracerWithCustomSource1Display.SelectionCellLabelFontFile = ''
    streamTracerWithCustomSource1Display.SelectionPointLabelFontFile = ''
    streamTracerWithCustomSource1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    streamTracerWithCustomSource1Display.ScaleTransferFunction.Points = [-265.9908051551271, 0.0, 0.5, 0.0, 294.75470870808783, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    streamTracerWithCustomSource1Display.OpacityTransferFunction.Points = [-265.9908051551271, 0.0, 0.5, 0.0, 294.75470870808783, 1.0, 0.5, 0.0]
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    streamTracerWithCustomSource1Display.DataAxesGrid.XTitle = 'X'
    streamTracerWithCustomSource1Display.DataAxesGrid.YTitle = 'Y'
    streamTracerWithCustomSource1Display.DataAxesGrid.ZTitle = 'Z'
    streamTracerWithCustomSource1Display.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
    streamTracerWithCustomSource1Display.DataAxesGrid.XTitleFontFamily = 'Times'
    streamTracerWithCustomSource1Display.DataAxesGrid.XTitleFontFile = ''
    streamTracerWithCustomSource1Display.DataAxesGrid.XTitleFontSize = 14
    streamTracerWithCustomSource1Display.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
    streamTracerWithCustomSource1Display.DataAxesGrid.YTitleFontFile = ''
    streamTracerWithCustomSource1Display.DataAxesGrid.YTitleFontSize = 14
    streamTracerWithCustomSource1Display.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
    streamTracerWithCustomSource1Display.DataAxesGrid.ZTitleFontFile = ''
    streamTracerWithCustomSource1Display.DataAxesGrid.ZTitleFontSize = 14
    streamTracerWithCustomSource1Display.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
    streamTracerWithCustomSource1Display.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
    streamTracerWithCustomSource1Display.DataAxesGrid.XLabelFontFile = ''
    streamTracerWithCustomSource1Display.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
    streamTracerWithCustomSource1Display.DataAxesGrid.YLabelFontFile = ''
    streamTracerWithCustomSource1Display.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
    streamTracerWithCustomSource1Display.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    streamTracerWithCustomSource1Display.PolarAxes.PolarAxisTitleFontFile = ''
    streamTracerWithCustomSource1Display.PolarAxes.PolarAxisLabelFontFile = ''
    streamTracerWithCustomSource1Display.PolarAxes.LastRadialAxisTextFontFile = ''
    streamTracerWithCustomSource1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

    # show data from calculator1
    calculator1Display = Show(calculator1, renderView1)
    
    # get color transfer function/color map for 'Totalpressure'
    totalpressureLUT = GetColorTransferFunction('Totalpressure')
    totalpressureLUT.ApplyPreset('Inferno (matplotlib)', True)
    totalpressureLUT.NanColor = [0.0, 1.0, 0.0]
    totalpressureLUT.ScalarRangeInitialized = 1.0
    
    # trace defaults for the display properties.
    calculator1Display.Representation = 'Surface'
    calculator1Display.ColorArrayName = ['POINTS', 'Total pressure']
    calculator1Display.LookupTable = totalpressureLUT
    calculator1Display.Ambient = 0.22
    calculator1Display.OSPRayScaleArray = 'Total pressure'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'None'
    calculator1Display.ScaleFactor = 1.3492461110358438
    calculator1Display.SelectScaleArray = 'Total pressure'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.GlyphTableIndexArray = 'Total pressure'
    calculator1Display.GaussianRadius = 0.0674623055517922
    calculator1Display.SetScaleArray = ['POINTS', 'Total pressure']
    calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    calculator1Display.OpacityArray = ['POINTS', 'Total pressure']
    calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
    calculator1Display.SelectionCellLabelFontFile = ''
    calculator1Display.SelectionPointLabelFontFile = ''
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    pressureMax = 35.0
    calculator1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, pressureMax, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    calculator1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, pressureMax, 1.0, 0.5, 0.0]
    if opacity: calculator1Display.Opacity = 0.8
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    calculator1Display.DataAxesGrid.XTitleFontFile = ''
    calculator1Display.DataAxesGrid.YTitleFontFile = ''
    calculator1Display.DataAxesGrid.ZTitleFontFile = ''
    calculator1Display.DataAxesGrid.XLabelFontFile = ''
    calculator1Display.DataAxesGrid.YLabelFontFile = ''
    calculator1Display.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    calculator1Display.PolarAxes.PolarAxisTitleFontFile = ''
    calculator1Display.PolarAxes.PolarAxisLabelFontFile = ''
    calculator1Display.PolarAxes.LastRadialAxisTextFontFile = ''
    calculator1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    
    # setup the color legend parameters for each legend in this view
    
    # get color legend/bar for bLUT in view renderView1
    bLUTColorBar = GetScalarBar(bLUT, renderView1)
    bLUTColorBar.AutoOrient = 0
    bLUTColorBar.Orientation = 'Horizontal'
    bLUTColorBar.WindowLocation = 'AnyLocation'
    bLUTColorBar.Position = [0.6, 0.2]
    bLUTColorBar.Title = '|B|'
    bLUTColorBar.ComponentTitle = '[nT]'
    bLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
    bLUTColorBar.TitleFontFile = ''
    bLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
    bLUTColorBar.LabelFontFile = ''
    bLUTColorBar.ScalarBarThickness = 20
    
    # set color bar visibility
    bLUTColorBar.Visibility = 1
    
    # get color legend/bar for totalpressureLUT in view renderView1
    totalpressureLUTColorBar = GetScalarBar(totalpressureLUT, renderView1)
    totalpressureLUTColorBar.AutoOrient = 0
    totalpressureLUTColorBar.Orientation = 'Horizontal'
    totalpressureLUTColorBar.WindowLocation = 'AnyLocation'
    totalpressureLUTColorBar.Position = [0.6, 0.725]
    totalpressureLUTColorBar.Title = 'Total pressure'
    totalpressureLUTColorBar.ComponentTitle = '(nPa)'
    totalpressureLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
    totalpressureLUTColorBar.TitleFontFile = ''
    totalpressureLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
    totalpressureLUTColorBar.LabelFontFile = ''
    totalpressureLUTColorBar.ScalarBarThickness = 20
    
    # set color bar visibility
    totalpressureLUTColorBar.Visibility = 1
    
    # ----------------------------------------------------------------
    # setup color maps and opacity maps used in the visualization
    # note: the Get..() functions create a new object, if needed
    # ----------------------------------------------------------------
    
    # get opacity transfer function/opacity map for 'B'
    bPWF = GetOpacityTransferFunction('B')
    bPWF.Points = [0.5143360840384497, 0.0, 0.5, 0.0, 4800.0, 1.0, 0.5, 0.0]
    bPWF.ScalarRangeInitialized = 1
    
    # Rescale transfer function
    bPWF.RescaleTransferFunction(Bmin, Bmax)
    
    # get opacity transfer function/opacity map for 'Totalpressure'
    totalpressurePWF = GetOpacityTransferFunction('Totalpressure')
    totalpressurePWF.Points = [0.0, 0.0, 0.5, 0.0, pressureMax, 1.0, 0.5, 0.0]
    totalpressurePWF.ScalarRangeInitialized = 1
    
    # ----------------------------------------------------------------
    # finally, restore active source
    #SetActiveSource(calculator1)
    # ----------------------------------------------------------------

    # Add a clip filter so we're not looking at the full domain
    # create a new 'Clip'
    clip1 = Clip(Input=streamTracerWithCustomSource1)
    clip1.ClipType = 'Plane'
    clip1.Scalars = ['POINTS', 'AngularVelocity']
    clip1.Value = -18860188.130234756
    
    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [0.0, 0.0, 0.0]
    
    # Properties modified on clip1.ClipType
    clip1.ClipType.Normal = [0.0, -1.0, 0.0]
    
    # show data in view
    clip1Display = Show(clip1, renderView1)
    
    # trace defaults for the display properties.
    clip1Display.Representation = 'Surface'
    clip1Display.ColorArrayName = ['POINTS', 'B']
    clip1Display.LookupTable = bLUT
    clip1Display.OSPRayScaleArray = 'AngularVelocity'
    clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    clip1Display.SelectOrientationVectors = 'Normals'
    clip1Display.ScaleFactor = 1.0107592344284058
    clip1Display.SelectScaleArray = 'AngularVelocity'
    clip1Display.GlyphType = 'Arrow'
    clip1Display.GlyphTableIndexArray = 'AngularVelocity'
    clip1Display.GaussianRadius = 0.05053796172142029
    clip1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
    clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
    clip1Display.OpacityArray = ['POINTS', 'AngularVelocity']
    clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
    clip1Display.DataAxesGrid = 'GridAxesRepresentation'
    clip1Display.SelectionCellLabelFontFile = ''
    clip1Display.SelectionPointLabelFontFile = ''
    clip1Display.PolarAxes = 'PolarAxesRepresentation'
    clip1Display.ScalarOpacityFunction = bPWF
    clip1Display.ScalarOpacityUnitDistance = 0.4975569674139967
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    clip1Display.ScaleTransferFunction.Points = [-8426125.789073758, 0.0, 0.5, 0.0, 1476579.976734953, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    clip1Display.OpacityTransferFunction.Points = [-8426125.789073758, 0.0, 0.5, 0.0, 1476579.976734953, 1.0, 0.5, 0.0]
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    clip1Display.DataAxesGrid.XTitleFontFile = ''
    clip1Display.DataAxesGrid.YTitleFontFile = ''
    clip1Display.DataAxesGrid.ZTitleFontFile = ''
    clip1Display.DataAxesGrid.XLabelFontFile = ''
    clip1Display.DataAxesGrid.YLabelFontFile = ''
    clip1Display.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
    clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
    clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
    clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    
    # hide streamline data in view, as we only need to see the output of the clip
    Hide(streamTracerWithCustomSource1, renderView1)
    
    # show color bar/color legend(s)
    clip1Display.SetScalarBarVisibility(renderView1, True)
    calculator1Display.SetScalarBarVisibility(renderView1, True)
    
    # Properties modified on clip1Display
    clip1Display.RenderLinesAsTubes = 1
    clip1Display.LineWidth = 3.0
    if opacity: clip1Display.Opacity = 0.75

    # Fix for colorbar location. Paraview seems to only let you set the location of one colorbar in the python...
    # That is, when one colorbar location is set the other colorbar jumps to the upper right corner
    bLUTColorBar.WindowLocation = 'UpperLeftCorner'
    bLUT.RescaleTransferFunction(Bmin, Bmax)
    totalpressureLUTColorBar.WindowLocation = 'UpperRightCorner'
    totalpressureLUT.RescaleTransferFunction(0, 80)

    # save image and state to file
    #-----------------------------------------Saving the image/video--------------------------------------------
    servermanager.SaveState('lastFig.pvsm')
    if properties['Movie'] == 'no':
        outname = os.path.splitext(os.path.split(pressurefile)[-1])[0]
        outname = os.path.join('images', outname.replace('_pressure', ''))
        SaveScreenshot(outname+'.png'. format(outname), magnification=1.0, quality=100, view=renderView1)
    else:
        SaveAnimation(filename = 'images/' + properties['Movie'] + '_movie.avi', FrameRate=2)

#========================================================================================================
if __name__=='__main__':
    read_config()

    #check for output directory
    figpath = 'images'
    if not os.path.isdir(figpath):
        os.mkdir(figpath)
    #temporary hardcode of input files - this needs to be moved back out to config
    vtxpath = os.path.abspath('vtk_files')
    pressurefile = os.path.join(vtxpath, 'restart_d20130317_t085100_pressure.vtp')
    fieldfile = os.path.join(vtxpath, 'restart_d20130317_t085100_field.vtu')
    pointsfile = os.path.join(vtxpath, 'sphere.vtp')
    pointsfile = os.path.join(vtxpath, 'disc.vtp')

    #generate either images or animations
    if properties['Movie'] == 'no':
        generateVisualization(pressurefile, fieldfile, pointsfile, opacity=False)
    else:
        raise NotImplementedError('Animation capability temporarily removed')

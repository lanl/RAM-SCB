## Visualization tools instructions

The visualization tools here generate VTK XML files from the RAM-SCB restart files.
They further provide a scripting capability for generating output images using Paraview's
`pvpython` shell (requires Paraview). Other VTK-based tools can be used with the
generated XML files.

### Making VTK XML files

The current visualization pipeline uses pressure data from RAM in the equatorial plane,
field data from SCB in 3D, and requires a set of input (seed) locations to trace field lines
through the SCB domain. The RAM plane and the seed locations are expected as VTK PolyData.
The SCB data are expected in VTK UnstructuredGrid format. All VTK XML files are expected to
be in a subdirectory called `vtk_files`
- Convenience routines to make VTK files with the seed locations are in `makeCustomSource.py`
- `convertRAMrestart.py` batch converts NetCDF restart files to `.vtp` and `.vtu` files
- `visualizeRAM.py` should be run in `pvpython` and will generate output PNG files

### Dependencies
The script to convert RAM-SCB restart files to VTK require `spacepy` and `netCDF4`.
These can both be obtained via `pip`.
As noted above, `pvpython` requires Paraview.

### Scripting visualizations

1. Put all netCDF restart files for which visulaizations are to be generated in the nc_files directory.
2. Modify config.txt as desired.
3. Execute ```source ram.sh``` (in bash shell) or ```source ram.csh``` (in csh or tcsh)
4. Look for the generated images in the images directory


#### Notes regarding the configuration file

Format: Anything before the colon is the property name, sometimes followed by possible values in parentheses. Anything after the colon is the property value. Property values are case-sensitive

Streamline source type: A sphere generates lesser streamlines. A slice of the magnetic field is a more appropriate choice
Plasma pressure display: Choose a species whose plasma pressure will be displayed
Camera Position: where the camera is
Camera Focal Point: where the camera is looking
Camera View Up: TODO
Scale: Indicate if you choose to see the scale
Movie: To generate images, set it to 'no'. To generate a movie, specify the groupname (first few characters that are common to all the restart files) whose movie is to be generated.
SaveState: Setting this to True will save the state in the images directory, so that the state can be loaded in the Paraview GUI later

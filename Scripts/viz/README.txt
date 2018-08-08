Execution instructions:
1) Put all netCDF restart files for which visulaizations are to be generated in the nc_files directory.
2) Modify config.txt as desired.
3) Execute 'bash ram.sh'
4) Look for the generated images in the images directory

Notes regarding the configuration file:

Format: Anything before the colon is the property name, sometimes followed by possible values in parantheses. Anything after the colon is the property value. Property values are case-sensitive


Streamline source type: A sphere generates lesser streamlines. A slice of the magnetic field is a more appropriate choice
Plasma pressure display: Choose a species whose plasma pressure will be displayed
Camera Position: where the camera is
Camera Focal Point: where the camera is looking
Camera View Up: TODO
Scale: Indicate if you choose to see the scale
Movie: To generate images, set it to 'no'. To generate a movie, specify the groupname (first few characters that are common to all the restart files) whose movie is to be generated.
SaveState: Setting this to True will save the state in the images directory, so that the state can be loaded in the Paraview GUI later
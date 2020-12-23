#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
heatfluxPointCloudcsv = CSVReader(FileName=['/home/tlooby/source/HEAT/rev2/data/nstx_116313/000851/HeatfluxPointCloud.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1490, 752]

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints9 = TableToPoints(Input=heatfluxPointCloudcsv)
tableToPoints9.XColumn = '# X'
tableToPoints9.YColumn = 'Y'
tableToPoints9.ZColumn = 'Z'

# show data in view
tableToPoints9Display = Show(tableToPoints9, renderView1)
# trace defaults for the display properties.
tableToPoints9Display.ColorArrayName = [None, '']

# set scalar coloring
ColorBy(tableToPoints9Display, ('POINTS', 'HeatFlux'))

# rescale color and/or opacity maps used to include current data range
tableToPoints9Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
tableToPoints9Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'HeatFlux'
heatFluxLUT = GetColorTransferFunction('HeatFlux')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
heatFluxLUT.ApplyPreset('Black-Body Radiation', True)

# get opacity transfer function/opacity map for 'HeatFlux'
heatFluxPWF = GetOpacityTransferFunction('HeatFlux')

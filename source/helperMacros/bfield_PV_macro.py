#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
b_pointcloudcsv = CSVReader(FileName=['/home/tlooby/source/HEAT/rev2/data/nstx_116313/000851/B_pointcloud.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1492, 752]

# get layout
viewLayout1 = GetLayout()

# Create a new 'SpreadSheet View'
#spreadSheetView1 = CreateView('SpreadSheetView')
#spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
#viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
#b_pointcloudcsvDisplay = Show(b_pointcloudcsv, spreadSheetView1)
# trace defaults for the display properties.
#b_pointcloudcsvDisplay.FieldAssociation = 'Row Data'

# destroy spreadSheetView1
#Delete(spreadSheetView1)
#del spreadSheetView1

# close an empty frame
#viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints6 = TableToPoints(Input=b_pointcloudcsv)
tableToPoints6.XColumn = '# X'
tableToPoints6.YColumn = 'Y'
tableToPoints6.ZColumn = 'Z'

# show data in view
tableToPoints6Display = Show(tableToPoints6, renderView1)
# trace defaults for the display properties.
tableToPoints6Display.ColorArrayName = [None, '']

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Calculator'
calculator3 = Calculator(Input=tableToPoints6)
# Properties modified on calculator3
calculator3.Function = '(iHat*Bx) + (jHat*By) + (kHat*Bz)'

# show data in view
calculator3Display = Show(calculator3, renderView1)
# trace defaults for the display properties.
calculator3Display.ColorArrayName = [None, '']

# hide data in view
Hide(tableToPoints6, renderView1)

# create a new 'Glyph'
glyph3 = Glyph(Input=calculator3,
    GlyphType='Arrow')
glyph3.Scalars = ['POINTS', 'Bx']
glyph3.Vectors = ['POINTS', 'Result']
glyph3.ScaleFactor = 35.50919392904
glyph3.GlyphTransform = 'Transform2'

# Properties modified on glyph3
glyph3.GlyphMode = 'Every Nth Point'
glyph3.Stride = 20

# get color transfer function/color map for 'Bx'
bxLUT = GetColorTransferFunction('GlyphVector')

# show data in view
glyph3Display = Show(glyph3, renderView1)
# trace defaults for the display properties.
glyph3Display.ColorArrayName = ['POINTS', 'GlyphVector']
glyph3Display.LookupTable = bxLUT

# show color bar/color legend
glyph3Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Bx'
bxPWF = GetOpacityTransferFunction('GlyphVector')

# set scalar coloring
ColorBy(glyph3Display, ('POINTS', 'GlyphVector'))

# rescale color and/or opacity maps used to include current data range
glyph3Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
glyph3Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'GlyphVector'
glyphVectorLUT = GetColorTransferFunction('GlyphVector')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
glyphVectorLUT.ApplyPreset('blue2yellow', True)

# get opacity transfer function/opacity map for 'GlyphVector'
glyphVectorPWF = GetOpacityTransferFunction('GlyphVector')

# hide data in view
Hide(calculator3, renderView1)

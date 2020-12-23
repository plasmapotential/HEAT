#fieldTrace_macro.py
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
structcsv = CSVReader(FileName=['/home/tlooby/source/HEAT/rev2/data/nstx_204118/001004/struct.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1492, 752]

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints7 = TableToPoints(Input=structcsv)
tableToPoints7.XColumn = '# X[m]'
tableToPoints7.YColumn = 'Y[m]'
tableToPoints7.ZColumn = 'Z[m]'

# show data in view
tableToPoints7Display = Show(tableToPoints7, renderView1)
# trace defaults for the display properties.
tableToPoints7Display.ColorArrayName = [None, '']

# create a new 'Programmable Filter'
programmableFilter2 = ProgrammableFilter(Input=tableToPoints7)
# Properties modified on programmableFilter2
programmableFilter2.Script = 'pdi = self.GetPolyDataInput()\npdo =  self.GetPolyDataOutput()\nnumPoints = pdi.GetNumberOfPoints()\npdo.Allocate()\nfor i in range(0, numPoints-1):\n    points = [i, i+1]\n    # VTK_LINE is 3\n    pdo.InsertNextCell(3, 2, points)'
programmableFilter2.RequestInformationScript = ''
programmableFilter2.RequestUpdateExtentScript = ''
programmableFilter2.PythonPath = ''

# show data in view
programmableFilter2Display = Show(programmableFilter2, renderView1)
# trace defaults for the display properties.
programmableFilter2Display.ColorArrayName = [None, '']

# hide data in view
Hide(tableToPoints7, renderView1)

# change solid color
programmableFilter2Display.DiffuseColor = [0.3333333333333333, 1.0, 0.0]

# Properties modified on programmableFilter2Display
programmableFilter2Display.LineWidth = 2.0

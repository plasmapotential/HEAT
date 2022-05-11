#this macro is to load multiple .foam files that all exist in a single
#heatFoam directory.  When HEAT creates .foam files for multiple PFCs
#it can be a pain to import all of them.  This script loops through
#all the subdirectories in heatFoamPath and opens the .foam files therein
#edit heatFoamPath, then import the Macro in paraview and run
#
# trace generated using paraview version 5.10.1

#### import the simple module from the paraview
from paraview.simple import *
import os

#===user edits these variables:
#peak temperature in colorbar (you can edit later in GUI)
maxT = 1500
#directory where heatFoam has output to
heatFoamPath = '/home/tom/HEAT/data/sparc_000001_sweep7/openFoam/heatFoam/'
PFCnames = [f.name for f in os.scandir(heatFoamPath) if f.is_dir()]


for PFC in PFCnames:
    regName = PFC + '.foam'
    f = heatFoamPath + PFC + '/' + regName


    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'OpenFOAMReader'
    foam = OpenFOAMReader(registrationName=regName, FileName=f)
    foam.MeshRegions = ['internalMesh']
    foam.CellArrays = ['DT', 'HF', 'T', 'fluxMag', 'gradTx', 'gradTy', 'gradTz', 'qFVM', 'qINV', 'thermCond']

    # get animation scene
    animationScene1 = GetAnimationScene()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    foamDisplay = Show(foam, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    foamDisplay.Representation = 'Surface'
    foamDisplay.ColorArrayName = [None, '']
    foamDisplay.SelectTCoordArray = 'None'
    foamDisplay.SelectNormalArray = 'None'
    foamDisplay.SelectTangentArray = 'None'
    foamDisplay.OSPRayScaleArray = 'T'
    foamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    foamDisplay.SelectOrientationVectors = 'None'
    foamDisplay.ScaleFactor = 21.775000000000002
    foamDisplay.SelectScaleArray = 'None'
    foamDisplay.GlyphType = 'Arrow'
    foamDisplay.GlyphTableIndexArray = 'None'
    foamDisplay.GaussianRadius = 1.08875
    foamDisplay.SetScaleArray = ['POINTS', 'T']
    foamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    foamDisplay.OpacityArray = ['POINTS', 'T']
    foamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    foamDisplay.DataAxesGrid = 'GridAxesRepresentation'
    foamDisplay.PolarAxes = 'PolarAxesRepresentation'
    foamDisplay.ScalarOpacityUnitDistance = 4.139798151157404
    foamDisplay.OpacityArrayName = ['POINTS', 'T']

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    foamDisplay.ScaleTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, maxT, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    foamDisplay.OpacityTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, maxT, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(foamDisplay, ('POINTS', 'T'))

    # rescale color and/or opacity maps used to include current data range
    foamDisplay.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    foamDisplay.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')

    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')

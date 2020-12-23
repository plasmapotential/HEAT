#csv2vtk.py
#Description:   csv vector data converter for various HEAT outputs
#Engineer:      T Looby
#Date:          20200205

import sys
import textwrap
import csv
#this step is completed in dashGUI.py
#path = '/opt/paraview/ParaView-5.7.0-MPI-Linux-Python3.7-64bit/lib/python3.7/site-packages'
#sys.path.append(path)


def convert2Glyph(pcfile, prefix):
    """
    Converts csv file to vtk file using paraview.simple for use in HEAT
    Then converts this vtkPolyData object into glyph vectors that correspond
    to the magnitudes of the vector components.
    """
    import os
    dir, file = os.path.split(pcfile)

    #Make a new directory specifically for Paraview VTK results, if it doesn't
    #already exist
    PVdir = dir+"/paraview"
    if not os.path.exists(PVdir):
        os.makedirs(PVdir)

    #Read first line from file and get the vector component labels
    with open(pcfile, newline='') as f:
        reader = csv.reader(f)
        CSVheader = next(reader)


    from paraview.simple import CSVReader, TableToPoints, CreateWriter, Calculator, Glyph
    #Generate Table
    pc = CSVReader(FileName=[pcfile])

    #Convert To Points
    t2p = TableToPoints(Input=pc)
    t2p.XColumn = '# X'
    t2p.YColumn = 'Y'
    t2p.ZColumn = 'Z'

    #Calculate Vector Magnitudes
    calculator = Calculator(Input=t2p)
    calculator.Function = '(iHat*'+CSVheader[3]+') + (jHat*'+CSVheader[4]+') + (kHat*'+CSVheader[5]+')'
    #calculator.Function = '(iHat*Bx) + (jHat*By) + (kHat*Bz)'

    #Create Glyphs
    glyph = Glyph(Input=calculator, GlyphType='Arrow')
    glyph.ScaleFactor = 20.0
    glyph.GlyphTransform = 'Transform2'
    glyph.GlyphMode = 'Every Nth Point'
    glyph.Stride = 20
    glyph.OrientationArray = ['POINTS', 'Result']
    glyph.ScaleArray = ['POINTS', 'Result']


    #Write to file
    writer = CreateWriter(PVdir+"/"+prefix+".vtk", glyph)
    writer.UpdatePipeline()
    del writer

    return


def convert2Trace(pcfile, prefix):
    """
    Converts csv file to vtk file using paraview.simple for use in HEAT
    Then converts this vtkPolyData object into a field line trace
    """
    import os
    dir, file = os.path.split(pcfile)

    #Make a new directory specifically for Paraview VTK results, if it doesn't
    #already exist
    PVdir = dir+"/paraview"
    if not os.path.exists(PVdir):
        os.makedirs(PVdir)

    from paraview.simple import CSVReader, TableToPoints, CreateWriter, ProgrammableFilter

    #Generate Table
    pc = CSVReader(FileName=[pcfile])

    #Convert To Points
    t2p = TableToPoints(Input=pc)
    t2p.XColumn = '# X[mm]'
    t2p.YColumn = 'Y[mm]'
    t2p.ZColumn = 'Z[mm]'

    #Programmable Filter
    progFil = ProgrammableFilter(Input=t2p)
    progFil.Script = ( 'pdi = self.GetPolyDataInput()\npdo =  self.GetPolyDataOutput()\nnumPoints = pdi.GetNumberOfPoints()\npdo.Allocate()\nfor i in range(0, numPoints-1):\n    points = [i, i+1]\n    # VTK_LINE is 3\n    pdo.InsertNextCell(3, 2, points)' )

    #Write to file
    writer = CreateWriter(PVdir+"/"+prefix+".vtk", progFil)
    writer.UpdatePipeline()
    del writer
    return

def convert2PointCloud(pcfile, prefix):
    """
    Converts csv file to vtk file using paraview.simple for use in HEAT
    """
    import os
    dir, file = os.path.split(pcfile)

    #Make a new directory specifically for Paraview VTK results, if it doesn't
    #already exist
    PVdir = dir+"/paraview"
    if not os.path.exists(PVdir):
        os.makedirs(PVdir)

    from paraview import simple
    #Import PC
    pc = simple.ParticlesReader(FileName=[pcfile])

    #Write to file
    writer = simple.CreateWriter(PVdir+"/"+prefix+".vtk", pc)
    writer.UpdatePipeline()


    #from paraview import simple
    ## create a new 'Table To Points'
    #hfPC = simple.CSVReader(FileName=[pcfile])
    #t2p = simple.TableToPoints(Input=hfPC)
    #t2p.XColumn = '# X'
    #t2p.YColumn = 'Y'
    #t2p.ZColumn = 'Z'
    #renderView1 = simple.GetActiveViewOrCreate('RenderView')
    #simple.SetActiveView(renderView1)
    #display = simple.Show(t2p, renderView1)
    #display.ColorArrayName = [None, '']
    #simple.ColorBy(display, ['POINTS', 'ShadowMask'])
    #display.RescaleTransferFunctionToDataRange(True)
    #display.SetScalarBarVisibility(renderView1, True)
    #heatFluxLUT = simple.GetColorTransferFunction('ShadowMask')
    #heatFluxLUT.ApplyPreset('Black-Body Radiation', True)
    #simple.SaveState(PVdir+"/"+prefix+".py")

    return



if __name__=='__main__':
    import sys
    import os
    from subprocess import check_output
    import signal
    print('Converting CSV to VTK format')
    if sys.argv[2] == 'glyph':
        convert2Glyph(sys.argv[1], sys.argv[3])
    elif sys.argv[2] == 'trace':
        convert2Trace(sys.argv[1], sys.argv[3])
    else:
        convert2PointCloud(sys.argv[1], sys.argv[3])
    print("Finished")

    #sometimes pvpython hangs, so we kill by pid after our work is done
    #pid = int(check_output(["pidof","pvpython"]))
    #os.kill(pid, signal.SIGKILL)

#vtkAlgorithms.py
#Description:   Various VTK algorithms for use with paraviewweb visualizer
#Engineer:      T Looby
#Date:          20200205

#Note:  This script should be run in combination with a PV visualization launcher

from paraview.util.vtkAlgorithm import *
#==============================================================================
#                   Read vtk files only
#==============================================================================
#@smproxy.reader(name="vtkReader3D", label="vtk TEST",
#                extensions="vtk", file_description="Vtk files")
#class CSVReader3D(VTKPythonAlgorithmBase):
#    def __init__(self):
#        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkPolyData')
#        self._filename = None
#        self._filename = '/u/tlooby/source/HEAT/rev4/data/204118/001004/B_pointcloud.vtk'
#
#        from vtkmodules.vtkIOLegacy import vtkGenericDataObjectReader
#        from vtkmodules.vtkFiltersProgrammable import vtkProgrammableFilter
#        from vtkmodules.vtkFiltersCore import vtkGlyph3D
#        from vtkmodules.vtkFiltersSources import vtkConeSource
#
#        self._reader = vtkGenericDataObjectReader()
#        self._reader.SetFileName(self._filename)
#        self._reader.Update()
#
#        sources = vtkConeSource()
#        sources.SetResolution(6)
#        sources.Update()
#
#        self._glyph = vtkGlyph3D()
#        self._glyph.SetInputConnection(self._reader.GetOutputPort())
#        self._glyph.SetSourceConnection(sources.GetOutputPort())
#        self._glyph.ScalingOn()
#        self._glyph.SetScaleModeToScaleByScalar()
#        self._glyph.SetVectorModeToUseVector()
#        self._glyph.OrientOn()
#        # Set the overall (multiplicative) scaling factor
#        self._glyph.SetScaleFactor(1)
#
#
#
#
#    @smproperty.stringvector(name="FileName")
#    @smdomain.filelist()
#    @smhint.filechooser(extensions="csv", file_description="Numpy CSV files")
#    def SetFileName(self, name):
#        """Specify filename for the file to read."""
#        if self._filename != name:
#            self._filename = name
#            self.Modified()
#
#    def RequestData(self, request, inInfo, outInfo):
#        from vtkmodules.vtkCommonDataModel import vtkPolyData
#        from vtkmodules.vtkFiltersCore import vtkGlyph3D
#        self._glyph.Update()
#        #output = vtkPolyData.GetData(outInfo, 0)
#        #output.ShallowCopy(self._glyph.GetOutput())
#        #output.Update()
#        output = vtkGlyph3D.GetData(outInfo, 0)
#        output.ShallowCopy(self._glyph.GetOutput())
#        output.Update()
#        return 1

#==============================================================================
#                   Read csv with vtkDelimitedTextReader
#==============================================================================
@smproxy.reader(name="CSVReader3D", label="CSV TEST",
                extensions="csv", file_description="CSV files")
class CSVReader3D(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkPolyData')
        self._filename = None
        self._filename = '/u/tlooby/source/HEAT/rev4/data/204118/001004/B_pointcloud.csv'
        from vtkmodules.vtkIOInfovis import vtkDelimitedTextReader
        from vtkmodules.vtkFiltersGeneral import vtkTableToPolyData

        # Create a Delimited Text Reader object
        self._csv_source = vtkDelimitedTextReader()
        self._csv_source.SetFieldDelimiterCharacters(",")
        self._csv_source.SetHaveHeaders(True)
        self._csv_source.SetDetectNumericColumns(True)
        self._csv_source.SetFileName(self._filename)

        #Convert Text Delimited Object to PolyData (points)
        self._realAlgorithm = vtkTableToPolyData()
        self._realAlgorithm.AddInputConnection(self._csv_source.GetOutputPort())
        self._realAlgorithm.SetXColumnIndex(0)
        self._realAlgorithm.SetYColumnIndex(1)
        self._realAlgorithm.SetZColumnIndex(2)
        self._realAlgorithm.SetXComponent(3)
        self._realAlgorithm.SetYComponent(4)
        self._realAlgorithm.SetZComponent(5)
        self._realAlgorithm.SetXColumn('# X')
        self._realAlgorithm.SetYColumn('Y')
        self._realAlgorithm.SetZColumn('Z')


    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="csv", file_description="Numpy CSV files")
    def SetFileName(self, name):
        """Specify filename for the file to read."""
        if self._filename != name:
            self._filename = name
            self.Modified()

    def RequestData(self, request, inInfo, outInfo):
        from vtkmodules.vtkCommonDataModel import vtkPolyData
        self._realAlgorithm.Update()
        output = vtkPolyData.GetData(outInfo, 0)
        output.ShallowCopy(self._realAlgorithm.GetOutput())
        output.Update()
        return 1

#==============================================================================
#                   Read csv with numpy
#==============================================================================
#@smproxy.reader(name="CSVReader3D", label="CSV TEST",
#                extensions="csv", file_description="CSV files")
#class CSVReader3D(VTKPythonAlgorithmBase):
#    def __init__(self):
#        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkPolyData')
#        self._filename = None
#        self._filename = '/u/tlooby/source/HEAT/rev4/data/204118/001004/B_pointcloud.csv'
#
#        import numpy as np
#        from paraview.vtk.util import numpy_support
#        from vtkmodules.vtkCommonDataModel import vtkTable
#        from paraview.vtk import VTK_DOUBLE
#        from vtkmodules.vtkFiltersGeneral import vtkTableToPolyData
#
#        arr = np.genfromtxt(self._filename, dtype=np.dtype(np.float64), skip_header=0, names=True, delimiter=',', autostrip=True)
#        tab = vtkTable()
#        for name in arr.dtype.names:
#            vtkarr = numpy_support.numpy_to_vtk( arr[name], deep=True, array_type=VTK_DOUBLE )
#            vtkarr.SetName(name)
#            tab.AddColumn(vtkarr)
#
#        #Convert Text Delimited Object to PolyData (points)
#        self._realAlgorithm = vtkTableToPolyData()
#        self._realAlgorithm.AddInputDataObject(0, tab)
#        self._realAlgorithm.SetXColumnIndex(0)
#        self._realAlgorithm.SetYColumnIndex(1)
#        self._realAlgorithm.SetZColumnIndex(2)
#        self._realAlgorithm.SetXComponent(3)
#        self._realAlgorithm.SetYComponent(4)
#        self._realAlgorithm.SetZComponent(5)
#        self._realAlgorithm.SetXColumn('X')
#        self._realAlgorithm.SetYColumn('Y')
#        self._realAlgorithm.SetZColumn('Z')
#
#    @smproperty.stringvector(name="FileName")
#    @smdomain.filelist()
#    @smhint.filechooser(extensions="csv", file_description="Numpy CSV files")
#    def SetFileName(self, name):
#        """Specify filename for the file to read."""
#        if self._filename != name:
#            self._filename = name
#            self.Modified()
#
#    def RequestData(self, request, inInfo, outInfo):
#        from vtkmodules.vtkCommonDataModel import vtkPolyData
#        self._realAlgorithm.Update()
#        output = vtkPolyData.GetData(outInfo, 0)
#        output.ShallowCopy(self._realAlgorithm.GetOutput())
#        return 1
#

#vtkOpsClass.py
#Description:   HEAT module to create vtk objects (.vtk,.vtp, etc.)
#Engineer:      T Looby
#Date:          20221108
import vtk
from vtk.util.numpy_support import vtk_to_numpy as vtk2np
import sys
import numpy as np

import logging
log = logging.getLogger(__name__)

class VTKops:
    def __init__(self):
        return

    def initializeMeshScalar(self, mesh, scalar, label):
        """
        initializes a HEAT vtp mesh object which requires a
        freecad mesh object and a scalar numpy array defined on the mesh

        scalar should be the length of mesh triangles in the mesh

        label is the name of the scalar that will be in the colorbar in paraview
        """
        self.mesh = mesh
        self.scalar = scalar
        self.Npts = len(scalar)
        self.label = label
        return

    def initializePointCloudScalar(self, ctrs, scalar, label):
        """
        initializes a HEAT vtp point cloud object which requires a
        freecad mesh object and a scalar numpy array defined on the mesh

        scalar should be the length of mesh triangles in the mesh

        label is the name of the scalar that will be in the colorbar in paraview
        """
        self.ctrs = ctrs
        self.scalar = scalar
        self.Npts = len(scalar)
        self.label = label
        return

    def initializeVectorField(self, ctrs, vecs, label):
        """
        initializes a HEAT vtp glyph object which requires a
        freecad mesh object and a vector numpy array defined on the mesh

        vec is a vector field of the format Nx,Ny,Nz, and is Nx3 in shape

        label is the name of the scalar that will be in the colorbar in paraview
        """
        self.ctrs = ctrs
        self.vecs = vecs
        self.Npts = len(vecs)
        self.label = label
        return

    def writeMeshVTP(self, outFile):
        """
        writes a .vtp file that contains the mesh triangle data, and has the
        scalar field data stored in the Colors object.
        """
        # setup colors
        Colors = vtk.vtkFloatArray()
        #Colors.SetNumberOfComponents(3)
        Colors.SetNumberOfTuples(self.Npts)
        Colors.SetName(self.label) #can change to any string

        #points
        vtkPts = vtk.vtkPoints()

        #build points and colors
        for i,facet in enumerate(self.mesh.Facets):
            for j in range(3):
                x = facet.Points[j][0]
                y = facet.Points[j][1]
                z = facet.Points[j][2]
                vtkPts.InsertNextPoint(x,y,z)
                #        Colors.InsertTuple( i*3+j, (arr[i],arr[i],arr[i]) )
                Colors.InsertTuple( i*3+j, [self.scalar[i]] )

        #build vtp triangular mesh
        Triangles = vtk.vtkCellArray()
        for i in range(self.Npts):
            Triangle = vtk.vtkTriangle()
            Triangle.GetPointIds().SetId(0, i*3+0)
            Triangle.GetPointIds().SetId(1, i*3+1)
            Triangle.GetPointIds().SetId(2, i*3+2)
            Triangles.InsertNextCell(Triangle)

        #build final vtp object for writing
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtkPts)
        polydata.SetPolys(Triangles)
        polydata.GetPointData().SetScalars(Colors)
        polydata.Modified()
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(outFile)
        writer.SetInputData(polydata)
        #writer.SetDataModeToBinary()
        writer.Write()

        return

    def writePointCloudVTP(self, outFile):
        """
        writes a point cloud vtk object.

        point cloud has scalar array defined as the color
        point cloud is defined on a set of x,y,z coordinates (ctrs)
        """
        #points
        vtkPts = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        # setup colors
        Colors = vtk.vtkFloatArray()
        #Colors.SetNumberOfComponents(3)
        Colors.SetNumberOfTuples(self.Npts)
        Colors.SetName(self.label) #can change to any string

        for i in range(self.Npts):
            x = self.ctrs[i,0]
            y = self.ctrs[i,1]
            z = self.ctrs[i,2]
            id = vtkPts.InsertNextPoint(x,y,z)
            cells.InsertNextCell(1)
            cells.InsertCellPoint(id)
            Colors.InsertTuple( i, [self.scalar[i]] )


        #build final vtp object for writing
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtkPts)
        polydata.SetVerts(cells)
        polydata.GetPointData().SetScalars(Colors)
        polydata.Modified()

        writer = vtk.vtkXMLPolyDataWriter()
        writer.DebugOn()
        writer.SetFileName(outFile)
        writer.SetInputData(polydata)
        #writer.SetDataModeToBinary()
        writer.Write()

        return

    def writeGlyphVTP(self, f, scale=0.005):
        """
        """
        vtkPts = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        vec = vtk.vtkDoubleArray()
        vec.SetName("vec")
        vec.SetNumberOfComponents(3)
        vec.SetNumberOfTuples(self.Npts)

        mag = vtk.vtkDoubleArray()
        mag.SetNumberOfValues(self.Npts)
        mag.SetName(self.label)

        m = np.linalg.norm(self.vecs, axis=1)

        for i in range(self.Npts):
            vec.SetTuple(i, list(self.vecs[i]))
            mag.SetValue(i, m[i])
            id = vtkPts.InsertNextPoint(self.ctrs[i,0],self.ctrs[i,1],self.ctrs[i,2])
            cells.InsertNextCell(1)
            cells.InsertCellPoint(id)

        # Add to point data array.
        poly = vtk.vtkPolyData()
        poly.SetPoints(vtkPts)
        poly.SetVerts(cells)
        poly.GetPointData().AddArray(vec)
        poly.GetPointData().AddArray(mag)
        poly.GetPointData().SetActiveScalars(self.label)
        poly.GetPointData().SetActiveVectors("vec")
        poly.Modified()

        # Create glyph
        arrow = vtk.vtkArrowSource()
        arrow.Update()
        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(poly)
        glyph.SetSourceConnection(arrow.GetOutputPort())
        glyph.SetScaleFactor(scale)
        glyph.OrientOn()
        glyph.SetVectorModeToUseVector()
        glyph.SetColorModeToColorByScalar()
        glyph.Update()

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(f)
        writer.SetInputData(glyph.GetOutput())
        #writer.SetDataModeToBinary()
        writer.Write()


        return

    def writeTraceVTP(self, xyz, f):
        """
        writes a field line trace for paraview where xyz are the coordinates
        and f is the file to write
        """
        points = vtk.vtkPoints()
        for x, y, z in xyz:
            points.InsertNextPoint(x, y, z)

        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(xyz))
        for i in range(len(xyz)):
            polyline.GetPointIds().SetId(i, i)

        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polyline)

        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(cells)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.DebugOn()
        writer.SetFileName(f)
        writer.SetInputData(polyData)
        #writer.SetDataModeToBinary()
        writer.Write()

        return
    
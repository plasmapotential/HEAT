#vtkOpsClass.py
#Description:   HEAT module to create vtk objects (.vtk,.vtp, etc.)
#Engineer:      T Looby
#Date:          20221108
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
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

    def writeMeshVTP(self, outFile, binary=True, compress=False):
        """
        writes a .vtp file that contains the mesh triangle data, and has the
        scalar field data stored in the Colors object.
        """
        facets = self.mesh.Facets
        ntri = len(facets)
        if ntri == 0:
            poly = vtk.vtkPolyData()
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(outFile)
            writer.SetInputData(poly)
            if binary: writer.SetDataModeToBinary()
            if not compress: writer.SetCompressorTypeToNone()
            writer.Write()
            return

        pts = np.empty((ntri * 3, 3), dtype=np.float32)
        for i, f in enumerate(facets):
            pts[3*i:3*i+3, :] = np.asarray(f.Points, dtype=np.float32)

        conn = np.arange(ntri * 3, dtype=np.int64)

        cells = vtk.vtkCellArray()
        if vtk.VTK_MAJOR_VERSION >= 9:
            offsets = np.arange(0, ntri * 3, 3, dtype=np.int64)
            cells.SetData(
                numpy_to_vtkIdTypeArray(offsets, deep=False),
                numpy_to_vtkIdTypeArray(conn,   deep=False)
            )
        else:
            ids = np.empty(ntri * 4, dtype=np.int64)
            ids[0::4] = 3
            ids[1::4] = conn[0::3]
            ids[2::4] = conn[1::3]
            ids[3::4] = conn[2::3]
            cells.SetCells(ntri, numpy_to_vtkIdTypeArray(ids, deep=False))

        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(numpy_to_vtk(pts, deep=False))

        # --- scalar per triangle (cell data) ---
        cell_scal = numpy_to_vtk(np.asarray(self.scalar, dtype=np.float32), deep=False)
        cell_scal.SetName(self.label)

        poly = vtk.vtkPolyData()
        poly.SetPoints(vtk_points)
        poly.SetPolys(cells)
        poly.GetCellData().SetScalars(cell_scal)   # <— key change (CellData)
        poly.Modified()
        # -------- Write (binary + no compression = fastest) -------
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(outFile)
        writer.SetInputData(poly)
        if binary:
            writer.SetDataModeToBinary()
        if not compress:
            writer.SetCompressorTypeToNone()
        writer.Write()


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
    
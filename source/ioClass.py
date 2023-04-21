#ioClass.py
#Description:   class for handling input / output pipelines
#Engineer:      T Looby
#Date:          20221108
import os
import numpy as np
import GUIscripts.vtkOpsClass as vtkOpsClass

import logging
log = logging.getLogger(__name__)

class IO_HEAT:
    def __init__(self):
        return

    def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
        """
        sets up pythonic string number formats for shot and timesteps
        """
        self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
        self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
        return

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization

        Here is a list of variables with description:
        testvar         dummy for testing

        """

        self.allowed_vars = [
                            'vtpPCOut',
                            'vtpMeshOut',
                            'csvOut',
                            ]
        return


    def setTypes(self):
        trueList = ['t','T','true','True','TRUE','Tru','TRU']
        try:
            if self.vtpPCOut in trueList:
                self.vtpPCMask = True
            else:
                self.vtpPCMask = False

            if self.vtpMeshOut in trueList:
                self.vtpMeshMask = True
            else:
                self.vtpMeshMask = False

            if self.csvOut in trueList:
                self.csvMask = True
            else:
                self.csvMask = False
        except:
            print("Error in input file.  No output filetypes provided.")
            print("Defaulting to csv and vtp")
            log.info("Error in input file.  No output filetypes provided.")
            log.info("Defaulting to csv and vtp")
            self.csvMask = True
            self.vtpMeshMask = True
            self.vtpPCMask = False
        return

    def outputMasks(self, list):
        """
        sets output filestream masks depending upon what is in list

        list comes from gui and corresponds to each flag below
        """
        self.vtpMeshMask = False
        self.vtpPCMask = False
        self.csvMask = False

        if "vtpMesh" in list:
            self.vtpMeshMask = True
        if "vtpPC" in list:
            self.vtpPCMask = True
        if "csv" in list:
            self.csvMask = True

        print("\n---Output file settings---")
        print("vtpMeshMask: " + str(self.vtpMeshMask))
        print("vtpPCMask: " + str(self.vtpPCMask))
        print("csvMask: " + str(self.csvMask))
        log.info("\n---Output file settings---")
        log.info("vtpMeshMask: " + str(self.vtpMeshMask))
        log.info("vtpPCMask: " + str(self.vtpPCMask))
        log.info("csvMask: " + str(self.csvMask))

        return

    def writeMeshVTP(self, mesh, scalar, label, prefix, path, tag=None, PClabel=True):
        """
        writes a vtp mesh file
        output file contains the PFC mesh, and an array of heat flux values,
        scalar, that is indexed to match the mesh triangles.

        mesh - freecad mesh object
        scalar - array of scalar values that are indexed to match mesh triangles
        label - units/label that will be displayed in paraview colorbar
        tag - name tag for file (<tag>.vtp)
        path - file path where paraview folder lives
        """
        if tag is None:
            if PClabel == True:
                fName = prefix + '_mesh.vtp'
            else:
                fName = prefix + '.vtp'
        else:
            if PClabel == True:
                fName = prefix + '_'+tag+'_mesh.vtp'
            else:
                fName = prefix + '_' + tag + '.vtp'

        VTKops = vtkOpsClass.VTKops()
        VTKops.initializeMeshScalar(mesh, scalar, label)

        PVdir = path + "paraview/"
        f = PVdir+fName
        if not os.path.exists(PVdir):
            os.makedirs(PVdir)

        VTKops.writeMeshVTP(f)
        return

    def writePointCloudVTP(self, ctrs, scalar, label, prefix, path, tag=None, PClabel=True):
        """
        writes a vtk point cloud file
        output file contains the PFC mesh centers, and an array of heat flux values,
        scalar, that is indexed to match the mesh triangles.

        ctrs - centroids of freecad mesh object triangles in [m]
        scalar - array of scalar values that are indexed to match mesh triangles
        label - units/label that will be displayed in paraview colorbar
        tag - name tag for file (<tag>.vtp)
        path - file path where paraview folder lives
        """
        if tag is None:
            if PClabel == True:
                fName = prefix + '_PC.vtp'
            else:
                fName = prefix + '.vtp'
        else:
            if PClabel == True:
                fName = prefix + '_'+tag+'_PC.vtp'
            else:
                fName = prefix + '_' + tag + '.vtp'

        VTKops = vtkOpsClass.VTKops()
        VTKops.initializePointCloudScalar(ctrs*1000.0, scalar, label) #scale to mm

        PVdir = path + "paraview/"
        f = PVdir+fName
        if not os.path.exists(PVdir):
            os.makedirs(PVdir)

        VTKops.writePointCloudVTP(f)
        return

    def writeGlyphVTP(self, ctrs, vecs, label, prefix, path, tag=None):
        """
        writes a vtk glyph file
        output file contains the PFC mesh centers, and an array of heat flux values,
        scalar, that is indexed to match the mesh triangles.

        ctrs - centroids of freecad mesh object triangles
        scalar - array of scalar values that are indexed to match mesh triangles
        label - units/label that will be displayed in paraview colorbar
        tag - name tag for file (<tag>.vtp)
        path - file path where paraview folder lives
        """
        if tag is None:
            fName = prefix + '.vtp'
        else:
            fName = prefix + '_'+tag+'.vtp'

        VTKops = vtkOpsClass.VTKops()
        VTKops.initializeVectorField(ctrs, vecs, label)

        PVdir = path + "paraview/"
        f = PVdir+fName
        if not os.path.exists(PVdir):
            os.makedirs(PVdir)

        VTKops.writeGlyphVTP(f)
        return


    def writePointCloudCSV(self,centers,Scalar,dataPath,label,tag=None,prefix='Undefined'):
        """
        writes a point cloud CSV file
        """
        print("Writing point cloud CSV for "+prefix)
        log.info("Writing point cloud CSV for "+prefix)
        if tag is None:
            pcfile = dataPath + prefix + '.csv'
        else:
            pcfile = dataPath + prefix + '_'+tag+'.csv'
        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = Scalar
        head = "X,Y,Z,"+label
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)
        return

    def writeGlyphCSV(self,centers,vecs,path,prefix,head,tag=None):
        """
        In paraview use TableToPoints => Calculator => Glyph
        Calculator should have this formula:
        (iHat*Nx) + (jHat*Ny) + (kHat*Nz)
        where Nx,Ny,Nz are labels in header
        """
        print("Creating Glyph "+prefix)
        log.info("Creating Glyph "+prefix)
        if tag is None:
            pcfile = path + prefix +'.csv'
        else:
            pcfile = path + prefix + '_' +tag+ '.csv'

        data = np.zeros((len(centers), 6))
        data[:,0] = centers[:,0]*1000.0
        data[:,1] = centers[:,1]*1000.0
        data[:,2] = centers[:,2]*1000.0
        data[:,3] = vecs[:,0]
        data[:,4] = vecs[:,1]
        data[:,5] = vecs[:,2]
        np.savetxt(pcfile, data, delimiter=',',fmt='%.10f', header=head)
        return

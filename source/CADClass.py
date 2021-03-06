#CADClass.py
#Description:   Base HEAT CAD module
#Engineer:      T Looby
#Date:          20191107

import sys
import os

#this happens in dashGUI.py now
#FREECADPATH = '/opt/freecad/appImage/squashfs-root/usr/lib'
#oldpath = sys.path
#sys.path.append(FREECADPATH)
#sys.path = [FREECADPATH]
import FreeCAD
import Part
import Mesh

import MeshPart
#sys.path = oldpath
import Import
import stl


import time
import numpy as np
import pandas as pd

import toolsClass
tools = toolsClass.tools()

import logging
log = logging.getLogger(__name__)

class CAD:
    """
    General CAD class

    Inputs:
    mode    GUI or CMD
    infile  input file in CSV format with each variable on a new line:
            variable, value

    Generally, there are part objects, which correspond to parts from an STEP
    (ISO 10303-21) file.  There are also mesh objects, which correspond to a
    mesh from an STL file.  The region of interest (ROI) is a list of part
    numbers that we are interested in.  When we get a ROI, we then build
    corresponding part (self.ROIparts), mesh (self.ROImeshes), face normal
    (self.ROInorms), and face center (self.ROIctrs) objects, which correspond
    to the ROI list by index.

    Maybe this picture will help:

                            ROI
                    _________|_______________
                   /     /       \     \     \
                parts  meshes   ctrs  norms  areas


    """
    def __init__(self, rootDir, dataPath):
        """
        dataPath is the location where we write all output to
        """
        self.rootDir = rootDir
        tools.rootDir = self.rootDir
        self.dataPath = dataPath
        tools.dataPath = self.dataPath
        return

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization

        Here is a list of variables with description:
        testvar         dummy for testing
        STPfile         STEP file location (ISO 10303-21)
        STLpath         location to save STLs
        rootDir         Root of this HEAT run
        ROIGridRes      max length of mesh triangle edge in [mm] on ROI
        gridRes         max length of mesh triangle edge in [mm] not in ROI
        FreeCADPath     path to freecad .so files on local machine
        permute_mask    boolean indicating if CAD file had z as 'upward' direction
                        NSTXU engineers have y as vertical coordinate for some reason
        unitConvert     Scalar that the input vectors from STL/STP can be multiplied
                        with to get units in meters.  If input is in inches, then
                        unitConvert = 0.0254, if input is in meters, then
                        unitConvert = 1.0
        """


        self.allowed_vars = [
                            'rootDir',
                            'STPfile',
                            #'STLpath',
                            'ROIGridRes',
                            'gridRes',
                            'FreeCADPath',
                            'PVPath',
                            'permute_mask',
                            'unitConvert',
                            'assembly_mask']
        return

    def setTypes(self):
        """
        Nothing to do for this class
        """
        pass


    def getROI(self, timestepMap):
        """
        Writes ROI as list to CAD object.  Input is timestepMap dataframe
        which is read by function in PFCClass.
        """
        self.ROI = timestepMap['PFCname'].values
        self.ROIList = list(set(self.ROI))
        self.ROIparts = ['None' for i in range(len(self.ROI))]
        self.ROImeshes = ['None' for i in range(len(self.ROI))]
        self.ROIctrs = ['None' for i in range(len(self.ROI))]
        self.ROInorms = ['None' for i in range(len(self.ROI))]
        return

    def getIntersectsFromFile(self, timestepMap):
        """
        Writes intersections to CAD object
        """
        self.ROImapDirections = []
        self.ROIintersects = []

        #parse intersects column and make indexed list of intersection parts
        # for each ROI part
        for row in timestepMap['intersectName']:
            self.ROIintersects.append(row.split(':'))

        #list of unique intersect parts
        self.intersectList = list(set( [j for row in self.ROIintersects for j in row] ))

        #initialize intersect variables
        self.intersectParts = ['None' for i in range(len(self.intersectList))]
        self.intersectMeshes = ['None' for i in range(len(self.intersectList))]
        self.intersectCtrs = ['None' for i in range(len(self.intersectList))]
        self.intersectNorms = ['None' for i in range(len(self.intersectList))]

        return


    def getROImeshes(self, resolution=None):
        """
        Checks to see if STLs at desired resolution exist.  If they do, load em.
        If they don't, create them.
        """
        if resolution == None:  resolution=self.ROIGridRes
        for idx,partnum in enumerate(self.ROI):
            name = self.STLpath + partnum + "___" + resolution +"mm.stl"
            if os.path.exists(name):
                print("Mesh exists, loading...")
                self.loadROIMesh(name,idx)
            else:
                print("New mesh.  Creating...")
                self.ROIobjFromPartnum(partnum,idx)

        #Now get face centers, normals, areas
        self.ROInorms,self.ROIctrs,self.ROIareas = self.normsCentersAreas(self.ROImeshes)
        return

    def getIntersectMeshes(self, resolution=None):
        """
        Checks to see if STLs at desired resolution exist.  If they do, load em.
        If they don't, create them.
        """
        if resolution == None:  resolution=self.gridRes
        print("TEST")
        print(resolution)
        for partnum in self.intersectList:
            name = self.STLpath + partnum + "___" + resolution +"mm.stl"
            if os.path.exists(name):
                print("Mesh exists, loading...")
                self.loadIntersectMesh(name)
            else:
                print("New mesh.  Creating...")
                self.intersectObjFromPartnum(partnum, resolution)

        #Now get face centers, normals, areas
        self.intersectNorms,self.intersectCtrs,self.intersectAreas = self.normsCentersAreas(self.intersectMeshes)
        return

    def ROIobjFromPartnum(self, partslist, idx):
        """
        Generates ROI objects from list of part numbers.
        """
        #Check if this is a single file or list and make it a list
        if type(partslist) == str:
            partslist = [partslist]
        #Build a list of parts CAD objects
        parts = []
        for part in partslist:
#            idx = np.where(np.asarray(self.ROI) == part)[0][0]
            count = 0
            for i in range(len(self.CADobjs)):
                if part == self.CADobjs[i].Label:
                    count += 1
                    self.ROIparts[idx] = self.CADobjs[i]
                    self.ROImeshes[idx] = self.part2mesh(self.ROIparts[idx])[0]

            if count == 0:
                print("Part "+part+" not found in CAD.  Cannot Mesh!")
                log.info("Part "+part+" not found in CAD.  Cannot Mesh!")
        return

    def intersectObjFromPartnum(self, partslist, resolution):
        """
        Generates intersect objects from list of part numbers.

        if resolution is 'standard' then generates mesh using FreeCAD
        Standard algorithm
        """
        #Check if this is a single file or list and make it a list
        if type(partslist) == str:
            partslist = [partslist]
        #Build a list of parts CAD objects
        parts = []
        for part in partslist:
            count = 0
            idx = np.where(np.asarray(self.intersectList) == part)[0][0]
            for i in range(len(self.CADobjs)):
                if part == self.CADobjs[i].Label:
                    count += 1
                    self.intersectParts[idx] = self.CADobjs[i]
                    if resolution=="standard":
                        self.intersectMeshes[idx] = self.part2meshStandard(self.intersectParts[idx])[0]
                    else:
                        self.intersectMeshes[idx] = self.part2mesh(self.intersectParts[idx], resolution=resolution)[0]

            if count == 0:
                print("Part "+part+" not found in CAD.  Cannot Mesh!")
                log.info("Part "+part+" not found in CAD.  Cannot Mesh!")

        return



    def loadSTEP(self):
        """
        Loads CAD STEP (ISO 10303-21) file into object
        """
        self.CAD = Import.open(self.STPfile)
        self.CADdoc = FreeCAD.ActiveDocument

        #Coordinate permutation if necessary
        if self.permute_mask=='True':
            self.permuteSTEP()
            #self.permuteSTEPAssy()
            self.permute_mask = False

        #Save all parts/objects
        self.CADobjs = self.CADdoc.Objects
        self.CADparts = []
        for obj in self.CADobjs:
            if type(obj) == Part.Feature:
                self.CADparts.append(obj)

        print("Loaded STEP file: " + self.STPfile)
        log.info("Loaded STEP file: " + self.STPfile)
        return

    def saveSTEP(self, file, objs):
        """
        Saves CAD STEP (ISO 10303-21) file

        objs    CAD objects.  If you want to preserve Assembly architecture
                then you need to only save [FreeCAD.ActiveDocument.ASSEMBLY].
                Otherwise each part should be in a list
        file    filename to save into
        """
        try:
            ImportGui.export(objs, file)
            print("Saved new STEP file: " + file)
        except:
            print("Error saving STEP file.  Aborting.")
        return


    def permuteSTEPAssy(self):
        """
        Cyclic permutation on STPfile assembly preserving right hand rule
        works on assemblies
        """
        ang = 90.0
        #here we assume we received CAD from NSTXU engineers who have y-axis vertical
        axis = FreeCAD.Vector(1,0,0)
        rot = FreeCAD.Rotation(axis,ang)
        if self.assembly_mask:
            self.CADdoc.ASSEMBLY.Placement = FreeCAD.Placement(axis,rot)
        print("CAD Permutation Complete")
        log.info("CAD Permutation Complete")
        return

    def permuteSTEP(self):
        """
        Cyclic permutation on STPfile assembly preserving right hand rule
        Works on individual parts
        """
        rot = FreeCAD.Placement( FreeCAD.Vector(0,0,0), FreeCAD.Rotation(0,0,90) )
        for obj in FreeCAD.ActiveDocument.Objects:
            if type(obj) == Part.Feature:
                obj.Placement = rot.multiply(obj.Placement)
        print("CAD Permutation Complete")
        log.info("CAD Permutation Complete")
        return

    def getLabels(self, parts):
        """
        Gets labels from a list of CAD parts and returns a list of labels
        """
        #Check if this is a single part or list and make it a list
        if type(parts) != list:
            parts = [parts]
        labels = []
        for part in parts:
            labels.append(part.Label)
        return labels

    def part2mesh(self, part, resolution=None, mode='fine'):
        """
        Converts CAD object to mesh object, and adds mesh object to CAD document
        if part is a list of objects, returns a list of meshes.
        If part isn't a list, freecad throws an error.  Use this to determine
        if part is a list of parts or a single object, and handle each case
        correctly.  Returns a list of mesh objects

        This function uses the FreeCAD Mefisto algorithm, and defines mesh
        by maximum edge length (resolution)
        """
        if resolution == None: resolution = float(self.ROIGridRes)
        else: resolution = float(resolution)
        #Check if this is a single file or list and make it a list
        if type(part) != list:
            part = [part]
        meshes = []
        for i in range(len(part)):
            shape = part[i].Shape.copy(False)
            shape.Placement = part[i].getGlobalPlacement()
            print('Meshing part ' + part[i].Label)
            log.info('Meshing part ' + part[i].Label)
            mesh = MeshPart.meshFromShape(shape, MaxLength=resolution)
            meshes.append(mesh)
        print("Converted parts to mesh objects at resolution: {:f}".format(resolution))
        log.info("Converted parts to mesh objects at resolution: {:f}".format(resolution))
        return meshes

    def part2meshStandard(self, part, surfDev=1.0, angDev=0.523599):
        """
        Converts CAD object to mesh object, and adds mesh object to CAD document
        if part is a list of objects, returns a list of meshes.
        If part isn't a list, freecad throws an error.  Use this to determine
        if part is a list of parts or a single object, and handle each case
        correctly.  Returns a list of mesh objects

        This function uses the FreeCAD Standard algorithm, and defines mesh
        by surface and angular deviation.  Default surface deviation is 0.1mm,
        and default angular deviation is 0.523599rad (30deg)
        """
        #Check if this is a single file or list and make it a list
        if type(part) != list:
            part = [part]
        meshes = []
        for i in range(len(part)):
            shape = part[i].Shape.copy(False)
            shape.Placement = part[i].getGlobalPlacement()
            print('Meshing part ' + part[i].Label)
            log.info('Meshing part ' + part[i].Label)
            mesh = MeshPart.meshFromShape(Shape=shape,
                                          LinearDeflection=surfDev,
                                          AngularDeflection=angDev,
                                          Relative=False)
            meshes.append(mesh)
        print("Converted parts to mesh objects using Standard algorithm.")
        log.info("Converted parts to mesh objects using Standard algorithm.")
        return meshes


    def writeMesh2file(self, mesh, label, path='./', resolution=None):
        """
        Writes a mesh object to STL file named by part number.
        If mesh is a list of mesh objects, then write a separate file for
        each mesh object in the list.  Does not overwrite / clobber
        """
        if resolution == None:
            resolution = self.ROIGridRes

        if type(resolution) != str:
            resolution=str(resolution)
        #Check if this is a single file or list and make it a list
        if type(mesh) != list:
            mesh = [mesh]
        if type(label)!= np.ndarray:
            if type(label) != list:
                label = [label]


        #Recursively make dirs for STLs
        print("making STL directory")
        log.info("making STL directory: "+path)
        os.makedirs(path, exist_ok=True)

        for i in range(len(mesh)):
            # ___ (3 underdashes) is the str we use to separate mesh name from resolution
            # this MATTERS when we read back in a mesh (see self.loadROIMesh and self.loadIntersectMesh)
            filename = path + label[i] + "___" + resolution +"mm.stl"
            print("Writing mesh file: " + filename)
            log.info("Writing mesh file: " + filename)
            if os.path.exists(filename):
                pass
            else:
                mesh[i].write(filename)
        print("\nWrote meshes to file at resolution: " + resolution)
        log.info("\nWrote meshes to file at resolution: " + resolution)
        return

    def loadROIMesh(self, filenames, idx):
        """
        Reads a previously generated STL file and saves object into class.  If
        filename is a list of filenames, then read each into a separate index
        of mesh variable in class.  filename should match a part number from the
        ROI
        """
        #Check if this is a single file or list and make it a list
        if type(filenames) == str:
            filenames = [filenames]
        for file in filenames:
            mesh = Mesh.Mesh(file)
            partnum = file.split('/')[-1].split('___')[0]
#            idx = np.where(np.asarray(self.ROI) == partnum)[0][0]
            #Find CAD object that matches this part number
            for i in range(len(self.CADobjs)):
                if partnum == self.CADobjs[i].Label:
                    self.ROIparts[idx] = self.CADobjs[i]
            self.ROImeshes[idx] = mesh
        print("Loaded STL files")
        log.info("Loaded STL files")
        return

    def loadIntersectMesh(self, filenames):
        """
        Reads a previously generated STL file and saves object into class.  If
        filename is a list of filenames, then read each into a separate index
        of mesh variable in class.
        """
        #Check if this is a single file or list and make it a list
        if type(filenames) == str:
            filenames = [filenames]
        for file in filenames:
            print('Loading ' + file)
            log.info('Loading ' + file)
            mesh = Mesh.Mesh(file)
            partnum = file.split('/')[-1].split('___')[0]
            idx = np.where(np.asarray(self.intersectList) == partnum)[0][0]
            #Find CAD object that matches this part number
            for i in range(len(self.CADobjs)):
                if partnum == self.CADobjs[i].Label:
                    self.intersectParts[idx] = self.CADobjs[i]
            self.intersectMeshes[idx] = mesh
        print("Loaded STL files")
        log.info("Loaded STL files")
        return



    def load1Mesh(self, filename):
        """
        Reads a previously generated STL file and generates 1 mesh object.
        """
        mesh = Mesh.Mesh(filename)
        return mesh

    def normsCentersAreas(self, meshes):
        """
        Gets face normals and face centers.  Both norms and centers are arrays
        of length mesh.CountFacets, consisting of three components (x,y,z) per
        facet
        """
        #Check if this is a single mesh or list and make it a list
        if type(meshes) != list:
            meshes = [meshes]

        norms = []
        centers = []
        areas = []
        for mesh in meshes:
            #mesh = obj.Mesh
            if (mesh == None) or (mesh=='None'):
                print("No Mesh for one of these objects.  Did you have a typo in input file?")
                print("Check HEAT output for Mesh Not Found errors")
                log.info("No Mesh for one of these objects.  Did you have a typo in input file?")
                log.info("Check HEAT output for Mesh Not Found errors")
            N_facets = mesh.CountFacets
            x = np.zeros((N_facets,3))
            y = np.zeros((N_facets,3))
            z = np.zeros((N_facets,3))

            for i,facet in enumerate(mesh.Facets):
                #mesh points
                for j in range(3):
                    x[i][j] = facet.Points[j][0]
                    y[i][j] = facet.Points[j][1]
                    z[i][j] = facet.Points[j][2]

            # scale and permute if necessary
            x,y,z = self.scale_and_permute(x,y,z)
            # get face normals and face centers
            norms.append(self.faceNormals(mesh))
            centers.append(self.faceCenters(x,y,z))
            areas.append(self.faceAreas(mesh))
        return norms,centers,areas


    def scale_and_permute(self, x_old, y_old, z_old, permute_mask=False, unitConvert=1.0):
        """
        Scales input mesh vectors if necessary (ie for unit conversion)
        Performs coordinate permutation on input mesh vectors if necessary
        (ie if CAD had y-axis as vertical axis)
        """
        if hasattr(self,'permute_mask'):
            permute_mask = self.permute_mask
        if hasattr(self,'unitConvert'):
            unitConvert = self.unitConvert

        #First handle coordinate permutations (preserve right hand rule)
        if permute_mask==True:
            x = z_old
            y = x_old
            z = y_old
        else:
            x = x_old
            y = y_old
            z = z_old
        #Scale inputs to get units in meters
        x *= float(unitConvert)
        y *= float(unitConvert)
        z *= float(unitConvert)
        return x, y, z

    def faceNormals(self, mesh):
        """
        returns normal vectors for single freecad mesh object in cartesian
        coordinates
        """
        #face normals
        normals = []
        for i, facet in enumerate(mesh.Facets):
            vec = np.zeros((3))
            for j in range(3):
                vec[j] = facet.Normal[j]
            normals.append(vec)
        return np.asarray(normals)

    def faceAreas(self, mesh):
        """
        returns face areas for mesh element
        """
        #face area
        areas = []
        for i, facet in enumerate(mesh.Facets):
            areas.append(facet.Area)
        return np.asarray(areas)


    def faceCenters(self, x, y, z):
        """
        returns centers of freecad mesh triangle in cartesian coordinates
        """
        #face centers
        centers = np.zeros((len(x), 3))
        centers[:,0] = np.sum(x,axis=1)/3.0
        centers[:,1] = np.sum(y,axis=1)/3.0
        centers[:,2] = np.sum(z,axis=1)/3.0
        return centers

    def stp2stl(self, resolution = None):
        """
        Reads in an STEP file (ISO 10303-21) and outputs an STL triangular mesh
        with maximum edge length defined by self.gridRes
        """
        if resolution == None: resolution = self.ROIGridRes

        t0 = time.time()
        shape = Part.Shape()
        shape.read(self.STPfile)
        print("CAD STEP read took {:f} seconds".format(time.time() - t0))
        mesh_shp = MeshPart.meshFromShape(shape, MaxLength=resolution)
        print("Part mesh took {:f} seconds".format(time.time() - t0))
        mesh_shp.write(self.STLfile)
        print("STP => STL conversion completed in {:f} seconds".format(time.time() - t0))

        return

    def getCOMs(self, parts):
        """
        Finds center of masses for a list of parts
        """
        #Check if this is a single part or list and make it a list
        if type(parts) != list:
            parts = [parts]

        #Find center of masses (COMs)
        COMs = []
        for part in parts:
            com = part.Shape.CenterOfMass
            COMs.append(com)

        COMs = np.asarray(COMs)
        return COMs

    def findRelativePhi(self, sourceParts, targetParts):
        """
        Finds relative phi between source part and all other parts in CADparts.
        Phi for each part is calculated at the center of mass (COM).  deltaPhis
        is an array with sources in one dimension and targets in other
        dimension, describing the relative angle between sources and targets.
        Source parts are all parts in the ROI.  Target parts are all the parts
        in the input STP file.
        """
        #Check if this is a single part or list and make it a list
        if type(sourceParts) != list:
            sourceParts = [sourceParts]
        if type(targetParts) != list:
            targetParts = [targetParts]

        deltaPhis = np.zeros((len(sourceParts),len(targetParts)))
        for i,source in enumerate(sourceParts):
            sourceCOM = self.getCOMs(source)[0]
            r,z,phi = tools.xyz2cyl(sourceCOM[0],sourceCOM[1],sourceCOM[2])
            for j,target in enumerate(targetParts):
                targetCOM = self.getCOMs(target)[0]
                r_t,z_t,phi_t = tools.xyz2cyl(targetCOM[0],targetCOM[1],targetCOM[2])
#                #Handle cases where we have negative angles
                if phi_t < 0:
                    phi_t = phi_t + 2*np.pi
                deltaPhis[i,j] = phi_t - phi
                if deltaPhis[i,j] > np.pi:
                    deltaPhis[i,j] = deltaPhis[i,j] - 2*np.pi

        return deltaPhis



    def findPotentialIntersectParts(self,deltaPhis,sourceBts,sourceParts,targetParts):
        """
        Uses center of masses (COMs) of the parts in the ROI and the toroidal
        field at each COM to determine which other parts need to be checked
        for intersections.
        intersect_mask is a bitmask (binary) array with sources in one dimension
        and targets in other dimension.  Value is 1 if target is potential
        intersection part and 0 if target is not potential intersection part.
        The toriodal component of magnetic field at source part is the
        parameter used to determine if an intersection is possible.
        """
        #Check if this is a single part or list and make it a list
        if type(sourceParts) != list:
            sourceParts = [sourceParts]
        if type(targetParts) != list:
            targetParts = [targetParts]

        intersect_mask = np.zeros((len(sourceParts),len(targetParts)))

        #Only use tiles who are upstream in toroidal field direction
        for i,source in enumerate(sourceParts):
            if sourceBts[i] > 0:
                intersect_mask[np.where(deltaPhis <= 0)] = 1
            else:
                intersect_mask[np.where(deltaPhis >= 0)] = 1

        #Now eliminate tiles that are over 0.25m from source COM in any direction
        sourceCOMs = self.getCOMs(sourceParts)
        targetCOMs = self.getCOMs(targetParts)
        for i in range(len(sourceCOMs)):
            for j in range(len(targetCOMs)):
                delta = np.abs(np.subtract(sourceCOMs[i],targetCOMs[j]))
                if (delta[0] > 500) or (delta[1] > 500) or (delta[2] > 500):
                     intersect_mask[i,j] = 0
        return intersect_mask

    def meshPotentialIntersects(self, intersect_mask, sourceParts, targetParts,
                                resolution):
        """
        Creates mesh objects for the parts set to 1 in intersect_mask.
        First checks to see if mesh at desired resolution exists.  If mesh
        already exists, loads it, otherwise creates new mesh object.
        """
        #Check if this is a single part or list and make it a list
        if type(sourceParts) != list:
            sourceParts = [sourceParts]
        if type(targetParts) != list:
            targetParts = [targetParts]

        #Make a list of all the target part numbers we need to check
        idx = np.where(intersect_mask == 1)
        parts = []
        labels = []
        for j,target in enumerate(targetParts):
            if intersect_mask[j] == 1:
                parts.append(target)
                labels.append(target.Label)

        parts = list(parts)
        labels = list(labels)
        meshes = self.part2mesh(parts,resolution)
        return meshes, labels

    def write_normal_pointcloud(self,centers,norms,dataPath, tag=None):
        """
        In paraview use TableToPoints => Calculator => Glyph
        Calculator should have this formula:
        (iHat*Nx) + (jHat*Ny) + (kHat*Nz)
        """
        print("Creating Normal Point Cloud")
        log.info("Creating Normal Point Cloud")
        if tag is None:
            pcfile = dataPath + 'NormPointCloud.csv'
        else:
            pcfile = dataPath + 'NormPointCloud_' +tag+ '.csv'

        pc = np.zeros((len(centers), 6))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = norms[:,0]
        pc[:,4] = norms[:,1]
        pc[:,5] = norms[:,2]
        head = """X,Y,Z,Nx,Ny,Nz"""
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)
        #np.savetxt('points.asc', pc[:,:-1], delimiter=' ',fmt='%.10f')
        #print("Wrote point cloud file: " + pcfile)

        #Now save a vtk file for paraviewweb
        tools.createVTKOutput(pcfile, 'glyph', 'Norm_pointcloud')
        return

    def stripSTPfile(self,partfile,rawSTP,outfile=None,partsOnly=True):
        """
        Strips an stpfile down to the parts defined in the parts csv input file
        STPfile is the output file defined in input file to HEAT, and is included
        in self variable.  rawSTP is the input CAD file that we want to strip.

        If partsOnly is True then we only copy parts (not assemblies), which
        technically means only objects with type=Part.Feature
        """
        print('Stripping STP file')
        t0 = time.time()
        with open(partfile) as f:
            part_list = f.read().splitlines()
        print("Read parts list...")
        print(part_list)

        #Input / Output STEP files
        infile = rawSTP
        if outfile is None:
            outfile = self.STPfile

        #If a shape has a label in part_list, keep it
        CAD = Import.open(infile)
        newobj = []
        count = 0
        for i,obj in enumerate(FreeCAD.ActiveDocument.Objects):
            if any(substring in obj.Label for substring in part_list):
                #conditional to check if item is part (Part.Feature) or
                # assembly (App.Part).
                # This could be adapted in future to exclude objects containing
                # specific string (like "_ASM") in Label that CAD engineer uses
                # for assemblies
                if partsOnly==True:
                    if type(obj)==Part.Feature:
                        count+=1
                        newobj.append(obj)
                        newobj[-1].Placement = obj.getGlobalPlacement()
                else:
                    count+=1
                    newobj.append(obj)
                    newobj[-1].Placement = obj.getGlobalPlacement()

        #Export to a new step file
        Import.export(newobj, outfile)
        print("Step file export complete.")
        print("Exported {:d} part objects".format(count))
        print("Execution took {:f} seconds".format(time.time() - t0))
        return

    def readIntersects(self, infile, targetParts):
        """
        read intersection parts by partname
        returns bitmask where 1 means that part is a potential intersection and
        0 means there is not potential intersection on that target
        """
        data = pd.read_csv(infile, sep=',', comment='#', names=['Part'], skipinitialspace=True)
        intersect = np.zeros((len(targetParts)))
        for j,target in enumerate(targetParts):
            for name in data['Part']:
                if target.Label == str(name):
                    intersect[j] = 1
        return intersect

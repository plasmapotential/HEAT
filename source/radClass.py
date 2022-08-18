#RADClass.py
#Description:   Base HEAT Radiated Power Module
#Engineer:      T Looby
#Date:          20220613
import os
import sys
import numpy as np
import pandas as pd
import logging
import multiprocessing
import time
import open3d as o3d
from scipy.spatial import ConvexHull

import toolsClass
log = logging.getLogger(__name__)
tools = toolsClass.tools()

class RAD:
    """
    HEAT Radiated Power Class.  For calculating photon heat loads

    """
    def __init__(self, rootDir=None, dataPath=None, chmod=0o774, UID=-1, GID=-1):
        """
        dataPath is the location where we write all output to
        """
        self.rootDir = rootDir
        tools.rootDir = self.rootDir
        self.dataPath = dataPath
        tools.dataPath = self.dataPath
        self.chmod = chmod
        self.GID = GID
        self.UID = UID
        return

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization

        Here is a list of variables with description:
        """


        self.allowed_vars = [
                             'radFile',
                             'Ntor',
                             'Nref',
                             'phiMin',
                             'phiMax',
                            ]
        return

    def setTypes(self):
        """
        set input variable types
        """
        self.Ntor = int(self.Ntor)
        self.Nref = int(self.Nref)
        self.phiMin = float(self.phiMin)
        self.phiMax = float(self.phiMax)
        return

    def read2DSourceFile(self, file):
        """
        reads a comma delimited file that describes the radiated power on an
        R,Z grid. This point cloud (pc) should represent the axisymmetric
        plasma radiated power across a single poloidal plane (ie sliced at phi).
        The power defined at each RZ point is the radiated photon power.

        It is assumed that the RZ points are located at the centroids of mesh
        elements
        """
        print("Reading 2D photon radiation source file: "+file)
        log.info("Reading 2D photon radiation source file: "+file)
        self.PC2D = pd.read_csv(file, header=0, names=['R','Z','P']).values #convert to m
        #self.PC2D /= 1000.0
        return

    def paraview2D(self):
        """
        creates a 2D .vtk representation of self.PC2D
        """
        return

    def paraview3D(self):
        """
        creates a 3D .vtk representation of self.PC3D
        """
        return

    def create3DFrom2D(self):
        """
        creates a 3D point cloud from a 2D point cloud.
        assumes the 2D point cloud is in a poloidal plane (ie sliced at a phi)
        and extrudes the 2D point cloud around the machine in toroidal direction.
        This extrusion results in 3D voxels (3D version of a 2D pixel), where
        each point is the centroid of the voxel.  User specifies resolution
        of voxels in toroidal direction.

        shape of output array is NR X 4 X Nphi
        """
        print('Building 3D radiated power profile')
        log.info('Building 3D radiated power profile')

        #build out source array (R,Z,phi)
        NR = self.PC2D.shape[0]
        self.Nphi = len(self.phis)
        PC2D = np.append( self.PC2D, np.zeros([NR,1]), 1 ) #add column for phi
        self.PC3D = np.repeat(PC2D[:,:,np.newaxis],self.Nphi,axis=2) #repeat Nphi times
        self.PC3D[:,-1,:] = np.repeat(self.phis[np.newaxis], NR, axis=0) #fill in phis
        #normalize Power to the number of toroidal steps
        self.PC3D[:,2,:] /= self.Nphi

        #RZphi = np.vstack([self.PC3D[:,0,:].flatten(),self.PC3D[:,1,:].flatten()]).T
        #RZphi = np.vstack([RZ.T, self.PC3D[:,3,:].flatten()]).T
        self.NradPts = len(self.PC3D[:,0,:].flatten())
        return

    def getPhis(self, Ntor, phiMin=0.0, phiMax=360.0):
        """
        Ntor:   number of discrete voxels in toroidal direction (360 would result
                in voxels of 1 degree toroidal extent)

        phiMin: minimum angle of phi in degrees
        phiMax: maximum angle of phi in degrees
        """
        self.phis = np.linspace(phiMin, phiMax, Ntor)
        return



    def preparePowerTransfer(self, PFC, CAD, mode=None):
        """
        builds the meshes and point clouds necessary for calculation.

        prepares the radiation sources and targets (ROI PFCs) for the
        radiated power calc.  Gets the source centers, target centers, and
        intersection mesh vertices.
        """
        #build out source centers
        x,y,z = tools.cyl2xyz(self.PC3D[:,0,:].flatten(),
                                     self.PC3D[:,1,:].flatten(),
                                     self.PC3D[:,3,:].flatten() )
        self.sources = np.vstack([x,y,z]).T
        #power on source voxels
        self.sourcePower = self.PC3D[:,2,:].flatten()

        #build out target array (ROI centers for this PFC)
        self.targetCtrs = PFC.centers
        self.targetNorms = PFC.norms
        self.targetAreas = PFC.areas
        targetMeshes = []

        #build out intersect array
        numTargetFaces = 0
        intersectPoints = []
        intersectNorms = []
        totalMeshCounter = 0

        #add ROI mesh
        for i,target in enumerate(CAD.ROImeshes):
            if CAD.ROIList[i] == PFC.name:
                targetMeshes.append(target)
                totalMeshCounter+=target.CountFacets
                numTargetFaces += target.CountFacets
                #append target data
                for face in target.Facets:
                    intersectPoints.append(face.Points)
                    intersectNorms.append(face.Normal)


        for i,intersect in enumerate(CAD.intersectMeshes):
            try:
                totalMeshCounter+=intersect.CountFacets
            except:
                print("Cannot count faces because "+CAD.intersectList[i]+" is not a mesh!")
                continue
            #check if this target is a potential intersection
            if CAD.intersectList[i] in PFC.intersects:
                #exclude self shadowing
                if CAD.intersectList[i] == PFC.name:
                    pass
                else:
                    targetMeshes.append(intersect)
                    numTargetFaces += intersect.CountFacets
                    #append target data
                    for face in intersect.Facets:
                        intersectPoints.append(face.Points)
                        intersectNorms.append(face.Normal)



        self.intersectPoints = np.asarray(intersectPoints)/1000.0 #scale to m
        intersectNorms = np.asarray(intersectNorms)

        self.p1 = self.intersectPoints[:,0,:]  #point 1 of mesh triangle
        self.p2 = self.intersectPoints[:,1,:]  #point 2 of mesh triangle
        self.p3 = self.intersectPoints[:,2,:]  #point 3 of mesh triangle

        self.Ni = len(self.sources)
        self.Nj = len(self.targetCtrs)

        #build single mesh for KDtree
        if mode=='kdtree':
            self.combinedMesh = CAD.createEmptyMesh()
            for m in targetMeshes:
                self.combinedMesh.addFacets(m.Facets)

        #build objects for open3D ray tracing
        elif mode=='open3d':
            combinedMesh = CAD.createEmptyMesh()
            for m in targetMeshes:
                combinedMesh.addFacets(m.Facets)
            oldMask = CAD.overWriteMask
            CAD.overWriteMask = True
            CAD.writeMesh2file(combinedMesh, 'combinedMesh', path=CAD.STLpath, resolution='standard')
            CAD.overWriteMask = oldMask
            self.meshFile = CAD.STLpath + 'combinedMesh' + "___standard.stl"
            #self.meshFile = '/home/tom/source/dummyOutput/SOLID843___5.00mm.stl'

        print("\nTotal Rad Intersection Faces: {:d}".format(totalMeshCounter))
        print("Rad Intersect Faces for this PFC: {:d}".format(numTargetFaces))
        print("# of radiated power point sources: {:d}".format(self.Ni))
        print("# of PFC ROI mesh elements: {:d}".format(self.Nj))
        print("Total number of source-target ray-tracing calculations: {:d}\n".format(self.Ni*self.Nj))

        return

    def calculatePowerTransferOpen3D(self, mode=None):
        """
        Maps power between sources and targets (ROI PFCs).  Uses Open3D to
        perform ray tracing.  Open3D can be optimized for CPU or GPU.

        Uses Open3D to accelerate the calculation:
        Zhou, Qian-Yi, Jaesik Park, and Vladlen Koltun. "Open3D: A modern
        library for 3D data processing." arXiv preprint arXiv:1801.09847 (2018).
        """
        t0 = time.time()
        powerFrac = np.zeros((self.Ni,self.Nj))
        Psum = np.zeros((self.Nj))
        self.hullPower = np.zeros((self.Ni))
        for i in range(self.Ni):
            r_ij = np.zeros((self.Nj,3))
            r_ij = self.targetCtrs - self.sources[i]
            #r_ij *= 1000.0
            rMag = np.linalg.norm(r_ij, axis=1)
            rNorm = r_ij / rMag.reshape((-1,1))
            rdotn = np.sum(rNorm*self.targetNorms, axis=1)
            q1 = np.tile(self.sources[i]*1000.0,(self.Nj,1))

            #build mesh and tensors for open3d
            mesh = o3d.io.read_triangle_mesh(self.meshFile)
            mesh.compute_vertex_normals()
            mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
            scene = o3d.t.geometry.RaycastingScene()
            mesh_id = scene.add_triangles(mesh)

            #calculate intersections
            rays = o3d.core.Tensor([np.hstack([q1,rNorm])],dtype=o3d.core.Dtype.Float32)
            hits = scene.cast_rays(rays)

            #convert open3d CPU tensors back to numpy
            hitMap = hits['primitive_ids'][0].numpy()
            distMap = hits['t_hit'][0].numpy()


            for j in range(self.Nj):
                powerFrac[i,j] = np.abs(rdotn[j])*self.targetAreas[j]/(4*np.pi*rMag[j]**2)
                #assign power
                #if hitMap[i,j] == j and distMap[i,j] >= rMag[i,j]:
                if hitMap[j] == j:
                    Psum[j] += self.sourcePower[i]*powerFrac[i,j]
                else:
                    Psum[j] += 0.0

                #for testing
                idxTest=None
                if j==idxTest:
                    print('\n===')
                    print(i)
                    print(self.targetCtrs[j])
                    print(rNorm[j]*rMag[j])
                    print(hitMap[j])
                    print(self.sourcePower[i])
                    print(powerFrac[i,j])

            #compute convex hull on unit sphere around point i and calculate
            #power balance via ratio of hull area to sphere area
            #calculate spherical coordinates for each target point, with i as (0,0,0)
            #note that we perform a coordinate permutation on rNorm, to prevent
            #angle wrap cases: (x,y,z) => (z,x,y)
            #If you are getting abnormally convex hull power, try a different permutation
            #first transform to 2D (phi,theta) => (x,y) using Jacobians
            theta = np.arccos( (rNorm[:,0]) ) #z after permutation
            phi = np.arctan2( (rNorm[:,2]), (rNorm[:,1])  ) # y,x after permutation
            #map (phi,theta) to cartesian plane using Jacobian
            points = np.vstack([phi, -np.cos(theta)]).T
            #calculate the convex hull in 2D
            hull = ConvexHull(points)
            #hull area for 2D points is volume, scale power to full solid angle
            self.hullPower[i] = hull.volume / (4*np.pi) * self.sourcePower[i]

        self.pdotn = rdotn
        self.powerFrac = powerFrac
        self.targetPower = Psum

        print("Photon tracing took {:f} seconds \n".format(time.time()-t0))
        log.info("Photon tracing took {:f} seconds \n".format(time.time()-t0))
        return

    def calculatePowerTransfer(self, mode=None):
        """
        Maps power between sources and targets (ROI PFCs).  Uses CPU
        multiprocessing without acceleration structures
        """
        self.pdotn = np.zeros((self.Ni,self.Nj))
        self.powerFrac = np.zeros((self.Ni,self.Nj))
        self.targetPower = np.zeros((self.Nj))

        #Prepare multiple cores for mesh builder
        Ncores = multiprocessing.cpu_count() - 2 #reserve 2 cores for overhead
        #in case we run on single core machine
        if Ncores <= 0:
            Ncores = 1
        print('Initializing intersection check across {:d} cores'.format(Ncores))
        print('Spawning tasks to multiprocessing workers')
        try:
            #manager can be used for locking, but shouldnt be necessary so long
            #as we dont overlap write locations between workers
            pool = multiprocessing.Pool(Ncores)
            if mode=='kdtree':
                tools.combinedMesh = self.combinedMesh
                output = np.asarray(pool.map(self.powerFracMapParallelKDtree, np.arange(self.Nj)))
            else:
                output = np.asarray(pool.map(self.powerFracMapParallelNoAcc, np.arange(self.Nj)))
        finally:
            pool.close()
            pool.join()
            del pool


        print("Multiprocessing complete")
        for j in range(self.Nj):
            self.pdotn[:,j] = output[j,0]
            self.powerFrac[:,j] = output[j,1]
            self.targetPower[j] = output[j,2]

        return

    def powerFracMapParallelNoAcc(self, j):
        """
        Calculates the fractional contribution of power from each self.sources
        point to the target mesh element (ROI PFC), j.  Includes intersection
        checking with the intersect mesh along the way.

        This function is meant to be called in parallel, where j represents
        the target (ROI PFC) mesh element we are calculating the power
        transfer for.  This function does not calculate the explicit power, but
        rather it builds the self.powerFrac (i,j) matrix, which can be scaled
        to any input power.

        does not include acceleration structures, so this is a brute force method
        """
        #source to target vectors
        p_ij = self.targetCtrs[j] - self.sources
        pMag = np.linalg.norm(p_ij, axis=1)
        pNorm = p_ij / pMag[:, np.newaxis]
        pdotn = np.sum(pNorm*self.targetNorms[j], axis=1)

        #backface culling
        shadowMask = np.zeros((self.Ni))
        backFaces = np.where( pdotn > 0 )[0]
        shadowMask[backFaces] = 1.0
        use0 = np.where(shadowMask == 0)[0]
        #use0 =np.arange(self.Ni)
        #should we add backface culling for intersects here?!

        #run in matrix or loop mode depending upon available memory?
        #availableRAM = psutil.virtual_memory().available #bytes
        #neededRAM = len(self.sources) * len(self.targets) * len(self.intersects) * 112 #bytes

        q1 = self.sources[use0,:]
        q2 = self.targetCtrs[j]
        Nt = len(self.p1)

        powerFrac = np.zeros((self.Ni))

        #get power fraction multiprocessing Array object and shape to be Ni X Nj
        #powerFrac = np.frombuffer(self.powerFrac.get_obj()).reshape((self.Ni, self.Nj))
        #powerFrac = np.array(self.powerFrac[:]).reshape(self.Ni, self.Nj)
        #loop thru sources checking if shadowed for this target
        Psum = 0.0
        for i in range(len(q1)):
            #Perform Intersection Test
            q13D = np.repeat(q1[i,np.newaxis], Nt, axis=0)
            q23D = np.repeat(q2[np.newaxis], Nt, axis=0)
            sign1 = np.sign(tools.signedVolume2(q13D,self.p1,self.p2,self.p3))
            sign2 = np.sign(tools.signedVolume2(q23D,self.p1,self.p2,self.p3))
            sign3 = np.sign(tools.signedVolume2(q13D,q23D,self.p1,self.p2))
            sign4 = np.sign(tools.signedVolume2(q13D,q23D,self.p2,self.p3))
            sign5 = np.sign(tools.signedVolume2(q13D,q23D,self.p3,self.p1))
            test1 = (sign1 != sign2)
            test2 = np.logical_and(sign3==sign4,sign3==sign5)

            #we intersected something, so face is shadowed
            if np.sum(np.logical_and(test1,test2)) > 0:
                powerFrac[i] = 0.0
            #we didn't intersect. assign fraction of power
            else:
                powerFrac[i] = np.abs(pdotn[i])*self.targetAreas[j]/(4*np.pi*pMag[i]**2)
                Psum += self.sourcePower[i]*powerFrac[i]

        return pdotn, powerFrac, Psum

    def calculateBRDF(self, material):
        """
        constructs the Bidirectional Reflection Distribution Function for a user
        defined material.  Once constructed, user can evaluate specific BRDF
        values for input and output angles.
        """
        return

    def savePowerFrac(self,PFC):
        """
        saves power fracs into a file
        """
        f =PFC.controlfilePath + 'radFracs_' + PFC.name + '.csv'
        np.savetxt(f, PFC.radPowerFracs, delimiter=',',fmt='%.10f')
        return

    def write_Prad_pointcloud(self,centers,Prad,dataPath,tag=None):
        print("Creating Heat Flux Point Cloud")
        log.info("Creating Heat Flux Point Cloud")
        prefix = 'P_radSource'
        if tag is None:
            pcfile = dataPath + prefix + '.csv'
        else:
            pcfile = dataPath + prefix + '_'+tag+'.csv'
        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = Prad
        head = "X,Y,Z,Prad"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

        #Now save a vtk file for paraviewweb
        if tag is None:
            tools.createVTKOutput(pcfile, 'points', prefix)
        else:
            name = prefix+'_'+tag
            tools.createVTKOutput(pcfile, 'points', name)
        return

    def writeRadFileData(self,radFile,radData,tmpDir):
        """
        writes data passed in string object (from GUI) to files in
        tmpDir directory for use later on in HEAT

        the self.tmpDir directory is accessible to the GUI users for uploading
        and downloading

        this function is called from GUI because objects are json / base64
        """
        import base64
        data = radData.encode("utf8").split(b";base64,")[1]
        path = tmpDir + radFile
        print("Writing local radFile: "+path)
        log.info("Writing local radFile: "+path)
        with open(path, 'wb') as f:
            f.write(base64.decodebytes(data))

        return path

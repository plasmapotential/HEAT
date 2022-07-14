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
        Nphi = len(self.phis)
        PC2D = np.append( self.PC2D, np.zeros([NR,1]), 1 ) #add column for phi
        self.PC3D = np.repeat(PC2D[:,:,np.newaxis],Nphi,axis=2) #repeat Nphi times
        self.PC3D[:,-1,:] = np.repeat(self.phis[np.newaxis], NR, axis=0) #fill in phis
        #RZphi = np.vstack([self.PC3D[:,0,:].flatten(),self.PC3D[:,1,:].flatten()]).T
        #RZphi = np.vstack([RZ.T, self.PC3D[:,3,:].flatten()]).T
        self.NradPts = len(self.RAD.PC3D[:,0,:].flatten())
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

    def preparePowerTransfer(self, PFC, CAD):
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

        #build out intersect array
        numTargetFaces = 0
        intersectPoints = []
        intersectNorms = []
        totalMeshCounter = 0
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

        print("\nTotal Rad Intersection Faces: {:d}".format(totalMeshCounter))
        print("Rad Intersect Faces for this PFC: {:d}".format(numTargetFaces))
        print("# of radiated power point sources: {:d}".format(self.Ni))
        print("# of PFC ROI mesh elements: {:d}".format(self.Nj))
        print("Total number of source-target ray-tracing calculations: {:d}\n".format(self.Ni*self.Nj))

        return

    def calculatePowerTransfer(self):
        """
        Maps power between sources and targets (ROI PFCs).  Uses CPU
        multiprocessing.
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
        #Do this try clause to kill any zombie threads that don't terminate
#        try:
#            #manager can be used for locking, but shouldnt be necessary so long
#            #as we dont overlap write locations between workers
#            manager = multiprocessing.Manager()
#            lock=manager.lock()
#            self.powerFrac = manager.Array('f',np.zeros((self.Ni*self.Nj)))
#            self.targetPower = manager.Array('f',np.zeros((self.Nj)))
#            #self.powerFrac = multiprocessing.Array('f',self.Ni*self.Nj)
#            #self.targetPower = multiprocessing.Array('f',self.Nj)
#            pool = multiprocessing.Pool(Ncores)
#            pool.map(self.powerFracMapParallel, np.arange(self.Nj))
#        finally:
#            pool.close()
#            pool.join()
#            del pool
#            del manager

        #Do this try clause to kill any zombie threads that don't terminate
        try:
            #manager can be used for locking, but shouldnt be necessary so long
            #as we dont overlap write locations between workers
            pool = multiprocessing.Pool(Ncores)
            output = np.asarray(pool.map(self.powerFracMapParallel, np.arange(self.Nj)))
        finally:
            pool.close()
            pool.join()
            del pool


        print("Multiprocessing complete")
        for j in range(self.Nj):
            self.pdotn[:,j] = output[j,0]
            self.powerFrac[:,j] = output[j,1]
            self.targetPower[j] = output[j,2]
        print(self.pdotn[:,0])
        print(self.powerFrac[:,0])
        print(self.targetPower[0])
        return

    def powerFracMapParallel(self, j):
        """
        Calculates the fractional contribution of power from each self.sources
        point to the target mesh element (ROI PFC), j.  Includes intersection
        checking with the intersect mesh along the way.

        This function is meant to be called in parallel, where j represents
        the target (ROI PFC) mesh element we are calculating the power
        transfer for.  This function does not calculate the explicit power, but
        rather it builds the self.powerFrac (i,j) matrix, which can be scaled
        to any input power.
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

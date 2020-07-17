#pfcClass.py
#Description:   PFC Module
#Engineer:      T Looby
#Date:          20200616
import numpy as np
import toolsClass
import heatfluxClass
import os
import sys
import time
import pandas as pd
import multiprocessing
tools = toolsClass.tools()

import logging
log = logging.getLogger(__name__)


class PFC:
    """
    this class is where we connect CAD to MHD to HF objects.
    Basically, we create a PFC object for each tile that contains all
    the CAD information regarding the geometry (face centers, face normals, etc.),
    the MHD information (timesteps, equilibrium objects, mapDirection, etc.)
    and the HF information (profile, lq, Psol, etc.)
    This class also builds out the directory tree for each simulation
    """
    def __init__(self, timestepMapRow):
        self.timeStr = timestepMapRow[0]
        self.name = timestepMapRow[1]
        self.mapDirection = timestepMapRow[2]
        self.mapDirectionStruct = self.mapDirection
        self.ionDirection = self.mapDirectionStruct
        self.timeLimits = np.asarray( self.timeStr.split(':') ).astype(int)
        deltat = self.timeLimits[1]-self.timeLimits[0]
        self.timesteps = np.linspace(self.timeLimits[0],self.timeLimits[1],deltat+1,dtype=int)
        self.intersects = timestepMapRow[3].split(':')

        return

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization

        Here is a list of variables with description:
        testvar         dummy for testing

        """


        self.allowed_vars = [
                            'name',
                            'mesh',
                            'centers',
                            'areas',
                            'norms',
                            'Nfaces',
                            'shadowed_mask',
                            't',
                            'ep'
                            ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        return


    def makePFC(self,MHD,CAD,ROIidx,mapDirection,clobberFlag=True):
        """
        Inputs are CAD object, MHD object (with all ep timesteps),
        and ROIidx, which is the index in the ROI that this pfc object
        corresponds to.

        Each PFC object has within it:
        timesteps
        ep objects (equilibrium) indexed to correspond to timesteps
        centers
        normals
        areas
        meshes
        shadowMasks
        """
        #get CAD objects corresponding to this PFC
        self.centers = CAD.ROIctrs[ROIidx] / (1000.0) #convert to meters
        self.norms = CAD.ROInorms[ROIidx]
        self.areas = CAD.ROIareas[ROIidx] / (1000.0**2) #convert to meters^2
        self.mesh = CAD.ROImeshes[ROIidx]
        self.Nfaces = self.mesh.CountFacets
        self.qDiv = np.zeros((len(self.centers)))

        #get MHD objects corresponding to those timesteps
        self.EPs = [] #containers for multiple ep
        self.shadowMasks = [] #container for multiple shadowed_mask
        self.powerSum = [] #container for multiple power summations
        self.tIndexes = [] #container for mapping between PC timesteps and MHD timesteps
        for t in self.timesteps:
            idx = np.where(t==MHD.timesteps)[0]
            if len(idx) > 0:
                self.EPs.append(MHD.ep[idx[0]])
                self.shadowMasks.append(self.backfaceCulling(self.centers,self.norms,MHD,MHD.ep[idx[0]],mapDirection))
                self.powerSum.append(0.0)
                self.tIndexes.append(idx[0])

        #set up file directory structure using first timestep as dummy
        # (it is changed later in GUIclass loop)
        self.t = MHD.timesteps[self.tIndexes[0]]
        self.controlfile = '_lamCTL.dat'
        self.controlfileStruct = '_struct_CTL.dat'
        self.controlfilePath = MHD.dataPath + '/' + '{:06d}/'.format(self.t) + self.name + '/'
        self.gridfile = MHD.dataPath + '/' + '{:06d}/'.format(self.t) + self.name + '/grid.dat'
        self.gridfileStruct = MHD.dataPath + '/' + '{:06d}/'.format(self.t) + self.name + '/struct_grid.dat'
        self.outputFile = self.controlfilePath + 'lam.dat'
        self.structOutfile = self.controlfilePath + 'struct.dat'

        return

    def resetPFCeps(self, MHD):
        """
        resets the PFC objects ep objects using MHD objects eps
        """
        self.EPs = [] #containers for multiple ep gets reset here
        for t in self.timesteps:
            idx = np.where(t==MHD.timesteps)[0]
            if len(idx) > 0:
                self.EPs.append(MHD.ep[idx[0]])
        return


    def backfaceCulling(self,centers,norms,MHD,ep,ionDir):
        """
        Use backface culling technique from computer graphics to eliminate
        faces that are on the 'back side' of the PFC tiles.  Mathematically,
        if: B_hat dot n_hat >= 0, then the face is shadowed.
        returns 1 if face is shadowed, 0 if face is not shadowed
        """
        xyz = centers
        r,z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        BNorms = MHD.Bfield_pointcloud(ep, r, z, phi, ionDir, normal=True)
        dot = np.multiply(norms, BNorms).sum(1)
        shadowed_mask = np.zeros((len(centers)))
        shadowed_mask[np.where(dot >= 0.0)] = 1
        return shadowed_mask

    def write_shadow_pointcloud(self,centers,scalar,dataPath,tag=None):
        print("Creating Shadow Point Cloud")
        log.info("Creating Shadow Point Cloud")

        if tag is None:
            pcfile = dataPath + 'ShadowPointCloud.csv'
        else:
            pcfile = dataPath + 'ShadowPointCloud_'+tag+'.csv'

        #print("Shadow point cloud filename: "+pcfile)
        #log.info("Shadow point cloud filename: "+pcfile)

        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = scalar
        head = "X,Y,Z,ShadowMask"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

        #Now save a vtk file for paraviewweb
        if tag is None:
            tools.createVTKOutput(pcfile, 'points', 'ShadowMask')
        else:
            name = 'ShadowMask_'+tag
            tools.createVTKOutput(pcfile, 'points', name)
        return

    def write_bdotn_pointcloud(self,centers,bdotn,dataPath,tag=None):
        print("Creating bdotn Point Cloud")
        log.info("Creating bdotn Point Cloud")
        if tag is None:
            pcfile = dataPath + 'bdotnPointCloud.csv'
        else:
            pcfile = dataPath + 'bdotnPointCloud_'+tag+'.csv'
        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = bdotn
        head = "X,Y,Z,bdotn"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

        #Now save a vtk file for paraviewweb
        if tag is None:
            tools.createVTKOutput(pcfile, 'points', 'bdotn')
        else:
            name = 'bdotn_'+tag
            tools.createVTKOutput(pcfile, 'points', name)
        return

    def findShadows_structure(self,MHD,CAD,verbose=False, shadowMaskClouds=False):
        """
        Find shadowed faces for a given heatFluxPart object using MAFOT structure.
        Traces field lines from PFC surface, looking for intersections with
        triangles from mesh
        """
        use = np.where(self.shadowed_mask == 0)[0]
        numTargetFaces = 0
        targetPoints = []
        targetNorms = []

        print("\nFinding intersections for {:d} faces".format(len(self.centers[use])))
        log.info("\nFinding intersections for {:d} faces".format(len(self.centers[use])))
        print('Number of target parts: {:f}'.format(len(self.intersects)))
        log.info('Number of target parts: {:f}'.format(len(self.intersects)))

        totalMeshCounter = 0
        for i,target in enumerate(CAD.intersectMeshes):
            totalMeshCounter+=target.CountFacets
            #check if this target is a potential intersection
            if CAD.intersectList[i] in self.intersects:
                numTargetFaces += target.CountFacets
                #append target data
                for face in target.Facets:
                    targetPoints.append(face.Points)
                    targetNorms.append(face.Normal)
        targetPoints = np.asarray(targetPoints)/1000.0 #scale to m
        targetNorms = np.asarray(targetNorms)
        print("TOTAL INTERSECTION FACES: {:d}".format(totalMeshCounter))
        #for debugging, save a shadowmask at each step up fieldline
        if shadowMaskClouds == True:
            self.write_shadow_pointcloud(self.centers,self.shadowed_mask,self.controlfilePath,tag='original')

        #MAFOT always returns ordered points in CCW (from top) direction,
        #so we need to check which direction we are running so we read the output
        #correctly
        if self.mapDirectionStruct == -1:
            print('Tracing with reversed Map Direction')
            log.info('Tracing with reversed Map Direction')
            startIdx = 0
        else:
            print('Tracing with forward Map Direction')
            log.info('Tracing with forward Map Direction')
            startIdx = 1

        #===INTERSECTION TEST 1 (tricky frontface culling / first step up field line)
        dphi = 1.0
        MHD.ittStruct = 1.0
        numSteps = MHD.nTrace #actual trace is (numSteps + 1)*dphi degrees
        #If numSteps = 0, dont do intersection checking
        if numSteps > 0:
            CTLfile = self.controlfilePath + self.controlfileStruct
            MHD.writeControlFile(CTLfile, self.t, self.mapDirectionStruct, mode='struct')

            #Perform first integration step
            MHD.writeMAFOTpointfile(self.centers[use],self.gridfileStruct)
            MHD.getMultipleFieldPaths(dphi, self.gridfileStruct, self.controlfilePath, self.controlfileStruct)
            structData = self.readStructOutput(self.structOutfile)
            os.remove(self.structOutfile) #clean up

            #this is for printing information about a specific mesh element
#            ptIdx = np.where(use==16886)[0]
            ptIdx = None

            #First we do basic intersection checking
            intersect_mask = self.intersectTestBasic(structData,self.norms[use],
                                                    targetPoints,targetNorms,
                                                    MHD,self.ep,self.mapDirectionStruct,
                                                    ptIdx)

            self.shadowed_mask[use] = intersect_mask

            #for debugging, save a shadowmask at each step up fieldline
            if shadowMaskClouds == True:
                self.write_shadow_pointcloud(self.centers,self.shadowed_mask,self.controlfilePath,tag='test0')

        #===INTERSECTION TEST 2 (multiple steps up field line)
        #Starts at second step up field line
        if numSteps > 1:
            MHD.writeControlFile(self.controlfileStruct, self.t, self.mapDirectionStruct, mode='struct')
            use = np.where(self.shadowed_mask == 0)[0]
            intersect_mask2 = np.zeros((len(use)))

            #Perform fist integration step but dont use for finding intersections
            MHD.writeMAFOTpointfile(self.centers[use],self.gridfileStruct)
            MHD.getMultipleFieldPaths(dphi, self.gridfileStruct, self.controlfilePath, self.controlfileStruct)
            structData = self.readStructOutput(self.structOutfile)
            os.remove(self.structOutfile) #clean up

            #Perform subsequent integration steps.  Use the point we left off at in
            #last loop iteration as the point we launch from in next loop iteration
            #This amounts to 'walking' up the field line looking for intersections,
            #which is important when field line curvature makes intersections happen
            #farther than 1-2 degrees from PFC surface.
            use2 = np.where(intersect_mask2 == 0)[0]
            for i in range(numSteps):
                print("\nIntersect Trace #2 Step {:d}".format(i))
                log.info("\nIntersect Trace #2 Step {:d}".format(i))
                useOld = use2
                use2 = np.where(intersect_mask2 == 0)[0]
                indexes = np.where([x==useOld for x in use2])[1] #map current steps's index back to intersect_mask2 index

                StartPoints = structData[startIdx::2,:][indexes] #odd indexes are second trace point
                MHD.writeMAFOTpointfile(StartPoints,self.gridfileStruct)
                MHD.getMultipleFieldPaths(dphi, self.gridfileStruct, self.controlfilePath, self.controlfileStruct)
                structData = self.readStructOutput(self.structOutfile)
                os.remove(self.structOutfile) #clean up
                intersect_mask2[use2] = self.intersectTest2(structData,targetPoints,self.mapDirectionStruct)

                #for debugging, save a shadowmask at each step up fieldline
                if shadowMaskClouds == True:
                    self.shadowed_mask[use] = intersect_mask2
                    self.write_shadow_pointcloud(self.centers,self.shadowed_mask,self.controlfilePath,tag='test{:d}'.format(i+1))

            #Now revise shadowed_mask taking intersections into account
            self.shadowed_mask[use] = intersect_mask2
        return

    def longRangeIntersectCheck(self,qDiv,thresh,distPhi,MHD,CAD):
        """
        Runs a long range intersection check.

        Because some PFC components may have STL faces that are in locations
        with very high Bt/(Br+Bz) ratio, it might take longer than 10-20 degrees
        for them to intersect with anything.

        This function requires a threshold power.  Any PFC face with power
        higher than the threshold is checked.

        The flagged faces are then checked for intersections over distPhi degrees
        """
        print("\nLong Range Trace Begins")
        log.info("\nLong Range Trace Begins")

        numTargetFaces = 0
        targetPoints = []
        print('Number of target parts: {:f}'.format(len(self.intersects)))
        for i,target in enumerate(CAD.intersectMeshes):
            #check if this target is a potential intersection
            if CAD.intersectList[i] in self.intersects:
                numTargetFaces += target.CountFacets
                #append target data
                for face in target.Facets:
                    targetPoints.append(face.Points)
        targetPoints = np.asarray(targetPoints)/1000.0 #scale to m

        use = np.where(np.abs(qDiv) > thresh)[0]
        found = 0
        print('Number of points above threshold: {:f}'.format(len(use)))
        log.info('Number of points above threshold: {:f}'.format(len(use)))
        for idx in use:
            startPoint = self.centers[idx]
#            print("Start Point")
#            print(startPoint)
            MHD.writeMAFOTpointfile(startPoint,self.gridfileStruct)
            MHD.getMultipleFieldPaths(distPhi, self.gridfileStruct, self.controlfilePath, self.controlfileStruct)
            structData = self.readStructOutput(self.structOutfile)
            os.remove(self.structOutfile) #clean up
            #create source format that intersectTest2 reads (odd start points, even end points)
            x = np.insert(structData[:,0], np.arange(len(structData[:,0])), structData[:,0])
            y = np.insert(structData[:,1], np.arange(len(structData[:,1])), structData[:,1])
            z = np.insert(structData[:,2], np.arange(len(structData[:,2])), structData[:,2])
            sources = np.array([x,y,z]).T[1:-1]
            intersects = self.intersectTest2(sources,targetPoints,self.mapDirectionStruct)
            if np.sum(intersects) > 0:
#                print("Found long range intersection")
#                log.info("Found long range intersection")
                qDiv[idx] = 0.0
                found +=1
        print("Long range trace found {:f} intersections".format(found))
        log.info("Long range trace found {:f} intersections".format(found))
        return qDiv

    def intersectTestBasic(self,sources,sourceNorms,
                        targets,targetNorms,MHD,ep,ionDir,ptIdx=None):
        """
        checks if any of the lines (field line traces) generated by MAFOT
        struct program intersect any of the target mesh faces.

        source represents the face where we want to find the heat flux on.
        here it is a 'source' for the field line launching that happens in
        MAFOT.  target represents the potential intersection faces, where a
        field line launched from the source face may intersect.

        Returns a boolean matrix (bitmask) of shape (N) that is true where
        sourceface i intersects with any targetface

        This function first checks for intersection with all of the faces that are
        labeled as 'shadowed' by the backface culling algorithm.  This will
        capture the majority of the intersections, but will miss intricate
        details like when the poloidal field reverses in a gap.

        The function then checks for intersections for the sources that the
        aforementioned check deemed to be heat loaded (no intersection), against
        all of the target faces.  It eliminates self intersections during this
        step by making sure that if an intersection occurs, it is with a face
        that is far away in any direction
        """
        q1 = sources[::2,:]  #even indexes are first trace point
        q2 = sources[1::2,:]  #odd indexes are second trace point
        p1 = targets[:,0,:]  #point 1 of mesh triangle
        p2 = targets[:,1,:]  #point 2 of mesh triangle
        p3 = targets[:,2,:]  #point 3 of mesh triangle

        N = len(q2)
        Nt = len(p1)
        print('{:d} Source Faces and {:d} Target Faces'.format(N,Nt))


        #Cull target front face (PFC front faces will only intersect with
        # other PFC shadowed faces)
        x = np.zeros((Nt,3))
        y = np.zeros((Nt,3))
        z = np.zeros((Nt,3))
        x[:,0] = p1[:,0]
        x[:,1] = p2[:,0]
        x[:,2] = p3[:,0]
        y[:,0] = p1[:,1]
        y[:,1] = p2[:,1]
        y[:,2] = p3[:,1]
        z[:,0] = p1[:,2]
        z[:,1] = p2[:,2]
        z[:,2] = p3[:,2]
        tools.targetCtrs = tools.faceCenters(x,y,z)
        targetBackfaceMask = self.backfaceCulling(tools.targetCtrs,targetNorms,MHD,ep,ionDir)

        #Potential intersection faces
        use = np.where(targetBackfaceMask == 1)[0]
#        use = np.arange(0,Nt)
        Nt_use = len(use)

        tools.q1 = q1
        tools.q2 = q2
        tools.p1 = p1[use]
        tools.p2 = p2[use]
        tools.p3 = p3[use]
        tools.Nt = Nt_use

        print('Entering basic intersection test for {:d} potential intersection faces'.format(Nt_use))
        log.info('Entering basic intersection test for {:d} potential intersection faces'.format(Nt_use))
        t0 = time.time()

        #Prepare intersectionTest across multiple cores
        Ncores = multiprocessing.cpu_count() -2 #reserve 2 cores for overhead
        print('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        log.info('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        #each worker receives a single start and end point (q1 and q2),
        #corresponding to one trace from the MAFOT structure output.
        #The worker accesses the many potential intersection triangle vertices
        #through tools class variables (p1,p2,p3).  Making p1,p2,p3 tools class
        #variables eliminates the overhead of transmitting these matrices to
        #each worker, which yields about an order of magnitude speedup.
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
        pool = multiprocessing.Pool(Ncores)
        mask = np.asarray(pool.map(tools.intersectTestParallel, np.arange(N)))
        pool.close()
        print('Found {:f} shadowed faces'.format(np.sum(mask)))
        log.info('Found {:f} shadowed faces'.format(np.sum(mask)))
        print('Time elapsed: {:f}'.format(time.time() - t0))
        log.info('Time elapsed: {:f}'.format(time.time() - t0))



        #Now remove trickier to identify intersections.  Test 2
        #Remove intersections that are outside of 5mm from us.
        #This is run with 'use' representing all faces that were identified
        #as not shadowed (mask = 0) above.
#        use = np.where( mask == 0 )[0]
#        tools.q1 = q1
#        tools.q2 = q2
#        tools.p1 = p1
#        tools.p2 = p2
#        tools.p3 = p3
#        tools.Nt = Nt
#        print('Entering intersection Test #2 for {:d} potential intersection faces'.format(Nt))
#        log.info('Entering intersection Test #2 for {:d} potential intersection faces'.format(Nt))
#        t0 = time.time()
#        print('Spawning tasks to workers')
#        log.info('Spawning tasks to workers')
#        pool = multiprocessing.Pool(Ncores)
#        mask[use] = np.asarray(pool.map(tools.intersectTestParallel_selfCheck, use))
#        pool.close()
#        print('Found {:f} shadowed faces'.format(np.sum(mask[use])))
#        log.info('Found {:f} shadowed faces'.format(np.sum(mask[use])))
#        print('Time elapsed: {:f}'.format(time.time() - t0))
#        log.info('Time elapsed: {:f}'.format(time.time() - t0))

        return mask




    def intersectTest2(self,sources,targets,mapDirection):
        """
        Run an intersection test against all possible source faces
        sources is endpoints of line
        targets is points that make up triangle faces we check for intersections
           against

        This test is called repeatedly from findShadows_structure, as we integrate
        up a field line.  Each time it is called, it checks if the sources line
        intersects a targets face
        """
        q1 = sources[::2,:] #even indexes are first trace point
        q2 = sources[1::2,:] #odd indexes are second trace point
        p1 = targets[:,0,:] #point 1 of mesh triangle
        p2 = targets[:,1,:] #point 2 of mesh triangle
        p3 = targets[:,2,:] #point 3 of mesh triangle

        tools.q1 = q1
        tools.q2 = q2
        tools.p1 = p1
        tools.p2 = p2
        tools.p3 = p3
        N = len(q2)
        Nt = len(p1)
        tools.Nt = Nt

        print('Entering intersection Test #1 for {:d} potential intersection faces'.format(Nt))
        log.info('Entering intersection Test #1 for {:d} potential intersection faces'.format(Nt))
        t0 = time.time()

        #Prepare intersectionTest across multiple cores
        Ncores = multiprocessing.cpu_count() -2 #reserve 2 cores for overhead
        print('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        log.info('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        #each worker receives a single start and end point (q1 and q2),
        #corresponding to one trace from the MAFOT structure output.
        #The worker accesses the many potential intersection triangle vertices
        #through tools class variables (p1,p2,p3).  Making p1,p2,p3 tools class
        #variables eliminates the overhead of transmitting these matrices to
        #each worker, which yields about an order of magnitude speedup.
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
        pool = multiprocessing.Pool(Ncores)
        mask = np.asarray(pool.map(tools.intersectTestParallel, np.arange(N)))
        pool.close()

        print('Found {:f} shadowed faces'.format(np.sum(mask)))
        log.info('Found {:f} shadowed faces'.format(np.sum(mask)))
        print('Time elapsed: {:f}'.format(time.time() - t0))
        log.info('Time elapsed: {:f}'.format(time.time() - t0))

        return mask

    def readStructOutput(self,file):
        """
        Reads output file from MAFOT structure program
        """
        structdata = np.genfromtxt(file,comments='#')
        xyz = np.zeros((len(structdata),3))
        xyz[:,0] = structdata[:,0]
        xyz[:,1] = structdata[:,1]
        xyz[:,2] = structdata[:,2]
        #remove rows with zeros (invalids)
        #xyzFinal = xyz[~(np.abs(xyz)<1e-99).any(axis=1)]
        print('Read structure output for {:d} points'.format(len(xyz)))
        log.info('Read structure output for {:d} points'.format(len(xyz)))
        return xyz

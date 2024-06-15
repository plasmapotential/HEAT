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
#from scipy.spatial import ConvexHull


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

    def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
        """
        sets up pythonic string number formats for shot and timesteps
        """
        self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
        self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
        return

    def allowed_class_vars(self):
        """
        .. Writes a list of recognized class variables to HEAT object
        .. Used for error checking input files and for initialization

        Photon Radiation Heat Flux Variables:
        -------------------------------------

        :radFile: CSV file path where the photon emission data is stored.  Photon emission
          data should be provided in a CSV file with columns R,Z,MW, corresonding to an
          axisymmetric photon emission profile given in [MW].  Path should be absolute path.
        :Ntor: Number of times the axisymmetric profile is repeated toroidally.  In the limit
          that Ntor goes to infinity, emission profile becomes toroidally continuous.
        :Nref: Number of reflections to trace for.  Currently reflections are not implemented
          in HEAT so this should be set to 1
        :phiMin: Minimum toroidal angle of emission extent [degrees].  The emission profile
          read from a file will be duplicated Ntor times between phiMin and phiMax.
        :phiMax: Maximum toroidal angle of emission extent [degrees].  The emission profile
          read from a file will be duplicated Ntor times between phiMin and phiMax.

        
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
        try:
            self.Ntor = int(self.Ntor)
            self.Nref = int(self.Nref)
            self.phiMin = float(self.phiMin)
            self.phiMax = float(self.phiMax)
        except:
            print("Could not initialize RadPower variables.  Bailing...")
            log.info("Could not initialize RadPower variables.  Bailing...")
        return

    def read2DSourceFile(self, file):
        """
        Reads a comma delimited file that describes the radiated power on an
        R,Z grid. This point cloud (pc) should represent the axisymmetric
        plasma radiated power across a single poloidal plane (ie sliced at phi).
        The power defined at each RZ point is the radiated photon power.

        It is assumed that the RZ points are located at the centroids of mesh
        elements

        Each row of this file corresponds to a separate source point.  The columns
        of the file are:

        :R: Radial coordinate of the source point [m]
        :Z: Z coordinate of the source point [m]
        :MW: The power [MW] associated with that mesh point.  This power represents the 
          toroidally integrated power.  In other words, the total power of a toroidally
          revolved mesh element centered at R,Z

        """
        print("Reading 2D photon radiation source file: "+file)
        log.info("Reading 2D photon radiation source file: "+file)
        self.PC2D = pd.read_csv(file, header=0, names=['R','Z','P']).values 
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
        PC2D = np.append( self.PC2D, np.zeros([NR,1]), 1 ) #add column for phi
        self.PC3D = np.repeat(PC2D[:,:,np.newaxis],self.Nphi,axis=2) #repeat Nphi times
        self.PC3D[:,-1,:] = np.repeat(self.phis[np.newaxis], NR, axis=0) #fill in phis
        #normalize power to the number of toroidal steps and toroidal section width, deltaPhi
        self.PC3D[:,2,:] = self.PC3D[:,2,:] * (self.deltaPhi / 360.0) / self.Nphi

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
        self.phis = np.linspace(phiMin, phiMax, Ntor+2)[1:-1]
        self.Nphi = len(self.phis)
        self.deltaPhi = phiMax - phiMin
        print("Radiated Power Source Angles:")
        print(self.phis)
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
                log.info("Cannot count faces because "+CAD.intersectList[i]+" is not a mesh!")
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
            #self.meshFile = '/home/tom/source/dummyOutput/SOLID843___5.00mm.stl' #for testing

        elif mode=='mitsuba':
            combinedMesh = CAD.createEmptyMesh()
            for m in targetMeshes:
                combinedMesh.addFacets(m.Facets)
            oldMask = CAD.overWriteMask
            CAD.overWriteMask = True
            #mitsuba only reads PLYs!!!
            CAD.writeMesh2file(combinedMesh, 'combinedMesh', path=CAD.STLpath, resolution='standard', fType='ply')
            CAD.overWriteMask = oldMask
            self.meshFile = CAD.STLpath + 'combinedMesh' + "___standard.ply"    


        print("\nTotal Rad Intersection Faces: {:d}".format(totalMeshCounter))
        print("Rad Intersect Faces for this PFC: {:d}".format(numTargetFaces))
        print("# of radiated power point sources: {:d}".format(self.Ni))
        print("# of PFC ROI mesh elements: {:d}".format(self.Nj))
        print("Total number of source-target ray-tracing calculations: {:d}\n".format(self.Ni*self.Nj))
        print("Running intersection check...")
        log.info("\nTotal Rad Intersection Faces: {:d}".format(totalMeshCounter))
        log.info("Rad Intersect Faces for this PFC: {:d}".format(numTargetFaces))
        log.info("# of radiated power point sources: {:d}".format(self.Ni))
        log.info("# of PFC ROI mesh elements: {:d}".format(self.Nj))
        log.info("Total number of source-target ray-tracing calculations: {:d}\n".format(self.Ni*self.Nj))
        log.info("Running intersection check...")
        return



    def calculatePowerTransferMitsubaLoop(self, mode=None, mitsubaMode='llvm', fType='ply', batch_size=1000):
        """
        Maps power between sources and targets (ROI PFCs).  Uses Mitsuba3 to
        perform ray tracing.  Mitsuba3 can be optimized for CPU or GPU.

        Code largely developed by A. Rosenthal (CFS)
        Adapted to HEAT by T. Looby

        Mitsuba:
        @software{jakob2022mitsuba3,
            title = {Mitsuba 3 renderer},
            author = {Wenzel Jakob and Sébastien Speierer and Nicolas Roussel and Merlin Nimier-David and Delio Vicini and Tizian Zeltner and Baptiste Nicolet and Miguel Crespo and Vincent Leroy and Ziyi Zhang},
            note = {https://mitsuba-renderer.org},
            version = {3.0.1},
            year = 2022,
        }

        """
        import drjit as dr
        import mitsuba as mi


        t0 = time.time()
        powerFrac = np.zeros((self.Ni,self.Nj))
        Psum = np.zeros((self.Nj))
        self.hullPower = np.zeros((self.Ni))
        self.powerFrac = np.zeros((self.Ni))

        print("Building radiation scene...")    

        # Initialize Mitsuba
        if mitsubaMode == 'cuda':
            mi.set_variant('cuda_ad_rgb')
        else:
            mi.set_variant('llvm_ad_rgb')

        scene = {'type': 'scene'}
        scene['path_int'] = {
            'type': 'path',
            'max_depth': 1  # Maximum number of bounces
        }

        name = 'allMesh' #name of entire mesh we saved

        scene[name] = {
            'type': fType,
            'filename': self.meshFile,
            'to_world': mi.ScalarTransform4f.scale(0.001),  # Convert from mm to m
            'face_normals': True,  # This prevents smoothing of sharp-corners by discarding surface-normals. Useful for engineering CAD.
            'sensor': {
                'type': 'irradiancemeter',
                'film': {
                    'type': 'hdrfilm',
                    'pixel_format': 'luminance',
                    'filter': {'type': 'box'},
                    'width': 1,
                    'height': 1,
                }
            }
        }

        t = time.time()
        print('Loading Scene....')
        scene = mi.load_dict(scene)
        print('Time to load scene: ' + str(time.time() - t))
        print('\n')

        t = time.time()
        # Convert the numpy arrays to Dr.Jit dynamic arrays
        #tx,tx,tz,tp are source values for x,y,z,power
        tx = dr.cuda.Float(self.sources[:,0]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,0])
        ty = dr.cuda.Float(self.sources[:,1]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,1])
        tz = dr.cuda.Float(self.sources[:,2]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,2])
        tp = dr.cuda.Float(self.sourcePower) if mode == 'cuda' else dr.llvm.Float(self.sourcePower)

        # Grab all the shape ids to correctly identify the sensor
        idL = [x.id() for x in scene.shapes()]

        print('Time to build sources: ' + str(time.time() - t))

        t = time.time()
        # Find the sensor
        sensInd = idL.index(name)
        sensS = scene.shapes()[sensInd]

        params = mi.traverse(sensS)
        verts = params['vertex_positions']
        norms = params['vertex_normals']

        # Find all the face and areas
        totFace = sensS.face_count()
        if mitsubaMode == 'cuda':
            indF = sensS.face_indices(dr.arange(dr.cuda.ad.UInt, totFace))
        else:
            indF = sensS.face_indices(dr.arange(dr.llvm.ad.UInt, totFace))

        # Label each vertex with its index
        vertNum = int(len(verts) / 3)
        vertA = np.array(verts)

        vertAR = vertA.reshape((vertNum, 3))
        indFA = np.array(indF)

        # Find the vertices of each face
        faceV = np.take(vertAR, indFA.flatten(), axis=0)

        # Find the area of each face and the mean point of each face
        faceV1 = faceV[::3, :]
        faceV2 = faceV[1::3, :]
        faceV3 = faceV[2::3, :]

        # Find the area of each face
        cross1 = faceV2 - faceV1
        cross2 = faceV3 - faceV1

        cross = np.cross(cross1, cross2)
        area = np.linalg.norm(cross, axis=1) / 2

        #find the mean point of each face
        #again can probably be sped up
        center = np.asarray([x.mean(0) for x in np.split(faceV,totFace,axis = 0)])

        #dynamically allocate batch size based upon available memory
        if batch_size == "auto":
            batch_size = int(2**32 / totFace  * 0.5)

        print("Using batch size of: {:d}".format(batch_size))
        
        N_sources = len(tx)
        N_targets = len(center) #all intersections, including those not in ROI
        targetPower = np.zeros((N_targets))
        print('Time to prepare mesh: ' + str(time.time() - t))

        t0 = time.time()
        
        self.powCount = np.zeros((totFace))
        

        # Initialize the loop state (listing all variables that are modified inside the loop)
        #loopFlag = mi.UInt32(0)
        #loop = mi.Loop(name="", state=lambda: (loopFlag))
        #i=0
        idx = mi.UInt32(0)
        i=0
        loop = mi.Loop(name="Example Loop", state=lambda: (idx))
        N = int(N_sources / batch_size)
        N_dr = mi.UInt32(N)

        while loop(idx<N_dr):
            self.mitsubaLoop(N_sources, batch_size, totFace, mitsubaMode,
                    center, tx, ty, tz, tp, area, scene, idx)
        
        self.targetPower = self.powCount[:self.Nj]
        print('Mitsuba Calc Time: '+str(time.time()-t0))

        return
    
    def mitsubaLoop(self, N_sources, batch_size, totFace, mitsubaMode,
                    center, tx, ty, tz, tp, area, scene, idx):

        import drjit as dr
        import mitsuba as mi

        i = int(idx[0])*batch_size

        tLoop = time.time()
        #dynamically adjust size of batch depending on how many
        #sources are left to calculate
        if batch_size == 1:
            size_i = batch_size
        elif (N_sources - (batch_size+i)) < 0:
            size_i = np.abs(N_sources - i)
            #if we change size_i, the powCount array will be a different 
            #dimension.  so we take the existing powCount and export it
            #to the targetPower array before changing the size for the 
            #last loop iteration where size_i changes
            targetPower += np.sum(np.array(powCount.copy()).reshape(size_i, totFace), axis=0)
            #now change the size of powCount for the last iteration
            powCount = dr.zeros(dr.llvm.ad.Float, totFace*size_i)
        else:
            size_i = batch_size
        #initialize the power tally if this is the first iteration
        if i==0:
            powCount = dr.zeros(dr.llvm.ad.Float, totFace*size_i)
        print("\n=== Running Sources: {:d} - {:d} ===".format(i, i+size_i))
        dummy = np.ones((size_i))
        #take all the x values from the sensor and pair with all the x values from the power distribution
        if mitsubaMode == 'cuda':
            sx_batch,tx_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,0]),dr.cuda.ad.Float(tx[i:i+size_i]),indexing = 'xy') #indexing ij makes it so that the first numpower points of sx are all the same (tx[0])
            sy_batch,ty_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,1]),dr.cuda.ad.Float(ty[i:i+size_i]),indexing = 'xy')
            sz_batch,tz_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,2]),dr.cuda.ad.Float(tz[i:i+size_i]),indexing = 'xy')
            _,tp_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,0]),dr.cuda.ad.Float(tp[i:i+size_i]),indexing = 'xy') #make the power meshed in the same way
            sA,_ = dr.meshgrid(dr.cuda.ad.Float(area),dr.cuda.ad.Float(dummy),indexing = 'xy')
        else:
            sx_batch,tx_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,0]),dr.llvm.ad.Float(tx[i:i+size_i]),indexing = 'xy') #indexing ij makes it so that the first numpower points of sx are all the same (tx[0])
            sy_batch,ty_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,1]),dr.llvm.ad.Float(ty[i:i+size_i]),indexing = 'xy')
            sz_batch,tz_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,2]),dr.llvm.ad.Float(tz[i:i+size_i]),indexing = 'xy')
            _,tp_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,0]),dr.llvm.ad.Float(tp[i:i+size_i]),indexing = 'xy') #make the power meshed in the same way
            sA,_ = dr.meshgrid(dr.llvm.ad.Float(area),dr.llvm.ad.Float(dummy),indexing = 'xy')
        #print('Time to meshgrid: ' + str(time.time() - t))
        #define the ray origins which start at the HEAT grid points
        origin = mi.Point3f(tx_batch,ty_batch,tz_batch)
        #and the target points which are the sensor faces
        target = mi.Vector3f(sx_batch,sy_batch,sz_batch)
        direction = dr.normalize(target-origin)
        ray = mi.Ray3f(o=origin, d=direction)
        #t = time.time()
        si = scene.ray_intersect(ray)
        #valid = si.is_valid()
        #active = mi.Bool(valid)
        #print('Time to run intersection: '+str(time.time()-t))
        
        #gives the index of the primitive triangle that was hit
        primI = si.prim_index
        pow = dr.abs(dr.dot(dr.normalize(si.n),direction)) * sA * tp_batch/(4*dr.pi*dr.power(si.t,2))
        
        #finitePow = dr.isfinite(pow)
        piA = dr.arange(dr.llvm.ad.UInt,totFace)
        piA = dr.tile(piA,size_i)
        mask = dr.eq(primI,piA)
        #mask2 = dr.eq(mask, active)
        #mask3 = dr.eq(mask2, finitePow)
        correctHit = dr.compress(mask)
        #print('t6: '+str(time.time()-t))
        powTmp = dr.zeros(dr.llvm.ad.Float, totFace*size_i)
        tmp = dr.gather(type(pow),source = pow, index = correctHit) 
        dr.scatter(powTmp, value=tmp, index=correctHit)
        ##for troubleshooting.  print data for a source -> target trace
        #srcIdx = 18289
        #roiIdx = 1700
        #if np.logical_and(srcIdx > i, srcIdx < i+size_i):
        #    idxBatch = srcIdx - i
        #    print(np.array(tp_batch).reshape(size_i, totFace)[idxBatch, roiIdx])
        #    print(np.array(pow).reshape(size_i, totFace)[idxBatch,roiIdx])
        #    print(np.array(powTmp).reshape(size_i, totFace)[idxBatch, roiIdx])
        #    input()
        powCount += powTmp
        targetPower = np.sum(np.array(powCount).reshape(size_i, totFace), axis=0)
        print("Batch Loop Time: {:f}[s]".format(time.time() - tLoop))
        self.powCount += targetPower

        idx += size_i
        return 


    def calculatePowerTransferMitsubaJIT(self, mode=None, mitsubaMode='llvm', fType='ply', batch_size=1000):
        """
        Maps power between sources and targets (ROI PFCs).  Uses Mitsuba3 to
        perform ray tracing.  Mitsuba3 can be optimized for CPU or GPU.

        Code largely developed by A. Rosenthal (CFS)
        Adapted to HEAT by T. Looby

        Mitsuba:
        @software{jakob2022mitsuba3,
            title = {Mitsuba 3 renderer},
            author = {Wenzel Jakob and Sébastien Speierer and Nicolas Roussel and Merlin Nimier-David and Delio Vicini and Tizian Zeltner and Baptiste Nicolet and Miguel Crespo and Vincent Leroy and Ziyi Zhang},
            note = {https://mitsuba-renderer.org},
            version = {3.0.1},
            year = 2022,
        }

        """
        import drjit as dr
        import mitsuba as mi

        t0 = time.time()
        powerFrac = np.zeros((self.Ni,self.Nj))
        Psum = np.zeros((self.Nj))
        self.hullPower = np.zeros((self.Ni))
        self.powerFrac = np.zeros((self.Ni))

        print("Building radiation scene...")    

        # Initialize Mitsuba
        if mitsubaMode == 'cuda':
            mi.set_variant('cuda_ad_rgb')
        else:
            mi.set_variant('llvm_ad_rgb')

        scene = {'type': 'scene'}
        scene['path_int'] = {
            'type': 'path',
            'max_depth': 1  # Maximum number of bounces
        }

        name = 'allMesh' #name of entire mesh we saved

        scene[name] = {
            'type': fType,
            'filename': self.meshFile,
            'to_world': mi.ScalarTransform4f.scale(0.001),  # Convert from mm to m
            'face_normals': True,  # This prevents smoothing of sharp-corners by discarding surface-normals. Useful for engineering CAD.
            'sensor': {
                'type': 'irradiancemeter',
                'film': {
                    'type': 'hdrfilm',
                    'pixel_format': 'luminance',
                    'filter': {'type': 'box'},
                    'width': 1,
                    'height': 1,
                }
            }
        }

        t = time.time()
        print('Loading Scene....')
        scene = mi.load_dict(scene)
        print('Time to load scene: ' + str(time.time() - t))
        print('\n')

        t = time.time()
        # Convert the numpy arrays to Dr.Jit dynamic arrays
        #tx,tx,tz,tp are source values for x,y,z,power
        tx = dr.cuda.Float(self.sources[:,0]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,0])
        ty = dr.cuda.Float(self.sources[:,1]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,1])
        tz = dr.cuda.Float(self.sources[:,2]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,2])
        tp = dr.cuda.Float(self.sourcePower) if mode == 'cuda' else dr.llvm.Float(self.sourcePower)

        # Grab all the shape ids to correctly identify the sensor
        idL = [x.id() for x in scene.shapes()]

        print('Time to build sources: ' + str(time.time() - t))

        t = time.time()
        # Find the sensor
        sensInd = idL.index(name)
        sensS = scene.shapes()[sensInd]

        params = mi.traverse(sensS)
        verts = params['vertex_positions']
        norms = params['vertex_normals']

        # Find all the face and areas
        totFace = sensS.face_count()
        if mitsubaMode == 'cuda':
            indF = sensS.face_indices(dr.arange(dr.cuda.ad.UInt, totFace))
        else:
            indF = sensS.face_indices(dr.arange(dr.llvm.ad.UInt, totFace))

        # Label each vertex with its index
        vertNum = int(len(verts) / 3)
        vertA = np.array(verts)

        vertAR = vertA.reshape((vertNum, 3))
        indFA = np.array(indF)

        # Find the vertices of each face
        faceV = np.take(vertAR, indFA.flatten(), axis=0)

        # Find the area of each face and the mean point of each face
        faceV1 = faceV[::3, :]
        faceV2 = faceV[1::3, :]
        faceV3 = faceV[2::3, :]

        # Find the area of each face
        cross1 = faceV2 - faceV1
        cross2 = faceV3 - faceV1

        cross = np.cross(cross1, cross2)
        area = np.linalg.norm(cross, axis=1) / 2

        #find the mean point of each face
        #again can probably be sped up
        center = np.asarray([x.mean(0) for x in np.split(faceV,totFace,axis = 0)])

        #dynamically allocate batch size based upon available memory
        if batch_size == "auto":
            batch_size = int(2**32 / totFace  * 0.5)

        print("Using batch size of: {:d}".format(batch_size))
        
        N_sources = len(tx)
        N_targets = len(center) #all intersections, including those not in ROI
        targetPower = np.zeros((N_targets))
        print('Time to prepare mesh: ' + str(time.time() - t))

        t0 = time.time()
        
        targetPower = np.zeros((totFace))

        # Initialize the loop state (listing all variables that are modified inside the loop)
        #loopFlag = mi.UInt32(0)
        #loop = mi.Loop(name="", state=lambda: (loopFlag))
        #i=0

        for i in range(0, N_sources, batch_size):
        #while loop(loopFlag == 0):
            tLoop = time.time()

            #dynamically adjust size of batch depending on how many
            #sources are left to calculate
            if batch_size == 1:
                size_i = batch_size
            elif (N_sources - (batch_size+i)) < 0:
                size_i = np.abs(N_sources - i)
                #if we change size_i, the powCount array will be a different 
                #dimension.  so we take the existing powCount and export it
                #to the targetPower array before changing the size for the 
                #last loop iteration where size_i changes
                targetPower += np.sum(np.array(powCount.copy()).reshape(size_i, totFace), axis=0)
                #now change the size of powCount for the last iteration
                if mitsubaMode == 'cuda':
                    powCount = dr.zeros(dr.cuda.ad.Float, totFace*size_i)
                else:
                    powCount = dr.zeros(dr.llvm.ad.Float, totFace*size_i)
            else:
                size_i = batch_size


            #initialize the power tally if this is the first iteration
            if i==0:
                if mitsubaMode == 'cuda':
                    powCount = dr.zeros(dr.llvm.ad.Float, totFace*size_i)
                else:
                    powCount = dr.zeros(dr.llvm.ad.Float, totFace*size_i)

            print("\n=== Running Sources: {:d} - {:d} ===".format(i, i+size_i))
            dummy = np.ones((size_i))

            #take all the x values from the sensor and pair with all the x values from the power distribution
            if mitsubaMode == 'cuda':
                sx_batch,tx_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,0]),dr.cuda.ad.Float(tx[i:i+size_i]),indexing = 'xy') #indexing ij makes it so that the first numpower points of sx are all the same (tx[0])
                sy_batch,ty_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,1]),dr.cuda.ad.Float(ty[i:i+size_i]),indexing = 'xy')
                sz_batch,tz_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,2]),dr.cuda.ad.Float(tz[i:i+size_i]),indexing = 'xy')
                _,tp_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,0]),dr.cuda.ad.Float(tp[i:i+size_i]),indexing = 'xy') #make the power meshed in the same way
                sA,_ = dr.meshgrid(dr.cuda.ad.Float(area),dr.cuda.ad.Float(dummy),indexing = 'xy')
            else:
                sx_batch,tx_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,0]),dr.llvm.ad.Float(tx[i:i+size_i]),indexing = 'xy') #indexing ij makes it so that the first numpower points of sx are all the same (tx[0])
                sy_batch,ty_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,1]),dr.llvm.ad.Float(ty[i:i+size_i]),indexing = 'xy')
                sz_batch,tz_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,2]),dr.llvm.ad.Float(tz[i:i+size_i]),indexing = 'xy')
                _,tp_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,0]),dr.llvm.ad.Float(tp[i:i+size_i]),indexing = 'xy') #make the power meshed in the same way
                sA,_ = dr.meshgrid(dr.llvm.ad.Float(area),dr.llvm.ad.Float(dummy),indexing = 'xy')

            #print('Time to meshgrid: ' + str(time.time() - t))

            #define the ray origins which start at the HEAT grid points
            origin = mi.Point3f(tx_batch,ty_batch,tz_batch)
            #and the target points which are the sensor faces
            target = mi.Vector3f(sx_batch,sy_batch,sz_batch)
            direction = dr.normalize(target-origin)
            ray = mi.Ray3f(o=origin, d=direction)

            #t = time.time()
            si = scene.ray_intersect(ray)
            #valid = si.is_valid()
            #active = mi.Bool(valid)
            #print('Time to run intersection: '+str(time.time()-t))
            
            #gives the index of the primitive triangle that was hit
            primI = si.prim_index
            pow = dr.abs(dr.dot(dr.normalize(si.n),direction)) * sA * tp_batch/(4*dr.pi*dr.power(si.t,2))
            
            #finitePow = dr.isfinite(pow)
            if mitsubaMode == 'cuda':
                piA = dr.arange(dr.cuda.ad.UInt,totFace)
                powTmp = dr.zeros(dr.cuda.ad.Float, totFace*size_i)
            else:
                piA = dr.arange(dr.llvm.ad.UInt,totFace)
                powTmp = dr.zeros(dr.llvm.ad.Float, totFace*size_i)
            piA = dr.tile(piA,size_i)
            mask = dr.eq(primI,piA)
            #mask2 = dr.eq(mask, active)
            #mask3 = dr.eq(mask2, finitePow)
            correctHit = dr.compress(mask)

            #print('t6: '+str(time.time()-t))
            tmp = dr.gather(type(pow),source = pow, index = correctHit) 
            dr.scatter(powTmp, value=tmp, index=correctHit)


            ##for troubleshooting.  print data for a source -> target trace
            #srcIdx = 18289
            #roiIdx = 1700
            #if np.logical_and(srcIdx > i, srcIdx < i+size_i):
            #    idxBatch = srcIdx - i
            #    print(np.array(tp_batch).reshape(size_i, totFace)[idxBatch, roiIdx])
            #    print(np.array(pow).reshape(size_i, totFace)[idxBatch,roiIdx])
            #    print(np.array(powTmp).reshape(size_i, totFace)[idxBatch, roiIdx])
            #    input()

            powCount += powTmp

            if i == (N_sources - batch_size):
                targetPower += np.sum(np.array(powCount).reshape(size_i, totFace), axis=0)
                #loopFlag = 1


            #take out the garbage
            del sx_batch,tx_batch, sy_batch,ty_batch, sz_batch,tz_batch, tp_batch, sA
            del origin, target, direction, ray, si
            del primI, pow, piA, mask, correctHit, powTmp

            print("Batch Loop Time: {:f}[s]".format(time.time() - tLoop))


        self.targetPower = targetPower[:self.Nj]
        print('Mitsuba Calc Time: '+str(time.time()-t0))
        
        return


    def calculatePowerTransferMitsubaNumpy(self, mode=None, mitsubaMode='cuda', fType='ply', batch_size=100):
        """
        Maps power between sources and targets (ROI PFCs).  Uses Mitsuba3 to
        perform ray tracing.  Mitsuba3 can be optimized for CPU or GPU.

        Code largely developed by A. Rosenthal (CFS)
        Adapted to HEAT by T. Looby

        Mitsuba:
        @software{jakob2022mitsuba3,
            title = {Mitsuba 3 renderer},
            author = {Wenzel Jakob and Sébastien Speierer and Nicolas Roussel and Merlin Nimier-David and Delio Vicini and Tizian Zeltner and Baptiste Nicolet and Miguel Crespo and Vincent Leroy and Ziyi Zhang},
            note = {https://mitsuba-renderer.org},
            version = {3.0.1},
            year = 2022,
        }

        """
        import drjit as dr
        import mitsuba as mi

        t0 = time.time()
        powerFrac = np.zeros((self.Ni,self.Nj))
        Psum = np.zeros((self.Nj))
        self.hullPower = np.zeros((self.Ni))
        self.powerFrac = np.zeros((self.Ni))

        print("Building radiation scene...")    

        # Initialize Mitsuba
        if mitsubaMode == 'cuda':
            mi.set_variant('cuda_ad_rgb')
        else:
            mi.set_variant('llvm_ad_rgb')

        scene = {'type': 'scene'}
        scene['path_int'] = {
            'type': 'path',
            'max_depth': 1  # Maximum number of bounces
        }

        name = 'allMesh' #name of entire mesh we saved

        scene[name] = {
            'type': fType,
            'filename': self.meshFile,
            'to_world': mi.ScalarTransform4f.scale(0.001),  # Convert from mm to m
            'face_normals': True,  # This prevents smoothing of sharp-corners by discarding surface-normals. Useful for engineering CAD.
            'sensor': {
                'type': 'irradiancemeter',
                'film': {
                    'type': 'hdrfilm',
                    'pixel_format': 'luminance',
                    'filter': {'type': 'box'},
                    'width': 1,
                    'height': 1,
                }
            }
        }

        t = time.time()
        print('Loading Scene....')
        scene = mi.load_dict(scene)
        print('Time to load scene: ' + str(time.time() - t))
        print('\n')

        t = time.time()
        # Convert the numpy arrays to Dr.Jit dynamic arrays
        tx = dr.cuda.Float(self.sources[:,0]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,0])
        ty = dr.cuda.Float(self.sources[:,1]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,1])
        tz = dr.cuda.Float(self.sources[:,2]) if mode == 'cuda' else dr.llvm.Float(self.sources[:,2])
        tp = dr.cuda.Float(self.sourcePower) if mode == 'cuda' else dr.llvm.Float(self.sourcePower)

        # Grab all the shape ids to correctly identify the sensor
        idL = [x.id() for x in scene.shapes()]

        print('Time to build sources: ' + str(time.time() - t))

        t = time.time()
        # Find the sensor
        sensInd = idL.index(name)
        sensS = scene.shapes()[sensInd]

        params = mi.traverse(sensS)
        verts = params['vertex_positions']
        norms = params['vertex_normals']

        # Find all the face and areas
        totFace = sensS.face_count()
        if mitsubaMode == 'cuda':
            indF = sensS.face_indices(dr.arange(dr.cuda.ad.UInt, totFace))
        else:
            indF = sensS.face_indices(dr.arange(dr.llvm.ad.UInt, totFace))

        # Label each vertex with its index
        vertNum = int(len(verts) / 3)
        vertA = np.array(verts)

        vertAR = vertA.reshape((vertNum, 3))
        indFA = np.array(indF)

        # Find the vertices of each face
        faceV = np.take(vertAR, indFA.flatten(), axis=0)

        # Find the area of each face and the mean point of each face
        faceV1 = faceV[::3, :]
        faceV2 = faceV[1::3, :]
        faceV3 = faceV[2::3, :]

        # Find the area of each face
        cross1 = faceV2 - faceV1
        cross2 = faceV3 - faceV1

        cross = np.cross(cross1, cross2)
        area = np.linalg.norm(cross, axis=1) / 2

        #find the mean point of each face
        #again can probably be sped up
        center = np.asarray([x.mean(0) for x in np.split(faceV,totFace,axis = 0)])

        #dynamically allocate batch size based upon available memory
        if batch_size == "auto":
            batch_size = int(2**32 / totFace  * 0.5)

        print("Calculated batch size of: {:d}".format(batch_size))
        
        N_sources = len(tx)
        N_targets = len(center)
        targetPower = np.zeros((N_targets))
        print('Time to prepare mesh: ' + str(time.time() - t))

        t0 = time.time()
        for i in range(0, N_sources, batch_size):
            tLoop = time.time()
            #dynamically adjust size of batch depending on how many
            #sources are left to calculate
            if (N_sources - (batch_size+i)) < 0:
                size_i = np.abs(N_sources - i)
            else:
                size_i = batch_size     
            
            print("\n=== Running Sources: {:d} - {:d} ===".format(i, i+size_i))
            t = time.time()
            dummy = np.ones((size_i))

            #take all the x values from the sensor and pair with all the x values from the power distribution
            if mitsubaMode == 'cuda':
                sx_batch,tx_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,0]),dr.cuda.ad.Float(tx[i:i+size_i]),indexing = 'xy') #indexing ij makes it so that the first numpower points of sx are all the same (tx[0])
                sy_batch,ty_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,1]),dr.cuda.ad.Float(ty[i:i+size_i]),indexing = 'xy')
                sz_batch,tz_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,2]),dr.cuda.ad.Float(tz[i:i+size_i]),indexing = 'xy')
                _,tp_batch = dr.meshgrid(dr.cuda.ad.Float(center[:,0]),dr.cuda.ad.Float(tp[i:i+size_i]),indexing = 'xy') #make the power meshed in the same way
                sA,_ = dr.meshgrid(dr.cuda.ad.Float(area),dr.cuda.ad.Float(dummy),indexing = 'xy')
            else:
                sx_batch,tx_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,0]),dr.llvm.ad.Float(tx[i:i+size_i]),indexing = 'xy') #indexing ij makes it so that the first numpower points of sx are all the same (tx[0])
                sy_batch,ty_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,1]),dr.llvm.ad.Float(ty[i:i+size_i]),indexing = 'xy')
                sz_batch,tz_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,2]),dr.llvm.ad.Float(tz[i:i+size_i]),indexing = 'xy')
                _,tp_batch = dr.meshgrid(dr.llvm.ad.Float(center[:,0]),dr.llvm.ad.Float(tp[i:i+size_i]),indexing = 'xy') #make the power meshed in the same way
                sA,_ = dr.meshgrid(dr.llvm.ad.Float(area),dr.llvm.ad.Float(dummy),indexing = 'xy')

            print('Time to meshgrid: ' + str(time.time() - t))

            #define the ray origins which start at the HEAT grid points
            origin = mi.Point3f(tx_batch,ty_batch,tz_batch)
            #and the target points which are the sensor faces
            target = mi.Vector3f(sx_batch,sy_batch,sz_batch)
            direction = dr.normalize(target-origin)
            ray = mi.Ray3f(o=origin, d=direction)

            t = time.time()
            si = scene.ray_intersect(ray)
            print('Time to run intersection: '+str(time.time()-t))
            
            t = time.time()
            #gives the index of the primitive triangle that was hit
            primI = si.prim_index
            pow = dr.abs(dr.dot(dr.normalize(si.n),direction)) * sA * tp_batch/(4*dr.pi*dr.power(si.t,2))
            
            lenX = size_i
            lenY = totFace
            #mitsuba arrays converted to numpy and reshaped to size_i X totFaces
            primArray = np.array(primI).reshape(lenX, lenY)
            distArray = np.array(si.t).reshape(lenX, lenY)
            powArray = np.array(pow).reshape(lenX, lenY)

            infMask = np.where(distArray==np.inf)
            powArray[infMask] == 0.0
            print('Time to build numpy: '+str(time.time()-t))

            t = time.time()
            #boolean where true if primI equals its j index
            #j index corresponds to target mesh elements, so this means the ray
            #traced from source i to target j terminated on target j (and not
            #on a different target face)
            mask = primArray == np.arange(totFace)[np.newaxis,:]

            #assign all the shadowed faces 0 power contribution
            shadow = np.where(mask==False)
            powArray[shadow] = 0.0

            sumP = np.sum(powArray, axis=0)
            targetPower += sumP
            print('Time to build shadows: '+str(time.time()-t))
            print("Batch Loop Time: {:f}[s]".format(time.time() - tLoop))

        self.targetPower = targetPower[:self.Nj]
        print('Mitsuba Calc Time: '+str(time.time()-t0))
        
        return


    def calculatePowerTransferOpen3D(self, mode=None, powFracSave=False):
        """
        Maps power between sources and targets (ROI PFCs).  Uses Open3D to
        perform ray tracing.  Open3D can be optimized for CPU or GPU.

        Uses Open3D to accelerate the calculation:
        Zhou, Qian-Yi, Jaesik Park, and Vladlen Koltun. "Open3D: A modern
        library for 3D data processing." arXiv preprint arXiv:1801.09847 (2018).
        """
        import psutil
        t0 = time.time()

        Psum = np.zeros((self.Nj))
        self.hullPower = np.zeros((self.Ni))

        print("Building radiation scene...")
        #build mesh and tensors for open3d
        mesh = o3d.io.read_triangle_mesh(self.meshFile)
        mesh.compute_vertex_normals()
        mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
        scene = o3d.t.geometry.RaycastingScene()
        mesh_id = scene.add_triangles(mesh)
        print("Scene building took {:f} [s]\n".format(time.time() - t0))

        #check how much memory we have
        #memNeeded = np.dtype(np.float64).itemsize * self.Ni * self.Nj
        #print("Required memory for this calculation: {:f} GB".format(memNeeded / (1024**3)))
        memAvail = psutil.virtual_memory().available
        print("Available memory: {:f} GB".format(memAvail / (1024**3)))  

        #setup netcdf (this is slow)
        if powFracSave:
            print("Saving power fractions.\n")
            from netCDF4 import Dataset
            # Create a powerFrac file on disk
            # Initialize the NetCDF file
            with Dataset(self.powFracFile, 'w', format='NETCDF4') as ncfile:
                # Define dimensions
                ncfile.createDimension('Ni', self.Ni)
                ncfile.createDimension('Nj', self.Nj)
                # Create a variable to store the data
                data_var = ncfile.createVariable('powerFrac', 'f4', ('Ni', 'Nj'))

        for i in range(self.Ni):
            if i%1000 == 0:
                if i==0:
                    t1 = time.time()
                process = psutil.Process()
                mem_info = process.memory_info()
                usage = mem_info.rss / (1024 * 1024)
                print("Source point {:d}.  1k time: {:f}. Memory Used: {:f} MB".format(i, time.time() - t1, usage))
                t1 = time.time()

            #r_ij = np.zeros((self.Nj,3))
            r_ij = self.targetCtrs - self.sources[i]
            #r_ij *= 1000.0
            rMag = np.linalg.norm(r_ij, axis=1)
            rNorm = r_ij / rMag.reshape((-1,1))
            rdotn = np.sum(rNorm*self.targetNorms, axis=1)
            q1 = np.tile(self.sources[i]*1000.0,(self.Nj,1))

            #calculate intersections
            rays = o3d.core.Tensor([np.hstack([q1,rNorm])],dtype=o3d.core.Dtype.Float32)
            hits = scene.cast_rays(rays)

            #convert open3d CPU tensors back to numpy
            hitMap = hits['primitive_ids'][0].numpy()
            distMap = hits['t_hit'][0].numpy()
            condition = hitMap == np.arange(self.Nj)
            powFrac = np.abs(rdotn)*self.targetAreas/(4*np.pi*rMag**2)

            if powFracSave == True:
                # Open the NetCDF file in append mode and write the column
                # this is slow!
                with Dataset(self.powFracFile, 'a', format='NETCDF4') as ncfile:
                    data_var = ncfile.variables['powerFrac']
                    # Write the column to the file. This operation writes directly to disk.
                    data_var[i, :] = powFrac

            Psum += condition * self.sourcePower[i] * powFrac

#            #compute convex hull on unit sphere around point i and calculate
#            #power balance via ratio of hull area to sphere area
#            #calculate spherical coordinates for each target point, with i as (0,0,0)
#            #note that we perform a coordinate permutation on rNorm, to prevent
#            #angle wrap cases: (x,y,z) => (z,x,y)
#            #If you are getting abnormally convex hull power, try a different permutation
#            #first transform to 2D (phi,theta) => (x,y) using Jacobians
#            theta = np.arccos( (rNorm[:,0]) ) #z after permutation
#            phi = np.arctan2( (rNorm[:,2]), (rNorm[:,1])  ) # y,x after permutation
#            #map (phi,theta) to cartesian plane using Jacobian
#            points = np.vstack([phi, -np.cos(theta)]).T
#            #calculate the convex hull in 2D
#            hull = ConvexHull(points)
#            #hull area for 2D points is volume, scale power to full solid angle
#            self.hullPower[i] = hull.volume / (4*np.pi) * self.sourcePower[i]

        self.pdotn = rdotn
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

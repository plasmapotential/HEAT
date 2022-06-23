#RADClass.py
#Description:   Base HEAT Radiated Power
#Engineer:      T Looby
#Date:          20220613
import numpy as np
import pandas as pd


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
                            ]
        return

    def setTypes(self):
        """
        Nothing to do for this class
        """
        return

    def read2DpcSourceFile(self, file):
        """
        reads a comma delimited file that describes the radiated power on an
        R,Z grid. This point cloud (pc) should represent the axisymmetric
        plasma radiated power across a single poloidal plane (ie sliced at phi).
        The power defined at each RZ point is the radiated photon power.

        It is assumed that the RZ points are located at the centroids of mesh
        elements
        """
        self.PC2D = pd.read_csv(file, header=0, names=['R','Z','P'])
        return

    def paraview2D(self):
        """
        """
        return

    def paraview3D(self):
        """
        """
        return

    def create3DpcFrom2D(self, Ntor):
        """
        creates a 3D point cloud from a 2D point cloud.
        assumes the 2D point cloud is in a poloidal plane (ie sliced at a phi)
        and extrudes the 2D point cloud around the machine in toroidal direction.
        This extrusion results in 3D voxels (3D version of a 2D pixel), where
        each point is the centroid of the voxel.  User specifies resolution
        of voxels in toroidal direction.

        Ntor:   number of discrete voxels in toroidal direction (360 would result
                in voxels of 1 degree toroidal extent)

        """
        return



    def source2TargetVectors(self, sourceXYZ, targetXYZ):
        """

        """
        return

    def vectors2Distances(self, vecs):
        """
        """
        return

    def normalizeVectors(self, vecs):
        """
        """
        return

    def backfaceCulling(self):
        """
        """
        return

    def findShadows(self):
        """
        """
        return

    def calculatePowerTransfer(self):
        """
        """
        return

    def calculateBRDF(self, material):
        """
        constructs the Bidirectional Reflection Distribution Function for a user
        defined material.  Once constructed, user can evaluate specific BRDF
        values for input and output angles.
        """
        return

    def mapPowerParallel(self, i):
        """
        maps power from source to target, looking for intersections along the way
        """
        #source to target vectors
        p_ij = self.source[i] - self.targetCtrs
        pMag = np.linalg.norm(r_ij, axis=1)
        pNorm = p_ij / pMag[:, np.newaxis]
        pdotn = np.sum(pNorm*self.targetNorms, axis=1)

        #backface culling
        backFaces = np.where( pdotn > 0 )[0]
        self.shadowMask[backFaces] = 1.0
        use0 = np.where(self.shadowMask == 0)[0]

        #should we add backface culling for intersects here?!

        #run in matrix or loop mode depending upon available memory
        #availableRAM = psutil.virtual_memory().available #bytes
        #neededRAM = len(self.sources) * len(self.targets) * len(self.intersects) * 112 #bytes

        q1 = self.source[i]
        q2 = self.targetCtrs[use0]
        Nt = len(self.intersects)

        #loop thru targets checking if shadowed
        for j in len(q2):
            #Perform Intersection Test
            q13D = np.repeat(self.q1[np.newaxis], Nt, axis=0)
            q23D = np.repeat(self.q2[j,np.newaxis], Nt, axis=0)
            sign1 = np.sign(self.signedVolume2(q13D,self.p1,self.p2,self.p3))
            sign2 = np.sign(self.signedVolume2(q23D,self.p1,self.p2,self.p3))
            sign3 = np.sign(self.signedVolume2(q13D,q23D,self.p1,self.p2))
            sign4 = np.sign(self.signedVolume2(q13D,q23D,self.p2,self.p3))
            sign5 = np.sign(self.signedVolume2(q13D,q23D,self.p3,self.p1))
            test1 = (sign1 != sign2)
            test2 = np.logical_and(sign3==sign4,sign3==sign5)

            #we intersected something, so face is shadowed
            if np.sum(np.logical_and(test1,test2)) > 0:
                P = 0.0
            #we didn't intersect. assign power
            else:
                P = pdotn[j]*self.Prad[i]*self.targetArea[j]/(4*np.pi*pMag**2)

        return P

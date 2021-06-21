import numpy as np
import pandas as pd
import plotly.graph_objects as go
import toolsClass
import multiprocessing
import time
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from tqdm.contrib.concurrent import process_map
tools = toolsClass.tools()

import logging
log = logging.getLogger(__name__)


class GYRO:

    def __init__(self, rootDir, dataPath):
        """
        rootDir is root location of python modules (where dashGUI.py lives)
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

        """


        self.allowed_vars = [
                    'N_gyroSteps',
                    'gyroDeg',
                    'gyroT_eV',
                    'N_vSlice',
                    'N_vPhase',
                    'N_gyroPhase',
                    'ionMassAMU',
                    'vMode',
                    'ionFrac'
                            ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        self.N_gyroSteps = int(self.N_gyroSteps)
        self.gyroDeg = int(self.gyroDeg)
        self.gyroT_eV = float(self.gyroT_eV)
        self.N_vSlice = int(self.N_vSlice)
        self.N_vPhase = int(self.N_vPhase)
        self.N_gyroPhase = int(self.N_gyroPhase)
        self.ionMassAMU = float(self.ionMassAMU)
        self.vMode = self.vMode
        self.ionFrac = float(self.ionFrac)
        return

    def setupConstants(self, ionMassAMU=2.014):
        """
        Sets up constants

        default mass is deuterium 2.014 MeV/c^2
        """
        #unit conversions
        self.kg2eV = 5.609e35 #1kg = 5.609e35 eV/c^2
        self.eV2K = 1.160e4 #1ev=1.160e4 K
        #constants
        self.AMU = 931.494e6 #ev/c^2
        self.kB = 8.617e-5 #ev/K
        self.e = 1.602e-19 # C
        self.c = 299792458 #m/s

        self.mass_eV = ionMassAMU * self.AMU
        self.Z=1 #assuming isotopes of hydrogen here

        return

    def temp2thermalVelocity(self, T_eV):
        """
        Calculates thermal velocity from a temperature

        T_eV is temperature in eV

        can also be found with: d/dv( v*f(v) ) = 0
        """
        return np.sqrt(T_eV/(self.mass_eV/self.c**2))


    def setupFreqs(self, B):
        """
        Calculates frequencies, periods, that are dependent upon B
        These definitions follow Freidberg Section 7.7.

        B is magnetic field magnitude

        """
        self.omegaGyro = self.Z * self.e * B / (self.mass_eV / self.kg2eV)
        if np.isscalar(self.omegaGyro):
            self.omegaGyro = np.array([self.omegaGyro])

        self.fGyro = self.omegaGyro/(2*np.pi)
        self.TGyro = 1.0/self.fGyro
        return

    def setupRadius(self, vPerp):
        """
        calculates gyro radius.

        rGyro has a column for each MC run (N_MC columns), and a
        row for each point on the PFC (N_pts), so it is a matrix
        of shape: N_pts X N_MC
        """
        N_pts = len(self.omegaGyro)
        #get number of vPerps
        if np.isscalar(vPerp):
            vPerp = np.array([vPerp])
            N_MC = 1
        else:
            N_MC = len(vPerp)

        self.rGyro = np.zeros((N_pts,N_MC))

        for i in range(N_MC):
            self.rGyro[:,i] = vPerp[i] / self.omegaGyro

        return

    def setupVelocities(self):
        """
        sets up velocities based upon vMode input from GUI

        len(self.t1) is number of points in divertor we are calculating HF on
        """
        #get velocity space phase angles
        self.uniformVelPhaseAngle()

        if self.vMode == 'single':
            print("Gyro orbit calculation from single plasma temperature")
            log.info("Gyro orbit calculation from single plasma temperature")
            self.T0 = np.ones((len(self.sourceCenters)))*self.gyroT_eV
            #get average velocity for each temperature point
            self.vThermal = self.temp2thermalVelocity(self.T0)
            #set upper bound of v*f(v) (note that this cuts off high energy particles)
            self.vMax = 4 * self.vThermal
            #get 100 points to initialize functional form of f(v) (note this is a 2D matrix cause vMax is 2D)
            self.vScan = np.linspace(0,self.vMax,100).T
            #get velocity slices for each T0
            self.pullEqualProbabilityVelocities()

        else:
            #TO ADD THIS YOU WILL NEED TO PASS IN XYZ COORDINATES OF CTRS AND INTERPOLATE
            print("3D plasma temperature interpolation from file not yet supported.  Run gyro orbits in single mode")
            log.info("3D plasma temperature interpolation from file not yet supported.  Run gyro orbits in single mode")

        return

    def pullEqualProbabilityVelocities(self):
        """
        creates vSlices: array of velocities indexed to match T0 array (or PFC.centers)

        each vSlice is positioned at a place in the PDF so it has an equal probability
        of occuring.  ie the area under the PDF curve between each vSlice is equal.
        """
        self.vSlices = np.ones((len(self.T0),self.N_vSlice))*np.nan
        self.energySlices = np.zeros((len(self.T0),self.N_vSlice))
        for i in range(len(self.T0)):
            #get velocity range for this T0
            v = self.vScan[i,:]
            #generate the (here maxwellian) PDF
            pdf = lambda x: (self.mass_eV/self.c**2) / (self.T0[i]) * np.exp(-(self.mass_eV/self.c**2 * x**2) / (2*self.T0[i]) )
            #speed pdf
            v_pdf = v * pdf(v)
            #generate the CDF
            v_cdf = np.cumsum(v_pdf[1:])*np.diff(v)
            v_cdf = np.insert(v_cdf, 0, 0)
            #create bspline interpolators for the cdf and cdf inverse
            inverseCDF = interp1d(v_cdf, v, kind='cubic')
            forwardCDF = interp1d(v, v_cdf, kind='cubic')
            #calculate N_vSlice velocities for each pdf each with equal area (probability)
            cdfMax = v_cdf[-1]
            cdfMin = v_cdf[0]
            sliceWidth = cdfMax / (self.N_vSlice+1)
            cdfSlices = np.linspace(0,1,self.N_vSlice+2)[1:-1]
            self.vSlices[i,:] = inverseCDF(cdfSlices)

            #Now find energies that correspond to these vSlices
            #we integrate: v**2 * f(v)
            #energy pdf (missing 1/2*mass but that gets divided out later anyways )
            energy = lambda x: x**2 * pdf(x)
            #if there is only 1 vSlice integrate entire pdf
            if len(self.vSlices[i]==1):
                vLo = 0.0
                vHi = self.vMax[i]
                self.energySlices[i] = integrate.quad(energy, vLo, vHi)[0]
            #if there are multiple vSlices use them as integral bounds
            else:
                for j in range(len(self.vSlices[i])-1):
                    if j==0:
                        vLo = 0.0
                        vHi = self.vMax[i,0]
                    elif j==len(self.vSlices[i])-2:
                        vLo = self.vSlices[i,-1]
                        vHi = self.vMax[i]
                    else:
                        vLo = self.vSlices[i,j-1]
                        vHi = self.vSlices[i,j]

                    self.energySlices[i,j] = integrate.quad(energy, vLo, vHi)[0]

        print("Found N_vPhase velocities of equal probability")
        log.info("Found N_vPhase velocities of equal probability")
        return

    def maxwellian(self, x):
        """
        returns a maxwellian distribution with x as independent variable.  x
        must either be scalar or len(T0)

        Uses mass, T0, c, from class
        """
        pdf = (self.mass_eV/self.c**2) / (self.T0) * np.exp(-(self.mass_eV/self.c**2 * x**2) / (2*self.T0) )
        return pdf

    def uniformGyroPhaseAngle(self):
        """
        Uniform sampling of a uniform distribution between 0 and 2pi

        returns angles in radians
        """
        self.gyroPhases = np.linspace(0,2*np.pi,self.N_gyroPhase+1)[:-1]
        return

    def uniformVelPhaseAngle(self):
        """
        Sampling of a uniform distribution between 0 and pi/2 (only forward velocities)
        vPerp is x-axis of velocity space
        vParallel is y-axis of velocity space

        returns angles in radians
        """
        self.vPhases = np.linspace(0.0,np.pi/2,self.N_vPhase+2)[1:-1]
        return


    def singleGyroTrace(self,vPerp,vParallel,gyroPhase,N_gyroSteps,
                        BtraceXYZ,controlfilePath,TGyro,rGyro,omegaGyro,
                        verbose=True):
        """
        Calculates the gyro-Orbit path and saves to .csv and .vtk

        vPerp and vParallel [m/s] are in velocities
        gyroPhase [degrees] is initial orbit phase angle
        N_gyroSteps is number of discrete line segments per gyro period
        BtraceXYZ is the points of the Bfield trace that we will gyrate about
        """
        print("Calculating gyro trace...")
        #Loop thru B field trace while tracing gyro orbit
        helixTrace = None
        for i in range(len(BtraceXYZ)-1):
            #points in this iteration
            p0 = BtraceXYZ[i,:]
            p1 = BtraceXYZ[i+1,:]
            #vector
            delP = p1 - p0
            #magnitude
            magP = np.sqrt(delP[0]**2 + delP[1]**2 + delP[2]**2)
            #time it takes to transit line segment
            delta_t = magP / (vParallel)
            #Number of steps in line segment
            Tsample = self.TGyro / N_gyroSteps
            Nsteps = int(delta_t / Tsample)
            #length (in time) along guiding center
            t = np.linspace(0,delta_t,Nsteps+1)
            #guiding center location
            xGC = np.linspace(p0[0],p1[0],Nsteps+1)
            yGC = np.linspace(p0[1],p1[1],Nsteps+1)
            zGC = np.linspace(p0[2],p1[2],Nsteps+1)
            # construct orthogonal system for coordinate transformation
            w = delP
            if np.all(w==[0,0,1]):
                u = np.cross(w,[0,1,0]) #prevent failure if bhat = [0,0,1]
            else:
                u = np.cross(w,[0,0,1]) #this would fail if bhat = [0,0,1] (rare)
            v = np.cross(w,u)
            #normalize
            u = u / np.sqrt(u.dot(u))
            v = v / np.sqrt(v.dot(v))
            w = w / np.sqrt(w.dot(w))
            xfm = np.vstack([u,v,w]).T
            #get helix path along (proxy) z axis reference frame
            x_helix = self.rGyro*np.cos(self.omegaGyro*t + gyroPhase)
            y_helix = self.rGyro*np.sin(self.omegaGyro*t + gyroPhase)
            z_helix = np.zeros((len(t)))
            #perform rotation to field line reference frame
            helix = np.vstack([x_helix,y_helix,z_helix]).T
            helix_rot = np.zeros((len(helix),3))
            for j,coord in enumerate(helix):
                helix_rot[j,:] = helix[j,0]*u + helix[j,1]*v + helix[j,2]*w
            #perform translation to field line reference frame
            helix_rot[:,0] += xGC
            helix_rot[:,1] += yGC
            helix_rot[:,2] += zGC
            #update gyroPhase variable so next iteration starts here
            gyroPhase = self.omegaGyro*t[-1] + gyroPhase
            #append to helix trace
            if helixTrace is None:
                helixTrace = helix_rot
            else:
                helixTrace = np.vstack([helixTrace,helix_rot])

        helixTrace*=1000.0 #scale for ParaView
        print("Saving data to CSV and VTK formats")
        #save data to csv format
        head = 'X[mm],Y[mm],Z[mm]'
        np.savetxt(controlfilePath+'helix.csv', helixTrace, delimiter=',', header=head)
        #save data to vtk format
        tools.createVTKOutput(controlfilePath+'helix.csv', 'trace', 'Gyro_trace')

        if verbose==True:
            print("V_perp = {:f} [m/s]".format(vPerp))
            print("V_parallel = {:f} [m/s]".format(vParallel))
            print("Cyclotron Freq = {:f} [rad/s]".format(self.omegaGyro[0]))
            print("Cyclotron Freq = {:f} [Hz]".format(self.fGyro[0]))
            print("Gyro Radius = {:f} [m]".format(self.rGyro[0][0]))
            print("Number of gyro points = {:f}".format(len(helixTrace)))
            print("Longitudinal dist between gyro points = {:f} [m]".format(magP/float(Nsteps)))
            print("Each line segment length ~ {:f} [m]".format(magP))
        return

    def gyroTraceParallel(self, i):
        """
        parallelized gyro trace.  called by multiprocessing.pool.map()

        i is index of parallel run from multiprocessing

        p0 and p1 are start/end points of field line trace

        writes helical trace to self.helixTrace[i] in 2D matrix format:
            columns = X,Y,Z
            rows = steps up helical trace

        also updates self.lastPhase for use in next iteration step
        """
        t1 = time.time()
        #vector
        delP = self.p1[i] - self.p0[i]
        #magnitude
        magP = np.sqrt(delP[0]**2 + delP[1]**2 + delP[2]**2)
        #time it takes to transit line segment
        delta_t = magP / (self.vParallelMC[i])
        #Number of steps in line segment
        Tsample = self.TGyro[i] / self.N_gyroSteps
        Nsteps = int(delta_t / Tsample)

        #length (in time) along guiding center
        t = np.linspace(0,delta_t,Nsteps+1)
        #guiding center location
        xGC = np.linspace(self.p0[i,0],self.p1[i,0],Nsteps+1)
        yGC = np.linspace(self.p0[i,1],self.p1[i,1],Nsteps+1)
        zGC = np.linspace(self.p0[i,2],self.p1[i,2],Nsteps+1)
        arrGC = np.vstack([xGC,yGC,zGC]).T
        # construct orthogonal system for coordinate transformation
        w = delP
        if np.all(w==[0,0,1]):
            u = np.cross(w,[0,1,0]) #prevent failure if bhat = [0,0,1]
        else:
            u = np.cross(w,[0,0,1]) #this would fail if bhat = [0,0,1] (rare)
        v = np.cross(w,u)
        #normalize
        u = u / np.sqrt(u.dot(u))
        v = v / np.sqrt(v.dot(v))
        w = w / np.sqrt(w.dot(w))
        xfm = np.vstack([u,v,w]).T
        #get helix path along (proxy) z axis reference frame
        x_helix = self.rGyroMC[i]*np.cos(self.omegaGyro[i]*t + self.lastPhase[i])
        y_helix = self.rGyroMC[i]*np.sin(self.omegaGyro[i]*t + self.lastPhase[i])
        z_helix = np.zeros((len(t)))
        #perform rotation to field line reference frame
        helix = np.vstack([x_helix,y_helix,z_helix]).T
        helix_rot = np.zeros((len(helix),3))
        for j,coord in enumerate(helix):
            helix_rot[j,:] = helix[j,0]*u + helix[j,1]*v + helix[j,2]*w
        #perform translation to field line reference frame
        helix_rot[:,0] += xGC
        helix_rot[:,1] += yGC
        helix_rot[:,2] += zGC

        #shift entire helix to ensure we capture intersections in p0 plane
        helix_rot[:,0] += w[0]*0.0003
        helix_rot[:,1] += w[1]*0.0003
        helix_rot[:,2] += w[2]*0.0003

        #update gyroPhase variable so next iteration starts here
        lastPhase = self.omegaGyro[i]*t[-1] + self.lastPhase[i]

        #print the trace for a specific index
        if self.traceIndex2 is not None:
            if self.traceIndex2 == i:
                #print("Saving Index data to CSV and VTK formats")
                #save data to csv format
                head = 'X[mm],Y[mm],Z[mm]'
                np.savetxt(self.controlfilePath+'helix{:d}.csv'.format(self.N_GCdeg), helix_rot*1000.0, delimiter=',', header=head)
                #save data to vtk format
                tools.createVTKOutput(self.controlfilePath+'helix{:d}.csv'.format(self.N_GCdeg),
                                        'trace', 'traceHelix{:d}'.format(self.N_GCdeg),verbose=False)
                #guiding center
                #np.savetxt(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg), arrGC*1000.0, delimiter=',', header=head)
                #save data to vtk format
                #tools.createVTKOutput(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg),
                #                        'trace', 'traceGC{:d}'.format(self.N_GCdeg),verbose=False)


        #=== intersection checking ===
        q1 = helix_rot[:-1,:]
        q2 = helix_rot[1:,:]

        #Filter by psi
        if self.psiFilterSwitch == True:
            use = np.where(self.psiMask[:,self.indexMap[i]] == 1)[0]
        else:
            use = np.arange(len(self.PFC_t1))
        Nt = len(use)

        #using full array (no for loop)
        RAMlimit = False
        if RAMlimit == False:
            q13D = np.repeat(q1[:,np.newaxis], Nt, axis=1)
            q23D = np.repeat(q2[:,np.newaxis], Nt, axis=1)
            sign1 = np.sign(tools.signedVolume2(q13D,self.PFC_t1,self.PFC_t2,self.PFC_t3,ax=2))
            sign2 = np.sign(tools.signedVolume2(q23D,self.PFC_t1,self.PFC_t2,self.PFC_t3,ax=2))
            sign3 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t1,self.PFC_t2,ax=2))
            sign4 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t2,self.PFC_t3,ax=2))
            sign5 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t3,self.PFC_t1,ax=2))
            test1 = (sign1 != sign2)
            test2 = np.logical_and(sign3==sign4,sign3==sign5)
            test3 = np.logical_and(test1,test2)
            loc = np.where(np.logical_and(test1,test2))

            #result=1 if we intersected, otherwise NaN
            if np.sum(np.logical_and(test1,test2)) > 0:
                #only take first index (ie first intersection location)
                loc = np.where(np.logical_and(test1,test2))
                loc = loc[0][0],loc[1][0]
                index = use[loc[1]]
                #if self.traceIndex2 == i:
                #    print("TEST!!!")
                #    print(np.where(np.logical_and(test1,test2))[0])
                #    print(index)
                vec = (q2[loc[0]] - q1[loc[0]]) / np.linalg.norm(q2[loc[0]]-q1[loc[0]])
                hdotn = np.dot(self.intersectNorms[index],vec)
            else:
                index = np.NaN
                hdotn = np.NaN
        #using loop
        else:
            #loop thru each step of helical path looking for intersections
            for j in range(len(helix_rot)-1):
                #Perform Intersection Test
                q13D = np.repeat(q1[j,np.newaxis], Nt, axis=0)
                q23D = np.repeat(q2[j,np.newaxis], Nt, axis=0)
                sign1 = np.sign(tools.signedVolume2(q13D,self.PFC_t1[use],self.PFC_t2[use],self.PFC_t3[use]))
                sign2 = np.sign(tools.signedVolume2(q23D,self.PFC_t1[use],self.PFC_t2[use],self.PFC_t3[use]))
                sign3 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t1[use],self.PFC_t2[use]))
                sign4 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t2[use],self.PFC_t3[use]))
                sign5 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t3[use],self.PFC_t1[use]))
                test1 = (sign1 != sign2)
                test2 = np.logical_and(sign3==sign4,sign3==sign5)

                #result=1 if we intersected, otherwise NaN
                if np.sum(np.logical_and(test1,test2)) > 0:
                    #only take first index (ie first intersection location)
                    #YOU SHOULD CHECK THIS TO MAKE SURE [0][0] is the first face along field line
                    index = use[ np.where(np.logical_and(test1,test2))[0][0] ]
                    #if self.traceIndex2 == i:
                    #    print("TEST!!!")
                    #    print(np.where(np.logical_and(test1,test2))[0])
                    #    print(index)
                    vec = (q2[j] - q1[j]) / np.linalg.norm(q2[j]-q1[j])
                    hdotn = np.dot(self.intersectNorms[index],vec)
                    break
                else:
                    index = np.NaN
                    hdotn = np.NaN

        print("TIME: {:f}".format(time.time() - t1))
        return lastPhase, index, hdotn

    def multipleGyroTrace(self):
        """
        Calculates the helical path for multiple points,
        each with different gyroRadii, using multiprocessing

        Btrace is one step of a field line trace for each point (MAFOT structure output)
        phase is phase angle (updated from last trace step)

        updates lastPhase variable and helixTrace
        """
        #magnetic field trace
        self.helixTrace = [None] * len(self.p0)
        N = len(self.p1)
        #Prepare helical trace across multiple cores
        Ncores = multiprocessing.cpu_count() -2 #reserve 2 cores for overhead
        print('Initializing parallel helix trace across {:d} cores'.format(Ncores))
        log.info('Initializing parallel helix trace across {:d} cores'.format(Ncores))
        #each worker receives a single start and end point (p0 and p1),
        #corresponding to one trace from the MAFOT structure output.
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
#OLD METHOD
#        #Do this try clause to kill any zombie threads that don't terminate
#        try:
#            pool = multiprocessing.Pool(Ncores)
#            output = np.asarray(pool.map(self.gyroTraceParallel, np.arange(N)))
#            self.lastPhase = output[:,0]
#            intersectRecord = output[:,1]
#        finally:
#            pool.close()
#            pool.join()
#            del pool

        #multiprocessing with status bar (equiv to multiprocessing.Pool.map())
        print("Multiprocessing gyro trace:")
        output = process_map(self.gyroTraceParallel, range(N), max_workers=Ncores, chunksize=1)
        output = np.asarray(output)

        intersectRecord = output[:,1]
        use = np.where(np.isnan(intersectRecord)==True)[0]
        self.lastPhase = output[:,0][use]
        hdotn = output[:,2]
        print('Parallel helix trace complete')
        log.info('Parallel helix trace complete')
        return intersectRecord, hdotn


    def writeIntersectRecord(self, gyroPhase, vPhase, vSlice, file, idx):
        """
        writes intersectRecord to CSV file

        1 file for each gyroPhase, vPhase, vSlice
        """
        print("Writing out intersectRecords")
        log.info("Writing out intersectRecords")
        #write the velocities for this run to a comment in file
        f = open(file, 'w')
        f.write('# gyroPhase: {:f} [radians]\n'.format(self.gyroPhases[gyroPhase]))
        f.write('# vPhase: {:f} [radians]\n'.format(self.vPhases[vPhase]))
        rec = self.intersectRecord[gyroPhase,vPhase,vSlice,:][idx]
        data = {
                'face': pd.Series(np.arange(len(rec))),
                'intersectFace': pd.Series(rec),
                'vPerp[m/s]': pd.Series(self.vPerpMC),
                'vParallel[m/s]': pd.Series(self.vParallelMC),
                'rGyro[m]': pd.Series(self.rGyroMC),
                }
        df = pd.DataFrame(data)
        df.to_csv(f,index=False)
        f.close()
        return

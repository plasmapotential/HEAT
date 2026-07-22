#gyroClass.py
#Description:   Gyro orbit module
#Engineer:      T Looby
#Date:          Fall 2021
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import toolsClass
import multiprocessing
import time
import open3d as o3d
from scipy.interpolate import interp1d
import scipy.integrate as integrate
#from tqdm.contrib.concurrent import process_map #for process bar.  very slow...
tools = toolsClass.tools()

import logging
log = logging.getLogger(__name__)

class GYRO:

    def __init__(self, rootDir, dataPath, chmod=0o774, UID=-1, GID=-1):
        """
        rootDir is root location of python modules (where dashGUI.py lives)
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

        Gyro Orbit Heat Flux Variables:
        -------------------------------------

        :N_gyroSteps: number of discrete line segments per helical gyro period [integer].  
          Higher values mean better approximation of the helical trajectory but come at 
          the cost of longer computation times.
        :gyroTraceLength: number of steps to trace for gyro orbit calculation [integer].  Step width
          is defined by the MHD EQ variable dpinit.  Should really change this name to
          gyroTraceLength as it is not directly related to degrees.  Total toroidal distance 
          of trace is gyroTraceLength * dpinit
        :gyroT_eV: Plasma ion temperature [eV].  This temperature corresponds to the mean
          total velocity in the ion velocity distribution function.
        :N_vSlice: Number of macroparticle samples to take from the velocity distribution
          function [integer].  Each sample defines the ion energies.
        :N_vPhase: Number of macroparticle samples to take from the velocity phase 
          distribution function [integer].  Each sample the ion vPerp and v|| components. 
        :N_gyroPhase: Number of macroparticle gyroPhase samples to take from a uniform
          2pi distribution [integer].  Each sample corresponds to birth phase angle of particles about
          guiding center.
        :ionMassAMU: Ion mass in atomic mass units [AMU].
        :vMode: determines if single temperature is defined for entire PFC or if each element
          on PFC mesh has a unique plasma temperature.  Can be single or mesh.  For now,
          only single works.
        :ionFrac: fraction of PSOL carried by ions (as opposed to electrons) [0-1].  Power
          carried by the ions will be P*ionFrac
        :gyroSources: name of CAD object to be used as the gyro source plane.  

        For more information on the gyro orbit module and corresponding physics, see [6].

        """


        self.allowed_vars = [
                    'N_gyroSteps',
                    'gyroTraceLength',
                    'gyroT_eV',
                    'N_vSlice',
                    'N_vPhase',
                    'N_gyroPhase',
                    'ionMassAMU',
                    'vMode',
                    'ionFrac',
                    'gyroSources',
                            ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        integers = [
                    'N_gyroSteps',
                    'gyroTraceLength',
                    'N_vSlice',
                    'N_vPhase',
                    'N_gyroPhase',
                    ]
        floats = [
                  'ionFrac',
                  'gyroT_eV',
                  'ionMassAMU',
                ]

        for var in integers:
            if (getattr(self, var) is not None) and (~np.isnan(float(getattr(self, var)))):
                try:
                    setattr(self, var, tools.makeInt(getattr(self, var)))
                except:
                    print("Error with input file var "+var+".  Perhaps you have invalid input values?")
                    log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")
        for var in floats:
            if var is not None:
                if (getattr(self, var) is not None) and (~np.isnan(float(getattr(self, var)))):
                    try:
                        setattr(self, var, tools.makeFloat(getattr(self, var)))
                    except:
                        print("Error with input file var "+var+".  Perhaps you have invalid input values?")
                        log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")



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
        self.diamag = -1 #diamagnetism = -1 for ions, 1 for electrons

        self.mass_eV = ionMassAMU * self.AMU
        self.Z=1 #assuming isotopes of hydrogen here

        return

    def temp2thermalVelocity(self, T_eV):
        """
        Calculates thermal velocity from a temperature, where thermal velocity
        is defined as the most probable speed

        T_eV is temperature in eV

        can also be found with: d/dv( v*f(v) ) = 0

        note that this is for v, not vPerp or v||
        """
        return np.sqrt(2.0*T_eV/(self.mass_eV/self.c**2))


    def setupFreqs(self, B):
        """
        Calculates frequencies, periods, that are dependent upon B
        These definitions follow Freidberg Section 7.7.

        B is magnetic field magnitude

        """
        self.omegaGyro = self.Z * self.e * B / (self.mass_eV / self.kg2eV)
        if np.isscalar(self.omegaGyro):
            self.omegaGyro = np.array([self.omegaGyro])

        self.fGyro = np.abs(self.omegaGyro)/(2*np.pi)
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
            self.rGyro[:,i] = vPerp[i] / np.abs(self.omegaGyro)

        return

    def setupVelocities(self, N):
        """
        sets up velocities based upon vMode input from GUI

        N is the number of source mesh elements (ie len(PFC.centers) )

        len(self.t1) is number of points in divertor we are calculating HF on
        """
        #get velocity space phase angles
        self.uniformVelPhaseAngle()

        if self.vMode == 'single':
            print("Gyro orbit calculation from single plasma temperature")
            log.info("Gyro orbit calculation from single plasma temperature")
            self.T0 = np.ones((N))*self.gyroT_eV
            #get average velocity for each temperature point
            self.vThermal = self.temp2thermalVelocity(self.T0)
            #set upper bound of v*f(v) (note that this cuts off high energy particles)
            self.vMax = 5 * self.vThermal
            #get 100 points to initialize functional form of f(v) (note this is a 2D matrix cause vMax is 2D)
            self.vScan = np.linspace(0,self.vMax,10000).T
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

        in loop, i is mesh element index
        """
        self.vSlices = np.ones((len(self.T0),self.N_vSlice))*np.nan
        self.energySlices = np.zeros((len(self.T0),self.N_vSlice))
        self.energyIntegrals = np.zeros((len(self.T0),self.N_vSlice))
        self.energyFracs = np.zeros((len(self.T0),self.N_vSlice))
        self.vBounds = np.zeros((len(self.T0),self.N_vSlice+1))
        for i in range(len(self.T0)):
            #get speed range for this T0
            v = self.vScan[i,:]
            #generate the (here maxwellian) velocity vector PDF
            #pdf = lambda x: (self.mass_eV/self.c**2) / (self.T0[i]) * np.exp(-(self.mass_eV/self.c**2 * x**2) / (2*self.T0[i]) )
            pdf = lambda x: ( (self.mass_eV/self.c**2) / (2 * np.pi * self.T0[i]) )**(3.0/2.0) * np.exp(-(self.mass_eV/self.c**2 * x**2) / (2*self.T0[i]) )
            #speed pdf (integrate over solid angle)
            v_pdf = 4*np.pi * v**2 * pdf(v)
            #generate the CDF
            v_cdf = np.cumsum(v_pdf[1:])*np.diff(v)
            v_cdf = np.insert(v_cdf, 0, 0)
            #create bspline interpolators for the cdf and cdf inverse
            inverseCDF = interp1d(v_cdf, v, kind='linear')
            forwardCDF = interp1d(v, v_cdf, kind='linear')
            #CDF location of vSlices and bin boundaries
            cdfBounds = np.linspace(0,v_cdf[-1],self.N_vSlice+1)
            #CDF location of velocity bin bounds omitting 0 and 1
            #old method does not make vSlices truly bin centers
            #cdfBounds = np.linspace(0,1,self.N_vSlice+1)[1:-1]

            #old method 2 spaces centers uniformly
#            #calculate N_vSlice velocities for each pdf each with equal area (probability)
#            cdfMax = v_cdf[-1]
#            cdfMin = v_cdf[0]
#            sliceWidth = cdfMax / (self.N_vSlice+1)
#            #CDF location of vSlices omitting 0 and 1
#            cdfSlices = np.linspace(0,1,self.N_vSlice+2)[1:-1]
#            #CDF location of velocity bin bounds omitting 0 and 1
#            #old method does not make vSlices truly bin centers
#            #cdfBounds = np.linspace(0,1,self.N_vSlice+1)[1:-1]
#            #new method makes vSlices bin centers, except for the end bins
#            cdfBounds = np.diff(cdfSlices)/2.0 + cdfSlices[:-1]
#            #vSlices are Maxwellian distribution sample locations (@ bin centers)
#            self.vSlices[i,:] = inverseCDF(cdfSlices)
#            vBounds = inverseCDF(cdfBounds)
#            vBounds = np.insert(vBounds,0,0)
#            vBounds = np.append(vBounds,self.vMax[i])

            #new method spaces bins uniformly, then makes vSlices center of these bins in CDF space
            cdfSlices = np.diff(cdfBounds)/2.0 + cdfBounds[:-1]
            #vSlices are Maxwellian distribution sample locations (@ bin centers)
            self.vSlices[i,:] = inverseCDF(cdfSlices)
            vBounds = inverseCDF(cdfBounds)
            self.vBounds[i,:] = vBounds
            #print(cdfBounds)
            #print(cdfSlices)
            #print(self.vBounds)
            #print(self.vSlices)

            #Now find energies that correspond to these vSlices
            #we integrate: v**2 * f(v)
            #energy pdf (missing 1/2*mass but that gets divided out later anyways )
            #EofV = lambda x: x**2 * pdf(x)
            #EofV = lambda x: 4*np.pi * x**4 * pdf(x)

            #old HEAT method
            #f_E = lambda x: 2 * np.sqrt(x / np.pi) * (self.T0[i])**(-3.0/2.0) * np.exp(-x / self.T0[i])
            #reviewer#2 method
            f_E = lambda x: x**2 * np.exp(-x / self.T0[i])

            #energy slices that correspond to velocity slices
            self.energySlices[i,:] = f_E(0.5 * (self.mass_eV/self.c**2) * self.vSlices[i,:]**2)
            #energy integrals
            for j in range(self.N_vSlice):
                Elo = 0.5 * (self.mass_eV/self.c**2) * vBounds[j]**2
                Ehi = 0.5 * (self.mass_eV/self.c**2) * vBounds[j+1]**2
                self.energyIntegrals[i,j] = integrate.quad(f_E, Elo, Ehi)[0]
            energyTotal = self.energyIntegrals[i,:].sum()
            ##for testing
            #if i==0:
            #    print("Integral Test===")
            #    print(energyTotal)
            #    print(integrate.quad(f_E, 0.0, self.vMax[i])[0])

            #energy fractions
            for j in range(self.N_vSlice):
                self.energyFracs[i,j] = self.energyIntegrals[i,j] / energyTotal

        print("Found N_vPhase velocities of equal probability")
        log.info("Found N_vPhase velocities of equal probability")
        return

    def uniformGyroPhaseAngle(self):
        """
        Uniform sampling between 0 and 2pi

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
        #original method
        #self.vPhases = np.linspace(0.0,np.pi/2,self.N_vPhase+2)[1:-1]

        #new method
        #thankyou to NF reviewer #2 for noticing we were not using the correct
        #sampling function:  we sample uniformly in cos(2*vP) here, which is
        #the CDF of the velocity space pdf
        #with bounds in vP between (0,pi/2)
        #self.vPhases = np.arccos(np.linspace(1.0,-1.0,self.N_vPhase+2)[1:-1])/2.0
        #self.vPhases = np.arccos(1-2*np.linspace(0,1,self.N_vPhase+2))/2.0

        #beta, the vPhase angle
        b = np.linspace(0,np.pi/2,10000).T
        #heat flux pdf for velocity phase
        pdf = lambda x: np.cos(x)*np.sin(x)*2
        b_pdf = pdf(b)
        #generate the CDF
        b_cdf = np.cumsum(b_pdf[1:])*np.diff(b)
        b_cdf = np.insert(b_cdf, 0, 0)
        #create bspline interpolators for the cdf and cdf inverse
        inverseCDF = interp1d(b_cdf, b, kind='linear')
        forwardCDF = interp1d(b, b_cdf, kind='linear')
        #CDF location of vSlices and bin boundaries
        cdfBounds = np.linspace(0,b_cdf[-1],self.N_vPhase+1)
        #spaces bins uniformly, then makes vSlices center of these bins in CDF space
        cdfSlices = np.diff(cdfBounds)/2.0 + cdfBounds[:-1]
        #vSlices are Maxwellian distribution sample locations (@ bin centers)
        self.vPhases= inverseCDF(cdfSlices)
        vBounds = inverseCDF(cdfBounds)

        return


    def singleGyroTrace(self,vPerp,vParallel,gyroPhase,N_gyroSteps,
                        BtraceXYZ,controlfilePath,TGyro,rGyro,omegaGyro,
                        tag=None, verbose=True):
        """
        Calculates the gyro-Orbit path and saves to .csv and .vtk

        vPerp and vParallel [m/s] are in velocities
        gyroPhase [degrees] is initial orbit phase angle
        N_gyroSteps is number of discrete line segments per gyro period
        BtraceXYZ is the points of the Bfield trace that we will gyrate about
        """
        print("Calculating gyro trace...")
        log.info("Calculating gyro trace...")
        log.info(omegaGyro)
        #Loop thru B field trace while tracing gyro orbit
        helixTrace = None
        for i in range(len(BtraceXYZ)-1):
            #points in this iteration
            p0 = BtraceXYZ[i,:]
            p1 = BtraceXYZ[i+1,:]
            #vector
            delP = p1 - p0
            #magnitude or length of line segment
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
            x_helix = rGyro*np.cos(omegaGyro*t + gyroPhase)
            y_helix = self.diamag*rGyro*np.sin(omegaGyro*t + gyroPhase)
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
            gyroPhase = omegaGyro*t[-1] + gyroPhase
            #append to helix trace
            if helixTrace is None:
                helixTrace = helix_rot
            else:
                helixTrace = np.vstack([helixTrace,helix_rot])

        helixTrace*=1000.0 #scale for ParaView
        print("Saving data to CSV and VTK formats")
        #save data to csv format
        head = 'X[mm],Y[mm],Z[mm]'
        if tag == None:
            name = 'helixTrace'
        else:
            name = 'helixTrace_'+tag

        tools.savetxt(controlfilePath+name+'.csv', helixTrace, delimiter=',', header=head)
        #save data to vtk format
        tools.createVTKOutput(controlfilePath+name+'.csv', 'trace', name)

        if verbose==True:
            print("\nV_perp = {:f} [m/s]".format(vPerp))
            print("V_parallel = {:f} [m/s]".format(vParallel))
            print("Cyclotron Freq = {:f} [rad/s]".format(self.omegaGyro[0]))
            print("Cyclotron Freq = {:f} [Hz]".format(self.fGyro[0]))
            print("Gyro Radius = {:f} [m]".format(self.rGyro[0][0]))
            print("Number of gyro points = {:f}".format(len(helixTrace)))
            print("Longitudinal dist between gyro points = {:f} [m]".format(magP/float(Nsteps)))
            print("Each line segment length ~ {:f} [m]".format(magP))
            log.info("\nV_perp = {:f} [m/s]".format(vPerp))
            log.info("V_parallel = {:f} [m/s]".format(vParallel))
            log.info("Cyclotron Freq = {:f} [rad/s]".format(self.omegaGyro[0]))
            log.info("Cyclotron Freq = {:f} [Hz]".format(self.fGyro[0]))
            log.info("Gyro Radius = {:f} [m]".format(self.rGyro[0][0]))
            log.info("Number of gyro points = {:f}".format(len(helixTrace)))
            log.info("Longitudinal dist between gyro points = {:f} [m]".format(magP/float(Nsteps)))
            log.info("Each line segment length ~ {:f} [m]".format(magP))
        return

    def buildHelixParallel(self, i):
        """
        builds the helical trajectory of a macroparticle for the entire toroidal
        extent of the trace.  Returns this helix
        """
        #vector linking Bfield points
        delP = self.p1[i,:,:] - self.p0[i,:,:]
        for j,dP in enumerate(delP):
            #magnitude
            magP = np.sqrt(dP[0]**2 + dP[1]**2 + dP[2]**2)
            #time it takes to transit line segment
            delta_t = magP / (self.vParallelMC[i])
            #Number of steps in line segment
            Tsample = self.TGyro[i] / self.N_gyroSteps
            Nsteps = int(delta_t / Tsample)
            #length (in time) along guiding center
            t = np.linspace(0,delta_t,Nsteps+1)
            #guiding center location
            xGC = np.linspace(self.p0[i,j,0],self.p1[i,j,0],Nsteps+1)
            yGC = np.linspace(self.p0[i,j,1],self.p1[i,j,1],Nsteps+1)
            zGC = np.linspace(self.p0[i,j,2],self.p1[i,j,2],Nsteps+1)
            arrGC = np.vstack([xGC,yGC,zGC]).T
            # construct orthogonal system for coordinate transformation
            w = dP
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
            rGyro = self.rGyroMC[i]
            omega = self.omegaGyro[i]
            theta = self.lastPhase[i]
            x_helix = rGyro*np.cos(omega*t + theta)
            y_helix = self.diamag*rGyro*np.sin(omega*t + theta)
            z_helix = np.zeros((len(t)))
            #perform rotation to field line reference frame
            helix = np.vstack([x_helix,y_helix,z_helix]).T
            helix_rot = np.zeros((len(helix),3))
            for k,coord in enumerate(helix):
                helix_rot[k,:] = helix[k,0]*u + helix[k,1]*v + helix[k,2]*w
            #perform translation to field line reference frame
            helix_rot[:,0] += xGC
            helix_rot[:,1] += yGC
            helix_rot[:,2] += zGC

            #shift entire helix to ensure we capture intersections in p0 plane
            helix_rot[:,0] += w[0]*0.0003
            helix_rot[:,1] += w[1]*0.0003
            helix_rot[:,2] += w[2]*0.0003

            #update gyroPhase variable so next iteration starts here
            lastPhase = omega*t[-1] + theta

            if j==0:
                helix_final = helix_rot
            else:
                helix_final = np.vstack([helix_final, helix_rot])
        return helix_final

    def multipleGyroTraceOpen3D(self):
        """
        Calculates the helical trajectory of ions, then checks if any mesh
        elements are intersected along the path.  This loops thru each
        mesh element in the source plane (i) and checks for ray-triangle
        intersections between the ray and all the mesh in the PFC's
        intersectList.

        Uses Open3D to accelerate the calculation:
        Zhou, Qian-Yi, Jaesik Park, and Vladlen Koltun. "Open3D: A modern
        library for 3D data processing." arXiv preprint arXiv:1801.09847 (2018).

        meant to be run once per field line step
        """
        N = len(self.p1)
        intersectRecord = np.ones((N))*np.nan
        hdotn = np.ones((N))*np.nan
        lastPhases = np.zeros((N))

        #cast variables to 32bit for C
        vertices = np.array(self.targets, dtype=np.float32)
        triangles = np.array(np.arange(self.PFC_Nt*3).reshape(self.PFC_Nt,3), dtype=np.uint32)
        #triangles = np.array(np.repeat(np.arange(self.PFC_Nt),3).reshape(self.PFC_Nt,3), dtype=np.uint32)

        #build intersection mesh and tensors for open3d ray-tri calcs
        mesh = o3d.t.geometry.TriangleMesh()
        scene = o3d.t.geometry.RaycastingScene()
        mesh_id = scene.add_triangles(vertices, triangles)

        #loop thru each gyroSource mesh face looking for intersections along helix
        for i in range(N):
            #magnetic field trace
            self.helixTrace = [None] * len(self.p0)
            #vector linking Bfield points
            delP = self.p1[i] - self.p0[i]
            #magnitude
            magP = np.sqrt(delP[0]**2 + delP[1]**2 + delP[2]**2)
            #time it takes to transit line segment
            delta_t = magP / (self.vParallelMC[self.GYRO_HLXmap][i])
            #Number of steps in line segment
            Tsample = self.TGyro[self.GYRO_HLXmap][i] / self.N_gyroSteps
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
            rGyro = self.rGyroMC[self.GYRO_HLXmap][i]
            omega = self.omegaGyro[self.GYRO_HLXmap][i]
            theta = self.lastPhase[self.GYRO_HLXmap][i]
            x_helix = rGyro*np.cos(omega*t + theta)
            y_helix = self.diamag*rGyro*np.sin(omega*t + theta)
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
            lastPhase = omega*t[-1] + theta
            lastPhases[i] = lastPhase

            #=== intersection checking begins here ===
            q1 = helix_rot[:-1,:]
            q2 = helix_rot[1:,:]

            r = q2-q1
            rMag = np.linalg.norm(r, axis=1)
            rNorm = r / rMag.reshape((-1,1))

            #calculate intersections
            mask = np.ones((len(q2)))
            rays = o3d.core.Tensor([np.hstack([np.float32(q1),np.float32(rNorm)])],dtype=o3d.core.Dtype.Float32)
            hits = scene.cast_rays(rays)
            #convert open3d CPU tensors back to numpy
            hitMap = hits['primitive_ids'][0].numpy()
            distMap = hits['t_hit'][0].numpy()

            #escapes occur where we have 32 bits all set: 0xFFFFFFFF = 4294967295 base10
            escapes = np.where(hitMap == 4294967295)[0]
            mask[escapes] = 0.0

            #when distances to target exceed the trace step, we exclude any hits
            tooLong = np.where(distMap > rMag)[0]
            mask[tooLong] = 0.0

            if np.sum(mask)>0:
                #first hit instance
                idxHit = np.where(mask==1)[0][0]
                intersectRecord[i] = self.PFCintersectMap[hitMap[idxHit]]
                hdotn[i] = np.dot(self.intersectNorms[int(intersectRecord[i])],rNorm[idxHit])
            else:
                #return nan if we didnt hit anything
                intersectRecord[i] = np.nan
                hdotn[i] = np.nan

            #print the trace for a specific index
            if self.traceIndex2 is not None:
                if self.traceIndex2 == i:
                    ##for testing
                    #if np.sum(mask)>0:
                    #    print(idxHit)
                    #    print(intersectRecord[i])
                    #    print(rNorm[idxHit])
                    #    print(self.intersectNorms[int(intersectRecord[i])])
                    #print(hitMap)
                    #print(distMap)
                    #print(rMag)
                    #print(mask)
                    #print(intersectRecord[i])

                    #print("Saving Index data to CSV and VTK formats")
                    #save data to csv format
                    head = 'X[mm],Y[mm],Z[mm]'
                    tools.savetxt(self.controlfilePath+'helix{:d}.csv'.format(self.N_GCdeg), helix_rot*1000.0, delimiter=',', header=head)
                    #save data to vtk format
                    tools.createVTKOutput(self.controlfilePath+'helix{:d}.csv'.format(self.N_GCdeg),
                                        'trace', 'traceHelix{:d}'.format(self.N_GCdeg),verbose=False)
                    #guiding center
                    #np.savetxt(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg), arrGC*1000.0, delimiter=',', header=head)
                    #save data to vtk format
                    #tools.createVTKOutput(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg),
                    #                        'trace', 'traceGC{:d}'.format(self.N_GCdeg),verbose=False)
                    print("Intersection Location for ROIidx: {:f}".format(intersectRecord[i]))



        self.lastPhase[self.GYRO_HLXmap] = lastPhases

        print('Parallel helix trace complete')
        log.info('Parallel helix trace complete')
        return intersectRecord, hdotn

    def multipleGyroTraceOpen3D_NoLoop(self):
        """
        Calculates the helical trajectory of ions, then checks if any mesh
        elements are intersected along the path.  This loops thru each
        mesh element in the source plane (i) and checks for ray-triangle
        intersections between the ray and all the mesh in the PFC's
        intersectList.

        Currently, this is called once every dpinit steps in MAFOT.  If this
        is very slow, you could run all the trace steps in a single calculation
        by removing the loop in PFCclass


        Uses Open3D to accelerate the calculation:
        Zhou, Qian-Yi, Jaesik Park, and Vladlen Koltun. "Open3D: A modern
        library for 3D data processing." arXiv preprint arXiv:1801.09847 (2018).

        meant to be run once for all field line steps
        """
        N = self.Nctrs
        intersectRecord = np.ones((N))*np.nan
        hdotn = np.ones((N))*np.nan


        #cast variables to 32bit for C
        vertices = np.array(self.targets, dtype=np.float32)
        triangles = np.array(np.arange(self.PFC_Nt*3).reshape(self.PFC_Nt,3), dtype=np.uint32)
        #triangles = np.array(np.repeat(np.arange(self.PFC_Nt),3).reshape(self.PFC_Nt,3), dtype=np.uint32)

        #build intersection mesh and tensors for open3d ray-tri calcs
        mesh = o3d.t.geometry.TriangleMesh()
        scene = o3d.t.geometry.RaycastingScene()
        mesh_id = scene.add_triangles(vertices, triangles)

        #Prepare helical trace across multiple cores
        Ncores = multiprocessing.cpu_count() - 2 #reserve 2 cores for overhead
        #in case we run on single core machine
        if Ncores <= 0:
            Ncores = 1
        #the overhead on this calculation can be high, so limit to 3 cores
        elif Ncores > 3:
            Ncores = 3
        print('Initializing parallel helix trace across {:d} cores'.format(Ncores))
        log.info('Initializing parallel helix trace across {:d} cores'.format(Ncores))
        #each worker receives a single start and end point (p0 and p1),
        #corresponding to one trace from the MAFOT structure output.
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
        #multiprocessing with normal methods
        #Do this try clause to kill any zombie threads that don't terminate
        t0 = time.time()
        try:
            pool = multiprocessing.Pool(3)
            output = np.array(pool.map(self.buildHelixParallel, np.arange(N)))
        finally:
            pool.close()
            pool.join()
            del pool

        helix = output
        print("Helix trace took {:f} seconds".format(time.time() - t0))
        log.info("Helix trace took {:f} seconds".format(time.time() - t0))

        print("Initializing open3D ray-triangle intersection checks")
        log.info("Initializing open3D ray-triangle intersection checks")
        for i in range(N):
            #=== intersection checking begins here ===
            q1 = helix[i][:-1,:]
            q2 = helix[i][1:,:]

            r = q2-q1
            rMag = np.linalg.norm(r, axis=1)
            rNorm = r / rMag.reshape((-1,1))

            #calculate intersections
            mask = np.ones((len(q2)))
            rays = o3d.core.Tensor([np.hstack([np.float32(q1),np.float32(rNorm)])],dtype=o3d.core.Dtype.Float32)
            hits = scene.cast_rays(rays)
            #convert open3d CPU tensors back to numpy
            hitMap = hits['primitive_ids'][0].numpy()
            distMap = hits['t_hit'][0].numpy()

            #escapes occur where we have 32 bits all set: 0xFFFFFFFF = 4294967295 base10
            escapes = np.where(hitMap == 4294967295)[0]
            mask[escapes] = 0.0

            #when distances to target exceed the trace step, we exclude any hits
            tooLong = np.where(distMap > rMag)[0]
            mask[tooLong] = 0.0

            if np.sum(mask)>0:
                #first hit instance
                idxHit = np.where(mask==1)[0][0]
                intersectRecord[i] = self.PFCintersectMap[hitMap[idxHit]]
                hdotn[i] = np.dot(self.intersectNorms[int(intersectRecord[i])],rNorm[idxHit])
            else:
                #return nan if we didnt hit anything
                intersectRecord[i] = np.nan
                hdotn[i] = np.nan

            #print the trace for a specific index
            if self.traceIndex2 is not None:
                if self.traceIndex2 == i:
                    ##for testing
                    #if np.sum(mask)>0:
                    #    print(idxHit)
                    #    print(intersectRecord[i])
                    #    print(rNorm[idxHit])
                    #    print(self.intersectNorms[int(intersectRecord[i])])
                    #print(hitMap)
                    #print(distMap)
                    #print(rMag)
                    #print(mask)
                    #print(intersectRecord[i])

                    #print("Saving Index data to CSV and VTK formats")
                    #save data to csv format
                    head = 'X[mm],Y[mm],Z[mm]'
                    tools.savetxt(self.controlfilePath+'helix{:d}.csv'.format(i), helix[i]*1000.0, delimiter=',', header=head)
                    #save data to vtk format
                    tools.createVTKOutput(self.controlfilePath+'helix{:d}.csv'.format(i),
                                        'trace', 'traceHelix{:d}'.format(i),verbose=False)
                    #guiding center
                    #np.savetxt(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg), arrGC*1000.0, delimiter=',', header=head)
                    #save data to vtk format
                    #tools.createVTKOutput(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg),
                    #                        'trace', 'traceGC{:d}'.format(self.N_GCdeg),verbose=False)
                    print("Intersection Location for ROIidx: {:f}".format(intersectRecord[i]))



        print('Parallel helix trace complete')
        log.info('Parallel helix trace complete')
        return intersectRecord, hdotn

    def gyroTraceParallel(self, i, mode='MT'):
        """
        parallelized gyro trace.  called by multiprocessing.pool.map()

        i is index of parallel run from multiprocessing, corresponds to a mesh face
        we are tracing in the ROI

        writes helical trace to self.helixTrace[i] in 2D matrix format:
            columns = X,Y,Z
            rows = steps up helical trace

        also updates self.lastPhase for use in next iteration step

        mode options are:
        -Signed Volume Loop: 'SigVolLoop'
        -Signed Volume Matrix:  'SigVolMat'
        -Moller-Trumbore Algorithm: 'MT'
        """
        #vector
        delP = self.p1[i] - self.p0[i]
        #magnitude
        magP = np.sqrt(delP[0]**2 + delP[1]**2 + delP[2]**2)
        #time it takes to transit line segment
        delta_t = magP / (self.vParallelMC[self.GYRO_HLXmap][i])
        #Number of steps in line segment
        Tsample = self.TGyro[self.GYRO_HLXmap][i] / self.N_gyroSteps
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
        rGyro = self.rGyroMC[self.GYRO_HLXmap][i]
        omega = self.omegaGyro[self.GYRO_HLXmap][i]
        theta = self.lastPhase[self.GYRO_HLXmap][i]
        x_helix = rGyro*np.cos(omega*t + theta)
        y_helix = self.diamag*rGyro*np.sin(omega*t + theta)
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
        lastPhase = omega*t[-1] + theta

        #=== intersection checking ===
        q1 = helix_rot[:-1,:]
        q2 = helix_rot[1:,:]

        #Filter by psi
        if self.psiFilterSwitch == True:
            psiP1 = self.PFC_psiP1
            psiP2 = self.PFC_psiP2
            psiP3 = self.PFC_psiP3
            psiMin = self.psiMin[i]
            psiMax = self.psiMax[i]

            #account for psi sign convention
            if psiMin > psiMax:
                pMin = psiMax
                pMax = psiMin
            else:
                pMin = psiMin
                pMax = psiMax

            #target faces outside of this toroidal slice
            test0 = np.logical_and(np.logical_and(psiP1 < pMin, psiP2 < pMin), psiP3 < pMin)
            test1 = np.logical_and(np.logical_and(psiP1 > pMax, psiP2 > pMax), psiP3 > pMax)
            test = np.logical_or(test0,test1)
            usePsi = np.where(test == False)[0]

        else:
            usePsi = np.arange(len(self.PFC_t1))

        #Filter by toroidal angle
        if self.phiFilterSwitch == True:
            phiP1 = self.PFC_phiP1
            phiP2 = self.PFC_phiP2
            phiP3 = self.PFC_phiP3
            phiMin = self.phiMin[i]
            phiMax = self.phiMax[i]

            #angle wrap cases (assumes we never trace in MAFOT steps larger than 10degrees)
            if np.abs(phiMin-phiMax) > np.radians(5):
                phiP1[phiP1<0] += 2*np.pi
                phiP2[phiP2<0] += 2*np.pi
                phiP3[phiP3<0] += 2*np.pi
                if phiMin < 0: phiMin+=2*np.pi
                if phiMax < 0: phiMax+=2*np.pi

            #account for toroidal sign convention
            if phiMin > phiMax:
                pMin = phiMax
                pMax = phiMin
            else:
                pMin = phiMin
                pMax = phiMax


            #target faces outside of this toroidal slice
            test0 = np.logical_and(np.logical_and(phiP1 < pMin, phiP2 < pMin), phiP3 < pMin)
            test1 = np.logical_and(np.logical_and(phiP1 > pMax, phiP2 > pMax), phiP3 > pMax)
            test = np.logical_or(test0,test1)
            usePhi = np.where(test == False)[0]

        else:
            usePhi = np.arange(len(self.PFC_t1))

        #combine filter algorithms
        use = np.intersect1d(usePsi,usePhi)

        Nt = len(use)

        t0 = time.time()
        #using full array (no for loop)
        if mode == 'SigVolMat':
            q13D = np.repeat(q1[:,np.newaxis], Nt, axis=1)
            q23D = np.repeat(q2[:,np.newaxis], Nt, axis=1)

            sign1 = np.sign(tools.signedVolume2(q13D,self.PFC_t1[use],self.PFC_t2[use],self.PFC_t3[use],ax=2))
            sign2 = np.sign(tools.signedVolume2(q23D,self.PFC_t1[use],self.PFC_t2[use],self.PFC_t3[use],ax=2))
            sign3 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t1[use],self.PFC_t2[use],ax=2))
            sign4 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t2[use],self.PFC_t3[use],ax=2))
            sign5 = np.sign(tools.signedVolume2(q13D,q23D,self.PFC_t3[use],self.PFC_t1[use],ax=2))

            test1 = (sign1 != sign2)
            test2 = np.logical_and(sign3==sign4,sign3==sign5)
            loc = np.where(np.logical_and(test1,test2))

            #result=1 if we intersected, otherwise NaN
            if np.sum(loc) > 0:
                #only take first index (ie first intersection location)
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
        elif mode=='SigVolLoop':
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

        # Intersection check using adapted version of Moller-Trumbore Algorithm:
        # Möller, Tomas; Trumbore, Ben (1997). "Fast, Minimum Storage Ray-Triangle Intersection".
        #      Journal of Graphics Tools. 2: 21–28. doi:10.1080/10867651.1997.10487468.
        else:
            E1 = (self.PFC_t2[use] - self.PFC_t1[use])
            E2 = (self.PFC_t3[use] - self.PFC_t1[use])
            D = (q2-q1)
            Dmag = np.linalg.norm(D, axis=1)
            eps = 0.0
            for j in range(len(helix_rot)-1):
                D[j] = D[j] / np.linalg.norm(D, axis=1)[j]
                h = np.cross(D[j], E2)
                a = np.sum(E1*h, axis=1)
                test1 = np.logical_and( a>-eps, a<eps) #ray parallel to triangle
                with np.errstate(divide='ignore', invalid='ignore'):
                    #test1 = a<eps #ray parallel to triangle
                    f=1.0/a
                    s = q1[j] - self.PFC_t1[use]
                    u = f * np.sum(s*h, axis=1)
                    test2 = np.logical_or(u<0.0, u>1.0) #ray inside triangle
                    q = np.cross(s,E1)
                    v = f*np.sum(D[j]*q, axis=1)
                    test3 =  np.logical_or(v<0.0, (u+v)>1.0) #ray inside triangle
                    l = f*np.sum(E2*q, axis=1)
                    test4 = np.logical_or(l<0.0, l>Dmag[j]) #ray long enough to intersect triangle
                if np.sum(~np.any([test1,test2,test3,test4], axis=0))>0:
                    #we assume first intersection in this array is the intersection
                    PFC_index = use[ np.where(np.any([test1,test2,test3,test4], axis=0)==False)[0][0] ]
                    #map this index (of self.PFC_tX) back to global index (of self.tX)
                    index = self.PFCintersectMap[PFC_index]
                    #gyro trace incident angle:
                    vec = (q2[j] - q1[j]) / np.linalg.norm(q2[j]-q1[j])
                    hdotn = np.dot(self.intersectNorms[index],vec)
                    break
                else:
                    PFC_index = np.NaN
                    index = np.NaN
                    hdotn = np.NaN


        #print the trace for a specific index
        if self.traceIndex2 is not None:
            if self.traceIndex2 == i:
                if np.sum(~np.any([test1,test2,test3,test4], axis=0))>0:
                    print("TEST====")
                    print(use[ np.where(np.any([test1,test2,test3,test4], axis=0)==False)[0] ])
                #print("Saving Index data to CSV and VTK formats")
                #save data to csv format
                head = 'X[mm],Y[mm],Z[mm]'
                tools.savetxt(self.controlfilePath+'helix{:d}.csv'.format(self.N_GCdeg), helix_rot*1000.0, delimiter=',', header=head)
                #save data to vtk format
                tools.createVTKOutput(self.controlfilePath+'helix{:d}.csv'.format(self.N_GCdeg),
                                        'trace', 'traceHelix{:d}'.format(self.N_GCdeg),verbose=False)
                #guiding center
                #np.savetxt(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg), arrGC*1000.0, delimiter=',', header=head)
                #save data to vtk format
                #tools.createVTKOutput(self.controlfilePath+'GC{:d}.csv'.format(self.N_GCdeg),
                #                        'trace', 'traceGC{:d}'.format(self.N_GCdeg),verbose=False)

                print("Intersection Index: {:f}".format(index))
                print("PFC Index: {:f}".format(PFC_index))



        t1 = time.time() - t0
        return lastPhase, index, hdotn, t1

    def multipleGyroTrace(self):
        """
        Calculates the helical path for multiple points,
        each with different gyroRadii, using multiprocessing

        Btrace is one step of a field line trace for each point (MAFOT structure output)
        phase is phase angle (updated from last trace step)

        updates lastPhase variable and helixTrace
        """
#        #include toroidal angle filtering
#        GYRO.phiFilterSwitch = False
        #magnetic field trace
        self.helixTrace = [None] * len(self.p0)
        N = len(self.p1)
        #Prepare helical trace across multiple cores
        Ncores = multiprocessing.cpu_count() - 2 #reserve 2 cores for overhead
        #in case we run on single core machine
        if Ncores <= 0:
            Ncores = 1
        print('Initializing parallel helix trace across {:d} cores'.format(Ncores))
        log.info('Initializing parallel helix trace across {:d} cores'.format(Ncores))
        #each worker receives a single start and end point (p0 and p1),
        #corresponding to one trace from the MAFOT structure output.
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
        #multiprocessing with normal methods
        #Do this try clause to kill any zombie threads that don't terminate
        try:
            pool = multiprocessing.Pool(Ncores)
            output = np.asarray(pool.map(self.gyroTraceParallel, np.arange(N)))
        finally:
            pool.close()
            pool.join()
            del pool

        #multiprocessing with status bar (equiv to multiprocessing.Pool.map())
#        print("Multiprocessing gyro trace:")
#        output = process_map(self.gyroTraceParallel, range(N), max_workers=Ncores, chunksize=1)
#        output = np.asarray(output)

        intersectRecord = output[:,1]
#        use = np.where(np.isnan(intersectRecord)==True)[0]
#        self.lastPhase = output[:,0][use]
        self.lastPhase[self.GYRO_HLXmap] = output[:,0]
        #uncomment for gyro trace incident angle:
        hdotn = output[:,2]

        #uncomment for avg time / calc
        print("Intersection Calc. Avg. time = {:f} [s]".format(np.sum(output[:,3]) / N))
        log.info("Intersection Calc. Avg. time = {:f} [s]".format(np.sum(output[:,3]) / N))
        print('Parallel helix trace complete')
        log.info('Parallel helix trace complete')
        return intersectRecord, hdotn

    def writeIntersectRecord(self, gyroPhase, vPhase, vSlice, faces, file):
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
        rec = self.intersectRecord[gyroPhase,vPhase,vSlice,:]
        data = {
                'face': pd.Series(faces),
                'intersectFace': pd.Series(rec),
                'vPerp[m/s]': pd.Series(self.vPerpMC),
                'vParallel[m/s]': pd.Series(self.vParallelMC),
                'rGyro[m]': pd.Series(self.rGyroMC),
                }
        df = pd.DataFrame(data)
        df.to_csv(f,index=False)
        f.close()
        return

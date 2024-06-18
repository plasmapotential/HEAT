#runawayClass.py
#Description:   Runaway Electron Module module
#Engineer:      A Feyrer
#Date:          April 2024

"""
Class for runaway electrons
"""

import os
import sys
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import scipy.interpolate as interp
from scipy.interpolate import interp1d
import scipy.integrate as integrate
import pandas as pd
import open3d as o3d
import time
import shutil


#HEAT classes
import toolsClass
tools = toolsClass.tools()
import ioClass
IO = ioClass.IO_HEAT()

import logging
log = logging.getLogger(__name__)


class Runaways:

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
        tools.physicsConstants(self)
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
        :RE_Eav: Average energy of the runaway electron distribution function
        :RE_I: Total current carried by runaway electrons
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
        :gyroSources: name of CAD object to be used as the gyro source plane.  

        For more information on the gyro orbit module and corresponding physics, see [6].

        """


        self.allowed_vars = [
                    'N_gyroSteps',
                    'gyroTraceLength',
                    'RE_Eav',
                    'RE_I',
                    'N_vSlice',
                    'N_vPhase',
                    'N_gyroPhase',
                    'vMode',
                    'gyroSources',
                            ] #check take out ionMassAMU?
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
                  'RE_Eav',
                  'RE_I',
                  
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
        self.me = 0.511e6 #eV

        self.mass_eV = ionMassAMU * self.AMU
        self.Z=1 #assuming isotopes of hydrogen here check

        return
        
    """
    Helpful functions for dealing with relativstic dynamics
    """
    
    def E2v_electron(self, RE_KE):
    	"""
    	Calculates the velocity from the energy of an electron
    	
    	E is in eV
    	
    	Assumes mass energy is not included in E
    	
    	   	
    	"""
    	return np.sqrt(self.c**2 * (1 - (self.me)**2/(RE_KE + self.me)**2))


    def calc_gamma(self, v):
        """
        Takes in the velocity and returns the relativistic gamma
        necessary when calculating helical trajectories
        """
        
        return (1/(1 - v**2/self.c**2))
    
    
    def readREFile(self, path:str):
        """
        reads a filament csv input file

        generates a data pd array
        """

        filFile = path + 'Runaways.csv'
        self.REData = pd.read_csv(filFile, sep=',', comment='#', skipinitialspace=True)
        print(self.REData.columns)
        print(self.REData.to_dict().keys())
        return

    """
    Add intialize_REdist() (not from dict)? 
    """    
    
    def setupRETime(self):
        """
        sets up filament timesteps using data from filament file
        """
        try:
            df = self.REData      
            df['N_dt'] = (df['tMax[s]'] - df['tMin[s]']) / df['dt[s]']
            df['N_dt'] = df['N_dt'].round(0)
            df['N_dt'] = df['N_dt'].astype(int)
            df['tNew'] = df.apply( lambda x:  np.linspace(x['tMin[s]'],x['tMax[s]'],x['N_dt']+1), axis=1)
            self.tsFil = df['tNew'].values

        except:
            print("\n\nError setting up filament time.  Check your filament file\n\n")
            log.info("\n\nError setting up filament time.  Check your filament file\n\n")
            sys.exit()
        return
    
    def initializeREdistFromDict(self, REDict:dict, id:int, ep:object):
        """
        initializes Runaway electron distribution from a dictionary (read from csv file in HEAT)

        dictionary is originally a pandas dataframe, so each row in the filament.csv
        HEAT file is a new id # in the dict.  each row can therefore be accessed by:
        filDict['<parameter>'][<row#>]
        ep is equilParams object
        """
        #initialize this filament
        try:
            self.rCtr = REDict['rCtr[m]'][id]
            self.zCtr = REDict['zCtr[m]'][id]
            self.phi = np.radians(REDict['phiCtr[deg]'][id])
            self.xCtr = self.rCtr * np.cos(self.phi)
            self.yCtr = self.rCtr * np.sin(self.phi)
            self.sig_b = REDict['sig_b[m]'][id]
            self.sig_r = REDict['sig_r[m]'][id]
            self.sig_p = REDict['sig_p[m]'][id]
            self.tMin = REDict['tMin[s]'][id]
            self.tMax = REDict['tMax[s]'][id]
            self.dt = REDict['dt[s]'][id]
            self.N_sig_r = REDict['N_sig_r'][id]
            self.N_sig_p = REDict['N_sig_p'][id]
            self.N_sig_b = REDict['N_sig_b'][id]
            self.N_r = REDict['N_r'][id]
            self.N_p = REDict['N_p'][id]
            self.N_b = REDict['N_b'][id]
            self.N_vS = REDict['N_vS'][id]
            self.E_av = REDict['E_av[eV]'][id]
            self.decay_t = REDict['decay_t[s]'][id]
            self.N_src_t = REDict['N_src_t'][id]
            self.E_max = REDict['E_max[eV]'][id]
            self.E_min = REDict['E_min[eV]'][id]
            self.IRE = REDict['IRE[A]'][id]
            self.I_dir = REDict['I_dir'][id] #Electrons will be moving in opposite direction of I_dir
            # self.Uniform_vel todo
            self.v_r = 0 #backdoor thing, should  probably change this

            self.ep = ep

        except:
            print("\n\nCould not initialize runaways.  Check your input file!")
            log.info("\n\nCould not initialize runaways.  Check your input file!")
            sys.exit()

        return
    
    def calcTotalE(self):
        vave = self.E2v_electron(self.E_av)
        Nre = self.IRE*2*np.pi*self.ep.g['RmAxis']/(vave* self.e)
        self.E0 = Nre * self.E_av * self.e # check e

    
    def RECtrBtrace(self, MHD:object, t:float):
        """
        traces magnetic field at RE dist center
        
        blatantly stolen from filament class

        PFC object used to get paths
        MHD object used for traces
        t used for name of geqdsk for MAFOT
        """
        #calculate field line for filament center
        xyz = np.array([self.xCtr,self.yCtr,self.zCtr])
        dphi = 1.0
        MHD.ittStruct = 360.0 #360 degree trace
        MHD.writeMAFOTpointfile(xyz,self.gridfileStruct)
        MHD.writeControlFile(self.controlfilePath+self.controlfileStruct, t, 0, mode='struct') #0 for both directions
        MHD.getFieldpath(dphi, self.gridfileStruct, self.controlfilePath, self.controlfileStruct, paraview_mask=False, tag=None)
        Btrace = tools.readStructOutput(self.structOutfile) 
        os.remove(self.structOutfile)
        return Btrace

    def findGuidingCenterPaths(self, pts:np.ndarray, MHD:object, traceDir:float):
        """
        creates MAFOT point files and control files, then runs heatstructure MAFOT program
        
        Again stolen from filament class
        
        todo: ask Alex or Tom about directional defintions 
        
        pts are pts we want to trace from
        MHD is a HEAT MHD object
        traceDir is not self.traceDir, but rather a variable set by filament tracer that 
           defines MAFOT mapDirection
           
        """
        MHD.ittStruct = 1.0
        #forward trace
        MHD.writeMAFOTpointfile(pts,self.gridfileStruct)
        MHD.writeControlFile(self.controlfilePath+self.controlfileStruct, self.tEQ, traceDir, mode='struct') #0 for both directions
        MHD.getMultipleFieldPaths(1.0, self.gridfileStruct, self.controlfilePath,
                                self.controlfileStruct)

        structData = tools.readStructOutput(self.structOutfile)
        os.remove(self.structOutfile) #clean up
        #traces are always saved in positive toroidal direction
        if traceDir == 1.0:
            q1 = structData[0::2,:] #even indexes are first trace point
            q2 = structData[1::2,:] #odd indexes are second trace point
        elif traceDir == -1.0:
            q1 = structData[1::2,:] #odd indexes are first trace point
            q2 = structData[0::2,:] #even indexes are second trace point

        return q1, q2
    
    def buildIntersectionMesh(self):
        """
        builds intersection mesh for open3d ray tracing
        """
        Nt = len(self.t1)
        targets = np.hstack([self.t1, self.t2, self.t3]).reshape(Nt*3,3)
        #cast variables to 32bit for C
        vertices = np.array(targets.reshape(Nt*3,3), dtype=np.float32)
        triangles = np.array(np.arange(Nt*3).reshape(Nt,3), dtype=np.uint32)
        #build intersection mesh and tensors for open3d
        self.mesh = o3d.t.geometry.TriangleMesh()
        self.scene = o3d.t.geometry.RaycastingScene()
        self.mesh_id = self.scene.add_triangles(vertices, triangles)
        return

    def traceREParticles(self, MHD: object, ts:np.ndarray , tIdx:int):
        """
        Traces RE macro-particles
        
        
        Based on the filament tracing code 
        
        todo: figure this out, think about it a bit harder
        
        Uses Open3D to accelerate the calculation:
        Zhou, Qian-Yi, Jaesik Park, and Vladlen Koltun. "Open3D: A modern
        library for 3D data processing." arXiv preprint arXiv:1801.09847 (2018).
        """
        pts = self.xyzPts.reshape(self.N_b*self.N_r*self.N_p, 3)
        print(pts)
        N_pts = len(pts)
        N_ts = len(ts)
        print(ts)
        print(self.N_vS)


        print('Number of target faces: {:f}'.format(self.Nt))
        log.info('Number of target faces: {:f}'.format(self.Nt))
        print('Number of filament source points: {:f}'.format(N_pts))
        log.info('Number of filament source points: {:f}'.format(N_pts))
        # x = 5/0
        #setup velocities
        self.setupParallelVelocities()
        #initialize trace step matrix
        self.xyzSteps = np.zeros((self.N_vS, N_pts, N_ts, 3))

        #build intersection mesh
        self.buildIntersectionMesh()
        intersectRecord = np.ones((self.N_vS, N_pts, N_ts))*np.nan
        #hdotn = np.ones((N_pts))*np.nan
        #initialize shadowMask matrix
        shadowMask = np.zeros((N_pts))

        for i in range(self.N_vS):
            print("\n---Tracing for velocity slice: {:f} [m/s]---\n".format(self.vSlices[0,i]))
            log.info("\n---Tracing for velocity slice: {:f} [m/s]---\n".format(self.vSlices[0,i]))
            #update time counting variables
            tsNext = np.ones((N_pts))*ts[tIdx]
            t_tot = np.ones((N_pts))*ts[tIdx]
            tIdxNext = np.zeros((N_pts), dtype=int) + tIdx
            vec = np.zeros((N_pts, 3))
            frac = np.zeros((N_pts))
            use = np.arange(N_pts)

            #save macroparticle coordinates at t0, increment next
            self.xyzSteps[i,:,tIdxNext,:] = pts
            tIdxNext += 1
            tsNext = ts[tIdxNext]
            launchPt = pts.copy()
            v_b = self.vSlices[:,i]
            sumCount = 0

            #check if v_b is positive or negative and trace accordingly
            # todo ask about this? 
            if v_b[0] > 0:
                traceDir = 1.0
            else:
                traceDir = -1.0

            #walk along field line, correcting for v_r, looking for intersections,
            #saving coordinates at ts
            while len(use) > 0:
                #calculate guiding center path for one step
                q1, q2 = self.findGuidingCenterPaths(launchPt[use], MHD, traceDir)
                #correct guiding center path using radial velocity
                d_b = np.linalg.norm((q2-q1), axis=1) #distance along field
                t_b = d_b / np.abs(v_b)[use] #time along field for this step
                t_tot[use] += t_b #cumulative time along field
                d_r = 0 # I think the radial distance will be approx 0 here? t_b * self.v_r #radial distance

                #generate local radial coordinate transformation
                R = np.sqrt(q2[:,0]**2 + q2[:,1]**2)
                Z = q2[:,2]
                phi = np.arctan2(q2[:,1], q2[:,0])
                rVec = self.fluxSurfNorms(self.ep, R, Z)
                rX, rY, rZ = tools.cyl2xyz(rVec[:,0], rVec[:,2], phi, degrees=False)
                rVecXYZ = np.vstack((rX,rY,rZ)).T
                rMag = np.linalg.norm(rVecXYZ, axis=1)
                rN = rVecXYZ / rMag[:,None]
                
                #magnetic field line trace corrected with radial (psi) velocity component
                q3 = q2 #+ d_r[:,np.newaxis] * rN
                vec[use] = q3-q1
                frac[use] = (tsNext[use] - t_tot[use] + t_b) / t_b

                if len(use) > 0:
                    launchPt[use] = q3

                #construct rays
                rMag = np.linalg.norm(vec[use], axis=1)
                rNorm = vec[use] / rMag.reshape((-1,1))

                #calculate intersections
                mask = np.ones((len(q3)))
                rays = o3d.core.Tensor([np.hstack([np.float32(q1),np.float32(rNorm)])],dtype=o3d.core.Dtype.Float32)
                hits = self.scene.cast_rays(rays)
                #convert open3d CPU tensors back to numpy
                hitMap = hits['primitive_ids'][0].numpy()
                distMap = hits['t_hit'][0].numpy()

                #escapes occur where we have 32 bits all set: 0xFFFFFFFF = 4294967295 base10
                escapes = np.where(hitMap == 4294967295)[0]
                mask[escapes] = 0.0
                #when distances to target exceed the trace step, we exclude any hits
                tooLong = np.where(distMap > rMag)[0]
                mask[tooLong] = 0.0

                #put hit index in intersectRecord for tIdxNext if we hit
                if np.sum(mask)>0:
                    #we found a hit, or multiple
                    idxHit = np.where(mask==1)[0]
                    intersectRecord[i, use[idxHit], tIdxNext[use[idxHit]]] = hitMap[idxHit]
                    #hdotn[i] = np.dot(self.intersectNorms[int(intersectRecord[i])],rNorm[idxHit])
                    sumCount += np.sum(mask)
                #else:
                #    #return nan if we didnt hit anything
                #    intersectRecord[i, use, tIdxNext[use]] = np.nan
                #    #hdotn[i] = np.nan

                #particles we need to keep tracing (didnt hit and less than tMax)
                test1 = np.where(np.isnan(intersectRecord[i, use, tIdxNext[use]]) == True)[0]
                test2 = np.where(t_tot[use] < self.tMax)[0]

                #use = np.where(np.logical_or(test1,test2)==True)[0]
                #use = np.intersect1d(test1,test2)
                #use = np.intersect1d(np.intersect1d(test1,test2),use)
                use = use[np.intersect1d(test1,test2)]


                #For testing
                
                #print(self.tMax)
                #print(d_b[0])
                #print(t_b[0])
                #print(t_tot)
                #print(d_r[0])
                #print(q1[0])
                #print(q2[0])
                #print(q3[0])
                #print(rN[0])
                #print(d_r.T.shape)
                #print(rN.shape)
                #print(tsNext)
                #print(self.xyzSteps)
                #input()

                #record particle locations for particles that exceeded timestep
                use2 = np.where(t_tot > tsNext)[0] #particles that crossed a timestep
                use3 = np.intersect1d(use,use2) #crossed a timestep and still not crossed tMax, indexed to N_pts
                use4 = np.where(t_tot[use] > tsNext[use])[0] #crossed a timestep and still not crossed tMax, indexed to use
                if len(use3) > 0:
                    self.xyzSteps[i,use3,tIdxNext[use3],:] = q1[use4] + vec[use3]*frac[use3,None]
                    tIdxNext[use3]+=1
                    tsNext[use3] = ts[tIdxNext[use3]]

            #record the final timestep       
            self.xyzSteps[i,:,-1,:] = launchPt + vec*frac[:,None]

        #record the final intersectRecord
        self.intersectRecord = intersectRecord
        return
    
    def setupParallelVelocities(self, pdf = None):
        """
        split each source point into N macro-particles, 
        using the Energy distribution function, then translating energy into velocity
        

        """
        self.vSlices = np.ones((self.N_b*self.N_r*self.N_p, self.N_vS))*np.nan #velocity for each macroparticle
        self.energySlices = np.ones((self.N_b*self.N_r*self.N_p, self.N_vS))*np.nan #Energy in eV for each macroparticle
        self.energyFracs = np.zeros((self.N_b*self.N_r*self.N_p, self.N_vS)) # Fraction of the total energy from that starting point
        self.vBounds = np.zeros((self.N_b*self.N_r*self.N_p, self.N_vS+1))

        E = np.linspace(self.E_min, self.E_max, 10000)
        #If a pdf is not put in, uses an exponentially decaying profile
       
        if pdf is None:
            Epdf =  1/self.E_av * np.exp(- E/self.E_av)

        Ecdf = np.cumsum(Epdf[1:])*np.diff(E)
        Ecdf = np.insert(Ecdf, 0, 0)
        #create bspline interpolators for the cdf and cdf inverse
        inverseCDF = interp1d(Ecdf, E, kind='linear')
        forwardCDF = interp1d(E, Ecdf, kind='linear')
        cdfBounds = np.linspace(0,Ecdf[-1],self.N_vS+1)
        #space bins uniformly, then makes vSlices center of these bins in CDF space
        cdfSlices = np.diff(cdfBounds)/2.0 + cdfBounds[:-1]

        
        
        for i in range(self.N_b*self.N_r*self.N_p):
            # If you only want 1 velocity slice, sets that energy to the average enrgy
            if self.N_vS == 1:
                self.energySlices[i,:] = np.array([self.E_av]) 
            else:
                self.energySlices[i,:] = inverseCDF(cdfSlices)
            self.vSlices[i,:] = self.E2v_electron(self.energySlices[i,:])

            #energy fracs (for energy fluxes) (could probably do this outside of a loop)
            energyTotal = np.sum(self.energySlices[i,:])
            for j in range(self.N_vS):
                self.energyFracs[i,j] = self.energySlices[i,j] / energyTotal 
        return
            
    def test_vel_splitting(self, E_av, E_min, E_max, N_vs, N_b = 1, N_r = 1, N_p = 1):
        self.E_av = E_av
        self.E_min = E_min
        self.E_max = E_max
        self.N_b = N_b
        self.N_r = N_r
        self.N_p = N_p
        self.N_vS = N_vs
        self.setupParallelVelocities()
        print(self.energySlices/1e6)
        print(self.vSlices)
        print(self.energyFracs)
        
    def createSource(self, t:float, Btrace:np.ndarray):
        """
        creates a filament at time t

        t is timestep at which we calculate filament
        so t in gaussian is t-tMin

        Btrace is magnetic field line trace from MAFOT for filament ctr

        """
        self.gridPsiThetaDistAtCtr(self.rCtr, self.zCtr, multR=10.0, multZ = 10.0)
        #discretize filament into macroparticle sources along field line
        self.discretizeRunaways(self.N_r,self.N_p,self.N_b, Btrace, self.N_sig_r, self.N_sig_p, self.N_sig_b)
        self.gaussianAtPts(self.ctrPts, self.xyzPts, t-self.tMin, self.v_r)
        return
               
            
        
    def fluxSurfNorms(self, ep: object, R:np.ndarray, Z:np.ndarray):
        
        """
        Calculates vectors normal to poloidal flux surfaces at R,Z coordinates
        
        stolen from filament class

        returns normal vector coordinates in (R, Phi, Z)

        phiNorm is returned as 0
        """       
        if np.isscalar(R):                   
            R = np.array([R])
        if np.isscalar(Z):                   
            Z = np.array([Z])

        Br = ep.BRFunc.ev(R,Z)
        Bz = ep.BZFunc.ev(R,Z)
        Bt = ep.BtFunc.ev(R,Z)
        Bp = np.sqrt(Br**2+Bz**2)

        Bmag = np.sqrt(Bt**2 + Bp**2)


        bt = np.zeros((len(Bt), 3)) #R,t,Z
        bp = np.zeros((len(Bp), 3)) #R,t,Z
    
        #bt[:,1] = Bt / Bmag
        #bp[:,0] = Br / Bmag
        #bp[:,2] = Bz / Bmag

        bt[:,1] = -1.0
        bp[:,0] = Br / Bp
        bp[:,2] = Bz / Bp


        norms = np.cross(bt,bp)


        mag = np.linalg.norm(norms, axis=1)
        r_hat = np.zeros((len(norms), 3))
        r_hat = norms / mag[:, None]

        return r_hat
    
    
    def poloidalVectors(self, ep: object, R: np.ndarray ,Z: np.ndarray):
        """
        Calculates vectors tangent to poloidal flux surfaces in RZ plane at R,Z coordinates
        
        returns normal vector coordinates in (R, Phi, Z)

        phiNorm is returned as 0
        """

        if np.isscalar(R):                   
            R = np.array([R])
        if np.isscalar(Z):                   
            Z = np.array([Z])


        Br = ep.BRFunc.ev(R,Z)
        Bz = ep.BZFunc.ev(R,Z)
        Bt = ep.BtFunc.ev(R,Z)
        Bp = np.sqrt(Br**2+Bz**2)

        Bmag = np.sqrt(Bt**2 + Bp**2)

        bp = np.zeros((len(Bp), 3)) #R,phi,Z

        bp[:,0] = Br / Bp
        bp[:,2] = Bz / Bp

        return bp 
    

    
    def interpolateTrace(self, traceData:np.ndarray, N:int, addRawData=False):
        """
        given a trace of coordinates, interpolates to N discrete points
        """
        #calculate distance along wall
        dist = self.distance(traceData)
        
        # Resolution we need given the inputs
        resolution = (dist[-1] - dist[0])/float(N)
        print("Resolution: {:f} m".format(resolution))
        # Calculate how many total grid points we need based upon resolution and
        # total distance around curve/wall
        numpoints = int(np.ceil((dist[-1] - dist[0])/resolution))
        # Spline Interpolation (linear) - Make higher resolution wall.
        interpolator = interp.interp1d(dist, traceData, kind='slinear', axis=0)
        alpha = np.linspace(dist[0], dist[-1], numpoints+2)[1:-1]
        #add in raw data points so that we preserve the corners
        if addRawData == True:
            alpha = np.sort(np.hstack([alpha, dist]))
        
        interpolated_points = interpolator(alpha)
        arr, idx = np.unique(interpolated_points, axis=0, return_index=True)
        interpolated_points = arr[np.argsort(idx)]
     
        return interpolated_points, alpha
    
    def distance(self, traceData:np.ndarray):
        """
        Calculate distance along curve/wall (also called S) of ND points:
        """
        distance = np.cumsum(np.sqrt(np.sum(np.diff(traceData,axis=0)**2,axis=1)))
        distance = np.insert(distance, 0, 0)
        return distance

    def gridPsiThetaDistAtCtr(self, rCtr:float, zCtr:float, multR=10.0, multZ = 10.0):
        """
        generates an RZ grid, and then calculates distances [m] on that grid
        in the psi and theta directions from (self.rCtr, self.zCtr)

        these distances can then be used by gaussian

        multR is multiplier to create R width of grid: R width = multR*sigma_r
        multZ is multiplier to create Z height of grid: Z height = multZ*sigma_z
        """
        ep = self.ep
        sigma_r = self.sig_r
        sigma_p = self.sig_p

        #build rectilinear coordinates around the ctr
        r = np.linspace(rCtr - multR*sigma_r, rCtr + multR*sigma_r, 100)
        z = np.linspace(zCtr - multZ*sigma_p, zCtr + multZ*sigma_p, 100)
        R,Z = np.meshgrid(r,z)

        psiCtr, distPsi, thetaCtr, distTheta = self.fluxCoordDistance(rCtr,zCtr,R,Z)

        #set class variables
        self.distPsi = distPsi
        self.distTheta = distTheta
        self.rGrid = r
        self.zGrid = z
        self.psiCtr = psiCtr
        self.thetaCtr = thetaCtr
        self.gridShape = R.shape
        return


    def gaussian1D(self, B:float, x:np.ndarray, x0=0.0):
        """
        returns a 1D gaussian
        """
        g = np.sqrt(B/np.pi) * np.exp(-B*(x-x0)**2)
        return g
        
    def gaussianAtPts(self, ctrPts: np.ndarray, xyzPts: np.ndarray, t: float, v_r:float):
        """
        calculates a gaussian, centered at ctrPts, evaluated at xyzPts, at time t
        v_r is radial velocity

        saves gaussian values on 3D grid, as well as distances (psi, theta) used to
        evaluate the gaussian        

        calculates fractional density, n/n0, at each point
        """
        gaussian = np.zeros((xyzPts.shape[:-1]))
        dPsi = np.zeros((xyzPts.shape[:-1]))
        dTheta = np.zeros((xyzPts.shape[:-1]))
        dB = np.zeros((xyzPts.shape[:-1]))

        for i,ctr in enumerate(ctrPts):
            rCtr, zCtr, phiCtr = tools.xyz2cyl(ctr[0],ctr[1],ctr[2])
            xyz = xyzPts[i]
            R,Z,phi = tools.xyz2cyl(xyz[:,:,0],xyz[:,:,1],xyz[:,:,2])
            psiCtr, distPsi, thetaCtr, distTheta = self.fluxCoordDistance(rCtr,zCtr,R,Z)

            #2D gaussian (not used)
            #g = self.gaussian2D(distPsi,distTheta,self.sig_r,self.sig_p,t,v_r,self.E0)

            #3D gaussian
            g = self.gaussian3D(distPsi,distTheta,self.distB[i],self.sig_r,self.sig_p,self.sig_b,t,v_r)
            gaussian[i,:,:] = g
            dPsi[i,:,:] = distPsi
            dTheta[i,:,:] = distTheta
            dB[i,:,:] = self.distB[i]

        #now build weights for gaussian  
        int_r = np.zeros((self.r_pts.shape))
        rBounds = np.diff(self.r_pts) / 2.0 + self.r_pts[:-1]
        rBounds = np.insert(rBounds, 0, -10*self.sig_r)
        rBounds = np.insert(rBounds, len(rBounds), 10*self.sig_r)
        for i in range(self.N_r):
            lo = rBounds[i]
            hi = rBounds[i+1]
            int_r[i] = integrate.quad(self.f_r, lo, hi)[0]

        int_p = np.zeros((self.p_pts.shape))
        pBounds = np.diff(self.p_pts) / 2.0 + self.p_pts[:-1]
        pBounds = np.insert(pBounds, 0, -10*self.sig_p)
        pBounds = np.insert(pBounds, len(pBounds), 10*self.sig_p)
        for i in range(self.N_p):
            lo = pBounds[i]
            hi = pBounds[i+1]
            int_p[i] = integrate.quad(self.f_p, lo, hi)[0]

        int_b = np.zeros((self.b_pts.shape))
        bBounds = np.diff(self.b_pts) / 2.0 + self.b_pts[:-1]
        bBounds = np.insert(bBounds, 0, -10*self.sig_b)
        bBounds = np.insert(bBounds, len(bBounds), 10*self.sig_b)
        for i in range(self.N_b):
            lo = bBounds[i]
            hi = bBounds[i+1]
            int_b[i] = integrate.quad(self.f_b, lo, hi)[0]

        #now build time weight, slices are also bin boundaries here
        self.tSrc = self.ts[:self.N_src_t]
        self.f_t = lambda t: (1.0 / self.decay_t) * np.exp(-(t-self.tSrc[0]) / self.decay_t)
        int_t = np.zeros((self.N_src_t))
        tBounds = np.zeros((self.N_src_t+1))
        tBounds[:-1] = self.tSrc.copy()
        tBounds[-1] = self.tSrc[-1] + self.dt
        for i in range(self.N_src_t):
            lo = tBounds[i]
            hi = tBounds[i+1]
            int_t[i] = integrate.quad(self.f_t, lo, hi)[0]

        #scale weights so that sum=1 (all particles/energy gets exhausted)
        tFracSum = np.sum(int_t)
        tFrac = int_t / tFracSum

        #build 4D energy density hypercube weight matrix
        B,R,P,T = np.meshgrid(int_b,int_r,int_p,tFrac, indexing='ij')
        weights = B*R*P*T

        self.density = weights
        self.distPsi = dPsi
        self.distTheta = dTheta
        return

    def fluxCoordDistance(self, r0:float, z0:float, R:np.ndarray, Z:np.ndarray):
        """
        calculates euclidean distance along flux coordinate surfaces (psi and poloidal)
        from rCtr,zCtr on an R,Z grid. 

        rCtr and zCtr:  coordinate to calculate flux coordinates at / distance from
        R,Z: meshgrid around rCtr,zCtr

        returns psiCtr, distPsi, thetaCtr, distTheta
        """
        ep = self.ep

        #flux coordinates of ctr
        psiCtr = ep.psiFunc(r0, z0)
        thetaCtr = self.thetaFromRZ(ep, r0, z0)

        #==Find radial coordinates
        #Find coordinate transformation at ctr
        psiaxis = ep.g['psiAxis']
        psiedge = ep.g['psiSep']
        deltaPsi = np.abs(psiedge - psiaxis)
        Bp = ep.BpFunc(r0,z0)
        gradPsi = Bp*r0
        xfm = gradPsi / deltaPsi
        #convert distance in psi to euclidean coordinates
        psiRZ = ep.psiFunc.ev(R,Z)
        distPsi = (psiRZ - psiCtr) / xfm #on grid


        #===Find poloidal coordinates
        thetaRZ = self.thetaFromRZ(ep, R.flatten(), Z.flatten())
        #surfR, surfZ = ep.flux_surface(psiCtr, Npts, thetaRZ[sortIdx])
        surfR = ep.g['lcfs'][:,0]
        surfZ = ep.g['lcfs'][:,1]
        thetaSurf = self.thetaFromRZ(ep, surfR, surfZ)
        #calculate distance along flux surface
        surface = np.vstack((surfR, surfZ)).T
        dSurf = self.distance(surface)
        #interpolator that maps theta back to distance along the flux surface
        interpolator = interp.interp1d(thetaSurf, dSurf, kind='slinear', axis=0)
        dCtr = interpolator(thetaCtr)
        distTheta = interpolator(thetaRZ) - dCtr
        distTheta = distTheta.reshape(R.shape) #on grid

        return psiCtr, distPsi, thetaCtr, distTheta
        
    def thetaFromRZ(self, ep:object, R:np.ndarray, Z:np.ndarray):
        """
        calculates poloidal coordinate, theta, at R,Z given an equilibrium (ep) object

        returns theta
        """
        r_hat = self.fluxSurfNorms(ep, R, Z)
        r_hat = np.delete(r_hat, 1, 1)  # delete second column (phi)
        R_hat = np.array([1.0, 0.0])
        theta = np.arccos(np.dot(r_hat, R_hat))
        zMid = ep.g['ZmAxis']
        idx = np.where(Z < zMid)[0]
        theta[idx] *= -1.0
        return theta
    
    
    def discretizeRunaways(self, N_r: int, N_p: int, N_b: int, Btrace: np.ndarray, 
                           N_sig_r: int, N_sig_p: int, N_sig_b: int):
        """
        takes user input parameters and builds a filament 
        
        N_r: number of discrete points in radial (psi) direction
        N_p: number of discrete points in poloidal direction
        N_b: number of discrete points along B field line
        N_sig_r : # of self.sig_r widths for points to cover in +/- radial direction
        N_sig_p : # of self.sig_p widths for points to cover in +/- poloidal direction
        N_sig_0 : # of self.sig_0 widths for points to cover in +/- parallel direction
        
        stolen from filament class
        """
        d_b = self.distance(Btrace)
        #find index of filament center
        ctr = np.array([self.xCtr, self.yCtr, self.zCtr])
        idx = np.argmin(np.sum(np.abs(Btrace-ctr), axis=1))      
        d_b = d_b - d_b[idx]
        use = np.where(np.abs(d_b) < N_sig_b*self.sig_b)[0]

        #calculate center coordinates along field line
        self.ctrPts, alpha = self.interpolateTrace(Btrace[use], N_b, addRawData=False)
        self.distB = alpha - np.max(alpha) / 2.0

        #calculate points in parallel (B) direction
        bMax = N_sig_b * self.sig_b
        b_pts = np.linspace(-bMax, bMax, N_b)

        #calculate points in radial (psi) direction
        rMax = N_sig_r * self.sig_r
        r_pts = np.linspace(-rMax, rMax, N_r+2)[1:-1]

        #calculate points in poloidal direction
        pMax = N_sig_p * self.sig_p
        p_pts = np.linspace(-pMax, pMax, N_p+2)[1:-1]

        #get radial (psi) vectors for each ctrPt
        R = np.sqrt(self.ctrPts[:,0]**2+self.ctrPts[:,1]**2)
        Z = self.ctrPts[:,2]
        phi = np.arctan2(self.ctrPts[:,1], self.ctrPts[:,0])
        rVec = self.fluxSurfNorms(self.ep, R, Z)

        #convert from R,phi,Z to xyz
        rX, rY, rZ = tools.cyl2xyz(rVec[:,0], rVec[:,2], phi, degrees=False)
        rVecXYZ = np.vstack((rX,rY,rZ)).T
        rMag = np.linalg.norm(rVecXYZ, axis=1)
        rN = rVecXYZ / rMag[:,None]

        #get poloidal vectors for each ctrPt
        pVec = self.poloidalVectors(self.ep, R, Z)
        #convert from R,phi,Z to xyz
        pX, pY, pZ = tools.cyl2xyz(pVec[:,0], pVec[:,2], phi, degrees=False)
        pVecXYZ = np.vstack((pX,pY,pZ)).T
        pMag = np.linalg.norm(pVecXYZ, axis=1)
        pN = pVecXYZ / pMag[:,None]


        #create xyz coordinates for each filament source pt
        xyzPts = np.ones((N_b, N_r, N_p,3))*np.nan
        for i in range(N_b):
            for j in range(N_r):
                for k in range(N_p):
                    xyzPts[i,j,k,:] = self.ctrPts[i] + r_pts[j]*rN[i] - p_pts[k]*pN[i]

        self.xyzPts = xyzPts
        self.rN = rN
        self.pN = pN
        self.b_pts = b_pts
        self.p_pts = p_pts
        self.r_pts = r_pts

        #project user defined toroidal rotation velocity along the field line at center point
        
        '''
        Br = self.ep.BRFunc.ev(self.rCtr,self.zCtr)
        Bz = self.ep.BZFunc.ev(self.rCtr,self.zCtr)
        Bt = self.ep.BtFunc.ev(self.rCtr,self.zCtr)
        Bp = np.sqrt(Br**2+Bz**2)
        Bmag = np.sqrt(Bt**2 + Bp**2)
        self.v_rot_b = self.v_t * Bmag / Bt
        '''

        return
    
    def gaussian3D(self, dx:np.ndarray, dy:np.ndarray, dz:np.ndarray, 
                   sigX:float, sigY:float, sigZ:float, 
                   t:float , v_r:float, A=1.0):
        """
        calculates gaussian function with three spatial dimensions and 1 time dimension

        corresponds to density in Fundamenski advective-diffusive model

        dx,dy,dz are distances from center point
        sigX,sigY,sigZ are standard deviations in each direction
        t is timestep
        v_r is radial velocity (other velocity components not implemented, but could be)
        """
        A = 1.0

        B_x = 1.0 / (2*sigX**2)
        B_y = 1.0 / (2*sigY**2)
        B_z = 1.0 / (2*sigZ**2)
        #1D gaussians
        self.f_r = lambda x: self.gaussian1D(B_x, x)
        self.f_p = lambda y: self.gaussian1D(B_y, y)
        self.f_b = lambda z: self.gaussian1D(B_z, z)
        self.f_G = lambda x,y,z: self.f_r(x) * self.f_p(y) * self.f_b(z)

        pdfX = self.f_r(dx - t*v_r)
        pdfY = self.f_p(dy)
        pdfZ = self.f_b(dz)

        g = A * pdfX * pdfY * pdfZ

#        print("TEST=====!!!") 
#        intX = integrate.quad(f_x, 0, sigX*10)[0]
#        intY = integrate.quad(f_y, 0, sigY*10)[0]
#        intZ = integrate.quad(f_z, 0, sigZ*10)[0]
#        print(self.E0)
#        print(intX)
#        print(intY)
#        print(intZ)
#        print(intX * intY * intZ)
#        input()

        return g
    
    
    """
    Below are calculations for helical trajectories
    have not implimented yet, will do at some later date
    """

    def setupFreqs(self, B, v):
        """
        Calculates frequencies, periods, that are dependent upon B
        These definitions follow Freidberg Section 7.7.

        B is magnetic field magnitude
        Requires that v, B have the same length

        """
        gamma = self.calc_gamma(v)
        self.omegaGyro = self.e * B / (self.me / self.kg2eV) / gamma
        if np.isscalar(self.omegaGyro):
            self.omegaGyro = np.array([self.omegaGyro])

        self.fGyro = np.abs(self.omegaGyro)/(2*np.pi)
        self.TGyro = 1.0/self.fGyro
        return

    def setupRadius(self, vPerp):
        """
        calculates gyro radius using equations from Boozer PofP 2015
	
	
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
            self.rGyro[:,i] =  vPerp[i] / np.abs(self.omegaGyro)

        return

    def setupVelocities(self, N):
	#todo
        return

    def pullEqualProbabilityVelocities(self):
	#todo
        return

    def uniformGyroPhaseAngle(self):
        """
        Uniform sampling between 0 and 2pi

        returns angles in radians
        """
        self.gyroPhases = np.linspace(0,2*np.pi,self.N_gyroPhase+1)[:-1]
        return

    def uniformVelPhaseAngle(self):


        return


   
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

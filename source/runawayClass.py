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
import time
import shutil


#HEAT classes
import toolsClass
tools = toolsClass.tools()
import ioClass
IO = ioClass.IO_HEAT()
from rayTracerClass import shadowKernels

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

        Runaway Electrons Heat Flux Variables:
        -------------------------------------

        :id:  a unique string tag (or name) that is assigned to the RE.  example: fil1pt
        :tMin[s]: time in seconds of runaway birth. example: 1000e-6
        :tMax[s]: time in seconds to stop tracking runaways. example: 1500e-6
        :dt[s]: timestep size in seconds
        :decay_t[s]: decay constant for runaway birth energy.  The RE beam can be born over a 
          series of timesteps.  If N_src_t is > 1, this variable describes the exponential decay.
          At each of the birth timesteps, HEAT sources particles with a decaying total energy,
          which is prescribed by this exponential decay constant.  See function gaussianAtPts() 
          for more information.
        :N_src_t: number of birth timesteps over which we source particles.  if set to 1,
          particles are all born instantaneously.
        :rCtr[m]: radial coordinate in meters of runaway beam centroid at birth
        :zCtr[m]: vertical coordinate in meters of runaway beam centroid at birth
        :phiCtr[deg]: toroidal coordinate in degrees of runaway beam centroid at birth
        :sig_r[m]: Gaussian width of runaway beam in radial direction in meters
        :N_sig_r: Width of runaway beam radial Gaussian in units of sig_r
        :N_r: Number of discrete birth locations in radial direction
        :sig_p[m]: Gaussian width of runaway beam in poloidal direction in meters
        :N_sig_p: Width of runaway beam poloidal Gaussian in units of sig_p
        :N_p: Number of discrete birth locations in poloidal direction
        :sig_b[m]: Gaussian width of runaway beam in along field line in meters
        :N_sig_b: Width of runaway beam Gaussian along field line in units of sig_b
        :N_b: Number of discrete birth locations along field line
        :N_vS: Number of samples from the parallel velocity distribution function. If set to 1 all particles will have energy E_av
        :E_av[eV]: Average energy of runaway electrons. Tends to be in MeV range
        :E_max[eV]: Maximum energy used for distribution function
        :E_min[eV]: minimum energy used for distribution function
        :IRE[A]: Current of runaway beam. Used to determine the total energy of particles
        :I_dir: direction of current. Runaways will flow in opposite direction.
        :Pitch_angle: Pitch angle of particles
        :Drift: True or false, whether or not drifts are included in tracking of runaway electron trajectories

        

        """


        self.allowed_vars = [
                    'id',
                    'tMin[s]',
                    'tMax[s]', 
                    'dt[s]', 
                    'decay_t[s]', 
                    'N_src_t',
                    'rCtr[m]',
                    'zCtr[m]',
                    'phiCtr[deg]',
                    'sig_r[m]',
                    'N_r',
                    'sig_p[m]',
                    'N_p',
                    'sig_v[m]',
                    'N_b',
                    'N_vs',
                    
                    'E_av[eV]',
                    'E_min[eV]',
                    'E_max[eV]',
                    'IRE[A]',
                    'I_dir',
                    'Pitch_angle',
                    'Drift',
                            ]
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

        KE is in eV

        """
        return np.sqrt(self.c**2 * (1 - (self.me)**2/(RE_KE + self.me)**2))


    def calc_gamma(self, v):
        """
        Takes in the velocity and returns the relativistic gamma (Lorentz factor)
        necessary when calculating helical trajectories
        """
        
        return (1/np.sqrt(1 - v**2/self.c**2))
    
    
    def readREFile(self, path:str):
        """
        reads a runaway csv input file

        generates a data pd array
        """

        filFile = path + 'Runaways.csv'
        self.REData = pd.read_csv(filFile, sep=',', comment='#', skipinitialspace=True)
        print(self.REData.columns)
        print(self.REData.to_dict().keys())
        return

    """
    TO DO: Add intialize_REdist() (not from dict)? 
    """    
    
    def setupRETime(self):
        """
        sets up RE timesteps using data from RE file
        """
        try:
            df = self.REData      
            df['N_dt'] = (df['tMax[s]'] - df['tMin[s]']) / df['dt[s]']
            df['N_dt'] = df['N_dt'].round(0)
            df['N_dt'] = df['N_dt'].astype(int)
            df['tNew'] = df.apply( lambda x:  np.linspace(x['tMin[s]'],x['tMax[s]'],x['N_dt']+1), axis=1)
            self.tsFil = df['tNew'].to_numpy()

        except:
            print("\n\nError setting up RE time.  Check your RE file\n\n")
            log.info("\n\nError setting up RE time.  Check your RE file\n\n")
            sys.exit()
        return
    
    def initializeREdistFromDict(self, REDict:dict, id:int, ep:object):
        """
        initializes Runaway electron distribution from a dictionary (read from csv file in HEAT)

        dictionary is originally a pandas dataframe, so each row in the runaway.csv
        HEAT file is a new id # in the dict.  each row can therefore be accessed by:
        filDict['<parameter>'][<row#>]
        ep is equilParams object
        """
        #initialize this RE
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
            self.pitch = REDict['Pitch_angle'][id]
            self.drift = REDict['Drift'][id]
            
            # self.Uniform_vel todo
            self.v_r = 0 #backdoor thing, should  probably change this

            self.ep = ep

        except:
            print("\n\nCould not initialize runaways.  Check your input file!")
            log.info("\n\nCould not initialize runaways.  Check your input file!")
            sys.exit()

        return
    
    def calcTotalE(self):
        """
        Calculates the total energy of the RE beam including magnetic and kinetic energy
        """
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
        dpinit = getattr(MHD, 'dpinit', None) or 1.0
        MHD.getFieldpath(dphi, dpinit, self.gridfileStruct, self.controlfilePath, self.controlfileStruct, paraview_mask=False, tag=None, bbox=MHD.mafot_bbox)
        Btrace = tools.readStructOutput(self.structOutfile) 
        os.remove(self.structOutfile)
        return Btrace

    def findGuidingCenterPaths(self, pts:np.ndarray, MHD:object, traceDir:float, verbose = False):
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
                                self.controlfileStruct, bbox=MHD.mafot_bbox)

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
        Ray–mesh scene via rayTracerClass.shadowKernels (Open3D or Mitsuba).

        Engine sets runawayRayTracer from HF.rayTracer.
        """
        tri = np.stack([self.t1, self.t2, self.t3], axis=1)
        raw = getattr(self, 'runawayRayTracer', None) or 'open3d'
        raw = str(raw).strip().lower()
        self._re_rt = shadowKernels()
        try:
            if raw in ('mitsuba_gpu', 'mitsuba_cuda'):
                self._re_rt.buildMitsubaScene(tri, mitsubaMode='cuda')
                self._re_rt_engine = 'mitsuba'
                self._re_mitsuba_mode = 'cuda'
            elif raw == 'mitsuba_cpu':
                self._re_rt.buildMitsubaScene(tri, mitsubaMode='cpu')
                self._re_rt_engine = 'mitsuba'
                self._re_mitsuba_mode = 'cpu'
            else:
                self._re_rt.buildOpen3Dscene(tri.astype(np.float32))
                self._re_rt_engine = 'open3d'
        except Exception as e:
            #e.g. requested Mitsuba variant not available in this environment
            log.warning("RE Mitsuba scene build failed (%s); falling back to open3d", e)
            print("RE Mitsuba scene build failed ({}); falling back to open3d".format(e))
            self._re_rt = shadowKernels()
            self._re_rt.buildOpen3Dscene(tri.astype(np.float32))
            self._re_rt_engine = 'open3d'
        self.scene = getattr(self._re_rt, 'scene', None)
        return

    def traceREParticles(self, MHD: object, ts:np.ndarray , tIdx:int):
        """
        Traces RE macro-particles
        
        
        Based on the filament tracing code 
        
        
        Ray–mesh tests use rayTracerClass.shadowKernels; engine sets runawayRayTracer
        from HF.rayTracer.
        """
        pts = self.xyzPts.reshape(self.N_b*self.N_r*self.N_p, 3)
        N_pts = len(pts)
        N_ts = len(ts)

        print('Number of target faces: {:f}'.format(self.Nt))
        log.info('Number of target faces: {:f}'.format(self.Nt))
        print('Number of runaway source points: {:f}'.format(N_pts))
        log.info('Number of runaway source points: {:f}'.format(N_pts))

        #setup velocities + output arrays
        self.setupParallelVelocities()
        self.xyzSteps = np.zeros((self.N_vS, N_pts, N_ts, 3))
        self.buildIntersectionMesh()
        intersectRecord = np.ones((self.N_vS, N_pts, N_ts))*np.nan

        #MHD-equilibrium segment boundaries.  The RE flight can cross MHD timesteps; a single MAFOT
        #trace uses one g-file, so we trace one batched one-shot per segment and carry the end
        #positions into the next equilibrium (this replaces the old step-by-step equilibrium reload).
        MHDtimes = np.append(MHD.timesteps, self.tMax)
        Rctr = np.mean(np.sqrt(pts[:,0]**2 + pts[:,1]**2))
        dpinit = float(getattr(MHD, 'dpinit', 1.0))
        ep0, tEQ0 = self.ep, self.tEQ          #source equilibrium (restored per slice and at the end)

        for i in range(self.N_vS):
            print("\n---Tracing for velocity slice: {:f} [m/s]---\n".format(self.vSlices[0,i]))
            log.info("\n---Tracing for velocity slice: {:f} [m/s]---\n".format(self.vSlices[0,i]))

            self.xyzSteps[i,:,tIdx,:] = pts
            v_par = abs(float(self.vSlices[0,i]))          #relativistic parallel speed (sets the time scale)
            if v_par < 1e-12:
                continue
            speed = v_par
            traceDir = 1.0 if self.vSlices[0,i] > 0 else -1.0

            #--- MAFOT particle mode: MAFOT integrates the relativistic guiding-center drift (CPU and
            #    GPU).  PHYSICS MAPPING TO VALIDATE against the old Python d_z drift:
            #      sigma  -> drift sign (electrons move opposite I_dir),
            #      Ekin   -> slice kinetic energy [keV],
            #      lambda -> pitch / perp-energy fraction (~ self.pitch**2, self.pitch = v_perp/v_total).
            if self.drift:
                MHD.ParticleDirection = -int(np.sign(self.I_dir))
                MHD.Ekin              = float(self.energySlices[0,i]) / 1000.0
                MHD.Lambda            = float(self.pitch)**2
                MHD.Mass              = 1
                MHD.ParticleCharge    = -1
            else:
                MHD.ParticleDirection = 0                  #no drift -> pure field line
            MHD.v_par = 0.0; MHD.v_radial = 0.0; MHD.v_tor = 0.0   #prescribed drifts off

            launchPt = pts.copy()
            use = np.arange(N_pts)
            current_time = ts[tIdx]
            i_MHDtNext = 1
            self.ep, self.tEQ = ep0, tEQ0                  #restart each slice at the source equilibrium

            #--- walk the MHD-equilibrium segments; each is one batched one-shot trace ---
            while len(use) > 0 and current_time < self.tMax - 1e-15:
                seg_end = min(MHDtimes[i_MHDtNext], self.tMax)
                seg_dur = seg_end - current_time
                if seg_dur > 1e-15:
                    nSteps = int(np.ceil(seg_dur * v_par * 180.0 / (np.pi * Rctr) / dpinit * 1.1)) + 2
                    MHD.ittStruct = float(nSteps)
                    MHD.writeMAFOTpointfile(launchPt[use], self.gridfileStruct)
                    MHD.writeControlFile(self.controlfilePath+self.controlfileStruct, self.tEQ, traceDir, mode='struct')
                    MHD.getMultipleFieldPaths(1.0, self.gridfileStruct, self.controlfilePath, self.controlfileStruct, bbox=MHD.mafot_bbox)
                    full = np.genfromtxt(self.structOutfile, comments='#')
                    os.remove(self.structOutfile)
                    if full.ndim == 1:
                        full = full.reshape(1, -1)
                    xyzAll = full[:, 0:3]

                    #split per active marker on the exact launch point (Lc is unreliable on the CPU)
                    use_pts = launchPt[use]
                    nUse = len(use)
                    launchIdx = []
                    pos = 0
                    for k in range(nUse):
                        rel = np.where(np.abs(xyzAll[pos:pos+nSteps+2] - use_pts[k]).max(axis=1) < 1e-7)[0]
                        if len(rel) == 0:
                            rel = np.where(np.abs(xyzAll[pos:] - use_pts[k]).max(axis=1) < 1e-7)[0]
                        if len(rel) == 0:
                            raise ValueError("RE one-shot: could not locate launch point of active marker {:d}".format(k))
                        launchIdx.append(pos + int(rel[0])); pos = launchIdx[-1] + 1
                    if traceDir == 1.0:
                        bnd = launchIdx + [len(xyzAll)]
                        trajs = [xyzAll[bnd[k]:bnd[k+1]] for k in range(nUse)]
                    else:
                        bnd = [-1] + launchIdx
                        trajs = [xyzAll[bnd[k]+1:bnd[k+1]+1][::-1] for k in range(nUse)]

                    #batch ray-mesh intersection over every segment of every active marker
                    q1L, q2L, ownerL = [], [], []
                    for k in range(nUse):
                        tm = trajs[k]
                        if len(tm) >= 2:
                            q1L.append(tm[:-1]); q2L.append(tm[1:])
                            ownerL.append(np.full(len(tm)-1, k, dtype=np.int64))
                    if len(q1L) > 0:
                        q1 = np.concatenate(q1L); q2 = np.concatenate(q2L); owner = np.concatenate(ownerL)
                        d = q2 - q1; segMag = np.linalg.norm(d, axis=1)
                        with np.errstate(invalid='ignore', divide='ignore'):
                            rNorm = np.where(segMag[:,None] > 0, d/segMag[:,None], 0.0)
                        if getattr(self, '_re_rt_engine', 'open3d') == 'mitsuba':
                            prim = np.asarray(self._re_rt.intersect_segments_mitsuba(q1, q2, self._re_rt.scene, mitsubaMode=self._re_mitsuba_mode))
                            allHit = prim >= 0; allFace = prim.astype(np.int64); allDist = np.zeros(len(prim))
                        else:
                            hitMap, distMap = self._re_rt.cast_rays_open3d_segments(q1, rNorm, segMag)
                            allHit = (hitMap != 4294967295) & (distMap <= segMag + 1e-9)
                            allFace = hitMap.astype(np.int64); allDist = np.clip(distMap, 0.0, segMag)
                    else:
                        owner = np.zeros(0, dtype=np.int64)

                    #per active marker: record timesteps in this segment, carry the end position
                    newLaunch = launchPt.copy()
                    stillActive = np.ones(nUse, dtype=bool)
                    for k in range(nUse):
                        m = int(use[k])
                        tm = trajs[k]
                        if len(tm) < 2:
                            stillActive[k] = False             #left the domain immediately
                            continue
                        arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(tm, axis=0), axis=1))])
                        tRel = arc / speed                     #time since this segment's start
                        sel = np.where(owner == k)[0]
                        face = -1; tHit = tRel[-1]
                        if len(sel) > 0:
                            h = np.where(allHit[sel])[0]
                            if len(h) > 0:
                                kk = int(h[0]); tHit = (arc[kk] + allDist[sel[kk]]) / speed; face = int(allFace[sel[kk]])
                        tEndMarker = min(tHit, seg_dur)
                        #record positions at timesteps falling in (current_time, current_time+tEndMarker]
                        for j in range(tIdx, N_ts):
                            rel_t = ts[j] - current_time
                            if rel_t < 0:
                                continue
                            if rel_t > tEndMarker + 1e-12:
                                break
                            self.xyzSteps[i, m, j, :] = [np.interp(rel_t, tRel, tm[:,c]) for c in range(3)]
                        if face >= 0:
                            jhit = int(np.searchsorted(ts, current_time + tHit))
                            if jhit < N_ts:
                                intersectRecord[i, m, jhit] = face
                            stillActive[k] = False             #hit geometry -> stop tracing
                        else:
                            newLaunch[m] = [np.interp(seg_dur, tRel, tm[:,c]) for c in range(3)]   #carry to next equilibrium
                    launchPt = newLaunch
                    use = use[stillActive]

                current_time = seg_end
                if current_time < self.tMax - 1e-15:
                    self.tEQ = MHDtimes[i_MHDtNext]
                    self.ep  = MHD.ep[i_MHDtNext]
                    i_MHDtNext += 1

        self.intersectRecord = intersectRecord
        self.ep, self.tEQ = ep0, tEQ0          #restore the source equilibrium
        return
    
    def setupParallelVelocities(self, pdf = None):
        """
        split each source point into N macro-particles, 
        using the Energy distribution function, then translating energy into velocity
        

        """
        self.vSlices = np.ones((self.N_b*self.N_r*self.N_p, self.N_vS))*np.nan #velocity for each macroparticle
        self.gammas = np.ones((self.N_b*self.N_r*self.N_p, self.N_vS))*np.nan
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
            vtot = self.E2v_electron(self.energySlices[i,:])
            self.vSlices[i,:] = vtot *  np.sqrt(1 - self.pitch**2)
            self.gammas[i,:] = self.calc_gamma(vtot)


            #energy fracs (for energy fluxes) (could probably do this outside of a loop)
            energyTotal = np.sum(self.energySlices[i,:])
            for j in range(self.N_vS):
                self.energyFracs[i,j] = self.energySlices[i,j] / energyTotal 
        #print(self.vSlices)
        self.vPerps = self.pitch * self.vSlices / np.sqrt(1 - self.pitch **2)


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
        creates a RE beam at time t

        t is timestep at which we calculate RE beam
        so t in gaussian is t-tMin

        Btrace is magnetic field line trace from MAFOT for RE ctr

        """
        self.gridPsiThetaDistAtCtr(self.rCtr, self.zCtr, multR=10.0, multZ = 10.0)
        #discretize filament into macroparticle sources along field line
        # Using gaussian function for macroparticle energy TODO: update with more accurate density function
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

    def gridPsiThetaDistAtCtr(self, rCtr:float, zCtr:float, multR= 10.0, multZ = 10.0): 
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
        surfR = np.append(surfR, surfR[0])
        surfZ = np.append(surfZ, surfZ[0])
        thetaSurf = self.thetaFromRZ(ep, surfR, surfZ) #thetas at the LCFS
        #thetaSurf[-1] += 2 * np.pi
        #calculate distance along flux surface
        surface = np.vstack((surfR, surfZ)).T
        dSurf = self.distance(surface)
        
        #Fixing 2pi issue
        pi_ind = np.argmin(np.abs(thetaSurf - np.pi))
        negpi_ind = np.argmin(np.abs(thetaSurf + np.pi))
        thetaSurf = np.append(thetaSurf, np.array([thetaSurf[negpi_ind] + 2 * np.pi, thetaSurf[pi_ind] - 2 * np.pi]))
        dSurf = np.append(dSurf, np.array([dSurf[negpi_ind], dSurf[pi_ind]]))
        
        #interpolator that maps theta back to distance along the flux surface
        #interpolator = interp.interp1d(thetaSurf, dSurf, kind='slinear', axis=0)
        #dCtr = interpolator(thetaCtr)
        #distTheta = interpolator(thetaRZ) - dCtr
        
        dCtr = np.interp(thetaCtr, thetaSurf, dSurf)
        distTheta = np.interp(thetaRZ, thetaSurf, dSurf) - dCtr
        
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
        if type(Z) != np.ndarray:
            if Z < zMid:
                theta*=-1.0 
        else:
            if (Z < zMid).any():
                idx = np.where(Z < zMid)[0]
                theta[idx] *= -1.0
        return theta
    
    
    def discretizeRunaways(self, N_r: int, N_p: int, N_b: int, Btrace: np.ndarray, 
                           N_sig_r: int, N_sig_p: int, N_sig_b: int):
        """
        takes user input parameters and builds a RE beam 
        
        N_r: number of discrete points in radial (psi) direction
        N_p: number of discrete points in poloidal direction
        N_b: number of discrete points along B field line
        N_sig_r : # of self.sig_r widths for points to cover in +/- radial direction
        N_sig_p : # of self.sig_p widths for points to cover in +/- poloidal direction
        N_sig_0 : # of self.sig_0 widths for points to cover in +/- parallel direction
        
        stolen from filament class
        """
        d_b = self.distance(Btrace)
        #find index of RE center
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


        #create xyz coordinates for each RE source pt
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

#engineClass.py
#Description:   HEAT engine
#Engineer:      T Looby
#Date:          20220131
"""
HEAT Engine

Connects GUI/TUI to other HEAT classes.  Steps thru time solving for
HF, T, etc.
"""
import CADClass
import MHDClass
import toolsClass
import heatfluxClass
import openFOAMclass
import pfcClass
import gyroClass
import runawayClass
import filamentClass
import radClass
import ioClass
import plasma3DClass
import elmerClass
import time
import numpy as np
import logging
import os
import sys
import pandas as pd
import shutil
import errno
import copy
import EFIT.equilParams_class as EP
import GUIscripts.plotlyGUIplots as pgp
import trimesh
import multiprocessing

log = logging.getLogger(__name__)
tools = toolsClass.tools()

class engineObj():
    def __init__(self, rootDir, dataPath, OFbashrc, chmod, UID, GID, tsSigFigs=9, shotSigFigs=6):
        #number of significant figures after radix for timesteps
        self.tsSigFigs=tsSigFigs
        self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
        
        #number of significant figures for shot numbers
        self.shotSigFigs = shotSigFigs
        self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
        #file read/write/execute permissions
        self.chmod = chmod
        self.GID = GID
        self.UID = UID
        #where python source code is located (dashGUI.py)
        self.rootDir = rootDir
        #where we are saving data / HEAT output
        self.dataPath = dataPath
        #initialize all the HEAT python submodules (subclasses)
        self.initializeEveryone()
        #Set timestepMap to nothing
        self.timestepMap = None
        #Make a tmp dir that we will use for loading/unloading files from GUI
        self.makeTmpDir(dataPath)
        #initialize bashrc for OF
        self.OF.OFbashrc = OFbashrc
        #list of input files.  will stay None if using single file (see self.runHEAT)
        self.inputFileList = None
        return

    def makeTmpDir(self,dataPath):
        """
        makes a temp directory in rootDir path for user uploaded gfiles

        the self.tmpDir directory is accessible to the GUI users for uploading
        and downloading
        """
        tempDir = dataPath + '/tmpDir/'
        self.tmpDir = tempDir
        self.MHD.tmpDir = tempDir
        tools.makeDir(tempDir, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)
        return

    def machineSelect(self, MachFlag, machineList):
        """
        Select a machine and set the necessary paths
        """
        self.machineList = machineList
        self.MHD.machineList = machineList
        self.MachFlag = MachFlag
        self.MHD.MachFlag = MachFlag
        self.setInitialFiles()
        self.setHiddenInputs()
        return

    def initializeEveryone(self):
        """
        Create objects that we can reference later on
        """
        self.MHD = MHDClass.MHD(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.CAD = CADClass.CAD(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.HF = heatfluxClass.heatFlux(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.OF = openFOAMclass.OpenFOAM(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.GYRO = gyroClass.GYRO(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.RE = runawayClass.Runaways(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.RAD = radClass.RAD(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.FIL = filamentClass.filament(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)
        self.IO = ioClass.IO_HEAT(self.chmod, self.UID, self.GID)
        self.plasma3D = plasma3DClass.plasma3D()
        self.hf3D = plasma3DClass.heatflux3D()
        self.FEM = elmerClass.FEM(self.rootDir, self.dataPath, self.chmod, self.UID, self.GID)

        #set up class variables for each object
        self.MHD.allowed_class_vars()
        self.CAD.allowed_class_vars()
        self.HF.allowed_class_vars()
        self.OF.allowed_class_vars()
        self.GYRO.allowed_class_vars()
        self.RE.allowed_class_vars()
        self.RAD.allowed_class_vars()
        self.IO.allowed_class_vars()
        self.FIL.allowed_class_vars()
        self.FEM.allowed_class_vars()

        #setup number formats for each object
        self.MHD.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.CAD.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.HF.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.OF.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.GYRO.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.RE.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.RAD.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.IO.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.FIL.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.plasma3D.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.hf3D.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        self.FEM.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
        tools.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)

        return

    def refreshSubclasses(self):
        """
        re-initializes variables in subclasses
        """
        #initialize all the HEAT python submodules (subclasses)
        self.initializeEveryone()
        #select machine specific variables
        self.machineSelect(self.MachFlag, self.machineList)
        return


    def setInitialFiles(self):
        """
        sets files back to default settings
        infile is path to input file with all HEAT parameters
        PartsFile is path to file with parts we will calculate HF on
        IntersectFile is path to file with parts we will check for intersections on
        """
        self.infile = self.rootDir + '/inputs/default_input.csv'
        if self.MachFlag == 'nstx':
            print('Loading NSTX-U Input Filestream')
            log.info('Loading NSTX-U Input Filestream')
            self.CAD.machPath = self.dataPath + '/NSTX'
            self.OF.meshDir = self.dataPath + '/NSTX/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/NSTX/STLs/'
            self.CAD.STPpath = self.dataPath + '/NSTX/STPs/'

        elif self.MachFlag == 'st40':
            print('Loading ST40 Input Filestream')
            log.info('Loading ST40 Input Filestream')
            self.CAD.machPath = self.dataPath + '/ST40'
            self.OF.meshDir = self.dataPath + '/ST40/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/ST40/STLs/'
            self.CAD.STPpath = self.dataPath + '/ST40/STPs/'

        elif self.MachFlag == 'd3d':
            print('Loading DIII-D Input Filestream')
            log.info('Loading DIII-D Input Filestream')
            self.CAD.machPath = self.dataPath + '/D3D'
            self.OF.meshDir = self.dataPath + '/D3D/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/D3D/STLs/'
            self.CAD.STPpath = self.dataPath + '/D3D/STPs/'

        elif self.MachFlag == 'step':
            print('Loading STEP Input Filestream')
            log.info('Loading STEP Input Filestream')
            self.CAD.machPath = self.dataPath + '/STEP'
            self.OF.meshDir = self.dataPath + '/STEP/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/STEP/STLs/'
            self.CAD.STPpath = self.dataPath + '/STEP/STPs/'

        elif self.MachFlag == 'sparc':
            print('Loading SPARC Input Filestream')
            log.info('Loading SPARC Input Filestream')
            self.CAD.machPath = self.dataPath + '/SPARC'
            self.OF.meshDir = self.dataPath + '/SPARC/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/SPARC/STLs/'
            self.CAD.STPpath = self.dataPath + '/SPARC/STPs/'

        elif self.MachFlag == 'arc':
            print('Loading ARC Input Filestream')
            log.info('Loading ARC Input Filestream')
            self.CAD.machPath = self.dataPath + '/ARC'
            self.OF.meshDir = self.dataPath + '/ARC/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/ARC/STLs/'
            self.CAD.STPpath = self.dataPath + '/ARC/STPs/'

        elif self.MachFlag == 'cmod':
            print('Loading CMOD Input Filestream')
            log.info('Loading CMOD Input Filestream')
            self.CAD.machPath = self.dataPath + '/CMOD'
            self.OF.meshDir = self.dataPath + '/CMOD/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/CMOD/STLs/'
            self.CAD.STPpath = self.dataPath + '/CMOD/STPs/'

        elif self.MachFlag == 'west':
            print('Loading WEST Input Filestream')
            log.info('Loading WEST Input Filestream')
            self.CAD.machPath = self.dataPath + '/WEST'
            self.OF.meshDir = self.dataPath + '/WEST/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/WEST/STLs/'
            self.CAD.STPpath = self.dataPath + '/WEST/STPs/'

        elif self.MachFlag == 'kstar':
            print('Loading K-STAR Input Filestream')
            log.info('Loading K-STAR Input Filestream')
            self.CAD.machPath = self.dataPath + '/KSTAR'
            self.OF.meshDir = self.dataPath + '/KSTAR/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/KSTAR/STLs/'
            self.CAD.STPpath = self.dataPath + '/KSTAR/STPs/'

        elif self.MachFlag == 'aug':
            print('Loading AUG Input Filestream')
            log.info('Loading AUG Input Filestream')
            self.CAD.machPath = self.dataPath + '/AUG'
            self.OF.meshDir = self.dataPath + '/AUG/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/AUG/STLs/'
            self.CAD.STPpath = self.dataPath + '/AUG/STPs/'

        elif self.MachFlag == 'nstx':
            print('Loading NSTX Input Filestream')
            log.info('Loading NSTX Input Filestream')
            self.CAD.machPath = self.dataPath + '/NSTX'
            self.OF.meshDir = self.dataPath + '/NSTX/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/NSTX/STLs/'
            self.CAD.STPpath = self.dataPath + '/NSTX/STPs/'

        elif self.MachFlag == 'tcv':
            print('Loading TCV Input Filestream')
            log.info('Loading TCV Input Filestream')
            self.CAD.machPath = self.dataPath + '/TCV'
            self.OF.meshDir = self.dataPath + '/TCV/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/TCV/STLs/'
            self.CAD.STPpath = self.dataPath + '/TCV/STPs/'


        else:
            print("INVALID MACHINE SELECTION!  Defaulting to OTHER!")
            log.info("INVALID MACHINE SELECTION!  Defaulting to OTHER!")
            self.CAD.machPath = self.dataPath + '/OTHER'
            self.OF.meshDir = self.dataPath + '/OTHER/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/OTHER/STLs/'
            self.CAD.STPpath = self.dataPath + '/OTHER/STPs/'

        self.OF.templateCase = self.rootDir + '/openFoamTemplates/heatFoamTemplate'
        self.OF.templateDir = self.rootDir + '/openFoamTemplates/templateDicts'
        self.OF.materialDir = self.rootDir + '/openFoamTemplates/materials'

        return

    def setHiddenInputs(self, file = None):
        """
        sets up inputs based upon MachFlag that are needed for HEAT runs but
        are not in GUI
        """
        #MAFOT variables
        self.MHD.Nphi = 1
        self.MHD.ittLaminar = 10.0
        self.MHD.Nswall = 1
        self.MHD.phistart = 1 #(deg)
        self.MHD.PlasmaResponse = 0 #(0=no,>1=yes)
        self.MHD.Field = -1 #(-3=VMEC,-2=SIESTA,-1=gfile,M3DC1:0=Eq,1=I-coil,2=both)
        self.MHD.target = 0 #(0=useSwall)
        self.MHD.createPoints = 2 #(2=target)
        self.MHD.useFilament = 0 #(0=no)
        self.MHD.useTe_profile = 0 #(0=no)
        self.MHD.ParticleDirection = 0 #(1=co-pass,-1=ctr-pass,0=field-lines)
        self.MHD.ParticleCharge = 1 #(-1=electrons,>=1=ions)
        self.MHD.Ekin = 100.0 #[keV]
        self.MHD.Lambda = 0.1
        self.MHD.Mass = 2
        self.MHD.useM3DC1 = 0.0
        self.MHD.Smin = 0.0
        self.MHD.Smax = 5.0
        self.MHD.phimin = 0.0
        self.MHD.phimax = 6.28
        self.MHD.Zmin = -5.0
        self.MHD.Zmax = 5.0
        self.MHD.Rmin = 0.01
        self.MHD.Rmax = 13.0
        self.MHD.useECcoil = 0
        self.MHD.useIcoil = 0
        self.MHD.useCcoil = 0
        self.MHD.useFcoil = 0
        self.MHD.useBcoil = 0
        self.MHD.useBus = 0

        if self.MachFlag == 'nstx':
            self.CAD.permute_mask = True
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = True

        elif self.MachFlag == 'st40':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'd3d':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'step':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'sparc':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'arc':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'west':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'kstar':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'aug':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'tcv':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        elif self.MachFlag == 'cmod':
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        else:
            print("INVALID MACHINE SELECTION!  Defaulting to OTHER!")
            log.info("INVALID MACHINE SELECTION!  Defaulting to OTHER!")
            self.CAD.permute_mask = False
            self.CAD.unitConvert = 1.0
            self.CAD.assembly_mask = False

        return

    def setupTime(self, timesteps:np.ndarray, shot:int, tag=None, clobberFlag=True):
        """
        sets up timesteps.  timesteps includes timesteps from batchfile and from
        filament file.  

        creates self.timesteps, which has all timesteps (FIL.tsFil and MHD.timesteps)

        """
        if tag is None:
            tag = ''
        else:
            tag = '_'+tag

        if self.dataPath[-1]!='/':
            self.shotPath = self.dataPath + '/' + self.MachFlag +"_"+self.shotFmt.format(shot) +tag + "/"
        else:
            self.shotPath = self.dataPath + self.MachFlag +"_"+self.shotFmt.format(shot) +tag + "/"


        #make tree for this shot
        tools.makeDir(self.shotPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        #make tree for this timestep
        if self.shotPath[-1] != '/': self.shotPath += '/'

        #make file for each timestep
        for t in timesteps:
            #clobber and make time directory
            timeDir = self.shotPath + self.tsFmt.format(t) +'/'
            tools.makeDir(timeDir, clobberFlag=clobberFlag, mode=self.chmod, UID=self.UID, GID=self.GID)

        self.timesteps = timesteps
        self.shot = shot
        return


    def getMHDInputsForGUI(self,shot=None,traceLength=None,dpinit=None,
                     eqList=None,eqData=None,plasma3Dmask=None,
                     psiMult=1.0, BtMult=1.0, IpMult=1.0,
                     ):
        """
        Get the mhd inputs.  only used in GUI mode
        """
        tools.vars2None(self.MHD)
        tools.read_input_file(self.MHD, infile=self.infile)
        self.MHD.MachFlag=self.MachFlag #override input file with machine user selected
        self.MHD.setTypes()


        if shot is not None:
            self.MHD.shot = shot

        if traceLength is not None:
            self.MHD.traceLength = traceLength
            self.MHD.nTrace = int(traceLength / dpinit)
        if dpinit is not None:
            self.MHD.dpinit = dpinit
        self.MHD.eqList = eqList
        if eqList is not None:
            self.MHD.writeEQdata(eqList, eqData)
            self.MHD.eqList = eqList

        if plasma3Dmask is not None:
            self.plasma3D.plasma3Dmask = plasma3Dmask

        #multipliers for EQ variables
        self.MHD.psiMult = float(psiMult)
        self.MHD.BtMult = float(BtMult)
        self.MHD.IpMult = float(IpMult)

        #determine if the EQ are GEQDSKs or IMAS formatted JSON / NetCDF
        #and get the timesteps accordingly
        self.MHD.EQmode = self.MHD.determineEQFiletype(eqList[0]) #assume we are not mixing and matching EQ types
        if self.MHD.EQmode == 'geqdsk':
            self.timesteps = self.MHD.getGEQDSKtimesteps(eqList)
        elif self.MHD.EQmode == 'json':
            self.timesteps = np.array([])
            for eq in eqList:
                data = self.IO.readJSON(self.tmpDir + eq)
                self.timesteps = np.append( self.timesteps, np.round(np.array(data['equilibrium']['time']), 8) )

        #make tree branch for this shot
        self.setupTime(self.timesteps, shot, clobberFlag=False)
        self.MHD.shotPath = self.shotPath

        self.MHD.getGEQDSK(self.timesteps,eqList)
        if self.MHD.EQmode != 'geqdsk':
            self.MHD.makeEFITobjects()
        else:
            self.MHD.makeEFITobjects()
        self.NCPUs = multiprocessing.cpu_count() - 2 #reserve 2 cores for overhead
        self.MHD.psiSepLimiter = None

        self.MHD.setTypes()

        print('psiSep0 = {:f}'.format(self.MHD.ep[0].g['psiSep']))
        print('psiAxis0 = {:f}'.format(self.MHD.ep[0].g['psiAxis']))
        print('Nlcfs0: {:f}'.format(self.MHD.ep[0].g['Nlcfs']))
        print('length Rlcfs0: {:f}'.format(len(self.MHD.ep[0].g['lcfs'][:,0])))
        log.info('psiSep0 = {:f}'.format(self.MHD.ep[0].g['psiSep']))
        log.info('psiAxis0 = {:f}'.format(self.MHD.ep[0].g['psiAxis']))
        log.info('Nlcfs0: {:f}'.format(self.MHD.ep[0].g['Nlcfs']))
        if self.plasma3D.plasma3Dmask:
            print('Solving for 3D plasmas with MAFOT')
            log.info('Solving for 3D plasmas with MAFOT')
        else:
            print('Solving for 2D plasmas with EFIT (no MAFOT)')
            log.info('Solving for 2D plasmas with EFIT (no MAFOT)')
        return


    def gfileClean(self, psiRZMult,psiSepMult,psiAxisMult,FpolMult,
                   psiRZAdd,psiSepAdd,psiAxisAdd,FpolAdd,
                   Bt0Mult,Bt0Add,IpMult,IpAdd,t,correctAllts):
        """
        multiplies values in MHD ep object with scalars defined by user in html gui
        """
        print("psiRZ Multiplier = {:f}".format(psiRZMult))
        log.info("psiRZ Multiplier = {:f}".format(psiRZMult))
        print("psiRZ Addition = {:f}".format(psiRZAdd))
        log.info("psiRZ Addition = {:f}".format(psiRZAdd))

        print("psiSep Multipltier = {:f}".format(psiSepMult))
        log.info("psiSep Multipltier = {:f}".format(psiSepMult))
        print("psiSep Addition = {:f}".format(psiSepAdd))
        log.info("psiSep Addition = {:f}".format(psiSepAdd))

        print("psiAxis Multipltier = {:f}".format(psiAxisMult))
        log.info("psiAxis Multipltier = {:f}".format(psiAxisMult))
        print("psiAxis Addition = {:f}".format(psiAxisAdd))
        log.info("psiAxis Addition = {:f}".format(psiAxisAdd))

        print("Fpol Multiplier = {:f}".format(FpolMult))
        log.info("Fpol Multiplier = {:f}".format(FpolMult))
        print("Fpol Addition = {:f}".format(FpolAdd))
        log.info("Fpol Addition = {:f}".format(FpolAdd))

        print("Bt0 Multiplier = {:f}".format(Bt0Mult))
        log.info("Bt0 Multiplier = {:f}".format(Bt0Mult))
        print("Bt0 Addition = {:f}".format(Bt0Add))
        log.info("Bt0 Addition = {:f}".format(Bt0Add))

        print("Ip Multiplier = {:f}".format(IpMult))
        log.info("Ip Multiplier = {:f}".format(IpMult))
        print("Ip Addition = {:f}".format(IpAdd))
        log.info("Ip Addition = {:f}".format(IpAdd))

        #all timesteps or just the timestep viewed in GUI?
        if correctAllts is None:
            allMask = False
            idx = np.where(t==self.MHD.timesteps)[0]
        elif 'all' in correctAllts:
            allMask = True
            idx = self.MHD.timesteps

        print("Apply changes to all timesteps?: "+str(allMask))
        log.info("Apply changes to all timesteps?: "+str(allMask))


        for i in range(len(idx)):
            self.MHD.ep[i].g['psiRZ'] *= psiRZMult
            self.MHD.ep[i].g['psiSep'] *= psiSepMult
            self.MHD.ep[i].g['psiAxis'] *= psiAxisMult
            self.MHD.ep[i].g['Fpol'] *= FpolMult
            self.MHD.ep[i].g['Bt0'] *= Bt0Mult
            self.MHD.ep[i].g['Ip'] *= IpMult

            self.MHD.ep[i].g['psiRZ'] += psiRZAdd
            self.MHD.ep[i].g['psiSep'] += psiSepAdd
            self.MHD.ep[i].g['psiAxis'] += psiAxisAdd
            self.MHD.ep[i].g['Fpol'] += FpolAdd
            self.MHD.ep[i].g['Bt0'] += Bt0Add
            self.MHD.ep[i].g['Ip'] += IpAdd

            psi = self.MHD.ep[i].g['psiRZ']
            psiSep = self.MHD.ep[i].g['psiSep']
            psiAxis = self.MHD.ep[i].g['psiAxis']
            self.MHD.ep[i].g['psiRZn'] = (psi - psiAxis) / (psiSep - psiAxis)
        return


    def findPsiSepfromEQ(self,t, rNew=None):
        """
        finds psiSep by stepping to/from core and calculating
        minimum psiN along Z-plane at each R location.  Increments in um
        """
        tIdx = np.where(t==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[tIdx]
        #redefine LCFS to be tangent to CAD maximum R (because rNew=None)
        self.newLCFS(t, rNew=rNew, zNew=None, psiSep=None)
        print("CAD rTangent: {:f}".format(self.MHD.rTangent))
        if rNew is None:
            rSep = self.MHD.rTangent
        else:
            rSep = float(rNew)

        zMin = ep.g['ZmAxis'] - 0.25
        zMax = ep.g['ZmAxis'] + 0.25
        zWall = np.linspace(zMin, zMax, 100000)
        while(True):
            print("rSep = {:f}".format(rSep))
            rWall = np.ones((len(zWall)))*rSep
            psiN = ep.psiFunc.ev(rWall, zWall).min()
            print("psiN Minimum = {:f}".format(psiN))
            if psiN < 1.0:
                rSep -= 1e-6 #step 1um away from core

            else:
                break
        self.MHD.rTangent = rSep
        #now that we have found rTangent and the psiSep that makes this CAD
        #truly limited, write this psiSep to all gfiles in HEAT tree
        self.newLCFSallTimesteps(rNew=self.MHD.rTangent, zNew=None, psiSep=None)
        #self.MHD.makeEFITobjects()
        for PFC in self.PFCs:
            PFC.resetPFCeps(self.MHD)
        
        return

    def findPsiSepfromPFCs(self, t, rNew=None):
        """
        finds psiSep for limiters by incrementally increasing R_{psiSep, IMP},
        at a specific time.  Then rewrite all gfiles in MHD.timesteps with this
        new LCFS.

        Both MHD and PFC objects must be defined before running this function
        """
        tIdx = np.where(t==self.MHD.timesteps)[0][0]
        gfile = self.MHD.shotPath + self.tsFmt.format(t) + '/g'+self.shotFmt.format(self.MHD.shot)+'_'+self.tsFmt.format(t)
        #redefine LCFS to be tangent to CAD maximum R (because rNew=None)
        self.newLCFS(t, rNew=rNew, zNew=None, psiSep=None)
        print("CAD rTangent: {:f}".format(self.MHD.rTangent))

        #run 2D equilibrium for all points in PFCs and determine if this rTangent
        #actually resulted in all psiN > 1.0.  If not, add 1um to rTangent and
        #test it again, looping until all points in CAD are in SOL
        while(True):
            privateN = 0
            for i,PFC in enumerate(self.PFCs):
                PFC.shadowed_mask = np.zeros((len(PFC.centers)))
                PFC.ep = EP.equilParams(gfile)
                self.MHD.psi2DfromEQ(PFC)
                psiMinimum = min(PFC.psimin)
                if psiMinimum < 1.0:
                    print(PFC.ep.g['psiSep'])
                    print(psiMinimum)
                    privateN += 1
                #psiIdx = np.argmin(PFC.psimin)
                #x[i] = PFC.centers[psiIdx][0]
                #y[i] = PFC.centers[psiIdx][1]
                #z[i] = PFC.centers[psiIdx][2]
                #R,Z,phi = tools.xyz2cyl(x,y,z)
            if privateN > 0:
                print("Number of PFCs with points in core: {:d}".format(privateN))
                print("psiSep: {:f}".format(self.MHD.psiSepLimiter))
                print("Incrementing rTangent by 10 um...")
                self.MHD.rTangent += 1e-6 #add 1 um
                print("New rTangent: {:f}".format(self.MHD.rTangent))
                self.newLCFS(t, rNew=self.MHD.rTangent, zNew=None, psiSep=None)

            else:
                break

        #now that we have found rTangent and the psiSep that makes this CAD
        #truly limited, write this psiSep to all gfiles in HEAT tree
        self.newLCFSallTimesteps(rNew=self.MHD.rTangent, zNew=None, psiSep=None)
        #self.MHD.makeEFITobjects()
        return

    def newLCFSallTimesteps(self, rNew=None, zNew=None, psiSep=None):
        """
        Loops through the timesteps in MHD object and overwrites new gfiles
        with user defined psiSep
        """
        for t in self.MHD.timesteps:
            self.newLCFS(t,rNew,zNew,psiSep)
        return

    def newLCFS(self, t, rNew=None, zNew=None, psiSep=None):
        """
        resets the lcfs so that it is defined as psi value at rminNew, zminNew,
        or psiSep defines new LCFS. For use with limited discharges

        overwrites existing gfile
        """
        print("redefining LCFS")
        log.info("redefining LCFS")
        idx = np.where(t==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[idx]

        if zNew in ['None', 'none', 'na', 'NA', 'N/A', 'n/a', '', ' ', None]:
            zNew = None
        else:
            zNew = float(zNew)

        #find maximum value of R in all of the PFCs
        if rNew in ['None', 'none', 'na', 'NA', 'N/A', 'n/a', '', ' ', None]:
            print('Using CAD rNew tangent point')
            rNew = 0.0
            for PFC in self.PFCs:
                if PFC.Rmax > rNew:
                    rNew = PFC.Rmax
        else:
            rNew = float(rNew)

        g = self.MHD.renormalizeLCFS(self.MHD.ep[idx], rNew, zNew, psiSep)
        ep.g['psiSep'] = g['psiSep']
        ep.g['lcfs'] = g['lcfs']
        ep.g['Nlcfs'] = g['Nlcfs']

        #set psi and R for LCFS tangent point
        self.MHD.psiSepLimiter = g['psiSep']
        self.MHD.rTangent = rNew

        #overwrite existing gfile
        gfile = self.MHD.shotPath + self.tsFmt.format(t) + '/g'+self.shotFmt.format(self.MHD.shot)+'_'+self.tsFmt.format(t)
        self.MHD.writeGfile(gfile, shot=self.MHD.shot, time=t, ep=ep)
        self.MHD.ep[idx] = EP.equilParams(gfile)
        for PFC in self.PFCs:
            PFC.resetPFCeps(self.MHD)
        return

    def writeGfile(self, newGfile=None, shot=None, t=None):
        """
        writes a new gfile from EP object already in MHD object
        """
        idx = np.where(t==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[idx]
        if os.path.isabs(newGfile) is False:
            #save to tmpDir if path is not absolute
            newGfile = self.tmpDir + newGfile
        else:
            print("Please enter a filename that is not absolute (no directories)")
            log.info("Please enter a filename that is not absolute (no directories)")
            return
        print("Writing new gFile: " + newGfile)
        log.info("Writing new gFile: " + newGfile)
        self.MHD.writeGfile(newGfile, shot, t, ep)
        return

    def createGfileZip(self):
        """
        creates a zip file with all the GEQDSKs in self.MHD.timesteps
        """
        sendBack = self.tmpDir + 'savedGeqdsks.zip'
        #stash files in a subdir in tmpDir
        dir = self.tmpDir+'geqdsks/'
        tools.makeDir(dir)

        #write all new geqdsks to tmpDir
        for i,t in enumerate(self.MHD.timesteps):
            f = dir + 'g'+self.shotFmt.format(self.MHD.shot) +'_'+ self.tsFmt.format(t)
            self.MHD.writeGfile(f, self.MHD.shot, t, self.MHD.ep[i])

        #combine geqdsks into zip file for user download
        shutil.make_archive(self.tmpDir + 'savedGeqdsks', 'zip', dir)

        return sendBack

    def interpolateGfile(self, t):
        """
        finds values of all gfile parameters for user defined time using the
        gfiles in the self.MHD.ep object

        interpolates at given t then writes to file in same gFileInterpolate
        directory

        returns the name of the new gFile
        """
        print("Interpolating gFile")
        log.info("Interpolating gFile")
        t = int(t)
        ep = self.MHD.gFileInterpolate(t)
        g = 'g'+self.shotFmt.format(self.MHD.shot) +'_'+ self.tsFmt.format(t)
        gFileName = self.tmpDir + g
        self.MHD.writeGfile(gFileName,self.MHD.shot,t,ep)
        print("gFile Interpolated")
        log.info("gFile Interpolated")
        return gFileName

    def interpolateNsteps(self, gfiles, timesteps, N):
        """
        interpolates N steps between gfiles arranged at user defined timesteps
        Saves resulting gfiles into tmpDir, then creates zip file containing
        them all

        timesteps and gfiles should be sorted so that they are increasing in
        chronological order

        gfiles should be named following d3d convention by Lau:
        g<XXXXXX>.<YYYYY>
        where <XXXXXX> is shot number and <YYYYY> is timestep[ms]
        """
        #change filename to HEAT GEQDSK naming convention
        self.MHD.tmax = max(timesteps)
        self.MHD.tmin = min(timesteps)
        shot = self.MHD.shot
        newGfiles = []
        for i,f in enumerate(gfiles):
            g = 'g'+self.shotFmt.format(shot) +'_'+ self.tsFmt.format(timesteps[i])
            newGfiles.append(g)
            old = self.tmpDir + gfiles[i]
            new = self.tmpDir + newGfiles[i]
            try:
                shutil.copyfile(old, new)
            except:
                print("Could not copy timestep {:d}.  Skipping.".format(timesteps[i]))
                log.info("Could not copy timestep {:d}.  Skipping.".format(timesteps[i]))


        #rebuild eq objects
        self.MHD.getGEQDSK(self.MachFlag,newGfiles)
        self.MHD.makeEFITobjects()


        #interpolate between existing gfiles
        newNames = []
        nTime = len(timesteps)
        for i in range(nTime-1):
            times = np.linspace(timesteps[i],timesteps[i+1],N)
            for t in times:
                newName = self.interpolateGfile(t)
                newNames.append(newName)

        #now zip all these new gFiles into a single file that the user may
        #download from GUI
        from zipfile import ZipFile
        from os.path import basename
        zipFile = self.tmpDir + 'InterpolatedGfiles.zip'
        zipObj = ZipFile(zipFile, 'w')
        for f in newNames:
            zipObj.write(f, basename(f))
        zipObj.close()

        return

    def getIOInputs(self, list=None):
        """
        Loads input/output filestream inputs

        if list is None (terminal mode), load IO parameters from file, else
        list should be a list which corresponds to IO.allowed_vars
        """
        if list == None:
            #file saving
            tools.initializeInput(self.IO, infile=self.infile)
        else:
            self.IO.outputMasks(list)

        return


    def getCADinputs(self,gridRes=None,gTx=None,gTy=None,gTz=None, mode='tui'):
        """
        Loads CAD inputs
        """
        if mode == 'gui':
            if gridRes is not None:
                #check if intersection grid resolution string is a number,
                #if not use standard mesh algorithms
                if tools.is_number(gridRes):
                    self.CAD.gridRes = gridRes
                else:
                    self.CAD.gridRes = "standard"

            self.CAD.xT = gTx
            self.CAD.yT = gTy
            self.CAD.zT = gTz
        else:
            tools.initializeInput(self.CAD, infile=self.infile)
        
        return

    def getCAD(self,STPfile=None,STPdata=None, ts=None):
        """
        Loads CAD file for GUI users
        """
        falseList = [False, 'F', 'f', 'false', 'False', 'FALSE']

        if hasattr(self.CAD, 'gridRes'):
            pass
        else:
            tools.initializeInput(self.CAD, infile=self.infile)
        self.CAD.rootDir = self.rootDir #set HEAT rootDir

        #make machine specific CAD path with correct permissions
        tools.makeDir(self.CAD.machPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)

        if STPfile is not None:
            #make STP path if it doesnt exist
            tools.makeDir(self.CAD.STPpath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
            newSTPpath = self.CAD.STPpath + STPfile
            #check to see if this STP file exists and write data to the file
            if os.path.isfile(newSTPpath) == False:
                print("New STP file.  Writing")
                with open(newSTPpath, 'wb') as f:
                    f.write(STPdata)
                try:
                    atime = os.stat(newSTPpath).st_atime
                    os.utime(newSTPpath, (atime, ts))
                except OSError as e:
                    print(f"Could not set file timestamps: {e}")
                self.CAD.overWriteMask = True #we need to also overwrite meshes
                try:
                    os.chmod(newSTPpath, self.chmod)
                    os.chown(newSTPpath, self.UID, self.GID)
                except OSError as e:
                    print(f"Could not set file permissions: {e}")
            else:
                #if file was modified, overwrite
                if ts != os.stat(newSTPpath).st_mtime:
                    print("File was modified since last HEAT upload.  Overwriting...")
                    with open(newSTPpath, 'wb') as f:
                        f.write(STPdata)
                    try:
                        atime = os.stat(newSTPpath).st_atime
                        os.utime(newSTPpath, (atime, ts))
                    except OSError as e:
                        print(f"Could not set file timestamps: {e}")
                    self.CAD.overWriteMask = True #we need to also overwrite meshes
                    try:
                        os.chmod(newSTPpath, self.chmod)
                        os.chown(newSTPpath, self.UID, self.GID)
                    except OSError as e:
                        print(f"Could not set file permissions: {e}")
                #meshes are already up to date, only overwrite if user requests
                else:
                    if self.CAD.overWriteMask in falseList:
                        self.CAD.overWriteMask = False
                        print("STP file is already in the HEAT database.  Not overwriting...")
                    else:
                        self.CAD.overWriteMask = True
                        print("Overwriting CAD data per user input file overWriteMask variable.")

            self.CAD.STPfile = newSTPpath
            #load STP file using FreeCAD
            self.CAD.loadSTEP()

        return

    def getCADfromTUI(self,CADfile):
        """
        Loads CAD file for terminal users.  Here we call this file an STPfile,
        but it could be other formats (ie BREP, FCStd, IGES, etc.)
        """
        tools.makeDir(self.CAD.STPpath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        #if CAD file is set to None, do not load
        #(this is true when user brings their own meshes and no STEP file)
        if type(CADfile) != str:
            print("CADfile column in batchFile set to None.  Skipping CAD load.")
            print("This is ok if you are bringing your own mesh (BYOM)")
            log.info("CADfile column in batchFile set to None.  Skipping CAD load.")
            log.info("This is ok if you are bringing your own mesh (BYOM)")

        else:
            #get file name
            stpName = os.path.basename(CADfile)
            #time last modified
            mtime_orig = os.stat(CADfile).st_mtime
            #time last read
            atime_orig = os.stat(CADfile).st_atime
            #we will copy to this cache directory for future use
            newSTPpath = self.CAD.STPpath + stpName
            #check to see if this STP file exists and write data to the file
            if os.path.isfile(newSTPpath) == False:
                print("New STP file.  Writing")
                shutil.copyfile(CADfile, newSTPpath)
                #set modified timestamps to match original
                try:
                    os.utime(newSTPpath, (atime_orig, mtime_orig))
                except OSError as e:
                    print(f"Could not set file timestamps: {e}")
                self.CAD.overWriteMask = True #we need to also overwrite meshes
            else:
                #time last modified
                mtime_new = os.stat(newSTPpath).st_mtime
                #time last read
                atime_new = os.stat(newSTPpath).st_atime
                #if file was modified, overwrite
                if mtime_orig != mtime_new:
                    print("File was modified since last HEAT upload.  Overwriting...")
                    shutil.copyfile(CADfile, newSTPpath)
                    print(atime_orig)
                    print(mtime_orig)
                    try:
                        os.utime(newSTPpath, (atime_orig, mtime_orig))
                    except OSError as e:
                        print(f"Could not set file timestamps: {e}")
                    self.CAD.overWriteMask = True #we need to also overwrite meshes
                    try:
                        os.chmod(newSTPpath, self.chmod)
                        os.chown(newSTPpath, self.UID, self.GID)
                    except OSError as e:
                        print(f"Could not set file permissions: {e}")
                else:
                    falseList = [False, 'F', 'f', 'false', 'False', 'FALSE']
                    if self.CAD.overWriteMask in falseList:
                        self.CAD.overWriteMask = False
                        print("STP file is already in the HEAT database.  Not overwriting...")
                    else:
                        self.CAD.overWriteMask = True
                        print("Overwriting CAD data per user input file overWriteMask variable.")
            self.CAD.STPfile = newSTPpath
            #load STP file using FreeCAD
            self.CAD.loadSTEP()
        return

    def getPFCdataFromGUI(self, data):
        """
        initializes timestepMap from GUI rather than from file
        """
        self.timestepMap = pd.DataFrame.from_dict(data)[list (data[0].keys())]
        self.timestepMap = self.timestepMap.rename(columns=lambda x: x.strip())
        self.timestepMap['PFCname'] = self.timestepMap['PFCname'].str.strip()
        self.timestepMap['intersectName'] = self.timestepMap['intersectName'].str.strip()
        self.timestepMap['excludeName'] = self.timestepMap['excludeName'].str.strip()
        self.timestepMap['DivCode'] = self.timestepMap['DivCode'].str.strip()
        return

    def readPFCfile(self, infile):
        """
        The HEAT PFC File:
        ------------------

        The PFC file defines which CAD objects comprise the region of interest (ROI),
        as well as various parameters for each ROI object.  The PFC file is a CSV
        file in which each row corresponds to a separate ROI object.

        A user can either use the HEAT algorithms to read a parametric CAD file (ie .step)
        and generate meshes or Bring Your Own Mesh (BYOM).  When using the HEAT algorithms 
        for meshing, the PFC file contains information about the part objects in the CAD file.  
        When using BYOM, the PFC file contains the name of the user's mesh file, and the 
        resolution is ignored.  The columns in the PFC file are as follows:

        :timesteps: the timesteps during which we should calculate quantities on this
          ROI object
        :PFCname: the name of the CAD object as it appears in the CAD file.  If the 
          user Brings Your Own Mesh (BYOM), then the PFC name should be the location
          name of the mesh (.stl) file.
        :resolution: the maximum length [mm] of any triangular mesh element for this
          ROI object.  This is a proxy for the resolution.  If the user Brings Your
          Own Mesh (BYOM), then this parameter is ignored.
        :DivCode: divertor code.  This can be: LO, LI, UO, UI, which correspond to:
          Lower Outer, Lower Inner, Upper Outer, Upper Inner.  These codes
          are how each PFC in the ROI get flagged as belonging to a specific
          divertor, and will affect their power later (Psol * powerFrac).  The fractions
          for each region are defined in the X_input.csv file
        :intersectName: name of the PFCs that may cast a magnetic shadow upon the ROI object
          we are calculating power on.  If we launch field lines from the
          PFCname part and follow them up the SOL, we may hit one of these
          intersectName PFCs.  If the user is unsure of this value, "all"
          can be specified to check against all parts in the CAD.
          Multiple PFCs can be specified by using ":" between part names.
        :excludeName: name of PFCs to exclude in the intersection check.  This can be
         useful when we use the "all" switch for intersectName and want to exclude 
         some obvious PFCs (for example a gyroSourcePlane)


        """
        if infile == None:
            print("No timesteps input file.  Please provide input file")
            log.info("No timesteps input file.  Please provide input file")
            sys.exit()

        self.timestepMap = pd.read_csv(infile, comment='#')
        self.timestepMap = self.timestepMap.rename(columns=lambda x: x.strip())
        self.timestepMap['PFCname'] = self.timestepMap['PFCname'].str.strip()
        self.timestepMap['intersectName'] = self.timestepMap['intersectName'].str.strip()
        self.timestepMap['excludeName'] = self.timestepMap['excludeName'].str.strip()
        self.timestepMap['DivCode'] = self.timestepMap['DivCode'].str.strip()
        return

    def getPFCinputs(self, defaultMask=True):
        """
        Load heat flux PFCs and intersection calculation PFCs from input files
        Generate meshes for each
        Generate file for HEAT calculation for each heat flux PFC
        """
        #if defaultMask is true read default pfc file
        if defaultMask is True:
            self.readPFCfile(self.pfcFile)
        self.CAD.getROI(self.timestepMap)
        self.CAD.getROImeshes()
        self.CAD.writeMesh2file(self.CAD.ROImeshes, self.CAD.ROI, self.CAD.ROIresolutions, path=self.CAD.STLpath)
        print("Calculating HF on these tiles:")
        log.info("Calculating HF on these tiles:")
        print(self.CAD.ROIList)
        log.info(self.CAD.ROIList)

        #Find potential intersections by file as they correspond to ROI PFCs,
        # then mesh them using FreeCAD Standard mesh algorithm
        # or use user supplied meshes (BYOM=True)
        self.CAD.getIntersectsFromFile(self.timestepMap)
        self.CAD.getIntersectMeshes(resolution=self.CAD.gridRes)
        self.CAD.writeMesh2file(self.CAD.intersectMeshes,
                                self.CAD.intersectList,
                                path=self.CAD.STLpath,
                                resolution=self.CAD.gridRes
                                )
        print("All potential intersects on these tiles:")
        log.info("All potential intersects on these tiles:")
        print(self.CAD.intersectList)
        log.info(self.CAD.intersectList)
        print(self.MHD.timesteps)

        self.PFCs = []
        #initialize PFC objects for each ROI part
        for i,row in self.timestepMap.iterrows():
            PFC = pfcClass.PFC(row, self.rootDir, self.dataPath, self.CAD.intersectList, self.timesteps)
            PFC.setupNumberFormats(self.tsSigFigs, self.shotSigFigs)
            PFC.makePFC(self.MHD, self.CAD, i, clobberFlag=True)
            self.PFCs.append(PFC)

        #Build HEAT file tree
        tools.buildDirectories(self.CAD.ROIList,
                               self.MHD.timesteps,
                               self.MHD.shotPath,
                               clobberFlag=True,
                               chmod=self.chmod,
                               UID=self.UID,
                               GID=self.GID
                              )

        for PFC in self.PFCs:
            ctrs = PFC.centers
            R,Z,phi = tools.xyz2cyl(ctrs[:,0],ctrs[:,1],ctrs[:,2])
            PFC.Rmax = max(R)
            print("PFC: "+PFC.name+" Maximum R: {:f}".format(max(R)))
        return

    def savePFCfile(self):
        """
        saves a default PFC file (no PFC data, just the template)
        to the tmpDir directory so the user can download
        """
        tools.saveDefaultPFCfile(self.tmpDir)
        return

    def getHFInputs(self,hfMode,
                    lqCN,lqCF,lqPN,lqPF,S,
                    fracCN,fracCF,fracPN,fracPF,
                    fracUI,fracUO,fracLI,fracLO,
                    lqCNmode,lqCFmode,lqPNmode,lqPFmode,SMode,
                    qBG,P,radFrac,fG,
                    qFilePath,qFileTag,
                    rzqFile, rayTracer, rzqFiledata=None, tIdx=0):
        """
        get heat flux inputs from gui or input file
        """
        self.HF.hfMode = hfMode
        self.HF.lqCN = lqCN
        self.HF.lqCF = lqCF
        self.HF.lqPN = lqPN
        self.HF.lqPF = lqPF
        self.HF.S = S
        self.HF.P = P
        self.HF.radFrac = radFrac
        self.HF.qBG = qBG
        self.HF.fracCN = fracCN
        self.HF.fracCF = fracCF
        self.HF.fracPN = fracPN
        self.HF.fracPF = fracPF
        self.HF.fracUI = fracUI
        self.HF.fracUO = fracUO
        self.HF.fracLI = fracLI
        self.HF.fracLO = fracLO
        self.HF.lqCNmode = lqCNmode
        self.HF.lqCFmode = lqCFmode
        self.HF.lqPNmode = lqPNmode
        self.HF.lqPFmode = lqPFmode
        self.HF.SMode = SMode
        self.HF.fG = fG
        self.HF.rayTracer = str(rayTracer).strip()

        allowed_qTags = [None, 'none', 'NA', 'None', 'N']
        if qFileTag in allowed_qTags:
            qFileTag = None
        if qFilePath in allowed_qTags:
            qFilePath = None
        self.HF.qFilePath = qFilePath
        self.HF.qFileTag = qFileTag

        if rzqFiledata != None:
            print("Saving rzq file to tmpDir.")
            self.HF.rzqFile = self.HF.writerzqFileData(rzqFile, rzqFiledata, self.tmpDir)
        else:
            self.HF.rzqFile = rzqFile

        self.HF.readrzqprofile(self.HF.rzqFile)

        self.HF.setTypes()

        print("Calculating Psol:")
        #fraction of power conducted to PFC surfaces
        self.HF.Psol = (1-self.HF.radFrac)*self.HF.P
        print("Psol value is:", self.HF.Psol)


        print("HF Mode = "+hfMode)
        log.info("Hf Mode = "+hfMode)
        if hfMode != 'qFile' or hfMode != 'rzqprofile':
            print("P = {:f}".format(self.HF.P))
            log.info("P = {:f}".format(self.HF.P))
            print("Fraction of P Radiated from Core = {:f}".format(self.HF.radFrac))
            log.info("Fraction of P Radiated from Core = {:f}".format(self.HF.radFrac))
            print("Psol = {:f}".format(self.HF.Psol))
            log.info("Psol = {:f}".format(self.HF.Psol))
            print("Upper Inner Div Power Fraction: {:f}".format(self.HF.fracUI))
            log.info("Upper Inner Div Power Fraction: {:f}".format(self.HF.fracUI))
            print("Upper Outer Div Power Fraction: {:f}".format(self.HF.fracUO))
            log.info("Upper Outer Div Power Fraction: {:f}".format(self.HF.fracUO))
            print("Lower Inner Div Power Fraction: {:f}".format(self.HF.fracLI))
            log.info("Lower Inner Div Power Fraction: {:f}".format(self.HF.fracLI))
            print("Lower Outer Div Power Fraction: {:f}".format(self.HF.fracLO))
            log.info("Lower Outer Div Power Fraction: {:f}".format(self.HF.fracLO))

            #get regression parameters from MHD EQ
            self.HF.getRegressionParams(self.MHD.ep[tIdx])

            if hasattr(self.MHD, 'ep'):
                self.HF.HFdataDict = self.HF.getHFtableData(self.MHD.ep[0])
            else:
                if self.HF.lqCNmode == 'eich':
                    print("You need to load MHD to load Eich profiles")
                    log.info("You need to load MHD to load Eich profiles")
                    self.HF.HFdataDict = {}
                    raise ValueError("MHD EQ must be loaded before HF Settings to calculate Eich lq")
                else:
                    self.HF.HFdataDict = self.HF.getHFtableData(ep=None)

                if self.HF.lqCF == 'horacek':
                    print("You need to load MHD to load Horacek profiles")
                    log.info("You need to load MHD to load Horacek profiles")
                    self.HF.HFdataDict = {}
                    raise ValueError("MHD EQ must be loaded before HF Settings to calculate Horacek lqCF")
                else:
                    self.HF.HFdataDict = self.HF.getHFtableData(ep=None)

                if self.HF.SMode == 'makowski':
                    print("You need to load MHD to load Makowski profiles")
                    log.info("You need to load MHD to load Makowski profiles")
                    self.HF.HFdataDict = {}
                    raise ValueError("MHD EQ must be loaded before HF Settings to calculate Makowski S")
                else:
                    self.HF.HFdataDict = self.HF.getHFtableData(ep=None)

        else: #if hfMode is qFile
            self.HF.HFdataDict = self.HF.getHFtableData(ep=None)

        if hfMode=='eich':
            print("lqCN = {:f}".format(self.HF.lqCN))
            print("S = {:f}".format(self.HF.S))
            print("qBG = {:f}".format(self.HF.qBG))
            print("fG = {:f}".format(self.HF.fG))
            log.info("lqEich = {:f}".format(self.HF.lqCN))
            log.info("S = {:f}".format(self.HF.S))
            log.info("qBG = {:f}".format(self.HF.qBG))
            log.info("fG = {:f}".format(self.HF.fG))
        elif hfMode=='limiter':
            print("lqCN = {:f}".format(self.HF.lqCN))
            print("lqCF = {:f}".format(self.HF.lqCF))
            print("fracCN = {:f}".format(self.HF.fracCN))
            print("fracCF = {:f}".format(self.HF.fracCF))
            log.info("lqCN = {:f}".format(self.HF.lqCN))
            log.info("lqCF = {:f}".format(self.HF.lqCF))
            log.info("fracCN = {:f}".format(self.HF.fracCN))
            log.info("fracCF = {:f}".format(self.HF.fracCF))
        elif hfMode=='multiExp':
            print("lqCN = {:f}".format(self.HF.lqCN))
            print("lqCF = {:f}".format(self.HF.lqCF))
            print("lqPN = {:f}".format(self.HF.lqPN))
            print("lqPF = {:f}".format(self.HF.lqPF))
            print("fracCN = {:f}".format(self.HF.fracCN))
            print("fracCF = {:f}".format(self.HF.fracCF))
            print("fracPN = {:f}".format(self.HF.fracPN))
            print("fracPF = {:f}".format(self.HF.fracPF))
            log.info("lqCN = {:f}".format(self.HF.lqCN))
            log.info("lqCF = {:f}".format(self.HF.lqCF))
            log.info("lqPN = {:f}".format(self.HF.lqPN))
            log.info("lqPF = {:f}".format(self.HF.lqPF))
            log.info("fracCN = {:f}".format(self.HF.fracCN))
            log.info("fracCF = {:f}".format(self.HF.fracCF))
            log.info("fracPN = {:f}".format(self.HF.fracPN))
            log.info("fracPF = {:f}".format(self.HF.fracPF))
        elif hfMode=='qFile':
            print("qFilePath = "+qFilePath)
            print("qFileTag = "+qFileTag)
            log.info("qFilePath = "+qFilePath)
            log.info("qFileTag = "+qFileTag)
        elif hfMode == 'tophat':
            print("lqCN = {:f}".format(self.HF.lqCN))
            log.info("lqCN = {:f}".format(self.HF.lqCN))
        elif hfMode=='rzqprofile': 
            print("rzqFile = ", rzqFile)
            log.info("rzqFile = %s",rzqFile)
        return

    def loadHFParams(self, infile=None, rzqFiledata=None, tIdx=0):
        """
        function for loading HF parameters on the fly (ie in the time loop)
        """
        if infile==None:
            self.initializeHF()
        else:
            self.initializeHF(infile)
        #initialize optical HF data from input file
        self.getHFInputs(self.HF.hfMode,
                         self.HF.lqCN,
                         self.HF.lqCF,
                         self.HF.lqPN,
                         self.HF.lqPF,
                         self.HF.S,
                         self.HF.fracCN,
                         self.HF.fracCF,
                         self.HF.fracPN,
                         self.HF.fracPF,
                         self.HF.fracUI,
                         self.HF.fracUO,
                         self.HF.fracLI,
                         self.HF.fracLO,
                         self.HF.lqCNmode,
                         self.HF.lqCFmode,
                         self.HF.lqPNmode,
                         self.HF.lqPFmode,
                         self.HF.SMode,
                         self.HF.qBG,
                         self.HF.P,
                         self.HF.radFrac,
                         self.HF.fG,
                         self.HF.qFilePath,
                         self.HF.qFileTag,
                         self.HF.rzqFile,
                         self.HF.rayTracer,
                         rzqFiledata,
                         tIdx)
        return


    def loadGYROParams(self, infile=None):
        """
        function for loading GYRO parameters on the fly (ie in the time loop)
        """
        if infile==None:
            self.initializeGYRO()
        else:
            self.initializeGYRO(infile)

        #initialize optical HF data from input file
        self.getGyroInputs(
                         self.GYRO.N_gyroSteps,
                         self.GYRO.N_gyroPhase,
                         self.GYRO.gyroTraceLength,
                         self.GYRO.ionMassAMU,
                         self.GYRO.vMode,
                         self.GYRO.gyroT_eV,
                         self.GYRO.N_vPhase,
                         self.GYRO.N_vSlice,
                         self.GYRO.ionFrac,
                         self.GYRO.gyroSources,
                         )
        return

    def getGyroInputs(self,N_gyroSteps,N_gyroPhase,gyroTraceLength,ionMassAMU,vMode,gyroT_eV,
                      N_vPhase, N_vSlice, ionFrac, gyroSources):
        """
        Sets up the gyro module
        """
        self.GYRO.N_gyroSteps = int(N_gyroSteps)
        self.GYRO.gyroTraceLength = int(gyroTraceLength)
        self.GYRO.gyroT_eV = float(gyroT_eV)
        self.GYRO.N_vSlice = int(N_vSlice)
        self.GYRO.N_vPhase = int(N_vPhase)
        self.GYRO.N_gyroPhase = int(N_gyroPhase)
        self.GYRO.N_MC = self.GYRO.N_gyroPhase*self.GYRO.N_vSlice*self.GYRO.N_vPhase
        self.GYRO.ionMassAMU = float(ionMassAMU)
        self.GYRO.vMode = vMode
        self.GYRO.ionFrac = float(ionFrac)
        #set up power source
        self.GYRO.gyroSourceTag = str(gyroSources)
        if 'allROI' in self.GYRO.gyroSourceTag:
            try:
                self.GYRO.gyroSources = self.CAD.ROIList
            except:
                print("NO CAD Loaded.  Cannot initialize ROIList.")
        else:
            if type(gyroSources) == list: #GUI mode
                self.GYRO.gyroSources = gyroSources
            else: #terminal mode
                self.GYRO.gyroSources = [x.split(":") for x in [self.GYRO.gyroSources]][0]
            #self.CAD.getGyroSources(gyroSources)
            #self.CAD.getGyroSourceMeshes()
            #self.CAD.writeMesh2file(self.CAD.gyroMeshes, self.CAD.gyroSources, path=self.CAD.STLpath)
            print("Gyro Traces Launched from these tiles:")
            log.info("Gyro Traces Launched from these tiles:")
            print(self.GYRO.gyroSources)
            log.info(self.GYRO.gyroSources)

        #toroidal step size taken from MHD object
        self.GYRO.dpinit = self.MHD.dpinit

        #set up GYRO object
        self.GYRO.setupConstants(self.GYRO.ionMassAMU)
        print('Loaded Gyro Orbit Settings')
        print('# Steps per helix period = {:f}'.format(float(N_gyroSteps)))
        print('Gyro tracing distance [degrees] = {:f}'.format(float(gyroTraceLength)))
        print('Plasma Temperature Mode = ' + vMode)
        print('Number of Monte Carlo runs per point = {:f}'.format(float(self.GYRO.N_MC)))
        print("Source of gyro orbit power = "+self.GYRO.gyroSourceTag)
        return
        
    def loadREParams(self, infile=None):
        """
        function for loading RE parameters on the fly (ie in the time loop)
        """
        if infile==None:
            self.initializeRE()
        else:
            self.initializeRE(infile)

        #initialize optical HF data from input file
        self.getREInputs(
                         self.RE.N_gyroSteps,
                         self.RE.gyroTraceLength,
                         self.RE.RE_Eav,
                       	 self.RE.RE_I,
                         self.RE.N_vSlice,
                         self.RE.N_vPhase,
                         self.RE.ionMassAMU,
                         self.RE.vMode,
                         self.RE.gyroSources,
                         )
        return

    def getREInputs(self,N_gyroSteps, gyroTraceLength, RE_Eav, RE_I,
                      N_vPhase, N_vSlice, ionMassAMU, vMode, gyroSources):
        """
        Sets up the RE module
        """
        self.RE.N_gyroSteps = int(N_gyroSteps)
        self.RE.gyroTraceLength = int(gyroTraceLength)
        self.RE.RE_Eav = float(RE_Eav)
        self.RE.RE_I = float(Re_I)
        self.RE.N_vSlice = int(N_vSlice)
        self.RE.N_vPhase = int(N_vPhase)
        self.RE.N_gyroPhase = int(N_gyroPhase)
        self.RE.N_MC = self.GYRO.N_gyroPhase*self.GYRO.N_vSlice*self.GYRO.N_vPhase
        self.RE.ionMassAMU = float(ionMassAMU)
        self.RE.vMode = vMode

        #set up power source
        self.RE.gyroSourceTag = str(gyroSources)
        if 'allROI' in self.RE.gyroSourceTag:
            try:
                self.RE.gyroSources = self.CAD.ROIList
            except:
                print("NO CAD Loaded.  Cannot initialize ROIList.")
        else:
            if type(gyroSources) == list: #GUI mode
                self.RE.gyroSources = gyroSources
            else: #terminal mode
                self.RE.gyroSources = [x.split(":") for x in [self.RE.gyroSources]][0]

            print("Runaway Traces Launched from these tiles:")
            log.info("Runaway Traces Launched from these tiles:")
            print(self.RE.gyroSources)
            log.info(self.RE.gyroSources)

        #toroidal step size taken from MHD object
        self.RE.dpinit = self.MHD.dpinit

        #set up RE object
        self.RE.setupConstants(self.RE.ionMassAMU)
        print('Loaded RE Settings')

        return

    def loadRADParams(self, infile=None):
        """
        function for loading RAD parameters on the fly (ie in the time loop)
        """
        if infile==None:
            self.initializeRAD()
        else:
            self.initializeRAD(infile)

        #initialize optical HF data from input file
        self.getRADInputs(
                          self.RAD.radFile,
                          self.RAD.Ntor,
                          self.RAD.Nref,
                          self.RAD.phiMin,
                          self.RAD.phiMax,
                          self.RAD.rayTracer,
                          self.RAD.Prad_mult,
                          self.RAD.saveRadFrac,
                         )
        return

    def getRADInputs(self, radFile, Ntor, Nref, phiMin, phiMax, 
                     rayTracer, Prad_mult, saveRadFrac, radData=None):
        """
        Sets up the RAD module

        """
        #radData is None when running in terminal mode
        if radData != None:
            self.RAD.radFile = self.RAD.writeRadFileData(radFile, radData, self.tmpDir)
        else:
            self.RAD.radFile = radFile

        #in TUI mode we already set these types, but in GUI mode we do it here
        self.RAD.Ntor = int(Ntor)
        self.RAD.Nref = int(Nref)
        self.RAD.phiMin = float(phiMin)
        self.RAD.phiMax = float(phiMax)
        self.RAD.rayTracer = str(rayTracer)
        self.RAD.Prad_mult = float(Prad_mult)
        self.RAD.saveRadFrac = bool(saveRadFrac)

        #read (R,Z,P) photon radiation source file (csv format)
        self.RAD.read2DSourceFile(self.RAD.radFile)

        #get phi that RZ profiles will be placed at in 3D
        self.RAD.getPhis(Ntor, phiMin, phiMax)

        #extrude 2D radiation profile to 3D profile
        self.RAD.create3DFrom2D()

        print('Number of radiation source points: {:d}'.format(self.RAD.NradPts))
        log.info('Number of radiation source points: {:d}'.format(self.RAD.NradPts))
        print("Radiation module loaded")
        log.info("Radiation module loaded")
        return


    def getFilamentInputs(self):
        """
        sets up the filament module
        """
        return


    def bfieldAtSurface(self, PFC, paraview=False):
        """
        Calculate the B field at tile surface

        """
        ctrs = PFC.centers
        R,Z,phi = tools.xyz2cyl(ctrs[:,0],ctrs[:,1],ctrs[:,2])
        #B vector field
        PFC.Bxyz = self.MHD.Bfield_pointcloud(PFC.ep, R, Z, phi)

        if PFC.mergedPFCs == True:
            print("Scattering back to PFC objects")
            PFC.scatter_back('B', self.IO, self.MHD.shotPath)
        else:
            prefix='BfieldGlyph'
            header = "X,Y,Z,Bx,By,Bz"
            path = PFC.controlfilePath
            label = 'B [T]'
            tag = None
            if self.IO.csvMask == True:
                self.IO.writeGlyphCSV(ctrs,PFC.Bxyz,path,prefix,header,tag)
            if self.IO.vtpPCMask == True:
                self.IO.writeGlyphVTP(ctrs,PFC.Bxyz,label,prefix,path,tag)
            if self.IO.csvMask==False and self.IO.vtpPCMask == False:
                print("To write glyphs (Normal and Bfield) you must choose a PC option")
                log.info("To write glyphs (Normal and Bfield) you must choose a PC option")


        #B scalar point clouds
        #PFC.Bp, PFC.Bt, PFC.Br, PFC.Bz = self.MHD.B_pointclouds(PFC.ep, R, Z)
        #if paraview == True:
        #    self.MHD.write_B_pointclouds(ctrs,PFC.Bp,PFC.Bt,PFC.Br,PFC.Bz,PFC.controlfilePath)

        return

    def bfieldMagnitude(self, PFC):
        """
        Calculate B field magnitude point cloud for Bfield Vectors at tile surface
        """
        PFC.Bmag = np.zeros((len(PFC.Bxyz), 4))
        PFC.Bmag[:,0] = PFC.centers[:,0] # X
        PFC.Bmag[:,1] = PFC.centers[:,1] # Y
        PFC.Bmag[:,2] = PFC.centers[:,2] # Z
        PFC.Bmag[:,3] = np.sqrt(PFC.Bxyz[:,0]**2+PFC.Bxyz[:,1]**2+PFC.Bxyz[:,2]**2)
        PFC.Bsign = np.sign(PFC.ep.g['Bt0'])
        return

    def BtraceMultiple(self, t, data=None, boundbox=True):
        """
        Run a MAFOT structure trace from multiple points defined in the gui/tui

        if boundbox is True, will not check for intersections against (R,Z) PFC contour
        """
        if data is not None:
            data = pd.DataFrame.from_dict(data)[list (data[0].keys())]
            data = data.rename(columns=lambda x: x.strip())
        else:
            try:
                #read the Btrace data from a csv file
                data = self.MHD.readBtraceFile()
            except:
                print("===== Could not read Btrace CSV file!  Skipping Btraces... =====")
                return

        data = data.astype({"x[mm]": float, "y[mm]": float, "z[mm]": float, "traceDirection": int, "Length[deg]":float, "stepSize[deg]":float})

        tIdx = np.where(float(t)==self.MHD.timesteps)[0][0]
        x = data['x[mm]'] / 1000.0
        y = data['y[mm]'] / 1000.0
        z = data['z[mm]'] / 1000.0

        xyz = np.array([x,y,z]).T
        controlfile = '_structCTL.dat'

        if len(xyz.shape) > 1:
            R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        else:
            R,Z,phi = tools.xyz2cyl(xyz[0],xyz[1],xyz[2])


        Bt = self.MHD.ep[tIdx].BtFunc.ev(R,Z)
        BR = self.MHD.ep[tIdx].BRFunc.ev(R,Z)
        BZ = self.MHD.ep[tIdx].BZFunc.ev(R,Z)

        gridfile = self.MHD.shotPath + self.tsFmt.format(t) + '/struct_grid.dat'
        controlfilePath = self.MHD.shotPath + self.tsFmt.format(t) + '/'
        #structOutfile = controlfilePath + 'struct.dat'

        for i in range(len(xyz)):
            print("\n --- Tracing pt {:d} for {:f} degrees ---".format(i, data['Length[deg]'][i]))
            tag = 'pt{:03d}'.format(i)
            self.MHD.ittStruct = data['Length[deg]'][i] / data['stepSize[deg]'][i]
            self.MHD.dpinit = data['stepSize[deg]'][i]
            self.MHD.writeControlFile(controlfile, t, data['traceDirection'][i], mode='struct')
            self.MHD.writeMAFOTpointfile(xyz[i,:],gridfile)
            self.MHD.getFieldpath(1.0, data['stepSize[deg]'][i], gridfile, controlfilePath, controlfile, paraview_mask=True, tag=tag, bbox=boundbox)
            #os.remove(structOutfile)

            outfile = controlfilePath+'struct_'+tag+'.csv'
            self.IO.writeTraceVTP(outfile, 'Field_trace_' + tag, controlfilePath)

            print('Converted file to ParaView formatted CSV.')
            log.info('Converted file to ParaView formatted CSV.')
        return


    def Btrace(self,x,y,z,t,direction,traceDeg,dpinit,tag=None):
        """
        Run a MAFOT structure trace from a point defined in the gui

        assumes xyz is in meters
        """
#        idx = np.where(t==self.MHD.timesteps)[0][0]
#        ep = self.MHD.ep[idx]
        t = int(t)
        direction=int(direction)
        x = float(x)
        y = float(y)
        z = float(z)

        xyz = np.array([x,y,z])
        controlfile = '_structCTL.dat'
        dphi = 1.0

        self.MHD.ittStruct = float(traceDeg)
        gridfile = self.MHD.shotPath + self.tsFmt.format(t) + '/struct_grid.dat'
        controlfilePath = self.MHD.shotPath + self.tsFmt.format(t) + '/'

        self.MHD.ittStruct = traceDeg / dpinit
        self.MHD.dpinit = dpinit
        self.MHD.writeControlFile(controlfile, t, direction, mode='struct')
        self.MHD.writeMAFOTpointfile(xyz,gridfile)
        self.MHD.getFieldpath(dphi, gridfile, controlfilePath, controlfile, paraview_mask=True, tag=tag)
        return

    def gyroTrace(self,x,y,z,t,gPhase,vPhase,gyroTraceLength,dpinit,N_helix,traceDirection,gyroT_eV,tag=None):
        """
        performs a single gyro orbit trace

        (x,y,z) are locations where we launch trace from
        gyroPhase is initial phase angle of orbit in degrees
        gyroTraceLength is the number of degrees we will trace for
        N_gyroSteps is the number of discrete lines we approximate helical path by
        """
        print("\n========Gyro Trace Initialized========")
        #get bField trace from this point
        self.Btrace(x,y,z,t,traceDirection,gyroTraceLength,dpinit,tag)
        #read bField trace csv output

        structOutfile = self.MHD.shotPath + self.tsFmt.format(t) + '/struct.dat'
        controlfilePath = self.MHD.shotPath + self.tsFmt.format(t) + '/'

        BtraceXYZ = tools.readStructOutput(structOutfile) #[m]
        #Setup gyro orbit trace constants and velocities
        self.GYRO.setupConstants()
        v = self.GYRO.temp2thermalVelocity(float(gyroT_eV))
        vPerp = v*np.cos(np.radians(vPhase))
        vParallel = v*np.sin(np.radians(vPhase))
        # Evaluate B
        R,Z,phi = tools.xyz2cyl(float(x),float(y),float(z))#mm => m
        tIdx = np.where(float(t)==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[tIdx]
        Bt = ep.BtFunc.ev(R,Z)
        BR = ep.BRFunc.ev(R,Z)
        BZ = ep.BZFunc.ev(R,Z)
        B = np.sqrt(BR**2 + Bt**2 + BZ**2)
        Bsign = np.sign(ep.g['Bt0'])
        #Calculate frequencies and gyro radius
        self.GYRO.setupFreqs(B*Bsign)
        self.GYRO.setupRadius(vPerp)
        #trace helix and save to CSV and VTK formats
        self.GYRO.singleGyroTrace(vPerp,vParallel,float(gPhase),float(N_helix),BtraceXYZ,controlfilePath,
                                  self.GYRO.TGyro[0],self.GYRO.rGyro[0],self.GYRO.omegaGyro[0],tag=tag)
        print("Perpendicular velocity of {:f} m/s".format(float(vPerp)))
        print("Parallel velocity of {:f} m/s".format(float(vParallel)))
        print("B magnitude = {:f}".format(B))
        return

    def gyroTraceMultiple(self,data,t):
        """
        Run a MAFOT structure trace from multiple points defined in the gui,
        then calculate helical trajectory around trace
        """
        data = pd.DataFrame.from_dict(data)[list (data[0].keys())]
        data = data.rename(columns=lambda x: x.strip())
        data = data.astype({"x[mm]": float, "y[mm]": float, "z[mm]": float,
                            "T[eV]":float,"gPhase[deg]":int,"vPhase[deg]":int,
                            "N_helix":int, "traceDirection": int,
                            "Length[deg]":float, "stepSize[deg]":float})

        t = int(t)
        tIdx = np.where(float(t)==self.MHD.timesteps)[0][0]

        traceDirection=data['traceDirection']
        dpinit = data['stepSize[deg]']
        x = data['x[mm]'] / 1000.0
        y = data['y[mm]'] / 1000.0
        z = data['z[mm]'] / 1000.0
        gPhase = data['gPhase[deg]']
        vPhase = data['vPhase[deg]']
        gyroT_eV = data['T[eV]']
        N_helix = data['N_helix']
        gyroTraceLength = data['Length[deg]']

        xyz = np.array([x,y,z]).T
        controlfile = '_structCTL.dat'
        dphi = 1.0

        if len(xyz.shape) > 1:
            R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
            Ntraces = len(xyz)
        else:
            R,Z,phi = tools.xyz2cyl(xyz[0],xyz[1],xyz[2])
            Ntraces = 1

        controlfilePath = self.MHD.shotPath + self.tsFmt.format(t) + '/'
        structOutfile = controlfilePath + 'struct.dat'

        for i in range(Ntraces):
            self.gyroTrace(x[i],y[i],z[i],t,gPhase[i],vPhase[i],gyroTraceLength[i],
                           dpinit[i],N_helix[i],traceDirection[i],gyroT_eV[i],
                           tag='pt{:03d}'.format(i))
            os.remove(structOutfile)
        return
        
    def REgyroTrace(self,x,y,z,t,gPhase,vPhase,gyroTraceLength,dpinit,N_helix,traceDirection,RE_KE,tag=None):
        """
        todo
        performs a single gyro orbit trace for a runaway electron macro particle

        (x,y,z) are locations where we launch trace from
        gyroPhase is initial phase angle of orbit in degrees
        gyroTraceLength is the number of degrees we will trace for
        N_gyroSteps is the number of discrete lines we approximate helical path by
        """
        print("\n========RE Gyro Trace Initialized========")
        #get bField trace from this point
        self.Btrace(x,y,z,t,traceDirection,gyroTraceLength,dpinit,tag)
        #read bField trace csv output

        structOutfile = self.MHD.shotPath + self.tsFmt.format(t) + '/struct.dat'
        controlfilePath = self.MHD.shotPath + self.tsFmt.format(t) + '/'

        BtraceXYZ = tools.readStructOutput(structOutfile) #[m]
        #Setup gyro orbit trace constants and velocities
        self.RE.setupConstants()
        v = self.RE.E2v_electron(float(RE_KE))
        vPerp = v*np.cos(np.radians(vPhase))
        vParallel = v*np.sin(np.radians(vPhase))
        # Evaluate B
        R,Z,phi = tools.xyz2cyl(float(x),float(y),float(z))#mm => m
        tIdx = np.where(float(t)==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[tIdx]
        Bt = ep.BtFunc.ev(R,Z)
        BR = ep.BRFunc.ev(R,Z)
        BZ = ep.BZFunc.ev(R,Z)
        B = np.sqrt(BR**2 + Bt**2 + BZ**2)
        Bsign = np.sign(ep.g['Bt0'])
        #Calculate frequencies and gyro radius
        self.RE.setupFreqs(B*Bsign, v)
        self.RE.setupRadius(vPerp)
        #trace helix and save to CSV and VTK formats
        self.RE.singleGyroTrace(vPerp,vParallel,float(gPhase),float(N_helix),BtraceXYZ,controlfilePath,
                                  self.GYRO.TGyro[0],self.GYRO.rGyro[0],self.GYRO.omegaGyro[0],tag=tag)
        print("Perpendicular velocity of {:f} m/s".format(float(vPerp)))
        print("Parallel velocity of {:f} m/s".format(float(vParallel)))
        print("B magnitude = {:f}".format(B))
        return

    def REgyroTraceMultiple(self,data,t):
        """
        Run a MAFOT structure trace from multiple points defined in the gui,
        then calculate helical trajectory around trace
        """
        data = pd.DataFrame.from_dict(data)[list (data[0].keys())]
        data = data.rename(columns=lambda x: x.strip())
        data = data.astype({"x[mm]": float, "y[mm]": float, "z[mm]": float,
                            "T[eV]":float,"gPhase[deg]":int,"vPhase[deg]":int,
                            "N_helix":int, "traceDirection": int,
                            "Length[deg]":float, "stepSize[deg]":float})

        t = int(t)
        tIdx = np.where(float(t)==self.MHD.timesteps)[0][0]

        traceDirection=data['traceDirection']
        dpinit = data['stepSize[deg]']
        x = data['x[mm]'] / 1000.0
        y = data['y[mm]'] / 1000.0
        z = data['z[mm]'] / 1000.0
        gPhase = data['gPhase[deg]']
        vPhase = data['vPhase[deg]']
        gyroT_eV = data['T[eV]']
        N_helix = data['N_helix']
        gyroTraceLength = data['Length[deg]']

        xyz = np.array([x,y,z]).T
        controlfile = '_structCTL.dat'
        dphi = 1.0

        if len(xyz.shape) > 1:
            R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
            Ntraces = len(xyz)
        else:
            R,Z,phi = tools.xyz2cyl(xyz[0],xyz[1],xyz[2])
            Ntraces = 1

        controlfilePath = self.MHD.shotPath + self.tsFmt.format(t) + '/'
        structOutfile = controlfilePath + 'struct.dat'

        for i in range(Ntraces):
            self.gyroTrace(x[i],y[i],z[i],t,gPhase[i],vPhase[i],gyroTraceLength[i],
                           dpinit[i],N_helix[i],traceDirection[i],gyroT_eV[i],
                           tag='pt{:03d}'.format(i))
            os.remove(structOutfile)
        return


    def NormPC(self, PFC):
        """
        create a normal vector point cloud for mesh centers on tile surface
        """
        if PFC.mergedPFCs == True:
            print("Scattering back to PFC objects")
            PFC.scatter_back('norm', self.IO, self.MHD.shotPath)
        else:
            prefix='NormGlyph'
            header = "X,Y,Z,Nx,Ny,Nz"
            path = PFC.controlfilePath
            tag = None
            label = 'N'
            if self.IO.csvMask == True:
                self.IO.writeGlyphCSV(PFC.centers,PFC.norms,path,prefix,header,tag)
            if self.IO.vtpPCMask == True:
                self.IO.writeGlyphVTP(PFC.centers,PFC.norms,label,prefix,path,tag)
        return

    def shadowPC(self, PFC):
        """
        create a pointcloud for mesh center locations where 1=shadowed, 0=not shadowed
        """
        prefix = 'shadowMask'
        label = 'shadowMask'
        path = PFC.controlfilePath
        tag = None
        if self.IO.csvMask == True:
            self.IO.writePointCloudCSV(PFC.centers,PFC.shadowed_mask,path,label,tag,prefix)
        if self.IO.vtpPCMask == True:
            self.IO.writePointCloudVTP(PFC.centers,PFC.shadowed_mask,label,prefix,path,tag)
        if self.IO.vtpMeshMask == True:
            self.IO.writeMeshVTP(PFC.mesh, PFC.shadowed_mask, label, prefix, path, tag)
        if self.IO.glbMeshMask == True:
            self.IO.writeMeshGLB(PFC.mesh, PFC.shadowed_mask, label, prefix, path, tag)
        return


    def powerDirPC(self, PFC):
        """
        create a pointcloud for mesh center locations for power direction
        """
        if PFC.mergedPFCs == True:
            print("Scattering back to PFC objects")
            PFC.scatter_back('pwrDir', self.IO, self.MHD.shotPath)
        else:
            prefix = 'powerDir'
            label = 'powerDir'
            path = PFC.controlfilePath
            tag = None
            if self.IO.csvMask == True:
                self.IO.writePointCloudCSV(PFC.centers,PFC.powerDir,path,label,tag,prefix)
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(PFC.centers,PFC.powerDir,label,prefix,path,tag)
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(PFC.mesh, PFC.powerDir, label, prefix, path, tag)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(PFC.mesh, PFC.powerDir, label, prefix, path, tag)
        return

    def bdotnPC(self, PFC):
        """
        makes a cos(alpha) point cloud:b_hat dot n_hat
        where b_hat is normalized magnetic field vector and n_hat is mesh surface
        normal vector
        """
        self.HF.HFincidentAngle(PFC,self.MHD)

        if PFC.mergedPFCs == True:
            print("Scattering back to PFC objects")
            PFC.scatter_back('bdotn', self.IO, self.MHD.shotPath)
        else:
            prefix = 'bdotn'
            label = '$\hat{b} \cdot \hat{n}$'
            path = PFC.controlfilePath
            tag = None
            if self.IO.csvMask == True:
                self.IO.writePointCloudCSV(PFC.centers,PFC.bdotn,path,label,tag,prefix)
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(PFC.centers,PFC.bdotn,label,prefix,path,tag)
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(PFC.mesh, PFC.bdotn, label, prefix, path, tag)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(PFC.mesh, PFC.bdotn, label, prefix, path, tag)
        return

    def initializeHF(self, infile=None):
        """
        Initialize heat flux variables
        """
        print("-"*70)
        print("Heat flux parameters read from file")
        log.info("Heat flux parameters read from file")
        #Initialize HF Object
        if infile == None:
            tools.initializeInput(self.HF, infile=self.infile)
        else:
            tools.initializeInput(self.HF, infile=infile)
        return

    def initializeGYRO(self, infile=None):
        """
        Initialize gyro orbit heat flux variables
        """
        print("-"*70)
        print("Gyro orbit parameters read from file")
        log.info("Gyro orbit parameters read from file")
        #Initialize HF Object
        if infile == None:
            tools.initializeInput(self.GYRO, infile=self.infile)
        else:
            tools.initializeInput(self.GYRO, infile=infile)
        return
        
    def initializeRE(self, infile=None):
        """
        Initialize runaway electron heat flux variables
        """
        print("-"*70)
        print("Runaway Electron parameters read from file")
        log.info("Runaway Electron parameters read from file")
        #Initialize RE Object
        if infile == None:
            tools.initializeInput(self.RE, infile=self.infile)
        else:
            tools.initializeInput(self.RE, infile=infile)
        return

    def initializeRAD(self, infile=None):
        """
        Initialize radiated power heat flux variables
        """
        print("-"*70)
        print("Radiated power parameters read from file")
        log.info("Radiated power parameters read from file")
        #Initialize RAD Object
        if infile == None:
            tools.initializeInput(self.RAD, infile=self.infile)
        else:
            tools.initializeInput(self.RAD, infile=infile)
        return

    def initializeFIL(self, infile=None):
        """
        Initialize filament heat flux variables
        """
        print("-"*70)
        print("Filament power parameters read from file")
        log.info("Filament power parameters read from file")
        #Initialize FIL Object
        if infile == None:
            tools.initializeInput(self.FIL, infile=self.infile)
        else:
            tools.initializeInput(self.FIL, infile=infile)
        return


    def runHEAT(self, runList):
        """
        Run a HEAT calculation.  This is called from gui/tui by user.
        Creates point clouds at each mesh center of PFC objects in CAD ROI

        Steps forward in time, and then solves everything for each PFC at that
        timestep

        runList options are:
        B               Bfield glyphs
        psiN            normalized psi 
        pwrDir          power direction 
        norm            normal glyphs
        bdotn           bdotn
        hfOpt           optical heat flux 
        hfGyro          gyro orbit heat flux 
        hfRad           photon radiation heat flux
        hfFil           filament heat flux
        hfRE            runaway electron trace
        Btrace          magnetic field line trace from file
        """
        print('\n')
        print("-"*70)
        print("HEAT RUN INITIALIZED")
        log.info("HEAT RUN INITIALIZED")
        t0 = time.time()


        #make sure that something in runList can be run in this function, else return
        allowedOptions = ['hfOpt', 'pwrDir', 'bdotn', 'B', 'psiN', 'norm', 'hfGyro', 'hfRad', 'hfFil', 'hfRE', 'Btrace']
        if len([i for i in runList if i in allowedOptions]) < 1:
            self.runList = runList
            print("No HEAT runList option to run.  Breaking out of engineClass runHEAT loop.")
            log.info("No HEAT runList option to run.  Breaking out of engineClass runHEAT loop.")
            return
        else:
            self.runList = runList


        #paraview movie dir
        PVdir = self.MHD.shotPath + "paraview/"
        tools.makeDir(PVdir, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)

        #=========================
        # ===== Steady State =====    
        #=========================

        #set up variables for power balance calculation
        powerTesselate = np.zeros((len(self.MHD.timesteps)))
        powerTrue = np.zeros((len(self.MHD.timesteps)))
        #set up electron frac if not in gyro mode
        if 'hfGyro' not in runList:
            self.HF.elecFrac = 1.0
        else:
            self.HF.elecFrac = 1.0 - self.GYRO.ionFrac

        #list of dictionaries for time varying inputs
        self.inputDicts = []

        # Time Loop 1: HF, bdotn, B, psi, Norms
        for tIdx,t in enumerate(self.MHD.timesteps):
            print('\n')
            print("-"*80)
            log.info("-"*80)
            print("Timestep: "+self.tsFmt.format(t))
            log.info("Timestep: "+self.tsFmt.format(t))
            print("-"*80)
            log.info("-"*80)
            #if user supplied multiple input files in TUI, parse at each timestep
            #note that if running openfoam, only the last timestep's input file
            #will be used for the openFOAM settings.
            if self.inputFileList is not None:
                print("Received a list of input files.")
                log.info("Received a list of input files.")
                self.inputDicts.append(self.loadInputs(inFile=self.inputFileList[tIdx]))
            else:
                self.inputDicts.append(self.getCurrentInputs())
                
            #run B field tracer
            if 'Btrace' in runList:
                self.BtraceMultiple(t, boundbox=False)
                #if this is only a Btrace, skip the PFC dependent steps
                if 'Btrace' in runList and len(runList) == 1:
                    continue

            # 3Dplasma general setup
            if self.plasma3D.plasma3Dmask:
                gFile = self.MHD.shotPath + self.tsFmt.format(t) + '/' + self.MHD.gFiles[tIdx]
                self.plasma3D.initializePlasma3D(self.MHD.shot, t, gFile, self.MHD.tmpDir[0:-1])   # remove / at the end of paths   Also: this no longer reads the input file. This is now done by self.loadInputs
                self.plasma3D.setBoundaryBox(self.MHD, self.CAD)
                self.hf3D.initializeHF3D(self.MHD.tmpDir[0:-1])     # this no longer reads the input file. This is now done by self.loadInputs
                self.plasma3D.print_settings()
                self.hf3D.print_settings()

            #build list of PFCs for this timestep
            PFClist = []
            for PFC in self.PFCs:
                if t not in PFC.timesteps:
                    pass
                else:
                    PFClist.append(PFC)

            #if we are doing all PFCs in one pass
            if self.CAD.mergedPFCs == True:
                if tIdx == 0:
                    #build a merged PFC from all PFCs
                    PFCs = [pfcClass.mergedPFCs(self.PFCs, self.MHD, self.tsSigFigs, self.shotSigFigs, self.chmod, self.UID, self.GID)]
            #if we are doing PFCs one by one
            else:
                PFCs = PFClist

            # Loop through all PFCs
            for PFC in PFCs:
                #set up file directory structure
                PFC.controlfile = '_lamCTL.dat'
                PFC.controlfileStruct = '_struct_CTL.dat'
                PFC.controlfilePath = self.MHD.shotPath + self.tsFmt.format(t) +'/'+PFC.name+'/'
                PFC.gridfile = PFC.controlfilePath + 'grid.dat'
                PFC.gridfileStruct = PFC.controlfilePath + 'struct_grid.dat'
                PFC.outputFile = PFC.controlfilePath + 'lam.dat'
                PFC.structOutfile = PFC.controlfilePath + 'struct.dat'
                #set up time and equilibrium
                PFC.t = t
                PFC.ep = PFC.EPs[tIdx]
                PFC.shadowed_mask = PFC.shadowMasks[tIdx].copy()
                #bfield info for this timestep
                r,z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])
                PFC.BNorms = self.MHD.Bfield_pointcloud(PFC.ep, r, z, phi, powerDir=None, normal=True)
                PFC.bdotn = np.multiply(PFC.norms, PFC.BNorms).sum(1)
                PFC.powerDir = np.sign(PFC.bdotn)*np.sign(PFC.ep.g['Bt0'])*-1.0
                #powerDir can also be calculated using dot product of phi
                #PFC.bdotphi = np.multiply(PFC.BNorms, PFC.phiVec).sum(1)
                #PFC.powerDir = np.sign(PFC.bdotn)*np.sign(PFC.bdotphi)*-1.0
                
                print('\n')
                print("*"*80)
                print('PFC Name: '+ PFC.name+', timestep: '+self.tsFmt.format(t))
                print("*"*80)
                print('\n')
                log.info('\n')
                log.info("*"*80)
                log.info('PFC Name: '+ PFC.name+', timestep: '+self.tsFmt.format(t))
                log.info("*"*80)
                log.info('\n')
                
                # 3Dplasma PFC specific setup
                if self.plasma3D.plasma3Dmask:
                    self.plasma3D.updatePFCdata(PFC.controlfilePath[0:-1])  # remove / at the end of paths
                    self.hf3D.updatePFCdata(PFC.ep, PFC.controlfilePath[0:-1])
                
                if 'hfOpt' in runList:
                    #load HF settings for this timestep if applicable (terminal mode)
                    try:
                        self.loadHFParams(infile=self.inputFileList[tIdx], tIdx=tIdx)
                    except Exception as e:
                        print("Could not load HF parameters.  Expected for GUI.  Check error message:")
                        print(e)
                    
                    #build list of divertor power sharing fracs and broadcast to PFC mesh elements
                    divFracs = [self.HF.fracUI, self.HF.fracUO, self.HF.fracLI, self.HF.fracLO]
                    PFC.broadcastDivFracs(divFracs)

                    #check if this timestep contains an MHD EQ we already traced
                    repeatIdx = self.MHD.check4repeatedEQ(PFC.ep, PFC.EPs[:tIdx])
                    #if the inputs are different, flag it to prevent copying
                    #shadowMask of repeatIdx
                    self.newInputsFlag = True
                    if repeatIdx != None:
                        if self.inputDicts[-1] == self.inputDicts[repeatIdx]:
                            self.newInputsFlag = False
                    #get the optical heat flux
                    if self.HF.qFileTag is None:
                        self.HF_PFC(PFC, repeatIdx, PFC.tag)
                    else:
                        #try to read HF from file
                        val = self.HF.readqFile(PFC, t)
                        #read from file failed, run regular q calculation
                        if val == -1:
                            self.HF_PFC(PFC, repeatIdx, PFC.tag)
                    PFC.shadowMasks[tIdx] = PFC.shadowed_mask.copy()
                    PFC.powerSumOptical[tIdx] = self.HF.power_sum_mesh(PFC, mode='optical')
                    print('\nMaximum optical heat load on tile: {:f}'.format(max(PFC.qDiv)))
                    print('Theoretical optical power to this divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac*self.HF.elecFrac))
                    print('Tessellated divertor power to this PFC = {:f}'.format(PFC.powerSumOptical[tIdx]))
                    log.info('\nMaximum heat load on tile: {:f}'.format(max(PFC.qDiv)))
                    log.info('Theoretical optical power to this divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac*self.HF.elecFrac))
                    log.info('Tessellated Total Power = {:f}'.format(PFC.powerSumOptical[tIdx]))
                    print("Optical Calculation Time Elapsed: {:f}".format(time.time() - t0))
                    log.info("Optical Calculation Time Elapsed: {:f}\n".format(time.time() - t0))
                    powerTesselate[tIdx] += PFC.powerSumOptical[tIdx]
                    
                if 'hfRad' in runList:
                    #load RAD settings for this timestep if applicable (terminal mode)
                    try:
                        self.loadRADParams(infile=self.inputFileList[tIdx])
                    except Exception as e:
                        print("Could not load RAD parameters.  Expected for GUI.  Check error message:")
                        print(e)
                    #location where we will save a memmap if necessary
                    self.RAD.powFracFile = self.MHD.shotPath + self.tsFmt.format(t) +'/photonPowerFrac.nc'
                    #calculate the radiated power on the PFC mesh
                    self.radPower(PFC)
                    #save output files
                    self.radPowerOutput(PFC)
                    PFC.powerSumRad[tIdx] = np.sum(PFC.Prad)
                    print('\nSummation radiated power to this PFC = {:0.10f}'.format(PFC.powerSumRad[tIdx]))
                    log.info('\nSummation radiated power to this PFC = {:0.10f}'.format(PFC.powerSumRad[tIdx]))
                    print('Peak qRad to this PFC: {:0.10f}'.format(np.max(PFC.qRad)))
                    log.info('Peak qRad to this PFC: {:0.10f}'.format(np.max(PFC.qRad)))
                if 'B' in runList:
                    print('Writing Bfield Glyphs')
                    self.bfieldAtSurface(PFC,paraview=True)
                if 'psiN' in runList:
                    self.psiPC(PFC)
                if 'norm' in runList:
                    self.NormPC(PFC)
                if 'pwrDir' in runList:
                    self.powerDirPC(PFC)
                if 'bdotn' in runList:
                    self.bdotnPC(PFC)

            totalPowPow = 0
            totalPowPowCirc = 0
            for PFC in PFCs:
                if 'hfOpt' in runList:
                    print("\n=== Final PFC tallies: Optical ===")
                    tmpPow = self.HF.power_sum_mesh(PFC, scale2circ=False, verbose=False)
                    tmpPowCirc = self.HF.power_sum_mesh(PFC, scale2circ=True, verbose=False)
                    totalPowPow += tmpPow
                    totalPowPowCirc += tmpPowCirc
                    print(PFC.name + ":\t{:.6f}".format(tmpPow))
                    log.info(PFC.name + ":\t{:.6f}".format(tmpPow))
                    print("PFC array sum: {:.6f}".format(totalPowPow))
                    log.info("PFC array sum: {:.6f}".format(totalPowPow))
                    print("scale2circ sum:\t{:.6f}".format(totalPowPowCirc))
                    log.info("scale2circ sum:\t{:.6f}".format(totalPowPowCirc))
                if 'hfRad' in runList:
                    print("=== Final PFC tallies: Photon radiation ===")
                    try:
                        print(PFC.name + ":\t{:.6f}".format(np.sum(PFC.Prad)))
                        log.info(PFC.name + ":\t{:.6f}".format(np.sum(PFC.Prad)))
                    except:
                        print("No radiated power available")
                        log.info("No radiated power available")

            if 'hfGyro' in runList:
                print("\n===+++ GYRO ORBIT CALCULATION +++===")
                log.info("\n===+++ GYRO ORBIT CALCULATION +++===")
                tGyro = time.time()

                #mergedPFCs not supported currently for gyro orbits (could be changed in future)
                if self.CAD.mergedPFCs == True:
                    print("Cannot use mergedPFC mode for gyro orbit calculation. Set mergedPFCs to False in input file.  Aborting")
                    log.info("Cannot use mergedPFC mode for gyro orbit calculation. Set mergedPFCs to False in input file.  Aborting")
                    sys.exit()

                #load GYRO settings for this timestep if applicable (terminal mode)
                try:
                    self.loadGYROParams(infile=self.inputFileList[tIdx])
                except:
                    print('Could not load gyro-orbit parameters.  Expected for GUI.')
                    pass

                self.getGyroMeshes()

                #generate index maps
                self.prepareGyroMaps(tIdx)
                countGyroSources=0
                for PFC in self.PFCs:
                    PFC.powerSumGyro = np.zeros((len(self.MHD.timesteps)))
                    if t not in PFC.timesteps:
                        pass
                    else:
                        if PFC.name in self.GYRO.gyroSources:
                            countGyroSources+=1
                            print("\n------PFC: "+PFC.name+"------")
                            self.gyroOrbitIntersects(PFC)

                if countGyroSources == 0:
                    print("\nYOU DO NOT HAVE ANY GYRO ORBIT SOURCE CAD DEFINED!")
                    print("Cannot complete gyro-orbit calculation.  You must")
                    print("properly define gyro sources in input file and in")
                    print("PFC file.  Exiting...")
                    return

                #redistribute ion optical power and build intersectRecord
                self.gyroOrbitHF()
                #path for this timestep
                tPath = self.MHD.shotPath + self.tsFmt.format(t) + '/'
                #write intersectRecord to CSV file
                self.intersectRecordCSV(tPath)

                #contribute to the tallies
                for PFC in self.PFCs:
                    if t not in PFC.timesteps:
                        pass
                    else:
                        Psum = self.HF.power_sum_mesh(PFC, mode='gyro', scale2circ=True)
                        PFC.powerSumGyro[tIdx] += Psum
                        powerTesselate[tIdx] += Psum
                        print(PFC.name + ' tessellated Gyro Power = {:f}'.format(Psum))
                        print('Theoretical gyro power to this divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac*self.GYRO.ionFrac))
                        print('Gyro Tessellated Power = {:f}'.format(np.sum(PFC.powerSumGyro)))
                        print('Escaped Gyro Power = {:f}'.format(self.GYRO.gyroNanPower))
                        print("Gyro orbit calculation took: {:f} [s]\n".format(time.time() - tGyro))


                #print some stats to terminal after all runs completed
                print("\n=== Final PFC Tallies: Gyro ===")
                log.info("\n=== Final PFC Tallies: Gyro ===")
                totalPowPow = 0
                for PFC in self.PFCs:
                    if t not in PFC.timesteps:
                        pass
                    else:
                        tmpPow = self.HF.power_sum_mesh(PFC, mode='gyro', scale2circ=False, verbose=False)
                        totalPowPow += tmpPow
                        print(PFC.name + ":\t{:.6f}".format(tmpPow))
                        log.info(PFC.name + ":\t{:.6f}".format(tmpPow))
                        print("PFC array sum: {:.6f}".format(totalPowPow))
                        log.info("PFC array sum: {:.6f}".format(totalPowPow))
                        print("PFC Max HF: {:.6f}".format(max(PFC.qGyro)))
                        log.info("PFC Max HF: {:.6f}".format(max(PFC.qGyro)))
                print("=== Final PFC Tallies: Optical ===")
                log.info("=== Final PFC Tallies: Optical ===")
                totalPowPow = 0
                for PFC in self.PFCs:
                    if t not in PFC.timesteps:
                        pass
                    else:
                        tmpPow = self.HF.power_sum_mesh(PFC, scale2circ=False, verbose=False)
                        totalPowPow += tmpPow
                        print(PFC.name + ":\t{:.6f}".format(tmpPow))
                        log.info(PFC.name + ":\t{:.6f}".format(tmpPow))
                        print("PFC array sum: {:.6f}".format(totalPowPow))
                        log.info("PFC array sum: {:.6f}".format(totalPowPow))
                if 'hfRad' in runList:
                    print("=== Final PFC Tallies: Photon Radiation ===")
                    log.info("=== Final PFC Tallies: Photon Radiation ===")
                    totalPowPow = 0
                    for PFC in PFCs:
                        if t not in PFC.timesteps:
                            pass
                        else:
                            tmpPow = self.HF.power_sum_mesh(PFC, mode='rad', scale2circ=False, verbose=False)
                            totalPowPow += tmpPow
                            print(PFC.name + ":\t{:.6f}".format(tmpPow))
                            log.info(PFC.name + ":\t{:.6f}".format(tmpPow))
                            print("PFC array sum: {:.6f}".format(totalPowPow))
                            log.info("PFC array sum: {:.6f}".format(totalPowPow))

            print("Completed all steady state heat flux calculations\n")
            log.info("Completed all steady state heat flux calculations\n")

            #generating allSources heat fluxes
            test1 = 'hfGyro' in runList
            test2 = 'hfOpt' in runList
            test3 = 'hfRad' in runList
            test4 = 'T' in runList #openfoam and elmer read allSources.csv, so we write it
            test5 = 'elmer' in runList
            if (test1 + test2 + test3 + test4 + test5) > 1:
                #set up time and equilibrium
                PFC.t = t
                for PFC in self.PFCs: #using self.PFCs here, not PFCs
                    if t not in PFC.timesteps:
                        pass
                    else:
                        print("Creating allSources CSV files")
                        q = np.zeros((len(PFC.centers)))
                        if 'hfOpt' in runList:
                            q += PFC.qDiv
                        if 'hfGyro' in runList:
                            q += PFC.qGyro
                        if 'hfRad' in runList:
                            q += PFC.qRad
                        #write hf files
                        prefix = 'HF_allSources'
                        label = '$MW/m^2$'
                        path = self.MHD.shotPath + self.tsFmt.format(t) +'/'+PFC.name+'/'
                        if self.IO.csvMask == True:
                            self.IO.writePointCloudCSV(PFC.centers,q,path,label,PFC.tag,prefix)
                        if self.IO.vtpPCMask == True:
                            self.IO.writePointCloudVTP(PFC.centers,q,label,prefix,path,PFC.tag)
                        if self.IO.vtpMeshMask == True:
                            self.IO.writeMeshVTP(PFC.mesh, q, label, prefix, path, PFC.tag)
                        if self.IO.glbMeshMask == True:
                            self.IO.writeMeshGLB(PFC.mesh, q, label, prefix, path, PFC.tag)

        #if we used a mergedPFCs object, revert back to the old PFCs object
        if self.CAD.mergedPFCs == True:
           self.PFCs = PFCs[0].PFClist 

        #if we are only running a Btrace, skip PFC specific operations
        if 'Btrace' in runList and len(runList) == 1:
            pass
        #otherwise postprocess the data into combined pointclouds
        else:

            # Time Loop: postprocessing for steady state heat loads
            for tIdx,t in enumerate(self.MHD.timesteps):
                #path for this timestep
                tPath = self.MHD.shotPath + self.tsFmt.format(t) + '/'
                #merge multiple pointclouds into one single pointcloud for visualization
                self.combinePFCpointcloud(runList, tPath, tIdx)
                #copy each timestep's composite point clouds to central location for
                #paraview postprocessing (movies)
                self.combineTimeSteps(runList, t)




        #=========================
        # ===== Transients =======
        #=========================
        #transient filament heat flux calculation
        if 'hfFil' in runList:
            filDict = self.FIL.filData.to_dict()

            #build filament meshes
            self.getFilMeshes()

            #loop thru ROI PFCs initializing filament HF matrix
            for PFC in self.PFCs:
                #initialize all self.timesteps, even if not in PFC.timesteps
                PFC.qFil = np.zeros((len(PFC.centers), len(self.timesteps)))
                PFC.Edep = np.zeros((len(PFC.centers), len(self.timesteps)))
                PFC.ptclDep = np.zeros((len(PFC.centers), len(self.timesteps)))
                PFC.filTimesteps = self.timesteps



            #loop through each filament
            for idx,ts in enumerate(self.FIL.tsFil):
                print('\n')
                print("-"*80)
                log.info("-"*80)

                id = filDict['id'][idx]
                N_src_t = filDict['N_src_t'][idx]

                print("Filament ID: "+id)
                log.info("Filament ID: "+id)
                print("-"*80)
                log.info("-"*80)

                EtotROI = 0.0
                pTotROI = 0.0
                EtotAll = 0.0
                pTotAll = 0.0

                #get the steady state timestep that precedes this transient timestep
                both = np.intersect1d(ts, self.MHD.timesteps)
                epIdx = np.where( np.min(both)==self.MHD.timesteps )[0][0]
                self.FIL.initializeFilamentFromDict(filDict, idx, self.MHD.ep[epIdx])
                self.FIL.ts = ts


                #loop through source timesteps for this filament
                tCount = 0
                for tIdx,t in enumerate(ts):
                    print('\n')
                    print("-"*30)
                    print("Filament Timestep: "+self.tsFmt.format(t))
                    log.info("\nFilament Timestep: "+self.tsFmt.format(t))

                    #set up file directory structure
                    timeDir = self.MHD.shotPath + self.tsFmt.format(t) + '/'  
                    self.FIL = self.MHD.setupMAFOTdirectory(timeDir, self.FIL)

                    if tIdx == 0:                       
                        #trace magnetic field at filament center at t0
                        Btrace = self.FIL.filamentCtrBtrace(self.MHD, t)

                    #build source for this timestep
                    if tCount < N_src_t:
                        self.getFilamentSource(t, id, Btrace, tIdx)
                    #else:
                    #    break

                    #trace macroparticles from source at this timestep
                    if tIdx < self.FIL.N_src_t:
                        self.FIL.tEQ = ts[0]
                        self.FIL.traceFilamentParticles(self.MHD, ts, tIdx)
                        #loop thru ROI PFCs, mapping power to targets
                        for PFC in self.PFCs:
                            if t not in PFC.timesteps:
                                pass
                            else:
                                print("*"*20)
                                print('PFC Name: '+ PFC.name+', timestep: '+self.tsFmt.format(t))
                                log.info("*"*20)
                                log.info('PFC Name: '+ PFC.name+', timestep: '+self.tsFmt.format(t))

                                pfcDir = self.MHD.shotPath + self.tsFmt.format(t) +'/'+PFC.name+'/'
                                tools.makeDir(pfcDir, clobberFlag=False)
                                self.HF.filamentHeatFlux(self.FIL, PFC, ts, tIdx)
                                self.HF.filamentParticleFlux(self.FIL, PFC, ts, tIdx)
                                pTotROI += np.sum(PFC.ptclDep)
                                EtotROI += np.sum(PFC.Edep)

                        #energy balance calculation
                        energy, particles = self.filDepositedEnergyParticles(tIdx)
                        EtotAll += energy
                        pTotAll += particles

                        print("Generating trace output")
                        log.info("Generating trace output")
                        self.filamentTraceOutput(id,t,tIdx)

                    else:
                        print("No more source timesteps to trace.  Breaking loop.")
                        log.info("No more source timesteps to trace.  Breaking loop.")
                        break

                    tCount += 1


                #print energy balance stats
                print("\n\nTotal Energy Deposited on ROI PFCs: {:f}".format(EtotROI))
                print("Total Energy Deposited on All PFCs: {:f}".format(EtotAll))
                print("Theoretical total energy: {:f}".format(self.FIL.E0))
                print("Energy balance: {:0.3f}%".format(EtotAll / self.FIL.E0 * 100.0))
                log.info("\n\nTotal Energy Deposited on ROI PFCs: {:f}".format(EtotROI))
                log.info("Total Energy Deposited on All PFCs: {:f}".format(EtotAll))
                log.info("Theoretical total energy: {:f}".format(self.FIL.E0))
                log.info("Energy balance: {:0.3f}%".format(EtotAll / self.FIL.E0 * 100.0))

                #print particle balance stats
                Nptcls = self.FIL.N_b*self.FIL.N_r*self.FIL.N_p*self.FIL.N_vS
                print("N Particles Deposited on All PFCs: {:f}".format(np.sum(pTotAll)))
                print("Theoretical N particles: {:f}".format(Nptcls))
                print("Particle balance: {:0.3f}%\n".format(pTotAll / Nptcls * 100.0))
                log.info("N Particles Deposited on All PFCs: {:f}".format(np.sum(pTotAll)))
                log.info("Theoretical N particles: {:f}".format(Nptcls))
                log.info("Particle balance: {:0.3f}%\n".format(pTotAll / Nptcls * 100.0))


                #copy heat fluxes to the paraview movie directory
                print("Copying HF to PV movieDir")
                log.info("Copying HF to PV movieDir")
                self.saveFilamentHFOutput(ts, id)
                self.saveFilamentParticleOutput(ts, id)


            #copy filament sources at this timestep to the paraview movie directory
            print("\nBuilding paraview movie directory...can take some time")
            for i,ts in enumerate(self.FIL.tsFil):
                id = filDict['id'][i]
                N_src_t = filDict['N_src_t'][i]
                tCount = 0
                for t in ts:
                    if tCount < N_src_t: 
                        oldPath = self.MHD.shotPath + self.tsFmt.format(t) + '/paraview/'
                        newPath = self.MHD.shotPath + '/paraview/'
                        name = 'filamentSource_'+id+'_' + self.tsFmt.format(t)
                        self.combineFilTimesteps(name, oldPath, newPath)
                    tCount +=1

        if 'hfRE' in runList:
            REDict = self.RE.REData.to_dict()

            #build filament meshes
            self.getFilMeshes(filtype = 'RE')

            #loop thru ROI PFCs initializing filament HF matrix
            for PFC in self.PFCs:
                #initialize all self.timesteps, even if not in PFC.timesteps
                PFC.qFil = np.zeros((len(PFC.centers), len(self.timesteps)))
                PFC.Edep = np.zeros((len(PFC.centers), len(self.timesteps)))
                PFC.ptclDep = np.zeros((len(PFC.centers), len(self.timesteps)))
                PFC.FilTimesteps = self.timesteps



            #loop through each filament
            for idx,ts in enumerate(self.RE.tsFil):
                print('\n')
                print("-"*80)
                log.info("-"*80)

                id = REDict['id'][idx]
                N_src_t = REDict['N_src_t'][idx]

                print("Runaway ID: "+id)
                log.info("Runaway ID: "+id)
                print("-"*80)
                log.info("-"*80)

                EtotROI = 0.0
                pTotROI = 0.0
                EtotAll = 0.0
                pTotAll = 0.0

                #get the steady state timestep that precedes this transient timestep
                both = np.intersect1d(ts, self.MHD.timesteps)
                print('intersection', both)
                print('ts', ts)
                print('MHD times', self.MHD.timesteps)
                epIdx = np.where( np.min(both)==self.MHD.timesteps )[0][0]


                print(REDict)
                self.RE.initializeREdistFromDict(REDict, idx, self.MHD.ep[epIdx])
                self.RE.ts = ts


                #loop through source timesteps for this filament
                tCount = 0
                for tIdx,t in enumerate(ts):
                    print('\n')
                    print("-"*30)
                    print("Filament Timestep: "+self.tsFmt.format(t))
                    log.info("\nFilament Timestep: "+self.tsFmt.format(t))

                    #set up file directory structure
                    timeDir = self.MHD.shotPath + self.tsFmt.format(t) + '/'  
                    self.RE = self.MHD.setupMAFOTdirectory(timeDir, self.RE)

                    if tIdx == 0:                       
                        #trace magnetic field at filament center at t0
                        Btrace = self.RE.RECtrBtrace(self.MHD, t)

                    #build source for this timestep
                    if tCount < N_src_t:
                        self.RE.createSource(t, Btrace)
                        self.filamentSourceOutput(id, t, tIdx, filtype = 'RE') #akf
                    #else:
                    #    break

                    #trace macroparticles from source at this timestep
                    #If N_src_t = 1 this only runs the one time, this is where all of your particles are created and followed
                    if tIdx < self.RE.N_src_t:
                        self.RE.tEQ = ts[0]
                        self.RE.traceREParticles(self.MHD, ts, tIdx)
                        #loop thru ROI PFCs, mapping power to targets
                        for PFC in self.PFCs:
                            if t not in PFC.timesteps:
                                pass
                            else:
                                print("*"*20)
                                print('PFC Name: '+ PFC.name+', timestep: '+self.tsFmt.format(t))
                                log.info("*"*20)
                                log.info('PFC Name: '+ PFC.name+', timestep: '+self.tsFmt.format(t))

                                pfcDir = self.MHD.shotPath + self.tsFmt.format(t) +'/'+PFC.name+'/'
                                tools.makeDir(pfcDir, clobberFlag=False)
                                self.HF.REHeatFlux(self.RE, PFC, ts, tIdx)
                                #self.HF.REParticleFlux(self.RE, PFC, ts, tIdx)
                                pTotROI += np.sum(PFC.ptclDep)
                                EtotROI += np.sum(PFC.Edep)

                        #energy balance calculation
                        energy, particles = self.filDepositedEnergyParticles(tIdx, filtype = 'RE') #akf
                        EtotAll += energy
                        pTotAll += particles

                        print("Generating trace output")
                        log.info("Generating trace output")
                        self.filamentTraceOutput(id,t,tIdx, filtype = 'RE')

                    else:
                        print("No more source timesteps to trace.  Breaking loop.")
                        log.info("No more source timesteps to trace.  Breaking loop.")
                        break

                    tCount += 1


                #print energy balance stats
                print("\n\nTotal Energy Deposited on ROI PFCs: {:f}".format(EtotROI))
                print("Total Energy Deposited on All PFCs: {:f}".format(EtotAll))
                print("Theoretical total energy: {:f}".format(self.RE.E0))
                print("Energy balance: {:0.3f}%".format(EtotAll / self.RE.E0 * 100.0))
                log.info("\n\nTotal Energy Deposited on ROI PFCs: {:f}".format(EtotROI))
                log.info("Total Energy Deposited on All PFCs: {:f}".format(EtotAll))
                log.info("Theoretical total energy: {:f}".format(self.RE.E0))
                log.info("Energy balance: {:0.3f}%".format(EtotAll / self.RE.E0 * 100.0))

                #print particle balance stats
                Nptcls = self.RE.N_b*self.RE.N_r*self.RE.N_p*self.RE.N_vS
                print("N Particles Deposited on All PFCs: {:f}".format(np.sum(pTotAll)))
                print("Theoretical N particles: {:f}".format(Nptcls))
                print("Particle balance: {:0.3f}%\n".format(pTotAll / Nptcls * 100.0))
                log.info("N Particles Deposited on All PFCs: {:f}".format(np.sum(pTotAll)))
                log.info("Theoretical N particles: {:f}".format(Nptcls))
                log.info("Particle balance: {:0.3f}%\n".format(pTotAll / Nptcls * 100.0))


                #copy heat fluxes to the paraview movie directory
                print("Copying HF to PV movieDir")
                log.info("Copying HF to PV movieDir")
                self.saveFilamentHFOutput(ts, id)
                self.saveFilamentParticleOutput(ts, id)


            #copy filament sources at this timestep to the paraview movie directory
            print("\nBuilding paraview movie directory...can take some time")
            for i,ts in enumerate(self.RE.tsFil):
                id = REDict['id'][i]
                N_src_t = REDict['N_src_t'][i]
                tCount = 0
                for t in ts:
                    if tCount < N_src_t: 
                        oldPath = self.MHD.shotPath + self.tsFmt.format(t) + '/paraview/'
                        newPath = self.MHD.shotPath + '/paraview/'
                        name = 'filamentSource_'+id+'_' + self.tsFmt.format(t)
                        self.combineFilTimesteps(name, oldPath, newPath)
                    tCount +=1

        #set tree permissions
        tools.recursivePermissions(self.MHD.shotPath, self.UID, self.GID, self.chmod)
            
        print("Total Time Elapsed: {:f}".format(time.time() - t0))
        log.info("Total Time Elapsed: {:f}".format(time.time() - t0))
        print("\nCompleted HEAT run\n")
        log.info("\nCompleted HEAT run\n")

        #make a sound when complete
#        os.system('spd-say -t female2 "HEAT run complete"')

        #end of runHEAT
        return


    #==========================================================
    #      runHEAT helper functions
    #==========================================================

    #--- Filaments ---

    def filDepositedEnergyParticles(self, tIdx:int, filtype = 'fil'):
        """
        loops through intersectRecord and calculates the sum of all deposited energy
        and particles on any PFC, including PFCs outside of the ROI
        """
        
        if filtype == 'RE':
            obj = self.RE
        else:
            obj = self.FIL
        density = obj.density[:,:,:,tIdx].reshape(obj.N_b*obj.N_r*obj.N_p)

        energy = np.zeros(density.shape)
        ptcls = 0.0
        for i in range(obj.N_vS):
            hits = np.any(~np.isnan(obj.intersectRecord[i,:,:]), axis=1)
            #ptcls +=  density * self.FIL.velocityFracs[:,i] * hits
            #energy += density * self.FIL.energyFracs[:,i] * hits
            #energy += E * self.FIL.velocityFracs[:,i] * hits
            #ptcls += np.sum(hits)

            ptcls += np.sum(hits)
            energy += obj.E0 * density * obj.energyFracs[:,i] * hits
            
        return np.sum(energy), ptcls

    def getFilamentSource(self, t:float, id:str, Btrace:np.ndarray, tIdx: int):
        """
        generates a filament source object along a flux coordinate grid
        saves output in csv and vtp point cloud formats

        only creates N_src_t sources
        """
        #set up time and equilibrium
        self.FIL.createSource(t, Btrace)
        self.filamentSourceOutput(id, t, tIdx)
        return


    def filamentSourceOutput(self,id: str, t: float, tIdx: int, filtype = 'fil'):
        """
        saves filament source profile
        """
        
        if filtype == 'RE':
            obj = self.RE
        else:
            obj = self.FIL

        #save filament data to file
        tag = self.tsFmt.format(t)
        prefix = 'filamentSource_'+id
        label = 'Filament Source'
        xyzData = obj.xyzPts.reshape(obj.N_b*obj.N_r*obj.N_p, 3)
        scalarData = obj.density[:,:,:,tIdx].reshape(obj.N_b*obj.N_r*obj.N_p)
        path = obj.controlfilePath
        if self.IO.csvMask == True:
            self.IO.writePointCloudCSV(xyzData,scalarData,path,label,tag,prefix)
        if self.IO.vtpPCMask == True:
            self.IO.writePointCloudVTP(xyzData,scalarData,label,prefix,path,tag)

        return

    def filamentTraceOutput(self, id: str, t_source: float, tIdx: int, colorbar=True, filtype = 'fil'):
        """
        saves filament traces in CSV or VTP format

        """ 
        if filtype == 'RE':
            obj = self.RE
        else:
            obj = self.FIL
        
        N_ts = int((obj.tMax - obj.tMin) / obj.dt)+1
        ts = np.linspace(obj.tMin, obj.tMax, N_ts)
        path = self.MHD.shotPath 
        for i in range(obj.N_vS):
            for j,t in enumerate(ts):
                #save filament trajectory data to file
                tag = self.tsFmt.format(t)
                prefix = 'filamentTrace_'+id+'_vS{:03d}_tsSrc'.format(i)+self.tsFmt.format(t_source)
                label = 'Filament Trace'
                xyzData = obj.xyzSteps[i,:,j,:].reshape(obj.N_b*obj.N_r*obj.N_p, 3)
                if colorbar == True:
                    scalarData = obj.density[:,:,:,tIdx].reshape(obj.N_b*obj.N_r*obj.N_p) * obj.energyFracs[:,i]
                else:
                    scalarData = np.ones((obj.N_b*obj.N_r*obj.N_p))
                if self.IO.csvMask == True:
                    self.IO.writePointCloudCSV(xyzData,scalarData,path,label,tag,prefix)
                if self.IO.vtpPCMask == True:
                    self.IO.writePointCloudVTP(xyzData,scalarData,label,prefix,path,tag, PClabel=False)

        return

    def saveFilamentHFOutput(self, ts:np.ndarray, id:str, EdepMask=True):
        """
        saves heat fluxes calculated on PFC in VTP or CSV format
        """
        path = self.MHD.shotPath
        
        #create mesh
        mesh = self.CAD.createEmptyMesh()
        #update mesh placement to reflect global translations
        mesh = self.CAD.globalMeshTranslation(mesh)
        [mesh.addMesh(PFC.mesh) for PFC in self.PFCs]
        c = [PFC.centers for PFC in self.PFCs]
        ctrs = np.concatenate(c)

        #create heat fluxes
        for i,t in enumerate(ts):
            #write hf files
            tag = self.tsFmt.format(t)
            print("Adding PV movieDir timestep: " + tag)
            log.info("Adding PV movieDir timestep: " + tag)

            prefix = 'HF_filaments_all_'+id+'_'
            label = '$W/m^2$'
            prefixE = 'Edep_filaments_all_'+id+'_'
            labelE = '$J$'
            q = np.array([])
            #for PFC in self.PFCs:
            #    q = np.append(q, PFC.qFil[:,i])
            #    mesh.addMesh(PFC.mesh)
            #trying to be faster
            qArr = [PFC.qFil[:,i] for PFC in self.PFCs]
            q = np.concatenate(qArr)
            EArr = [PFC.Edep[:,i] for PFC in self.PFCs]
            Edep = np.concatenate(EArr)

            #save output
            if self.IO.csvMask == True:
                 self.IO.writePointCloudCSV(ctrs,q,path+'paraview/',label,tag,prefix) #fluxes 
                 if EdepMask==True:
                     self.IO.writePointCloudCSV(ctrs,Edep,path+'paraview/',labelE,tag,prefixE) #energies
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(mesh, q, label, prefix, path, tag, PClabel=False)
                if EdepMask==True:
                    self.IO.writeMeshVTP(mesh, Edep, labelE, prefix, path, tag, PClabel=False)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(mesh, q, label, prefix, path, tag, PClabel=False)
                if EdepMask==True:
                    self.IO.writeMeshGLB(mesh, Edep, labelE, prefix, path, tag, PClabel=False)
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(ctrs,q,label,prefix,path+'paraview/',tag, PClabel=True)
                if EdepMask==True:
                    self.IO.writePointCloudVTP(ctrs,Edep,labelE,prefix,path+'paraview/',tag, PClabel=True)

        return


    def saveFilamentParticleOutput(self, ts:np.ndarray, id:str):
        """
        saves particle fluxes calculated on PFC
        """
        path = self.MHD.shotPath
        
        #create mesh
        mesh = self.CAD.createEmptyMesh()
        #update mesh placement to reflect global translations
        mesh = self.CAD.globalMeshTranslation(mesh)
        [mesh.addMesh(PFC.mesh) for PFC in self.PFCs]
        c = [PFC.centers for PFC in self.PFCs]
        ctrs = np.concatenate(c)

        #particle fluxes
        for i,t in enumerate(ts):
            #write hf files
            tag = self.tsFmt.format(t)
            print("Adding PV movieDir timestep: " + tag)
            log.info("Adding PV movieDir timestep: " + tag)

            prefix = 'Particles_filaments_all_'+id+'_'
            label = '$Particles$'
            p = np.array([])
            #for PFC in self.PFCs:
            #    q = np.append(q, PFC.qFil[:,i])
            #    mesh.addMesh(PFC.mesh)
            #trying to be faster
            pArr = [PFC.ptclDep[:,i] for PFC in self.PFCs]
            p = np.concatenate(pArr)

            #particle flux
            if self.IO.csvMask == True:
                 self.IO.writePointCloudCSV(ctrs,p,path+'paraview/',label,tag,prefix) #fluxes 
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(mesh, p, label, prefix, path, tag, PClabel=False)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(mesh, p, label, prefix, path, tag, PClabel=False)
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(ctrs,p,label,prefix,path,tag, PClabel=True)

        return




    def combineFilTimesteps(self, name, oldPath, newPath):
        """
        combines timesteps into a single directory for movie-making in paraview

        uses IO flags to determine which files to copy
        """
        old = oldPath + 'PC_' + name
        new = newPath + name
        tools.makeDir(newPath, clobberFlag=False)

        if self.IO.vtpPCMask == True:
            shutil.copyfile(old+'.vtp', new+'.vtp')
        return


    def getFilMeshes(self, filtype = 'fil'):
        """
        sets up filament meshes, independent of timestep
        
        very similar to getGYROMeshes, so could probably be consolidated into single function
        one day...
        """
        
        if filtype == 'RE':
            print("Using Runaway electrons")
            obj = self.RE
        else:
            obj = self.FIL
        
        print("\nBuilding filament meshes and mappings")
        log.info("\nBuilding filament meshes and mappings")
        totalMeshCounter = 0
        numTargetFaces = 0
        numROIFaces = 0
        targetPoints = []
        targetNorms = []
        obj.CADtargetNames = []
        obj.CADROIindexes = []
        obj.CADROINames = []

      
        #build arrays for intersections
        #first include the PFCs in the ROI
        print("CAD ROI List:")
        print(self.CAD.ROIList)

        for i,target in enumerate(self.CAD.ROImeshes):
            totalMeshCounter+=target.CountFacets
            numTargetFaces += target.CountFacets
            numROIFaces += target.CountFacets
            #append target data
            for face in target.Facets:
                obj.CADtargetNames.append(self.CAD.ROIList[i]) #do this for future HF reassignment
                obj.CADROIindexes.append(i)
                obj.CADROINames.append(self.CAD.ROIList[i])
                targetPoints.append(face.Points)
                targetNorms.append(face.Normal)

        #now include PFCs in the intersection list not in ROI
        for i,target in enumerate(self.CAD.intersectMeshes):
            totalMeshCounter+=target.CountFacets
            #we already have the ROI version of this PFC
            if self.CAD.intersectList[i] in self.CAD.ROIList:
                pass
            else:
                print("Adding target "+self.CAD.intersectList[i]+" to intersects with {:f} faces".format(target.CountFacets))
                numTargetFaces += target.CountFacets
                #append target data
                for face in target.Facets:
                    obj.CADtargetNames.append(self.CAD.intersectList[i]) #do this for future HF reassignment
                    targetPoints.append(face.Points)
                    targetNorms.append(face.Normal)

        #targets
        targetPoints = np.asarray(targetPoints)/1000.0 #scale to m
        targetNorms = np.asarray(targetNorms)
        obj.targetPoints = targetPoints
        obj.targetNorms = targetNorms
        obj.t1 = targetPoints[:,0,:] #target point 1 of mesh triangle
        obj.t2 = targetPoints[:,1,:] #target point 2 of mesh triangle
        obj.t3 = targetPoints[:,2,:] #target point 3 of mesh triangle
        obj.Nt = len(obj.t1)
        obj.intersectCenters = tools.getTargetCenters(targetPoints)
        obj.intersectNorms = np.zeros(targetNorms.shape)
        mag = np.linalg.norm(targetNorms,axis=1)
        for i in range(len(targetNorms)):
            obj.intersectNorms[i,:] = targetNorms[i,:] / mag[i]
        print("Total FIL Intersect Faces: {:d}".format(obj.Nt))

        obj.N_CADROI = len(obj.CADROINames)
        #maps from Targets to ROI
        obj.CADTGT_CADROImap = np.arange(obj.N_CADROI)

        return

    #--- Optical approximation ---
    def HF_PFC(self, PFC, repeatIdx=None, tag=None):
        """
        meat and potatoes of the HF calculation.  Called in loop or by parallel
        processes for each PFC object.  Only calculates optical heat flux
        """
        #Check for intersections with MAFOT struct
        t0 = time.time()
        val = -1
        if self.plasma3D.loadHF:
            f = self.plasma3D.loadBasePath + '/' + self.HF.tsFmt.format(PFC.t) + '/' + PFC.name + '/shadowMask.csv'
            val = plasma3DClass.readShadowFile(f, PFC)
        if val == -1:
            #check if this is a repeated MHD EQ
            #and that the inputs have not changed
            if (repeatIdx == None) or (self.newInputsFlag == True):
                #---shadowMask calculation
                PFC.findOpticalShadows(self.MHD,self.CAD,self.HF.rayTracer,batchMode=False) #walk up field line
                #PFC.findOpticalShadows(self.MHD,self.CAD,self.HF.rayTracer,batchMode=True) #entire trace at once
                #original HEAT homebrew MT ray-triangle method
                #PFC.findShadows_structure(self.MHD, self.CAD)
            else:
                PFC.shadowed_mask = PFC.shadowMasks[repeatIdx].copy()

        #PFC.findIntersectionFreeCADKDTree(self.MHD,self.CAD)
        print("Intersection calculation took {:f} [s]\n".format(time.time() - t0))

        #Run MAFOT laminar for 3D plasmas
        if self.plasma3D.plasma3Dmask:
            print('-'*80)
            print('\n----Solving for 3D plasmas with MAFOT----')
            log.info('\n----Solving for 3D plasmas with MAFOT----')
            use = np.where(PFC.shadowed_mask == 0)[0]
            #self.plasma3D.updatePointsFromCenters(PFC.centers[use])
            self.plasma3D.updatePointsFromVertices(PFC.vertices['x'][use,:], PFC.vertices['y'][use,:], PFC.vertices['z'][use,:], PFC.centers[use])
            if self.plasma3D.loadHF:
                f = self.plasma3D.loadBasePath + '/' + self.HF.tsFmt.format(PFC.t) + '/' + PFC.name
                self.plasma3D.copyAndRead(path = f, tag = 'opticalHF')
            else: 
                self.plasma3D.launchLaminar(self.NCPUs, tag = 'opticalHF')   # use MapDirection = 0. If problem, then we need to split here into fwd and bwd direction separately
                self.plasma3D.cleanUp(tag = 'opticalHF')      # removes the MAFOT log files
            
            # check for invalid points (psimin = 10) and remove; there should be none, but just in case
            invalid = self.plasma3D.checkValidOutput()    # this does NOT change self.plasma3D.psimin
            PFC.shadowed_mask[use[invalid]] = 1
            use = np.where(PFC.shadowed_mask == 0)[0]
            if (len(PFC.centers[use]) != len(self.plasma3D.psimin[~invalid])): 
                raise ValueError('psimin array does not match PFC centers. Abort!')
            
            # define and update psimin and Lc in PFC class
            PFC.psimin = np.zeros(PFC.centers[:,0].shape)
            PFC.Lc = np.zeros(PFC.centers[:,0].shape)
            PFC.psimin[use] = self.plasma3D.psimin[~invalid]
            PFC.Lc[use] = self.plasma3D.Lc[~invalid]
            
            print('\n' + '-'*80)
            print('\n----Calculating 3D Heat Flux Profile----')
            log.info('\n----Calculating 3D Heat Flux Profile----')
            self.hf3D.updateLaminarData(PFC.psimin[use],PFC.Lc[use])

            #TO DO: make this hf3D function compatible with divFracs as a mesh quantity
            PFC.powerFrac = self.HF.getDivertorPowerFraction(PFC.DivCode)
            self.hf3D.heatflux(PFC.DivCode, PFC.powerFrac)                # heat flux is scaled by power fraction here
            print("PFC "+PFC.name+" has {:.2f}% of the total power".format(PFC.powerFrac*100.0))
            log.info("PFC "+PFC.name+" has {:.2f}% of the total power".format(PFC.powerFrac*100.0))

            # assign and smooth the parallel heat flux
            q = np.zeros(PFC.centers[:,0].shape)
            q[use] = self.hf3D.q                                          # this is the parallel heat flux q||
            q = self.hf3D.smoothq(PFC.allNeighbours, PFC.centers, PFC.shadowed_mask, q = q)	# give full grid variables, not 'use' subset
            
            qDiv = np.zeros(PFC.centers[:,0].shape)
            qDiv[use] = q[use] * np.abs(PFC.bdotn[use]) * self.HF.elecFrac        # this is the incident heat flux

        #get psi from gfile for 2D plasmas
        else:
            print('\n----Solving for 2D plasmas with EFIT----')
            log.info('\n----Solving for 2D plasmas with EFIT----')
            self.MHD.psi2DfromEQ(PFC)

            #Create Heat Flux Profile
            print('\n----Calculating Heat Flux Profile----')
            log.info('\n----Calculating Heat Flux Profile----')
            q = self.HF.getHFprofile(PFC)   # this is q||
        
            # get the incident heat flux
            qDiv = self.HF.q_div(PFC, self.MHD, q) * self.HF.elecFrac

        #Save data to class variable for future use
        PFC.q = q
        PFC.qDiv = qDiv
        PFC.qOpticalList.append(PFC.qDiv)
        #Create pointclouds for paraview
        print('\n----Creating Output Files----')
        log.info('\n----Creating Output Files----')
        R,Z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])

        if PFC.mergedPFCs == True:
            print("Scattering back to PFC objects")
            PFC.scatter_back('hfOpt', self.IO, self.MHD.shotPath)
        else:
            #write all the files
            prefix = 'HF_optical'
            label = '$MW/m^2$'
            path = PFC.controlfilePath
            if self.IO.csvMask == True:
                self.IO.writePointCloudCSV(PFC.centers,qDiv,path,label,PFC.tag,prefix)
                self.IO.writePointCloudCSV(PFC.centers,PFC.shadowed_mask,path,'shadowMask',PFC.tag,'shadowMask')
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(PFC.centers,qDiv,label,prefix,path,PFC.tag)
                self.IO.writePointCloudVTP(PFC.centers,PFC.shadowed_mask,'shadowMask','shadowMask',path,PFC.tag)
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(PFC.mesh, qDiv, label, prefix, path, PFC.tag)
                self.IO.writeMeshVTP(PFC.mesh, PFC.shadowed_mask, 'shadowMask','shadowMask', path, PFC.tag)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(PFC.mesh, qDiv, label, prefix, path, PFC.tag)
                self.IO.writeMeshGLB(PFC.mesh, PFC.shadowed_mask, 'shadowMask','shadowMask', path, PFC.tag)


        #structOutfile = MHD.shotPath + self.tsFmt.format(t) +'/struct.csv'
        #HF.PointCloudfromStructOutput(structOutfile)
        return

    #--- Radiated power (photons) ---
    def radPower(self,PFC):
        """
        runs the radiated power calculation
        """
        if self.RAD.rayTracer == 'mitsuba_gpu':
            mode = 'cuda'
        elif self.RAD.rayTracer == 'mitsuba_cuda':
            mode = 'cuda'
        else:
            mode = 'cpu'

        #setup the radiated power calculation
        self.RAD.preparePowerTransfer(PFC, self.CAD)
        #trace rays
        if 'mitsuba' in self.RAD.rayTracer:
            print("Using Mitsuba")
            #calculate photon load on PFC using mitsuba JIT or numpy
            self.RAD.calculatePowerTransferMitsubaJIT(mitsubaMode=mode, fType='ply') 
            #self.RAD.calculatePowerTransferMitsubaNumpy(mitsubaMode='cpu', fType='ply')     
        elif self.RAD.rayTracer=="heat":
            #calculate photon load on PFC using legacy methods (brute force)
            self.RAD.calculatePowerTransfer()
        else:
            print("Using Open3D")
            #calculate photon load on PFC using open3d
            self.RAD.calculatePowerTransferOpen3D()
        #assign variables to the PFC itself
        PFC.Prad = self.RAD.targetPower
        PFC.qRad = PFC.Prad / PFC.areas
        PFC.qRadList.append(PFC.qRad)

        #calculate photon radiation shadowMask
        shadowMask = np.ones((self.RAD.Nj))
        loaded = np.where(PFC.qRad != 0.0)
        shadowMask[loaded] = 0.0
        PFC.radShadowMaskList.append(shadowMask)
        PFC.radShadowMask = shadowMask

        return

    def radPowerOutput(self,PFC):
        """
        saves radiated power output
        """
        #if saveFracs==True:
        #    self.RAD.savePowerFrac(PFC)

        if PFC.mergedPFCs == True:
            print("Scattering back to PFC objects")
            PFC.scatter_back('hfRad', self.IO, self.MHD.shotPath)
        else:
            prefix = 'HF_rad'
            label = '$MW/m^2$'
            path = PFC.controlfilePath
            if self.IO.csvMask == True:
                self.IO.writePointCloudCSV(PFC.centers,PFC.qRad,path,label,PFC.tag,prefix)
                self.IO.writePointCloudCSV(self.RAD.sources,self.RAD.sourcePower,path,'$MW$',PFC.tag,'Prad')
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(PFC.centers,PFC.qRad,label,prefix,path,PFC.tag)
                self.IO.writePointCloudVTP(self.RAD.sources,self.RAD.sourcePower,'$MW$','Prad',path,PFC.tag)
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(PFC.mesh, PFC.qRad, label, prefix, path, PFC.tag)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(PFC.mesh, PFC.qRad, label, prefix, path, PFC.tag)

        return



    #--- Gyro Orbits ---

    def gyroOrbitIntersects(self, PFC, mode='open3d'):
        """
        Calculates the gyro orbit intersects for a PFC object

        overwrites shadowMask and psimin
        """
        t0 = time.time()
        #Get B field vector on PFC surface
        self.bfieldAtSurface(PFC)
        #Get B field magnitude on surface
        self.bfieldMagnitude(PFC)
        #Trace B field upstream from PFC surface
        PFC.findGuidingCenterPaths(self.MHD, self.GYRO)
        #Trace helical path downstream, checking for intersections
        if mode == 'open3d':
            #use Open3D ray tracing (100X faster)
            PFC.findHelicalPathsOpen3D(self.GYRO)
        else:
            #use HEAT homebrew ray tracing
            PFC.findHelicalPaths(self.GYRO)

        return

    def gyroOrbitHF(self):
        """
        loop thru PFCs reassigning power based upon intersectRecord
        """
        #find gyro shadowMask
        self.GYRO.shadowMask = np.ones((self.GYRO.N_CADROI))
        shadows = np.unique(self.GYRO.intersectRecord.flatten())
        shadows = shadows[~np.isnan(shadows)]
        shadows = shadows.astype(int)
        shadows = np.where(shadows <= self.GYRO.N_CADROI)[0]
        self.GYRO.shadowMask[shadows] = 0.0

        #redistribute power
        for PFC in self.PFCs:
            if PFC.name in self.GYRO.gyroSources:
                #setup velocities and velocity phase angles
                self.GYRO.setupVelocities(PFC.N_gyroCenters)
                #setup gyroPhase angle
                self.GYRO.uniformGyroPhaseAngle()
                #setup frequencies
                self.GYRO.setupFreqs(PFC.Bmag[PFC.PFC_GYROmap,-1]*PFC.Bsign)
                #self.HF.gyroHF(self.GYRO, PFC)
                self.HF.gyroHF2(self.GYRO, PFC)

        for PFC in self.PFCs:
            #print("PFC NAME: "+PFC.name)
            #print(PFC.CADTGT_PFCmap)
            #assign gyro power to PFC
            PFC.Pgyro = self.GYRO.gyroPowMatrix[PFC.CADTGT_PFCmap]
            PFC.qGyro = PFC.Pgyro / PFC.areas
            PFC.qGyroList.append(PFC.qGyro)
            #gyro shadowMask = 1 if a face is shadowed
            PFC.gyroShadowMask = self.GYRO.shadowMask[PFC.CADROI_PFCmap]
            PFC.gyroShadowMaskList.append(PFC.gyroShadowMask)
            #Create pointcloud for paraview
            R,Z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])

            prefix = 'HF_gyro'
            label = '$MW/m^2$'
            path = PFC.controlfilePath

            if 'hfRad' in self.runList:
                qAll = PFC.qDiv + PFC.qGyro + PFC.qRad
            else:
                qAll = PFC.qDiv + PFC.qGyro

            if self.IO.csvMask == True:
                self.IO.writePointCloudCSV(PFC.centers,PFC.qGyro,path,label,PFC.tag,prefix)
                self.IO.writePointCloudCSV(PFC.centers,qAll,path,label,PFC.tag,'HF_allSources')
                self.IO.writePointCloudCSV(PFC.centers,PFC.qGyro+PFC.qDiv,path,label,PFC.tag,prefix)
                self.IO.writePointCloudCSV(PFC.centers,PFC.gyroShadowMask,path,'shadowMask',PFC.tag,'shadowMaskGyro')
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(PFC.centers,PFC.qGyro,label,prefix,path,PFC.tag)
                self.IO.writePointCloudVTP(PFC.centers,qAll,label,'HF_allSources',path,PFC.tag)
                self.IO.writePointCloudVTP(PFC.centers,PFC.gyroShadowMask,'shadowMask','shadowMaskGyro',path,PFC.tag)
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(PFC.mesh, PFC.qGyro, label, prefix, path, PFC.tag)
                self.IO.writeMeshVTP(PFC.mesh, qAll, label, 'HF_allSources', path, PFC.tag)
                self.IO.writeMeshVTP(PFC.mesh, PFC.gyroShadowMask, 'shadowMask','shadowMaskGyro', path, PFC.tag)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(PFC.mesh, PFC.qGyro, label, prefix, path, PFC.tag)
                self.IO.writeMeshGLB(PFC.mesh, qAll, label, 'HF_allSources', path, PFC.tag)
                self.IO.writeMeshGLB(PFC.mesh, PFC.gyroShadowMask, 'shadowMask','shadowMaskGyro', path, PFC.tag)


        return

    def getGyroMeshes(self):
        """
        set up gyro meshes, independent of timestep
        """
        print("\nBuilding gyro-orbit meshes and mappings")
        log.info("\nBuilding gyro-orbit meshes and mappings")
        totalMeshCounter = 0
        numTargetFaces = 0
        numROIFaces = 0
        numGyroFaces = 0
        targetPoints = []
        targetNorms = []
        gyroPoints = []
        gyroNorms = []
        self.GYRO.CADtargetNames = []
        self.GYRO.CADROINames = []
        self.GYRO.CADROIindexes = []
        self.GYRO.CADROICenters = []


        #build arrays for intersections
        #first include the PFCs in the ROI
        print("CAD ROI List:")
        print(self.CAD.ROIList)
        print("GYRO Source List")
        print(self.GYRO.gyroSources)

        for i,target in enumerate(self.CAD.ROImeshes):
            totalMeshCounter+=target.CountFacets
            numTargetFaces += target.CountFacets
            numROIFaces += target.CountFacets
            #append target data
            for face in target.Facets:
                self.GYRO.CADtargetNames.append(self.CAD.ROIList[i]) #do this for future HF reassignment
                self.GYRO.CADROIindexes.append(i)
                self.GYRO.CADROINames.append(self.CAD.ROIList[i])
                targetPoints.append(face.Points)
                targetNorms.append(face.Normal)

        #now include PFCs in the intersection list not in ROI
        for i,target in enumerate(self.CAD.intersectMeshes):
            totalMeshCounter+=target.CountFacets
            #we already have the ROI version of this PFC
            if self.CAD.intersectList[i] in self.CAD.ROIList:
                pass
            else:
                print("Adding target "+self.CAD.intersectList[i]+" to intersects with {:f} faces".format(target.CountFacets))
                numTargetFaces += target.CountFacets
                #append target data
                for face in target.Facets:
                    self.GYRO.CADtargetNames.append(self.CAD.intersectList[i]) #do this for future HF reassignment
                    targetPoints.append(face.Points)
                    targetNorms.append(face.Normal)

        #PFCs that are gyroSources
        #for i,target in enumerate(self.CAD.gyroMeshes):
        #    numGyroFaces += target.CountFacets
        #    #append target data
        #    for face in target.Facets:
        #        gyroPoints.append(face.Points)
        #        gyroNorms.append(face.Normal)
        #gyro sources
        #gyroPoints = np.asarray(gyroPoints)/1000.0 #scale to m
        #gyroNorms = np.asarray(gyroNorms)
        #self.GYRO.s1 = gyroPoints[:,0,:] #source point 1 of mesh triangle
        #self.GYRO.s2 = gyroPoints[:,1,:] #source point 2 of mesh triangle
        #self.GYRO.s3 = gyroPoints[:,2,:] #source point 3 of mesh triangle
        #self.GYRO.Ns = len(self.GYRO.s1)
        #self.GYRO.gyroCenters = tools.getTargetCenters(gyroPoints)

        #targets
        targetPoints = np.asarray(targetPoints)/1000.0 #scale to m
        targetNorms = np.asarray(targetNorms)
        self.GYRO.t1 = targetPoints[:,0,:] #target point 1 of mesh triangle
        self.GYRO.t2 = targetPoints[:,1,:] #target point 2 of mesh triangle
        self.GYRO.t3 = targetPoints[:,2,:] #target point 3 of mesh triangle
        self.GYRO.Nt = len(self.GYRO.t1)
        self.GYRO.intersectCenters = tools.getTargetCenters(targetPoints)
        self.GYRO.intersectNorms = np.zeros(targetNorms.shape)
        mag = np.linalg.norm(targetNorms,axis=1)
        for i in range(len(targetNorms)):
            self.GYRO.intersectNorms[i,:] = targetNorms[i,:] / mag[i]
        print("Total Gyro Intersect Faces: {:d}".format(self.GYRO.Nt))

        self.GYRO.N_CADROI = len(self.GYRO.CADROINames)
        #maps from Targets to ROI
        self.GYRO.CADTGT_CADROImap = np.arange(self.GYRO.N_CADROI)

        return

    def prepareGyroMaps(self, tIdx):
        """
        timestep dependent gyro calculation preparations

        there are mappings between various abstract containers, that are defined below.
        The mapping between containers is assigned variables with a <FROM|TO>map
        format.  Each container has centers, names, etc. associated with it.  Nested
        maps are represented in comments elsewhere in the code with (flawed) Dirac notation:
        ie: <ROI|PFC><PFC|GYRO> = <ROI|GYRO> = ROIPFCmap[PFCGYROmap]

        CADTGT:  all target faces, at various resolutions, from CAD
        CADROI:  all ROI faces, at single HF resolution, in order of CAD
        PFCROI:  all ROI faces, at single HF resolution, in order of PFCs
        HOT:  CADROI faces that are not optically shadowed ie hot
        PFC:  PFCROI faces on a specific PFC
        GYRO: PFC faces that are not shadowed
        HLX:  Faces that we are calculating helical trajectories on

        In the crappy picture below, non-parentheses variables are containers
        and parentheses variables are maps.  The lines signify the nesting

                          CADTGT
                             |
                             |
                    (CADTGT_CADROImap)
                             |
                             |
            ---------------CADROI-----(ROI_HOTmap)-------
            |                |                          |
            |                |                          |
            |         (CADROI_CADPFCmap)                |
            |                |                          |
     (CADROI_PFCmap)         |                          |
            |              PFCROI                       |
            |                |                          |
            |                |                          |
            |         (PFCROI_PFCmap)                  HOT-----IntersectRecord
            |                |                          |
            |                |                          |
            ----------------PFC                         |
                             |                          |
                             |                          |
                       (PFC_GYROmap)                    |
                             |                          |
                             |                          |
                            GYRO------(HOT_GYROmap)------
                             |
                             |
                       (GYRO_HLXmap)
                             |
                             |
                            HLX

        """

        Npoints = 0
        self.GYRO.N_HOT = 0
        #make arrays for the faces we care about (not optically shadowed)
        for i,PFC in enumerate(self.PFCs):
            if self.MHD.timesteps[tIdx] not in PFC.timesteps:
                pass
            else:
                #maps from CADROI to this PFC: <CADROI|PFC>
                #PFC.CADROI_PFCmap = np.where(np.array(self.GYRO.CADROINames)==PFC.name)[0]
                test = np.logical_and(np.array(self.GYRO.CADROINames)==PFC.name, np.array(self.GYRO.CADROIindexes)==i)
                PFC.CADROI_PFCmap = np.where(test==True)[0]
                #maps from PFC to GYRO (not optically shadowed) faces: <PFC|GYRO>
                PFC.PFC_GYROmap = np.where(PFC.shadowMasks[tIdx] == 0)[0]
                #maps from the CADROI to this PFCs gyro: <CADROI|PFC><PFC|GYRO>
                PFC.CADROI_GYROmap = PFC.CADROI_PFCmap[PFC.PFC_GYROmap]
                PFC.gyroCenters = PFC.centers[PFC.PFC_GYROmap]
                PFC.N_gyroCenters = len(PFC.gyroCenters)

                if i==0:
                    #maps from CADROI to HOT (all PFCs not shadowed) faces: <CADROI|HOT>
                    self.GYRO.CADROI_HOTmap = PFC.CADROI_GYROmap
                    self.GYRO.PFCROI_HOTmap = PFC.PFC_GYROmap
                    self.GYRO.PFCROINames = [PFC.name]*len(PFC.centers)
                    self.GYRO.PFCROIindexes = np.ones((len(PFC.centers)))*i
                    #maps from CAD indexes (not always the same) to PFC indexes (always from PFC file)
                    #<CADROI|PFCROI>
                    self.GYRO.CADROI_PFCROImap = PFC.CADROI_PFCmap
                else:
                    self.GYRO.CADROI_HOTmap = np.append(self.GYRO.CADROI_HOTmap, PFC.CADROI_GYROmap)
                    self.GYRO.PFCROI_HOTmap = np.append(self.GYRO.PFCROI_HOTmap, PFC.PFC_GYROmap+Npoints)
                    self.GYRO.PFCROINames += [PFC.name]*len(PFC.centers)
                    self.GYRO.PFCROIindexes = np.append(self.GYRO.PFCROIindexes, np.ones((len(PFC.centers)))*i)
                    self.GYRO.CADROI_PFCROImap = np.append(self.GYRO.CADROI_PFCROImap, PFC.CADROI_PFCmap)

                self.GYRO.N_HOT += int(PFC.N_gyroCenters)
                Npoints += len(PFC.centers)
                #maps from HOT to GYRO: <CAD|HOT|GYRO> or <PFC|HOT|GYRO>
                #PFC.CADHOT_GYROmap = np.where(np.array(self.GYRO.CADROINames)[self.GYRO.CADROI_HOTmap]==PFC.name)[0]
                #PFC.PFCHOT_GYROmap = np.where(np.array(self.GYRO.PFCROINames)[self.GYRO.PFCROI_HOTmap]==PFC.name)[0]
                test = np.logical_and(np.array(self.GYRO.CADROINames)[self.GYRO.CADROI_HOTmap]==PFC.name, np.array(self.GYRO.CADROIindexes)[self.GYRO.CADROI_HOTmap]==i)
                PFC.CADHOT_GYROmap = np.where(test==True)[0]
                test = np.logical_and(np.array(self.GYRO.PFCROINames)[self.GYRO.PFCROI_HOTmap]==PFC.name, np.array(self.GYRO.PFCROIindexes)[self.GYRO.CADROI_HOTmap]==i)
                PFC.PFCHOT_GYROmap = np.where(test==True)[0]

        #create mapping between ALLPFCs and ROI (because CAD list order can be diff
        #from PFC list order)
        for i,PFC in enumerate(self.PFCs):
            if self.MHD.timesteps[tIdx] not in PFC.timesteps:
                pass
            else:
                #maps from PFCROI to this PFC: <PFCROI|PFC>
                #PFC.PFCROI_PFCmap = np.where(np.array(self.GYRO.PFCROINames)==PFC.name)[0]
                test = np.logical_and(np.array(self.GYRO.PFCROINames)==PFC.name, np.array(self.GYRO.PFCROIindexes)==i)
                PFC.PFCROI_PFCmap = np.where(test==True)[0]
                #maps from targets to this PFC: <CADTGT|CADROI><CADROI|PFCROI><PFCROI|PFC>
                PFC.CADTGT_PFCmap = self.GYRO.CADTGT_CADROImap[PFC.CADROI_PFCmap]


        print("Total Gyro ROI Faces: {:d}".format(self.GYRO.N_CADROI))
        print("Gyro faces not optically shadowed: {:d}".format(self.GYRO.N_HOT))
        self.GYRO.gyroPowMatrix = np.zeros((self.GYRO.Nt))
        self.GYRO.gyroHFMatrix = np.zeros((self.GYRO.Nt))
        self.GYRO.gyroNanPower = 0.0
        #set up intersectRecord, which records index of intersection face
        #this 4D array has:
        #one dimension for gyroPhase angles,
        #one dimension for vPhases,
        #one dimension for vSlices,
        #and one dimesion for the PFC.centers points we are tracing
        #each element is the face that was intersected for that macroparticle
        self.GYRO.intersectRecord = np.ones((self.GYRO.N_gyroPhase,
                                            self.GYRO.N_vPhase,
                                            self.GYRO.N_vSlice,
                                            self.GYRO.N_HOT), dtype=int)*np.nan

        self.GYRO.hdotn = np.ones((self.GYRO.N_gyroPhase,
                                            self.GYRO.N_vPhase,
                                            self.GYRO.N_vSlice,
                                            self.GYRO.N_HOT), dtype=int)*np.nan
        return


    def intersectRecordCSV(self, tPath):
        """
        writes intersectRecord to CSV file
        """
        for gyroPhase in range(self.GYRO.N_gyroPhase):
            for vPhase in range(self.GYRO.N_vPhase):
                for vSlice in range(self.GYRO.N_vSlice):
                    file = tPath+'intersectRecord_All_{:d}_{:d}_{:d}.dat'.format(gyroPhase,vPhase,vSlice)
                    #self.GYRO.writeIntersectRecord(gyroPhase,vPhase,vSlice,self.GYRO.CADROI_HOTmap,file)
                    self.GYRO.writeIntersectRecord(gyroPhase,vPhase,vSlice,self.GYRO.PFCROI_HOTmap,file)

        return


    #--- Generic Outputs + File Saving ---
    def combinePFCpointcloud(self, runList, tPath, tIdx):
        """
        Combines multiple pointclouds into a single pointcloud then saves to file.
        Optimized to avoid quadratic growth from np.append and to cut copies.
        """
        vtpMeshMask = self.IO.vtpMeshMask
        vtpPCMask   = self.IO.vtpPCMask
        csvMask     = self.IO.csvMask
        glbMeshMask = self.IO.glbMeshMask

        want_opt  = 'hfOpt'  in runList
        want_gyro = 'hfGyro' in runList
        want_rad  = 'hfRad'  in runList
        want_all  = want_opt or want_gyro or want_rad
        want_shPC = 'shadowPC' in runList
        want_pdir = 'pwrDir' in runList
        want_bdn  = 'bdotn'  in runList
        want_psi  = 'psiN'   in runList
        want_norm = 'norm'   in runList
        want_B    = 'B'      in runList

        # -------- Pass 1: count points and merge mesh --------
        mesh = self.CAD.createEmptyMesh()
        mesh = self.CAD.globalMeshTranslation(mesh)

        counts = []
        valid_pfcs = []
        for PFC in self.PFCs:
            if self.MHD.timesteps[tIdx] not in PFC.timesteps:
                print("This timestep not in PFC.timesteps. Skipping PFC…")
                continue
            mesh.addMesh(PFC.mesh)
            n = len(PFC.centers)
            counts.append(n)
            valid_pfcs.append(PFC)

        if not counts:
            print("No PFCs at this timestep; nothing to combine.")
            return

        Npoints = int(sum(counts))

        # -------- Preallocate outputs we actually need --------
        centers = np.empty((Npoints, 3), dtype=np.float32)

        hfOpticalNumpy   = np.empty(Npoints, dtype=np.float32) if want_opt  else None
        hfGyroNumpy      = np.empty(Npoints, dtype=np.float32) if want_gyro else None
        hfRadNumpy       = np.empty(Npoints, dtype=np.float32) if want_rad  else None
        hfAllNumpy       = np.empty(Npoints, dtype=np.float32) if want_all  else None

        shadowNumpy      = np.empty(Npoints, dtype=np.uint8)   if want_opt or want_shPC else None
        shadowGyroNumpy  = np.empty(Npoints, dtype=np.uint8)   if want_gyro else None
        shadowRadNumpy   = np.empty(Npoints, dtype=np.uint8)   if want_rad  else None

        powerDirNumpy    = np.empty(Npoints, dtype=np.int8)    if want_pdir else None
        bdotnNumpy       = np.empty(Npoints, dtype=np.float32) if want_bdn  else None
        psiNumpy         = np.empty(Npoints, dtype=np.float32) if want_psi  else None
        normNumpy        = np.empty((Npoints, 3), dtype=np.float32) if want_norm else None
        bFieldNumpy      = np.empty((Npoints, 3), dtype=np.float32) if want_B    else None

        # -------- Pass 2: fill by slices --------
        off = 0
        for PFC, n in zip(valid_pfcs, counts):
            s, e = off, off + n

            # centers
            # (cast once to float32 to avoid hidden copies later)
            centers[s:e, :] = PFC.centers.astype(np.float32, copy=False)

            # optical
            if want_opt:
                qopt = PFC.qOpticalList[tIdx].astype(np.float32, copy=False)
                hfOpticalNumpy[s:e] = qopt
                if shadowNumpy is not None:
                    shadowNumpy[s:e] = PFC.shadowMasks[tIdx].astype(np.uint8, copy=False)

            # gyro
            if want_gyro:
                qg = PFC.qGyroList[tIdx].astype(np.float32, copy=False)
                hfGyroNumpy[s:e] = qg
                if shadowGyroNumpy is not None:
                    shadowGyroNumpy[s:e] = PFC.gyroShadowMaskList[tIdx].astype(np.uint8, copy=False)

            # rad
            if want_rad:
                qr = PFC.qRadList[tIdx].astype(np.float32, copy=False)
                hfRadNumpy[s:e] = qr
                if shadowRadNumpy is not None:
                    shadowRadNumpy[s:e] = PFC.radShadowMaskList[tIdx].astype(np.uint8, copy=False)

            # all sources (sum as we go)
            if want_all:
                total = np.zeros(n, dtype=np.float32)
                if want_opt:  total += hfOpticalNumpy[s:e]
                if want_gyro: total += hfGyroNumpy[s:e]
                if want_rad:  total += hfRadNumpy[s:e]
                hfAllNumpy[s:e] = total

            if want_pdir:
                powerDirNumpy[s:e] = PFC.powerDir.astype(np.int8, copy=False)

            if want_bdn:
                bdotnNumpy[s:e] = PFC.bdotn.astype(np.float32, copy=False)

            if want_psi:
                psiNumpy[s:e] = PFC.psimin.astype(np.float32, copy=False)

            if want_norm:
                # Nx,Ny,Nz per point
                normNumpy[s:e, :] = PFC.norms.astype(np.float32, copy=False)

            if want_B:
                bFieldNumpy[s:e, :] = PFC.Bxyz.astype(np.float32, copy=False)

            off = e

        tag = 'all'

        # -------- Write files --------
        if want_opt:
            prefix = 'HF_optical'; label = '$MW/m^2$'
            if csvMask:   self.IO.writePointCloudCSV(centers, hfOpticalNumpy, tPath, label, tag, prefix)
            if vtpPCMask: self.IO.writePointCloudVTP(centers, hfOpticalNumpy, label, prefix, tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, hfOpticalNumpy, label, prefix, tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, hfOpticalNumpy, label, prefix, tPath, tag)
            # shadow paired with optical
            if csvMask:   self.IO.writePointCloudCSV(centers, shadowNumpy, tPath, 'shadowMask', tag, 'shadowMask')
            if vtpPCMask: self.IO.writePointCloudVTP(centers, shadowNumpy, 'shadowMask', 'shadowMask', tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, shadowNumpy, 'shadowMask', 'shadowMask', tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, shadowNumpy, 'shadowMask', 'shadowMask', tPath, tag)

        if want_gyro:
            prefix = 'HF_gyro'; label = '$MW/m^2$'
            if csvMask:   self.IO.writePointCloudCSV(centers, hfGyroNumpy, tPath, label, tag, prefix)
            if vtpPCMask: self.IO.writePointCloudVTP(centers, hfGyroNumpy, label, prefix, tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, hfGyroNumpy, label, prefix, tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, hfGyroNumpy, label, prefix, tPath, tag)
            # mask
            if csvMask:   self.IO.writePointCloudCSV(centers, shadowGyroNumpy, tPath, 'shadowMask', tag, 'shadowMaskGyro')
            if vtpPCMask: self.IO.writePointCloudVTP(centers, shadowGyroNumpy, 'shadowMask', 'shadowMaskGyro', tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, shadowGyroNumpy, 'shadowMask', 'shadowMaskGyro', tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, shadowGyroNumpy, 'shadowMask', 'shadowMaskGyro', tPath, tag)

        if want_rad:
            prefix = 'HF_rad'; label = '$MW/m^2$'
            if csvMask:
                self.IO.writePointCloudCSV(centers, hfRadNumpy, tPath, label, tag, prefix)
                self.IO.writePointCloudCSV(self.RAD.sources, self.RAD.sourcePower, tPath, '$MW$', tag, 'Prad')
            if vtpPCMask:
                self.IO.writePointCloudVTP(centers, hfRadNumpy, label, prefix, tPath, tag)
                self.IO.writePointCloudVTP(self.RAD.sources, self.RAD.sourcePower, '$MW$', 'Prad', tPath, tag)
            if vtpMeshMask: self.IO.writeMeshVTP(mesh, hfRadNumpy, label, prefix, tPath, tag)
            if glbMeshMask: self.IO.writeMeshGLB(mesh, hfRadNumpy, label, prefix, tPath, tag)
            # mask
            if csvMask:   self.IO.writePointCloudCSV(centers, shadowRadNumpy, tPath, 'shadowMask', tag, 'shadowMaskRad')
            if vtpPCMask: self.IO.writePointCloudVTP(centers, shadowRadNumpy, 'shadowMask', 'shadowMaskRad', tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, shadowRadNumpy, 'shadowMask', 'shadowMaskRad', tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, shadowRadNumpy, 'shadowMask', 'shadowMaskRad', tPath, tag)

        if want_all:
            prefix = 'HF_allSources'; label = '$MW/m^2$'
            if csvMask:   self.IO.writePointCloudCSV(centers, hfAllNumpy, tPath, label, tag, prefix)
            if vtpPCMask: self.IO.writePointCloudVTP(centers, hfAllNumpy, label, prefix, tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, hfAllNumpy, label, prefix, tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, hfAllNumpy, label, prefix, tPath, tag)

        if want_shPC and not want_opt:
            # if user asked only for a shadow PC export
            prefix = 'shadowMask'; label = 'shadowMask'
            if csvMask:   self.IO.writePointCloudCSV(centers, shadowNumpy, tPath, label, tag, prefix)
            if vtpPCMask: self.IO.writePointCloudVTP(centers, shadowNumpy, label, prefix, tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, shadowNumpy, label, prefix, tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, shadowNumpy, label, prefix, tPath, tag)

        if want_pdir:
            prefix = 'powerDir'; label = 'powerDir'
            if csvMask:   self.IO.writePointCloudCSV(centers, powerDirNumpy, tPath, label, tag, prefix)
            if vtpPCMask: self.IO.writePointCloudVTP(centers, powerDirNumpy, label, prefix, tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, powerDirNumpy, label, prefix, tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, powerDirNumpy, label, prefix, tPath, tag)

        if want_bdn:
            prefix = 'bdotn'; label = '$\\hat{b} \\cdot \\hat{n}$'
            if csvMask:   self.IO.writePointCloudCSV(centers, bdotnNumpy, tPath, label, tag, prefix)
            if vtpPCMask: self.IO.writePointCloudVTP(centers, bdotnNumpy, label, prefix, tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, bdotnNumpy, label, prefix, tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, bdotnNumpy, label, prefix, tPath, tag)

        if want_psi:
            prefix = 'psiN'; label = '$\\psi_N$'
            if csvMask:   self.IO.writePointCloudCSV(centers, psiNumpy, tPath, label, tag, prefix)
            if vtpPCMask: self.IO.writePointCloudVTP(centers, psiNumpy, label, prefix, tPath, tag)
            if vtpMeshMask:self.IO.writeMeshVTP(mesh, psiNumpy, label, prefix, tPath, tag)
            if glbMeshMask:self.IO.writeMeshGLB(mesh, psiNumpy, label, prefix, tPath, tag)

        if want_norm:
            prefix='NormGlyph'; header="X,Y,Z,Nx,Ny,Nz"; label="N"
            if csvMask:   self.IO.writeGlyphCSV(centers, normNumpy, tPath, prefix, header, tag)
            if vtpPCMask: self.IO.writeGlyphVTP(centers, normNumpy, label, prefix, tPath, tag)

        if want_B:
            prefix='BfieldGlyph'; header="X,Y,Z,Bx,By,Bz"; label="B [T]"
            if csvMask:   self.IO.writeGlyphCSV(centers, bFieldNumpy, tPath, prefix, header, tag)
            if vtpPCMask: self.IO.writeGlyphVTP(centers, bFieldNumpy, label, prefix, tPath, tag)

        print("Wrote combined pointclouds"); log.info("Wrote combined pointclouds")
        return


    def combineTimeSteps(self, runList, t):
        """
        save composite csv from each timestep into a single directory for
        making paraview movies
        """
        movieDir = self.MHD.shotPath + 'paraview/'
        tPath = self.MHD.shotPath + self.tsFmt.format(t) + '/'

        vtpMeshMask = self.IO.vtpMeshMask
        vtpPCMask = self.IO.vtpPCMask
        csvMask = self.IO.csvMask
        glbMeshMask = self.IO.glbMeshMask


        #first try to make new directory
        tools.makeDir(movieDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        tStr = self.tsFmt.format(t)
        if 'hfOpt' in runList:
            #TO DO - make all the below runList options use the masks
            if csvMask == True:
                src = tPath + 'HF_optical_all.csv'
                dest = movieDir + 'hfOptical_'+tStr+'.csv'
                shutil.copyfile(src,dest)
            if vtpMeshMask == True:
                src = tPath + '/paraview/HF_optical_all_mesh.vtp'
                dest = movieDir + 'hfOptical_mesh_'+tStr+'.vtp'.format(t)
                shutil.copyfile(src,dest)
                src = tPath + '/paraview/shadowMask_all_mesh.vtp'
                dest = movieDir + 'shadowMask_mesh_'+tStr+'.vtp'
                shutil.copyfile(src,dest)
            if glbMeshMask == True:
                src = tPath + '/paraview/HF_optical_all_mesh.glb'
                dest = movieDir + 'hfOptical_mesh_'+tStr+'.glb'.format(t)
                shutil.copyfile(src,dest)
                src = tPath + '/paraview/shadowMask_all_mesh.glb'
                dest = movieDir + 'shadowMask_mesh_'+tStr+'.glb'
                shutil.copyfile(src,dest)

        if 'shadowPC' in runList:
            src = tPath + 'shadowMask_all.csv'
            dest = movieDir + 'shadowMask_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
        if 'pwrDir' in runList:
            src = tPath + 'powerDir_all.csv'
            dest = movieDir + 'powerDir_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
        if 'bdotn' in runList:
            src = tPath + 'bdotn_all.csv'
            dest = movieDir + 'bdotn_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
        if 'psiN' in runList:
            src = tPath + 'psiN_all.csv'
            dest = movieDir + 'psiN_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
        if 'norm' in runList:
            src = tPath + 'NormGlyph_all.csv'
            dest = movieDir + 'NormGlyph_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
        if 'B' in runList:
            src = tPath + 'BfieldGlyph_all.csv'
            dest = movieDir + 'BfieldGlyph_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
        if 'hfGyro' in runList:
            src = tPath + 'HF_gyro_all.csv'
            dest = movieDir + 'hfGyro_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
            src = tPath + 'HF_allSources_all.csv'
            dest = movieDir + 'hfAll_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)
            src = tPath + 'shadowMaskGyro_all.csv'
            dest = movieDir + 'shadowMaskGyro_'+tStr+'.csv'.format(t)
            shutil.copyfile(src,dest)

        #set tree permissions
        tools.recursivePermissions(movieDir, self.UID, self.GID, self.chmod)

        return

    def psiPC(self, PFC):
        """
        creates a poloidal flux, psi, point cloud on tile surface
        you need to have run the cad, mhd, and hf initialization processes
        before running this function
        """
        #PFC.shadowed_mask = np.zeros((len(PFC.shadowed_mask)))
        self.getPsiEverywhere(PFC, PFC.tag)

        print("Completed psiN calculation")
        log.info("Completed psiN calculation")

        return

    def getPsiEverywhere(self, PFC, save2File=True, tag=None):
        """
        get psi all over the PFC (including shadowed regions).
        """
        #Run MAFOT laminar for 3D plasmas
        if self.plasma3D.plasma3Dmask:
            print('Solving for 3D plasmas with MAFOT')
            log.info('Solving for 3D plasmas with MAFOT')
            self.plasma3D.updatePointsFromCenters(PFC.centers)
            self.plasma3D.launchLaminar(self.NCPUs, tag = 'psiOnly')
            self.plasma3D.cleanUp(tag = 'psiOnly')      # removes the MAFOT log files
            invalid = self.plasma3D.checkValidOutput()    # this does not update self.plasma3D.psimin
            if(np.sum(invalid) > 0): 
                print('****** WARNING *******')
                print('psimin could not be computed for all points.')
                print('Failed points will have psimin = 10.')
                print('Reason: they are most likely outside the M3D-C1 simulation domain.\n')
            PFC.psimin = self.plasma3D.psimin     # this defines and declares PFC.psimin
            
        #get psi from gfile for 2D plasmas
        else:
            print('Solving for 2D plasmas with EFIT')
            log.info('Solving for 2D plasmas with EFIT')
            PFC.shadowed_mask = np.zeros((len(PFC.shadowed_mask)))
            self.MHD.psi2DfromEQ(PFC)

        if PFC.mergedPFCs == True:
            print("Scattering back to PFC objects")
            PFC.scatter_back('psiN', self.IO, self.MHD.shotPath)
        else:
            prefix = 'psiN'
            label = '$\psi_N$'
            path = PFC.controlfilePath
            tag=None
            if self.IO.csvMask == True:
                self.IO.writePointCloudCSV(PFC.centers,PFC.psimin,path,label,tag,prefix)
            if self.IO.vtpPCMask == True:
                self.IO.writePointCloudVTP(PFC.centers,PFC.psimin,label,prefix,path,tag)
            if self.IO.vtpMeshMask == True:
                self.IO.writeMeshVTP(PFC.mesh, PFC.psimin, label, prefix, path, tag)
            if self.IO.glbMeshMask == True:
                self.IO.writeMeshGLB(PFC.mesh, PFC.psimin, label, prefix, path, tag)

        return


    #--- Default I/O and dicts ---

    def getDefaultDict(self):
        """
        returns an empty dict with each inputFile parameter
        """
        emptyDict = {
                    'gridRes': None,
                    'shot':None,
                    'tmin':None,
                    'tmax':None,
                    'traceLength': None,
                    'dpinit': None,
                    'psiMult': None,
                    'BtMult': None,
                    'IpMult':None,
                    'dataPath': None,
                    'hfMode': None,
                    'lqCN': None,
                    'lqCF': None,
                    'lqPN': None,
                    'lqPF': None,
                    'lqCNmode': None,
                    'lqCFmode': None,
                    'lqPNmode': None,
                    'lqPFmode': None,
                    'S': None,
                    'fracCN': None,
                    'fracCF': None,
                    'fracPN': None,
                    'fracPF': None,
                    'fracUI': None,
                    'fracUO': None,
                    'fracLI': None,
                    'fracLO': None,
                    'P': None,
                    'radFrac' : None,
                    'qBG' : None,
                    'fG' : None,
                    'qFilePath' : None,
                    'qFileTag' : None,
                    'rzqFile' : None,
                    'N_gyroSteps': None,
                    'gyroTraceLength': None,
                    'gyroT_eV': None,
                    'N_vSlice': None,
                    'N_vPhase': None,
                    'N_gyroPhase': None,
                    'ionMassAMU': None,
                    'vMode': None,
                    'ionFrac': None,
                    'gyroSources': None,
                    'OFtMin': None,
                    'OFtMax': None,
                    'deltaT': None,
                    'writeDeltaT': None,
                    'STLscale': None,
                    'meshMinLev': None,
                    'meshMaxLev': None,
                    'material': None,
                    'N_gyroSteps': None,
                    'gyroTraceLength': None,
                    'gyroT_eV': None,
                    'N_vSlice': None,
                    'N_vPhase': None,
                    'N_gyroPhase': None,
                    'ionMassAMU': None,
                    'vMode': None,
                    'ionFrac': None,
                    'gyroSources': None,
                    'radFile':None,
                    'phiMin':None,
                    'phiMax':None,
                    'Ntor':None,
                    'Nref':None,
                    'plasma3Dmask':None,
                    'itt':None,
                    'response':None,
                    'selectField':None,
                    'useIcoil':None,
                    'sigma':None,
                    'charge':None,
                    'Ekin':None,
                    'Lambda':None,
                    'Mass':None,
                    'loadHF':None,
                    'loadBasePath':None,
                    'NCPUs':None,
                    'Lcmin':None,
                    'lcfs':None,
                    'teProfileData':None,
                    'neProfileData':None,
                    'kappa':None,
                    'model':None,
                    }
        return emptyDict

    def loadInputs(self, inFile=None):
        """
        loads defaults from file rather than from GUI
        """
        self.setInitialFiles()

        if inFile is not None:
            self.infile = inFile

        tools.initializeInput(self.MHD, self.infile)
        tools.initializeInput(self.CAD, self.infile)
        tools.initializeInput(self.HF, self.infile)
        tools.initializeInput(self.GYRO, self.infile)
        tools.initializeInput(self.FIL, self.infile)
        tools.initializeInput(self.RAD, self.infile)
        tools.initializeInput(self.OF, self.infile)
        tools.initializeInput(self.plasma3D, self.infile)
        tools.initializeInput(self.hf3D, self.infile)
        tools.initializeInput(self.FEM, self.infile)

        inputDict = {
                    'shot': self.MHD.shot,
                    'tmin': self.MHD.tmin,
                    'tmax': self.MHD.tmax,
                    'traceLength': self.MHD.traceLength,
                    'dpinit': self.MHD.dpinit,
                    'psiMult':self.MHD.psiMult,
                    'BtMult':self.MHD.BtMult,
                    'IpMult':self.MHD.IpMult,                    
                    'gridRes': self.CAD.gridRes,
                    'hfMode': self.HF.hfMode,
                    'lqEich': self.HF.lqCN,
                    'S': self.HF.S,
                    'SMode': self.HF.SMode,
                    'lqCN': self.HF.lqCN,
                    'lqCF': self.HF.lqCF,
                    'lqPN': self.HF.lqPN,
                    'lqPF': self.HF.lqPF,
                    'lqTopHat': self.HF.lqCN,
                    'lqCNmode': self.HF.lqCNmode,
                    'lqCFmode': self.HF.lqCFmode,
                    'lqPNmode': self.HF.lqPNmode,
                    'lqPFmode': self.HF.lqPFmode,
                    'fracCN': self.HF.fracCN,
                    'fracCF': self.HF.fracCF,
                    'fracPN': self.HF.fracPN,
                    'fracPF': self.HF.fracPF,
                    'P': self.HF.P,
                    'radFrac' : self.HF.radFrac,
                    'fracUI':self.HF.fracUI,
                    'fracUO':self.HF.fracUO,
                    'fracLI':self.HF.fracLI,
                    'fracLO':self.HF.fracLO,
                    'qBG' : self.HF.qBG,
                    'fG' : self.HF.fG,
                    'qFilePath': self.HF.qFilePath,
                    'qFileTag': self.HF.qFileTag,
                    'rzqFile' : self.HF.rzqFile,
                    'OFtMin': self.OF.OFtMin,
                    'OFtMax': self.OF.OFtMax,
                    'deltaT': self.OF.deltaT,
                    'writeDeltaT': self.OF.writeDeltaT,
                    'STLscale': self.OF.STLscale,
                    'meshMinLevel': self.OF.meshMinLevel,
                    'meshMaxLevel': self.OF.meshMaxLevel,
                    'material': self.OF.material,
                    'N_gyroSteps': self.GYRO.N_gyroSteps,
                    'gyroTraceLength': self.GYRO.gyroTraceLength,
                    'gyroT_eV': self.GYRO.gyroT_eV,
                    'N_vSlice': self.GYRO.N_vSlice,
                    'N_vPhase': self.GYRO.N_vPhase,
                    'N_gyroPhase': self.GYRO.N_gyroPhase,
                    'ionMassAMU': self.GYRO.ionMassAMU,
                    'vMode': self.GYRO.vMode,
                    'ionFrac': self.GYRO.ionFrac,
                    'gyroSources': self.GYRO.gyroSources,
                    'radFile': self.RAD.radFile,
                    'phiMin':self.RAD.phiMin,
                    'phiMax':self.RAD.phiMax,
                    'Ntor':self.RAD.Ntor,
                    'Nref':self.RAD.Nref,
                    'plasma3Dmask':self.plasma3D.plasma3Dmask,
                    'itt':self.plasma3D.itt,
                    'response':self.plasma3D.response,
                    'selectField':self.plasma3D.selectField,
                    'useIcoil':self.plasma3D.useIcoil,
                    'sigma':self.plasma3D.sigma,
                    'charge':self.plasma3D.charge,
                    'Ekin':self.plasma3D.Ekin,
                    'Lambda':self.plasma3D.Lambda,
                    'Mass':self.plasma3D.Mass,
                    'loadHF':self.plasma3D.loadHF,
                    'loadBasePath':self.plasma3D.loadBasePath,
                    'NCPUs':self.plasma3D.NCPUs,
                    'Lcmin':self.hf3D.Lcmin,
                    'lcfs':self.hf3D.lcfs, 
                    'teProfileData':self.hf3D.teProfileData,
                    'neProfileData':self.hf3D.neProfileData,
                    'kappa':self.hf3D.kappa,
                    'model':self.hf3D.model
                    }
        print("Loaded inputs")

        return inputDict

    def getCurrentInputs(self):
        """
        loads current values for input file variables
        """
        inputDict = {
                    'shot': self.MHD.shot,
                    'tmin': self.MHD.tmin,
                    'tmax': self.MHD.tmax,
                    'traceLength': self.MHD.traceLength,
                    'dpinit': self.MHD.dpinit,
                    'psiMult':self.MHD.psiMult,
                    'BtMult':self.MHD.BtMult,
                    'IpMult':self.MHD.IpMult,
                    'gridRes': self.CAD.gridRes,
                    'hfMode': self.HF.hfMode,
                    'lqEich': self.HF.lqCN,
                    'S': self.HF.S,
                    'SMode': self.HF.SMode,
                    'lqCN': self.HF.lqCN,
                    'lqCF': self.HF.lqCF,
                    'lqPN': self.HF.lqPN,
                    'lqPF': self.HF.lqPF,
                    'lqTopHat': self.HF.lqCN,
                    'lqCNmode': self.HF.lqCNmode,
                    'lqCFmode': self.HF.lqCFmode,
                    'lqPNmode': self.HF.lqPNmode,
                    'lqPFmode': self.HF.lqPFmode,
                    'fracCN': self.HF.fracCN,
                    'fracCF': self.HF.fracCF,
                    'fracPN': self.HF.fracPN,
                    'fracPF': self.HF.fracPF,
                    'P': self.HF.P,
                    'radFrac' : self.HF.radFrac,
                    'fracUI':self.HF.fracUI,
                    'fracUO':self.HF.fracUO,
                    'fracLI':self.HF.fracLI,
                    'fracLO':self.HF.fracLO,
                    'qBG' : self.HF.qBG,
                    'fG' : self.HF.fG,
                    'qFilePath': self.HF.qFilePath,
                    'qFileTag': self.HF.qFileTag,
                    'rzqFile': self.HF.rzqFile,
                    'OFtMin': self.OF.OFtMin,
                    'OFtMax': self.OF.OFtMax,
                    'deltaT': self.OF.deltaT,
                    'writeDeltaT': self.OF.writeDeltaT,
                    'STLscale': self.OF.STLscale,
                    'meshMinLevel': self.OF.meshMinLevel,
                    'meshMaxLevel': self.OF.meshMaxLevel,
                    'material': self.OF.material,
                    'N_gyroSteps': self.GYRO.N_gyroSteps,
                    'gyroTraceLength': self.GYRO.gyroTraceLength,
                    'gyroT_eV': self.GYRO.gyroT_eV,
                    'N_vSlice': self.GYRO.N_vSlice,
                    'N_vPhase': self.GYRO.N_vPhase,
                    'N_gyroPhase': self.GYRO.N_gyroPhase,
                    'ionMassAMU': self.GYRO.ionMassAMU,
                    'vMode': self.GYRO.vMode,
                    'ionFrac': self.GYRO.ionFrac,
                    'gyroSources': self.GYRO.gyroSources,
                    'radFile': self.RAD.radFile,
                    'phiMin':self.RAD.phiMin,
                    'phiMax':self.RAD.phiMax,
                    'Ntor':self.RAD.Ntor,
                    'Nref':self.RAD.Nref,
                    'plasma3Dmask':self.plasma3D.plasma3Dmask,
                    'itt':self.plasma3D.itt,
                    'response':self.plasma3D.response,
                    'selectField':self.plasma3D.selectField,
                    'useIcoil':self.plasma3D.useIcoil,
                    'sigma':self.plasma3D.sigma,
                    'charge':self.plasma3D.charge,
                    'Ekin':self.plasma3D.Ekin,
                    'Lambda':self.plasma3D.Lambda,
                    'Mass':self.plasma3D.Mass,
                    'loadHF':self.plasma3D.loadHF,
                    'loadBasePath':self.plasma3D.loadBasePath,
                    'NCPUs':self.plasma3D.NCPUs,
                    'Lcmin':self.hf3D.Lcmin,
                    'lcfs':self.hf3D.lcfs, 
                    'teProfileData':self.hf3D.teProfileData,
                    'neProfileData':self.hf3D.neProfileData,
                    'kappa':self.hf3D.kappa,
                    'model':self.hf3D.model,
                    'meshFEMminRes':self.FEM.meshFEMminRes,
                    'meshFEMmaxRes':self.FEM.meshFEMmaxRes,
                    'elmerDir':self.FEM.elmerDir,
                    'elmerFile':self.FEM.elmerFile,
                    'elmerHEATlib':self.FEM.elmerHEATlib
                    }
        print("Loaded current inputs")

        return inputDict

    def loadPFCParts(self):
        """
        loads parts list from an input file
        """
        with open(self.IntersectFile) as f:
            lines = f.readlines()

        for idx, line in enumerate(lines):
            lines[idx] = line.replace('#','').rstrip()

        return lines[:]

    def loadPFCDefaults(self):
        """
        returns two lists
        1) list of parts to calculate HF on
        2) list of parts to check for intersections with
        """
        parts = pd.read_csv(self.PartsFile, sep=',', comment='#', names=['parts'], skipinitialspace=True)
        intersects = pd.read_csv(self.IntersectFile, sep=',', comment='#', names=['intersects'], skipinitialspace=True)
        return parts['parts'].to_list(), intersects['intersects'].to_list()

    def writePFCs(self,parts,intersects):
        """
        writes a PartsFile and an IntersectFile.  Updates the respective class
        variables so that subsequent function calls reflect these new files.
        """
        allParts = self.loadPFCParts()
        path, default = os.path.split(self.PartsFile)
        self.PartsFile = path + '/userParts.csv'
        self.IntersectFile = path + '/userIntersects.csv'

        with open(self.PartsFile, 'w') as f:
            for name in allParts:
                if name in parts:
                    f.write(name+'\n')
                else:
                    f.write('#'+name+'\n')

        with open(self.IntersectFile, 'w') as f:
            for name in allParts:
                if name in intersects:
                    f.write(name+'\n')
                else:
                    f.write('#'+name+'\n')

        return


    #--- OpenFOAM ---

    def loadOF(self,OFtMin,OFtMax,OFminMeshLevel,OFmaxMeshLevel,
                      OFSTLscale, OFbashrc, OFdeltaT, OFwriteDeltaT, materialSelect):
        """
        loads user OF GUI settings

        OFtMin is when we start OF simulation (can be before HF)
        OFtMax is when we stop OF simulation (can be after HF)
        OFminMeshLev is minimum refinement level for snappyhexmesh (default is 1)
        OFminMeshLev is maximum refinement level for snappyhexmesh (default is 3)
        OFSTLscale is scalar for unit conversion (default is 1)
        OFbashrc is file location on system to source OF binaries / libs
        OFdeltaT is timestep size for FVM simulation.  Defaults to 1ms (0.001s)
        """
        self.OF.OFtMin = float(OFtMin) #[s] for openfoam
        self.OF.OFtMax = float(OFtMax) #[s] for openfoam
        self.OF.meshMinLevel = int(OFminMeshLevel)
        self.OF.meshMaxLevel = int(OFmaxMeshLevel)
        self.OF.STLscale = OFSTLscale
        self.OF.cmd3Dmesh = 'meshAndPatch'
        #the OF bashrc must be sourced before running OF
        self.OF.cmdSourceOF = 'source ' + OFbashrc
        self.OF.OFbashrc = OFbashrc
        self.OF.cmdThermal = 'runThermal'
        self.OF.cmdTprobe = 'runTprobe'
        self.OF.deltaT = float(OFdeltaT) #this comes in [sec]
        self.OF.writeDeltaT = float(OFwriteDeltaT) #this comes in [sec]
        self.OF.material = materialSelect

        print("Material Selection: "+self.OF.material)

        print("Loaded OF data")
        log.info("Loaded OF data")

        return

    def runOpenFOAM(self):
        """
        sets up and runs OpenFOAM finite volume simulation to determine
        temperature distribution in a PFC.  For use with time varying multiple
        tile cases

        Note that HEAT must have been run BEFORE running this function
        The user cannot change PFC, MHD, CAD, or HF, settings in GUI in between
        running HEAT to generate HF and this function, because it depends on all
        of them. Reads heat flux csv from HEAT tree for each timestep and
        applies as boundary condition for openFOAM FVM simulation, so you need
        the heat flux csvs to be in the tree already

        HEAT calculates heat fluxes at discrete timesteps.  Because these may
        not always be uniform in time (ie every millisecond) we need a method
        for resolving the discrepancies between the various timesteps we have in
        HEAT.  We have three different timestep variables that we need to juggle:
            - MHD.timesteps: timesteps where we have equilibrium data
            - PFC.timesteps: subset of MHD timesteps for each PFC tile
            - OF.timesteps: evenly spaced timesteps for FVM analyses

        In order to couple these timestep variables, we need to have conditionals
        that tests whether each OF.timestep is:
            - in MHD.timesteps?
            - in PFC.timesteps?

        Additionally, for each OF.timestep, there are three potential cases:
        1)  If the OF.timestep aligns perfectly with a PFC timestep the boundary
            condition is easy, we just assign the heat flux HEAT calculated for that
            timestep as the boundary condition.
        2)  If, however, the OF.timestep is in between PFC.timesteps, then we
            do not write anything and openFOAM linearly interpolates (in time)
            the heat flux between PFC.timesteps
        3) There is one more case, when we are outside of the MHD.timesteps.  This
            can happen when you want to run your FVM simulation for a duration longer
            than the plasma discharge (ie to capture the exponential decay of
            temperature within a PFC for a few minutes after a shot ends).  Here
            I just assign a heat flux of 0 MW/m^2 to the boundary
        """
        #uncomment this line when you need to patch the HF before T calc
        #input("Patch files then press enter to continue...")

        print('Setting Up OF run')
        log.info('Setting Up OF run')

        #set up base OF directory for this discharge
        self.OF.OFDir = self.MHD.shotPath + 'openFoam'
        self.OF.caseDir = self.MHD.shotPath + 'openFoam/heatFoam'
        tools.makeDir(self.OF.caseDir, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)

        #set up directory for all .foam files
        #self.OF.allFoamsDir = self.OF.caseDir + '/allFoams'
        #tools.makeDir(self.OF.allFoamsDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        #set up OF parts for each PFC part
        for PFC in self.PFCs:
            #check if PFC is a gyroSource Plane
            if hasattr(self.GYRO, 'gyroSourceTag'):
                #check if PFC is a gyroSource Plane
                if self.GYRO.gyroSourceTag == 'allROI':
                    pass
                else:
                    if PFC.name in self.GYRO.gyroSources:
                        print("Not including "+PFC.name+" in thermal analysis")
                        continue

            print("Running openFOAM for PFC: "+PFC.name)
            log.info("Running openFOAM for PFC: "+PFC.name)
            partDir = self.OF.caseDir + '/' + PFC.name.replace(" ", "_")
            self.OF.partDir = partDir
            #replace spaces with underscores to prevent pyfoam from reading args as
            # space delimited list
            self.OF.partName = PFC.name.replace(" ", "_")

            #tools.makeDir(partDir, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)
            #copy heatFoam template directory to this location
            try:
                tools.copytree(self.OF.templateCase, partDir)
                if self.OF.OFtMin != 0:
                    t0 = partDir+'/0'
                    t0new = partDir+'/{:f}'.format(self.OF.OFtMin).rstrip('0').rstrip('.')
                    tools.copytree(t0, t0new)
                    shutil.rmtree(t0)
            except OSError as e:
                print('COULD NOT COPY TEMPLATE DIRECTORY!  Aborting!')
                print(e)
                return

            #Set up partDir with user selected material properties
            print("Material Selection: "+self.OF.material)
            matDst = partDir + '/constant/'
            if self.OF.material == "ATJ":
                print("Using ATJ Material Template")
                log.info("Using ATJ Material Template")
                matSrc = self.OF.materialDir + '/ATJ/'
            elif self.OF.material == "MOLY":
                print("Using MOLY Material Template")
                log.info("Using MOLY Material Template")
                matSrc = self.OF.materialDir + '/MOLY/'
            elif self.OF.material == "TUNG":
                print("Using TUNG Material Template")
                log.info("Using TUNG Material Template")
                matSrc = self.OF.materialDir + '/TUNG/'
            elif self.OF.material == "TUNG_SPARC":
                print("Using TUNG_SPARC Material Template")
                log.info("Using TUNG_SPARC Material Template")
                matSrc = self.OF.materialDir + '/TUNG_SPARC/'
            elif self.OF.material == "WHA":
                print("Using WHA Material Template")
                log.info("Using WHA Material Template")                
                matSrc = self.OF.materialDir + '/WHA/'
            elif self.OF.material == "SGL":
                print("Using SGL Material Template")
                log.info("Using SGL Material Template")                
                matSrc = self.OF.materialDir + '/SGLR6510/'
            else: #SGL6510
                print("COULD NOT FIND MATERIAL!")
                log.info("COULD NOT FIND MATERIAL!")
                print("Using USER Material Template")
                log.info("Using USER Material Template")                
                matSrc = self.OF.materialDir + '/USER/'
            try:
                shutil.copyfile(matSrc+'DT', matDst+'DT')
                shutil.copyfile(matSrc+'thermCond', matDst+'thermCond')
            except OSError as e:
                print('COULD NOT COPY MATERIAL PROPERTIES!  Aborting!')
                print(e)

            #set up timesteps
            N_t = int((self.OF.OFtMax - self.OF.OFtMin) / self.OF.deltaT)
            OFtimesteps = np.round(np.linspace(self.OF.OFtMin, self.OF.OFtMax, N_t+1), self.tsSigFigs)

            #create symbolic link to STL file
            print("Creating openFOAM symlink to STL")
            log.info("Creating openFOAM symlink to STL")
            #old method used ROIGridRes, but this can sometimes create non-watertight mesh
            #new method uses standard FreeCAD mesher (not mefisto) to make better volume mesh
            #standard meshing algorithm
            PFC.OFpart = PFC.name + "___" + "standard"

            if self.CAD.STLpath[-1] == '/':
                stlfile = self.CAD.STLpath + PFC.OFpart +".stl"
            else:
                stlfile = self.CAD.STLpath +'/'+ PFC.OFpart +".stl"

            #if the standard mesh for this part doesn't exist, create it using freecad
            #in the HEAT CAD module
            if os.path.exists(stlfile)==True:
                print("Mesh already exists.")
                log.info("Mesh already exists.")
            else:
                partIdx = np.where(self.CAD.ROI == PFC.name)[0][0]
                part = self.CAD.ROIparts[partIdx]
                meshSTL = self.CAD.part2meshStandard(part)
                self.CAD.writeMesh2file(meshSTL, PFC.name, path=self.CAD.STLpath, resolution='standard')

            #replace spaces with underscores to prevent pyfoam from reading args as
            # space delimited list
            PFC.OFpart = PFC.OFpart.replace(" ", "_")
            triSurfaceLocation = partDir+'/constant/triSurface/' + PFC.OFpart +".stl"

            #copy STL to triSurface location (copyfile for S3 compat)
            shutil.copyfile(stlfile, triSurfaceLocation)

            #Test STL to make sure it is watertight (open3d)
            #meshInQuestion = o3d.io.read_triangle_mesh(stlfile)
            #validMesh = meshInQuestion.is_watertight()
            #Test STL to make sure it is watertight (freecad)
            partIdx = np.where(self.CAD.ROI == PFC.name)[0][0]
            part = self.CAD.ROIparts[partIdx]
            meshInQuestion = self.CAD.readMesh(stlfile)
            validMesh = meshInQuestion.isSolid()

            if validMesh == True:
                print("STL file for "+PFC.name+" is watertight")
                log.info("STL file for "+PFC.name+" is watertight")
            else:
                #try to repair
                print("Trying to repair leaky mesh...")
                log.info("Trying to repair leaky mesh...")
                p = partDir+'/constant/triSurface/'
                meshName = PFC.OFpart +".stl"
                #validMesh = self.CAD.repairMeshOpen3D(p, meshName, meshName, part, 'dummySTP.step')
                validMesh = self.CAD.repairMeshFreeCAD(p, meshName, PFC.name, meshInQuestion)
                if validMesh == True:
                    print("repaired mesh...")
                    log.info("repaired mesh...")                   

                else:
                    print("\n=====================================================")
                    print("WARNING!!!  STL file for "+PFC.name+" is NOT watertight.")
                    print("OpenFOAM thermal analysis will probably fail")
                    print("=====================================================\n")
                    log.info("\n=====================================================")
                    log.info("WARNING!!!  STL file for "+PFC.name+" is NOT watertight.")
                    log.info("OpenFOAM thermal analysis will probably fail")
                    print("=====================================================\n")

            #update blockmesh bounds for each PFC and
            # give 10mm of clearance on each side
            self.OF.xMin = (PFC.centers[:,0].min() - 0.01)*1000.0
            self.OF.xMax = (PFC.centers[:,0].max() + 0.01)*1000.0
            self.OF.yMin = (PFC.centers[:,1].min() - 0.01)*1000.0
            self.OF.yMax = (PFC.centers[:,1].max() + 0.01)*1000.0
            self.OF.zMin = (PFC.centers[:,2].min() - 0.01)*1000.0
            self.OF.zMax = (PFC.centers[:,2].max() + 0.01)*1000.0

            #find a point inside the mesh
            #old method:  only works for solid objects (if COM is hollow, fails)
            #self.OF.xMid = (self.OF.xMax-self.OF.xMin)/2.0 + self.OF.xMin
            #self.OF.yMid = (self.OF.yMax-self.OF.yMin)/2.0 + self.OF.yMin
            #self.OF.zMid = (self.OF.zMax-self.OF.zMin)/2.0 + self.OF.zMin

            #another old method. steps just inside a mesh face
            #find point 2mm inside of face 0 (this will fail if part width is < 2mm)
            #smallStep = 0.002 #2mm
            #self.OF.xMid = (PFC.centers[0,0] - smallStep*PFC.norms[0,0])*1000.0
            #self.OF.yMid = (PFC.centers[0,1] - smallStep*PFC.norms[0,1])*1000.0
            #self.OF.zMid = (PFC.centers[0,2] - smallStep*PFC.norms[0,2])*1000.0

            #find a point inside the mesh via ray tracing to other mesh faces
            #in the PFC.
            #set up other faces in PFC (targets)
            print("Calculating volume mesh internal coordinate")
            log.info("Calculating volume mesh internal coordinate")
            targetPoints = []
            partIdx = np.where(self.CAD.ROI == PFC.name)[0][0]
            part = self.CAD.ROImeshes[partIdx]
            norms = PFC.norms
            ctrs = PFC.centers

            #partIdx = np.where(np.array(self.CAD.intersectList) == PFC.name)[0][0]
            #part = self.CAD.intersectMeshes[partIdx]
            #norms = self.CAD.intersectNorms[partIdx]
            #ctrs = self.CAD.intersectCtrs[partIdx]/1000.0

            for face in part.Facets:
                targetPoints.append(face.Points)
            targetPoints = np.asarray(targetPoints) / 1000.0 #scale to m
            targetPoints = targetPoints[1:,:] #omit launch point face
            #find all intersections with self along inverse normal
            tools.q1 = ctrs[0,:]
            tools.q2 = ctrs[0,:] -10.0 * norms[0,:]            
            tools.p1 = targetPoints[:,0,:] #point 1 of mesh triangle
            tools.p2 = targetPoints[:,1,:] #point 2 of mesh triangle
            tools.p3 = targetPoints[:,2,:] #point 3 of mesh triangle
            tools.Nt = len(tools.p1)
            mask = tools.intersectTestSingleRay()
            use = np.where(mask == True)[0] + 1 #account for face we deleted above (index 0)
            #find which one of these intersection faces is closest
            dist = np.linalg.norm(ctrs[0,:] - ctrs[use,:], axis=1)           
            #smallStep = dist[closest] / 2.0
            closest = use[np.argmin(dist)]
            #find the distance that face is along the inverse normal
            distAlongNorm = np.dot((ctrs[0,:] - ctrs[closest,:]), norms[0,:])
            #put a point between those two planes
            smallStep = distAlongNorm / 2.0

            #place midpoint between index 0 and the closest intersection face along inverse normal
            self.OF.xMid = (ctrs[0,0] - smallStep*norms[0,0])*1000.0
            self.OF.yMid = (ctrs[0,1] - smallStep*norms[0,1])*1000.0
            self.OF.zMid = (ctrs[0,2] - smallStep*norms[0,2])*1000.0

            #dynamically write template variables to templateVarFile
            print("Building openFOAM templates and shell scripts")
            log.info("Building openFOAM templates and shell scripts")
            templateVarFile = partDir + '/system/templateVariables'
            STLpart = PFC.OFpart +".stl"
            self.OF.writeOFtemplateVarFile(templateVarFile, STLpart)
            
            self.OF.writeShellScript()

            #create OF dictionaries from templates
            self.OF.createDictionaries(self.OF.templateDir,
                                       partDir,
                                       templateVarFile,
                                       #stlfile)
                                       triSurfaceLocation)

            #generate 3D volume mesh or copy from file
            print("Generating volume mesh")
            log.info("Generating volume mesh")
            self.OF.generate3Dmesh(PFC.OFpart, self.CAD.overWriteMask)

            qDiv = np.zeros((len(PFC.centers)))
            ctrs = copy.copy(PFC.centers)*1000.0


            #in the future if you want adaptive meshing for each timestep this whole
            #section needs to be revised.  for now, you can only have 1 boundary
            #condition set of points.  We take timestep 0 from the PFC
            #and get the points from a file.  This enables user to swap out
            #HF_allSources.csv in between HF and Temp runs, to create non-uniform mesh
            #get the mesh centers from one of the timesteps
            t = self.MHD.timesteps[0]
            HFcsv = self.MHD.shotPath + self.tsFmt.format(t) + '/' + PFC.name + '/HF_allSources.csv'
            OFcenters = pd.read_csv(HFcsv).iloc[:,0:3].values

            #cycle through timesteps and get HF data from HEAT tree
            for t in OFtimesteps:
                print("openFOAM timestep: {:f}".format(t))
                log.info("openFOAM timestep: {:f}".format(t))

                #make timeDir folder and copy HF
                timeDir = partDir + '/{:f}'.format(t).rstrip('0').rstrip('.')
                #timestep field prescriptions
                HFt0 = partDir + '/{:f}'.format(self.OF.OFtMin).rstrip('0').rstrip('.')+'/HF'
                HFtStep = partDir + '/{:f}'.format(t).rstrip('0').rstrip('.')+'/HF'

                if self.OF.OFtMin != t:
                    try:
                        #shutil.copytree(t0new,timeDir
                        tools.makeDir(timeDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
                        #shutil.copy(HFt0, HFtStep)
                        shutil.copyfile(HFt0, os.path.join(timeDir, os.path.basename(HFt0)))

                    except:
                        print("***")
                        print("Could not create OF directory for timestep {:d}".format(t))
                        print("(this is expected for t=0)")
                        print("***")
                        pass


                # determine heat flux boundary condition
                if (t in PFC.timesteps) and (t in self.MHD.timesteps):
                    #we explicitly calculated HF for this timestep
                    print("OF.timestep: {:f} in PFC.timesteps".format(t))
                    log.info("OF.timestep: {:f} in PFC.timesteps".format(t))
                    HFcsv = self.MHD.shotPath + self.tsFmt.format(t) + '/' + PFC.name + '/HF_allSources.csv'
                    qDiv = pd.read_csv(HFcsv)['$MW/m^2$'].values #this is the HF column header in the CSV file
                    #OFcenters = pd.read_csv(HFcsv).iloc[:,0:3].values
                    #write boundary condition
                    print("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    log.info("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    self.HF.write_openFOAM_boundary(OFcenters,qDiv,partDir,t)
                elif (t < self.MHD.timesteps.min()) or (t > self.MHD.timesteps.max()):
                    #apply zero HF outside of discharge domain (ie tiles cooling)
                    print("OF.timestep: {:f} outside MHD domain".format(t))
                    log.info("OF.timestep: {:f} outside MHD domain".format(t))

                    qDiv = np.zeros((len(OFcenters)))
                    #write boundary condition
                    print("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    log.info("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    self.HF.write_openFOAM_boundary(OFcenters,qDiv,partDir,t)
                else:
                    #boundary using last timestep that we calculated a HF for
                    #(basically a heaviside function in time)
                    #print("OF.timestep: {:d} using heaviside from last PFC.timestep".format(t))
                    #log.info("OF.timestep: {:d} using heaviside from last PFC.timestep".format(t))
                    #write boundary condition
                    #print("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    #self.HF.write_openFOAM_boundary(ctrs,qDiv,partDir,t)

                    #openFOAM linear interpolation in time using timeVaryingMappedFixedValue
                    print("HF being linearly interpolated by OF at this t")
                    pass


            #run openfoam thermal analysis using heatFoam solver
            self.OF.runThermalAnalysis()
            print("thermal analysis complete...")
            log.info("thermal analysis complete...")

            #THIS DOES NOT WORK.  FUTURE WORK TO CREATE SINGLE .foam FILE
            ##build single .foam directory structure
            #print("Building allFoams directory")
            #log.info("Building allFoams directory")
            #for t in OFtimesteps:
            #    #make timeDir folder and copy HF
            #    timeDir = partDir + '/{:f}'.format(t).rstrip('0').rstrip('.')
            #    #make the allFoam timestep directory for this timestep
            #    allFoamTimeDir = self.OF.allFoamsDir + '/{:f}'.format(t).rstrip('0').rstrip('.')
            #    tools.makeDir(allFoamTimeDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)

        #set tree permissions
        tools.recursivePermissions(self.OF.OFDir, self.UID, self.GID, self.chmod)
        print("openFOAM run completed.")
        log.info("openFOAM run completed.")
        return


    def loadElmer(self):
        """
        loads an elmer FEM file and prepares for an elmer simulation
        """ 
        self.FEM.loadElmerFile()        
        return
    
    def runElmerFEM(self, meshAlg='gmsh'):
        """
        runs an Elmer FEM simulation
        """
        #build Elmer FEM output directory
        self.FEM.elmerOutDir = self.MHD.shotPath + 'elmer/'

        HFvars = ['hfOpt', 'hfRad', 'hfGyro', 'hfFil']
        #only clobber if we created a new HF
        if len(set(HFvars).intersection(self.runList)) > 0:
            tools.makeDir(self.FEM.elmerOutDir, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)
        else:
            tools.makeDir(self.FEM.elmerOutDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        
        #initialize list of variables we accept as none in elmerFile
        noneArray = ['None', 'NA', 'none', 'NONE', 'na', '']

        #loop thru PFCs, building meshes
        for PFC in self.PFCs:
            #parameters from elmerFile
            params = self.FEM.elmerData[PFC.name]

            #user did not supply a mesh file
            if params['meshFile'] in noneArray:
                print("Creating new FEM mesh.")
                log.info("Creating new FEM mesh.")
                partIdx = self.CAD.ROIList.index(PFC.name)
                part = self.CAD.ROIparts[partIdx]

                if meshAlg == 'gmsh':
                    print("Using GMSH mesh algorithm")
                    log.info("Using GMSH mesh algorithm")
                    meshName = PFC.name + '_GMSH_{:0.3f}mm_{:0.3f}mm'.format(self.FEM.meshFEMminRes, self.FEM.meshFEMmaxRes)
                    mesh = self.CAD.createFEMmeshGmsh(part, minLength=self.FEM.meshFEMminRes, maxLength=self.FEM.meshFEMmaxRes, name=meshName)
                else:
                    print("Using NETGEN mesh algorithm")
                    log.info("Using NETGEN mesh algorithm")
                    meshName = PFC.name + '_NETGEN_{:0.3f}mm'.format(self.FEM.meshFEMmaxRes)
                    mesh = self.CAD.createFEMmeshNetgen(part, MaxSize=self.FEM.meshFEMmaxRes, name=meshName)

                meshFile = self.FEM.elmerOutDir + meshName + '.unv'
                self.CAD.exportFEMmesh(mesh, meshFile)

            #user supplied a mesh file
            else:
                #check if mesh file from user exists in .unv format
                meshExists = os.path.isfile(self.FEM.elmerDir + params['meshFile'])
                #if mesh .unv does not exist then make it
                if meshExists != True:
                    print("User supplied FEM mesh does not exist.  Creating new mesh.")
                    log.info("User supplied FEM mesh does not exist.  Creating new mesh.")
                    partIdx = self.CAD.ROIList.index(PFC.name)
                    part = self.CAD.ROIparts[partIdx]
                    if meshAlg == 'gmsh':
                        print("Using GMSH mesh algorithm")
                        log.info("Using GMSH mesh algorithm")
                        meshName = PFC.name + '_GMSH_{:0.3f}mm_{:0.3f}mm'.format(self.FEM.meshFEMminRes, self.FEM.meshFEMmaxRes)
                        mesh = self.CAD.createFEMmeshGmsh(part, minLength=self.FEM.meshFEMminRes, maxLength=self.FEM.meshFEMmaxRes, name=meshName)
                    else:
                        print("Using NETGEN mesh algorithm")
                        log.info("Using NETGEN mesh algorithm")
                        meshName = PFC.name + '_NETGEN_{:0.3f}mm'.format(self.FEM.meshFEMmaxRes)
                        mesh = self.CAD.createFEMmeshNetgen(part, MaxSize=self.FEM.meshFEMmaxRes, name=meshName)
                    meshFile = self.FEM.elmerOutDir + meshName + '.unv'
                    self.CAD.exportFEMmesh(mesh, meshFile)
                #mesh exists
                else:
                    print("Using user supplied FEM mesh.")
                    log.info("Using user supplied FEM mesh.")
                    meshName = params['meshFile']  
                    meshFile = self.FEM.elmerDir + meshName         

            #location where we will save the Elmer grid mesh (its a directory)
            meshDir = self.FEM.elmerOutDir + PFC.name
            self.FEM.buildElmerMesh(meshDir, meshFile)

            PFC.meshDir = meshDir
            PFC.meshFile = meshFile
            PFC.meshName = meshName
            PFC.SIFfile = params['SIF']

        #build a timestep array from the SIF       
        self.FEM.buildTimesteps(PFC.SIFfile)

        #loop through Elmer timesteps, assigning HF from MHD timesteps
        #to the Elmer timesteps as necessary.  Interpolates the heat 
        #flux when the Elmer timesteps fall between MHD timesteps.
        #If the Elmer timesteps fall outside of the MHD timestep domain,
        #assigns 0.
        #only assign HF values if we ran a heat flux calculation,
        #otherwise we assume the (node, HF) .dat files are already in the ElmerDir
        if len(set(HFvars).intersection(self.runList)) > 0:
            for tIdx,t in enumerate(self.FEM.ts):
                for PFC in self.PFCs:
                    if t in self.MHD.timesteps:
                        tMHD = t
                    else:
                        tMHD = self.MHD.timesteps[0]
                    tMin = np.min(self.MHD.timesteps)
                    tMax = np.max(self.MHD.timesteps)
                    if t < tMin:
                        #timesteps outside of PFC domain are assigned 0 HF on surface
                        self.FEM.interpolateHFtoMesh(PFC, t, tMin, hfFile=None)
                    elif t > tMax:
                        #timesteps outside of PFC domain are assigned 0 HF on surface
                        self.FEM.interpolateHFtoMesh(PFC, t, tMin, hfFile=None)
                    elif t not in PFC.timesteps:
                        #timesteps within the MHD domain but not an MHD timestep get linear
                        #interpolated to the Elmer timestep
                        idx = np.where(tMHD == self.MHD.timesteps)[0][0]
                        tNext = self.MHD.timesteps[idx+1]
                        tprev = self.MHD.timesteps[idx]
                        pfcDir = self.MHD.shotPath + self.tsFmt.format(tprev) +'/'+PFC.name+'/'
                        hfFile = pfcDir + "HF_allSources.csv"
                        pfcDirNext = self.MHD.shotPath + self.tsFmt.format(tNext) +'/'+PFC.name+'/'
                        hfFileNext= pfcDirNext + "HF_allSources.csv"
                        hfFileNew = self.FEM.interpolateHFinTime(hfFile, hfFileNext, tMHD, tNext, t)
                        self.FEM.interpolateHFtoMesh(PFC, t, tMin, hfFileNew) 
                    else:
                        #Elmer timesteps align with MHD timesteps
                        pfcDir = self.MHD.shotPath + self.tsFmt.format(tMHD) +'/'+PFC.name+'/'
                        hfFile = pfcDir + "HF_allSources.csv"
                        self.FEM.interpolateHFtoMesh(PFC, t, tMin, hfFile)
   

        #loop through PFCs, running Elmer Solvers
        for PFC in self.PFCs:
            #parameters from elmerFile
            params = self.FEM.elmerData[PFC.name]
            #copy the ReX init file
            self.FEM.copyReXinit(PFC)
            #solve the Elmer system
            self.FEM.runElmerSolve(params['SIF'], PFC.name)

        #set tree permissions
        tools.recursivePermissions(self.FEM.elmerOutDir, self.UID, self.GID, self.chmod)
        return

    #--- Plots ---

    def getOFMinMaxPlots(self):
        """
        returns plotly figure that has data from openFOAM postProcessing
        fieldMinMax.dat files for each PFC

        The function that gets the data during the run is located in the
        controlDict for that run, as the fieldMinMax1 function.  It can
        be found in the openFOAMTemplates/templateDicts/ dir for HEAT
        """
        data = []
        pfcNames = []
        for PFC in self.PFCs:
            #check if PFC is a gyroSource Plane
            if hasattr(self.GYRO, 'gyroSourceTag'):
                if self.GYRO.gyroSourceTag == 'allROI' :
                    pass
                else:
                    if PFC.name in self.GYRO.gyroSources:
                        print("Not including "+PFC.name+" in MinMax plots")
                        continue

            partDir = self.MHD.shotPath + 'openFoam/heatFoam/'+PFC.name

            #replace spaces with underscores to prevent pyfoam from reading args as
            # space delimited list
            partDir = partDir.replace(" ", "_")
            file = (partDir +
                    '/postProcessing/fieldMinMax1/{:f}'.format(self.OF.OFtMin).rstrip('0').rstrip('.')
                    +'/fieldMinMax.dat')
            data.append(self.OF.getMinMaxData(file))
            pfcNames.append(PFC.name)

        fig = pgp.plotlyOpenFOAMplot(data,pfcNames)

        #save interactive plotly plot in shotPath/plotly/OFminmax.html
        plotlyDir = self.MHD.shotPath + 'plotly'
        tools.makeDir(plotlyDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        plotPath = plotlyDir + '/OFminmax.html'
        fig.write_html(plotPath)

        #set tree permissions
        tools.recursivePermissions(plotlyDir, self.UID, self.GID, self.chmod)

        return fig

    def getHFdistPlots(self):
        """
        returns plotly figure with qDiv PFC surface distributions
        """
        heatFluxes = []
        labels = []
        for PFC in self.PFCs:
            heatFluxes.append(PFC.qDiv)
            labels.append(PFC.name)

        fig = pgp.plotlyqDivPlot(heatFluxes, labels, logPlot=True)

        #save interactive plotly plot in shotPath/plotly/HFdist.html
        plotlyDir = self.MHD.shotPath + 'plotly'
        tools.makeDir(plotlyDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        plotPath = plotlyDir + '/HFdist.html'
        fig.write_html(plotPath)

        #set tree permissions
        tools.recursivePermissions(plotlyDir, self.UID, self.GID, self.chmod)

        return fig

    def TprobeOF(self,x,y,z):
        """
        run temperature probe OF function: postProcess -func "probes"
        returns a dash figure for use in GUI
        """
        print("Solving for Temperature Probe")
        log.info("Solving for Temperature Probe")
        tData = []
        Tdata = []
        names = []
        for PFC in self.PFCs:
            #PFC boundary with 1mm buffer
            xMin = (PFC.centers[:,0].min() - 0.001)*1000.0
            xMax = (PFC.centers[:,0].max() + 0.001)*1000.0
            yMin = (PFC.centers[:,1].min() - 0.001)*1000.0
            yMax = (PFC.centers[:,1].max() + 0.001)*1000.0
            zMin = (PFC.centers[:,2].min() - 0.001)*1000.0
            zMax = (PFC.centers[:,2].max() + 0.001)*1000.0

            case1 = (x<xMax) and (x>xMin)
            case2 = (y<yMax) and (y>yMin)
            case3 = (z<zMax) and (z>zMin)

            #make sure the Tprobe is inside this tile
            if case1 and case2 and case3:
                print("Found Tprobe in PFC: "+PFC.name)
                log.info("Found Tprobe in PFC: "+PFC.name)
                partDir = self.OF.caseDir+'/'+PFC.name+'/'
                #replace spaces with underscores to prevent pyfoam from reading args as
                # space delimited list
                partDir = partDir.replace(" ", "_")
                self.OF.runTprobe(x,y,z,partDir)
                file = (partDir +
                    'postProcessing/probes/{:f}'.format(self.OF.OFtMin).rstrip('0').rstrip('.')
                    +'/T')
                data = np.genfromtxt(file,comments="#", autostrip=True)
                tData.append(data[:,0])
                Tdata.append(data[:,1])
                names.append(PFC.name)

            #if Tprobe not in this tile, dont make a figure
            else:
                print("Tprobe outside of PFC: "+PFC.name)
                log.info("Tprobe outside of PFC: "+PFC.name)
                pass

        fig = pgp.plotlyTprobes(tData,Tdata,names)

        #save interactive plotly plot in shotPath/plotly/Tprobes.html
        plotlyDir = self.MHD.shotPath + 'plotly'
        tools.makeDir(plotlyDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        plotPath = plotlyDir + '/Tprobes.html'
        fig.write_html(plotPath)

        #set tree permissions
        tools.recursivePermissions(plotlyDir, self.UID, self.GID, self.chmod)

        return fig

    def gyroPhasePlot(self):
        """
        return a gyrophase figure with gyro orbit phase angles
        """
        #setup gyroPhase angle
        self.GYRO.uniformGyroPhaseAngle()

        fig = pgp.plotlyGyroPhasePlot(np.degrees(self.GYRO.gyroPhases))
        return fig

    def vPhasePlot(self):
        """
        return a vPhase figure with velocity phase angles
        """
        #setup velocity phase angles
        self.GYRO.uniformVelPhaseAngle()
        fig = pgp.plotlyVPhasePlot(np.degrees(self.GYRO.vPhases))
        return fig

    def vSlicePlot(self):
        """
        return a vSlice figure with vSlices
        """
        #setup velocities and velocity phase angles
        self.GYRO.setupVelocities(1)
        fig = pgp.plotlyVSlicePlot(self.GYRO.mass_eV,
                                   self.GYRO.c,
                                   self.GYRO.T0[0],
                                   self.GYRO.vSlices[0,:],
                                   self.GYRO.vBounds[0,:],
                                   self.GYRO.vScan[0])
        return fig

    def cdfSlicePlot(self):
        """
        return a vSlice figure with vSlices
        """
        #setup velocities and velocity phase angles
        self.GYRO.setupVelocities(1)
        fig = pgp.plotlyallSlicePlot(self.GYRO.mass_eV,
                                     self.GYRO.c,
                                     self.GYRO.T0[0],
                                     self.GYRO.vSlices[0,:],
                                     self.GYRO.vBounds[0,:],
                                     self.GYRO.vScan[0],
                                     self.GYRO.N_vSlice,
                                     )
        return fig



    def writeInputTable(self, data):
        path = self.MHD.shotPath + '/inputData.csv'
        with open(path, 'w') as f:
            for line in data:
                p = line['Parameter']
                if "\u03bb" in p:
                    p = p.replace("\u03bb","lq")
                v = line['Value']
                f.write("{}:\t{}\n".format(p,v))
        return

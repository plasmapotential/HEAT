#GUIclass.py
#Description:   GUI Module
#Engineer:      T Looby
#Date:          20200115
"""
This is the launch point for the HEAT code, when accessing from the HTML GUI.
It calls HEAT functions and generates heat fluxes based upon user input.  These
functions are usually called from the dashGUI.py script, which is the FLASK
binding to html.
"""

import CADClass
import MHDClass
import toolsClass
import heatfluxClass
import openFOAMclass
import pfcClass
import gyroClass
import time
import numpy as np
import logging
import os
import pandas as pd
import shutil
import errno
import copy
import EFIT.equilParams_class as EP
import GUIscripts.plotlyGUIplots as pgp
import trimesh

log = logging.getLogger(__name__)
tools = toolsClass.tools()

def create_app(GUIobj):
    """
    Creates flask app with GUI object so that we can refer to it later.
    """
    from flask import Flask
    app = Flask(__name__)
    app.config['GUI'] = GUIobj
    return app

def create_DASH_app(GUIobj):
    import dash
    app = dash.Dash(__name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}])
    app.config['GUI'] = GUIobj
    return app

class GUIobj():
    def __init__(self, logFile, rootDir, dataPath, OFbashrc):
        #where HEAT log is written
        self.logFile = logFile
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
        try:
            os.mkdir(tempDir)
            print("Directory " , tempDir ,  " Created ")
        except FileExistsError:
            try: shutil.rmtree(tempDir)
            except OSError as e:
                print ("Error: %s - %s." % (e.filename, e.strerror))
                sys.exit()

            os.mkdir(tempDir)
            print("Directory " , tempDir ,  " Created ")
        return

    def machineSelect(self, MachFlag):
        """
        Select a machine and set the necessary paths
        """
        self.MachFlag = MachFlag
        self.MHD.MachFlag = MachFlag
        self.setInitialFiles()
        self.setHiddenInputs()
        return

    def initializeEveryone(self):
        """
        Create objects that we can reference later on
        """
        self.MHD = MHDClass.MHD(self.rootDir, self.dataPath)
        self.CAD = CADClass.CAD(self.rootDir, self.dataPath)
        self.HF = heatfluxClass.heatFlux(self.rootDir, self.dataPath)
        self.OF = openFOAMclass.OpenFOAM(self.rootDir, self.dataPath)
        self.GYRO = gyroClass.GYRO(self.rootDir, self.dataPath)

        #set up class variables for each object
        self.MHD.allowed_class_vars()
        self.CAD.allowed_class_vars()
        self.HF.allowed_class_vars()
        self.OF.allowed_class_vars()
        self.GYRO.allowed_class_vars()
        return


    def setInitialFiles(self):
        """
        sets files back to default settings
        infile is path to input file with all HEAT parameters
        PartsFile is path to file with parts we will calculate HF on
        IntersectFile is path to file with parts we will check for intersections on
        """
        if self.MachFlag == 'nstx':
            print('Loading NSTX-U Input Filestream')
            log.info('Loading NSTX-U Input Filestream')
            self.infile = self.rootDir + '/inputs/NSTXU/NSTXU_input.csv'
            self.pfcFile = self.rootDir + '/inputs/NSTXU/NSTXUpfcs.csv'
            self.OF.meshDir = self.dataPath + '/NSTX/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/NSTX/STLs/'
            self.CAD.STPpath = self.dataPath + '/NSTX/STPs/'

        elif self.MachFlag == 'st40':
            print('Loading ST40 Input Filestream')
            log.info('Loading ST40 Input Filestream')
            self.infile = self.rootDir + '/inputs/ST40/ST40_input.csv'
            self.pfcFile = self.rootDir + '/inputs/ST40/ST40pfcs.csv'
            self.OF.meshDir = self.dataPath + '/ST40/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/ST40/STLs/'
            self.CAD.STPpath = self.dataPath + '/ST40/STPs/'

        elif self.MachFlag == 'd3d':
            print('Loading DIII-D Input Filestream')
            log.info('Loading DIII-D Input Filestream')
            self.infile = self.rootDir + '/inputs/D3D/D3D_input.csv'
            self.pfcFile = self.rootDir + '/inputs/D3D/D3Dpfcs.csv'
            self.OF.meshDir = self.dataPath + '/D3D/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/D3D/STLs/'
            self.CAD.STPpath = self.dataPath + '/D3D/STPs/'

        elif self.MachFlag == 'step':
            print('Loading STEP Input Filestream')
            log.info('Loading STEP Input Filestream')
            self.infile = self.rootDir + '/inputs/STEP/STEP_input.csv'
            self.pfcFile = self.rootDir + '/inputs/STEP/STEPpfcs.csv'
            self.OF.meshDir = self.dataPath + '/STEP/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/STEP/STLs/'
            self.CAD.STPpath = self.dataPath + '/STEP/STPs/'

        elif self.MachFlag == 'sparc':
            print('Loading SPARC Input Filestream')
            log.info('Loading SPARC Input Filestream')
            self.infile = self.rootDir + '/inputs/SPARC/SPARC_input.csv'
            self.pfcFile = self.rootDir + '/inputs/SPARC/SPARCpfcs.csv'
            self.OF.meshDir = self.dataPath + '/SPARC/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/SPARC/STLs/'
            self.CAD.STPpath = self.dataPath + '/SPARC/STPs/'

        else:
            print("INVALID MACHINE SELECTION!  Defaulting to NSTX-U!")
            log.info("INVALID MACHINE SELECTION!  Defaulting to NSTX-U!")
            self.infile = self.rootDir + '/inputs/NSTXU/NSTXU_input.csv'
            self.pfcFile = self.rootDir + '/inputs/NSTXU/NSTXUpfcs.csv'
            self.OF.meshDir = self.dataPath + '/NSTX/3Dmeshes'
            self.CAD.STLpath = self.dataPath + '/NSTX/STLs/'
            self.CAD.STPpath = self.dataPath + '/NSTX/STPs/'

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
        self.MHD.Rmax = 5.0
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

        else:
            print("INVALID MACHINE SELECTION!  Defaulting to NSTX-U!")
            log.info("INVALID MACHINE SELECTION!  Defaulting to NSTX-U!")

        return


    def getMHDInputs(self,shot=None,tmin=None,tmax=None,nTrace=None,
                     gFileList=None,gFileData=None,plasma3Dmask=None,
                     ):
        """
        Get the mhd inputs from the gui or input file
        """
        tools.vars2None(self.MHD)
        tools.read_input_file(self.MHD, infile=self.infile)
        self.MHD.MachFlag=self.MachFlag #override input file with machine user selected
        self.MHD.setTypes()

        if shot is not None:
            self.MHD.shot = shot
        if tmin is not None:
            self.MHD.tmin = tmin
        if tmax is not None:
            self.MHD.tmax = tmax
        if nTrace is not None:
            self.MHD.nTrace = nTrace
        self.MHD.gFileList = gFileList
        if gFileList is not None:
            self.MHD.writeGfileData(gFileList, gFileData)

        if plasma3Dmask is not None:
            self.MHD.plasma3Dmask = plasma3Dmask

        self.MHD.tree = 'EFIT02'

        if self.dataPath[-1]!='/':
            self.MHD.shotPath = self.dataPath + '/' + self.MHD.MachFlag +"_{:06d}".format(self.MHD.shot)
        else:
            self.MHD.shotPath = self.dataPath + self.MHD.MachFlag +"_{:06d}".format(self.MHD.shot)

        self.MHD.get_mhd_inputs('nstx',self.MHD.gFileList)

#        self.t = self.MHD.timesteps[0]
        self.MHD.makeEFITobjects()
        self.NCPUs = 4
        self.MHD.psiSepLimiter = None

        self.MHD.setTypes()

        print('psiSep0 = {:f}'.format(self.MHD.ep[0].g['psiSep']))
        print('psiAxis0 = {:f}'.format(self.MHD.ep[0].g['psiAxis']))
        print('Nlcfs0: {:f}'.format(self.MHD.ep[0].g['Nlcfs']))
        print('length Rlcfs0: {:f}'.format(len(self.MHD.ep[0].g['lcfs'][:,0])))
        log.info('psiSep0 = {:f}'.format(self.MHD.ep[0].g['psiSep']))
        log.info('psiAxis0 = {:f}'.format(self.MHD.ep[0].g['psiAxis']))
        log.info('Nlcfs0: {:f}'.format(self.MHD.ep[0].g['Nlcfs']))
        if self.MHD.plasma3Dmask==1:
            print('Solving for 3D plasmas with MAFOT')
            log.info('Solving for 3D plasmas with MAFOT')
        else:
            print('Solving for 2D plasmas with EFIT (no MAFOT)')
            log.info('Solving for 2D plasmas with EFIT (no MAFOT)')
        return

    def gfileClean(self, psiRZMult,psiSepMult,psiAxisMult,FpolMult,
                   psiRZAdd,psiSepAdd,psiAxisAdd,FpolAdd, t):
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

        idx = np.where(t==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[idx]
        ep.g['psiRZ'] *= psiRZMult
        ep.g['psiSep'] *= psiSepMult
        ep.g['psiAxis'] *= psiAxisMult
        ep.g['Fpol'] *= FpolMult

        ep.g['psiRZ'] += psiRZAdd
        ep.g['psiSep'] += psiSepAdd
        ep.g['psiAxis'] += psiAxisAdd
        ep.g['Fpol'] += FpolAdd

        psi = ep.g['psiRZ']
        psiSep = ep.g['psiSep']
        psiAxis = ep.g['psiAxis']
        ep.g['psiRZn'] = (psi - psiAxis) / (psiSep - psiAxis)
        return


    def findPsiSepfromEQ(self,t, rNew=None):
        """
        finds psiSep by stepping to/from core and calculating
        minimum psiN along Z-plane at each R location.  Increments in um
        """
        tIdx = np.where(t==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[tIdx]
        gfile = self.MHD.shotPath + '/' + '{:06d}/'.format(t) + 'g{:6d}.{:05d}'.format(self.MHD.shot, t)
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

    def findPsiSepfromPFCs(self, t, rNew=None):
        """
        finds psiSep for limiters by incrementally increasing R_{psiSep, IMP},
        at a specific time.  Then rewrite all gfiles in MHD.timesteps with this
        new LCFS.

        Both MHD and PFC objects must be defined before running this function
        """
        tIdx = np.where(t==self.MHD.timesteps)[0][0]
        gfile = self.MHD.shotPath + '{:06d}/'.format(t) + 'g{:06d}.{:05d}'.format(self.MHD.shot, t)
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
        gfile = self.MHD.shotPath + '{:06d}/'.format(t) + 'g{:06d}.{:05d}'.format(self.MHD.shot, t)
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
        gFileName = self.tmpDir + 'g{:06d}.{:05d}'.format(self.MHD.shot,t)
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
        #change filename to d3d convention
        self.MHD.tmax = int(max(timesteps))
        self.MHD.tmin = int(min(timesteps))
        shot = self.MHD.shot
        newGfiles = []
        for i,f in enumerate(gfiles):
            newGfiles.append('g{:06d}.{:05d}'.format(shot,timesteps[i]))
            old = self.tmpDir + gfiles[i]
            new = self.tmpDir + newGfiles[i]
            try:
                shutil.copyfile(old, new)
            except:
                print("Could not copy timestep {:d}.  Skipping.".format(timesteps[i]))
                log.info("Could not copy timestep {:d}.  Skipping.".format(timesteps[i]))


        #rebuild eq objects
        self.MHD.get_mhd_inputs(self.MachFlag,newGfiles)
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

    def getCADResInputs(self,ROIGridRes=None,gridRes=None):
        """
        Loads CAD inputs
        """
        import numpy as np
        tools.initializeInput(self.CAD, infile=self.infile)
        self.CAD.rootDir = self.rootDir #set HEAT rootDir from HEATgui.py
        if ROIGridRes is not None:
            self.CAD.ROIGridRes = ROIGridRes
        if gridRes is not None:
            #check if intersection grid resolution string is a number,
            #if not use standard mesh algorithms
            if tools.is_number(gridRes):
                self.CAD.gridRes = gridRes
            else:
                self.CAD.gridRes = "standard"
        return

    def getCAD(self,STPfile=None,STPdata=None, ts=None):
        """
        Loads CAD file
        """
        import numpy as np
        if hasattr(self.CAD, 'gridRes'):
            pass
        else:
            tools.initializeInput(self.CAD, infile=self.infile)
        self.CAD.rootDir = self.rootDir #set HEAT rootDir from HEATgui.py

        if STPfile is not None:
            #make STP path if it doesnt exist
            try:
                os.makedirs(self.CAD.STPpath)
            except:
                print("Did not create STPpath.  It probably exists already (HEAT already ran on this machine).")
            newSTPpath = self.CAD.STPpath + STPfile
            #check to see if this STP file exists and write data to the file
            if os.path.isfile(newSTPpath) == False:
                print("New STP file.  Writing")
                with open(newSTPpath, 'wb') as f:
                    f.write(STPdata)
                atime = os.stat(newSTPpath).st_atime
                os.utime(newSTPpath, (atime, ts))
                self.CAD.overWriteMask = True #we need to also overwrite meshes
            else:
                #if file was modified, overwrite
                if ts != os.path.getmtime(newSTPpath):
                    print("File was modified since last HEAT upload.  Overwriting...")
                    with open(newSTPpath, 'wb') as f:
                        f.write(STPdata)
                    atime = os.stat(newSTPpath).st_atime
                    os.utime(newSTPpath, (atime, ts))
                    self.CAD.overWriteMask = True #we need to also overwrite meshes
                else:
                    self.CAD.overWriteMask = False #meshes are already up to date
                    print("STP file is already in the HEAT database.  Not overwriting...")
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
        self.timestepMap['DivCode'] = self.timestepMap['DivCode'].str.strip()
        return

    def readPFCfile(self, infile):
        """
        Reads PFC names, timesteps, intersect names, and mapdirection from file
        """
        if infile == None:
            print("No timesteps input file.  Please provide input file")
            log.info("No timesteps input file.  Please provide input file")
            sys.exit()

        self.timestepMap = pd.read_csv(infile, comment='#')
        self.timestepMap = self.timestepMap.rename(columns=lambda x: x.strip())
        self.timestepMap['PFCname'] = self.timestepMap['PFCname'].str.strip()
        self.timestepMap['intersectName'] = self.timestepMap['intersectName'].str.strip()
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
        self.CAD.writeMesh2file(self.CAD.ROImeshes, self.CAD.ROI, path=self.CAD.STLpath)
        print("Calculating HF on these tiles:")
        log.info("Calculating HF on these tiles:")
        print(self.CAD.ROIList)
        log.info(self.CAD.ROIList)

        self.PFCs = []
        #initialize PFC objects for each ROI part
        for i,row in enumerate(self.timestepMap.values):
            PFC = pfcClass.PFC(row, self.rootDir, self.dataPath)
            PFC.makePFC(self.MHD, self.CAD, i, self.timestepMap['powerDirection'][i], clobberFlag=True)
            self.PFCs.append(PFC)

        #Find potential intersections by file as they correspond to ROI PFCs,
        # then mesh 'em using FreeCAD Standard mesh algorithm
        self.CAD.getIntersectsFromFile(self.timestepMap)
        self.CAD.getIntersectMeshes(resolution=self.CAD.gridRes)
        self.CAD.writeMesh2file(self.CAD.intersectMeshes,
                                self.CAD.intersectList,
                                path=self.CAD.STLpath,
                                resolution=self.CAD.gridRes
                                )
        print("Potential intersects on these tiles:")
        log.info("Potential intersects on these tiles:")
        print(self.CAD.intersectList)
        log.info(self.CAD.intersectList)

        #Build HEAT file tree
        tools.buildDirectories(self.CAD.ROIList,
                                         self.MHD.timesteps,
                                         self.MHD.shotPath,
                                         clobberFlag=True
                                         )

        #assign tag if PFC is run in multiple directions (multiple lines in XXXpfc.csv)
        for PFC in self.PFCs:
            bool = np.where(self.timestepMap['PFCname'] == PFC.name)[0]
            if len(bool) > 1:
                if PFC.mapDirection > 0:
                    PFC.tag = 'forward'
                else:
                    PFC.tag = 'reverse'
            else:
                PFC.tag = None


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

    def getHFInputs(self,lqEich,S,P,qBG,lqPN,lqPF,lqCN,lqCF,
                        fracPN, fracPF, fracCN, fracCF,
                        fracUI,fracUO,fracLI,fracLO,
                        mode, LRmask, LRpower,
                        lqCNmode,lqCFmode,SMode,fG):
        """
        get heat flux inputs from gui or input file
        """
        self.initializeHF()
        self.HF.hfMode = mode
        self.HF.lqEich = float(lqEich)
        self.HF.lqCN = float(lqCN)
        self.HF.lqCF = float(lqCF)
        self.HF.lqPN = float(lqPN)
        self.HF.lqPF = float(lqPF)
        self.HF.S = float(S)
        self.HF.Psol = float(P)
        self.HF.qBG = float(qBG)
        self.HF.mode = mode
        self.HF.fracCN = float(fracCN)
        self.HF.fracCF = float(fracCF)
        self.HF.fracPN = float(fracPN)
        self.HF.fracPF = float(fracPF)
        self.HF.fracUI = float(fracUI)
        self.HF.fracUO = float(fracUO)
        self.HF.fracLI = float(fracLI)
        self.HF.fracLO = float(fracLO)
        if 'yes' in LRmask:
            self.HF.LRmask = True
            self.HF.LRpower = float(LRpower)
        else:
            self.HF.LRmask = False

        self.HF.lqCNmode = lqCNmode
        self.HF.lqCFmode = lqCFmode
        self.HF.SMode = SMode
        self.HF.fG = float(fG)

        print("Mode = "+mode)
        log.info("Mode = "+mode)
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
        print("Long range intersection checking: "+LRmask)

        if mode=='eich':
            print("lqEich = {:f}".format(self.HF.lqEich))
            print("S = {:f}".format(self.HF.S))
            print("qBG = {:f}".format(self.HF.qBG))
            print("fG = {:f}".format(self.HF.fG))
            log.info("lqEich = {:f}".format(self.HF.lqEich))
            log.info("S = {:f}".format(self.HF.S))
            log.info("qBG = {:f}".format(self.HF.qBG))
            log.info("fG = {:f}".format(self.HF.fG))
        elif mode=='limiter':
            print("lqCN = {:f}".format(self.HF.lqCN))
            print("lqCF = {:f}".format(self.HF.lqCF))
            print("fracCN = {:f}".format(self.HF.fracCN))
            print("fracCF = {:f}".format(self.HF.fracCF))
            log.info("lqCN = {:f}".format(self.HF.lqCN))
            log.info("lqCF = {:f}".format(self.HF.lqCF))
            log.info("fracCN = {:f}".format(self.HF.fracCN))
            log.info("fracCF = {:f}".format(self.HF.fracCF))
        elif mode=='multiExp':
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
        if hasattr(self.MHD, 'ep'):
            self.HF.getHFtableData(self.MHD.ep[0])
        return

    def getGyroInputs(self,N_gyroSteps,N_gyroPhase,gyroDeg,ionMassAMU,vMode,gyroT_eV,
                      N_vPhase, N_vSlice, ionFrac):
        """
        Sets up the gyro module
        """
        self.GYRO.N_gyroSteps = int(N_gyroSteps)
        self.GYRO.gyroDeg = int(gyroDeg)
        self.GYRO.gyroT_eV = float(gyroT_eV)
        self.GYRO.N_vSlice = int(N_vSlice)
        self.GYRO.N_vPhase = int(N_vPhase)
        self.GYRO.N_gyroPhase = int(N_gyroPhase)
        self.GYRO.N_MC = self.GYRO.N_gyroPhase*self.GYRO.N_vSlice*self.GYRO.N_vPhase
        self.GYRO.ionMassAMU = float(ionMassAMU)
        self.GYRO.vMode = vMode
        self.GYRO.ionFrac = float(ionFrac)
        #set up GYRO object
        self.GYRO.setupConstants(self.GYRO.ionMassAMU)
        print('Loaded Gyro Orbit Settings')
        print('# Steps per helix period = {:f}'.format(float(N_gyroSteps)))
        print('Gyro tracing distance [degrees] = {:f}'.format(float(gyroDeg)))
        print('Plasma Temperature Mode = ' + vMode)
        print('Number of Monte Carlo runs per point = {:f}'.format(float(self.GYRO.N_MC)))
        return

    def bfieldAtSurface(self, PFC):
        """
        Calculate the B field at tile surface

        """
        ctrs = PFC.centers
        R,Z,phi = tools.xyz2cyl(ctrs[:,0],ctrs[:,1],ctrs[:,2])
        PFC.Bxyz = self.MHD.Bfield_pointcloud(PFC.ep, R, Z, phi, PFC.powerDirection)
        self.MHD.write_B_pointcloud(ctrs,PFC.Bxyz,PFC.controlfilePath)
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
        return

    def BtraceMultiple(self,data,t):
        """
        Run a MAFOT structure trace from multiple points defined in the gui
        """
        data = pd.DataFrame.from_dict(data)[list (data[0].keys())]
        data = data.rename(columns=lambda x: x.strip())
        data = data.astype({"x[mm]": float, "y[mm]": float, "z[mm]": float, "traceDirection": int, "Length[deg]":float})

        t = int(t)
        traceDirection=data['traceDirection']
        x = data['x[mm]'] / 1000.0
        y = data['y[mm]'] / 1000.0
        z = data['z[mm]'] / 1000.0

        xyz = np.array([x,y,z]).T
        controlfile = '_structCTL.dat'
        dphi = 1.0

        if self.MHD.shotPath[-1]=='/':
            gridfile = self.MHD.shotPath + '{:06d}/struct_grid.dat'.format(t)
            controlfilePath =  self.MHD.shotPath + '{:06d}/'.format(t)
        else:
            gridfile = self.MHD.shotPath + '/' + '{:06d}/struct_grid.dat'.format(t)
            controlfilePath =  self.MHD.shotPath + '/' + '{:06d}/'.format(t)
        structOutfile = controlfilePath + 'struct.dat'
        for i in range(len(xyz)):
            self.MHD.ittStruct = data['Length[deg]'][i]
            self.MHD.writeControlFile(controlfile, t, data['traceDirection'][i], mode='struct')
            self.MHD.writeMAFOTpointfile(xyz[i,:],gridfile)
            self.MHD.getFieldpath(dphi, gridfile, controlfilePath, controlfile, paraview_mask=True, tag='pt{:03d}'.format(i))
            os.remove(structOutfile)
        return


    def Btrace(self,x,y,z,t,mapDirection, traceDeg):
        """
        Run a MAFOT structure trace from a point defined in the gui
        """
#        idx = np.where(t==self.MHD.timesteps)[0][0]
#        ep = self.MHD.ep[idx]
        t = int(t)
        mapDirection=int(mapDirection)
        x = float(x)/1000.0
        y = float(y)/1000.0
        z = float(z)/1000.0

        xyz = np.array([x,y,z])
        controlfile = '_structCTL.dat'
        dphi = 1.0

        self.MHD.ittStruct = float(traceDeg)
        if self.MHD.shotPath[-1]=='/':
            gridfile = self.MHD.shotPath + '{:06d}/struct_grid.dat'.format(t)
            controlfilePath =  self.MHD.shotPath + '{:06d}/'.format(t)
        else:
            gridfile = self.MHD.shotPath + '/' + '{:06d}/struct_grid.dat'.format(t)
            controlfilePath =  self.MHD.shotPath + '/' + '{:06d}/'.format(t)
        self.MHD.writeControlFile(controlfile, t, mapDirection, mode='struct')
        self.MHD.writeMAFOTpointfile(xyz,gridfile)
        self.MHD.getFieldpath(dphi, gridfile, controlfilePath, controlfile, paraview_mask=True)
        return

    def gyroTrace(self,x,y,z,t,gyroPhase,gyroDeg,N_gyroSteps,gyroDir,gyroT_eV):
        """
        performs a gyro orbit trace from GUI defined location

        (x,y,z) are locations where we launch trace from
        gyroPhase is initial phase angle of orbit in degrees
        gyroDeg is the number of degrees we will trace for
        N_gyroSteps is the number of discrete lines we approximate helical path by
        """
        print("\n========Gyro Trace Initialized========")
        #get bField trace from this point
        self.Btrace(x,y,z,t,gyroDir,gyroDeg)
        #read bField trace csv output
        if self.MHD.shotPath[-1]=='/':
            structOutfile = self.MHD.shotPath + '{:06d}/struct.dat'.format(t)
            controlfilePath =  self.MHD.shotPath + '{:06d}/'.format(t)
        else:
            structOutfile = self.MHD.shotPath + '/' + '{:06d}/struct.dat'.format(t)
            controlfilePath =  self.MHD.shotPath + '/' + '{:06d}/'.format(t)
        BtraceXYZ = tools.readStructOutput(structOutfile) #[m]
        #Setup gyro orbit trace constants and velocities
        self.GYRO.setupConstants()
        v = self.GYRO.temp2thermalVelocity(float(gyroT_eV))
        vPerp = v*np.cos(np.pi/4)
        vParallel = v*np.sin(np.pi/4)
        # Evaluate B
        R,Z,phi = tools.xyz2cyl(float(x)/1000.0,float(y)/1000.0,float(z)/1000.0)#mm => m
        tIdx = np.where(float(t)==self.MHD.timesteps)[0][0]
        ep = self.MHD.ep[tIdx]
        Bt = ep.BtFunc.ev(R,Z)
        BR = ep.BRFunc.ev(R,Z)
        BZ = ep.BZFunc.ev(R,Z)
        B = np.sqrt(BR**2 + Bt**2 + BZ**2)
        #Calculate frequencies and gyro radius
        self.GYRO.setupFreqs(B)
        self.GYRO.setupRadius(vPerp)
        #trace helix and save to CSV and VTK formats
        self.GYRO.singleGyroTrace(vPerp,vParallel,float(gyroPhase),float(N_gyroSteps),BtraceXYZ,controlfilePath,
                                  self.GYRO.TGyro[0],self.GYRO.rGyro[0],self.GYRO.omegaGyro[0],)
        print("Temp of {:f} eV yields thermal velocity of {:f} m/s".format(float(gyroT_eV), float(vPerp)))
        print("B magntidude = {:f}".format(B))
        return




    def NormPC(self, PFC):
        """
        create a normal vector point cloud for mesh centers on tile surface
        """
        self.CAD.write_normal_pointcloud(PFC.centers,PFC.norms,PFC.controlfilePath)
        return

    def shadowPC(self, PFC):
        """
        create a pointcloud for mesh center locations where 1=shadowed, 0=not shadowed
        """
        PFC.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,PFC.controlfilePath)
        return


    def backfacePC(self, PFC):
        """
        create a pointcloud for mesh center locations where 1=shadowed, 0=not shadowed
        using the backface culling algorithm
        """

        PFC.backfaceMask = PFC.backfaceCulling(PFC.centers,PFC.norms,self.MHD,PFC.ep,PFC.powerDirection)
        PFC.write_backface_pointcloud(PFC.centers,PFC.backfaceMask,PFC.controlfilePath)
        return

    def bdotnPC(self, PFC):
        """
        makes a cos(alpha) point cloud:b_hat dot n_hat
        where b_hat is normalized magnetic field vector and n_hat is mesh surface
        normal vector
        """
        self.HF.HFincidentAngle(PFC,self.MHD)
        PFC.write_bdotn_pointcloud(PFC.centers, PFC.bdotn, PFC.controlfilePath)


    def initializeHF(self):
        """
        Initialize heat flux variables
        """
        print("-"*70)
        print("HEAT FLUX MODULE INITIALIZED")
        log.info("HEAT FLUX MODULE INITIALIZED")
        #Initialize HF Object
        tools.initializeInput(self.HF, infile=self.infile)
        return

    def runHEAT(self, runList):
        """
        Run a HEAT calculation.  This is called from gui by user.
        Creates point clouds at each mesh center of PFC objects in CAD ROI

        Steps forward in time, and then solves everything for each PFC at that
        timestep

        runList options are:
        Bpc             Bfield point cloud
        psiPC           Normalized psi point cloud
        backfacePC      backface culling point cloud
        NormPC          Normal Field point cloud
        bdotnPC         bdotn poitn cloud
        HFpc            heat flux point cloud
        GyroPC          gyro orbit heat flux point cloud
        """
        print('\n')
        print("-"*70)
        print("HEAT RUN INITIALIZED")
        log.info("HEAT RUN INITIALIZED")
        t0 = time.time()
        allowedOptions = ['HFpc', 'backfacePC', 'bdotnPC', 'Bpc', 'psiPC', 'NormPC', 'GyroPC']
        #make sure that something in runList can be run in this function, else return
        if len([i for i in runList if i in allowedOptions]) < 1:
            print("No HEAT point cloud option to run")
            return
        #set up variables for power balance calculation
        powerTesselate = np.zeros((len(self.MHD.timesteps)))
        powerTrue = np.zeros((len(self.MHD.timesteps)))
        powerByTile = np.zeros((len(self.PFCs)))
        divCodes = []
        #set up electron frac if not in gyro mode
        if 'GyroPC' not in runList:
            self.HF.elecFrac = 1.0
        else:
            self.HF.elecFrac = 1.0 - self.GYRO.ionFrac

        # Time Loop 1: HF, bdotn, B, psi, Norms
        for tIdx,t in enumerate(self.MHD.timesteps):
            print('\n')
            print("-"*80)
            log.info("-"*80)
            print("Timestep: {:d}".format(t))
            log.info("Timestep: {:d}\n".format(t))
            print("-"*80)
            log.info("-"*80)
            for PFC in self.PFCs:
                if t not in PFC.timesteps:
                    pass
                else:
                    #set up file directory structure
                    PFC.controlfile = '_lamCTL.dat'
                    PFC.controlfileStruct = '_struct_CTL.dat'
                    if self.MHD.shotPath[-1]=='/':
                        PFC.controlfilePath = self.MHD.shotPath + '{:06d}/'.format(t) + PFC.name + '/'
                        PFC.gridfile = self.MHD.shotPath + '{:06d}/'.format(t) + PFC.name + '/grid.dat'
                        PFC.gridfileStruct = self.MHD.shotPath + '{:06d}/'.format(t) + PFC.name + '/struct_grid.dat'
                    else:
                        PFC.controlfilePath = self.MHD.shotPath + '/' + '{:06d}/'.format(t) + PFC.name + '/'
                        PFC.gridfile = self.MHD.shotPath + '/' + '{:06d}/'.format(t) + PFC.name + '/grid.dat'
                        PFC.gridfileStruct = self.MHD.shotPath + '/' + '{:06d}/'.format(t) + PFC.name + '/struct_grid.dat'
                    PFC.outputFile = PFC.controlfilePath + 'lam.dat'
                    PFC.structOutfile = PFC.controlfilePath + 'struct.dat'
                    #set up time and equilibrium
                    PFC.t = t
                    PFC.ep = PFC.EPs[tIdx]
                    PFC.shadowed_mask = PFC.shadowMasks[tIdx]
                    print('\n')
                    print("*"*20)
                    print('PFC Name: '+ PFC.name)
                    log.info('PFC Name: '+ PFC.name)
                    if 'HFpc' in runList:
                        self.HF_PFC(PFC, PFC.tag)
                        PFC.shadowMasks[tIdx] = PFC.shadowed_mask
                        PFC.powerSumOptical[tIdx] = self.HF.power_sum_mesh(PFC, mode='optical')
                        print('Maximum optical heat load on tile: {:f}'.format(max(PFC.qDiv)))
                        print('Theoretical optical power to this divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac*self.HF.elecFrac))
                        print('Tessellated divertor power to this PFC = {:f}'.format(PFC.powerSumOptical[tIdx]))
                        log.info('Maximum heat load on tile: {:f}'.format(max(PFC.qDiv)))
                        log.info('Power to this Divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac))
                        log.info('Tessellated Total Power = {:f}'.format(PFC.powerSumOptical[tIdx]))
                        print("\nTime Elapsed: {:f}".format(time.time() - t0))
                        log.info("\nTime Elapsed: {:f}".format(time.time() - t0))
                        powerTesselate[tIdx] += PFC.powerSumOptical[tIdx]
                        #Add ground truth power for all the PFCs, but not if we
                        #already counted this divertor
                        if PFC.DivCode not in divCodes:
                            powerTrue[tIdx] += self.HF.Psol*PFC.powerFrac
                        divCodes.append(PFC.DivCode)
                    if 'Bpc' in runList:
                        self.bfieldAtSurface(PFC)
                    if 'psiPC' in runList:
                        self.psiPC(PFC)
                    if 'NormPC' in runList:
                        self.NormPC(PFC)
                    if 'backfacePC' in runList:
                        self.backfacePC(PFC)
                    if 'bdotnPC' in runList:
                        self.bdotnPC(PFC)

        if 'HFpc' in runList:
            print('Total Input Power = {:f}'.format(np.sum(powerTrue)))
            print('Theoretical power to this divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac))
            print('Optical Tessellated Total Power = {:f}'.format(np.sum(PFC.powerSumOptical)))
            log.info('Total Input Power = {:f}'.format(np.sum(powerTrue)))
            log.info('Theoretical power to this divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac))
            log.info('Optical Tessellated Total Power = {:f}'.format(np.sum(PFC.powerSumOptical)))

            totalPowPow = 0
            for PFC in self.PFCs:
                print("=== Last timestep's PFC arrays: Optical ===")
                tmpPow = self.HF.power_sum_mesh(PFC, scale2circ=False, verbose=False)
                totalPowPow += tmpPow
                print(PFC.name + ":\t{:.6f}".format(tmpPow))
                log.info(PFC.name + ":\t{:.6f}".format(tmpPow))
                print("PFC array sum: {:.6f}".format(totalPowPow))
                log.info("PFC array sum: {:.6f}".format(totalPowPow))

        if 'GyroPC' in runList:
            print("\n===+++ GYRO ORBIT CALCULATION +++===")
            log.info("\n===+++ GYRO ORBIT CALCULATION +++===")
            PFC.powerSumGyro = np.zeros((len(self.MHD.timesteps)))
            tGyro = time.time()
            self.getGyroMeshes()

            # Time Loop: gyro orbits runs after optical HF calculated
            for tIdx,t in enumerate(self.MHD.timesteps):
                #path for this timestep
                if self.MHD.shotPath[-1]=='/':
                    PFC.controlfilePath = self.MHD.shotPath + '{:06d}/'.format(t) + PFC.name + '/'
                else:
                    PFC.controlfilePath = self.MHD.shotPath + '/' + '{:06d}/'.format(t) + PFC.name + '/'
                #generate index maps
                self.prepareGyroMaps(tIdx)
                for PFC in self.PFCs:
                    if t not in PFC.timesteps:
                        pass
                    else:
                        print("\n------PFC: "+PFC.name+"------")
                        self.gyroOrbitIntersects(PFC)
                self.gyroOrbitHF()
                Psum = self.HF.power_sum_mesh(PFC, mode='gyro', scale2circ=True)
                PFC.powerSumGyro[tIdx] += Psum
                powerTesselate[tIdx] += Psum
                print(PFC.name + ' tessellated Gyro Power = {:f}'.format(Psum))

                #path for this timestep
                tPath = self.MHD.shotPath + '/' + '{:06d}/'.format(t)
                #write intersectRecord to CSV file
                self.intersectRecordCSV(tPath)

                print('Theoretical gyro power to this divertor: {:f}'.format(self.HF.Psol*PFC.powerFrac*self.GYRO.ionFrac))
                print('Gyro Tessellated Power = {:f}'.format(np.sum(PFC.powerSumGyro)))
                print('Escaped Gyro Power = {:f}'.format(self.GYRO.gyroNanPower))
                print("Gyro orbit calculation took: {:f} [s]\n".format(time.time() - tGyro))

                print("=== Last timestep's PFC arrays: Gyro ===")
                totalPowPow = 0
                for PFC in self.PFCs:
                    tmpPow = self.HF.power_sum_mesh(PFC, mode='gyro', scale2circ=False, verbose=False)
                    totalPowPow += tmpPow
                    print(PFC.name + ":\t{:.6f}".format(tmpPow))
                    log.info(PFC.name + ":\t{:.6f}".format(tmpPow))
                    print("PFC array sum: {:.6f}".format(totalPowPow))
                    log.info("PFC array sum: {:.6f}".format(totalPowPow))
                    print("PFC Max HF: {:.6f}".format(max(PFC.qGyro)))
                    log.info("PFC Max HF: {:.6f}".format(max(PFC.qGyro)))

                totalPowPow = 0
                for PFC in self.PFCs:
                    print("=== Last timestep's PFC arrays: Optical ===")
                    tmpPow = self.HF.power_sum_mesh(PFC, scale2circ=False, verbose=False)
                    totalPowPow += tmpPow
                    print(PFC.name + ":\t{:.6f}".format(tmpPow))
                    log.info(PFC.name + ":\t{:.6f}".format(tmpPow))
                    print("PFC array sum: {:.6f}".format(totalPowPow))
                    log.info("PFC array sum: {:.6f}".format(totalPowPow))


        #time/PFC loop for generating allSources heat fluxes
        if 'GyroPC' in runList or 'HFpc' in runList:
            for tIdx,t in enumerate(self.MHD.timesteps):
                #path for this timestep
                if self.MHD.shotPath[-1]=='/':
                    PFC.controlfilePath = self.MHD.shotPath + '{:06d}/'.format(t) + PFC.name + '/'
                else:
                    PFC.controlfilePath = self.MHD.shotPath + '/' + '{:06d}/'.format(t) + PFC.name + '/'
                #set up time and equilibrium
                PFC.t = t
                for PFC in self.PFCs:
                    if t not in PFC.timesteps:
                        pass
                    else:
                        print("Creating allSources CSV files")
                        R,Z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])
                        if 'GyroPC' in runList:
                            self.HF.write_heatflux_pointcloud(PFC.centers,PFC.qGyro+PFC.qDiv,PFC.controlfilePath,tag=PFC.tag,mode='all')
                        else:
                            self.HF.write_heatflux_pointcloud(PFC.centers,PFC.qDiv,PFC.controlfilePath,tag=PFC.tag,mode='all')

        # Time Loop: postprocessing
        for tIdx,t in enumerate(self.MHD.timesteps):
            #path for this timestep
            tPath = self.MHD.shotPath + '/' + '{:06d}/'.format(t)
            #merge multiple pointclouds into one single pointcloud for visualization
            self.combinePFCpointcloud(runList, tPath)
            #copy each timestep's composite point clouds to central location for
            #paraview postprocessing (movies)
            self.combineTimeSteps(runList, t)


        print("Total Time Elapsed: {:f}".format(time.time() - t0))
        log.info("Total Time Elapsed: {:f}".format(time.time() - t0))
        print("\nCompleted HEAT run\n")
        log.info("\nCompleted HEAT run\n")
        #make a sound when complete
#        os.system('spd-say -t female2 "HEAT run complete"')

        return

    def HF_PFC(self, PFC, tag=None):
        """
        meat and potatoes of the HF calculation.  Called in loop or by parallel
        processes for each PFC object.  Only calculates optical heat flux
        """
        #Check for intersections with MAFOT struct
        t0 = time.time()
        PFC.findShadows_structure(self.MHD, self.CAD)
        print("Intersection calculation took {:f} s".format(time.time() - t0))

        #Run MAFOT laminar for 3D plasmas
        if self.MHD.plasma3Dmask==True:
#            print('\n')
#            print("-"*70)
#            print("MAFOT LAMINAR MODULE INITIALIZED")
#            log.info("MAFOT LAMINAR MODULE INITIALIZED")
            CTLfile=PFC.controlfilePath + PFC.controlfile
            self.MHD.writeControlFile(CTLfile, PFC.t, PFC.mapDirection, mode='laminar')
            use = np.where(PFC.shadowed_mask != 1)[0]
            self.MHD.writeMAFOTpointfile(PFC.centers[use],PFC.gridfile)
            self.MHD.runMAFOTlaminar(PFC.gridfile,PFC.controlfilePath,PFC.controlfile,self.NCPUs)
            self.HF.readMAFOTLaminarOutput(PFC,PFC.outputFile)
            os.remove(PFC.outputFile)
        #get psi from gfile for 2D plasmas
        else:
            self.MHD.psi2DfromEQ(PFC)

        #Create Heat Flux Profile
        q = self.HF.getHFprofile(PFC)
        qDiv = self.HF.q_div(PFC, self.MHD, q) * self.HF.elecFrac

        #Points over threshold power are likely errors, so check them
        if self.HF.LRmask == True:
            distPhi = 370.0 #degrees
            qDiv = PFC.longRangeIntersectCheck(qDiv,
                                               self.HF.LRpower,
                                               distPhi,
                                               self.MHD,
                                               self.CAD)

        #Save data to class variable for future use
        PFC.q = q
        PFC.qDiv = qDiv

        #Create pointclouds for paraview
        R,Z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])
        PFC.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,PFC.controlfilePath,PFC.tag)
        self.HF.write_heatflux_pointcloud(PFC.centers,qDiv,PFC.controlfilePath,PFC.tag)
        PFC.write_bdotn_pointcloud(PFC.centers, PFC.bdotn, PFC.controlfilePath,PFC.tag)
        #structOutfile = MHD.shotPath + '/' + '{:06d}/struct.csv'.format(PFC.t)
        #HF.PointCloudfromStructOutput(structOutfile)

        return

    def gyroOrbitIntersects(self, PFC):
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
            #setup velocities and velocity phase angles
            self.GYRO.setupVelocities(PFC.N_gyroCenters)
            #setup gyroPhase angle
            self.GYRO.uniformGyroPhaseAngle()
            #setup frequencies
            self.GYRO.setupFreqs(PFC.Bmag[PFC.PFC_GYROmap,-1])
            self.HF.gyroHF(self.GYRO, PFC)

        for PFC in self.PFCs:
            #print("PFC NAME: "+PFC.name)
            #print(PFC.CADTGT_PFCmap)
            #assign gyro power to PFC
            PFC.Pgyro = self.GYRO.gyroPowMatrix[PFC.CADTGT_PFCmap]
            PFC.qGyro = PFC.Pgyro / PFC.areas
            #gyro shadowMask = 1 if a face is shadowed
            PFC.gyroShadowMask = self.GYRO.shadowMask[PFC.CADROI_PFCmap]
            #Create pointcloud for paraview
            R,Z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])
            self.HF.write_heatflux_pointcloud(PFC.centers,PFC.qGyro,PFC.controlfilePath,tag=PFC.tag,mode='gyro')
            self.HF.write_heatflux_pointcloud(PFC.centers,PFC.qGyro+PFC.qDiv,PFC.controlfilePath,tag=PFC.tag,mode='all')
            PFC.write_shadow_pointcloud(PFC.centers,PFC.gyroShadowMask,PFC.controlfilePath,PFC.tag,mode='gyro')

        return

    def getGyroMeshes(self):
        """
        set up gyro meshes, independent of timestep
        """
        totalMeshCounter = 0
        numTargetFaces = 0
        numROIFaces = 0
        targetPoints = []
        targetNorms = []
        self.GYRO.CADtargetNames = []
        self.GYRO.CADROINames = []
        self.GYRO.CADROIindexes = []
        self.GYRO.CADROICenters = []


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
                numTargetFaces += target.CountFacets
                #append target data
                for face in target.Facets:
                    self.GYRO.CADtargetNames.append(self.CAD.intersectList[i]) #do this for future HF reassignment
                    targetPoints.append(face.Points)
                    targetNorms.append(face.Normal)

        targetPoints = np.asarray(targetPoints)/1000.0 #scale to m
        targetNorms = np.asarray(targetNorms)
        self.GYRO.t1 = targetPoints[:,0,:] #point 1 of mesh triangle
        self.GYRO.t2 = targetPoints[:,1,:] #point 2 of mesh triangle
        self.GYRO.t3 = targetPoints[:,2,:] #point 3 of mesh triangle
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
        self.GYRO.gyroNanPower = 0.0
        #set up intersectRecord, which records index of intersection face
        #this 4D array has:
        #one dimension for gyroPhase angles,
        #one dimension for vPhases,
        #one dimension for vSlices,
        #and one dimesion for the PFC.centers points we are tracing
        #each element is the face that was intersected for that MC run
        self.GYRO.intersectRecord = np.ones((self.GYRO.N_gyroPhase,
                                            self.GYRO.N_vPhase,
                                            self.GYRO.N_vSlice,
                                            self.GYRO.N_HOT), dtype=int)*np.NaN

        self.GYRO.hdotn = np.ones((self.GYRO.N_gyroPhase,
                                            self.GYRO.N_vPhase,
                                            self.GYRO.N_vSlice,
                                            self.GYRO.N_HOT), dtype=int)*np.NaN
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

    def combinePFCpointcloud(self, runList, tPath):
        """
        Combines multiple pointclouds into a single pointcloud then saves to file
        pcName is the name of the point cloud we are dealing with.  Run this
        function once for each timestep
        """
        hfOptical = []
        hfGyro = []
        hfAll = []
        shadow =[]
        backface =[]
        shadowGyro = []
        bdotn = []
        psi = []
        norm = np.array([])
        bField = np.array([])

        names = []
        centers = np.array([])
        Npoints = 0
        for PFC in self.PFCs:
            #for tiles that are run in multiple mapDirections
            #add quantities that are superposition of multiple mapDirection runs
            if PFC.name in names:
                idx = names.index(PFC.name)
                if 'HFpc' in runList:
                    hfOptical[idx]+=PFC.qDiv
                    hfAll[idx]+=PFC.qDiv
                    #HFPC always runs shadowMask too
                    shadow[idx]+=PFC.shadowed_mask
                    shadow[idx].astype(bool)
                elif 'GyroPC' in runList:
                    hfGyro[idx]+=PFC.qGyro
                    hfAll[idx]+=PFC.qGyro
                    shadowGyro[idx]+=PFC.gyroShadowMask
                    shadowGyro[idx].astype(bool)
                elif 'backfacePC' in runList:
                    backface[idx]+=PFC.backfaceMask
                    backface[idx].astype(bool)
                #some quantities cant be added (like psi)
                #just use one mapDirection for these cases
                else:
                    pass
            else:
                if 'HFpc' in runList:
                    hfOptical.append(PFC.qDiv)
                    shadow.append(PFC.shadowed_mask)
                    hfAll.append(PFC.qDiv)
                if 'GyroPC' in runList:
                    hfGyro.append(PFC.qGyro)
                    shadowGyro.append(PFC.gyroShadowMask)
                    hfAll[-1]+=PFC.qGyro
                if 'backfacePC' in runList:
                    backface.append(PFC.backfaceMask)
                if 'bdotnPC' in runList:
                    bdotn.append(PFC.bdotn)
                if 'psiPC' in runList:
                    psi.append(PFC.psimin)
                if 'NormPC' in runList:
                    norm = np.append(norm, PFC.norms)
                if 'Bpc' in runList:
                    bField = np.append(bField, PFC.Bxyz)
                #note that I don't do the normal vector NormPC (its same every timestep)
                #user can just get NormPCs individually for each tile
                Npoints += len(PFC.centers)
                centers = np.append(centers,PFC.centers)
                names.append(PFC.name)

        #now build something we can write to csv (ie numpy)
        hfOpticalNumpy = np.array([])
        hfGyroNumpy = np.array([])
        hfAllNumpy = np.array([])
        shadowNumpy = np.array([])
        shadowGyroNumpy = np.array([])
        backfaceNumpy = np.array([])
        bdotnNumpy = np.array([])
        psiNumpy = np.array([])
        normNumpy = np.array([])
        for arr in hfOptical:
            hfOpticalNumpy = np.append(hfOpticalNumpy, arr)
        for arr in hfGyro:
            hfGyroNumpy = np.append(hfGyroNumpy, arr)
        for arr in hfAll:
            hfAllNumpy = np.append(hfAllNumpy, arr)
        for arr in shadow:
            shadowNumpy = np.append(shadowNumpy, arr)
        for arr in shadowGyro:
            shadowGyroNumpy = np.append(shadowGyroNumpy, arr)
        for arr in backface:
            backfaceNumpy = np.append(backfaceNumpy, arr)
        for arr in bdotn:
            bdotnNumpy = np.append(bdotnNumpy, arr)
        for arr in psi:
            psiNumpy = np.append(psiNumpy, arr)

        tag='all'
        centers = centers.reshape(Npoints,3)
        if 'HFpc' in runList:
            self.HF.write_heatflux_pointcloud(centers,hfOpticalNumpy,tPath,tag)
            PFC.write_shadow_pointcloud(centers,shadowNumpy,tPath,tag)
            self.HF.write_heatflux_pointcloud(centers,hfAllNumpy,tPath,tag,'all')
        if 'GyroPC' in runList:
            self.HF.write_heatflux_pointcloud(centers,hfGyroNumpy,tPath,tag,'gyro')
            PFC.write_shadow_pointcloud(centers,shadowGyroNumpy,tPath,tag,'gyro')
        if 'shadowPC' in runList:
            PFC.write_shadow_pointcloud(centers,shadowNumpy,tPath,tag)
        if 'backfacePC' in runList:
            PFC.write_backface_pointcloud(centers,backfaceNumpy,tPath,tag)
        if 'bdotnPC' in runList:
            PFC.write_bdotn_pointcloud(centers,bdotnNumpy, tPath,tag)
        if 'psiPC' in runList:
            self.HF.write_psiN_pointcloud(centers,psiNumpy,tPath,tag)
            #self.HF.write_psiNThresh_pointcloud(centers,psiNumpy,tPath,tag)
        if 'NormPC' in runList:
            norm = norm.reshape(Npoints,3)
            self.CAD.write_normal_pointcloud(centers,norm,tPath,tag)
        if 'Bpc' in runList:
            bField = bField.reshape(Npoints,3)
            self.MHD.write_B_pointcloud(centers, bField, tPath, tag)

        print("Wrote combined pointclouds")
        log.info("Wrote combined pointclouds")


    def combineTimeSteps(self, runList, t):
        """
        save composite csv from each timestep into a single directory for
        making paraview movies
        """
        movieDir = self.MHD.shotPath + '/paraview/'
        tPath = self.MHD.shotPath + '/' + '{:06d}/'.format(t)
        #first try to make new directory
        try:
            os.mkdir(movieDir)
        except FileExistsError:
            pass
        if 'HFpc' in runList:
            src = tPath + 'HF_optical_all.csv'
            dest = movieDir + 'hfOptical_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
        if 'shadowPC' in runList:
            src = tPath + 'shadowMask_optical_all.csv'
            dest = movieDir + 'shadowMask_optical_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
        if 'backfacePC' in runList:
            src = tPath + 'backfaceMask_all.csv'
            dest = movieDir + 'backfaceMask_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
        if 'bdotnPC' in runList:
            src = tPath + 'bdotnPointCloud_all.csv'
            dest = movieDir + 'bdotn_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
        if 'psiPC' in runList:
            src = tPath + 'psiPointCloud_all.csv'
            dest = movieDir + 'psiN_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
        if 'NormPC' in runList:
            src = tPath + 'NormPointCloud_all.csv'
            dest = movieDir + 'normals_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
        if 'Bpc' in runList:
            src = tPath + 'B_pointcloud_all.csv'
            dest = movieDir + 'Bfield_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
        if 'GyroPC' in runList:
            src = tPath + 'HF_gyro_all.csv'
            dest = movieDir + 'hfGyro_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
            src = tPath + 'HF_allSources_all.csv'
            dest = movieDir + 'hfAll_{:06d}.csv'.format(t)
            shutil.copy(src,dest)
            src = tPath + 'shadowMask_gyro_all.csv'
            dest = movieDir + 'shadowMask_gyro_{:06d}.csv'.format(t)
            shutil.copy(src,dest)

        return

    def psiPC(self, PFC):
        """
        creates a poloidal flux, psi, point cloud on tile surface
        you need to have run the cad, mhd, and hf initialization processes
        before running this function
        """
        PFC.shadowed_mask = np.zeros((len(PFC.shadowed_mask)))
        self.getPsiEverywhere(PFC, PFC.tag)

        print("Completed psiPC calculation")
        log.info("Completed psiPC calculation")

        return

    def getPsiEverywhere(self, PFC, save2File=True, tag=None):
        """
        get psi all over the PFC (including shadowed regions).
        """
        #Run MAFOT laminar for 3D plasmas
        if self.MHD.plasma3Dmask==True:
            CTLfile=PFC.controlfilePath + PFC.controlfile
            self.MHD.writeControlFile(CTLfile, PFC.t, PFC.mapDirection, mode='laminar')
            self.MHD.writeMAFOTpointfile(PFC.centers,PFC.gridfile)
            self.MHD.runMAFOTlaminar(PFC.gridfile,PFC.controlfilePath,PFC.controlfile,self.NCPUs)
            self.HF.readMAFOTLaminarOutput(PFC,PFC.outputFile)
            use = np.where(PFC.psimin < 10)[0]
            self.HF.write_psiN_pointcloud(PFC.centers,PFC.psimin,PFC.controlfilePath, PFC.tag)
            os.remove(PFC.outputFile)
        #get psi from gfile for 2D plasmas
        else:
            self.MHD.psi2DfromEQ(PFC)
            self.HF.write_psiN_pointcloud(PFC.centers,PFC.psimin,PFC.controlfilePath, PFC.tag)
        return


    def getDefaultDict(self):
        """
        returns an empty dict with each GUI parameter
        """
        emptyDict = {
                    'shot':None,
                    'tmin':None,
                    'tmax':None,
                    'nTrace': None,
                    'ROIGridRes': None,
                    'gridRes': None,
                    'STPfile': None,
                    'hfMode': None,
                    'lqEich': None,
                    'S': None,
                    'lqCN': None,
                    'lqCF': None,
                    'lqPN': None,
                    'lqPF': None,
                    'fracCN': None,
                    'fracCF': None,
                    'fracPN': None,
                    'fracPF': None,
                    'Psol': None,
                    'fracUI': None,
                    'fracUO': None,
                    'fracLI': None,
                    'fracLO': None,
                    'qBG' : None,
                    'OFtMin': None,
                    'OFtMax': None,
                    'deltaT': None,
                    'writeDeltaT': None,
                    'STLscale': None,
                    'meshMinLev': None,
                    'meshMaxLev': None,
                    'N_gyroSteps': None,
                    'gyroDeg': None,
                    'gyroT_eV': None,
                    'N_vSlice': None,
                    'N_vPhase': None,
                    'N_gyroPhase': None,
                    'ionMassAMU': None,
                    #'vMode': None,
                    'ionFrac': None
                    }
        return emptyDict

    def loadDefaults(self, inFile=None):
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
        tools.initializeInput(self.OF, self.infile)

        inputDict = {
                    'shot': self.MHD.shot,
                    'tmin': self.MHD.tmin,
                    'tmax': self.MHD.tmax,
                    'nTrace': self.MHD.nTrace,
                    'ROIGridRes': self.CAD.ROIGridRes,
                    'gridRes': self.CAD.gridRes,
                    'lqEich': self.HF.lqEich,
                    'S': self.HF.S,
                    'lqCN': self.HF.lqCN,
                    'lqCF': self.HF.lqCF,
                    'lqPN': self.HF.lqPN,
                    'lqPF': self.HF.lqPF,
                    'fracCN': self.HF.fracCN,
                    'fracCF': self.HF.fracCF,
                    'fracPN': self.HF.fracPN,
                    'fracPF': self.HF.fracPF,
                    'Psol': self.HF.Psol,
                    'fracUI':self.HF.fracUI,
                    'fracUO':self.HF.fracUO,
                    'fracLI':self.HF.fracLI,
                    'fracLO':self.HF.fracLO,
                    'qBG' : self.HF.qBG,
                    'OFtMin': self.OF.OFtMin,
                    'OFtMax': self.OF.OFtMax,
                    'deltaT': self.OF.deltaT,
                    'writeDeltaT': self.OF.writeDeltaT,
                    'STLscale': self.OF.STLscale,
                    'meshMinLev': self.OF.meshMinLev,
                    'meshMaxLev': self.OF.meshMaxLev,
                    'N_gyroSteps': self.GYRO.N_gyroSteps,
                    'gyroDeg': self.GYRO.gyroDeg,
                    'gyroT_eV': self.GYRO.gyroT_eV,
                    'N_vSlice': self.GYRO.N_vSlice,
                    'N_vPhase': self.GYRO.N_vPhase,
                    'N_gyroPhase': self.GYRO.N_gyroPhase,
                    'ionMassAMU': self.GYRO.ionMassAMU,
                    #'hfMode': self.HF.hfMode,
                    #'vMode': self.GYRO.vMode,
                    'ionFrac': self.GYRO.ionFrac
                    }
        print("Loaded defaults")

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

    def loadOF(self,OFtMin,OFtMax,OFminMeshLev,OFmaxMeshLev,
                      OFSTLscale, OFbashrc, OFdeltaT, OFwriteDeltaT, materialSelect):
        """
        loads user OF GUI settings

        OFstartTime is when we start OF simulation (can be before HF)
        OFstopTime is when we start OF simulation (can be after HF)
        OFminMeshLev is minimum refinement level for snappyhexmesh (default is 1)
        OFminMeshLev is maximum refinement level for snappyhexmesh (default is 3)
        OFSTLscale is scalar for unit conversion (default is 1)
        OFbashrc is file location on system to source OF binaries / libs
        OFdeltaT is timestep size for FVM simulation.  Defaults to 1ms (0.001s)
        """
        self.OF.OFtMin = float(OFtMin)/1000.0 #to [s] for openfoam
        self.OF.OFtMax = float(OFtMax)/1000.0 #to [s] for openfoam
        self.OF.meshMinLevel = OFminMeshLev
        self.OF.meshMaxLevel = OFmaxMeshLev
        self.OF.STLscale = OFSTLscale
        self.OF.cmd3Dmesh = 'meshAndPatch'
        #the OF bashrc must be sourced before running OF
        self.OF.cmdSourceOF = 'source ' + OFbashrc
        self.OF.OFbashrc = OFbashrc
        self.OF.cmdThermal = 'runThermal'
        self.OF.cmdTprobe = 'runTprobe'
        self.OF.deltaT = float(OFdeltaT) #this comes in [sec]
        self.OF.writeDeltaT = float(OFwriteDeltaT) #this comes in [sec]
        self.OF.materialSelect = materialSelect

        print("Material Selection: "+self.OF.materialSelect)

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
        print('Setting Up OF run')
        log.info('Setting Up OF run')
        #set up base OF directory for this discharge
        if self.MHD.shotPath[-1]=='/':
            self.OF.caseDir = self.MHD.shotPath + 'openFoam/heatFoam'
        else:
            self.OF.caseDir = self.MHD.shotPath + '/openFoam/heatFoam'
        tools.makeDir(self.OF.caseDir)
        #set up OF parts for each PFC part
        for PFC in self.PFCs:
            print("Running openFOAM for PFC: "+PFC.name)
            log.info("Running openFOAM for PFC: "+PFC.name)
            partDir = self.OF.caseDir + '/' + PFC.name
            self.OF.partDir = partDir
            self.OF.partName = PFC.name

            # tools.makeDir(partDir)
            #copy heatFoam template directory to this location
            try:
                shutil.copytree(self.OF.templateCase, partDir)
                if self.OF.OFtMin != 0:
                    t0 = partDir+'/0'
                    t0new = partDir+'/{:f}'.format(self.OF.OFtMin).rstrip('0').rstrip('.')
                    os.rename(t0, t0new)
            except OSError as e:
                print('COULD NOT COPY TEMPLATE DIRECTORY!  Aborting!')
                print(e)
                return

            #Set up partDir with user selected material properties
            print("Material Selection: "+self.OF.materialSelect)
            matDst = partDir + '/constant/'
            if self.OF.materialSelect == "ATJ":
                matSrc = self.OF.materialDir + '/ATJ/'
            elif self.OF.materialSelect == "MOLY":
                matSrc = self.OF.materialDir + '/MOLY/'
            else: #SGL6510
                matSrc = self.OF.materialDir + '/SGLR6510/'
            try:
                shutil.copyfile(matSrc+'DT', matDst+'DT')
                shutil.copyfile(matSrc+'thermCond', matDst+'thermCond')
            except OSError as e:
                print('COULD NOT COPY MATERIAL PROPERTIES!  Aborting!')
                print(e)

            #set up timesteps
            tMin = self.OF.OFtMin*1000.0 #in [ms] for HEAT
            tMax = self.OF.OFtMax*1000.0 #in [ms] for HEAT
            arr = np.linspace(int(tMin), int(tMax), int((tMax-tMin)+1), dtype=int)
            OFtimesteps = arr[0::int(self.OF.deltaT*1000.0)]

            #create symbolic link to STL file
            print("Creating openFOAM symlink to STL")
            log.info("Creating openFOAM symlink to STL")
            #old method used ROIGridRes, but this can sometimes create non-watertight mesh
            #new method uses standard FreeCAD mesher (not mefisto) to make better volume mesh
            #PFC.OFpart = PFC.name + "___" + self.CAD.ROIGridRes
            PFC.OFpart = PFC.name + "___" + "standard"
            triSurfaceLocation = partDir+'/constant/triSurface/' + PFC.OFpart +"mm.stl"

            if self.CAD.STLpath[-1] == '/':
                stlfile = self.CAD.STLpath + PFC.OFpart +"mm.stl"
            else:
                stlfile = self.CAD.STLpath +'/'+ PFC.OFpart +"mm.stl"

            #if the standard mesh for this part doesn't exist, create it using freecad
            #in the HEAT CAD module
            if os.path.exists(stlfile)==True:
                pass
            else:
                partIdx = np.where(self.CAD.ROI == PFC.name)[0][0]
                part = self.CAD.ROIparts[partIdx]
                meshSTL = self.CAD.part2meshStandard(part)
                self.CAD.writeMesh2file(meshSTL, PFC.name, path=self.CAD.STLpath, resolution='standard')

            #create hard link to STL
            os.link(stlfile,triSurfaceLocation)

            #Test STL to make sure it is watertight
            meshInQuestion = trimesh.load(stlfile)
            validMesh = meshInQuestion.is_watertight
            if validMesh == True:
                print("STL file for "+PFC.name+" is watertight")
                log.info("STL file for "+PFC.name+" is watertight")
            else:
                print("\n=====================================================")
                print("WARNING!!!  STL file for "+PFC.name+"is NOT watertight.")
                print("OpenFOAM thermal analysis will probably fail")
                print("=====================================================\n")
                log.info("\n=====================================================")
                log.info("WARNING!!!  STL file for "+PFC.name+"is NOT watertight.")
                log.info("OpenFOAM thermal analysis will probably fail")
                print("=====================================================\n")

            #update middle points and blockmesh bounds for each PFC and
            # give 10mm of clearance on each side
            self.OF.xMin = (PFC.centers[:,0].min() - 0.01)*1000.0
            self.OF.xMax = (PFC.centers[:,0].max() + 0.01)*1000.0
            self.OF.yMin = (PFC.centers[:,1].min() - 0.01)*1000.0
            self.OF.yMax = (PFC.centers[:,1].max() + 0.01)*1000.0
            self.OF.zMin = (PFC.centers[:,2].min() - 0.01)*1000.0
            self.OF.zMax = (PFC.centers[:,2].max() + 0.01)*1000.0
            self.OF.xMid = (self.OF.xMax-self.OF.xMin)/2.0 + self.OF.xMin
            self.OF.yMid = (self.OF.yMax-self.OF.yMin)/2.0 + self.OF.yMin
            self.OF.zMid = (self.OF.zMax-self.OF.zMin)/2.0 + self.OF.zMin

#            #setup openfoam environment
#            try:
#                AppImage = os.environ["APPIMAGE"]
#                AppDir = os.environ["APPDIR"]
#                inAppImage = True
#            except:
#                inAppImage = False
#                AppDir = ''
#            if inAppImage == True:
#                print("Setting up OF apppImage environment")
#                log.info("Setting up OF apppImage environment")
#                OFplatformDir = AppDir + '/usr/share/openfoam/platforms/linux64GccDPInt32pt'
#                OFbinDir = OFplatformDir + '/bin'
#                OFlibDir = OFplatformDir + '/lib'
#                #make openfoam platform directory
#                os.mkdirs(OFplatformDir, exist_ok=True)
#                #symlink AppDir/usr/bin to openfoam bin (same for lib)
#                try:
#                    os.symlink('/usr/bin', OFbinDir)
#                    os.symlink('/usr/lib', OFlibDir)
#                except:
#                    print("could not link openfoam libs and bins")
#                    log.info("could not link openfoam libs and bins")



            #dynamically write template variables to templateVarFile
            print("Building openFOAM templates and shell scripts")
            log.info("Building openFOAM templates and shell scripts")
            templateVarFile = partDir + '/system/templateVariables'
            STLpart = PFC.OFpart +"mm.stl"
            self.OF.writeOFtemplateVarFile(templateVarFile, STLpart)
            self.OF.writeShellScript(self.logFile)

            #create OF dictionaries from templates
            self.OF.createDictionaries(self.OF.templateDir,
                                       partDir,
                                       templateVarFile,
                                       stlfile)

            #generate 3D volume mesh or coy from file
            print("Generating volume mesh")
            log.info("Generating volume mesh")
            self.OF.generate3Dmesh(PFC.OFpart)

            qDiv = np.zeros((len(PFC.centers)))
            ctrs = copy.copy(PFC.centers)*1000.0

            #cycle through timesteps and get HF data from HEAT tree
            for t in OFtimesteps:
                print("openFOAM timestep: {:d}".format(t))
                log.info("openFOAM timestep: {:d}".format(t))
                OFt = t/1000.0 #in [s] for openFOAM

                #copy variable files into each timestep for solver
                timeDir = partDir + '/{:f}'.format(OFt).rstrip('0').rstrip('.')
                #timestep field prescriptions
                HFt0 = partDir+'/{:f}'.format(self.OF.OFtMin).rstrip('0').rstrip('.')+'/HF'
                HFtStep = partDir+'/{:f}'.format(OFt).rstrip('0').rstrip('.')+'/HF'
                try:
                    #shutil.copytree(t0new,timeDir
                    os.mkdir(timeDir)
                    #shutil.copy(HFt0, HFtStep)
                    shutil.copy(HFt0, timeDir)

                except:
                    print("***")
                    print("Could not create OF directory for timestep {:d}".format(t))
                    print("(this is expected for t=0)")
                    print("***")
                    pass


                # determine heat flux boundary condition
                if (t in PFC.timesteps) and (t in self.MHD.timesteps):
                    #we explicitly calculated HF for this timestep
                    print("OF.timestep: {:d} in PFC.timesteps".format(t))
                    log.info("OF.timestep: {:d} in PFC.timesteps".format(t))
                    if self.MHD.shotPath[-1]=='/':
                        HFcsv = self.MHD.shotPath + '{:06d}/'.format(t) + PFC.name + '/HF_allSources.csv'
                    else:
                        HFcsv = self.MHD.shotPath + '/' + '{:06d}/'.format(t) + PFC.name + '/HF.allSources.csv'
                    qDiv = pd.read_csv(HFcsv)['HeatFlux'].values
                    #write boundary condition
                    print("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    self.HF.write_openFOAM_boundary(ctrs,qDiv,partDir,OFt)
                elif (t < self.MHD.timesteps.min()) or (t > self.MHD.timesteps.max()):
                    #apply zero HF outside of discharge domain (ie tiles cooling)
                    print("OF.timestep: {:d} outside MHD domain".format(t))
                    log.info("OF.timestep: {:d} outside MHD domain".format(t))
                    qDiv = np.zeros((len(ctrs)))
                    #write boundary condition
                    print("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    self.HF.write_openFOAM_boundary(ctrs,qDiv,partDir,OFt)
                else:
                    #boundary using last timestep that we calculated a HF for
                    #(basically a heaviside function in time)
                    #print("OF.timestep: {:d} using heaviside from last PFC.timestep".format(t))
                    #log.info("OF.timestep: {:d} using heaviside from last PFC.timestep".format(t))
                    #write boundary condition
                    #print("Maximum qDiv for this PFC and time: {:f}".format(qDiv.max()))
                    #self.HF.write_openFOAM_boundary(ctrs,qDiv,partDir,OFt)

                    #openFOAM linear interpolation in time using timeVaryingMappedFixedValue
                    print("HF being linearly interpolated by OF at this t")
                    pass


            #run openfoam thermal analysis using heatFoam solver
            self.OF.runThermalAnalysis()
            print("thermal analysis complete...")
            log.info("thermal analysis complete...")

        print("openFOAM run completed.")
        log.info("openFOAM run completed.")
        return

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
            if self.MHD.shotPath[-1]=='/':
                partDir = self.MHD.shotPath + 'openFoam/heatFoam/'+PFC.name
            else:
                partDir = self.MHD.shotPath + '/openFoam/heatFoam/'+PFC.name
            file = (partDir +
                    '/postProcessing/fieldMinMax1/{:f}'.format(self.OF.OFtMin).rstrip('0').rstrip('.')
                    +'/fieldMinMax.dat')
            data.append(self.OF.getMinMaxData(file))
            pfcNames.append(PFC.name)

        fig = pgp.plotlyOpenFOAMplot(data,pfcNames)

        #save interactive plotly plot in shotPath/plotly/OFminmax.html
        if self.MHD.shotPath[-1]=='/':
            plotlyDir = self.MHD.shotPath + 'plotly'
        else:
            plotlyDir = self.MHD.shotPath + '/plotly'
        try:
            os.mkdir(plotlyDir)
        except:
            print("Could not make plotly directory")
            log.info("Could not make plotly directory")

        plotPath = plotlyDir + '/OFminmax.html'
        fig.write_html(plotPath)
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
        if self.MHD.shotPath[-1]=='/':
            plotlyDir = self.MHD.shotPath + 'plotly'
        else:
            plotlyDir = self.MHD.shotPath + '/plotly'
        try:
            os.mkdir(plotlyDir)
        except:
            print("Could not make plotly directory")
            log.info("Could not make plotly directory")

        plotPath = plotlyDir + '/HFdist.html'
        fig.write_html(plotPath)

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
        if self.MHD.shotPath[-1]=='/':
            plotlyDir = self.MHD.shotPath + 'plotly'
        else:
            plotlyDir = self.MHD.shotPath + '/plotly'
        try:
            os.mkdir(plotlyDir)
        except:
            print("Could not make plotly directory")
            log.info("Could not make plotly directory")

        plotPath = plotlyDir + '/Tprobes.html'
        fig.write_html(plotPath)
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
                                   self.GYRO.vScan[0])
        return fig

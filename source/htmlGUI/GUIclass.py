#GUIclass.py
#Description:   GUI Module
#Engineer:      T Looby
#Date:          20200115
"""
This is the launch point for the HEAT code, when accessing from the HTML GUI.
It calls HEAT functions and generates heat fluxes based upon user input.  These
functions are usually called from the HEATgui.py script, which is the FLASK
binding to html.
"""

import CADClass
import MHDClass
import toolsClass
import heatfluxClass
import openFOAMclass
import time
import numpy as np
import logging
import os
import pandas as pd
import shutil
import errno
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

class GUIobj():
    def __init__(self, logFile, rootDir):
        self.logFile = logFile
        self.rootDir = rootDir
        self.initializeEveryone()


    def machineSelect(self, MachFlag):
        """
        Select a machine and set the necessary paths
        """
        self.MachFlag = MachFlag
        self.setInitialFiles()

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
            self.infile = '../inputs/NSTXU_input.csv'
            self.PartsFile = '../inputs/NSTXUparts.csv'
            self.IntersectFile = '../inputs/NSTXUintersects.csv'
            self.OF.meshDir = '/u/tlooby/NSTX/CAD/HEAT/3Dmeshes'

        elif self.MachFlag == 'st40':
            print('Loading ST40 Input Filestream')
            log.info('Loading ST40 Input Filestream')
            self.infile = '../inputs/ST40_input.csv'
            self.PartsFile = '../inputs/ST40parts.csv'
            self.IntersectFile = '../inputs/ST40intersects.csv'
            self.OF.meshDir = '/u/tlooby/ST40/CAD/HEAT/3Dmeshes'

        elif self.MachFlag == 'd3d':
            print('Loading DIII-D Input Filestream')
            log.info('Loading DIII-D Input Filestream')
            self.infile = '../inputs/D3D_input.csv'
            self.PartsFile = '../inputs/D3Dparts.csv'
            self.IntersectFile = '../inputs/D3Dintersects.csv'
            self.OF.meshDir = '/u/tlooby/D3D/CAD/HEAT/3Dmeshes'

        self.OF.templateCase = '../openFoamTemplates/heatFoamTemplate'
        self.OF.templateDir = '../openFoamTemplates/templateDicts'



    def initializeEveryone(self):
        """
        Create objects that we can reference later on
        """
        self.MHD = MHDClass.MHD()
        self.CAD = CADClass.CAD()
        self.HF = heatfluxClass.heatFlux()
        self.OF = openFOAMclass.OpenFOAM()
        return

    def getMHDInputs(self,shot=None,tmin=None,tmax=None,nTrace=None,gfile=None):
        """
        Get the mhd inputs from the gui or input file
        """
        self.MHD.allowed_class_vars()
        tools.vars2None(self.MHD)
        self.MHD.MachFlag=self.MachFlag
        tools.read_input_file(self.MHD, infile=self.infile)
        self.MHD.setTypes()

        if shot is not None:
            self.MHD.shot = shot
        if tmin is not None:
            self.MHD.tmin = tmin
        if tmax is not None:
            self.MHD.tmax = tmax
        if nTrace is not None:
            self.MHD.nTrace = nTrace
        self.MHD.gfile = gfile
        if gfile is not None:
            path, name = os.path.split(gfile)
            gshot, time = name[:13].split('.')
            shot = gshot[1:]
            self.MHD.shot = int(shot)
            #self.MHD.time = time



        self.MHD.tree = 'EFIT02'
        self.MHD.dataPath = self.MHD.dataPath + self.MHD.MachFlag +"_{:06d}".format(self.MHD.shot)
        self.MHD.get_mhd_inputs('nstx',self.MHD.gfile)
        if type(self.MHD.timesteps) == int:
            self.t = self.MHD.timesteps  #NEED TO ADD MULTIPLE TIMESTEPS
        else:
            self.t = self.MHD.timesteps[0]
        self.MHD.make1EFITobject(self.MHD.shot, self.t, self.MHD.gfile)
        self.NCPUs = 4
        self.controlfile = '_lamCTL.dat'
        self.controlfilePath = self.MHD.dataPath + '/' + '{:06d}/'.format(self.t)
        self.gridfile = self.MHD.dataPath + '/' + '{:06d}/grid.dat'.format(self.t)
        self.outputFile = self.controlfilePath + 'lam.dat'
        self.MHD.setTypes()

    def gfileClean(self, psiRZMult,psiSepMult,psiAxisMult,FpolMult):
        """
        multiplies values in MHD ep object with scalars defined by user in html gui
        """
        print("psiRZ Multiplier = {:f}".format(psiRZMult))
        log.info("psiRZ Multiplier = {:f}".format(psiRZMult))
        print("psiSep Multipltier = {:f}".format(psiSepMult))
        log.info("psiSep Multipltier = {:f}".format(psiSepMult))
        print("psiAxis Multipltier = {:f}".format(psiAxisMult))
        log.info("psiAxis Multipltier = {:f}".format(psiAxisMult))
        print("Fpol Multiplie = {:f}".format(FpolMult))
        log.info("Fpol Multiplie = {:f}".format(FpolMult))

        self.MHD.ep.g['psiRZ'] *= psiRZMult
        self.MHD.ep.g['psiSep'] *= psiSepMult
        self.MHD.ep.g['psiAxis'] *= psiAxisMult
        self.MHD.ep.g['Fpol'] *= FpolMult

        psi = self.MHD.ep.g['psiRZ']
        psiSep = self.MHD.ep.g['psiSep']
        psiAxis = self.MHD.ep.g['psiAxis']

        self.MHD.ep.g['psiRZn'] = (psi - psiAxis) / (psiSep - psiAxis)
        return

    def writeGfile(self, newShot, newTime, newGfile):
        """
        writes a new gfile
        """
        self.MHD.writeGfile(newGfile, newShot, newTime)
        return

    def getCADInputs(self,ROIGridRes=None,gridRes=None,STPfile=None):
        #Load CAD and generate ROI meshes
        import numpy as np
        tools.initializeInput(self.CAD, infile=self.infile)
        self.CAD.rootDir = self.rootDir #set HEAT rootDir from HEATgui.py
        if ROIGridRes is not None:
            self.CAD.ROIGridRes = ROIGridRes
        if gridRes is not None:
            self.CAD.gridRes = gridRes
        if STPfile is not None:
            self.CAD.STPfile = STPfile
        self.CAD.loadSTEP()
        self.CAD.getROIfromfile(self.PartsFile)
        self.CAD.getROImeshes()
        self.CAD.writeMesh2file(self.CAD.ROImeshes, self.CAD.ROI, path=self.CAD.STLpath)

        #Find potential intersections by file and generate corresponding meshes
        intersects = self.CAD.readIntersects(self.IntersectFile, self.CAD.CADparts)
        print("Potential intersects on these tiles:")
        log.info("Potential intersects on these tiles:")
        for i in range(len(self.CAD.CADparts)):
            if intersects[i] == 1.0:
                print(self.CAD.CADparts[i].Label)
                log.info(self.CAD.CADparts[i].Label)
        self.CAD.gridRes = gridRes
        self.targetMeshes, self.targetLabels = self.CAD.meshPotentialIntersects(intersects, self.CAD.ROIparts,self.CAD.CADparts, resolution=gridRes)
        print('Meshed intersections')
        self.CAD.writeMesh2file(self.targetMeshes, self.targetLabels, path=self.CAD.STLpath,resolution=gridRes)
        print('Wrote intersect meshes to file')
        return

    def getHFInputs(self,lq=None,S=None,Psol=None,qBG=None):
        """
        get heat flux inputs from gui or input file
        """
        self.initializeHF()

        if lq is not None:
            self.HF.lq = lq
        if S is not None:
            self.HF.S = S
        if Psol is not None:
            self.HF.Psol = Psol
        if qBG is not None:
            self.HF.qBG = qBG
        return

    def bfieldAtSurface(self):
        """
        Calculate the B field at tile surface
        """
        ctrs = self.CAD.ROIctrs[0]/1000.0
        R,Z,phi = tools.xyz2cyl(ctrs[:,0],ctrs[:,1],ctrs[:,2])
        Bxyz = self.MHD.Bfield_pointcloud(self.MHD.ep, R, Z, phi)
        self.MHD.write_B_pointcloud(ctrs,Bxyz,self.controlfilePath)

    def Btrace(self,x,y,z):
        """
        Run a MAFOT structure trace from a point defined in the gui
        """
        xyz = np.array([x,y,z])
        controlfile = '_structCTL.dat'
        dphi = 1.0
        #self.MHD.MapDirection = 0
        #self.MHD.ittStruct = 100.0
        self.MHD.ittStruct = self.MHD.nTrace+1
        gridfile = self.MHD.dataPath + '/' + '{:06d}/struct_grid.dat'.format(self.t)
        self.MHD.writeControlFile(controlfile, self.t, mode='struct')
        self.MHD.writeMAFOTpointfile(xyz,gridfile)
        self.MHD.getFieldpath(dphi, gridfile, self.controlfilePath, controlfile, paraview_mask=True)

    def NormPC(self):
        """
        create a normal vector point cloud for mesh centers on tile surface
        """
        self.HF.makePFCs(self.MHD, self.CAD.ROImeshes, self.CAD.ROIctrs,
                        self.CAD.ROInorms, self.CAD.ROIareas, self.CAD.ROIparts)
        PFC = self.HF.PFCs[0]
        self.CAD.write_normal_pointcloud(PFC.centers,PFC.norms,self.controlfilePath)

    def shadowPC(self, HFpc_flag):
        """
        create a pointcloud for mesh center locations where 1=shadowed, 0=not shadowed
        """
        self.HF.makePFCs(self.MHD, self.CAD.ROImeshes, self.CAD.ROIctrs,
                        self.CAD.ROInorms, self.CAD.ROIareas, self.CAD.ROIparts)

        PFC = self.HF.PFCs[0]

        #Check for intersections with MAFOT struct if it hasnt been done by HFPC already
        print("INTERSECTION CHECK INITIALIZED")
        log.info("INTERSECTION CHECK INITIALIZED")
        #self.HF.findShadows_structure(self.MHD,PFC,self.targetMeshes)

        #Write Point Cloud
        self.HF.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,self.controlfilePath)


    def initializeHF(self):
        """
        Initialize heat flux variables
        """
        print("-"*70)
        print("HEAT FLUX MODULE INITIALIZED")
        log.info("HEAT FLUX MODULE INITIALIZED")
        #Initialize HF Object
        tools.initializeInput(self.HF, infile=self.infile)


    def HFPC(self):
        """
        Run a heat flux calculation.  This is called from gui by user.
        Calculates eich profile and creates heat flux point clouds where
        heat flux is defined at each mesh center
        """
        t0 = time.time()
        self.HF.makePFCs(self.MHD, self.CAD.ROImeshes, self.CAD.ROIctrs,
                        self.CAD.ROInorms, self.CAD.ROIareas, self.CAD.ROIparts)

        #PFC Object
        PFC = self.HF.PFCs[0]
        #self.HF.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,self.controlfilePath)

        #Check for intersections with MAFOT struct
        print("INTERSECTION CHECK INITIALIZED")
        log.info("INTERSECTION CHECK INITIALIZED")
        self.HF.findShadows_structure(self.MHD,PFC,self.targetMeshes)

        #Run MAFOT laminar
        print('\n')
        print("-"*70)
        print("MAFOT LAMINAR MODULE INITIALIZED")
        log.info("MAFOT LAMINAR MODULE INITIALIZED")
        self.MHD.writeControlFile(self.controlfile, self.t, mode='laminar')
        use = np.where(PFC.shadowed_mask != 1)[0]
        self.MHD.writeMAFOTpointfile(PFC.centers[use],self.gridfile)
        self.MHD.runMAFOTlaminar(self.gridfile,self.controlfilePath,self.controlfile,self.NCPUs)

        #Create Heat Flux Profile
        print('\n')
        print("-"*70)
        print("HEAT FLUX CALCULATION")
        log.info("HEAT FLUX CALCULATION")
        self.HF.readMAFOTLaminarOutput(PFC,self.outputFile)
        q = self.HF.getHFprofile(PFC, self.HF.Psol, self.HF.lq, self.HF.S, self.HF.qBG)
        qDiv = self.HF.q_div(PFC, self.MHD, q)
        os.remove(self.outputFile)

        print('Maximum heat load on tile: {:f}'.format(max(qDiv)))
        print('Input Power = {:f}'.format(self.HF.Psol))
        print('Tessellated Total Power = {:f}'.format(self.HF.power_sum_mesh(PFC)))
        log.info('Maximum heat load on tile: {:f}'.format(max(qDiv)))
        log.info('Input Power = {:f}'.format(self.HF.Psol))
        log.info('Tessellated Total Power = {:f}'.format(self.HF.power_sum_mesh(PFC)))

        #Create pointclouds for paraview
        R,Z,phi = tools.xyz2cyl(PFC.centers[:,0],PFC.centers[:,1],PFC.centers[:,2])
        self.HF.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,self.controlfilePath)
        self.HF.write_heatflux_pointcloud(PFC.centers,qDiv,self.controlfilePath)
        #structOutfile = MHD.dataPath + '/' + '{:06d}/struct.csv'.format(PFC.t)
        #HF.PointCloudfromStructOutput(structOutfile)

        #Save data to class variable for future use
        PFC.q = q
        PFC.qDiv = qDiv
        self.HF.PFCs[0] = PFC
        print("\nTime Elapsed: {:f}".format(time.time() - t0))
        log.info("\nTime Elapsed: {:f}".format(time.time() - t0))

        print("Completed HFPC run")
        log.info("Completed HFPC run")

    def psiPC(self):
        """
        creates a poloidal flux, psi, point cloud on tile surface
        you need to have run the cad, mhd, and hf initialization processes
        before running this function
        """
        self.HF.makePFCs(self.MHD, self.CAD.ROImeshes, self.CAD.ROIctrs,
                        self.CAD.ROInorms, self.CAD.ROIareas, self.CAD.ROIparts)
        PFC = self.HF.PFCs[0]
        PFC.shadowed_mask = np.zeros((len(PFC.shadowed_mask)))
        #Run MAFOT laminar
        print('\n')
        print("-"*70)
        print("Writing psi point cloud")
        log.info("Writing psi point cloud")
        self.MHD.writeControlFile(self.controlfile, self.t, mode='laminar')
        self.MHD.writeMAFOTpointfile(PFC.centers,self.gridfile)
        self.MHD.runMAFOTlaminar(self.gridfile,self.controlfilePath,self.controlfile,self.NCPUs)
        self.HF.readMAFOTLaminarOutput(PFC,self.outputFile)
        use = np.where(PFC.psimin < 10)[0]
        self.HF.write_psiN_pointcloud(PFC.centers[use],PFC.psimin[use],self.controlfilePath)
        os.remove(self.outputFile)
        print("Completed psiPC calculation")
        log.info("Completed psiPC calculation")


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
        tools.initializeInput(self.OF, self.infile)

        inputDict = {
                    'shot': self.MHD.shot,
                    'tmin': self.MHD.tmin,
                    'tmax': self.MHD.tmax,
                    'nTrace': self.MHD.nTrace,
                    'ROIGridRes': self.CAD.ROIGridRes,
                    'gridRes': self.CAD.gridRes,
                    'STPfile': self.CAD.STPfile,
                    'lq': self.HF.lq,
                    'S': self.HF.S,
                    'Psol': self.HF.Psol,
                    'qBG' : self.HF.qBG,
                    'xMin': self.OF.xMin,
                    'xMax': self.OF.xMax,
                    'yMin': self.OF.yMin,
                    'yMax': self.OF.yMax,
                    'zMin': self.OF.zMin,
                    'zMax': self.OF.zMax,
                    'tMin': self.OF.tMin,
                    'tMax': self.OF.tMax,
                    'deltaT': self.OF.deltaT,
                    'writeDeltaT': self.OF.writeDeltaT,
                    'xProbe': self.OF.xProbe,
                    'yProbe': self.OF.yProbe,
                    'zProbe': self.OF.zProbe,
                    'STLscale': self.OF.STLscale,
                    'meshMinLevel': self.OF.meshMinLevel,
                    'meshMaxLevel': self.OF.meshMaxLevel,
                    'xMid': self.OF.xMid,
                    'yMid': self.OF.yMid,
                    'zMid': self.OF.zMid,
                    'STLfileName': self.OF.STLfileName,
                    'STLlayerName': self.OF.STLlayerName
                    }


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

    def getOpenFOAMinputs(self,OFsettingsDict=None):
        """
        Reads openfoam input variables from HEAT input file or from HEAT gui
        if OFsettingsDict=None, then we assume we are reading from a file (not gui)
                OFsettingsDict is a dictionary that contains all variables for OF run
                as described in (openFOAMclass.OpenFOAM.setTypes)
        """
        print('\n')
        print("-"*55)
        print("OPENFOAM MODULE INITIALIZED")
        log.info('\n')
        log.info("-"*55)
        log.info("OPENFOAM MODULE INITIALIZED")

        self.OF.allowed_class_vars()
        tools.vars2None(self.OF)
        if OFsettingsDict==None:
            tools.read_input_file(self.OF, infile=self.infile)
            self.OF.dict = tools.createDict(self.OF)
        else:
            tools.inputs_from_dict(self.OF, OFsettingsDict)
            self.OF.dict = tools.createDict(self.OF)
        self.OF.setTypes()


    def setUpOpenFOAM(self,STLpart,clobberFlag=True):
        """
        sets up OpenFOAM finite element simulation to determine temperature distribution
        in a PFC.  If you are running OF for multiple PFC tiles then you must
        call this function once for each STLpart.

        STLpart is name of STL file (ie E-ED1434-011_5.0mm)

        caseDir is solver case directory
        partDir is <caseDir>/<part> where <part> is the CAD part.  Each solver
            caseDir has multiple partDir inside, corresponding to running that
            solver for each part
        """
        print('Setting Up OF run')
        log.info('Setting Up OF run')
        #strip .stl extension for part name
        if STLpart[-4:] == '.stl':
            part = STLpart[:-4]
        else:
            part = STLpart
            STLpart = STLpart+'.stl'

        #Variables we will use for OF bash scripts
        self.OF.caseDir = self.controlfilePath + 'openFoam/heatFoam'
        self.OF.partDir = self.controlfilePath + 'openFoam/heatFoam/'+part
        self.OF.cmd3Dmesh = 'meshAndPatch'
        self.OF.cmdSourceOF = 'source /opt/OpenFOAM/OpenFOAM-v1912/etc/bashrc'
        self.OF.cmdThermal = 'runThermal'
        self.OF.cmdTprobe = 'runTprobe'
        self.OF.paraviewDir = self.controlfilePath+'paraview'
        self.OF.TprobeFile = self.OF.partDir+'/postProcessing/probes/0/T'

        #Try to make new OFcaseDir root directory for heatFoam
        try:
            os.mkdir(self.OF.caseDir)
        except:
            print('Could not make OF case directory')
            if clobberFlag == True:
                print('Overwriting because clobberFlag=True.')
                os.makedirs(self.OF.caseDir, exist_ok=True)

        #Try to make a caseDir for this specific part inside the directory we just made
#        try:
#            os.mkdir(self.OF.partDir)
#        except:
#            print('COULD NOT MAKE OpenFOAM PART DIR')

        #copy OF heatFoam template case to new case directory location
        try:
            shutil.copytree(self.OF.templateCase, self.OF.partDir)
        except:
            print('COULD NOT COPY TEMPLATE DIRECTORY')

        #copy STLfile from gui.CAD.STLpath to OFpartDir/constant/triSurface/
        triSurfaceLocation = self.OF.partDir+'/constant/triSurface'
        if self.CAD.STLpath[-1] == '/':
            STLlocation = self.CAD.STLpath + STLpart
        else:
            STLlocation = self.CAD.STLpath + '/' + STLpart
        print('Copying STL file: '+STLlocation)
        log.info('Copying STL file: '+STLlocation)
        try:
            shutil.copy(STLlocation, triSurfaceLocation)
        except:
            print('Could not copy STL file to openfoam directory')

        #dynamically write template variables to templateVarFile
        templateVarFile = self.OF.partDir + '/system/templateVariables'
        self.OF.writeOFtemplateVarFile(templateVarFile, STLpart)

        self.OF.writeShellScript(self.logFile)

        #get list of STL parts
        #parts = pd.read_csv(self.PartsFile, sep=',', comment='#', names=['parts'], skipinitialspace=True)
        #STLfiles = parts['parts'].to_list()

        #create OF dictionaries from templates
        self.OF.createDictionaries(self.OF.templateDir,
                                   self.OF.partDir,
                                   templateVarFile,
                                   STLpart)

        print('OF case directory ready for execution')
        log.info('OF case directory ready for execution')


    def runOpenFOAM(self, STLpart):
        """
        Runs OpenFOAM case that was setup with self.setUpOpenFOAM
        Currently configured to solve heat diffusion equation using temperature
        dependent thermal diffusion constant for single region
        (no multiregion analysis)

        dict is the dictionary of
        """
        #strip .stl extension for part name
        if STLpart[-4:] == '.stl':
            part = STLpart[:-4]
        else:
            part = STLpart
            STLpart = STLpart+'.stl'

        #Create Boundary Condition for openFOAM from HEAT heatflux data
        self.HF.write_openFOAM_boundary(self.HF.PFCs[0].centers,self.HF.PFCs[0].qDiv,self.OF.partDir)

        #generate 3D mesh if necessary
        self.OF.generate3Dmesh(part)

        #run openfoam thermal analysis using heatFoam solver

        self.OF.runThermalAnalysis()

    def TprobeOF(self,x,y,z):
        """
        if Tprobe_mask is true run temperature probe: postProcess -func "probes"
        """
        self.OF.runTprobe(x,y,z)
        self.OF.plotTprobes(self.OF.TprobeFile)

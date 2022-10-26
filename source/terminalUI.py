#terminalUI.py
#Description:   Terminal user interface (UI) to HEAT
#Engineer:      T Looby
#Date:          20220131
"""
provides an interface that can be run from a linux or python terminal.  connects
to the HEAT engine (engineClass.py).  Allows for single runs or batch mode.
"""
import os
import sys
import shutil
import numpy as np
import pandas as pd
import EFIT.equilParams_class as EP
import toolsClass
tools = toolsClass.tools()
import logging
log = logging.getLogger(__name__)
import multiprocessing
import time


#get relevant environment variables.  Needed for containers
logFile = os.environ["logFile"]
rootDir = os.environ["rootDir"]
dataPath = os.environ["dataPath"]
OFbashrc = os.environ["OFbashrc"]
FreeCADPath = os.environ["FreeCADPath"]
PVPath = os.environ["PVPath"]
pvpythonCMD = os.environ["pvpythonCMD"]
try:
    AppDir = os.environ["APPDIR"]
except:
    AppDir = 'Not in appImage'

try:
    chmod = int(os.environ["HEATchmod"], 8) #convert from base 8
except:
    chmod = int(0o774, 8)

try:
    GID = int(os.environ["dockerGID"]) #group ID
except:
    GID = -1
try:
    UID = int(os.environ["dockerUID"]) #user ID
except:
    UID = -1


print("CHMOD: " + (oct(chmod)))
print("UID: {:d}".format(UID))
print("GID: {:d}".format(GID))

#Import HEAT engine class
from engineClass import engineObj

class TUI():
    def __init__(self):
        """
        intialize terminal user interface (TUI) object
        """
        self.ENG = engineObj(logFile, rootDir, dataPath, OFbashrc, chmod, UID, GID)
        self.ENG.NCPUs = multiprocessing.cpu_count() - 2 #reserve 2 cores for overhead
        self.chmod = chmod
        self.GID = GID
        self.UID = UID
        #if data directory doesn't exist, create it
        tools.makeDir(dataPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        return


    def simulationSchedule(self, batchFile):
        """
        determines what simulations need to be run and their correct grouping
        via a batchFile
        """
        print("Reading batch file")
        log.info("Reading batch file")

        self.caseDir= os.path.dirname(batchFile)

        #read batch file
        data = pd.read_csv(batchFile, sep=',', comment='#', skipinitialspace=True)

        #determine simulation schedule
        machines = np.unique(data['MachFlag'].values)
        for mach in machines:
            if mach not in self.machineList:
                print("\n\nMachFlag was not properly set in batchFile!")
                print("You provided "+mach+", which is not in the machineList.")
                print("Machines must be one of the following: ")
                print(self.machineList)
                sys.exit()

        machData = []
        machTags = []
        #split data by machine
        for mach in machines:
            machData.append(data[data["MachFlag"]==mach])
            #split data by tags
            machTags.append(np.unique(machData[-1]['Tag']))

        self.machines = machines
        self.machData = machData
        self.machTags = machTags
        self.Nsim = sum(len(l) for l in machTags)
        self.batchData = data

        print("Number of simulations to be scheduled from batchFile: {:d}".format(self.Nsim))
        log.info("Number of simulations to be scheduled from batchFile: {:d}".format(self.Nsim))
        return


    def runSimulations(self):
        """
        run simulations in schedule
        """
        for i,mach in enumerate(self.machines):
            self.ENG.machineSelect(mach, self.machineList)
            data = self.machData[i]
            machInDir = self.caseDir + '/' + mach + '/'
            for j,tag in enumerate(self.machTags[i]):
                print('\n')
                print("-"*70)
                print(" "*20 + "Machine: "+mach+"   Tag: "+tag)
                print("-"*70)
                tagData = data[data["Tag"]==tag]
                N_thisTag = len(tagData)
                print("# Timesteps for this machine + tag combo: {:d}".format(N_thisTag))

                #get file paths associated with this tag from batchFile
                gFileNames = tagData['GEQDSK'].values
                gFilePaths = machInDir + gFileNames
                CADfiles = machInDir + tagData['CAD'].values
                PFCfiles = machInDir + tagData['PFC'].values
                inputFiles = machInDir + tagData['Input'].values
                runList = [x.split(":") for x in tagData['Output'].values]

                #refresh all subclasses
                self.ENG.refreshSubclasses()

                #read input file 0
                inputData = self.ENG.loadInputs(inFile=inputFiles[0])

                #build the HEAT tree for this tag
                self.prepareDirectories(mach,tag)

                #read GEQDSK and load into MHD object
                self.loadMHD(mach, machInDir, gFileNames)

                #read CAD and initialize CAD objects
                #note: current version of HEAT only supports single CAD file
                #per tag
                self.loadCAD(CADfiles[0])

                #read PFC file and initialize PFC objects
                #note: current version of HEAT only supports single CAD file
                #per tag
                self.loadPFCs(PFCfiles[0])

                #note that we load HF settings (optical, gyro, rad) dynamically
                #from input file in the self.ENG.runHEAT loop

                #run HEAT
                #note: current version of HEAT only supports single runList
                #per tag
                self.runHEAT(inputFiles, runList[0])

                print("Completed all HEAT runs\n")
                log.info("Completed all HEAT runs\n")



        return

    def runHEAT(self,inputFiles,runList):
        """
        runs HEAT engine.  Steps through time solving for vars in runList
        """
        t0 = time.time()
        #if user supplied multiple input files in TUI, they are parsed at each timestep
        #note that if running openfoam, only the last timestep's input file
        #will be used for the openFOAM settings.
        self.ENG.inputFileList = inputFiles
        self.ENG.runHEAT(runList)
        #run openFOAM
        if 'T' in runList:
            self.loadOF()
            #run openFOAM analysis
            self.ENG.runOpenFOAM()
        print("Total time: {:f}".format(time.time() - t0))
        return

    def prepareDirectories(self,mach,tag, clobber='y'):
        """
        build HEAT tree for mach + tag combo
        """
        if dataPath[-1]!='/':
            self.shotPath = dataPath + '/' + mach +"_{:06d}".format(self.ENG.MHD.shot)+"_"+tag+"/"
        else:
            self.shotPath = dataPath + mach +"_{:06d}".format(self.ENG.MHD.shot)+"_"+tag+"/"

        self.ENG.MHD.shotPath = self.shotPath

        #make tree branch for this shot
        tools.makeDir(self.shotPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)

        return


    def loadMHD(self, mach, tmpDir, gFiles):
        """
        loads GEQDSK file into HEAT tree and MHD object
        """
        #2D plasmas for now
        self.ENG.MHD.plasma3Dmask = 0

        if self.ENG.MHD.plasma3Dmask == 0:
            #initialize MHD
            self.ENG.MHD.tmpDir = tmpDir
            self.ENG.MHD.tree = 'EFIT02'
            self.ENG.MHD.getGEQDSK(machine=mach,gFileList=gFiles)
            self.ENG.MHD.makeEFITobjects()
            self.ENG.MHD.psiSepLimiter = None
            self.ENG.MHD.setTypes()
            self.ENG.MHD.nTrace = int(self.ENG.MHD.traceLength / self.ENG.MHD.dpinit)
        else:
            print('3D plasmas not yet available in HEAT!')
            sys.exit()
        return

    def loadCAD(self, STPfile):
        """
        loads CAD files into CAD object
        """
        self.ENG.CAD.rootDir = rootDir #set HEAT rootDir
        self.ENG.getCADfromTUI(STPfile)
        return

    def loadPFCs(self, PFCfile):
        """
        loads PFC file into PFC objects
        """
        self.ENG.readPFCfile(PFCfile)
        self.ENG.getPFCinputs(defaultMask=False)
        return

    def loadGYRO(self):
        """
        loads gyro orbit settings
        """
        return

    def loadOF(self):
        """
        Loads OF parameters
        """
        self.ENG.loadOF(
                    self.ENG.OF.OFtMin,
                    self.ENG.OF.OFtMax,
                    self.ENG.OF.meshMinLevel,
                    self.ENG.OF.meshMaxLevel,
                    self.ENG.OF.STLscale,
                    OFbashrc,
                    self.ENG.OF.deltaT,
                    self.ENG.OF.writeDeltaT,
                    self.ENG.OF.material
                    )
        return

    def saveBatchFile(self, path=None):
        """
        Saves a batchFile to path, otherwise $HOME
        """
        if path is None:
            file = '~/batchFile.dat'
        else:
            file = path

        print("Saving batchFile template to "+path)

        text = """
#HEAT batchFile
#For use when running HEAT in terminal / batch mode.  Each line is a new entry.
#
# The fist line of every batchFile should be (uncommented):
# MachFlag, Tag, GEQDSK, CAD, PFC, Input, Output
#
#===Column variables are defined as follows
# MachFlag: machine specific flag.
#           can be 'd3d','nstx','st40','step','sparc','west','kstar'
#
# Tag:  user specified tag to label the simulation by.  Tags represent
#       independent HEAT runs.  For time varying discharges with multiple
#       GEQDSK files, tag should be repeated on multiple lines with the GEQDSK
#       for each timestep in each line.
#
# GEQDSK:  magnetic equilibrium file (ie EFIT) in GEQDSK format
#          naming convention is g<shot>.<timestep> where <shot> is the integer
#          shot number (6 digits) and timestep is the timestep in ms (5 digits).
#          For example, shot 204118 timestep 50ms would be g204118.00050
#
# CAD: CAD file for the tag.  Note that HEAT will use the first CAD file provided
#      in for each tag.  Subsequent lines in that tag are ignored.  In other words,
#      there can only be one CAD file per tag.
#
# PFC: PFC file for the tag.  Note that HEAT will use the first PFC file provided
#      in for each tag.  Subsequent lines in that tag are ignored.  In other words,
#      there can only be one PFC file per tag.
#
# INPUT: Input file for the tag.  Input files can be time varying, but only the
#        HF Variables will be read at each timestep.
#
# Output: Defines what output HEAT should calculate.  Options are:
#         -hfOpt   optical heat flux point cloud
#         -hfGyro  gyro orbit heat flux point cloud
#         -B       magnetic field glyph cloud
#         -psiN    normalized poloidal flux point cloud
#         -pwrDir  powerDir point cloud
#         -bdotn   bdotn point cloud
#         -norm    normal vector glyph cloud
#         -T       temperature
#
#       for multiple outputs, separate options with : (ie hfOpt:psi:T).  Note
#       that HEAT will use the first options list provided for each tag.
#       Subsequent lines in that tag are ignored.  In other words, there can
#       only be one set of options per tag.
#
#
# Once you have a batchFile, you need to save all input files in the following
# directory structure, where <path> is wherever the batchFile is:
# <path>/batchFile.dat
# <path>/MachFlag/GEQDSK
# <path>/MachFlag/CAD
# <path>/MachFlag/PFC
# <path>/MachFlag/Input
#
#  Example line for an NSTX-U run:
#MachFlag, Tag, GEQDSK, CAD, PFC, Input, Output
#nstx,run1, g204118.00004, IBDH_2tiles.step, PFCs_run1.csv, NSTXU_input.csv, B:hfOpt
#
# And the directory structure would look like this
# <path>/batchFile.dat
# <path>/nstx/g204118.00004
# <path>/nstx/IBDH_2tiles.step
# <path>/nstx/PFCs_run1.csv
# <path>/nstx/NSTXU_input.csv
#
#
#
#
MachFlag, Tag, GEQDSK, CAD, PFC, Input, Output
            """
        with open(file,'w') as f:
            f.write(text)




        return

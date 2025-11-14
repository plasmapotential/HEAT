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
try:
    logFile = os.environ["logFile"]
    rootDir = os.environ["rootDir"]
    dataPath = os.environ["dataPath"]
    OFbashrc = os.environ["OFbashrc"]
    FreeCADPath = os.environ["FreeCADPath"]
except:
    print("Could not properly set up environment.  Proceed with caution!")

try:
    chmod = int(os.environ["HEATchmod"], 8) #convert from base 8
except:
    chmod = 0o774

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
import engineClass

class TUI():
    def __init__(self):
        """
        intialize terminal user interface (TUI) object
        """
        self.ENG = engineClass.engineObj(logFile, rootDir, dataPath, OFbashrc, chmod, UID, GID)
        self.ENG.NCPUs = multiprocessing.cpu_count() - 2 #reserve 2 cores for overhead
        self.chmod = chmod
        self.GID = GID
        self.UID = UID
        #if data directory doesn't exist, create it
        tools.makeDir(dataPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        #set up number formats for file IO and printing
        self.setupNumberFormats(self.ENG.tsSigFigs, self.ENG.shotSigFigs)
        return

    def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
        """
        sets up pythonic string number formats for shot and timesteps
        """
        self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
        self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
        return

    def simulationSchedule(self, batchFile):
        """
        determines what simulations need to be run and their correct grouping
        via a batchFile

        The fist line of every batchFile should be:
        MachFlag, Tag, EQ, CAD, PFC, Input, Output

        Column variables are defined as follows

        :MachFlag: machine specific flag.
          can be 'd3d','nstx','st40','step','sparc', 'arc', 'cmod, 'west','kstar','aug','tcv','other'

        :Tag:  user specified tag to label the simulation by.  Tags represent
          independent HEAT runs.  For time varying discharges with multiple
          GEQDSK files, tag should be repeated on multiple lines with the 
          TimeStep column changing.

        :Shot:  pulse number to use for saving HEAT output (MachFlag_Shot).  Default is to
          allow for 6 sig figs of shot numbers.  SigFigs can be changed in engineClass
          object initialization function call.
        
        :TimeStep: Timestep that the equilibrium defined in "GEQDSK" column corresponds
          to.  GEQDSKs can be named in multiple formats, but the TimeStep in 
          this row is what is used by HEAT.  For time varying discharges, tag 
          should be repeated on multiple lines with the TimeStep column changing.
          Timestep units are [s].  Default is to allow for 9 SigFigs after the
          radix, which corresponds to nanosecond resolution (defined in engineClass
          initialization function call)
        
        :EQ:  magnetic equilibrium file 
          can be in GEQDSK format or in a JSON following the IMAS schema
          Psi should be in units of Wb/rad (divide by 2pi) and the Bt0, Fpol, Psi
          and Ip values should reflect COCOS and match machine helicity.  This column
          can also be labeled GEQDSK for backwards compatability
          
        :CAD: CAD file for the tag.  Note that HEAT will use the first CAD file provided
          in for each tag.  Subsequent lines in that tag are ignored.  In other words,
          there can only be one CAD file per tag.  Can be STEP, IGES, or FCStd formats.

        :PFC: PFC file for the tag.  Note that HEAT will use the first PFC file provided
          in for each tag.  Subsequent lines in that tag are ignored.  In other words,
          there can only be one PFC file per tag.

        :INPUT: Input file for the tag.  Input files can be time varying, but only the
          HF Variables will be read at each timestep.

        :Output: Defines what output HEAT should calculate.  Options are:
                 -hfOpt   optical heat flux 
                 -hfGyro  gyro orbit heat flux 
                 -hfFil   filament heat flux 
                 -hfRE    runaway electron heat flux
                 -hfRad   radiated power heat 
                 -B       magnetic field glyphs
                 -psiN    normalized poloidal flux 
                 -pwrDir  powerDir point 
                 -bdotn   bdotn point
                 -norm    normal vector glyphs
                 -T       temperature calculation using openFOAM
                 -elmer   runs an Elmer FEM simulation
                 -Btrace  Bfield tracing from a trace csv file

          for multiple outputs, separate options with : (ie hfOpt:psi:T).  Note
          that HEAT will use the first options list provided for each tag.
          Subsequent lines in that tag are ignored.  In other words, there can
          only be one set of options per tag.


        Once you have a batchFile, you need to save all input files in the following
        directory structure, where <path> is wherever the batchFile is:

        - <path>/batchFile.dat
        - <path>/MachFlag/GEQDSK
        - <path>/MachFlag/CAD
        - <path>/MachFlag/PFC
        - <path>/MachFlag/Input


        Example line for an NSTX-U run for pulse 204118 timestep 4ms:

        MachFlag, Tag, Shot, TimeStep, EQ, CAD, PFC, Input, Output

        nstx,run1, 204118, 0.004, geqdsk.00004, IBDH_2tiles.step, PFCs_run1.csv, NSTXU_input.csv, B:hfOpt
        
        And the directory structure would look like this

        - <path>/batchFile.dat
        - <path>/nstx/g204118.00004
        - <path>/nstx/IBDH_2tiles.step
        - <path>/nstx/PFCs_run1.csv
        - <path>/nstx/NSTXU_input.csv

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
                log.info('\n')
                log.info("-"*70)
                log.info(" "*20 + "Machine: "+mach+"   Tag: "+tag)
                log.info("-"*70)

                tagData = data[data["Tag"]==tag]
                N_thisTag = len(tagData)
                print("# Timesteps for this machine + tag combo: {:d}".format(N_thisTag))

                #get file paths associated with this tag from batchFile
                try:

                    shots = tagData['Shot'].values #only 1 shot per tag allowed
                    timesteps = tagData['TimeStep'].values

                    #this conditional allows for various EQ formats
                    if 'GEQDSK' in tagData.keys():
                        eqFileNames = tagData['GEQDSK'].values
                    else:
                        eqFileNames = tagData['EQ'].values
                    
                    eqFilePaths = machInDir + eqFileNames

                    #if user is bringing their own meshes then the CAD file can be set to None
                    #NEED TO FIX THIS SO THAT THERE CAN BE NONE AND strings in same batchFile
                    if True in pd.isna(tagData['CAD'].values):
                        print("No CAD file provided, assuming you are providing STLs")
                        log.info("No CAD file provided, assuming you are providing STLs")
                        CADfiles = tagData['CAD'].values
                    #normal cases, user supplies CAD file of some kind
                    else:
                        CADfiles = machInDir + tagData['CAD'].values
                    
                    #CADfiles = machInDir + tagData['CAD'].values
                    PFCfiles = machInDir + tagData['PFC'].values
                    inputFiles = machInDir + tagData['Input'].values
                    runList = [x.strip().split(":") for x in tagData['Output'].values]
                    runList = np.unique([x for y in runList for x in y])
                except Exception as e:
                    print("\n\nSomething is wrong with your batchFile!  Error Trace:\n")
                    print(e)
                    sys.exit()

                #refresh all subclasses
                self.ENG.refreshSubclasses()

                #read input file 0
                inputData = self.ENG.loadInputs(inFile=inputFiles[0])

                #build the HEAT tree for this tag
                self.prepareDirectories(mach,tag)

                #load filament data 
                self.loadFilaments(runList, machInDir)

                #build timesteps
                self.loadTimeSteps(timesteps, shots[0], tag, self.ENG.FIL.tsFil)

                #read GEQDSK and load into MHD object
                self.loadMHD(machInDir, eqFileNames, timesteps)

                #if we are only running a Bfield trace, no need to load CAD
                if 'Btrace' in runList and len(runList) == 1:
                    pass
                else:
                    #read CAD and initialize CAD objects
                    #note: current version of HEAT only supports single CAD file
                    #per tag
                    self.loadCAD(CADfiles[0], machInDir)

                    #read PFC file and initialize PFC objects
                    #note: current version of HEAT only supports single CAD file
                    #per tag
                    self.loadPFCs(PFCfiles[0])

                #note that we load HF settings (optical, gyro, rad) dynamically
                #from input file in the self.ENG.runHEAT loop

                #set up output file stream
                self.ENG.getIOInputs()

                #run HEAT
                #note: current version of HEAT only supports single runList
                #per tag
                self.runHEAT(inputFiles, runList)

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
        #run thermal analysis
        if 'T' in runList:
            self.loadOF()
            if self.ENG.IO.csvMask == True:
                #run openFOAM analysis
                self.ENG.runOpenFOAM()
            else:
                print("Error!  Cannot run openFOAM unless you save CSV files!")
                print("Please turn on csv output file switch and re-run")
                log.info("Error!  Cannot run openFOAM unless you save CSV files!")
                log.info("Please turn on csv output file switch and re-run")

        if 'elmer' in runList:
            self.loadElmer()
            self.ENG.runElmerFEM()

        print("HEAT run completed...")
        log.info("HEAT run completed...")
        print("Total time: {:f}".format(time.time() - t0))
        log.info("Total time: {:f}".format(time.time() - t0))
        return

    def prepareDirectories(self,mach,tag):
        """
        build HEAT tree for mach + tag combo
        """
        if dataPath[-1]!='/':
            self.shotPath = dataPath + '/' + mach +"_"+self.shotFmt.format(self.ENG.MHD.shot)+"_"+tag+"/"
        else:
            self.shotPath = dataPath + mach +"_"+self.shotFmt.format(self.ENG.MHD.shot)+"_"+tag+"/"

        self.ENG.MHD.shotPath = self.shotPath

        #make tree branch for this shot
        tools.makeDir(self.shotPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)

        #make unique logfile for this tag
        logFile = self.shotPath + 'HEATlog.txt'
        self.ENG.logFile = logFile
        from logConfig import setup_logging
        log.info("Changing log file to new path for batchMode:")
        log.info(logFile)
        setup_logging(logfile_path=logFile)

        return

    def loadFilaments(self, runList, path):
        """
        loads a filament file and build timestep array for filament tracing

        runList is HEAT runlist of unique values
        path is path to filament file
        """
        if 'hfFil' in runList:
            self.ENG.FIL.readFilamentFile(path)
            self.ENG.FIL.setupFilamentTime()
        elif 'hfRE' in runList:
            self.ENG.RE.readREFile(path)
            self.ENG.RE.setupRETime()
            self.ENG.FIL.tsFil = self.ENG.RE.tsFil
        else:
            self.ENG.FIL.tsFil = None
        return

    def loadTimeSteps(self, ts, shot, tag, tsFil=None):
        """
        ts is timesteps in batchfile
        tsFil is 2D list of timesteps for each filament (can be finer resolution than batchfile)

        loads timesteps data from batchFile and filament file.  
        """
        if tsFil is not None:
            for row in tsFil:
                ts = np.append(ts, row)

        ts = np.unique(np.sort(ts))
        self.ENG.setupTime(ts, shot, tag)

        return


    def loadMHD(self, tmpDir, eqFiles, ts):
        """
        loads EQ files into HEAT tree and MHD object
        """
        #initialize MHD
        self.ENG.MHD.tmpDir = tmpDir
        self.ENG.MHD.EQmode = self.ENG.MHD.determineEQFiletype(eqFiles[0])
        self.ENG.MHD.getGEQDSK(ts, eqFiles)
        self.ENG.MHD.makeEFITobjects()
        self.ENG.MHD.psiSepLimiter = None
        self.ENG.MHD.setTypes()
        self.ENG.MHD.nTrace = int(self.ENG.MHD.traceLength / self.ENG.MHD.dpinit)
        return

    def loadCAD(self, CADfile, machInDir):
        """
        loads CAD files into CAD object
        """
        self.ENG.CAD.rootDir = rootDir #set HEAT rootDir
        self.ENG.CAD.machInDir = machInDir #HEATrun directory
        self.ENG.getCADinputs()
        self.ENG.getCADfromTUI(CADfile)
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
    
    def loadElmer(self):
        """
        Loads Elmer FEM parameters into engine object
        """
        print("Elmer FEM")
        self.ENG.FEM.readElmerFile()
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
# MachFlag, Tag, Shot, TimeStep, GEQDSK, CAD, PFC, Input, Output
#
#===Column variables are defined as follows
# MachFlag: machine specific flag.
#           can be 'd3d','nstx','st40','step','sparc', 'arc', 'cmod, 'west','kstar','aug','tcv','other'
#
# Tag:  user specified tag to label the simulation by.  Tags represent
#       independent HEAT runs.  For time varying discharges with multiple
#       GEQDSK files, tag should be repeated on multiple lines with the 
#       TimeStep column changing.
#
# Shot:  pulse number to use for saving HEAT output (MachFlag_Shot).  Default is to
#        allow for 6 sig figs of shot numbers.  SigFigs can be changed in engineClass
#        object initialization function call.
#
# TimeStep: Timestep that the equilibrium defined in "GEQDSK" column corresponds
#           to.  GEQDSKs can be named in multiple formats, but the TimeStep in 
#           this row is what is used by HEAT.  For time varying discharges, tag 
#           should be repeated on multiple lines with the TimeStep column changing.
#           Timestep units are [s].  Default is to allow for 9 SigFigs after the
#           radix, which corresponds to nanosecond resolution (defined in engineClass
#           initialization function call)
#
# EQ:  magnetic equilibrium file 
#          can be in GEQDSK format or in a JSON following the IMAS schema
#          Psi should be in units of Wb/rad (divide by 2pi) and the Bt0, Fpol, Psi
#          and Ip values should reflect COCOS and match machine helicity
#
# CAD: CAD file for the tag.  Note that HEAT will use the first CAD file provided
#      in for each tag.  Subsequent lines in that tag are ignored.  In other words,
#      there can only be one CAD file per tag.  Can be STEP, IGES, or FCStd formats.
#      User can also supply a mesh (.STL) file and set this value to None.  In that
#      case, the mesh file should be listed in the PFC CSV file for each ROI or
#      intersection part
#
# PFC: PFC file for the tag.  Note that HEAT will use the first PFC file provided
#      in for each tag.  Subsequent lines in that tag are ignored.  In other words,
#      there can only be one PFC file per tag.
#
# INPUT: Input file for the tag.  Input files can be time varying, but only the
#        HF Variables will be read at each timestep.
#
# Output: Defines what output HEAT should calculate.  Options are:
#         -hfOpt   optical heat flux 
#         -hfGyro  gyro orbit heat flux 
#         -hfFil   filament heat flux 
#         -hfRE    runaway electron heat flux
#         -hfRad   radiated power heat 
#         -B       magnetic field glyphs
#         -psiN    normalized poloidal flux 
#         -pwrDir  powerDir point 
#         -bdotn   bdotn point
#         -norm    normal vector glyphs
#         -T       temperature calculation using openFOAM
#         -elmer   runs an Elmer FEM simulation
#
#       for multiple outputs, separate options with : (ie hfOpt:psi:T).  Note
#       that HEAT will use the first options list provided for each tag.
#       Subsequent lines in that tag are ignored.  In other words, there can
#       only be one set of options per tag.
#
#
# Once you have a batchFile, you need to save all input files in the following
# directory structure, where <path> is wherever the batchFile is and <MachFlag>
# is the machine name (ie nstx):
# <path>/batchFile.dat
# <path>/<MachFlag>/GEQDSK
# <path>/<MachFlag>/CAD
# <path>/<MachFlag>/PFC
# <path>/<MachFlag>/Input
#
#  Example line for an NSTX-U run for pulse 204118 timestep 4ms:
#MachFlag, Tag, Shot, TimeStep, EQ, CAD, PFC, Input, Output
#nstx,run1, 204118, 0.004, geqdsk.00004, IBDH_2tiles.step, PFCs_run1.csv, NSTXU_input.csv, B:hfOpt
#  Example line for an NSTX-U run for pulse 204118 using JSON and user supplied mesh
#nstx,run2, 204118, 0.004, nstx_eq_204118.json, None, PFCs_run2.csv, NSTXU_input.csv, B:hfOpt
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
MachFlag, Tag, Shot, TimeStep, EQ, CAD, PFC, Input, Output
            """
        with open(file,'w') as f:
            f.write(text)




        return

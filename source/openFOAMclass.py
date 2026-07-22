#OpenFOAMclass.py
#Description:   OpenFOAM Module
#Engineer:      T Looby
#Date:          20200327

import PyFoam
import sys
import os
import shutil
import toolsClass
import subprocess
import numpy as np
import logging
import pandas as pd
import psutil
log = logging.getLogger(__name__)
tools = toolsClass.tools()
import logConfig

class OpenFOAM():
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
        #self.buildheatFoam()
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
          
        OpenFOAM Variables:
        -------------------

        :OFtMin: minimum timestep of openFOAM simulation [s]
        :OFtMax: maximum timestep of openFOAM simulation [s]
        :deltaT: timestep size for simulation [s]
        :writeDeltaT: timestep size for writing simulation [s].  Currently must be set to 
          same value as deltaT.
        :STLscale: scales points in STL mesh up/down by scalar value [float].  Can be used for
          unit conversions on meshes.
        :meshMinLevel: minimum number of refinement iterations to be completed when snapping
          a volume mesh to a surface mesh [integrer].  Increasing this value creates a finer
          volume mesh but takes longer.  Cannot be larger than meshMaxLevel.
        :meshMaxLevel: maximum number of refinement iterations to be completed when snapping
          a volume mesh to a surface mesh [integrer].  Increasing this value creates a finer
          volume mesh but takes longer.  Cannot be smaller than meshMinLevel.        

                
        """
        self.allowed_vars = [
                             'OFtMin',
                             'OFtMax',
                             'deltaT',
                             'writeDeltaT',
                             'STLscale',
                             'meshMinLevel',
                             'meshMaxLevel',
                             'material'
                             ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        integers = [
                    'meshMinLevel',
                    'meshMaxLevel',
                    ]
        floats = [
                  'OFtMin',
                  'OFtMax',
                  'deltaT',
                  'writeDeltaT',
                  'STLscale',
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

    def buildheatFoam(self):
        """
        build heatFoam if we are in appImage mode.  This is useful if you
        including the binary heatFoam + libs doesn't work when constructing
        the appImage.  This function compiles heatFoam inside the container.
        """
        try:
            print("Compiling heatFoam solver for openFOAM inside container")
            log.info("Compiling heatFoam solver for openFOAM inside container")
            #Copy the current environment
            current_env = os.environ.copy()
            #run blockMesh, snappyHexMesh, topoSet, createPatch
            homeDir = os.environ["homeDir"]
            OFversion = os.environ["OFversion"]
            heatFoamDir = homeDir+'/source/openfoam/OpenFOAM-v'+OFversion+'/applications/solvers/custom/heatFoam'
            buildFile = heatFoamDir + '/buildHEATfoam'
            with open(buildFile, 'w') as f:
                f.write('#!/bin/bash\n')
                f.write(self.cmdSourceOF + '\n')
                f.write('wclean\n')
                f.write('wmake\n')
            subprocess.run([buildFile], env=current_env, cwd=heatFoamDir)
            print("Compiled heatFoam")
            log.info("Compiled heatFoam")
        except:
            print("Did not build heatFoam")
            log.info("Did not build heatFoam")
        return

    def createDictionaries(self, templateDir, partDir, templateVarFile, STLfile):
        """
        Creates openfoam dictionaries from a template dictionary for a specific
        user specified HEAT configuration.

        Input Arguments:
        templateDir     location of template dictionary files
        partDir         location where we will run openFOAM
        STLfile         name of STLfile that we are using.  Should be located in STLpath (from CAD module)
        """

        if templateDir[-1] != '/':
            templateDir += '/'
        if partDir[-1] != '/':
            partDir += '/'

        STLlayerName = STLfile + '_firstSolid'

        #list of files to copy
        templateList = [
                        'blockMeshDict',
                        'controlDict',
                        #'surfaceFeatureExtractDict', #use for tricky meshes
                        'snappyHexMeshDict',
                        'topoSetDict',
                        'meshQualityDict',
                        'createPatchDict',
                        'decomposeParDict'
                        #, 'probes'
                        ]
        #generate dictionaries from template files
        for file in templateList:
            inFile = templateDir + file + '.template'
            outFile = partDir + 'system/' + file
            #cmd = ('pyFoamFromTemplate.py'
            #       ' --template-file='+inFile+
            #       ' --output-file='+outFile+
            #       ' --values-dictionary='+templateVarFile)
            #build pyFoam command
            args=[]
            #args.append('pyFoamFromTemplate.py')
            args.append('--template-file='+inFile)
            args.append('--output-file='+outFile)
            args.append('--values-dictionary='+templateVarFile)

            from PyFoam.Applications.FromTemplate import FromTemplate
            FromTemplate(args)
        return

    def writeOFtemplateVarFile(self, file, STLpart):
        """
        writes a file that pyfoam uses to populate the OF dictionary templates
        with real variables.

        file        path to the file we will write
        STLpart     name of STL that openfoam will create 3D mesh from
        rootDir     path to HEAT root
        """
        #determine number of logical cpus and leave 1 for overhead
        if psutil.cpu_count(logical=True) < 2:
            self.NCPU = 1
        else:
            self.NCPU = psutil.cpu_count(logical=True) - 1
        with open(file, 'w') as f:
            f.write('xMin {:f};\n'.format(self.xMin))
            f.write('xMax {:f};\n'.format(self.xMax))
            f.write('yMin {:f};\n'.format(self.yMin))
            f.write('yMax {:f};\n'.format(self.yMax))
            f.write('zMin {:f};\n'.format(self.zMin))
            f.write('zMax {:f};\n'.format(self.zMax))
            f.write('tMin {:f};\n'.format(self.OFtMin))
            f.write('tMax {:f};\n'.format(self.OFtMax))
            f.write('deltaT {:f};\n'.format(self.deltaT))
            f.write('writeDeltaT {:f};\n'.format(self.writeDeltaT))
            if hasattr(self, 'xProbe'):
                f.write('xProbe {:f};\n'.format(self.xProbe))
                f.write('yProbe {:f};\n'.format(self.yProbe))
                f.write('zProbe {:f};\n'.format(self.zProbe))
            f.write('STLscale {:f};\n'.format(self.STLscale))
            f.write('meshMinLevel {:d};\n'.format(self.meshMinLevel))
            f.write('meshMaxLevel {:d};\n'.format(self.meshMaxLevel))
            f.write('xMid {:f};\n'.format(self.xMid))
            f.write('yMid {:f};\n'.format(self.yMid))
            f.write('zMid {:f};\n'.format(self.zMid))
            f.write('STLfileName '+STLpart+';\n')
            f.write('STLlayerName '+STLpart+'_firstSolid;\n')
            f.write('NCPU '+str(self.NCPU)+';\n')
        return

    def writeShellScript(self, parallel=False):
        """
        Writes shell script that runs OF from the template file.  This file will
            write to the HEATlog.txt file for viewing in the GUI

        logFile is absolute location of HEAT log file
        if parallel is set to True, will mesh across multiple cores
            (this is buggy and sometimes fails)

        """
        #get current log file path
        logFile = logConfig.get_current_log_path()

        #check if we are in appImage mode to get correct bash location
        shebang = '#!/bin/bash'

        #get part directory
        if self.partDir[-1] != '/':
            self.partDir += '/'

        #determine number of logical cpus and leave 1 for overhead
        #now done in template file creation method
        #if psutil.cpu_count(logical=True) < 2:
        #    self.NCPU = 1
        #else:
        #    self.NCPU = psutil.cpu_count(logical=True) - 1


        #Write 3D meshing script
        file = self.partDir + self.cmd3Dmesh
        with open(file, 'w') as f:
            f.write(shebang + '\n')
            #source OF bashrc
            f.write(self.cmdSourceOF + '\n')

            #for tricky meshes, extract features
            #f.write('surfaceFeatureExtract | tee -a ' + logFile + '\n')

            f.write('blockMesh | tee -a ' + logFile + '\n')
            #f.write('blockMesh > ' + logFile + '\n')
            #single core meshing
            if parallel==False:
                f.write('snappyHexMesh -overwrite | tee -a '+logFile+ '\n')
                #f.write('snappyHexMesh -overwrite > '+logFile+ '\n')
            #parallel meshing
            else:
                #Run snappyHexMesh across multiple processors
                #this feature should be added in one day...sigh...for now its a placeholder
                #you need to download libscotch-6.0 from ubuntu repo,
                #then build scotchDecomp from source from src/parallel/decompose/scotchDecomp
                f.write('decomposePar | tee -a ' + logFile + '\n')
                f.write('mpirun --use-hwthread-cpus -np '+str(self.NCPU)+' snappyHexMesh -parallel -overwrite | tee -a '+logFile+ '\n')
                f.write('reconstructParMesh -mergeTol 1e-6 -latestTime -constant | tee -a '+logFile+ '\n')

        try:
            os.chmod(file, self.chmod)
        except OSError:
            pass

        #Write thermal analysis script
        file = self.partDir + self.cmdThermal
        foamFile = self.partDir + '*.foam'
        with open(file, 'w') as f:
            f.write(shebang + '\n')
            #source OF bashrc
            f.write(self.cmdSourceOF + '\n')
            f.write('topoSet | tee -a '+logFile+ '\n')
            f.write('createPatch -overwrite | tee -a '+logFile+ '\n')
            f.write('heatFoam | tee -a '+logFile+ '\n')
            #f.write('topoSet > '+logFile+ '\n')
            #f.write('createPatch -overwrite > '+logFile+ '\n')
            #f.write('heatFoam > '+logFile+ '\n')
            #f.write('paraFoam -touchAll\n')
            #f.write('touch '+ self.partName +'.foam\n')
        try:
            os.chmod(file, self.chmod)
        except OSError:
            pass

        #Write temperature probe script
        file = self.partDir + self.cmdTprobe
        with open(file, 'w') as f:
            f.write(shebang + '\n')
            #source OF bashrc
            f.write(self.cmdSourceOF + '\n')
            f.write('topoSet | tee -a '+logFile+ '\n')
            f.write('createPatch -overwrite | tee -a '+logFile+ '\n')
            f.write('postProcess -func "probes" | tee -a '+logFile+ '\n')
        try:
            os.chmod(file, self.chmod)
        except OSError:
            pass
        return


    def generate3Dmesh(self, part, overWriteMask=False):
        """
        Run OpenFOAM blockmesh and snappyhexmesh to generate a 3D mesh.  Then
        run topoSet and createPatch to create patches

        part is the part name for the part we are meshing.  Before this command
        is run there should be an STL file in the constant/triSurface directory
        of the openfoam case
        """
        print('Generating 3D Mesh in OpenFOAM')
        log.info('Generating 3D Mesh in OpenFOAM')

        #Check to see if 3D mesh exists for these settings
        minLev = int(self.meshMinLevel)
        maxLev = int(self.meshMaxLevel)

        if self.meshDir[-1] == '/':
            file = self.meshDir + part + '_{:d}_{:d}'.format(minLev,maxLev)
        else:
            file = self.meshDir + '/' + part + '_{:d}_{:d}'.format(minLev,maxLev)

        if self.partDir[-1] != '/':
            self.partDir += '/'

        #if 3D mesh doesn't exist, create it.  Otherwise copy it to OF.partDir/constant/polyMesh
        newFile = self.partDir+'constant/polyMesh'
        if os.path.isdir(file) == True: #polyMesh already exists
            #overWriteMask is True if STP file was modified since last use
            if overWriteMask == False:
                try:
                    tools.copytree(file, newFile)
                    print('Copied existing 3D mesh')
                except:
                    print('\nProblem copying 3D mesh...possibly doesnt exist...creating')
                    print('tried copying: '+file)
                    print('to this location: '+newFile+'\n')
            else:
                shutil.rmtree(file)
                self.createNewMesh(newFile,file)


        else:
            print("Creating new polymesh as none existed")
            self.createNewMesh(newFile,file)
        return

    def createNewMesh(self, newFile, file):
        """
        creates new FVM mesh and copies into a HEAT tree
        """
        #Copy the current environment
        current_env = os.environ.copy()
        #point to correct path for bash
        bashExec = '/bin/bash'
        #run blockMesh, snappyHexMesh
        meshCMD = self.partDir+self.cmd3Dmesh
        print(self.cmd3Dmesh)
        print()
        try:
            p = subprocess.run([meshCMD], env=current_env, cwd=self.partDir, shell=True, executable=bashExec)
            retcode = p.returncode
            if retcode < 0:
                print("OF mesh child was terminated by signal", -retcode, file=sys.stderr)
            else:
                print("OF mesh child returned", retcode, file=sys.stderr)
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)

        #Now copy this mesh to the 3D meshes folder for future use
        #tools.makeDir(newFile, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        tools.copytree(newFile, file)
        #set tree permissions
        tools.recursivePermissions(self.meshDir, self.UID, self.GID, self.chmod)
# This method left here for reference:
#            from PyFoam.Execution.BasicRunner import BasicRunner
#            solvers = ['blockMesh','snappyHexMesh','topoSet','createPatch']
#            dir = self.OF.partDir
#            case = self.OF.partDir
#            for solver in solvers:
#                run = BasicRunner(argv=[solver,dir,case])
#                run.start()
        return


    def runThermalAnalysis(self):
        """
        runs an openfoam thermal analysis in the specified location
        """
        print('Running Thermal Analysis')
        log.info('Running Thermal Analysis')
        print('See HEAT LogFile for Status')
        log.info('See HEAT LogFile  for Status')
        #Copy the current environment
        current_env = os.environ.copy()
        #point to correct path for bash
        bashExec = '/bin/bash'

        #run topoSet, createPatch, heatFoam, paraFoam -touchAll
        thermalCMD = self.partDir+self.cmdThermal
        try:
            p = subprocess.run([thermalCMD], env=current_env, cwd=self.partDir, shell=True, executable=bashExec)
            retcode = p.returncode
            if retcode < 0:
                print("thermal analysis child was terminated by signal", -retcode, file=sys.stderr)
            else:
                print("thermal analysis child returned", retcode, file=sys.stderr)

        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)

        from pathlib import Path
        Path(self.partDir+self.partName+'.foam').touch()

        #set tree permissions
        tools.recursivePermissions(self.partDir, self.UID, self.GID, self.chmod)
        return

    def runTprobe(self,x,y,z, partDir):
        """
        runs openfoam postprocessing function to find temperature at specified
        x,y,z probe location
        """
        print('Simulating Temperature Probe')
        log.info('Simulating Temperature Probe')

        from subprocess import run
        #Copy the current environment
        current_env = os.environ.copy()

        #write probe file to partDir
        file = partDir+'system/probes'
        with open(file,'w') as f:
            f.write('type probes;\n')
            f.write('libs ("libsampling.so");\n')
            f.write('fields (T);\n')
            f.write('probeLocations ( ({:.6f} {:.6f} {:.6f}) );\n'.format(x,y,z))

        #run topoSet, createPatch, postProcess -func "probes"
        TprobeCMD = partDir+self.cmdTprobe
        try:
            p = run([TprobeCMD], env=current_env, cwd=self.partDir, shell=True, executable='/bin/bash')
            retcode = p.returncode
            if retcode < 0:
                print("thermal probe child was terminated by signal", -retcode, file=sys.stderr)
            else:
                print("thermal probe child returned", retcode, file=sys.stderr)
        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)

        print("Completed Tprobe run")
        log.info("Completed Tprobe run")
        return

    def getMinMaxData(self, file):
        """
        gets min / max data that was assembled in openfoam postprocessing function

        file is location of openfoam output in postProcessing/fieldMinMax1/<timestep>
        returns dataframe
        """
        path, name = os.path.split(file)
        outfile = path + '/minMaxTnoTab.dat'
        if os.path.isfile(file):
            #first clean the tabs out of csv
            #Equivalent of: sed 's/\t/,/g' file > outfile
            with open(file, 'r') as fin:
                with open(outfile, 'w') as fout:
                    for line in fin:
                        fout.write(line.replace('\t',','))
            #outfile permissions
            try:
                os.chown(outfile,self.UID,self.GID)
                os.chmod(outfile,self.chmod)
            except OSError:
                pass
            #read data file
            data = pd.read_csv(outfile, header=1)
            data.columns = data.columns.str.strip()
            data = data.sort_values('field')
            data['field'] = data['field'].str.strip()


        else:
            print("minMaxTnoTab.dat does not exist!  This is probably because openFOAM failed silently.")
            log.info("minMaxTnoTab.dat does not exist!  This is probably because openFOAM failed silently.")
            raise ValueError("minMaxTnoTab.dat does not exist!  This is probably because openFOAM failed silently.")

        return data

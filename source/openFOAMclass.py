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

class OpenFOAM():
    def __init__(self, rootDir, dataPath):
        """
        rootDir is root location of python modules (where dashGUI.py lives)
        dataPath is the location where we write all output to
        """
        self.rootDir = rootDir
        tools.rootDir = self.rootDir
        self.dataPath = dataPath
        tools.dataPath = self.dataPath
        #self.buildheatFoam()
        return

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization
        """
        self.allowed_vars = ['xMin',
                             'xMax',
                             'yMin',
                             'yMax',
                             'zMin',
                             'zMax',
                             'tMin',
                             'tMax',
                             'deltaT',
                             'writeDeltaT',
                             'xProbe',
                             'yProbe',
                             'zProbe',
                             'STLscale',
                             'meshMinLevel',
                             'meshMaxLevel',
                             'xMid',
                             'yMid',
                             'zMid',
                             'STLfileName',
                             'STLlayerName',
                             ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file

        Template Variables:
        xMin            minimum X coordinate for blockMesh box
        xMax            maximum X coordinate for blockMesh box
        yMin            minimum Y coordinate for blockMesh box
        yMax            maximum Y coordinate for blockMesh box
        zMin            minimum Z coordinate for blockMesh box
        zMax            maximum Z coordinate for blockMesh box
        tMin            minimum timestep for controlDict
        tMax            maximum timestep for controlDict
        deltaT          solution time resolution for controlDict
        writeDeltaT     write time resolution for controlDict
        xProbe          x coordinate for temperature probe for postProcess
        yProbe          y coordinate for temperature probe for postProcess
        zProbe          z coordinate for temperature probe for postProcess
        STLfileName     STL file name for use with snappyHexMesh
        STLlayerName    STL layer name for snappyHexMesh
        STLscale        scale for unit conversion (1.0 = keep STL units)
        meshMinLevel    minimum level for snappyHexMesh refinement
        meshMaxLevel    maximum level for snappyHexMesh refinement
        xMid            x coordinate in center of tile of interest for snappyHexMesh
        yMid            x coordinate in center of tile of interest for snappyHexMesh
        zMid            x coordinate in center of tile of interest for snappyHexMesh
        OFbashrc        location of OpenFOAM installation bashrc file
        """
        self.xMin = float(self.xMin)
        self.xMax = float(self.xMax)
        self.yMin = float(self.yMin)
        self.yMax = float(self.yMax)
        self.zMin = float(self.zMin)
        self.zMax = float(self.zMax)
        self.tMin = float(self.tMin)
        self.tMax = float(self.tMax)
        self.deltaT = float(self.deltaT)
        self.writeDeltaT = float(self.writeDeltaT)
        self.STLscale = float(self.STLscale)
        self.meshMinLevel = int(self.meshMinLevel)
        self.meshMaxLevel = int(self.meshMaxLevel)
        self.xMid = float(self.xMid)
        self.yMid = float(self.yMid)
        self.zMid = float(self.zMid)
        if self.xProbe:
            self.xProbe = float(self.xProbe)
        if self.yProbe:
            self.yProbe = float(self.yProbe)
        if self.zProbe:
            self.zProbe = float(self.zProbe)
        return

    def buildheatFoam(self):
        """
        build heatFoam if we are in appImage mode.  This is useful if you
        including the binary heatFoam + libs doesn't work when constructing
        the appImage.  This function compiles heatFoam inside the appImage.
        """
        try:
            AppImage = os.environ["APPIMAGE"]
            inAppImage = True
            print("Compiling heatFoam solver for openFOAM inside appImage")
            log.info("Compiling heatFoam solver for openFOAM inside appImage")
            #Copy the current environment
            current_env = os.environ.copy()
            #run blockMesh, snappyHexMesh, topoSet, createPatch
            AppDir = os.environ["APPDIR"]
            heatFoamDir = AppDir+'/usr/lib/openfoam/openfoam1912/applications/solvers/custom/heatFoam'
            buildFile = heatFoamDir + '/buildHEATinAppImage'
            with open(buildFile, 'w') as f:
                f.write('#!/bin/bash\n')
                f.write(self.cmdSourceOF + '\n')
                f.write('wclean\n')
                f.write('wmake\n')
            subprocess.run([buildFile], env=current_env, cwd=heatFoamDir)
            print("Compiled heatFoam inside appImage")
            log.info("Compiled heatFoam inside appImage")
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
            NCPU = 1
        else:
            NCPU = psutil.cpu_count(logical=True) - 1
        with open(file, 'w') as f:
            for var in self.allowed_vars:
                f.write(var+' '+str(getattr(self,var))+';\n')
            f.write('STLfileName '+STLpart+';\n')
            f.write('STLlayerName '+STLpart+'_firstSolid;\n')
            f.write('NCPU '+str(NCPU)+';\n')
        return

    def writeShellScript(self,logFile, parallel=False):
        """
        Writes shell script that runs OF from the template file.  This file will
            write to the HEATlog.txt file for viewing in the GUI

        logFile is absolute location of HEAT log file
        if parallel is set to True, will mesh across multiple cores
            (this is buggy and sometimes fails)

        """
        #check if we are in appImage mode to get correct bash location
        try:
            AppImage = os.environ["APPIMAGE"]
            shebang = '#!' + os.environ["APPDIR"] + '/bin/bash'
            AppDir = os.environ["APPDIR"]
            inAppImage = True

        except:
            inAppImage = False
            shebang = '#!/bin/bash'
            AppDir = ''

        #get part directory
        if self.partDir[-1] != '/':
            self.partDir += '/'

        #determine number of logical cpus and leave 1 for overhead
        if psutil.cpu_count(logical=True) < 2:
            NCPU = 1
        else:
            NCPU = psutil.cpu_count(logical=True) - 1


#        #write openfoam sourcing script to generate environment variables
#        #basically replaces the default OF bashrc if in appimage
#        sourceFile = self.partDir + 'sourceOF'
#        with open(sourceFile, 'w') as f:
#            f.write(shebang + '\n')
#            if inAppImage == True:
#                f.write('export WM_PROJECT=OpenFOAM\n')
#                f.write('export WM_PROJECT_VERSION=v1912\n')
#                f.write('export WM_COMPILER_TYPE=system\n')
#                f.write('export WM_COMPILER=Gcc\n')
#                f.write('export WM_PRECISION_OPTION=DP\n')
#                f.write('export WM_LABEL_SIZE=32\n')
#                f.write('export WM_COMPILE_OPTION=Opt\n')
#                f.write('export WM_MPLIB=SYSTEMOPENMPI\n')
#                f.write('export WM_PROJECT_DIR='+AppDir+'/usr/lib/openfoam/openfoam1912\n')
#                f.write('source $WM_PROJECT_DIR/etc/config.sh/setup')
#
#            else:
#                f.write(self.cmdSourceOF + '\n')
#
#        os.chmod(sourceFile, 0o775)

        #Write 3D meshing script
        file = self.partDir + self.cmd3Dmesh
        with open(file, 'w') as f:
            f.write(shebang + '\n')
            #source OF bashrc if in dev mode (already sourced in appImage)
            if inAppImage == False:
                f.write(self.cmdSourceOF + '\n')
            f.write('blockMesh > ' + logFile + '\n')
            #single core meshing
            if parallel==False:
                f.write('snappyHexMesh -overwrite > '+logFile+ '\n')
            #parallel meshing
            else:
                #Run snappyHexMesh across multiple processors
                #this feature should be added in one day...sigh...for now its a placeholder
                #you'll have to set up blocks in blockMeshDict or this will not work.
                #problem is you need to allocate processors in advance.
                #see: "2.3.11 Running in parallel"
                #from https://cfd.direct/openfoam/user-guide/v8-damBreak/#x7-610002.3.11
                f.write('decomposePar > ' + logFile + '\n')
                f.write('mpirun --use-hwthread-cpus -np '+str(NCPU)+' snappyHexMesh -parallel -overwrite > '+logFile+ '\n')
                f.write('reconstructParMesh -mergeTol 1e-6 -latestTime -constant > '+logFile+ '\n')

        os.chmod(file, 0o744)

        #Write thermal analysis script
        file = self.partDir + self.cmdThermal
        foamFile = self.partDir + '*.foam'
        with open(file, 'w') as f:
            f.write(shebang + '\n')
            #source OF bashrc if in dev mode (already sourced in appImage)
            if inAppImage == False:
                f.write(self.cmdSourceOF + '\n')
            f.write('topoSet > '+logFile+ '\n')
            f.write('createPatch -overwrite > '+logFile+ '\n')
            f.write('heatFoam > '+logFile+ '\n')
            #f.write('paraFoam -touchAll\n')
            #f.write('touch '+ self.partName +'.foam\n')
        os.chmod(file, 0o775)

        #Write temperature probe script
        file = self.partDir + self.cmdTprobe
        with open(file, 'w') as f:
            f.write(shebang + '\n')
            #source OF bashrc if in dev mode (already sourced in appImage)
            if inAppImage == False:
                f.write(self.cmdSourceOF + '\n')
            f.write('topoSet > '+logFile+ '\n')
            f.write('createPatch -overwrite > '+logFile+ '\n')
            f.write('postProcess -func "probes" > '+logFile+ '\n')
        os.chmod(file, 0o744)
        return


    def generate3Dmesh(self, part):
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
        minLev = self.meshMinLevel
        maxLev = self.meshMaxLevel
        if self.meshDir[-1] == '/':
            file = self.meshDir + part + '_{:d}_{:d}'.format(minLev,maxLev)
        else:
            file = self.meshDir + '/' + part + '_{:d}_{:d}'.format(minLev,maxLev)

        if self.partDir[-1] != '/':
            self.partDir += '/'

        #if 3D mesh doesn't exist, create it.  Otherwise copy it to OF.partDir/constant/polyMesh
        newFile = self.partDir+'constant/polyMesh'
        if os.path.isdir(file) == True: #polyMesh already exists
            try:
                shutil.copytree(file, newFile)
                print('Copied existing 3D mesh')
            except:
                print('\nProblem copying 3D mesh...possibly doesnt exist...creating')
                print('tried copying: '+file)
                print('to this location: '+newFile+'\n')
        else:
            from subprocess import run
            #Copy the current environment
            current_env = os.environ.copy()
            try:
                appDir = os.environ["APPDIR"]
                inAppImage = True
            except:
                appDir = ''
                inAppImage = False
            #run blockMesh, snappyHexMesh
            meshCMD = self.partDir+self.cmd3Dmesh
            try:
                p = run([meshCMD], env=current_env, cwd=self.partDir, shell=True, executable=appDir+'/bin/bash')
                retcode = p.returncode
                if retcode < 0:
                    print("OF mesh child was terminated by signal", -retcode, file=sys.stderr)
                else:
                    print("OF mesh child returned", retcode, file=sys.stderr)
            except OSError as e:
                print("Execution failed:", e, file=sys.stderr)

            #Now copy this mesh to the 3D meshes folder for future use
            try:
                os.mkdir(newFile)
            except:
                print("Could not create polyMesh directory in OF folder")
                log.info("Could not create polyMesh directory in OF folder")
            shutil.copytree(newFile, file)
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

        from subprocess import run
        #Copy the current environment
        current_env = os.environ.copy()
        try:
            appDir = os.environ["APPDIR"]
            inAppImage = True
        except:
            appDir = ''
            inAppImage = False

        #run topoSet, createPatch, heatFoam, paraFoam -touchAll
        thermalCMD = self.partDir+self.cmdThermal
        try:
            p = run([thermalCMD], env=current_env, cwd=self.partDir, shell=True, executable=appDir+'/bin/bash')
            retcode = p.returncode
            if retcode < 0:
                print("thermal analysis child was terminated by signal", -retcode, file=sys.stderr)
            else:
                print("thermal analysis child returned", retcode, file=sys.stderr)

        except OSError as e:
            print("Execution failed:", e, file=sys.stderr)

        from pathlib import Path
        Path(self.partDir+self.partName+'.foam').touch()

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
        appDir = os.environ["APPDIR"]

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
            p = run([TprobeCMD], env=current_env, cwd=self.partDir, shell=True, executable=appDir+'/bin/bash')
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
        #first clean the tabs out of csv
        #Equivalent of: sed 's/\t/,/g' file > outfile
        with open(file, 'r') as fin:
            with open(outfile, 'w') as fout:
                for line in fin:
                    fout.write(line.replace('\t',','))
        #read data file
        data = pd.read_csv(outfile, header=1)
        data.columns = data.columns.str.strip()
        data = data.sort_values('field')
        data['field'] = data['field'].str.strip()

        return data

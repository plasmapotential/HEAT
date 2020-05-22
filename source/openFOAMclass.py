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
import matplotlib.pyplot as plt
import logging
log = logging.getLogger(__name__)

class OpenFOAM():
    def __init__(self):
        tools = toolsClass.tools()

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
                             'STLlayerName'
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
                        'createPatchDict' #, 'probes'
                        ]
        #generate dictionaries from template files
        for file in templateList:
            inFile = templateDir + file + '.template'
            outFile = partDir + 'system/' + file
            cmd = ('pyFoamFromTemplate.py'
                   ' --template-file='+inFile+
                   ' --output-file='+outFile+
                   ' --values-dictionary='+templateVarFile)

            #run pyFoamFromTemplate shell command to generate openfoam dictionaries
            os.system(cmd)


    def writeOFtemplateVarFile(self, file, STLpart):
        """
        writes a file that pyfoam uses to populate the OF dictionary templates
        with real variables.

        file        path to the file we will write
        STLpart     name of STL that openfoam will create 3D mesh from
        rootDir     path to HEAT root
        """
        with open(file, 'w') as f:
            for var in self.allowed_vars:
                f.write(var+' '+str(getattr(self,var))+';\n')
            f.write('STLfileName '+STLpart+';\n')
            f.write('STLlayerName '+STLpart+'_firstSolid;\n')

    def writeShellScript(self,logFile):
        """
        Writes shell script that runs OF from the template file.  This file will
            write to the HEATlog.txt file for viewing in the GUI

        logFile is absolute location of HEAT log file

        """
        if self.partDir[-1] != '/':
            self.partDir += '/'
        file = self.partDir + self.cmd3Dmesh

        #Write 3D meshing script
        with open(file, 'w') as f:
            f.write(self.cmdSourceOF + '\n')
            f.write('blockMesh > ' + logFile + '\n')
            #Run snappyHexMesh across multiple processors
            f.write('decomposePar > ' + logFile + '\n')
            #f.write('snappyHexMesh -overwrite > '+logFile+ '\n')
            f.write('mpirun -np 8 snappyHexMesh -parallel -overwrite > '+logFile+ '\n')
            f.write('reconstructParMesh -mergeTol 1e-6 -latestTime > '+logFile+ '\n')

        os.chmod(file, 0o744)

        #Write thermal analysis script
        file = self.partDir + self.cmdThermal
        foamFile = self.partDir + '*.foam'
        with open(file, 'w') as f:
            f.write(self.cmdSourceOF + '\n')
            f.write('topoSet > '+logFile+ '\n')
            f.write('createPatch -overwrite > '+logFile+ '\n')
            f.write('heatFoam > '+logFile+ '\n')
            f.write('paraFoam -touchAll'+ '\n')
        os.chmod(file, 0o744)

        #Write temperature probe script
        file = self.partDir + self.cmdTprobe
        with open(file, 'w') as f:
            f.write(self.cmdSourceOF + '\n')
            f.write('topoSet > '+logFile+ '\n')
            f.write('createPatch -overwrite > '+logFile+ '\n')
            f.write('postProcess -func "probes" > '+logFile+ '\n')
        os.chmod(file, 0o744)


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
        minLev = self.dict['meshMinLevel']
        maxLev = self.dict['meshMaxLevel']
        if self.meshDir[-1] == '/':
            file = self.meshDir + part + '_' + minLev + '_' + maxLev
        else:
            file = self.meshDir + '/' + part + '_' + minLev + '_' + maxLev

        if self.partDir[-1] != '/':
            self.partDir += '/'

        #if 3D mesh doesn't exist, create it.  Otherwise copy it to OF.partDir/constant/polyMesh
        newFile = self.partDir+'constant/polyMesh'
        if os.path.isdir(file) == True: #polyMesh already exists
            try:
                shutil.copytree(file, newFile)
            except:
                print('Problem copying 3D mesh...')
        else:
            #run blockMesh, snappyHexMesh, topoSet, createPatch
            subprocess.call(self.partDir+self.cmd3Dmesh, shell=True, cwd=self.partDir)
            #Now copy this mesh to the 3D meshes folder for future use
            shutil.copytree(newFile, file)

#            from PyFoam.Execution.BasicRunner import BasicRunner
#            solvers = ['blockMesh','snappyHexMesh','topoSet','createPatch']
#            dir = self.OF.partDir
#            case = self.OF.partDir
#            for solver in solvers:
#                run = BasicRunner(argv=[solver,dir,case])
#                run.start()

    def runThermalAnalysis(self):
        """
        runs an openfoam thermal analysis in the specified location
        """
        print('Running Thermal Analysis')
        log.info('Running Thermal Analysis')

        subprocess.call(self.partDir+self.cmdThermal, shell=True, cwd=self.partDir)




    def runTprobe(self,x,y,z):
        """
        runs openfoam postprocessing function to find temperature at specified
        x,y,z probe location
        """
        print('Simulating Temperature Probe')
        log.info('Simulating Temperature Probe')

        #write probe file to partDir
        file = self.partDir+'system/probes'
        with open(file,'w') as f:
            f.write('type probes;\n')
            f.write('libs ("libsampling.so");\n')
            f.write('fields (T);\n')
            f.write('probeLocations ( ({:.6f} {:.6f} {:.6f}) );\n'.format(x,y,z))


        subprocess.call(self.partDir+self.cmdTprobe, shell=True, cwd=self.partDir)

    def plotTprobes(self, file):
        """
        generates plot with all temperature probes on single plot

        file is filename to import

        uses genfromtxt to get data.
        first column is timesteps.  each additional column is an additional
            temperature probe
        """
        data = np.genfromtxt(file)
        #for idx in range(data.shape[1]-1):
        #    plt.plot(data[:,0], data[:,idx+1])

        plt.style.use('dark_background')
        plt.plot(data[:,0], data[:,1], 'o-')
        plt.title('Temperature Probe')
        plt.xlabel('Time [s]')
        plt.ylabel('Temperature [K]')

        path = os.path.dirname(file)
        plotFile = path+'/Tprobe.png'
        plt.savefig(plotFile,type="png",dpi=300)

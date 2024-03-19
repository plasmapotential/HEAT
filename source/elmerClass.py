#elmerClass.py
#Description:   HEAT interface to Elmer FEM
#Engineer:      T Looby
#Date:          20240206

"""
Class for connecting HEAT output to Elmer FEM

https://www.elmerfem.org/blog/

https://github.com/ElmerCSC/elmerfem

Elmer:  https://doi.org/10.5281/zenodo.7892181
"""
import sys
import os
import pandas as pd
import numpy as np
import shutil
from scipy.interpolate import griddata
import gmsh


#HEAT classes
import toolsClass
tools = toolsClass.tools()
import ioClass
IO = ioClass.IO_HEAT()

import logging
log = logging.getLogger(__name__)

class FEM:

    def __init__(self, rootDir:str, dataPath:str, chmod=0o774, UID=-1, GID=-1):
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
        return
    
    def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
        """
        sets up pythonic string number formats for shot and timesteps
        """
        self.tsSigFigs = tsSigFigs
        self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
        self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
        return

    def allowed_class_vars(self):
        """
        .. Writes a list of recognized class variables to HEAT object
        .. Used for error checking input files and for initialization
          
        Elmer Variables:
        -------------------

        :meshFEMmaxRes: float that indicated maximum mesh size in [mm].  0.0 means auto meshing
        :meshFEMminRes: float that indicated minimum mesh size in [mm].  0.0 means auto meshing
        :elmerDir: path to the Elmer case directory.  In container should be location that
          is bind mounted into container.  Requires full path.
        :elmerFile: CSV file that contains information about the Elmer run.  See
          function "readElmerFile()" for more information.  This file should be location
          in the elmerDir.
        :elmerHEATlib: Fortan shared object (.so) file that is called from the SIF.  This
          file should be located in the elmerDir.  We allow this file to be dynamically 
          loaded by the user so that they can employ any Elmer User Defined Function.
                
        """
        self.allowed_vars = [
                             'meshFEMminRes',
                             'meshFEMmaxRes',
                             'elmerDir',
                             'elmerFile',
                             'elmerHEATlib'
                             ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file

        """
        integers = [
                    ]
        floats = [
            'meshFEMminRes',
            'meshFEMmaxRes'
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
    
    def readElmerFile(self):
        """
        The HEAT Elmer File:
        ------------------

        The elmerFile is a file that describes the Elmer FEM analysis to run.
        It defines which PFCs to run through the FEM solver, as well as the
        Elmer FEM Solver Input Files (.SIF) for each PFC.  The columns are:

        :PFCname:  the name of the PFC to solve on.  Should match the CAD.
        :SIF:  the name of the .SIF file to use.  The SIF must be saved in 
          the Elmer Directory as declared in the HEAT input file.
        :meshFile: if the user wants to supply a FEM mesh in universal mesh
          (.unv) format, then they may include the name of the file here.
          The meshFile must be saved in the Elmer Directory as declared in
          the HEAT input file.  If the user wants HEAT to create a FEM mesh,
          then leave this field blank or write 'None'.  Keep in mind that if
          HEAT creates the mesh file autonomously, then there will not be 
          boundary groups used to constrain the model.


        """
        if self.elmerDir[-1] != '/': self.elmerDir +='/'
        if self.elmerFile is None:
            print("No elmerFile!  Cannot run Elmer without an elmerFile")
            log.info("No elmerFile!  Cannot run Elmer without an elmerFile")
            sys.exit()

        print("Reading elmerFile: " + self.elmerDir + self.elmerFile)
        log.info("Reading elmerFile: " + self.elmerDir + self.elmerFile)
        try:
            #self.elmerData = pd.read_csv(self.elmerDir + self.elmerFile, comment='#')
            #self.elmerData = self.elmerData.rename(columns=lambda x: x.strip())
            #self.elmerData['PFCname'] = self.elmerData['PFCname'].str.strip()
            #self.elmerData['SIF'] = self.elmerData['SIF'].str.strip()
            #self.elmerData['meshFile'] = self.elmerData['meshFile'].str.strip()
            import csv
            self.elmerData = {}
            with open(self.elmerDir + self.elmerFile, mode='r') as infile:
                reader = csv.reader(infile)
                for row in reader:
                    if 'PFCname' in row: #header row
                        header = row
                    else:
                        self.elmerData[row[0].strip()] = {header[1].strip():row[1].strip(), header[2].strip():row[2].strip()}
            print("Elmer Data:")
            print(self.elmerData)
            log.info("Elmer Data:")
            log.info(self.elmerData)

        except Exception as e:
            print(str(e))
            log.info(str(e))
            print("Could not read elmerFile!")
            log.info("Could not read elmerFile!")
        return
    
    def buildElmerMesh(self, meshDir, meshFile):
        """
        builds an Elmer FEM mesh
        """
        tools.makeDir(meshDir, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)
        #rhel is refinement level, autoclean reorders boundaries and removes unused elements
        #args = ['ElmerGrid', '8', '2', meshFile, '-autoclean', '-relh', '1.0', '-out', meshDir]
        #args = ['ElmerGrid', '8', '2', meshFile, '-relh', '1.0', '-out', meshDir]
        args = ['ElmerGrid', '8', '2', meshFile, '-autoclean', '-out', meshDir]

        current_env = os.environ.copy()
        #run ElmerGrid
        from subprocess import run
        run(args, env=current_env, cwd=meshDir)
        #update permissions
        tools.recursivePermissions(meshDir, self.UID, self.GID, self.chmod)        
        return
    
    def runElmerSolve(self, SIFfile, name):
        """
        runs the Elmer Solver
        """
        #copy SIF from elmerDir to elmerOutDir
        src = self.elmerDir + SIFfile
        dst = self.elmerOutDir + SIFfile
        #shutil.copyfile(src, dst)
        #overwrite Mesh DB line in SIF if necessary
        with open(src, 'r') as f:
            with open(dst, 'w') as f2:
                for line in f:
                    if 'Mesh DB' in line:
                        f2.write("  Mesh DB \""+self.elmerOutDir[:-1]+"\" \""+name+"\" \n")
                    else:
                        f2.write(line)

        #copy shared object file (.so) to elmerOutDir
        src = self.elmerDir + self.elmerHEATlib
        dst = self.elmerOutDir + self.elmerHEATlib
        shutil.copyfile(src, dst)


        args = ['ElmerSolver', SIFfile]
        current_env = os.environ.copy()
        #run Elmer Solver
        from subprocess import run
        run(args, env=current_env, cwd=self.elmerOutDir)
        return
    
    def interpolateHFtoMesh(self, PFC, t, tMin, hfFile=None):
        """
        Takes a PFC object, interpolates the qDiv result onto the surface node of
        a corresponding volume mesh, then saves a .csv file with columns nodeId, MW/m2.

        SIFfile must contain a Variable, nodalHFprefix, in BC section.  
        t is the timestep.        
        """
        #get the prefix
        src = self.elmerDir + PFC.SIFfile
        with open(src, 'r') as f:
            for line in f:
                if "nodalHFprefix" in line:
                    prefix = line.split("String ")[-1].strip()

        gmsh.initialize()
        gmsh.open(PFC.meshFile)
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        node_tags = np.array(node_tags)
        node_coords = np.array(node_coords).reshape(-1, 3)
        boundary_nodes = set()

        # Get the entities of the dimension you are interested in
        entities = gmsh.model.getEntities(2)  # 2 for surfaces

        for entity in entities:
            dim, tag = entity
            # Get the nodes of this surface entity
            _, _, surface_node_tags = gmsh.model.mesh.getElements(dim, tag)
            for node_tag_array in surface_node_tags:
                for node_tag in node_tag_array:
                    boundary_nodes.add(node_tag)  # Use add() for single elements

        # Extract the coordinates of boundary nodes
        boundaryPts = np.array([node_coords[np.where(node_tags == node_tag)[0][0]] for node_tag in boundary_nodes])

        # Perform the interpolation
        print("Interpolating HF to mesh surface...")
        #if None, then we assign zeros to the surface
        if hfFile == None:
            hfOnMesh = np.zeros((len(boundaryPts)))
        #if there is a hfFile, then we assign that HF to the surface
        else:
            hf_data = pd.read_csv(hfFile)
            # Extracting coordinates and heat flux from the CSV file
            points = hf_data[['# X', 'Y', 'Z']].values
            values = hf_data['$MW/m^2$'].values        
            hfOnMesh = griddata(points, values, boundaryPts, method='nearest', fill_value=0.0)
            #hfArray = np.hstack((boundaryPts, hfOnMesh[:,np.newaxis]))
            #nodes = np.array(list(boundary_nodes), dtype=int)

        #convert HF to Watts
        nodeArray = np.vstack([np.array(list(boundary_nodes), dtype=int), hfOnMesh*1e6]).T

        #for saving nodeId,HF
        name = prefix + '_' + self.tsFmt.format(t) + '.dat'
        np.savetxt(self.elmerOutDir + name, nodeArray, fmt="% .9E", delimiter=',')

        return
    
    def interpolateHFinTime(self, hfFile, hfFileNext, tLast, tNext, t):
        """
        given a heat flux file from last timestep, hfFile, and a heat flux
        file for the next timestep, hfFileNext, along with the times of
        the last timestep, next timestep, and current timestep, t,
        interpolate the heat flux at each xyz position in the hfFiles
        """
        hf_data_last = pd.read_csv(hfFile)
        # Extracting coordinates and heat flux from the CSV file
        xyz = hf_data_last[['# X', 'Y', 'Z']].values
        qLast = hf_data_last['$MW/m^2$'].values    

        hf_data_next = pd.read_csv(hfFileNext)
        qNext = hf_data_next['$MW/m^2$'].values   

        mult = (t-tLast) / (tNext-tLast)
        q = (qNext - qLast)*mult + qLast

        hfArray = np.hstack((xyz, q[:,np.newaxis]))
        newFile = self.elmerOutDir + 'HFtmp.csv'
        head = "X,Y,Z,$MW/m^2$"
        np.savetxt(newFile, hfArray, fmt="% .9E", delimiter=',', header=head)

        return newFile

    def buildTimesteps(self,SIFfile):
        """
        Parses an Elmer Solver Input File (SIF) and extracts the timestep
        intervals and sizes.  builds an Elmer timestep array, self.ts
        that contains all the timesteps Elmer will run for
        """
        #parse the SIF file and extract the timestep data
        SIF = self.elmerDir + SIFfile
        with open(SIF, 'r') as f:
            for line in f:
                if "Timestep intervals" in line:
                    tmp = line.split("(")[1]
                    lenInt = int(tmp.split(")")[0].strip())
                    tmp = line.split("=")[1].strip()
                    tmp2 = tmp.split(" ")
                    test = [x.isdigit() for x in tmp2]
                    idx = np.where(test)[0]
                    intervals = np.array(tmp2)[idx].astype(int)
                
                if "Timestep Sizes" in line:
                    tmp = line.split("=")[1].strip()
                    tmp2 = tmp.split(" ")
                    sizes = np.array(tmp2).astype(float)        

        #check that SIF is formatted correctly
        test1 = lenInt == len(intervals)
        test2 = lenInt == len(sizes)
        if np.logical_and(test1, test2) != True:
            print("You have a mismatch in your timestep rows in Elmer SIF. This will fail...")
            log.info("You have a mismatch in your timestep rows in Elmer SIF. This will fail...")

        #build a timestep array using this data that corresponds 
        #to what Elmer will solve
        ts = []
        for i in range(lenInt):
            itvl = intervals[i]
            s = sizes[i]
            for j in range(itvl+1):
                if len(ts) == 0:
                    val = 0.0
                else:
                    val = np.round(ts[-1] + s, self.tsSigFigs)
                ts.append(val)

        self.ts = ts
        self.lenInt = lenInt
        self.timeIntervals = intervals
        self.timeSizes = sizes
        
        return
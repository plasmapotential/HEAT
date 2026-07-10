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
from scipy.spatial import cKDTree
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
        self._hf_mapping_active = False
        self._pfc_boundary_cache = {}
        self._hf_csv_cache = {}
        self._hf_kdtree_cache = {}
        self._nodal_hf_prefix_cache = {}
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
                             'elmerHEATlib',
                             'numberpartitions'
                             ]
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file

        """
        integers = ['numberpartitions'
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
    
    def begin_elmer_hf_mapping(self, n_timesteps=None, n_pfcs=None):
        """
        Start HEAT -> Elmer nodal HF mapping: one GMSH session and per-run caches.
        """
        self._pfc_boundary_cache = {}
        self._hf_csv_cache = {}
        self._hf_kdtree_cache = {}
        self._nodal_hf_prefix_cache = {}
        if not self._hf_mapping_active:
            gmsh.initialize()
            self._hf_mapping_active = True
        if n_timesteps is not None and n_pfcs is not None:
            total = int(n_timesteps) * int(n_pfcs)
            msg = (
                "Starting Elmer HF mapping: {} Elmer timesteps x {} PFCs "
                "({} nodal HF maps)".format(n_timesteps, n_pfcs, total)
            )
            print(msg)
            log.info(msg)

    def log_hf_mapping_progress(self, step, total, t, PFC, mode):
        """Progress line for one PFC/timestep nodal HF map."""
        msg = "Elmer HF mapping [{}/{}]: t={} s, {}, {}".format(
            step, total, self.tsFmt.format(t), PFC.name, mode)
        print(msg)
        log.info(msg)

    def end_elmer_hf_mapping(self, n_maps=None):
        """Release GMSH and mapping caches after nodal HF files are written."""
        if self._hf_mapping_active:
            gmsh.finalize()
            self._hf_mapping_active = False
        self._pfc_boundary_cache = {}
        self._hf_csv_cache = {}
        self._hf_kdtree_cache = {}
        self._nodal_hf_prefix_cache = {}
        if n_maps is not None:
            msg = "Elmer HF mapping complete ({} nodal HF files written)".format(n_maps)
            print(msg)
            log.info(msg)

    def _read_nodal_hf_prefix(self, siffile):
        """Read nodalHFprefix from a SIF (cached by SIF path)."""
        if siffile in self._nodal_hf_prefix_cache:
            return self._nodal_hf_prefix_cache[siffile]
        src = self.elmerDir + siffile
        prefix = None
        with open(src, 'r') as f:
            for line in f:
                if "nodalHFprefix" in line:
                    prefix = line.split("String ")[-1].strip()
                    break
        if prefix is None:
            raise ValueError("nodalHFprefix not found in SIF: " + src)
        self._nodal_hf_prefix_cache[siffile] = prefix
        return prefix

    def _load_hf_csv(self, hfFile):
        """
        Load HEAT HF_allSources CSV as (xyz, q_MW_m2).

        hfFile may be a filesystem path or an in-memory (xyz, q) tuple from
        temporal interpolation.
        """
        if isinstance(hfFile, tuple):
            return hfFile
        if hfFile not in self._hf_csv_cache:
            hf_data = pd.read_csv(hfFile)
            xyz = hf_data[['# X', 'Y', 'Z']].to_numpy()
            q = hf_data['$MW/m^2$'].to_numpy()
            self._hf_csv_cache[hfFile] = (xyz, q)
        return self._hf_csv_cache[hfFile]

    def _get_pfc_boundary_cache(self, PFC):
        """
        Boundary nodes and coordinates for a PFC volume mesh (cached per meshFile).
        """
        key = PFC.meshFile
        if key in self._pfc_boundary_cache:
            return self._pfc_boundary_cache[key]

        if not self._hf_mapping_active:
            gmsh.initialize()
            self._hf_mapping_active = True

        gmsh.open(PFC.meshFile)
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        node_tags = np.asarray(node_tags, dtype=np.int64)
        node_coords = np.asarray(node_coords, dtype=np.float64).reshape(-1, 3)
        tag_to_idx = {int(tag): i for i, tag in enumerate(node_tags)}

        boundary_nodes = set()
        for dim, tag in gmsh.model.getEntities(2):
            _, _, surface_node_tags = gmsh.model.mesh.getElements(dim, tag)
            for node_tag_array in surface_node_tags:
                for node_tag in node_tag_array:
                    boundary_nodes.add(int(node_tag))

        boundary_nodes_arr = np.array(sorted(boundary_nodes), dtype=np.int64)
        boundary_indices = np.array([tag_to_idx[int(t)] for t in boundary_nodes_arr])
        boundaryPts = node_coords[boundary_indices]
        prefix = self._read_nodal_hf_prefix(PFC.SIFfile)

        self._pfc_boundary_cache[key] = {
            'boundary_nodes': boundary_nodes_arr,
            'boundaryPts': boundaryPts,
            'nodalHFprefix': prefix,
        }
        print("Cached Elmer boundary mesh for {} ({} nodes)".format(
            PFC.name, len(boundary_nodes_arr)))
        log.info("Cached Elmer boundary mesh for {} ({} nodes)".format(
            PFC.name, len(boundary_nodes_arr)))
        return self._pfc_boundary_cache[key]

    def _hf_on_boundary(self, boundaryPts, hfFile, tree_key=None):
        """Nearest-neighbor map of HEAT HF onto Elmer boundary nodes (cKDTree)."""
        xyz, values = self._load_hf_csv(hfFile)
        if tree_key is None:
            tree_key = hfFile if isinstance(hfFile, str) else ('inline', id(xyz))
        if tree_key not in self._hf_kdtree_cache:
            self._hf_kdtree_cache[tree_key] = cKDTree(xyz)
        idx = self._hf_kdtree_cache[tree_key].query(boundaryPts, k=1)[1]
        return values[idx]
    
    def buildElmerMesh(self, meshDir, meshFile):
        """
        builds an Elmer FEM mesh
        """
        tools.makeDir(meshDir, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)
        #rhel is refinement level, autoclean reorders boundaries and removes unused elements
        args = ['ElmerGrid', '8', '2', meshFile, '-autoclean', '-relh', '1.0', '-out', meshDir]
        #args = ['ElmerGrid', '8', '2', meshFile, '-relh', '1.0', '-out', meshDir]
        #args = ['ElmerGrid', '8', '2', meshFile, '-autoclean', '-out', meshDir]
        #args = ['ElmerGrid', '8', '2', meshFile, '-out', meshDir]
        #args = ['ElmerGrid', '8', '2', meshFile, '-boundorder', '-removeunused', '-out', meshDir]

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

        if self.numberpartitions > 0:
            args_mesh = ['ElmerGrid', '2', '2', name, '-metis', str(self.numberpartitions)]
            args = ['mpirun', '-np', str(self.numberpartitions), 'ElmerSolver', SIFfile]
            current_env = os.environ.copy()
		    #run Elmer Solver
            from subprocess import run
            run(args_mesh, env=current_env, cwd=self.elmerOutDir)		
            run(args, env=current_env, cwd=self.elmerOutDir)
            try:
                self.merge_Rex(name, SIFfile)
            except:
                print('no ReX calcs were done')
        else:
            args = ['ElmerSolver', SIFfile]
            current_env = os.environ.copy()
            #run Elmer Solver
            from subprocess import run
            run(args, env=current_env, cwd=self.elmerOutDir)
        return
    
    def interpolateHFtoMesh(self, PFC, t, tMin, hfFile=None, tree_key=None):
        """
        Takes a PFC object, interpolates the qDiv result onto the surface node of
        a corresponding volume mesh, then saves a .csv file with columns nodeId, MW/m2.

        SIFfile must contain a Variable, nodalHFprefix, in BC section.  
        t is the timestep.

        hfFile may be None, a path to HF_allSources.csv, or an in-memory (xyz, q) tuple.
        tree_key reuses the cKDTree built for a prior MHD CSV when xyz is unchanged.
        """
        cache = self._get_pfc_boundary_cache(PFC)
        boundary_nodes = cache['boundary_nodes']
        boundaryPts = cache['boundaryPts']
        prefix = cache['nodalHFprefix']

        if hfFile is None:
            hfOnMesh = np.zeros(len(boundaryPts), dtype=np.float64)
        else:
            hfOnMesh = self._hf_on_boundary(boundaryPts, hfFile, tree_key=tree_key)

        hfOnMesh = np.nan_to_num(hfOnMesh, nan=0.0, posinf=0.0, neginf=0.0)

        #convert HF to Watts
        nodeArray = np.column_stack([boundary_nodes, hfOnMesh * 1e6])

        name = prefix + '_' + self.tsFmt.format(t) + '.dat'
        out_path = self.elmerOutDir + name
        tools.savetxt(out_path, nodeArray, fmt="% .9E", delimiter=',')

        return

    def copyReXinit(self, PFC):
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
                if "nodalReXprefix" in line:
                    prefix = line.split("String ")[-1].strip()
        try:
            ReXfile = prefix + '.dat'
            #copy Rex init from elmerDir to elmerOutDir
            src = self.elmerDir + ReXfile
            dst = self.elmerOutDir + ReXfile
            shutil.copyfile(src, dst)
        except:
            print('no Rex init file provided')
        return

    def merge_Rex(self, name, SIFfile):
       import glob
	   # 1. Find all partition files created by your Fortran module and merge them to one
       #get the output name from the .sif file
       src = self.elmerDir + SIFfile
       with open(src, 'r') as f:
            for line in f:
                if "Filename" in line:
                    prefix = line.split("= ")[-1].strip().strip('"')
       search_pattern = self.elmerOutDir + prefix + '.p*'
       files = glob.glob(search_pattern)
       data_list = []
       for f in files:
          # Load each file (assuming they are comma-separated as we programmed)
          data_list.append(np.loadtxt(f, delimiter=','))
       merged_data = np.vstack(data_list)
       # This ensures your final file is perfectly ordered from Node 1 to Node N
       sorted_data = merged_data[merged_data[:, 0].argsort()]
       # 4. Save to a single, clean .dat file
       np.savetxt(self.elmerOutDir+prefix, sorted_data, delimiter=',', fmt=['%d', '%15.9E'])
       return
        
    def interpolateHFinTime(self, hfFile, hfFileNext, tLast, tNext, t):
        """
        given a heat flux file from last timestep, hfFile, and a heat flux
        file for the next timestep, hfFileNext, along with the times of
        the last timestep, next timestep, and current timestep, t,
        interpolate the heat flux at each xyz position in the hfFiles.

        Returns (xyz, q_MW_m2) for in-memory use (no HFtmp.csv disk write).
        """
        xyz, qLast = self._load_hf_csv(hfFile)
        _, qNext = self._load_hf_csv(hfFileNext)

        mult = (t - tLast) / (tNext - tLast)
        q = (qNext - qLast) * mult + qLast

        return (xyz, q)

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

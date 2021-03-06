#toolsClass.py
#Description:   Utility tools for HEAT
#Engineer:      T Looby
#Date:          20191117
import sys
import pandas as pd
import numpy as np
import os
import shutil
import logging
import subprocess
log = logging.getLogger(__name__)

class tools:
    """
    These are tools that are used by multiple modules in HEAT.  Stuff like
    setting up input file pipelines, generic coordinate transformations, etc.
    """

    def __init__(self):
        pass
        return

    def initializeInput(self, obj, infile=None):
        """
        Pulls input data for Command Line Interface (CLI)
        """
        #Get class variables
        obj.allowed_class_vars()

        #Set all variables to None
        self.vars2None(obj)

        if infile == None:
            print("You must provide an input file!")
            sys.exit()
        else:
            self.read_input_file(obj,infile)

        #After importing from file, set data types
        obj.setTypes()

        return

    def vars2None(self,obj):
        """
        Set all variables in allowed_vars to None
        """
        for var in obj.allowed_vars:
            setattr(obj, var, None)

    def read_input_file(self, obj, infile):
        """
        Reads input file and fills object with corresponding parameters

        Inputs:
        infile      path of input file

        Format for input file is comma delimited, # are comments.
        Example:
        #Important Comment
        variable_name,value

        """
        if infile == None:
            print("No input file.  Please provide input file")
            sys.exit()

        data = pd.read_csv(infile, sep=',', comment='#', names=['Var','Val'], skipinitialspace=True)

        #Dynamically initialize class variables
        for i in range(len(data)):
            if data['Var'][i] in obj.allowed_vars:
                setattr(obj, data['Var'][i], data['Val'][i])
            else:
                pass
                #print("Caution! Unrecognized variable in input file: "+data['Var'][i])

        return

    def saveInputFile(self, data, path, rootDir, dataPath):
        """
        saves a HEAT formatted data file from data, which is a dict
        file is saved to path location with the name HEATinput.csv
        """
        import CADClass
        import MHDClass
        import heatfluxClass
        import openFOAMclass
        MHD = MHDClass.MHD(rootDir, dataPath)
        CAD = CADClass.CAD(rootDir, dataPath)
        HF = heatfluxClass.heatFlux(rootDir, dataPath)
        OF = openFOAMclass.OpenFOAM(rootDir, dataPath)
        MHD.allowed_class_vars()
        CAD.allowed_class_vars()
        HF.allowed_class_vars()
        OF.allowed_class_vars()

        if path[-1] == '/':
            file = path + 'HEATinput.csv'
        else:
            file = path + '/HEATinput.csv'
        with open(file, 'w') as f:
            f.write("# Input file for HEAT\n")
            f.write("# Format is: Variable, Value\n")
            f.write("#=============================================================\n")
            f.write("#                CAD Variables\n")
            f.write("#=============================================================\n")
            for var in CAD.allowed_vars:
                if var in data:
                    f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ', \n')
            f.write("#=============================================================\n")
            f.write("#                MHD Variables\n")
            f.write("#=============================================================\n")
            for var in MHD.allowed_vars:
                if var in data:
                    f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ', \n')
            f.write("#=============================================================\n")
            f.write("#                HF Variables\n")
            f.write("#=============================================================\n")
            for var in HF.allowed_vars:
                if var in data:
                    f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ', \n')
            f.write("#=============================================================\n")
            f.write("#                OpenFOAM Variables\n")
            f.write("#=============================================================\n")
            for var in OF.allowed_vars:
                if var in data:
                    f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ', \n')

        print("Wrote HEAT input file")
        log.info("Wrote HEAT input file")
        return

    def saveDefaultPFCfile(self, path):
        """
        saves a default PFC input file to path
        file is named PFCinput.csv
        """
        if path[-1] == '/':
            file = path + 'PFCinput.csv'
        else:
            file = path + '/PFCinput.csv'
        with open(file, 'w') as f:
            f.write("# This is a HEAT input file responsible for: \n")
            f.write("# 1) assigning intersection PFCs to each heat flux calculation PFC \n")
            f.write("# 2) assigning MapDirection to each PFC \n")
            f.write("# 3) assigning the timesteps to use each PFC for \n")
            f.write("# \n")
            f.write("# you can use colons (:) to insert a range of timesteps or to insert multiple intersects \n")
            f.write("# but each PFC needs its own line (and sometimes two if you need to trace in both directions) \n")
            f.write("# commas are the delimiter here, so don't use them inside columns or you'll break pandas \n")
            f.write("# \n")
            f.write("# Example \n")
            f.write("# timesteps, PFCname, MapDirection, DivCode, intersectName\n")
            f.write("# 100:105, E-ED1450-013, -1, E-ED1450-013:E-ED1390-003:E-ED1450-071 \n")
            f.write("# \n")
            f.write("# this would calculate HF on E-ED1450-013 for all (available) timesteps between \n")
            f.write("# 100ms and 105ms, checking for intersections with E-ED1450-013, E-ED1390-003, \n")
            f.write("# and E-ED1450-071, with traces running in the reversed MapDirection \n")
            f.write("# \n")
            f.write("# \n")
            f.write("# \n")
            f.write("# 40-1783_2 TILE UPPER \n")
            f.write("# 40-1782_3 TILE LOWER TYPE A0 \n")
            f.write("timesteps, PFCname, MapDirection, DivCode, intersectName\n")

        print("Default PFC file saved")
        log.info("Default PFC file saved")
        return

    def createDict(self, obj):
        """
        Creates and returns a dictionary object from an <obj>.allowed_vars
        """
        dict = {}
        for var in obj.allowed_vars:
            dict[var]=getattr(obj,var)
        return dict


    def inputs_from_dict(self, obj, settingsDict):
        """
        Assigns the variables from a python dictionary to variables
        belonging to an object.
        """
        #Dynamically initialize class variables
        for var in settingsDict:
            setattr(obj, var, settingsDict[var])


    def xyz2cyl(self,x,y,z):
        """
        Converts x,y,z coordinates to r,z,phi
        """
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y,x)
        #phi = np.radians(phi)
        return r,z,phi

    def cyl2xyz(self,r,z,phi):
        """
        Converts r,z,phi coordinates to x,y,z
        phi will be converted to radians
        """
        r = np.asarray(r)
        z = np.asarray(z)
        phi = np.radians(np.asarray(phi))
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        return x,y,z


    def signedVolume(self,a,b,c,d):
        """
        Calculates signed volume
        NOTE: true signed volume would be multiplied by 1/6 but we don't
        do that here
        """

        return np.sum(np.cross(b-a, c-a) * (d-a))

    def signedVolume2(self,a,b,c,d):
        """
        Calculates signed volume
        NOTE: true signed volume would be multiplied by 1/6 but we don't
        do that here
        """
        #return (1.0/6.0)*np.diagonal(np.matmul(np.cross(b-a,c-a),(d-a).T))
        #return (1.0/6.0) * np.dot(np.cross(b-a,c-a), d-a)
        return np.sum(np.cross(b-a, c-a) * (d-a), axis=1)

    def faceCenters(self, x, y, z):
        """
        returns centers of triangle in cartesian coordinates
        """
        #face centers
        centers = np.zeros((len(x), 3))
        centers[:,0] = np.sum(x,axis=1)/3.0
        centers[:,1] = np.sum(y,axis=1)/3.0
        centers[:,2] = np.sum(z,axis=1)/3.0
        return centers

    def createVTKOutput(self, pcfile, outType, prefix):
        """
        Creates a vtk file from csv output
        The vtk file will have some paraview operation performed, such as
        converting vector components to a glyph, or calculating a trace from points
        and so it requires declaration of the intended output type (outType).
        spawn process in separte session
        outType is a string that can be one of the following:
        'glyph'     Create glyph from vector components
        'trace'     create a line trace between points

        prefix is what you want to output vtk file to be named

        ParaVIEW requires a binary file called pvpython
        Here I load the path to this file into the PATH environment variable,
        and also load the paraVIEW python libs into the PYTHONPATH env var.
        I do this here (instead of in main HEAT code), because this allows
        ParaVIEW's python version to be different from HEAT's python version.
        If you have two different python versions (ie 3.7 for PV and 3.8 for HEAT)
        then there will be conflicts unless you isolate the environments (like
        I do here).
        """
        import os
        current_env = os.environ.copy()
        pvpythonCMD = current_env["pvpythonCMD"]
#        #running in appImage (isolate PV environment from HEAT's)
#        try:
#            pvpythonCMD = current_env["pvpythonCMD"]
#        #running on dev machine
#        #(it is expected that you have set up env externally, perhaps in dashGUI.py)
#        except:
#            pvpythonCMD = 'pvpython'
        print("Spawning PVpython subprocess")
        log.info("Spawning PVpython subprocess")
        args = [pvpythonCMD, self.rootDir + '/GUIscripts/csv2vtk.py', pcfile, outType, prefix]
        from subprocess import run
        run(args, env=current_env)
        print("PVpython subprocess complete")
        log.info("PVpython subprocess complete")
        return

    def intersectTestParallel(self, i):
        """
        Intersection test that does not check for self intersections
        """
        q13D = np.repeat(self.q1[i,np.newaxis], self.Nt, axis=0)
        q23D = np.repeat(self.q2[i,np.newaxis], self.Nt, axis=0)
        sign1 = np.sign(self.signedVolume2(q13D,self.p1,self.p2,self.p3))
        sign2 = np.sign(self.signedVolume2(q23D,self.p1,self.p2,self.p3))
        sign3 = np.sign(self.signedVolume2(q13D,q23D,self.p1,self.p2))
        sign4 = np.sign(self.signedVolume2(q13D,q23D,self.p2,self.p3))
        sign5 = np.sign(self.signedVolume2(q13D,q23D,self.p3,self.p1))
        test1 = (sign1 != sign2)
        test2 = np.logical_and(sign3==sign4,sign3==sign5)
        targetIdx = np.where(np.logical_and(test1,test2))[0]

        #if we found an intersection return 1
        if np.sum(np.logical_and(test1,test2)) > 0:
            result = 1
        else:
            result = 0

        return result

    def intersectTestParallel_selfCheck(self,i):
        """
        Intersection test for first step up field line.  Not only checks for
        intersections, but also determines if these intersections are self
        intersections
        """
        q13D = np.repeat(self.q1[i,np.newaxis], self.Nt, axis=0)
        q23D = np.repeat(self.q2[i,np.newaxis], self.Nt, axis=0)
        sign1 = np.sign(self.signedVolume2(q13D,self.p1,self.p2,self.p3))
        sign2 = np.sign(self.signedVolume2(q23D,self.p1,self.p2,self.p3))
        sign3 = np.sign(self.signedVolume2(q13D,q23D,self.p1,self.p2))
        sign4 = np.sign(self.signedVolume2(q13D,q23D,self.p2,self.p3))
        sign5 = np.sign(self.signedVolume2(q13D,q23D,self.p3,self.p1))
        test1 = (sign1 != sign2)
        test2 = np.logical_and(sign3==sign4,sign3==sign5)
        targetIdx = np.where(np.logical_and(test1,test2))[0]

        result = 0
        #if we found an intersection
        if np.sum(np.logical_and(test1,test2)) > 0:
            #If there are multiple intersections, we arent worried about self
            # intersection because multiple faces are intersected
            if len(targetIdx) > 1:
                result = 1
            else:
                #Target centers for self intersection calculation
                delta2 = np.abs(self.targetCtrs[targetIdx]-self.q1[i])[0]
                #check if intersection face is over 10mm from start face
                #which implies that it is not a self intersection
#                if (delta2[0] > 0.02) or (delta2[1] > 0.02) or (delta2[2] > 0.2):
                if np.sqrt(delta2[0]**2 + delta2[1]**2 + delta2[2]**2) > 0.2:
                    result = 1
                else:
                    result = 0

        return result


    def readLaminarParallel(self,i):
        """
        Parallel version of MAFOT laminar read.
        Sometimes the MAFOT calculation returns an error and it discards the
        launch point instead of writing to the output file.  This results
        in less points coming out of MAFOT than we put in.  So we amend this
        by checking to make sure the points coming out are assigned to the correct
        location in the mesh centers.
        """
        for j in range(len(self.lamData[:,0])):
            if self.lamR[i] == self.lamData[j,0]:
                break

        return j

    def readStructParallel(self,i):
        """
        Parallel version of MAFOT structure read.
        Sometimes the MAFOT calculation returns an error and it discards the
        launch point instead of writing to the output file.  This results
        in less points coming out of MAFOT than we put in.  So we amend this
        by checking to make sure the points coming out are assigned to the correct
        location in the mesh centers.
        """
        test=1
        for j in range(len(self.structData[:,3])):
            if self.R[i] == self.structData[j,3]:
                test=0

        return test


    def buildDirectories(self, PFCnames, timesteps, dataPath, clobberFlag=True):
        """
        builds a HEAT tree from timesteps and PFC part names
        """
        for t in timesteps:
            #Build timestep directory
            if dataPath[-1]!='/':
                timeDir = dataPath + '/{:06d}/'.format(t)
            else:
                timeDir = dataPath + '{:06d}/'.format(t)
            try:
                os.mkdir(timeDir)
                print("Tree Directory " , timeDir ,  " Created ")
            except FileExistsError:
                clobberFlagTime = False #don't overwrite time directories, just PFC directories
                if clobberFlagTime is True:
                    try: shutil.rmtree(timeDir)
                    except OSError as e:
                        print ("Error: %s - %s." % (e.filename, e.strerror))
                        sys.exit()

                    os.mkdir(timeDir)
                    print("Tree Directory " , timeDir ,  " Created ")

            #build directory for each PFC partname
            for name in PFCnames:
                pfcDir = timeDir + name
                try:
                    os.mkdir(pfcDir)
                    print("Tree Directory " , pfcDir ,  " Created ")
                except FileExistsError:
                    if clobberFlag is True:
                        try: shutil.rmtree(pfcDir)
                        except OSError as e:
                            print ("Error: %s - %s." % (e.filename, e.strerror))
                            return

                        os.mkdir(pfcDir)
                        print("Tree Directory " , pfcDir ,  " Created ")
        return

    def makeDir(self, dir, clobberFlag=True):
        """
        builds a directory with clobber checking.  for general use.
        """
        try:
            os.makedirs(dir)
            print("Directory " , dir ,  " Created ")
        except:
            if clobberFlag is True:
                try:
                    shutil.rmtree(dir)
                    os.makedirs(dir)
                    print("Directory " , dir ,  " Created ")
                except OSError as e:
                    print ("Error: %s - %s." % (e.filename, e.strerror))
        return


    def is_number(self, s):
        """
        checks if s is a number
        s is a string
        returns boolean (true if number)
        """
        try:
            float(s)
            return True
        except ValueError:
            return False

#toolsClass.py
#Description:   Utility tools for HEAT
#Engineer:      T Looby
#Date:          20191117
#
#
#THIS MODULE IS FROM HEAT, AND AS SUCH IS PROTECTED UNDER THE MIT LICENSE!
#USERS MUST ATTRIBUTE THE SOURCE CODE
#See https://github.com/plasmapotential/HEAT for more information

import sys
import pandas as pd
import numpy as np
import os
import shutil
import logging
import subprocess
import psutil
log = logging.getLogger(__name__)

class tools:
    """
    These are tools that are used by multiple modules in HEAT.  Stuff like
    setting up input file pipelines, generic coordinate transformations, etc.
    """

    def __init__(self, chmod=0o774):
        #speed benchmarking tools
        self.testN = 0
        self.testTime = 0
        self.chmod = chmod
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
                noneList = ['none', None, 'None', 'NONE', '', 'na', 'NA', 'Na', 'nan', 'NaN', 'NAN', 'Nan']
                #handle Nones
                if type(data['Val'][i]) == str:
                    if (data['Val'][i].strip() in noneList):
                        setattr(obj, data['Var'][i], None)
                    else:
                        setattr(obj, data['Var'][i], data['Val'][i].strip())
                else:
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
        import gyroClass
        MHD = MHDClass.MHD(rootDir, dataPath)
        CAD = CADClass.CAD(rootDir, dataPath)
        HF = heatfluxClass.heatFlux(rootDir, dataPath)
        OF = openFOAMclass.OpenFOAM(rootDir, dataPath)
        GYRO = gyroClass.GYRO(rootDir, dataPath)
        MHD.allowed_class_vars()
        CAD.allowed_class_vars()
        HF.allowed_class_vars()
        OF.allowed_class_vars()
        GYRO.allowed_class_vars()

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
                    if (data[var] == None) or (data[var] == 'None'):
                        f.write(var + ', None \n')
                    else:
                        f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ', None \n')
            f.write("#=============================================================\n")
            f.write("#                MHD Variables\n")
            f.write("#=============================================================\n")
            for var in MHD.allowed_vars:
                if var in data:
                    if (data[var] == None) or (data[var] == 'None'):
                        f.write(var + ',  None \n')
                    else:
                        f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ',  None \n')
            f.write("#=============================================================\n")
            f.write("#                Optical HF Variables\n")
            f.write("#=============================================================\n")
            for var in HF.allowed_vars:
                if var in data:
                    if (data[var] == None) or (data[var] == 'None'):
                        f.write(var + ',  None \n')
                    else:
                        f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ',  None \n')
            f.write("#=============================================================\n")
            f.write("#                Ion Gyro Orbit HF Variables\n")
            f.write("#=============================================================\n")
            for var in GYRO.allowed_vars:
                if var in data:
                    if (data[var] == None) or (data[var] == 'None'):
                        f.write(var + ',  None \n')
                    else:
                        f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ',  None \n')
            f.write("#=============================================================\n")
            f.write("#                Radiated Power HF Variables\n")
            f.write("#=============================================================\n")
            for var in RAD.allowed_vars:
                if var in data:
                    if (data[var] == None) or (data[var] == 'None'):
                        f.write(var + ',  None \n')
                    else:
                        f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ',  None \n')
            f.write("#=============================================================\n")
            f.write("#                OpenFOAM Variables\n")
            f.write("#=============================================================\n")
            print(data)
            for var in OF.allowed_vars:
                print(var)
                if var in data:
                    if (data[var] == None) or (data[var] == 'None'):
                        f.write(var + ',  None \n')
                    else:
                        f.write(var + ', ' + str(data[var]) + '\n')
                else:
                    f.write(var + ',  None \n')

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
            f.write('#PFC CSV FILE\n')
            f.write('#\n')
            f.write('#This is a HEAT input file responsible for:\n')
            f.write('# 1) Defining the PFCs in the Region of Interest (ROI)\n')
            f.write('# 2) Defining the heat flux calculation mesh resolution for each PFC in ROI\n')
            f.write('# 3) Defining the intersection PFCs that the ROI PFCs could be shadowed by\n')
            f.write('# 4) Defining the timesteps we run HEAT for, w.r.t. each PFC\n')
            f.write('# 5) Defining which divertor each PFC is in\n')
            f.write('# 6) Defining any PFCs we want to exclude from the intersection calculation\n')
            f.write('#\n')
            f.write('#The first line in this file should always be:\n')
            f.write('#timesteps, PFCname, resolution, DivCode, intersectName, excludeName\n')
            f.write('#\n')
            f.write('#Each new row under this header line represents a new PFC in the ROI.\n')
            f.write('#Definitions of each column (comma separated variable) from this header line:\n')
            f.write('#  timesteps: timesteps we include the PFC in during the HEAT calculation\n')
            f.write('#  PFCname: name of the PFC we want to calculate heat loads on.  Name should\n')
            f.write('#           be taken directly from the part in the CAD STP file.\n')
            f.write('#  resolution: heat flux resolution on PFCname.\n')
            f.write('#  DivCode:  divertor code.  This can be: LO, LI, UO, UI, which correspond to:\n')
            f.write('#            Lower Outer, Lower Inner, Upper Outer, Upper Inner.  These codes\n')
            f.write('#            are how each PFC in the ROI get flagged as belonging to a specific\n')
            f.write('#            divertor, and will affect their power later (Psol * powerFrac)\n')
            f.write('#  intersectName:  name of the PFCs that may cast a magnetic shadow upon the PFC\n')
            f.write('#                  we are calculating power on.  If we launch field lines from the\n')
            f.write('#                  PFCname part and follow them up the SOL, we may hit one of these\n')
            f.write('#                  intersectName PFCs.  If the user is unsure of this value, "all"\n')
            f.write('#                  can be specified to check against all parts in the CAD.\n')
            f.write('#                  Multiple PFCs can be specified by using ":" between part names.\n')
            f.write('#  excludeName:  name of PFCs to exclude in the intersection check.  This can be\n')
            f.write('#                useful when we use the "all" switch for intersectName and want\n')
            f.write('#                to exclude some obvious PFCs.\n')
            f.write('#\n')
            f.write('#Example to calculate heat loads on PFC001, a PFC in the Lower Outer divertor,\n')
            f.write('# at 5.0 mm heat flux mesh resolution, for the first 100ms of the discharge, \n')
            f.write('# checking for shadows (intersections) with PFC001 (self), PFC002 and PFC003:\n')
            f.write('# timesteps, PFCname, resolution, DivCode, intersectName, excludeName\n')
            f.write('# 0:100, PFC001, 5.0, LO, PFC001:PFC002:PFC003, none\n')
            f.write('#\n')
            f.write('#Example to calculate heat loads on PFC001, a PFC in the Upper Inner divertor,\n')
            f.write('# at 1.0 mm heat flux mesh resolution, for timesteps 100:200ms of the discharge, \n')
            f.write('# checking for shadows (intersections), with all PFCs in the CAD file:\n')
            f.write('# timesteps, PFCname, resolution, DivCode, intersectName, excludeName\n')
            f.write('# 100:200, PFC001, 1.0, UI, all, none\n')
            f.write('#\n')
            f.write('#\n')
            f.write('#\n')
            f.write('timesteps, PFCname, resolution, DivCode, intersectName, excludeName\n')

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

    def signedVolume2(self,a,b,c,d, ax=1):
        """
        Calculates signed volume
        NOTE: true signed volume would be multiplied by 1/6 but we don't
        do that here
        """
        #return (1.0/6.0)*np.diagonal(np.matmul(np.cross(b-a,c-a),(d-a).T))
        #return (1.0/6.0) * np.dot(np.cross(b-a,c-a), d-a)
        return np.sum(np.cross(b-a, c-a) * (d-a), axis=ax)

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

    def createVTKOutput(self, pcfile, outType, prefix, verbose=False):
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
        if verbose==True:
            print("Spawning PVpython subprocess")
            log.info("Spawning PVpython subprocess")
        args = [pvpythonCMD, self.rootDir + '/GUIscripts/csv2vtk.py', pcfile, outType, prefix]
        from subprocess import run
        run(args, env=current_env)
        if verbose==True:
            print("PVpython subprocess complete")
            log.info("PVpython subprocess complete")
        return

    def intersectionTestParallelKdTree(self,i):
        """
        intersection test that uses FreeCAD's interal kd-tree method
        not used - left here for reference
        acceleration structure in HEAT is toroidal / poloidal flux filtering
        """
        rayOrig = q1[i]
        rayTerm = q2[i]
        rayVec = rayTerm - rayOrig
        rayDist = np.linalg.norm(rayVec)
        rayDir = rayVec / rayDist
        intersect = tools.mesh.nearestFacetOnRay((rayOrig[0],rayOrig[1],rayOrig[2]),(rayDir[0],rayDir[1],rayDir[2]))

        #we found an intersection
        if bool(intersect):
            #get the index of the intersection
            result = list(intersect.keys())[0]
            #get the location of the intersection
            #loc = list(intersect.values())[0]
        #no intersection on this mesh
        else:
            result  = None

        return result

    def intersectTestParallelMT(self, i):
        """
        Intersection uses adapted version of Moller-Trumbore Algorithm:
        Möller, Tomas; Trumbore, Ben (1997). "Fast, Minimum Storage Ray-Triangle Intersection".
           Journal of Graphics Tools. 2: 21–28. doi:10.1080/10867651.1997.10487468.

        also see this stack overflow thread:
        https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d/42752998#42752998

        i is index of parallel run from multiprocessing

        q1 and q2 are start/end points of trace
        p1,p2,p3 are points of potential intersection faces

        """
        #backface culling
        if self.bfCull == True:
            if self.powerDir[i] < 0:
                use0 = self.targetsFwdUse
                p1 = self.p1Fwd
                p2 = self.p2Fwd
                p3 = self.p3Fwd
                if self.phiFilterSwitch == True:
                    phiP1 = self.phiP1[self.targetsFwdUse]
                    phiP2 = self.phiP2[self.targetsFwdUse]
                    phiP3 = self.phiP3[self.targetsFwdUse]
                if self.psiFilterSwitch == True:
                    psiP1 = self.psiP1[self.targetsFwdUse]
                    psiP2 = self.psiP2[self.targetsFwdUse]
                    psiP3 = self.psiP3[self.targetsFwdUse]

            else:
                use0 = self.targetsRevUse
                p1 = self.p1Rev
                p2 = self.p2Rev
                p3 = self.p3Rev
                if self.phiFilterSwitch == True:
                    phiP1 = self.phiP1[self.targetsRevUse]
                    phiP2 = self.phiP2[self.targetsRevUse]
                    phiP3 = self.phiP3[self.targetsRevUse]
                if self.psiFilterSwitch == True:
                    psiP1 = self.psiP1[self.targetsRevUse]
                    psiP2 = self.psiP2[self.targetsRevUse]
                    psiP3 = self.psiP3[self.targetsRevUse]
        #check all faces (no culling)
        else:
            use0 = np.arange(len(self.p1))
            p1 = self.p1
            p2 = self.p2
            p3 = self.p3
            if self.phiFilterSwitch == True:
                phiP1 = self.phiP1
                phiP2 = self.phiP2
                phiP3 = self.phiP3
            if self.psiFilterSwitch == True:
                psiP1 = self.psiP1
                psiP2 = self.psiP2
                psiP3 = self.psiP3

        #Filter by psi
        if self.psiFilterSwitch == True:
            psiMin = self.psiMin[i]
            psiMax = self.psiMax[i]

            #account for psi sign convention
            if psiMin > psiMax:
                pMin = psiMax
                pMax = psiMin
            else:
                pMin = psiMin
                pMax = psiMax

            #target faces outside of this poloidal slice
            test0 = np.logical_and(np.logical_and(psiP1 < pMin, psiP2 < pMin), psiP3 < pMin)
            test1 = np.logical_and(np.logical_and(psiP1 > pMax, psiP2 > pMax), psiP3 > pMax)
            test = np.logical_or(test0,test1)
            usePsi = np.where(test == False)[0]

        else:
            usePsi = np.arange(len(p1))

        #Filter by phi (toroidal angle)
        if self.phiFilterSwitch == True:
            phiMin = self.phiMin[i]
            phiMax = self.phiMax[i]

            #angle wrap cases (assumes we never trace in MAFOT steps larger than 5degrees)
            if np.abs(phiMin-phiMax) > np.radians(5):
                phiP1[phiP1<0] += 2*np.pi
                phiP2[phiP2<0] += 2*np.pi
                phiP3[phiP3<0] += 2*np.pi
                if phiMin < 0: phiMin+=2*np.pi
                if phiMax < 0: phiMax+=2*np.pi

            #account for toroidal sign convention
            if phiMin > phiMax:
                pMin = phiMax
                pMax = phiMin
            else:
                pMin = phiMin
                pMax = phiMax

            #target faces outside of this toroidal slice
            test0 = np.logical_and(np.logical_and(phiP1 < pMin, phiP2 < pMin), phiP3 < pMin)
            test1 = np.logical_and(np.logical_and(phiP1 > pMax, phiP2 > pMax), phiP3 > pMax)
            test = np.logical_or(test0,test1)
            usePhi = np.where(test == False)[0]

        else:
            usePhi = np.arange(len(p1))

        #combine filter algorithms
        useF = np.intersect1d(usePsi,usePhi)

        Nt = len(use0[useF])

        #Perform Intersection Test
        D = self.D[i] / np.linalg.norm(self.D, axis=1)[i]
        eps = 0.0

        #old method.  depends upon triangle vertex permutations.  left for reference
        #h = np.cross(D, self.E2[use0[useF]])
        #a = np.sum(self.E1[use0[useF]]*h, axis=1)
        #test1 = np.logical_and( a>-eps, a<eps) #ray parallel to triangle
        #with np.errstate(divide='ignore', invalid='ignore'):
        #    f=1.0/a #sometimes results in div by 0 (python warning)
        #    s = self.q1[i] - self.p1[use0[useF]]
        #    u = f * np.sum(s*h, axis=1)
        #    test2 = np.logical_or(u<0.0, u>1.0) #ray inside triangle
        #    q = np.cross(s,self.E1[use0[useF]])
        #    v = f*np.sum(D*q, axis=1)
        #    test3 =  np.logical_or(v<0.0, (u+v)>1.0) #ray inside triangle
        #    l = f*np.sum(self.E2[use0[useF]]*q, axis=1) #sometimes results in invalid mult (python warning)
        #    test4 = np.logical_or(l<0.0, l>self.Dmag[i]) #ray long enough to intersect triangle

        #New method: impervious to triangle vertex permutations
        #rearranges matrix algebra as defined in MT paper, using vector triple product
        #so that self.N is precalculated
        T = self.q1[i] - self.p1[use0[useF]]
        TcrossD = np.cross(T,D)
        a = -np.sum(self.N[use0[useF]]*D, axis=1)
        with np.errstate(divide='ignore', invalid='ignore'):
            f=1.0/a #sometimes results in div by 0 (python warning)
            #barycentric coordinates
            u = f * np.sum(self.E2[use0[useF]]*TcrossD, axis=1)
            v = -f * np.sum(self.E1[use0[useF]]*TcrossD, axis=1)
            t = f * np.sum(T*self.N[use0[useF]], axis=1)
            #now test if intersection point is inside triangle (true=outside)
            test1 = np.logical_and( a>-eps, a<eps) #ray parallel to triangle
            test2 = np.logical_or(u<0.0, u>1.0) #ray inside triangle
            test3 = np.logical_or(v<0.0, v>1.0) #ray inside triangle
            test4 = (u+v)>1.0 #ray inside triangle
            test5 = np.logical_or(t<0.0, t>self.Dmag[i]) #ray long enough to intersect triangle
            if np.sum(~np.any([test1,test2,test3,test4,test5], axis=0))>0:
                result = 1
            else:
                result = 0

        #print which face we intersect with if ptIdx is this i
        if self.ptIdx is not None:
            if self.ptIdx == i:
                #we assume first intersection in this array is the intersection
                targetIdx = np.where(np.any([test1,test2,test3,test4], axis=0)==False)[0]
                if len(targetIdx)>0:
                    self.targetIdx = targetIdx
                    print("Found ptIdx's intersection target vertices:")
                    print(self.p1[use0[useF[self.targetIdx]]])
                    print(self.p2[use0[useF[self.targetIdx]]])
                    print(self.p3[use0[useF[self.targetIdx]]])
                print("ptIdx's last trace step:")
                print(self.q1[i])
                print(self.q2[i])

        return result


    def intersectTestParallel(self, i):
        """
        Intersection test that does not check for self intersections

        i is index of parallel run from multiprocessing

        q1 and q2 are start/end points of trace
        p1,p2,p3 are points of potential intersection faces

        using Signed Volume line + triangle intersection rule
        """
        #backface culling
        if self.bfCull == True:
            if self.powerDir[i] < 0:
                use0 = self.targetsFwdUse
                p1 = self.p1Fwd
                p2 = self.p2Fwd
                p3 = self.p3Fwd
                if self.phiFilterSwitch == True:
                    phiP1 = self.phiP1[self.targetsFwdUse]
                    phiP2 = self.phiP2[self.targetsFwdUse]
                    phiP3 = self.phiP3[self.targetsFwdUse]
                if self.psiFilterSwitch == True:
                    psiP1 = self.psiP1[self.targetsFwdUse]
                    psiP2 = self.psiP2[self.targetsFwdUse]
                    psiP3 = self.psiP3[self.targetsFwdUse]

            else:
                use0 = self.targetsRevUse
                p1 = self.p1Rev
                p2 = self.p2Rev
                p3 = self.p3Rev
                if self.phiFilterSwitch == True:
                    phiP1 = self.phiP1[self.targetsRevUse]
                    phiP2 = self.phiP2[self.targetsRevUse]
                    phiP3 = self.phiP3[self.targetsRevUse]
                if self.psiFilterSwitch == True:
                    psiP1 = self.psiP1
                    psiP2 = self.psiP2
                    psiP3 = self.psiP3
        #check all faces (no culling)
        else:
            use0 = np.arange(len(self.p1))
            p1 = self.p1
            p2 = self.p2
            p3 = self.p3
            if self.phiFilterSwitch == True:
                phiP1 = self.phiP1
                phiP2 = self.phiP2
                phiP3 = self.phiP3
            if self.psiFilterSwitch == True:
                psiP1 = self.psiP1
                psiP2 = self.psiP2
                psiP3 = self.psiP3


        #Filter by psi
        if self.psiFilterSwitch == True:
            psiMin = self.psiMin[i]
            psiMax = self.psiMax[i]

            #account for psi sign convention
            if psiMin > psiMax:
                pMin = psiMax
                pMax = psiMin
            else:
                pMin = psiMin
                pMax = psiMax

            #target faces outside of this toroidal slice
            test0 = np.logical_and(np.logical_and(psiP1 < pMin, psiP2 < pMin), psiP3 < pMin)
            test1 = np.logical_and(np.logical_and(psiP1 > pMax, psiP2 > pMax), psiP3 > pMax)
            test = np.logical_or(test0,test1)
            usePsi = np.where(test == False)[0]

        else:
            usePsi = np.arange(len(p1))

        #Filter by phi (toroidal angle)
        if self.phiFilterSwitch == True:
            phiMin = self.phiMin[i]
            phiMax = self.phiMax[i]

            #angle wrap cases (assumes we never trace in MAFOT steps larger than 5degrees)
            if np.abs(phiMin-phiMax) > np.radians(5):
                phiP1[phiP1<0] += 2*np.pi
                phiP2[phiP2<0] += 2*np.pi
                phiP3[phiP3<0] += 2*np.pi
                if phiMin < 0: phiMin+=2*np.pi
                if phiMax < 0: phiMax+=2*np.pi

            #account for toroidal sign convention
            if phiMin > phiMax:
                pMin = phiMax
                pMax = phiMin
            else:
                pMin = phiMin
                pMax = phiMax

            #target faces outside of this toroidal slice
            test0 = np.logical_and(np.logical_and(phiP1 < pMin, phiP2 < pMin), phiP3 < pMin)
            test1 = np.logical_and(np.logical_and(phiP1 > pMax, phiP2 > pMax), phiP3 > pMax)
            test = np.logical_or(test0,test1)
            usePhi = np.where(test == False)[0]

        else:
            usePhi = np.arange(len(p1))

        #combine filter algorithms
        useF = np.intersect1d(usePsi,usePhi)

        Nt = len(use0[useF])

        #Perform Intersection Test
        q13D = np.repeat(self.q1[i,np.newaxis], Nt, axis=0)
        q23D = np.repeat(self.q2[i,np.newaxis], Nt, axis=0)
        sign1 = np.sign(self.signedVolume2(q13D,self.p1[use0[useF]],self.p2[use0[useF]],self.p3[use0[useF]]))
        sign2 = np.sign(self.signedVolume2(q23D,self.p1[use0[useF]],self.p2[use0[useF]],self.p3[use0[useF]]))
        sign3 = np.sign(self.signedVolume2(q13D,q23D,self.p1[use0[useF]],self.p2[use0[useF]]))
        sign4 = np.sign(self.signedVolume2(q13D,q23D,self.p2[use0[useF]],self.p3[use0[useF]]))
        sign5 = np.sign(self.signedVolume2(q13D,q23D,self.p3[use0[useF]],self.p1[use0[useF]]))
        test1 = (sign1 != sign2)
        test2 = np.logical_and(sign3==sign4,sign3==sign5)

        #print which face we intersect with if ptIdx is this i
        if self.ptIdx is not None:
            if self.ptIdx == i:
                #we assume first intersection in this array is the intersection
                targetIdx = np.where(np.any([test1,test2,test3,test4], axis=0)==False)[0]
                if len(targetIdx)>0:
                    self.targetIdx = targetIdx
                    print("Found ptIdx's intersection target vertices:")
                    print(self.p1[use0[useF[self.targetIdx]]])
                    print(self.p2[use0[useF[self.targetIdx]]])
                    print(self.p3[use0[useF[self.targetIdx]]])
                print("ptIdx's last trace step:")
                print(self.q1[i])
                print(self.q2[i])

        #return 1 if we intersected
        if np.sum(np.logical_and(test1,test2)) > 0:
            result = 1
        else:
            result = 0

        return result


    def intersectTestNoLoop(self):
        """
        Intersection test that does not check for self intersections

        creates matrices instead of parallel processing, but this can
        easily saturate RAM.  only use for small checks.  otherwise use MT algorithm

        q1 and q2 are start/end points of trace
        p1,p2,p3 are points of potential intersection faces

        using sigvol line + triangle intersection rule
        """
        Nt = len(use)

        q13D = np.repeat(self.q1[:,np.newaxis], Nt, axis=1)
        q23D = np.repeat(self.q2[:,np.newaxis], Nt, axis=1)

        sign1 = np.sign(self.signedVolume2(q13D,self.p1[use],self.p2[use],self.p3[use],ax=2))
        sign2 = np.sign(self.signedVolume2(q23D,self.p1[use],self.p2[use],self.p3[use],ax=2))
        sign3 = np.sign(self.signedVolume2(q13D,q23D,self.p1[use],self.p2[use],ax=2))
        sign4 = np.sign(self.signedVolume2(q13D,q23D,self.p2[use],self.p3[use],ax=2))
        sign5 = np.sign(self.signedVolume2(q13D,q23D,self.p3[use],self.p1[use],ax=2))
        test1 = (sign1 != sign2)
        test2 = np.logical_and(sign3==sign4,sign3==sign5)
        return np.logical_and(test1,test2)


    def intersectTestSingleRay(self):
        """
        Intersection for single ray (q1 and q2 each are a single point)

        creates matrices instead of parallel processing, but this can
        easily saturate RAM.  only use for small checks.  otherwise use MT algorithm

        q1 and q2 are start/end points of trace
        p1,p2,p3 are points of potential intersection faces

        using sigvol line + triangle intersection rule
        """
        q13D = np.repeat(self.q1[np.newaxis,:], self.Nt, axis=0)
        q23D = np.repeat(self.q2[np.newaxis,:], self.Nt, axis=0)

        sign1 = np.sign(self.signedVolume2(q13D,self.p1,self.p2,self.p3))
        sign2 = np.sign(self.signedVolume2(q23D,self.p1,self.p2,self.p3))
        sign3 = np.sign(self.signedVolume2(q13D,q23D,self.p1,self.p2))
        sign4 = np.sign(self.signedVolume2(q13D,q23D,self.p2,self.p3))
        sign5 = np.sign(self.signedVolume2(q13D,q23D,self.p3,self.p1))
        test1 = (sign1 != sign2)
        test2 = np.logical_and(sign3==sign4,sign3==sign5)
        return np.logical_and(test1,test2)


    def intersectTestParallel_selfCheck(self,i):
        """
        LEGACY FUNCTION: left here for reference

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

    def buildMask(self, source, target, thresh=1):
        """
        builds a mask matrix of size (Ntarget, Nsource)
        elements are 1 if abs(source - target) < thresh
        otherwise 0
        """
        sourceMat = np.repeat(source[np.newaxis,:], len(target), axis=0)
        targetMat = np.repeat(target[:,np.newaxis], len(source), axis=1)
        mask = np.abs(sourceMat-targetMat)<thresh
        print("PSIMASK SHAPE: ")
        print(mask.shape)
        return mask

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


    def buildDirectories(self, PFCnames, timesteps, dataPath, clobberFlag=True, chmod=0o774, UID=-1, GID=-1):
        """
        builds a HEAT tree from timesteps and PFC part names
        """
        for t in timesteps:
            #Build timestep directory
            if dataPath[-1]!='/':
                timeDir = dataPath + '/{:06d}/'.format(t)
            else:
                timeDir = dataPath + '{:06d}/'.format(t)

            #don't overwrite time directories, just PFC directories
            self.makeDir(timeDir, clobberFlag=False, mode=chmod, UID=UID, GID=GID)

            #build directory for each PFC partname
            for name in PFCnames:
                pfcDir = timeDir + name
                #overwrite PFC directories
                self.makeDir(pfcDir, clobberFlag=True, mode=chmod, UID=UID, GID=GID)
            #set tree permissions
            self.recursivePermissions(timeDir, UID, GID, chmod)
        return

    def makeDir(self, path, clobberFlag=True, mode=None, UID=None, GID=None):
        """
        builds a directory with clobber checking.  for general use.
        user can pass an octal Linux permissions mode, and
        group ID #

        note: if this creates nested dirs only the path will have correct permissions!
              all dirs created above path will be 0o777 and $USER:$USER
        """
        #Make directory
        try:
            os.makedirs(path)
            print("Directory " , path ,  " Created ")
        except:
            if clobberFlag is True:
                try:
                    shutil.rmtree(path)
                    os.makedirs(path, mode=self.chmod)
                    print("Directory " , path ,  " clobbered and created ")
                except OSError as e:
                    print ("Error: %s - %s." % (e.filename, e.strerror))

        #change permissions
        if mode != None:
            try:
                os.chmod(path, mode)
            except OSError as e:
                print("Could not change directory permissions")
                print ("Error: %s - %s." % (e.filename, e.strerror))

        #change ownership
        if GID == None:
            GID = -1
        if UID == None:
            UID = -1
        try:
            os.chown(path, UID, GID)
        except OSError as e:
            print("Could not change directory ownership")
            print ("Error: %s - %s." % (e.filename, e.strerror))

        return

    def recursivePermissions(self, path, UID, GID, chmod):
        """
        recursively set permissions in a dir
        is identical to bash commands:
        chown -R <UID>:<GID> path
        chmod -R <chmod> path

        UID is user id
        GID is group id
        chmod is permissions in base 10 (or if in base 8 preface with 0o)
        """
        for dirpath, dirnames, filenames in os.walk(path):
            os.chown(dirpath, UID, GID)
            os.chmod(dirpath, chmod)
            for filename in filenames:
                os.chown(os.path.join(dirpath, filename), UID, GID)
                os.chmod(os.path.join(dirpath, filename), chmod)
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

    def readStructOutput(self,file):
        """
        Reads output file from MAFOT structure program
        """
        structdata = np.genfromtxt(file,comments='#')
        xyz = np.zeros((len(structdata),3))
        xyz[:,0] = structdata[:,0]
        xyz[:,1] = structdata[:,1]
        xyz[:,2] = structdata[:,2]
        #remove rows with zeros (invalids)
        #xyzFinal = xyz[~(np.abs(xyz)<1e-99).any(axis=1)]
        print('Read structure output for {:d} points'.format(len(xyz)))
        log.info('Read structure output for {:d} points'.format(len(xyz)))
        return xyz

    def getTargetCenters(self, targets):
        """
        returns target centers

        targets is 3 points [[p1],[p2],[p3]] that comprise a mesh triangle
        where [pN] = [xN,yN,zN]
        """
        p1 = targets[:,0,:]  #point 1 of mesh triangle
        p2 = targets[:,1,:]  #point 2 of mesh triangle
        p3 = targets[:,2,:]  #point 3 of mesh triangle
        Nt = len(p1)

        x = np.zeros((Nt,3))
        y = np.zeros((Nt,3))
        z = np.zeros((Nt,3))

        x[:,0] = p1[:,0]
        x[:,1] = p2[:,0]
        x[:,2] = p3[:,0]
        y[:,0] = p1[:,1]
        y[:,1] = p2[:,1]
        y[:,2] = p3[:,1]
        z[:,0] = p1[:,2]
        z[:,1] = p2[:,2]
        z[:,2] = p3[:,2]
        return self.faceCenters(x,y,z)

    def makeFloat(self, var):
        """
        converts var to float
        if var is None, returns None
        """
        if var == None:
            newVar = None
        else:
            try:
                newVar = float(var)
            except:
                newVar = var
        return newVar

    def makeInt(self, var):
        """
        converts var to int
        if var is None, returns None
        """
        if var == None:
            newVar = None
        else:
            try:
                newVar = int(float(var))
            except:
                newVar = var
        return newVar

    def VVdistortion(self, points):
        """
        distorts a mesh's xyz coordinates by applying a transform

        meshes is a list of meshes to distort
        deltaR is change in radius bounds
        deltaB is change in conical bounds
        N is toroidal mode number
        h is reference height
        R0 is reference radius
        """

        N = self.distortN
        R0 = self.distortR0
        deltaR = self.distortDeltaR
        deltaB = self.distortDeltaB
        h = self.distortH

        theta = np.arctan2(points[:,1], points[:,0])
        #distortion transforms
        xDist = np.sin(N*theta)*deltaR/R0 + deltaB*points[:,2]/h + 1
        yDist = np.cos(N*theta)*deltaR/R0 + deltaB*points[:,2]/h + 1
        zDist = np.ones((len(points)))

        #perform transform
        points[:,0] *= xDist
        points[:,1] *= yDist
        points[:,2] *= zDist
        return points

    def faceNormals(self, points):
        """
        returns an array of normal vectors from an array of triangle vertex points
        """
        norms = np.zeros((len(points),3))
        for i,pt in enumerate(points):
            A = pt[1] - pt[0]
            B = pt[2] - pt[0]
            norm = np.cross(A,B)
            norms[i,:] = norm / np.linalg.norm(norm)

        return norms

    def checkSignOfNorm(self, normDist, normTrue):
        """
        checks direction of Na against Nb and changes normDist to match normTrue

        this is for distorting the mesh where we know normTrue is the true surface
        normal and we dont want to accidentally flip the normals when we distort
        """

        dot = np.multiply(normDist,normTrue).sum(axis=1)
        flipLoc = np.where(dot < 0.0)[0]
        normDist[flipLoc,:] *= -1.0

        return normDist

    def calculatePowerDir(self, bdotn, Bt0):
        """
        calculate the power flow direction

        returns a vector of length bdotn than has the toroidal direction
        power flows (-1 is clockwise from above)
        """
        return np.sign(bdotn)*np.sign(Bt0)*-1.0

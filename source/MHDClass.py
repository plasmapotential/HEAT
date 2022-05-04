#MHDClass.py
#Description:   Base HEAT MHD Equilibrium (EQ) module
#Engineer:      T Looby
#Date:          20191111

import EFIT.equilParams_class as EP
import toolsClass
tools = toolsClass.tools()

import time
import numpy as np
import pandas as pd
import os
import shutil
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d

import logging
log = logging.getLogger(__name__)

def setupForTerminalUse(gFile=None, shot=None, time=None):
    """
    Sets up an MHD object so that it can be used from python console
    without running HEAT.  This is convenient when a user wants to load
    MHD EQ directly from the python console.

    To use, user needs to set pythonpath to include path of this file and EFIT:
    import sys
    EFITPath = '/home/tom/source'
    HEATPath = '/home/tom/source/HEAT/github/source'
    sys.path.append(EFITPath)
    sys.path.append(HEATPath)
    import MHDClass

    Then user needs to run this function to set up the MHD object:
    MHD = MHDClass.setupForTerminalUse()

    if user wants to load MHD EQ, they need to set gFile
    arguments when calling this function, otherwise it just returns an
    empty MHD object:
    MHD = MHDClass.setupForTerminalUse(gFile=gFilePath)

    if you want access some gFile parameters, use the ep.g:
    Example to see available parameters:
    MHD.ep.g.keys()
    Example to see Rlim, Zlim:
    MHD.ep.g['wall']



    returns an MHD class object ie: MHD = MHDClass.MHD()

    """
    rootDir = ''
    dataPath = ''

    import MHDClass
    MHD = MHDClass.MHD(rootDir,dataPath)


    if gFile!=None:
        print("Making MHD() object with ep included")
        MHD.ep = EP.equilParams(gFile)

    else:
        print("Not including ep in MHD() object")
    return MHD


class MHD:

    def __init__(self, rootDir, dataPath, chmod=0o774, UID=-1, GID=-1):
        """
        rootDir is root location of python modules (where dashGUI.py lives)
        dataPath is the location where we write all output to
        """
        self.rootDir = rootDir
        tools.rootDir = self.rootDir
        self.dataPath = dataPath
        self.chmod = chmod
        self.GID = GID
        self.UID = UID
        tools.dataPath = self.dataPath
        return

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization

        Here is a list of variables with description:
        testvar         dummy for testing
        MachFlag        Machine flag (nstx, d3d, cmod, etc.)
        shot            discharge number
        tree            EFIT tree (efit01, efit02, etc.)
        tmin            minimum timestep to consider
        tmax            maximum timestep to consider
        MHDpath         location where we will save / read gfiles
        """
        self.allowed_vars = [
                            #'MachFlag',
                            'shot',
                            'tmin',
                            'tmax',
                            'nTrace',
                            'dataPath',
                            'torFilt',
                            'psiFilt'
                            ]

        return


    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        integers = [
                    'shot',
                    'tmin',
                    'tmax',
                    'nTrace'
                    ]
        #floats = [
        #         'MapDirection',
        #         'MapDirectionStruct',
        #        ]
        bools = ['torFilt', 'psiFilt']

        for var in integers:
            if (getattr(self, var) is not None) and (~np.isnan(float(getattr(self, var)))):
                try:
                    setattr(self, var, tools.makeInt(getattr(self, var)))
                except:
                    print("Error with input file var "+var+".  Perhaps you have invalid input values?")
                    log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")

        trueStrings = ["T", "True", "t", "true", "TRUE", "y", "yes", "Yes", "YES", True]
        falseStrings =["F", "False", "f", "false", "FALSE", "n", "no", "No", "NO", False]
        for var in bools:
            if (getattr(self, var) is not None):
                try:
                    if getattr(self, var) in trueStrings:
                        setattr(self, var, True)
                    else:
                        setattr(self, var, False)
                except:
                    print("Error with input file var "+var+".  Perhaps you have invalid input values?")
                    log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")

        return


    def getGEQDSK(self,machine='nstx',gFileList=None):
        """
        get gfile from mds+ tree if gfileList is None, otherwise use file
        """
        #load gfile from MDS+
        if gFileList is None:
            self.timesteps = gfiles.GEQDSKFromMDS(
                                                  self.MachFlag,
                                                  self.shot,
                                                  self.tree,
                                                  self.tmin,
                                                  self.tmax,
                                                  self.shotPath,
                                                  clobber=True,
                                                  chmod=self.chmod,
                                                  GID=self.GID,
                                                  UID=self.UID
                                                  )
        #load from file in GUI or TUI
        else:
            self.timesteps = []
            if len(gFileList) > 1:
                self.N_gFiles = len(gFileList)
            for i,gfile in enumerate(gFileList):
                #test if gfile follows d3d naming convention
                try:
                    idx = gfile.find('.')
                    fmtstr = '0' + str(idx-1) + 'd'
                    shot, time = int(gfile[1:idx]), int(gfile[idx+1::])
                #otherwise define new timestep based upon index
                except:
                    shot=1
                    time=i+1
                ts = self.copyGfile2tree(gfile,shot,time)
                self.timesteps.append(ts)
            self.timesteps = np.array(self.timesteps)

        return


    def make1EFITobject(self,shot,t,gfile=None):
        """
        Creates an equilParams_class object for a SINGLE timesteps. equilParams
        is a class from the ORNL_Fusion github repo, and was developed by
        A. Wingen
        """
        if self.shotPath[-1]!='/':
            shotPath = self.shotPath+'/'
        else:
            shotPath = self.shotPath
        if gfile is None:
            gfile = shotPath+'{:06d}/g{:06d}.{:05d}'.format(t,shot,t)
        self.ep = EP.equilParams(gfile, shot, t)
        #Correct for weird helicity and field directions that are occasionally
        #in gfile.  Here we assume that Ip is the CCW direction as viewed from
        #above tokamak.
#        self.dsign = np.sign(self.ep.g['Ip'])
#        self.gsign = np.sign( (self.ep.g['psiAxis'] - self.ep.g['psiSep']) )
#        self.qsign = np.sign(self.ep.g['Fpol'][-1]) #F_edge sign
#        print('dsign = {:f}'.format(self.dsign))
#        print('gsign = {:f}'.format(self.gsign))
#        print('qsign = {:f}'.format(self.qsign))
#        self.ep.g['FFprime'] *= self.dsign*self.qsign*self.gsign
#        self.ep.g['Pprime'] *= self.dsign*self.qsign*self.gsign
#        self.ep.g['psiRZ'] *= self.dsign*self.qsign*self.gsign
#        self.ep.g['qpsi'] *= -self.qsign*self.dsign
#        self.ep.g['psiSep'] *= self.dsign*self.qsign*self.gsign
#        self.ep.g['psiAxis'] *= self.dsign*self.qsign*self.gsign
#        self.ep.g['Fpol'] *= self.dsign*self.qsign
        return

    def makeEFITobjects(self):
        """
        Creates an equilParams_class object for MULTIPLE timesteps. equilParams
        is a class from the ORNL_Fusion github repo, and was developed by
        A. Wingen

        gfiles should be placed in the dataPath before running this function
        """
        self.ep= ['None' for i in range(len(self.timesteps))]
        for idx,t in enumerate(self.timesteps):
            if self.shotPath[-1]!='/':
                shotPath = self.shotPath+'/'
            else:
                shotPath = self.shotPath
            gfile = shotPath+'{:06d}/g{:06d}.{:05d}'.format(t,self.shot,t)
            self.ep[idx] = EP.equilParams(gfile)
        return


    def Bfield_pointcloud(self, ep, R, Z, phi, powerDir=None, normal=False):
        """
        Creates a Bfield Pointcloud that can be saved in a .csv file
        This pointcloud is then used by ParaView to overlay the
        Bfield onto the CAD for analysis.  The format is:
        X,Y,Z,Bx,By,Bz
        where X,Y,Z are the spatial locations of the Bx,By,Bz vectors
        In paraview use TableToPoints => Calculator => Glyph
        Calculator should have this formula:
        (iHat*Bx) + (jHat*By) + (kHat*Bz)

        powerDirection is the direction power flows onto the PFC, and can be assigned by
        user if desired.  Otherwise, it is calculated in PFC class for each mesh
        element

        Note: this function returns the normal Bfield components if normal=True

        """
#        if powerDir != None:
#            Bt0Direction = np.sign(ep.g['Bt0']) #account for machine Bt direction
#            Bt = ep.BtFunc.ev(R,Z)*powerDir*Bt0Direction
#            BR = ep.BRFunc.ev(R,Z)*powerDir*Bt0Direction
#            BZ = ep.BZFunc.ev(R,Z)*powerDir*Bt0Direction
#        else:
        Bt = ep.BtFunc.ev(R,Z)
        BR = ep.BRFunc.ev(R,Z)
        BZ = ep.BZFunc.ev(R,Z)
        #print("TEST=====")
        #print(BR[0])
        #print(Bt[0])
        #print(BZ[0])
        #print(max(BR))
        #print(max(Bt))
        #print(max(BZ))
        try:
            Bxyz = np.zeros((len(BR),3))
        except:
            Bxyz = np.zeros((1,3))
        Bxyz[:,0] = (BR[:]*np.cos(phi[:]) - Bt[:]*np.sin(phi[:]))
        Bxyz[:,1] = (BR[:]*np.sin(phi[:]) + Bt[:]*np.cos(phi[:]))
        Bxyz[:,2] = BZ[:].copy()
        B = np.sqrt(Bxyz[:,0]**2 + Bxyz[:,1]**2 + Bxyz[:,2]**2)
        if normal:
            Bxyz[:,0] /= B
            Bxyz[:,1] /= B
            Bxyz[:,2] /= B
        return Bxyz

    def B_pointclouds(self, ep, R, Z):
        """
        returns a 1D vector of Bp, Bt, Br, Bz, values that correspond to
        input R,Z points
        """
        Br = ep.BRFunc.ev(R,Z)
        Bz = ep.BZFunc.ev(R,Z)
        Bt = ep.BtFunc.ev(R,Z)
        Bp = np.sqrt(Br**2+Bz**2)
        return Bp, Bt, Br, Bz


    def write_Bvec_pointcloud(self,centers,Bxyz, dataPath, tag=None):
        """
        In paraview use TableToPoints => Calculator => Glyph
        Calculator should have this formula:
        (iHat*Bx) + (jHat*By) + (kHat*Bz)
        """
        print("Creating B Field Point Cloud")
        log.info("Creating B Field Point Cloud")
        if tag is None:
            pcfile = dataPath + 'B_pointcloud.csv'
        else:
            pcfile = dataPath + 'B_pointcloud_'+tag+'.csv'
        pc = np.zeros((len(centers), 6))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = Bxyz[:,0]
        pc[:,4] = Bxyz[:,1]
        pc[:,5] = Bxyz[:,2]
        head = "X,Y,Z,Bx,By,Bz"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)
        #np.savetxt('points.asc', pc[:,:-1], delimiter=' ',fmt='%.10f')
        #print("Wrote point cloud file: " + pcfile)
        if tag is None:
            tools.createVTKOutput(pcfile, 'glyph', 'Bdotpower_pointcloud')
        else:
            name = 'B_pointcloud_'+tag
            tools.createVTKOutput(pcfile, 'glyph', name)
        return

    def write_B_pointclouds(self,centers,Bp,Bt,Br,Bz,dataPath, tag=None):
        """
        writes a point cloud of poloidal and toroidal field
        """
        print("Creating B Point Clouds")
        log.info("Creating B Point Clouds")

        prefixes = ['Bp','Bt','Br','Bz']
        arrays = [Bp, Bt, Br, Bz]

        for i,prefix in enumerate(prefixes):
            if tag is None:
                pcfile = dataPath + prefix + '.csv'
            else:
                pcfile = dataPath + prefix + '_'+tag+'.csv'
            print(pcfile)
            pc = np.zeros((len(centers), 4))
            pc[:,0] = centers[:,0]*1000.0
            pc[:,1] = centers[:,1]*1000.0
            pc[:,2] = centers[:,2]*1000.0
            pc[:,3] = arrays[i]
            head = "X,Y,Z,"+prefix
            np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

            #Now save a vtk file for paraviewweb
            if tag is None:
                tools.createVTKOutput(pcfile, 'points', prefix)
            else:
                name = prefix+'_'+tag
                tools.createVTKOutput(pcfile, 'points', name)


        return


    def writeControlFile(self, name, t, mapDirection, mode='laminar'):
        """
        Create a control file for MAFOT
        """
        if len(name.split('/')) > 1:
            save_location = name
        else:
            if self.shotPath[-1] == '/':
                save_location = self.shotPath + '{:06d}/'.format(t) + name
            else:
                save_location = self.shotPath + '/{:06d}/'.format(t) + name
        with open(save_location, 'w') as f:

            f.write('# Parameterfile for ' + self.MachFlag + ' Programs\n')
            f.write('# Shot: {:06d}\tTime: {:05d}ms\n'.format(int(self.shot), int(t)))
            if self.shotPath[-1] == '/':
                f.write('# Path: ' + self.shotPath + '{:06d}\n'.format(t))
            else:
                f.write('# Path: ' + self.shotPath + '/{:06d}\n'.format(t))

            f.write('Nphi=\t{:d}\n'.format(self.Nphi))

            #itt means different things depending on if we are tracing field line
            #or running full MAFOT laminar
            if mode=='laminar':
                f.write('itt=\t{:f}\n'.format(self.ittLaminar))
            elif mode=='gyro':
                f.write('itt=\t{:f}\n'.format(self.ittGyro))
            else:
                f.write('itt=\t{:f}\n'.format(self.ittStruct))
            #f.write('Smin=\t{:2f}\n'.format(self.Smin))
            #f.write('Smax=\t{:2f}\n'.format(self.Smax))
            f.write('Rmin=\t{:2f}\n'.format(self.Rmin))
            f.write('Rmax=\t{:2f}\n'.format(self.Rmax))
            f.write('Zmin=\t{:2f}\n'.format(self.Zmin))
            f.write('Zmax=\t{:2f}\n'.format(self.Zmax))
            #f.write('phimin=\t{:2f}\n'.format(self.phimin))
            #f.write('phimax=\t{:2f}\n'.format(self.phimax))
            f.write('Nswall=\t{:d}\n'.format(self.Nswall))

            f.write('phistart(deg)=\t{:2f}\n'.format(self.phistart))
            f.write('MapDirection=\t{:f}\n'.format(mapDirection))
            #We check here to see if we defined a multiplier for MAFOT trace direction
            #because MAFOT assumes increasing monotonic psiN (cant be decreasing)
#            if (self.structMapDirMultiply >= 0.0) or (self.structMapDirMultiply is None):
#                f.write('MapDirection=\t{:f}\n'.format(mapDirection))
#                print("Writing CTL file with mapDir = {:f}".format(mapDirection))
#                log.info("Writing CTL file with mapDir = {:f}".format(mapDirection))
#            else:
#                f.write('MapDirection=\t{:f}\n'.format(mapDirection*-1.0))
#                print("Writing CTL file with mapDir = {:f}".format(mapDirection*-1.0))
#                log.info("Writing CTL file with mapDir = {:f}".format(mapDirection*-1.0))
            f.write('PlasmaResponse(0=no,>1=yes)=\t{:d}\n'
                    .format(self.PlasmaResponse))
            f.write('Field(-3=VMEC,-2=SIESTA,-1=gfile,M3DC1:0=Eq,1=I-coil,2=both)=\t'
                    '{:d}\n'.format(self.Field))

            f.write('target(0=useSwall)=\t{:d}\n'.format(self.target))
            f.write('createPoints(2=target)=\t{:d}\n'.format(self.createPoints))

            if(self.MachFlag == 'iter'):
                f.write('useIcoil(0=no,1=yes)=\t{:d}\n'.format(self.useIcoil))
            elif(self.MachFlag == 'nstx'):
                f.write('useECcoil(0=no,1=yes)=\t{:d}\n'.format(self.useECcoil))
            elif(self.MachFlag == 'mast'):
                f.write('useCcoil(0=no,1=yes)=\t{:d}\n'.format(self.useCcoil))
                f.write('useIcoil(0=no,1=yes)=\t{:d}\n'.format(self.useIcoil))
            elif(self.MachFlag == 'd3d'):
                f.write('useFcoil(0=no,1=yes)=\t{:d}\n'.format(self.useFcoil))
                f.write('useCcoil(0=no,1=yes)=\t{:d}\n'.format(self.useCcoil))
                f.write('useIcoil(0=no,1=yes)=\t{:d}\n'.format(self.useIcoil))
            else:
                f.write('useECcoil(0=no,1=yes)=\t{:d}\n'.format(self.useECcoil))

            if self.MachFlag in self.machineList:
                f.write('useFilament(0=no)=\t{:d}\n'.format(self.useFilament))
                f.write('useTe_profile(0=no)=	{:d}\n'.format(self.useTe_profile))

            f.write('ParticleDirection(1=co-pass,-1=ctr-pass,0=field-lines)=\t{:d}\n'
                    .format(self.ParticleDirection))
            f.write('PartileCharge(-1=electrons,>=1=ions)=\t{:d}\n'
                    .format(self.ParticleCharge))
            f.write('Ekin[keV]=\t{:2f}\n'.format(self.Ekin))
            f.write('lambda=\t{:2f}\n'.format(self.Lambda))
            f.write('Mass=\t{:2f}\n'.format(self.Mass))

            if self.MachFlag in ['dt']:
                f.write('useFilament(0=no)=\t{:d}\n'.format(self.useFilament))
                f.write('useBusError(0=no,1=yes)=\t{:d}\n'.format(self.useBus))
                f.write('useBcoilError(0=no,1=yes)=\t{:d}\n'.format(self.useBcoil))

            f.write('pi=\t3.141592653589793\n')
            f.write('2*pi=\t6.283185307179586\n')
            return

    def psi2DfromEQ(self, PFC):
        """
        Returns psi from EFIT equilibrium rather than from 3D trace.  Does not
        use MAFOT
        """
        use = np.where(PFC.shadowed_mask == 0)[0]
        xyz = PFC.centers[use]
        R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        PFC.psimin = PFC.ep.psiFunc.ev(R,Z)
        return

    def psi2DfromEQandCtrs(self, xyz, ep):
        """
        Returns psi from EFIT equilibrium rather than from 3D trace.  Does not
        use MAFOT
        """
        R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        psi = ep.psiFunc.ev(R,Z)
        return psi

    def writeMAFOTpointfile(self,xyz,gridfile):
        """
        Creates a file that MAFOT's can use to trace a field line.  Points
        in this file are the starting points for traces, organized in columns
        of R,phi,Z

        input should be in meters
        """
        #Convert to meters
        #xyz/=1000.0
        if len(xyz.shape) > 1:
            R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        else:
            R,Z,phi = tools.xyz2cyl(xyz[0],xyz[1],xyz[2])


        phi = np.degrees(phi)
        #Save R,phi,Z into a special initial condition grid file for MAFOT
        array = np.column_stack((R,phi,Z))
        np.savetxt(gridfile, array, delimiter='\t',fmt='%.10f')
        return



    def getFieldpath(self, dphi, gridfile, controlfilePath,
                    controlfile, paraview_mask=False, tag=None):
        """
        Uses MAFOT's structure script to get the full magnetic field line
        starting from a user defined R,Z,phi location.  Integrates back up
        the field line in dphi (integer) increments.
        machine is a string name of machine 'nstx'.  Controlfile contains
        parameters for MAFOT. Scale is unit conversion scale (ie mm => m is 1000)
        If paraview=True, writes CSV output for paraview

        After this is completed you need to open this file in ParaView then
        do the following filters:
        1) Table to Points
        2) Programmable Filters (scale by 1000,1000,-1000 to account for [m]=>[mm])
        The programmable filter should run the following code:
            pdi = self.GetPolyDataInput()
            pdo =  self.GetPolyDataOutput()
            numPoints = pdi.GetNumberOfPoints()
            pdo.Allocate()
            for i in range(0, numPoints-1):
                points = [i, i+1]
                # VTK_LINE is 3
                pdo.InsertNextCell(3, 2, points)
        """
        t0 = time.time()
        #Create MAFOT shell command
        #MAFOT now runs a program called heatstructure that is generic for all machines
        args = []
        #MAFOT now runs a program called HEAT that is generic for all machines
        args.append('heatstructure')
        #args 1,2 are the number of degrees we want to run the trace for
        args.append('-d')
        args.append(str(dphi))
        #args 3,4 are the points that we launch traces from
        args.append('-P')
        args.append(gridfile)
        #args 5 is the MAFOT control file
        args.append(controlfile)
        #Copy the current environment (important when in appImage mode)
        current_env = os.environ.copy()
        #run MAFOT structure for points in gridfile
        from subprocess import run
        run(args, env=current_env, cwd=controlfilePath)
        run(['rm', 'log_*'], env=current_env, cwd=controlfilePath)

        if paraview_mask:
            #This file is not in the format we want, so we read it in, and then
            #write it out as a CSV.  Because this is a serial function that is only
            #called for special cases, we arent too concerned about speed or memory
            #so this works instead of making a more efficient pipeline
            readfile = controlfilePath + 'struct.dat'
            structdata = np.genfromtxt(readfile,comments='#')

            xyz = np.zeros((len(structdata),3))
            xyz[:,0] = structdata[:,0]
            xyz[:,1] = structdata[:,1]
            xyz[:,2] = structdata[:,2]
            xyz*=1000 #Scale for ParaVIEW
            xyz = xyz[~(np.abs(xyz[:,:2])<1e-50).any(axis=1)]
            head = 'X[mm],Y[mm],Z[mm]'
            #head = 'X[m],Y[m],Z[m],R[m],phi[rad]'
            if tag == None:
                outfile = controlfilePath+'struct.csv'
                vtkName = 'Field_trace'
            else:
                outfile = controlfilePath+'struct_'+tag+'.csv'
                vtkName = 'Field_trace_'+tag
            np.savetxt(outfile, xyz, delimiter=',', header=head)
            #Now save a vtk file for paraviewweb
            tools.createVTKOutput(outfile, 'trace', vtkName)
            print('Converted file to ParaView formatted CSV.')
            log.info('Converted file to ParaView formatted CSV.')

        return

    def getMultipleFieldPaths(self, dphi, gridfile, controlfilePath,
                                controlfile, paraview_mask=False):

        args = []
        #args 0 is MAFOT structure binary call
        #MAFOT now runs a program called HEAT that is generic for all machines
        args.append('heatstructure')
        #args 1,2 are the number of degrees we want to run the trace for
        args.append('-d')
        args.append(str(dphi))
        #args 3,4 are the points that we launch traces from
        args.append('-P')
        args.append(gridfile)
        #args 5 is the MAFOT control file
        args.append(controlfile)
        #Copy the current environment (important when in appImage mode)
        current_env = os.environ.copy()
        #run MAFOT structure for points in gridfile
        from subprocess import run
        run(args, env=current_env, cwd=controlfilePath)
        return

    def runMAFOTlaminar(self,gridfile,controlfilePath,controlfile,NCPUs):
        """
        Runs MAFOT laminar program to trace field lines and return connection
        length, psi, penetration depth, etc.
        """
        #Create MAFOT shell command
        #MAFOT now runs a program called HEAT that is generic for all machines
        tool = 'heatlaminar_mpi'
        cmd = 'mpirun'
        args = [cmd,'-n','{:d}'.format(NCPUs),tool,controlfile,'-P',gridfile]
        #Copy the current environment (important when in appImage mode)
        current_env = os.environ.copy()
        #run MAFOT structure for points in gridfile
        from subprocess import run
        run(args, env=current_env, cwd=controlfilePath)

        #you rewrote with the above subprocess.run method but never tested
        #this is what it was before:
        #from subprocess import call
        #flags = ' -n {:d} '.format(NCPUs) + tool + ' ' + controlfile + ' ' + '-P '+gridfile
        #call(cmd + flags, shell = True, cwd=controlfilePath)
        #call('rm log_*', shell = True, cwd=controlfilePath)	# clean up
        #os.remove(controlfilePath+'log_*') #clean up

        return

    def renormalizeLCFS(self, ep, rNew=None, zNew=None, psiSep=None):
        """
        script for changing the lcfs in gfile.  because we use CAD in HEAT, and
        this CAD is oftentimes at different r,z locations than the rlim,zlim from
        gfile, we sometimes need to adjust the gfile to reflect the CAD.  This
        is especially true for limited discharges, where if we use the rlim, zlim
        wall from the gfile and try to compute qDiv on a real CAD PFC the center
        stack CAD PFCs are sometimes in the core with psiN < 1.  This script
        updates the LCFS in the gfile so that only one point in the CAD touches
        the LCFS.

        User may supply either:
        -R coord with no Z coord, then function finds min psi for all Z at that R
        -R coord and Z coord, then function finds psi at that point
        -psiSep, then function finds flux surface at that psiSep
        """


        g = ep.g
        if zNew is not None:
            print("Redefine LCFS using point")
            log.info("Redefine LCFS using point")
            psiNew = ep.psiFunc_noN.ev(rNew, zNew)
        elif rNew is not None:
            print("Redefine LCFS using vector")
            log.info("Redefine LCFS using vector")
            zMin = g['ZmAxis'] - 0.25
            zMax = g['ZmAxis'] + 0.25
            zWall = np.linspace(zMin, zMax, 100000)
            rWall = np.ones((len(zWall)))*rNew
#            psiWall = ep.psiFunc_noN.ev(rWall, zWall)
#            psiNew = ep.psiFunc_noN.ev(rWall, zWall).min()
            psiaxis = g['psiAxis']
            psiedge = g['psiSep']
            psiNwall = ep.psiFunc.ev(rWall, zWall)
            psiWall = psiNwall*(psiedge - psiaxis) + psiaxis
            psiNew = psiWall.min()
            psiNewNormalized = ep.psiFunc.ev(rWall, zWall).min()
            idx = np.argmin(psiWall)
            idx2 = np.argmin(ep.psiFunc.ev(rWall, zWall))
            #print('R = {:f}, Z = {:f}'.format(rWall[idx2], zWall[idx2]))

#            print('New LCFS R: {:f}'.format(rWall[idx]))
#            print('New LCFS Z: {:f}'.format(zWall[idx]))
#            print('psi at RLCFS, ZLCFS: {:f}'.format(psiNew))
#            print('psiN at RLCFS, ZLCFS: {:f}'.format(psiNewNormalized))
#            print('Original psiSep: {:f}'.format(g['psiSep']))
#            print('Original psiAxis: {:f}'.format(g['psiAxis']))
#            print('Original Nlcfs: {:f}'.format(g['Nlcfs']))
#            log.info('New LCFS R: {:f}'.format(rWall[idx]))
#            log.info('New LCFS Z: {:f}'.format(zWall[idx]))
#            log.info('psi at RLCFS, ZLCFS: {:f}'.format(psiNew))
#            log.info('psiN at RLCFS, ZLCFS: {:f}'.format(psiNewNormalized))
#            log.info('Original psiSep: {:f}'.format(g['psiSep']))
#            log.info('Original psiAxis: {:f}'.format(g['psiAxis']))
#            log.info('Original Nlcfs: {:f}'.format(g['Nlcfs']))

        #this overrides any psiNew we got using R, Z
        if psiSep is not None:
            print("psi provided...overwriting R,Z")
            psiNewNormalized = (psiNew - g['psiAxis']) / (g['psiSep'] - g['psiAxis'])
            psiNew = psiSep
        #find contour
#        import matplotlib.pyplot as plt
#        CS = plt.contourf(g['R'],g['Z'],g['psiRZ'])
#        lcfsCS = plt.contour(CS, levels = [psiNew])
#        rlcfs = lcfsCS.allsegs[0][0][:,0]
#        zlcfs = lcfsCS.allsegs[0][0][:,1]
        surface = ep.getBs_FluxSur(psiNewNormalized)
        rlcfs = surface['Rs']
        zlcfs = surface['Zs']
        print("minimum LCFS R: {:f}".format(min(rlcfs)))
        print("psi at this R: {:f}".format(psiNewNormalized))
        print("Nlcfs")
        print(len(rlcfs))
        g['Nlcfs'] = len(rlcfs)
        g['lcfs'] = np.column_stack((rlcfs,zlcfs))
        g['psiSep'] = psiNew
        print('Renormalized LCFS to new (Rlcfs, Zlcfs)')
        log.info('Renormalized LCFS to new (Rlcfs, Zlcfs)')
        return g


    def check4repeatedEQ(self, ep, EPs):
        """
        checks if an ep object is repeated in a list of EPs
        """
        repeat = None
        for i,e in enumerate(EPs):
            test1 = ep.g['Ip'] == e.g['Ip']
            test2 = ep.g['Bt0'] == e.g['Bt0']
            test3 = np.all(ep.g['psiRZ'] == e.g['psiRZ'])
            test4 = np.all(ep.g['Fpol'] == e.g['Fpol'])
            if test1 and test2 and test3 and test4:
                repeat = i #found it
                print("MHD EQ repeated.  Using EQ from tIdx {:d}".format(i))
                break

        return repeat



    def copyGfile2tree(self,gFileName,shot,time,clobberflag=True):
        """
        Copies gfile to HEAT tree
        gFileName is name of gFile that is already located in self.tmpDir
        """
        oldgfile = self.tmpDir + gFileName
#        #try to make EP object if naming follows d3d gFile naming convention
#        try:
#            ep = EP.equilParams(oldgfile)
#            shot = ep.g['shot']
#            time = ep.g['time']
#
#        #if gfile doesn't follow naming convention define manually
#        except:
#            print("Couldn't open gFile with equilParams_class")
#            log.info("Couldn't open gFile with equilParams_class")
#            if self.shot is None:
#                shot = 1
#            else:
#                shot = self.shot
#            time = idx

        name = 'g{:06d}.{:05d}'.format(shot,time)
        #make tree for this shot
        tools.makeDir(self.shotPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        #make tree for this timestep
        if self.shotPath[-1] != '/': self.shotPath += '/'
        timeDir = self.shotPath + '{:06d}/'.format(time)
        newgfile = timeDir + name

        #clobber and make time directory
        tools.makeDir(timeDir, clobberFlag=clobberflag, mode=self.chmod, UID=self.UID, GID=self.GID)
        shutil.copyfile(oldgfile, newgfile)
        return time


    def writeGfileData(self,gFileList, gFileData):
        """
        writes data passed in string object (from GUI) to files in
        self.tmpDir directory for use later on in HEAT

        the self.tmpDir directory is accessible to the GUI users for uploading
        and downloading

        this function is called from GUI because objects are json / base64
        """
        import base64
        for i,gfile in enumerate(gFileList):
            data = gFileData[i].encode("utf8").split(b";base64,")[1]
            print("Writing gfile: "+gfile)
            log.info("Writing gfile: "+gfile)
            path = self.tmpDir + gfile
            with open(path, 'wb') as f:
                f.write(base64.decodebytes(data))

        return

    def writeGfile(self, file, shot=None, time=None, ep=None):
        """
        writes a new gfile.  for use with the cleaner script.  user must supply
        file: name of new gfile
        shot: new shot number
        time: new shot timestep [ms]
        """
        if ep==None:
            print("Warning no gFile provided.  Writing from gFile in memory.")
            log.info("Warning no gFile provided.  Writing from gFile in memory.")
            g = self.ep.g
        else:
            g = ep.g

        if shot==None:
            shot=1
        if time==None:
            time=1

        KVTOR = 0
        RVTOR = 1.7
        NMASS = 0
        RHOVN = np.zeros((g['NR']))

        print('Writing to path: ' +file)
        with open(file, 'w') as f:
            f.write('  EFIT    xx/xx/xxxx    #' + str(shot) + '  ' + str(time) + 'ms        ')
            f.write('   3 ' + str(g['NR']) + ' ' + str(g['NZ']) + '\n')
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['Xdim'], g['Zdim'], g['R0'], g['R1'], g['Zmid']))
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['RmAxis'], g['ZmAxis'], g['psiAxis'], g['psiSep'], g['Bt0']))
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['Ip'], 0, 0, 0, 0))
            f.write('% .9E% .9E% .9E% .9E% .9E\n'%(0,0,0,0,0))
            self._write_array(g['Fpol'], f)
            self._write_array(g['Pres'], f)
            self._write_array(g['FFprime'], f)
            self._write_array(g['Pprime'], f)
            self._write_array(g['psiRZ'].flatten(), f)
            self._write_array(g['qpsi'], f)
            f.write(str(g['Nlcfs']) + ' ' + str(g['Nwall']) + '\n')
            self._write_array(g['lcfs'].flatten(), f)
            self._write_array(g['wall'].flatten(), f)
            f.write(str(KVTOR) + ' ' + format(RVTOR, ' .9E') + ' ' + str(NMASS) + '\n')
            self._write_array(RHOVN, f)

        print('Wrote new gfile')


    def gFileInterpolate(self, newTime):
        """
        interpolates gfiles in time at newTime
        """
        #Set up all arrays
        ts = []
        RmAxisAll = []
        ZmAxisAll = []
        psiRZAll = []
        psiAxisAll = []
        psiSepAll = []
        Bt0All = []
        IpAll = []
        FpolAll = []
        PresAll = []
        FFprimeAll = []
        PprimeAll = []
        qpsiAll = []

        EPs = self.ep
        for ep in EPs:
            ts.append(ep.g['time'])
            RmAxisAll.append(ep.g['RmAxis'])
            ZmAxisAll.append(ep.g['ZmAxis'])
            psiRZAll.append(ep.g['psiRZ'])
            psiAxisAll.append(ep.g['psiAxis'])
            psiSepAll.append(ep.g['psiSep'])
            Bt0All.append(ep.g['Bt0'])
            IpAll.append(ep.g['Ip'])
            FpolAll.append(ep.g['Fpol'])
            PresAll.append(ep.g['Pres'])
            FFprimeAll.append(ep.g['FFprime'])
            PprimeAll.append(ep.g['Pprime'])
            qpsiAll.append(ep.g['qpsi'])


        R = EPs[0].g['R']
        Z = EPs[0].g['Z']
        ts = np.array(ts)
        RmAxisAll = np.array(RmAxisAll)
        ZmAxisAll = np.array(ZmAxisAll)
        psiRZAll = np.dstack(psiRZAll) #2D
        psiAxisAll = np.array(psiAxisAll)
        psiSepAll = np.array(psiSepAll)
        Bt0All = np.array(Bt0All)
        IpAll = np.array(IpAll)
        FpolAll = np.array(FpolAll).T
        PresAll = np.array(PresAll).T
        FFprimeAll = np.array(FFprimeAll).T
        PprimeAll = np.array(PprimeAll).T
        qpsiAll = np.array(qpsiAll).T
#       FpolAll = np.dstack(FpolAll)
#       FFprimeAll = np.dstack(FFprimeAll)
#       PprimeAll = np.dstack(PprimeAll)
#       qpsiAll = np.dstack(qpsiAll)

        #Set up interpolators
        RmAxisInterp = interp1d(ts,RmAxisAll)
        ZmAxisInterp = interp1d(ts,ZmAxisAll)
        psiRZInterp = RegularGridInterpolator((R, Z, ts), psiRZAll)
        psiAxisInterp = interp1d(ts, psiAxisAll)
        psiSepInterp = interp1d(ts, psiSepAll)
        Bt0Interp = interp1d(ts, Bt0All)
        IpInterp = interp1d(ts, IpAll)
        psiN = np.linspace(0,1,len(FpolAll[:,0]))
        FpolInterp = RegularGridInterpolator((psiN,ts), FpolAll)
        PresInterp = RegularGridInterpolator((psiN,ts), PresAll)
        FFprimeInterp = RegularGridInterpolator((psiN,ts), FFprimeAll)
        PprimeInterp = RegularGridInterpolator((psiN,ts), PprimeAll)
        qpsiInterp = RegularGridInterpolator((psiN,ts), qpsiAll)

        #Interpolate each parameter for the new timestep
        r,z = np.meshgrid(R,Z)
        RmAxis = RmAxisInterp(newTime)
        ZmAxis = ZmAxisInterp(newTime)
        psiRZ = psiRZInterp((r,z,newTime)).T
        psiAxis = psiAxisInterp(newTime)
        psiSep = psiSepInterp(newTime)
        Bt0 = Bt0Interp(newTime)
        Ip = IpInterp(newTime)
        Fpol = FpolInterp((psiN,newTime))
        Pres = PresInterp((psiN,newTime))
        FFprime = FFprimeInterp((psiN,newTime))
        Pprime = PprimeInterp((psiN,newTime))
        qpsi = qpsiInterp((psiN,newTime))

        #get new LCFS
        surface = ep.getBs_FluxSur(1.0)
        rlcfs = surface['Rs']
        zlcfs = surface['Zs']

        #with linear interpolation, infrequently RegularGridInterpolator cannot
        #resolve the points correctly and returns NANs.  This messes up future gFile
        #reading algorithms (although it shouldnt! come on!).  See RegularGridInterpolator
        #python webpage for more info
        #first do R
        mask = np.ones(len(rlcfs), dtype=bool)
        mask[np.argwhere(np.isnan(rlcfs))] = False
        rlcfs = rlcfs[mask]
        zlcfs = zlcfs[mask]
        #then do Z
        mask = np.ones(len(zlcfs), dtype=bool)
        mask[np.argwhere(np.isnan(zlcfs))] = False
        rlcfs = rlcfs[mask]
        zlcfs = zlcfs[mask]

        #make new dictionary with all this stuff
        newEP = lambda: None #empty object
        newEP.g = {
                    'RmAxis':RmAxis,
                    'ZmAxis':ZmAxis,
                    'psiRZ':psiRZ,
                    'psiAxis':psiAxis,
                    'psiSep':psiSep,
                    'Bt0':Bt0,
                    'Ip':Ip,
                    'Fpol':Fpol,
                    'Pres':Pres,
                    'FFprime':FFprime,
                    'Pprime':Pprime,
                    'qpsi':qpsi,
                    'NR':EPs[0].g['NR'],
                    'NZ':EPs[0].g['NZ'],
                    'Xdim':EPs[0].g['Xdim'],
                    'Zdim':EPs[0].g['Zdim'],
                    'R0':EPs[0].g['R0'],
                    'R1':EPs[0].g['R1'],
                    'Zmid':EPs[0].g['Zmid'],
                    'wall':EPs[0].g['wall'],
                    'Nwall':EPs[0].g['Nwall'],
                    'Nlcfs':len(rlcfs),
                    'lcfs':np.column_stack((rlcfs,zlcfs)),

                    }


        #return ep that can be written to file (note its not a real EP as defined by equilParams class)
        return newEP


    def GEQDSKFromMDS(machine, shot, tree='efit01', tmin=None, tmax=None,
                          rootDir=None, clobber=True, chmod=0o774, UID=-1, GID=-1 ):
        """
        Function for grabbing multiple gfiles associated with a single shot
        machine: string. either 'nstx' or 'd3d' or 'st40'
        you will need to add server info for others
        shot: integer.  Shot number of interest
        tmin: integer. initial time step in [ms]
        tmax: integer. final time step in [ms]
        note: if tmin and tmax are 'None', then we pull all timsteps
        if only tmin = 'None' then tmin = minimum timestep
        if only tmax = 'None' then tmax = maximum timestep
        """
        #Tune server address depending upon machine
        if machine == 'nstx': Server='skylark.pppl.gov'
        elif machine == 'd3d': Server='atlas.gat.com'
        elif machine == 'st40':
            Server='skylark2.pppl.gov'
            tree='st40'

        #Connect to server and get all timesteps.
        import MDSplus
        MDS = MDSplus.Connection(Server)
        MDS.openTree(tree, shot)
        base = 'RESULTS:GEQDSK:'
        # get time slice
        signal = 'GTIME'
        timesteps = MDS.get(base + signal).data()*1000
        MDS.closeTree(tree, shot)
        #Handle the multiple potential inputs discussed above
        if tmin==None and tmax==None:
            tmin = 0
            tmax = timesteps[-1]
        elif tmin==None and tmax!=None: tmin = 0
        elif tmin!=None and tmax==None: tmax = timesteps[-1]

        # find times between tmin and tmax
        use = np.where(np.logical_and(timesteps >=tmin, timesteps <=tmax))
        ts = timesteps[use].astype(int)

        #Directory where we will save these gfiles if one is not declared
        if rootDir==None: dirName = machine + '_{:06d}'.format(shot)
        else: dirName = rootDir

        tools.makeDir(dirName, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)

        #Pull gfiles and write into dirName/subdirectory where subdirectory
        #is exclusively for that timestep
        print('Pulling all gfiles between {:5d} and {:5d} [ms]'.format(int(tmin),int(tmax)))
        for t in ts:
            t_path = dirName + '/{:06d}'.format(t)
            tools.makeDir(t_path, clobberFlag=True, mode=self.chmod, UID=self.UID, GID=self.GID)
            self.MDSplusPull(shot,t,tree=tree,Server=Server,gpath=t_path)

        return ts


    def MDSplusPull(self, shot, time, tree='EFIT01', exact=False, Server='skylark.pppl.gov',
              gpath='.'):

        print('Reading shot =', shot, 'and time =', time, 'from MDS+ tree:', tree)
        # in case those are passed in as strings
        shot = int(shot)
        time = int(time)
        # Connect to server, open tree and go to g-file
        MDS = MDSplus.Connection(Server)
        MDS.openTree(tree, shot)
        base = 'RESULTS:GEQDSK:'
        # get time slice
        signal = 'GTIME'
        k = np.argmin(np.abs(MDS.get(base + signal).data()*1000 - time))
        time0 = int(MDS.get(base + signal).data()[k]*1000)
        if (time != time0):
            if exact:
                raise RuntimeError(tree + ' does not exactly contain time ' + str(time) + '  ->  Abort')
            else:
                print('Warning: ' + tree + ' does not exactly contain time ' + str(time) + ' the closest time is ' + str(time0))
                print('Fetching time slice ' + str(time0))
                #time = time0
        # store data in dictionary
        g = {'shot':shot, 'time':time}
        # get header line
        header = str(MDS.get(base + 'CASE').data()[k], 'utf-8')

        # get all signals, use same names as in read_g_file
        translate = {'MW': 'NR', 'MH': 'NZ', 'XDIM': 'Xdim', 'ZDIM': 'Zdim', 'RZERO': 'R0',
                    'RMAXIS': 'RmAxis', 'ZMAXIS': 'ZmAxis', 'SSIMAG': 'psiAxis',
                    'SSIBRY': 'psiSep', 'BCENTR': 'Bt0', 'CPASMA': 'Ip', 'FPOL': 'Fpol',
                    'PRES': 'Pres', 'FFPRIM': 'FFprime', 'PPRIME': 'Pprime', 'PSIRZ': 'psiRZ',
                    'QPSI': 'qpsi', 'NBDRY': 'Nlcfs', 'LIMITR': 'Nwall'}

        for signal in translate:
            try:
                g[translate[signal]] = MDS.get(base + signal).data()[k]
            except:
                raise ValueError(signal +' retrieval failed.')
        g['R1'] = MDS.get(base + 'RGRID1').data()[0]
        g['Zmid'] = 0.0
#       RLIM = MDS.get(base + 'RLIM').data()[:, 0]
#       ZLIM = MDS.get(base + 'ZLIM').data()[:, 1]
        RLIM = MDS.get(base + 'RLIM').data()[k]
        ZLIM = MDS.get(base + 'ZLIM').data()[k]
        g['wall'] = np.vstack((RLIM, ZLIM)).T
        RBBBS = MDS.get(base + 'RBDRY').data()[k]# [:g['Nlcfs']]
        ZBBBS = MDS.get(base + 'ZBDRY').data()[k]# [:g['Nlcfs']]
        g['lcfs'] = np.vstack((RBBBS, ZBBBS)).T
        KVTOR = 0
        RVTOR = 0#1.7
        NMASS = 0
        RHOVN = MDS.get(base + 'RHOVN').data()[k]

        # convert floats to integers
        for item in ['NR', 'NZ', 'Nlcfs', 'Nwall']:
            g[item] = int(g[item])
        # convert single (float32) to double (float64) and round
        for item in ['Xdim', 'Zdim', 'R0', 'R1', 'RmAxis', 'ZmAxis', 'psiAxis', 'psiSep',
                    'Bt0', 'Ip']:
            g[item] = np.round(np.float64(g[item]), 7)
        # convert single arrays (float32) to double arrays (float64)
        for item in ['Fpol', 'Pres', 'FFprime', 'Pprime', 'psiRZ', 'qpsi', 'lcfs', 'wall']:
            g[item] = np.array(g[item], dtype=np.float64)
        # write g-file to disk
        if not (gpath[-1] == '/'):
            gpath += '/'

        # Writing to match J. Menard's idl script
        # in /u/jmenard/idl/efit/efit_routines_jem.pro
        # Scale FFprime, Pprime, psiRZ, qpsi, Ip, psiSep, psiAxis,
        # per the AUTO keyword in J. Menard's function, 'rescale_gstructures'
        #       -From the Menard script line 395:
        # "This function re-scales the poloidal flux to be consistent with the
        # specified sign of Ip, where Ip is the toroidal component of the plasma
        # current in cylindrical coordinates, i.e. + --> C.C.W. as viewed from
        # above.  Note that in the re-definition of q, the theta flux coordinate
        # is C.C.W. in the poloidal plane, so that if Ip > 0 and Bt > 0, q will
        # be negative, since then B.grad(theta) < 0."

#       dsign = np.sign(g['Ip'])
#       gsign = np.sign( (g['psiAxis'] - g['psiSep']) )
#       qsign = np.sign(g['Fpol'][-1]) #F_edge sign
#       g['FFprime'] *= dsign*gsign
#       g['Pprime'] *= dsign*gsign
#       g['psiRZ'] *= dsign*gsign
#       g['qpsi'] *= -qsign*dsign
#       g['Ip'] *= dsign
#       g['psiSep'] *= dsign*gsign
#       g['psiAxis'] *= dsign*gsign

#       print('dsign = {:f}'.format(dsign))
#       print('gsign = {:f}'.format(gsign))
#       print('qsign = {:f}'.format(qsign))
#       g['FFprime'] *= dsign*qsign*gsign
#       g['Pprime'] *= dsign*qsign*gsign
#       g['psiRZ'] *= dsign*qsign*gsign
#       g['qpsi'] *= -qsign*dsign
#       g['psiSep'] *= dsign*qsign*gsign
#       g['psiAxis'] *= dsign*qsign*gsign
#       g['Fpol'] *= dsign*qsign

        # Now, write to file using same style as J. Menard script (listed above)
        # Using function in WRITE_GFILE for reference
        file = gpath + 'g' + format(shot, '06d') + '.' + format(time,'05d')
        with open(file, 'w') as f:
            if ('EFITD' in header or 'LRD' in header):
                f.write(header)
            else:
                f.write('  EFITD    xx/xx/xxxx    #' + str(shot) + '  ' + str(time) + 'ms        ')
                f.write('   3 ' + str(g['NR']) + ' ' + str(g['NZ']) + '\n')
                f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['Xdim'], g['Zdim'], g['R0'], g['R1'], g['Zmid']))
                f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['RmAxis'], g['ZmAxis'], g['psiAxis'], g['psiSep'], g['Bt0']))
                f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['Ip'], g['psiAxis'], 0, g['RmAxis'], 0))
                f.write('% .9E% .9E% .9E% .9E% .9E\n'%(g['ZmAxis'],0,g['psiSep'],0,0))
                write_array(g['Fpol'], f)
                write_array(g['Pres'], f)
                write_array(g['FFprime'], f)
                write_array(g['Pprime'], f)
                write_array(g['psiRZ'].flatten(), f)
                write_array(g['qpsi'], f)
                f.write(str(g['Nlcfs']) + ' ' + str(g['Nwall']) + '\n')
                write_array(g['lcfs'].flatten(), f)
                write_array(g['wall'].flatten(), f)
                f.write(str(KVTOR) + ' ' + format(RVTOR, ' .9E') + ' ' + str(NMASS) + '\n')
                write_array(RHOVN, f)
        os.chown(file, self.UID, self.GID)
        os.chmod(file, self.chmod)
        return time

    #=====================================================================
    #                       private functions
    #=====================================================================
    # --- _write_array -----------------------
    # write numpy array in format used in g-file:
    # 5 columns, 9 digit float with exponents and no spaces in front of negative numbers
    def _write_array(self, x, f):
        N = len(x)
        rows = int(N/5)  # integer division
        rest = N - 5*rows
        for i in range(rows):
            for j in range(5):
                    f.write('% .9E' % (x[i*5 + j]))
            f.write('\n')
        if(rest > 0):
            for j in range(rest):
                f.write('% .9E' % (x[rows*5 + j]))
            f.write('\n')

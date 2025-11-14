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
import matplotlib.pyplot as plt
import netCDF4

import logging
log = logging.getLogger(__name__)

def setupForTerminalUse(gFile=None, shot=None, time=0.0):
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
        #single geqdsk
        if type(gFile)==str:
            EQmode = MHD.determineEQFiletype(gFile)
            MHD.ep = EP.equilParams(gFile, EQmode=EQmode, time=time)
        #multiple geqdsks, filenames in a list
        elif type(gFile)==list:
            MHD.ep = []
            for f in gFile:
                MHD.ep.append(EP.equilParams(f))
        else:
            print("Unknown type for gFile argument.  Provide string path or list of paths")


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
        
        MHD EQ Variables:
        -----------------

        :shot: integer pulse number
        :tmin: minimum timestep of any MHD equilibrium in simulation [s]
        :tmax: maximum timestep of any MHD equilibrium in simulation [s]
        :traceLength: number of steps to trace along magnetic field lines looking for
          intersections
        :dpinit: toroidal length of each trace step up magnetic field line [degrees]
        :psiMult: multiplier to apply to all psi quantities (psiSep, psiAxis, psiRZ) for normalization (ie the dreaded 2pi) in EQ
        :BtMult: multiplier to apply to toroidal field quantities (Bt0, Fpol) in EQ
        :IpMult: multiplier to apply to plasma current, Ip in EQ
        :BtraceFile: path to file to be used for field line tracing

        """
        self.allowed_vars = [
                            'shot',
                            'tmin',
                            'tmax',
                            'traceLength',
                            'dpinit',
                            'psiMult',
                            'BtMult',
                            'IpMult',
                            'BtraceFile',
                            ]

        return


    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        integers = [
                    'shot',
                    'traceLength',
                    ]
        floats = [
                    'tmin',
                    'tmax',
                    'dpinit',
                    'psiMult',
                    'BtMult',
                    'IpMult',
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


    def determineEQFiletype(self,file):
        '''
        tests to see if an EQ file is a netcdf (.nc) or JSON (.json)

        if neither, we assume it is a GEQDSK
        '''
        extension = os.path.splitext(file)[1]
        if extension == '.nc':
            EQmode = 'netcdf'
        elif extension == '.json':
            EQmode = 'json'
        else:
            EQmode = 'geqdsk'
        return EQmode


    def getGEQDSKtimesteps(self, eqList):
        """
        gets timesteps for gfiles in tmpDir

        checks to see if they are in the correct format, otherwise just arbitrarily assigns time

        called from the GUI
        """
        #check if these GEQDSKs are named according to the HEAT convention g<shot>_<timestep>
        #where shot is an int and timestep can be any float (ie with or without radix)
        #test1 = np.all(np.array([len(x.split("_")) for x in eqList]) > 1)
        #test2 = np.all(np.array([type(g.split('_')[-1])==float for g in eqList]) > 1)
        useD3DtimeFmt = True
        useHeatTimeFmt = True
        for g in eqList:
            if ('.' in g) & ('_' not in g):
                time = g.split('.')[-1]
                try: time = float(time)
                except: useD3DtimeFmt = False
            else: useD3DtimeFmt = False
            if ('_' in g): 
                time = g.split('_')[-1]
                try: time = float(time)
                except: useHeatTimeFmt = False
            else: useHeatTimeFmt = False

        print("HEAT time format: " + str(useHeatTimeFmt))
        print("D3D time format: " + str(useD3DtimeFmt))
        log.info("HEAT time format: " + str(useHeatTimeFmt))
        log.info("D3D time format: " + str(useD3DtimeFmt))
 

        ts = []

        for i,g in enumerate(eqList):
            #GEQDSKs are named using timesteps
            if useHeatTimeFmt:
                ts.append(float(g.split('_')[-1]))
            #GEQDSKs are named by D3D naming convention
            elif useD3DtimeFmt:
                ts.append(float(g.split('.')[-1]))
            #GEQDSKs do not follow HEAT or D3D naming convention
            else:
                ts.append(float(i))
        return np.array(ts)

    def getGEQDSK(self, ts, eqList):
        """
        copies EQ into the HEAT output tree

        if netcdf EQ was provided, converts to GEQDSK (required for MAFOT)
        

        ts is list of timesteps
        eqList is list of names of eq files

        ts and eqList are indexed to match each other if GEQDSK
           if not GEQDSK, then eqList is single file with all eq inside

        geqdsk text file naming format for HEAT is:
        g = 'g'+self.shotFmt.format(self.shot) +'_'+ self.tsFmt.format(t)

        if netcdf, then should follow the IMAS standards

        psiMult is multiplier for psi quantities (ie psiRZ, psiSep, psiAxis)
        BtMult is multiplier for Bt quantities (Bt, Fpol)
        IpMult is multipler for Ip 
        """
        #setup multipliers if necessary
        if self.psiMult == None:
            self.psiMult = 1.0
        if self.BtMult == None:
            self.BtMult = 1.0
        if self.IpMult == None:
            self.IpMult = 1.0

        self.timesteps = ts
        self.gFiles = []
        self.eqFiles = []
        for i,t in enumerate(ts):
            #if all the EQ are in a single file, account for it here (GUI mode)
            #otherwise, use the EQ assigned to each timestep (TUI mode)
            if len(ts) > len(eqList):
                eq = eqList[0]
            else:
                eq = eqList[i]
            #in GUI tmpdir is upload dir.  in TUI, its the machDir
            oldeqfile = self.tmpDir + eq 
            timeDir = self.shotPath + self.tsFmt.format(t) +'/'
            if self.shotPath[-1] != '/': self.shotPath += '/'
            self.gFiles.append('g'+self.shotFmt.format(self.shot) + '_'+ self.tsFmt.format(t))
            newgfile = timeDir + self.gFiles[-1]
            
            #check if this is a netcdf and convert to MAFOT compatible file type (GEQDSK)
            EQmode = self.determineEQFiletype(oldeqfile)
            print("Equilibrium file for timestep "+str(t)+" in "+EQmode+" format")
            if EQmode != 'geqdsk':
                #convert to GEQDSK (for MAFOT)
                ep = EP.equilParams(oldeqfile, EQmode=EQmode, time=t, 
                                    psiMult=self.psiMult, BtMult=self.BtMult, IpMult=self.IpMult)
                self.writeGfile(newgfile, shot=self.shot, time=t, ep=ep)
            else:
                #file is already in MAFOT compatible format
                shutil.copyfile(oldeqfile, newgfile)

            #if not a geqdsk, still save the eq files into the HEAT tree
            if EQmode != 'geqdsk':
                #if there is one eq for all ts, save that to the shotDir
                if len(eqList) == 1:
                    #save the JSON/netcdf into the HEAT tree
                    self.singleEQfile = self.shotPath + eq
                    shutil.copyfile(oldeqfile, self.singleEQfile)
                #if there is an eq for each ts, save each to the corresponding timeDir
                else:
                    self.singleEQfile = None
                    neweqfile = timeDir + eq
                    self.eqFiles.append(neweqfile)
                    shutil.copyfile(oldeqfile, neweqfile)
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
            timeDir = self.shotPath + self.tsFmt.format(t) +'/'
            gfile = timeDir + 'g'+self.shotFmt.format(shot)+'_'+self.tsFmt.format(t)
        self.ep = EP.equilParams(gfile)#, shot, t)#, gtype='heat')
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
            timeDir = self.shotPath + self.tsFmt.format(t) +'/'
            #netcdf or json
            if self.EQmode != 'geqdsk':
                if self.singleEQfile == None:
                    eqfile = self.eqFiles[idx]
                else:
                    eqfile = self.singleEQfile
                self.ep[idx] = EP.equilParams(eqfile, EQmode=self.EQmode, time=t, 
                                              psiMult=self.psiMult, BtMult=self.BtMult, IpMult=self.IpMult)
            #geqdsk
            else:
                gfile = timeDir+self.gFiles[idx]
                self.ep[idx] = EP.equilParams(gfile, psiMult=self.psiMult, BtMult=self.BtMult, IpMult=self.IpMult)
        return




    def Bfield_pointcloud(self, ep, R, Z, phi, powerDir=None, normal=False,
                          helicityCheck=True):
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

        if helicityCheck is True, we update helicity based upon signs of variables
        in GEQDSK

        """

        Bt = ep.BtFunc.ev(R,Z)
        BR = ep.BRFunc.ev(R,Z)
        BZ = ep.BZFunc.ev(R,Z)
        BtMult = 1.0
        BpMult = 1.0

        if helicityCheck == True:
            #check for signs
            FSign = np.sign(ep.g['Fpol'][-1])
            Bt0Sign = np.sign(ep.g['Bt0'])
            IpSign = np.sign(ep.g['Ip'])
            psiSign = np.sign(ep.psiFunc.ev(R,Z))

            #correct field directions to match psi, Fpol, Ip, Bt0
            #the ep object (equilParams_class.py) uses Fpol to define Bt everywhere,
            #so we check to make sure the sign matches Bt0
            if FSign != Bt0Sign:
                BtMult = -1.0
            else:
                BtMult = 1.0
            #the ep object defines B by derivatives of psi, but doesnt take into
            #account whether or not psi is increasing / decreasing from axis to sep
            if ep.g['psiSep'] > ep.g['psiAxis']: #dPsi/dR > 0 at OMP
                BpMult = 1.0
            else:
                BpMult = -1.0
            #the ep object does not check for the sign of Ip when defining Bp,
            #so we check Bp at the omp to see if it follows the sign of Ip
            if np.sign(ep.g['Ip']) > 0:
                BpMult *= 1.0
            else:
                BpMult *= -1.0

            print("\n#====  Bfield helicity check ====")
            print("Fpol sign: {:f}".format(FSign))
            print("Bt0 sign: {:f}".format(Bt0Sign))
            print("Ip sign: {:f}".format(IpSign))
            print("psiRZ sign [0]: {:f}".format(psiSign[0]))
            print("BtMult: {:f}".format(BtMult))
            print("BpMult: {:f}".format(BpMult))
            log.info("\n#====  Bfield helicity check ====")
            log.info("Fpol sign: {:f}".format(FSign))
            log.info("Bt0 sign: {:f}".format(Bt0Sign))
            log.info("Ip sign: {:f}".format(IpSign))
            log.info("psiRZ sign [0]: {:f}".format(psiSign[0]))
            log.info("BtMult: {:f}".format(BtMult))
            log.info("BpMult: {:f}".format(BpMult))

        Bt *= BtMult
        BR *= BpMult
        BZ *= BpMult

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


    def setupMAFOTdirectory(self, controlfilePath:str, obj:object):
        """
        sets up necessary paths for MAFOT and saves them as object variables
        """
        if controlfilePath[-1]!='/': controlfilePath+='/'
        
        obj.controlfilePath = controlfilePath
        obj.controlfile = '_lamCTL.dat'
        obj.controlfileStruct = '_struct_CTL.dat'
        obj.gridfile = obj.controlfilePath + 'grid.dat'
        obj.gridfileStruct = obj.controlfilePath + 'struct_grid.dat'
        obj.outputFile = obj.controlfilePath + 'lam.dat'
        obj.structOutfile = obj.controlfilePath + 'struct.dat'
        return obj

    def writeControlFile(self, name, t, mapDirection, mode='laminar'):
        """
        Create a control file for MAFOT
        """
        if len(name.split('/')) > 1:
            save_location = name
        else:
            save_location = self.shotPath + self.tsFmt.format(t)+'/' + name


        with open(save_location, 'w') as f:

            f.write('# Parameterfile for HEAT Programs\n')
            f.write('# Shot: '+self.shotFmt.format(self.shot)+'\tTime: '+self.tsFmt.format(t)+'s\n')
            f.write('# Path: ' + self.shotPath  +self.tsFmt.format(t) + '/' + 'g'+self.shotFmt.format(self.shot) + '_'+ self.tsFmt.format(t) + '\n')
            #f.write('# Path: ' + self.shotPath  +self.tsFmt.format(t) + '/' + '\n')
            f.write('Nphi=\t{:d}\n'.format(self.Nphi))    # This must be entry index 0

            #itt means different things depending on if we are tracing field line
            #or running full MAFOT laminar
            if mode=='laminar':
                f.write('itt=\t{:f}\n'.format(self.ittLaminar))
            elif mode=='gyro':
                f.write('itt=\t{:f}\n'.format(self.ittGyro))
            else:
                f.write('itt=\t{:f}\n'.format(self.ittStruct))
                
            f.write('Rmin=\t{:2f}\n'.format(self.Rmin))    # This must be entry index 2
            f.write('Rmax=\t{:2f}\n'.format(self.Rmax))
            f.write('Zmin=\t{:2f}\n'.format(self.Zmin))
            f.write('Zmax=\t{:2f}\n'.format(self.Zmax))
            f.write('Nswall=\t{:d}\n'.format(self.Nswall))
            f.write('phistart(deg)=\t{:2f}\n'.format(self.phistart))
            f.write('MapDirection=\t{:f}\n'.format(mapDirection))
            f.write('PlasmaResponse(0=no,>1=yes)=\t{:d}\n'.format(self.PlasmaResponse))
            f.write('Field(-3=VMEC,-2=SIESTA,-1=gfile,M3DC1:0=Eq,1=I-coil,2=both)=\t{:d}\n'.format(self.Field))
            f.write('target(0=useSwall)=\t{:d}\n'.format(self.target))
            f.write('createPoints(2=target)=\t{:d}\n'.format(self.createPoints))    # This must be entry index 12
            #f.write('useFcoil(0=no,1=yes)=\t{:d}\n'.format(self.useFcoil))
            #f.write('useCcoil(0=no,1=yes)=\t{:d}\n'.format(self.useCcoil))
            #f.write('useIcoil(0=no,1=yes)=\t{:d}\n'.format(self.useIcoil))
            f.write('unused=\t0\n')
            f.write('unused=\t0\n')
            f.write('unused=\t0\n')
            f.write('ParticleDirection(1=co-pass,-1=ctr-pass,0=field-lines)=\t{:d}\n'.format(self.ParticleDirection))   # This must be entry index 16
            f.write('PartileCharge(-1=electrons,>=1=ions)=\t{:d}\n'.format(self.ParticleCharge))
            f.write('Ekin[keV]=\t{:2f}\n'.format(self.Ekin))
            f.write('lambda=\t{:2f}\n'.format(self.Lambda))
            f.write('Mass=\t{:2f}\n'.format(self.Mass))    # This must be entry index 20
            #f.write('useFilament(0=no)=\t{:d}\n'.format(self.useFilament))
            #f.write('useECcoil(0=no,1=yes)=\t{:d}\n'.format(self.useECcoil))
            f.write('unused=\t0\n')
            f.write('unused=\t0\n')
#            f.write('unused=\t0\n')
            f.write('dpinit=\t{:f}\n'.format(self.dpinit)) # This must be entry index 23
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
        PFC.psimin = PFC.ep.psiFunc.ev(R,Z) #psi_N
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


    def distance(self, traceData:np.ndarray):
        """
        Calculate distance along
        """
        distance = np.cumsum(np.sqrt(np.sum(np.diff(traceData,axis=0)**2,axis=1)))
        distance = np.insert(distance, 0, 0)
        return distance

    def getFieldpath(self, dphi, dpinit, gridfile, controlfilePath,
                    controlfile, paraview_mask=False, tag=None, bbox=True, distance_mask=True):
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
        #step size for output
        args.append('-d')
        args.append(str(dphi))
        #integrator step size
        args.append('-i')
        args.append(str(dpinit))
        #the points that we launch traces from
        args.append('-P')
        args.append(gridfile)
        #now use simple boundary instead of g-file wall as limiter
	    #This must be used with MAFOT version 5.7 or newer for the Shadow Mask calculation
        if bbox == True:
            args.append('-b')
	    #args 6 is the MAFOT control file
        args.append(controlfile)
        #args 6 is the tag, if not None
        if tag is not None: args.append(tag)
        #Copy the current environment (important when in appImage mode)
        current_env = os.environ.copy()
        #run MAFOT structure for points in gridfile
        print(args)
        from subprocess import run
        run(args, env=current_env, cwd=controlfilePath)
        try:
            print("Removing MAFOT logs")
            log.info("Removing MAFOT logs")
            run(['rm', 'log_*'], env=current_env, cwd=controlfilePath)
        except:
            print("Cannot delete MAFOT logs")
            log.info("Cannot delete MAFOT logs")

        if paraview_mask:
            #This file is not in the format we want, so we read it in, and then
            #write it out as a CSV.  Because this is a serial function that is only
            #called for special cases, we arent too concerned about speed or memory
            #so this works instead of making a more efficient pipeline
            if tag is None: readfile = controlfilePath + 'struct.dat'
            else: readfile = controlfilePath + 'struct_' + tag + '.dat'
            structdata = np.genfromtxt(readfile,comments='#')

            xyz = np.zeros((len(structdata),3))
            xyz[:,0] = structdata[:,0]
            xyz[:,1] = structdata[:,1]
            xyz[:,2] = structdata[:,2]
            xyz*=1000 #Scale for ParaVIEW
            xyz = xyz[~(np.abs(xyz[:,:2])<1e-50).any(axis=1)]
            head = 'X[mm],Y[mm],Z[mm],distance[mm]'
            #head = 'X[m],Y[m],Z[m],R[m],phi[rad]'
            if tag == None:
                outfile = controlfilePath+'struct.csv'
                vtkName = 'Field_trace'
            else:
                outfile = controlfilePath+'struct_'+tag+'.csv'
                vtkName = 'Field_trace_'+tag

            d = self.distance(xyz)
            xyzd = np.hstack((xyz, d.reshape(-1,1)))
            np.savetxt(outfile, xyzd, delimiter=',', header=head)
            
            #Now save a vtk file for paraviewweb
            #old method (new method happens in engineClass using ioClass)
            #delete this comment after v5.0
            #tools.createVTKOutput(outfile, 'trace', vtkName)

        return

    def getMultipleFieldPaths(self, dphi, gridfile, controlfilePath,
                                controlfile, bbox=True, paraview_mask=False):

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
        #args 5 use simple boundary instead of g-file wall as limiter
	    #This must be used with MAFOT version 5.7 or newer for the Shadow Mask calculation
        if bbox is True:
            args.append('-b')
	    #args 6 is the MAFOT control file
        args.append(controlfile)
        #Copy the current environment (important when in appImage mode)
        current_env = os.environ.copy()
        #run MAFOT structure for points in gridfile
        from subprocess import run
        run(args, env=current_env, cwd=controlfilePath)
        return


    def readBtraceFile(self):
        """
        .. Reads a Btrace file (file for tracing magnetic field lines from)
        .. first row file should be formatted like this, which defines the columns:
        .. x[mm], y[mm], z[mm], traceDirection, Length[deg], stepSize[deg]
        .. every new row represents a new trace that will be run in HEAT
        
        Columns in the file:
        --------------------

        :x: x coordinate of field line trace start locations in [mm]
        :y: y coordinate of field line trace start locations in [mm]
        :z: z coordinate of field line trace start locations in [mm]
        :traceDirection: can be 1 (forward toroidal angle), -1 (reverse toroidal angle), or 0 (both directions) to trace from point
        :Length[deg]: distance in degrees to trace from the point
        :stepSize[deg]: step size for magnetic field line integration steps in degrees
                
        """

        #read file
        data = pd.read_csv(self.BtraceFile, delimiter=",")
        data = data.astype({"x[mm]": float, "y[mm]": float, "z[mm]": float, "traceDirection": int, "Length[deg]":float, "stepSize[deg]":float})
        return data


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



    def writeEQdata(self,eqList, eqData):
        """
        writes data passed in string object (from GUI) to files in
        self.tmpDir directory for use later on in HEAT

        the self.tmpDir directory is accessible to the GUI users for uploading
        and downloading

        this function is called from GUI because objects are json / base64
        """
        import base64
        for i,gfile in enumerate(eqList):
            data = eqData[i].encode("utf8").split(b";base64,")[1]
            print("Writing gfile: "+gfile)
            log.info("Writing gfile: "+gfile)
            path = self.tmpDir + gfile
            with open(path, 'wb') as f:
                f.write(base64.decodebytes(data))

        return

    def writeNetCDF(self, data_dict, filename):
        """
        Write a dictionary containing numpy arrays to a NetCDF file.

        :param data_dict: Dictionary containing numpy arrays.
        :param filename: Name of the NetCDF file to write to.
        """
        import netCDF4 as nc
        with nc.Dataset(filename, 'w') as ds:
            for key, value in data_dict.items():
                if isinstance(value, np.ndarray):
                    # Create dimensions
                    for i, dim_size in enumerate(value.shape):
                        dim_name = f"{key}_dim_{i}"
                        ds.createDimension(dim_name, dim_size)

                    # Create the variable and assign data
                    var = ds.createVariable(key, value.dtype, tuple(f"{key}_dim_{i}" for i in range(value.ndim)))
                    var[:] = value

                elif isinstance(value, str):
                    # Handle string values
                    str_dim = f"{key}_strlen"
                    ds.createDimension(str_dim, len(value))
                    var = ds.createVariable(key, "S1", (str_dim,))
                    # Assign the string as a sequence of individual characters
                    var[:] = nc.stringtochar(np.array(list(value), 'S1'))

                elif isinstance(value, (int, float)):
                    # Handle scalar numeric values
                    var = ds.createVariable(key, type(value))
                    var[:] = value

        return


    def writeGfile(self, file, shot=None, time=None, ep=None):
        """
        writes a new gfile.  for use with the cleaner script.  user must supply
        file: name of new gfile
        shot: new shot number
        time: new shot timestep [ms]

        Note that this writes some data as 0 (ie rhovn, kvtor, etc.)
        """
        if ep==None:
            print("Warning no gFile provided for write operation.  Writing from gFile in memory.")
            log.info("Warning no gFile provided for write operation.  Writing from gFile in memory.")
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
            if ep.g['psiRZ'].shape[0] == ep.g['Xdim']:
                psiRZAll.append(ep.g['psiRZ'])
            else:
                psiRZAll.append(ep.g['psiRZ'].T)
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
        psiRZ = psiRZInterp((r,z,newTime))
        psiAxis = psiAxisInterp(newTime)
        psiSep = psiSepInterp(newTime)
        Bt0 = Bt0Interp(newTime)
        Ip = IpInterp(newTime)
        Fpol = FpolInterp((psiN,newTime))
        Pres = PresInterp((psiN,newTime))
        FFprime = FFprimeInterp((psiN,newTime))
        Pprime = PprimeInterp((psiN,newTime))
        qpsi = qpsiInterp((psiN,newTime))


        #get new lcfs
        levels = np.linspace(psiSep-0.05,psiSep+0.05,15)
        CS = plt.contourf(R,Z,psiRZ,levels,cmap=plt.cm.cividis)
        lcfsCS = plt.contour(CS, levels = [psiSep])
        #assuming plasma is approx. centered in machine here
        zMin = ZmAxis - 0.25
        zMax = ZmAxis + 0.25
        for i in range(len(lcfsCS.allsegs[0])):
            rlcfs = lcfsCS.allsegs[0][i][:,0]
            zlcfs = lcfsCS.allsegs[0][i][:,1]
            #this prevents us from getting locations not at midplane
            # (ie psiSep locations around coils outside wall)
            idx = np.where(np.logical_and(zlcfs>zMin,zlcfs<zMax))[0]
            if len(idx)>0: #found it
                break

        ##get new LCFS
        ##legacy method that works, but depends on ep (set in loop above)
        #surface = ep.getBs_FluxSur(1.0)
        #rlcfs = surface['Rs']
        #zlcfs = surface['Zs']
        ##with linear interpolation, infrequently RegularGridInterpolator cannot
        ##resolve the points correctly and returns NANs.  This messes up future gFile
        ##reading algorithms (although it shouldnt! come on!).  See RegularGridInterpolator
        ##python webpage for more info
        ##first do R
        #mask = np.ones(len(rlcfs), dtype=bool)
        #mask[np.argwhere(np.isnan(rlcfs))] = False
        #rlcfs = rlcfs[mask]
        #zlcfs = zlcfs[mask]
        ##then do Z
        #mask = np.ones(len(zlcfs), dtype=bool)
        #mask[np.argwhere(np.isnan(zlcfs))] = False
        #rlcfs = rlcfs[mask]
        #zlcfs = zlcfs[mask]
        #lcfsAll.append()

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

    def gFileInterpolateByS(self, newS, transposePsi=True, transposeFpol=True):
        """
        interpolates gfiles as a function of MHD.Spols at newS

        requires MHD object to have list of same length as MHD.ep, with the values
        of Spol for each ep.  this variable is MHD.Spols
        """
        #Set up all arrays
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
            if transposePsi==True:
                psiRZAll.append(ep.g['psiRZ'].T)
            else:
                psiRZAll.append(ep.g['psiRZ'])
                

            RmAxisAll.append(ep.g['RmAxis'])
            ZmAxisAll.append(ep.g['ZmAxis'])
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
        Spols = np.array(self.Spols)
        RmAxisAll = np.array(RmAxisAll)
        ZmAxisAll = np.array(ZmAxisAll)
        psiRZAll = np.dstack(psiRZAll) #2D
        psiAxisAll = np.array(psiAxisAll)
        psiSepAll = np.array(psiSepAll)
        Bt0All = np.array(Bt0All)
        IpAll = np.array(IpAll)
        if transposeFpol==True:
            FpolAll = np.array(FpolAll).T
            PresAll = np.array(PresAll).T
            FFprimeAll = np.array(FFprimeAll).T
            PprimeAll = np.array(PprimeAll).T
            qpsiAll = np.array(qpsiAll).T
        else:
            FpolAll = np.array(FpolAll)
            PresAll = np.array(PresAll)
            FFprimeAll = np.array(FFprimeAll)
            PprimeAll = np.array(PprimeAll)
            qpsiAll = np.array(qpsiAll)            
#       FpolAll = np.dstack(FpolAll)
#       FFprimeAll = np.dstack(FFprimeAll)
#       PprimeAll = np.dstack(PprimeAll)
#       qpsiAll = np.dstack(qpsiAll)

        #Set up interpolators
        RmAxisInterp = interp1d(Spols,RmAxisAll)
        ZmAxisInterp = interp1d(Spols,ZmAxisAll)
        psiRZInterp = RegularGridInterpolator((R, Z, Spols), psiRZAll)
        psiAxisInterp = interp1d(Spols, psiAxisAll)
        psiSepInterp = interp1d(Spols, psiSepAll)
        Bt0Interp = interp1d(Spols, Bt0All)
        IpInterp = interp1d(Spols, IpAll)
        psiN = np.linspace(0,1,len(FpolAll[:,0]))
        FpolInterp = RegularGridInterpolator((psiN,Spols), FpolAll)
        PresInterp = RegularGridInterpolator((psiN,Spols), PresAll)
        FFprimeInterp = RegularGridInterpolator((psiN,Spols), FFprimeAll)
        PprimeInterp = RegularGridInterpolator((psiN,Spols), PprimeAll)
        qpsiInterp = RegularGridInterpolator((psiN,Spols), qpsiAll)

        #Interpolate each parameter for the new timestep
        r,z = np.meshgrid(R,Z)
        RmAxis = RmAxisInterp(newS)
        ZmAxis = ZmAxisInterp(newS)
        if transposePsi==True:
            psiRZ = psiRZInterp((r,z,newS))
        else:
            psiRZ = psiRZInterp((r,z,newS)).T
        psiAxis = psiAxisInterp(newS)
        psiSep = psiSepInterp(newS)
        Bt0 = Bt0Interp(newS)
        Ip = IpInterp(newS)
        Fpol = FpolInterp((psiN,newS))
        Pres = PresInterp((psiN,newS))
        FFprime = FFprimeInterp((psiN,newS))
        Pprime = PprimeInterp((psiN,newS))
        qpsi = qpsiInterp((psiN,newS))

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

    #=====================================================================
    #                       LEGACY FUNCTIONS (not used now)
    #=====================================================================
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
        file = gpath + 'g' + self.shotFmt.format(shot) + '_' + self.tsFmt.format(time)
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

    def copyGfile2tree(self,gFileName,shot,time,clobberflag=True):
        """
        Copies gfile to HEAT tree
        gFileName is name of gFile that is already located in self.tmpDir
        """
        oldgfile = self.tmpDir + gFileName
        name = 'g'+self.shotFmt.format(shot)+'_'+self.tsFmt.format(time)
        #make tree for this shot
        tools.makeDir(self.shotPath, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        #make tree for this timestep
        timeDir = self.shotPath + self.tsFmt.format(time) + '/'
        newgfile = timeDir + name

        #clobber and make time directory
        tools.makeDir(timeDir, clobberFlag=False, mode=self.chmod, UID=self.UID, GID=self.GID)
        shutil.copyfile(oldgfile, newgfile)
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

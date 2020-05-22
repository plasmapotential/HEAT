#MHDClass.py
#Description:   Base HEAT MHD Equilibrium (EQ) module
#Engineer:      T Looby
#Date:          20191111

import EFIT.equilParams_class as EP
import gfiles
import toolsClass
tools = toolsClass.tools()

import time
import numpy as np
import pandas as pd
import os

import logging
log = logging.getLogger(__name__)

class MHD:

    def __init__(self):
        pass

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


        self.allowed_vars = ['testvar',
                            'shot',
                            'tree',
                            'tmin',
                            'tmax',
                            'nTrace',
                            'dataPath',
                            'MachFlag',
                            'Nphi',
                            'ittLaminar',
                            'ittStruct',
                            'Smin',
                            'Smax',
                            'phimin',
                            'phimax',
                            'Nswall',
                            'phistart',
                            'MapDirection',
                            'PlasmaResponse',
                            'Field',
                            'target',
                            'createPoints',
                            'useECcoil',
                            'useIcoil',
                            'useCcoil',
                            'useFcoil',
                            'useBcoil',
                            'useBus',
                            'useFilament',
                            'useTe_profile',
                            'ParticleDirection',
                            'ParticleCharge',
                            'Ekin',
                            'Lambda',
                            'Mass',
                            'Zmin',
                            'Zmax',
                            'Rmin',
                            'Rmax',
                            'dataPath']
        return


    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        self.shot = int(self.shot)
        self.tmin = int(self.tmin)
        self.tmax = int(self.tmax)
        self.Nphi = int(self.Nphi)
        self.Smin = float(self.Smin)
        self.Smax = float(self.Smax)
        self.phimin = float(self.phimin)
        self.phimax = float(self.phimax)
        self.Nswall = int(self.Nswall)
        self.ittLaminar = float(self.ittLaminar)
        self.ittStruct = float(self.ittStruct)
        self.phistart = float(self.phistart)
        self.MapDirection = int(self.MapDirection)
        self.PlasmaResponse = int(self.PlasmaResponse)
        self.Field = int(self.Field)
        self.target = int(self.target)
        self.createPoints = int(self.createPoints)
        self.Zmin = float(self.Zmin)
        self.Zmax = float(self.Zmax)
        self.Rmin = float(self.Rmin)
        self.Rmax = float(self.Rmax)
        if(self.MachFlag == 'iter'):
            self.useIcoil = int(self.useIcoil)
        elif(self.MachFlag == 'nstx'):
            self.useECcoil = int(self.useECcoil)
        elif(self.MachFlag == 'st40'):
            self.useECcoil = int(self.useECcoil)
        elif(self.MachFlag == 'mast'):
            self.useCcoil = int(self.useCcoil)
            self.useIcoil = int(self.useIcoil)
        elif(self.MachFlag == 'd3d'):
            self.useFcoil = int(self.useFcoil)
            self.useCcoil = int(self.useCcoil)
            self.useIcoil = int(self.useIcoil)
        self.useFilament = int(self.useFilament)
        self.useTe_profile = int(self.useTe_profile)
        self.ParticleDirection = int(self.ParticleDirection)
        self.ParticleCharge = int(self.ParticleCharge)
        self.Ekin = float(self.Ekin)
        self.Lambda = float(self.Lambda)
        self.Mass = float(self.Mass)
        if self.MachFlag in ['dt']:
            self.useBus = int(self.useBus)
            self.useBcoil = int(self.useBcoil)
        return

    # Pull fresh gfile data from MDS+ tree, create directory tree in the process
    def get_mhd_inputs(self,machine='nstx',gfile=None):
        """
        get gfile from mds+ tree if gfile is None, otherwise use file
        """
        if gfile==None:
            self.timesteps, self.preserveFlag = gfiles.write_multiple_gfiles(
                                                                self.MachFlag,
                                                                self.shot,
                                                                self.tree,
                                                                self.tmin,
                                                                self.tmax,
                                                                self.dataPath,
                                                                clobberwait=False)
        else:
            self.timesteps, self.preserveFlag = gfiles.loadgfile(machine,gfile,
                                        rootDir=self.dataPath, clobberwait=False)





    def make1EFITobject(self,shot,t,gfile=None):
        """
        Creates an equilParams_class object for a SINGLE timesteps. equilParams
        is a class from the ORNL_Fusion github repo, and was developed by
        A. Wingen
        """
        if self.dataPath[-1]!='/':
            dataPath = self.dataPath+'/'
        if gfile is None:
            gfile = dataPath+'{:06d}/g{:06d}.{:05d}'.format(t,shot,t)
        self.ep = EP.equilParams(gfile)
        #Correct for weird helicity and field directions that are occasionally
        #in gfile.  Here we assume that Ip is the CCW direction as viewed from
        #above tokamak.
        self.dsign = np.sign(self.ep.g['Ip'])
        self.gsign = np.sign( (self.ep.g['psiAxis'] - self.ep.g['psiSep']) )
        self.qsign = np.sign(self.ep.g['Fpol'][-1]) #F_edge sign
        print('dsign = {:f}'.format(self.dsign))
        print('gsign = {:f}'.format(self.gsign))
        print('qsign = {:f}'.format(self.qsign))
        self.ep.g['FFprime'] *= self.dsign*self.qsign*self.gsign
        self.ep.g['Pprime'] *= self.dsign*self.qsign*self.gsign
        self.ep.g['psiRZ'] *= self.dsign*self.qsign*self.gsign
        self.ep.g['qpsi'] *= -self.qsign*self.dsign
        self.ep.g['psiSep'] *= self.dsign*self.qsign*self.gsign
        self.ep.g['psiAxis'] *= self.dsign*self.qsign*self.gsign
        self.ep.g['Fpol'] *= self.dsign*self.qsign
        return

    def makeEFITobjects(self):
        """
        Creates an equilParams_class object for MULTIPLE timesteps. equilParams
        is a class from the ORNL_Fusion github repo, and was developed by
        A. Wingen
        """
        self.ep= ['None' for i in range(len(self.timesteps))]
        for idx,t in enumerate(self.timesteps):
            if self.dataPath[-1]!='/':
                dataPath = self.dataPath+'/'
            gfile = dataPath+'{:06d}/g{:06d}.{:05d}'.format(t,self.shot,t)
            self.ep[idx] = EP.equilParams(gfile)
        return

    def checkEFITobject(self):
        """
        checks to make sure that the EFIT object created from a gfile
        has the correct signs for poloidal flux, toroidal field, plasma current,
        etc., so that the helicity and field are in the correct orientation
        """






    def Bfield_pointcloud(self, ep, R, Z, phi, normal=False):
        """
        Creates a Bfield Pointcloud that can be saved in a .csv file
        This pointcloud is then used by ParaView to overlay the
        Bfield onto the CAD for analysis.  The format is:
        X,Y,Z,Bx,By,Bz
        where X,Y,Z are the spatial locations of the Bx,By,Bz vectors
        In paraview use TableToPoints => Calculator => Glyph
        Calculator should have this formula:
        (iHat*Bx) + (jHat*By) + (kHat*Bz)

        Note: this function returns the normal Bfield components if normal=True

        """
        #BR = ep.BRFunc.ev(R,Z)
        #BZ = ep.BZFunc.ev(R,Z)
#        Bt = ep.BtFunc.ev(R,Z)*self.MapDirection
#        BR = ep.BRFunc.ev(R,Z)
#        BZ = -ep.BZFunc.ev(R,Z)
        Bt = ep.BtFunc.ev(R,Z)#*self.MapDirection
        BR = ep.BRFunc.ev(R,Z)#*self.MapDirection
        BZ = ep.BZFunc.ev(R,Z)#*self.MapDirection

        #print(BR)
        #print(Bt)
        #print(BZ)
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

    def write_B_pointcloud(self,centers,Bxyz, dataPath):
        """
        In paraview use TableToPoints => Calculator => Glyph
        Calculator should have this formula:
        (iHat*Bx) + (jHat*By) + (kHat*Bz)
        """
        print("Creating B Field Point Cloud")
        log.info("Creating B Field Point Cloud")
        pcfile = dataPath + 'B_pointcloud.csv'
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
        tools.createVTKOutput(pcfile, 'glyph', 'B_pointcloud')
        return


    def writeControlFile(self, name, t, mode='laminar'):
        """
        Create a control file for MAFOT
        """
        save_location = self.dataPath + '/' + '{:06d}/'.format(t) + name
        with open(save_location, 'w') as f:

            f.write('# Parameterfile for ' + self.MachFlag + ' Programs\n')
            f.write('# Shot: {:06d}\tTime: {:05d}ms\n'.format(int(self.shot), int(t)))
            f.write('# Path: ' + self.dataPath + '/{:06d}\n'.format(t))

            f.write('Nphi=\t{:d}\n'.format(self.Nphi))

            #itt means different things depending on if we are tracing field line
            #or running full MAFOT laminar
            if mode=='laminar':
                f.write('itt=\t{:f}\n'.format(self.ittLaminar))
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
            f.write('MapDirection=\t{:d}\n'.format(self.MapDirection))


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
            elif(self.MachFlag == 'st40'):
                f.write('useECcoil(0=no,1=yes)=\t{:d}\n'.format(self.useECcoil))
            elif(self.MachFlag == 'mast'):
                f.write('useCcoil(0=no,1=yes)=\t{:d}\n'.format(self.useCcoil))
                f.write('useIcoil(0=no,1=yes)=\t{:d}\n'.format(self.useIcoil))
            elif(self.MachFlag == 'd3d'):
                f.write('useFcoil(0=no,1=yes)=\t{:d}\n'.format(self.useFcoil))
                f.write('useCcoil(0=no,1=yes)=\t{:d}\n'.format(self.useCcoil))
                f.write('useIcoil(0=no,1=yes)=\t{:d}\n'.format(self.useIcoil))

            if self.MachFlag in ['iter', 'nstx', 'mast', 'st40', 'd3d']:
                f.write('useFilament(0=no)=\t{:d}\n'.format(self.useFilament))
            if self.MachFlag in ['iter', 'nstx', 'st40', 'd3d']:
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


    def writeMAFOTpointfile(self,xyz,gridfile):
        """
        Creates a file that MAFOT's can use to trace a field line.  Points
        in this file are the starting points for traces, organized in columns
        of R,phi,Z
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
                    controlfile, paraview_mask=False):
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
        cmd = self.MachFlag + 'structure '
        flags = '-d {:f}'.format(int(dphi)) + ' '
        flags += '-P ' + gridfile + ' '
        #flags += '-Z '
        flags += controlfile
        args = [cmd,'-d','{:f}'.format(dphi),'-P',gridfile,'-Z',controlfile]
        #print("COMMAND:")
        #print(cmd+flags)
        #Run shell command.  Will create file with results
        from subprocess import call
        call(cmd + flags, shell = True, cwd=controlfilePath)
        call('rm log_*', shell = True, cwd=controlfilePath)	# clean up

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
            np.savetxt(controlfilePath+'struct.csv', xyz, delimiter=',', header=head)
            print('Converted file to ParaView formatted CSV.')
            log.info('Converted file to ParaView formatted CSV.')

            #Now save a vtk file for paraviewweb
            tools.createVTKOutput(controlfilePath+'struct.csv', 'trace', 'Field_trace')

        return

    def getMultipleFieldPaths(self, dphi, gridfile, controlfilePath,
                    controlfile, paraview_mask=False):

        if self.MachFlag == 'd3d':
            cmd = 'dtstructure '
        else:
            cmd = self.MachFlag + 'structure '
        flags = '-d {:f}'.format(float(dphi)) + ' '
        flags += '-P ' + gridfile + ' '
#        flags += '-Z '
        flags += controlfile
        #Run shell command.  Will create file with results
        from subprocess import call
        call(cmd + flags, shell = True, cwd=controlfilePath)
        #call('rm log_*', shell = True, cwd=controlfilePath)	# clean up

    def runMAFOTlaminar(self,gridfile,controlfilePath,controlfile,NCPUs):
        """
        Runs MAFOT laminar program to trace field lines and return connection
        length, psi, penetration depth, etc.
        """
        #Create MAFOT shell command
        if self.MachFlag == 'd3d':
            tool = 'dtlaminar_mpi'
        else:
            tool = self.MachFlag + 'laminar_mpi'
        cmd = 'mpirun'
        args = [cmd,'-n','{:d}'.format(NCPUs),tool,controlfile,'-P',gridfile]
        from subprocess import PIPE
        from subprocess import Popen
        #p = Popen(args,stdin=PIPE,stdout=PIPE,cwd=controlfilePath)
        #print("Finished MAFOT Laminar run")

        from subprocess import call
        flags = ' -n {:d} '.format(NCPUs) + tool + ' ' + controlfile + ' ' + '-P '+gridfile
        call(cmd + flags, shell = True, cwd=controlfilePath)
        #call('rm log_*', shell = True, cwd=controlfilePath)	# clean up
        #os.remove(controlfilePath+'log_*') #clean up

        return

    def writeGfile(self, file, shot, time):
        """
        writes a new gfile.  for use with the cleaner script.  user must supply
        file: name of new gfile
        shot: new shot number
        time: new shot timestep [ms]
        """
        g = self.ep.g
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

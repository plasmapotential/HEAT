import numpy as np
import plotly.graph_objects as go
import toolsClass
tools = toolsClass.tools()

import logging
log = logging.getLogger(__name__)


class GYRO:

    def __init__(self, rootDir, dataPath):
        """
        rootDir is root location of python modules (where dashGUI.py lives)
        dataPath is the location where we write all output to
        """
        self.rootDir = rootDir
        tools.rootDir = self.rootDir
        self.dataPath = dataPath
        tools.dataPath = self.dataPath
        return

    def setupConstants(self, species='D'):
        """
        Sets up constants

        if species is 'D', then it is set up as deuterium
        """
        #unit conversions
        self.kg2eV = 5.609e35 #1kg = 5.609e35 eV/c^2
        self.eV2K = 1.160e4 #1ev=1.160e4 K
        #constants
        self.AMU = 931.494e6 #ev/c^2
        self.kB = 8.617e-5 #ev/K
        self.e = 1.602e-19 # C
        self.c = 299792458 #m/s

        #deuterium
        if species == 'D':
            print("Using Deuterium Properties")
            self.mass_eV = 2.014*self.AMU  #eV/c^2
            self.Z = 1 #1 for deuterium
        #tritium
        elif species == 'T':
            print("Using Tritium Properties")
            self.mass_eV = 3.016*self.AMU  #eV/c^2
            self.Z = 1 #1 for deuterium
        #hydrogen
        else:
            print("Using Hydrogen Properties")
            self.mass_eV = 1.007*self.AMU #eV/c^2
            self.Z = 1 #1 for deuterium
        return

    def temp2thermalVelocity(self, T_eV):
        """
        Calculates thermal velocity from a temperature

        T_eV is temperature in eV
        """
        return np.sqrt(2*float(T_eV)/(self.mass_eV/self.c**2))


    def setupFreqs(self, B, v_perp0):
        """
        Calculates frequencies, periods, radii, that are dependent upon B and vPerp
        These definitions follow Freidberg Section 7.7

        B is magnetic field magnitude

        """
        #User defined but can be set here for tests
        #B = 5 #tesla
        #T0_eV = 15000 # eV

        self.omegaGyro = self.Z * self.e * B / (self.mass_eV / self.kg2eV)
        self.fGyro = self.omegaGyro/(2*np.pi)
        self.TGyro = 1.0/self.fGyro
        #gyro radius
        self.rGyro = v_perp0 / self.omegaGyro
        return


    def randomVelocity(self, T0, N=1):
        """
        Monte Carlo sampling of Maxwellian velocity distribution

        Returns N velocities distributed about T0 per the Maxwellian distribution
        """
        return

    def randomPhaseAngle(self, N=1):
        """
        Monte Carlo sampling of a uniform distribution between 0 and 2pi

        returns an array of N random values in radians
        """
        gyroPhase = np.random.uniform(0, 2*np.pi, N)
        return gyroPhase

    def singleGyroTrace(self,vPerp,vParallel,gyroPhase,N_gyroSteps,
                        BtraceXYZ,controlfilePath,verbose=True):
        """
        Calculates the gyro-Orbit path and saves to .csv and .vtk

        vPerp and vParallel [m/s] are in velocities
        gyroPhase [degrees] is initial orbit phase angle
        N_gyroSteps is number of discrete line segments per gyro period
        BtraceXYZ is the points of the Bfield trace that we will gyrate about
        """
        print("Calculating gyro trace...")
        #Loop thru B field trace while tracing gyro orbit
        helixTrace = None
        for i in range(len(BtraceXYZ)-1):
            #points in this iteration
            p0 = BtraceXYZ[i,:]
            p1 = BtraceXYZ[i+1,:]
            #vector
            delP = p1 - p0
            #magnitude
            magP = np.sqrt(delP[0]**2 + delP[1]**2 + delP[2]**2)
            #time it takes to transit line segment
            delta_t = magP / (vParallel)
            #Number of steps in line segment
            Tsample = self.TGyro / N_gyroSteps
            Nsteps = int(delta_t / Tsample)
            #length (in time) along guiding center
            t = np.linspace(0,delta_t,Nsteps+1)
            #guiding center location
            xGC = np.linspace(p0[0],p1[0],Nsteps+1)
            yGC = np.linspace(p0[1],p1[1],Nsteps+1)
            zGC = np.linspace(p0[2],p1[2],Nsteps+1)
            # construct orthogonal system for coordinate transformation
            w = delP
            if np.all(w==[0,0,1]):
                u = np.cross(w,[0,1,0]) #prevent failure if bhat = [0,0,1]
            else:
                u = np.cross(w,[0,0,1]) #this would fail if bhat = [0,0,1] (rare)
            v = np.cross(w,u)
            #normalize
            u = u / np.sqrt(u.dot(u))
            v = v / np.sqrt(v.dot(v))
            w = w / np.sqrt(w.dot(w))
            xfm = np.vstack([u,v,w]).T
            #get helix path along (proxy) z axis reference frame
            x_helix = self.rGyro*np.cos(self.omegaGyro*t + gyroPhase)
            y_helix = self.rGyro*np.sin(self.omegaGyro*t + gyroPhase)
            z_helix = np.zeros((len(t)))
            #perform rotation to field line reference frame
            helix = np.vstack([x_helix,y_helix,z_helix]).T
            helix_rot = np.zeros((len(helix),3))
            for j,coord in enumerate(helix):
                helix_rot[j,:] = helix[j,0]*u + helix[j,1]*v + helix[j,2]*w
            #perform translation to field line reference frame
            helix_rot[:,0] += xGC
            helix_rot[:,1] += yGC
            helix_rot[:,2] += zGC
            #update gyroPhase variable so next iteration starts here
            gyroPhase = self.omegaGyro*t[-1] + gyroPhase
            #append to helix trace
            if helixTrace is None:
                helixTrace = helix_rot
            else:
                helixTrace = np.vstack([helixTrace,helix_rot])

        helixTrace*=1000.0 #scale for ParaView
        print("Saving data to CSV and VTK formats")
        #save data to csv format
        head = 'X[mm],Y[mm],Z[mm]'
        np.savetxt(controlfilePath+'helix.csv', helixTrace, delimiter=',', header=head)
        #save data to vtk format
        tools.createVTKOutput(controlfilePath+'helix.csv', 'trace', 'Gyro_trace')

        if verbose==True:
            print("V_perp = {:f} [m/s]".format(vPerp))
            print("V_parallel = {:f} [m/s]".format(vParallel))
            print("Cyclotron Freq = {:f} [rad/s]".format(self.omegaGyro))
            print("Cyclotron Freq = {:f} [Hz]".format(self.fGyro))
            print("Gyro Radius = {:f} [m]".format(self.rGyro))
            print("Number of gyro points = {:f}".format(len(helixTrace)))
            print("Longitudinal dist between gyro points = {:f} [m]".format(magP/float(Nsteps)))
            print("Each line segment length ~ {:f} [m]".format(magP))
        return

#heatfluxClass.py
#Description:   Heat Flux Module
#Engineer:      T Looby
#Date:          20191117
import numpy as np
import toolsClass
import os
import sys
import time
tools = toolsClass.tools()

import EFIT.equilParams_class as EP
import scipy.interpolate as scinter
from scipy.optimize import bisect
from scipy.integrate import simps

import multiprocessing


import logging
log = logging.getLogger(__name__)

class heatFlux:

    def __init__(self):
        pass

    def allowed_class_vars(self):
        """
        Writes a list of recognized class variables to HEAT object
        Used for error checking input files and for initialization

        Here is a list of variables with description:
        testvar         dummy for testing

        """


        self.allowed_vars = ['testvar',
                            'profileType',
                            'lqEich',
                            'S',
                            'Psol',
                            'qBG',
                            'lqPN',
                            'lqPF',
                            'lqCN',
                            'lqCF',
                            'fracPN',
                            'fracPF',
                            'fracCN',
                            'fracCF']
        return

    def setTypes(self):
        """
        Set variable types for the stuff that isnt a string from the input file
        """
        self.lqEich = float(self.lqEich)
        self.S = float(self.S)
        self.Psol = float(self.Psol)
        self.qBG = float(self.qBG)
        self.lqPN = float(self.lqPN)
        self.lqPF = float(self.lqPF)
        self.lqCN = float(self.lqCN)
        self.lqCF = float(self.lqCF)
        self.fracPN = float(self.fracPN)
        self.fracPF = float(self.fracPF)
        self.fracCN = float(self.fracCN)
        self.fracCF = float(self.fracCF)
        return

    def makePFCs(self,MHDobj,sourceMeshes,sourceCenters,sourceNorms,sourceAreas,names):
        """
        Initialize a list of heat flux objects, corresponding to each part mesh
        passed to this function.  sourceMeshes should be a list of FreeCAD
        meshes that we want to find the heat fluxes on.
        """
        #Check if this is a single mesh or list and make it a list
        if type(sourceMeshes) != list:
            sourceMeshes = [sourceMeshes]
        if type(sourceCenters) != list:
            sourceCenters = [sourceCenters]
        if type(sourceNorms) != list:
            sourceNorms = [sourceNorms]
        if type(names) != list:
            names = [names]


        self.PFCs = []
        for i,source in enumerate(sourceMeshes):
            self.PFCs.append(self.PFC(MHDobj,
                                      source,
                                      sourceCenters[i],
                                      sourceNorms[i],
                                      sourceAreas[i],
                                      name=names[i]))

    def backfaceCulling(self,centers,norms,MHD):
        """
        Use backface culling technique from computer graphics to eliminate
        faces that are on the 'back side' of the PFC tiles.  Mathematically,
        if: B_hat dot n_hat >= 0, then the face is shadowed.
        returns 1 if face is shadowed, 0 if face is not shadowed
        """
        xyz = centers
        r,z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        BNorms = MHD.Bfield_pointcloud(MHD.ep, r, z, phi, normal=True)
        dot = np.multiply(norms, BNorms).sum(1)
        shadowed_mask = np.zeros((len(centers)))
        shadowed_mask[np.where(dot >= 0.0)] = 1
        return shadowed_mask

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

    def readMAFOTLaminarOutput(self,PFC,file):
        """
        Reads output from MAFOT laminar into PFC object
        """
        print('Reading laminar output')
        log.info('Reading laminar output')
        data = np.genfromtxt(file,comments='#')
        use = np.where(PFC.shadowed_mask != 1)[0]
        xyz = PFC.centers[use]
        r,z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        r = np.round(r,10) #number of decimal places in MAFOT output file

        #Sometimes the MAFOT calculation returns an error and it discards the
        #launch point instead of writing to the output file.  This results
        #in less points coming out of MAFOT than we put in.  So we amend this
        #by checking to make sure the points coming out are assigned to the correct
        #location in the mesh centers.
        #Could be parallelized in future

        #Create shared variables so that all multiprocessing cores can access
        #same variables without needless data transmission
        tools.lamData = data
        tools.lamR = r

        #Prepare intersectionTest across multiple cores
        t0 = time.time()
        Ncores = multiprocessing.cpu_count() -2 #reserve 2 cores for overhead
        print('Initializing parallel MAFOT laminar check across {:d} cores'.format(Ncores))
        log.info('Initializing parallel MAFOT laminar check across {:d} cores'.format(Ncores))
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
        pool = multiprocessing.Pool(Ncores)
        indexes = np.asarray(pool.map(tools.readLaminarParallel, np.arange(len(r))))
        pool.close()

        PFC.psimin = data[indexes,4]
        PFC.conLength = data[indexes,3]

        #R = data[:,0]
        #Z = data[:,1]
        #phi = np.radians(data[:,9])
        #x,y,z = tools.cyl2xyz(R,Z,phi)
        #PFC.centersLam = np.concatenate((x,y,z)).reshape((-1, 3), order='F')
        print('Laminar output read')
        print('Requested {:d} traces from MAFOT'.format(len(r)))
        print('MAFOT returned {:d} traces'.format(len(data[:,0])))
        log.info('Laminar output read')
        log.info('Requested {:d} traces from MAFOT'.format(len(r)))
        log.info('MAFOT returned {:d} traces'.format(len(data[:,0])))

        return

    def map_R_psi(self, psi, PFC):
        """
        Map normalized poloidal flux psi to R at midplane (Z = 0)
        psi is flat array
        return R(psi)
        """
        R = np.linspace(PFC.ep.g['RmAxis'], PFC.ep.g['R1'] + PFC.ep.g['Xdim'], 100)
        Z = np.zeros(len(R))
        p = PFC.ep.psiFunc.ev(R,Z)

        #In case of monotonically decreasing psi, sort R, p so that p is
        #monotonically increasing
        points = zip(R, p)
        points = sorted(points, key=lambda point: point[1]) # Sort list of tuples by p-value (psi)
        R, p = zip(*points)

        f = scinter.UnivariateSpline(p, R, s = 0, ext = 'const')	# psi outside of spline domain return the boundary value
        return f(psi)

    def getEichFromEQ(self, ep, MachFlag):
        """
        finds lqEich and S from equilibrium object using regression from
        Eich's paper:

        T. Eich et al., Nucl. Fusion, vol. 53, no. 9, p. 093031, Sep. 2013,
        doi: 10.1088/0029-5515/53/9/093031.

        if machine is NSTX, then use regression #15
        if machine is D3D, then use regression #10
        if machine is ST40, then use regression #15

        """
        # geometric quantities
        Rmax = ep.g['lcfs'][:,0].max()
        Rmin = ep.g['lcfs'][:,0].min()
        Rgeo = (Rmax + Rmin) / 2.0
        a = (Rmax - Rmin) / 2.0
        aspect = a/Rgeo

        if MachFlag == 'nstx':
            #Regression 15
            C = 1.35
            Cp = -0.02
            Cr = 0.04
            Cb = -0.92
            Ca = 0.42

            # Evaluate B at outboard midplane
            Z_omp_sol = 0.0
            Bp = PFC.ep.BpFunc.ev(Rmax,Z_omp_sol)

            #Evaluate lq, S
            self.lqEich = C * self.Psol**Cp * Rgeo**Cr * Bp**Cb * aspect**Ca # in mm
            self.S = 1.0 #in mm

        return




    def getHFprofile(self, PFC, MachFlag):
        """
        Calculates heat flux profile from psi.  Default is an Eich profile.

        mode can be 'eich', 'multiExp', 'limiter'

        """
        psi = PFC.psimin
        R_omp = self.map_R_psi(psi,PFC)
        Z_omp = np.zeros(R_omp.shape)
        # Evaluate B at midplane
        Bp_omp = PFC.ep.BpFunc.ev(R_omp,Z_omp)
        Bt_omp = PFC.ep.BtFunc.ev(R_omp,Z_omp)
        B_omp = np.sqrt(Bp_omp**2 + Bt_omp**2)
        xyz = PFC.centers
        R_div,Z_div,phi_div = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        print('phi_divMin = {:f}'.format(phi_div.min()))
        print('phi_divMax = {:f}'.format(phi_div.max()))
        # Evaluate B at Target Plate neglecting shadowed points
        Bp_div = PFC.ep.BpFunc.ev(R_div,Z_div)
        Bt_div = PFC.ep.BtFunc.ev(R_div,Z_div)
        B_div = np.sqrt(Bp_div**2 + Bt_div**2)
        #Calculate psi using gfile for scaling coefficient
        psi_EQ = PFC.ep.psiFunc.ev(R_div,Z_div)
        #Calculate poloidal flux expansion
        #fx = R_div*Bp_div / (R_omp*Bp_omp)
        q = np.zeros(PFC.centers[:,0].shape)
        use = np.where(PFC.shadowed_mask == 0)[0]

        #Multiple exponential profile (Brunner Profile)
        if self.mode=='multiExp' or self.mode=='limiter':
            q[use] = self.multiExp_profile_fluxspace(PFC, R_omp, Bp_omp, psi, self.mode)

        #Eich Profile
        else:
            #if Eich profile lq and S were not defined find them from equilibrium object
            if self.lqEich == None:
                self.getEichFromEQ(PFC.ep, MachFlag)
            q0 = self.scaleHF_fluxspace(PFC,self.lqEich,self.S,self.Psol)
            q[use] = self.eich_profile_fluxspace(PFC, self.lqEich, self.S, R_omp, Bp_omp, psi)
            q *= q0
            q += self.qBG

        return q

    def scaleHF_fluxspace(self, PFC, lqEich, S, P):
        """
        scales HF using Eich profile

        Get scale factor q||0 (q0) for heat flux via power balance:
        (input MW = output MW)
        Ignores wall psi and just creates a profile at OMP
        Creates a dense (1000pts) grid at the midplane to get higher resolution
        integral.  Integrates q_hat / B_omp with respect to psi.
        q||0 = P_div / ( 2*pi* integral(q_hat / B_omp)dPsi )

        return q0
        """
        # Get R and Z vectors at the midplane
        R_omp_sol = PFC.ep.g['lcfs'][:,0].max()
        R_omp_min = R_omp_sol - 5.0*lqEich*(1e-3) #in meters now
        R_omp_max = R_omp_sol + 20.0*lqEich*(1e-3) #in meters now
        R_omp = np.linspace(R_omp_min, R_omp_max, 1000)
        Z_omp = np.zeros(R_omp.shape)
        #Calculate flux at midplane using gfile
        psiN = PFC.ep.psiFunc.ev(R_omp,Z_omp)
        psi = PFC.ep.psiFunc_noN.ev(R_omp,Z_omp)
        PFC.psiMinLCFS = PFC.ep.psiFunc.ev(R_omp_sol,0.0)
        s_hat = psiN - PFC.psiMinLCFS
        # Evaluate B at outboard midplane
        Bp_omp = PFC.ep.BpFunc.ev(R_omp,Z_omp)
        Bt_omp = PFC.ep.BtFunc.ev(R_omp,Z_omp)
        B_omp = np.sqrt(Bp_omp**2 + Bt_omp**2)

        #Get q|| profile then integrate in Psi
        #Eich profile
        q_hat = self.eich_profile_fluxspace(PFC, lqEich, S, R_omp, Bp_omp, psiN)
        P0 = 2*np.pi * simps(q_hat / B_omp, psi)
        if P0 < 0: P0 = -P0
        #Scale to input power
        q0 = P/P0
        return q0



    #================== Eich Profile ============================
    def eich_profile_fluxspace(self, PFC, lq, S, R, Bp, psiN):
        """
        Based on the paper: T.Eich et al.,PRL 107, 215001 (2011)

        Here we adapt so that q is a function of normalized poloidal flux (psi)
        rather than distance from separatrix ie: q(psi).  This requires the
        coordinate transformation, xfm

        lq is heat flux width at midplane in mm
        lq_hat is heat flux width in flux coordinates
        s_hat is a flux coordinate
        S is the private flux region spreading in mm

        return array q1(psi)
        """
        from scipy.special import erfc
        # Convert to meters
        lq *= 1e-3
        S *= 1e-3
        psiaxis = PFC.ep.g['psiAxis']
        psiedge = PFC.ep.g['psiSep']
        deltaPsi = np.abs(psiedge - psiaxis)
        s_hat = psiN - PFC.psiMinLCFS
        # Gradient
        gradPsi = Bp*R
        xfm = gradPsi / deltaPsi
        # Decay width mapped to flux coordinates
        lq_hat = lq * xfm
        rho = s_hat/lq_hat
        rho_0 = S/(2.0*lq)
        #===Eich Profile as a function of psi
        q1 = 0.5 * np.exp(rho_0**2 - rho) * erfc(rho_0 - rho/(2*rho_0))
        nan_locations = np.isnan(q1)
        inf_locations = np.isinf(q1)
        q1[nan_locations] = 0.0
        q1[inf_locations] = 0.0
        return q1

    #================== Multiple Exponential Profile ============================
    def multiExp_profile_fluxspace(self, PFC, R, Bp, psiN, mode):
        """
        Multiple (4) exponential scaling for divertor plasmas is based upon
        Brunner's multiple exponential regression:
        D. Brunner, Nucl. Fusion, vol. 58, no. 9, p. 094002, Sep. 2018,
        doi: 10.1088/1741-4326/aad0d6.

        Double (2) exponential scaling for limiter plasmas follows same profile
        as Brunner's scaling, but without the private plasma exponential decay
        (because it's limited = no private plasma), as described in this paper:
        J. Horacek et al.,Plasma Phys. Control. Fusion, vol. 58, no. 7,
        p. 074005, Jul. 2016, doi: 10.1088/0741-3335/58/7/074005.

        Here we adapt so that q is a function of normalized poloidal flux (psi)
        rather than distance from separatrix ie: q(psi).  This requires the
        coordinate transformation, xfm

        A multiple exponential scaling is apparent in many tokamaks, with a
        'near' and 'far' SOL in both the common flux and private flux regions.
        We define each as follows:
        lqPN        near SOL decay length in private plasma
        lqPF        far SOL decay length in private plasma
        lqCN        near SOL decay length in common plasma
        lqCF        far SOL decay length in common plasma

        fracPN is fraction of power sharing between near exponential in private flux region
        fracPF is fraction of power sharing between far exponential in private flux region
        fracCN is fraction of power sharing between near exponential in common flux region
        fracCF is fraction of power sharing between far exponential in common flux region
        1 = fracNP + fracFP + fracNC + fracFC

        return array q1(psi)
        """
        # Convert to meters
        lqCN = self.lqCN*1e-3
        lqCF = self.lqCF*1e-3
        if mode=='multiExp':
            lqPN = self.lqPN*1e-3
            lqPF = self.lqPF*1e-3

        psiaxis = PFC.ep.g['psiAxis']
        psiedge = PFC.ep.g['psiSep']
        deltaPsi = np.abs(psiedge - psiaxis)
        # Gradient
        gradPsi = Bp*R
        xfm = gradPsi / deltaPsi

        # transform into flux space
        lqCN_hat = lqCN*xfm
        lqCF_hat = lqCF*xfm
        if mode=='multiExp':
            lqPN_hat = lqPN*xfm
            lqPF_hat = lqPF*xfm

        if mode=='multiExp':
            q0 = self.findScalingCoeffsMultiExp(PFC, lqCN, lqCF, lqPN, lqPF)
        else:
            q0 = self.findScalingCoeffsLimiter(PFC, lqCN, lqCF)

        s_hat = psiN - PFC.psiMinLCFS
        #find locations in Private vs Common flux regions
        useP = np.where(s_hat < 0.0)[0]
        useC = np.where(s_hat >= 0.0)[0]
        print('{:d} points in private flux region'.format(len(useP)) )
        print('{:d} points in common flux region'.format(len(useC)) )

        q = np.zeros(len(psiN))
        #===Brunner Profile as a function of psi
        if mode=='multiExp':
            if len(useC)>0:
                q[useC] = self.fracCN*np.exp(-s_hat[useC] / lqCN_hat[useC]) + self.fracCF*np.exp(-s_hat[useC] / lqCF_hat[useC])
            if len(useP)>0:
                q[useP] = self.fracPN*np.exp(s_hat[useP] / lqPN_hat[useP]) + self.fracPF*np.exp(s_hat[useP] / lqPF_hat[useP])
            q *= q0

        #===limiter Profile as a function of psi
        #technically there should only be one contact point between SOL and wall
        #but this may not be true in simulation all the time, so we just leave
        #any points inside the core (s_hat < 0) as q=0.
        else:
            q[useC] = self.fracCN*np.exp(-s_hat[useC] / lqCN_hat[useC]) + self.fracCF*np.exp(-s_hat[useC] / lqCF_hat[useC])
            q *= q0

        nan_locations = np.isnan(q)
        inf_locations = np.isinf(q)
        q[nan_locations] = 0.0
        q[inf_locations] = 0.0
        return q

    def findScalingCoeffsLimiter(self, PFC, lqCN, lqCF):
        """
        finds scaling coefficient for limiter heat flux profiles
        q0 = qn + qf = q0(fracNC + fracFC)
        where q0 is the peak HF

        fracNC is fraction of power sharing between near exponential in common flux region
        fracFC is fraction of power sharing between far exponential in common flux region
        """
        # Get R and Z vectors at the midplane
#        R_omp_sol = PFC.ep.g['lcfs'][:,0].max()
        R_omp_sol = self.map_R_psi(1.0,PFC)
        R_omp_min = R_omp_sol #this is a limited discharge so Rmin = Rlcfs
        R_omp_max = R_omp_sol + 20.0*(lqCN + lqCF)
        R_omp = np.linspace(R_omp_min, R_omp_max, 1000)
        Z_omp = np.zeros(R_omp.shape)

        # Evaluate B at outboard midplane
        Bp_omp = PFC.ep.BpFunc.ev(R_omp,Z_omp)
        Bt_omp = PFC.ep.BtFunc.ev(R_omp,Z_omp)
        B_omp = np.sqrt(Bp_omp**2 + Bt_omp**2)

        #Find coordinate transformation vector at midplane
        psiaxis = PFC.ep.g['psiAxis']
        psiedge = PFC.ep.g['psiSep']
        deltaPsi = np.abs(psiedge - psiaxis)
        gradPsi = Bp_omp*R_omp
        xfm = gradPsi / deltaPsi

        # transform hf width into flux space
        lqCN_hat = lqCN*xfm
        lqCF_hat = lqCF*xfm

        #Calculate flux at midplane using gfile
        psiN = PFC.ep.psiFunc.ev(R_omp,Z_omp)
        psi = PFC.ep.psiFunc_noN.ev(R_omp,Z_omp)
        PFC.psiMinLCFS = PFC.ep.psiFunc.ev(R_omp_sol,0.0)
        s_hat = psiN - PFC.psiMinLCFS
        print('psiMinLCFS: {:f}'.format(PFC.psiMinLCFS))
        print('un-normalized psiMinLCFS: {:f}'.format(PFC.ep.psiFunc_noN.ev(R_omp_sol,0.0)))
        print('Minimum s_hat: {:f}'.format(s_hat.min()))


        #integral in flux space
        qCN_hat = np.exp(-s_hat / lqCN_hat)
        qCF_hat = np.exp(-s_hat / lqCF_hat)
        intCN = simps(qCN_hat / B_omp, psi)
        intCF = simps(qCF_hat / B_omp, psi)

        q0 = (self.Psol/(2*np.pi)) / (intCN*self.fracCN + intCF*self.fracCF)

        print(q0)

        return q0

    def findScalingCoeffsMultiExp(self, PFC, lqCN, lqCF, lqPN, lqPF):
        """
        finds scaling coefficient for limiter heat flux profiles
        q0 = q0(fracNC + fracFC + fracNP + fracFP)
        where q0 is the peak HF

        fracPN is fraction of power sharing between near exponential in private flux region
        fracPF is fraction of power sharing between far exponential in private flux region
        fracCN is fraction of power sharing between near exponential in common flux region
        fracCF is fraction of power sharing between far exponential in common flux region
        1 = fracPN + fracPF + fracCN + fracCF

        """
        # Get R and Z vectors at the midplane
        R_omp_sol = PFC.ep.g['lcfs'][:,0].max()
        R_omp_min = R_omp_sol - 5.0*(lqPN + lqPF)
        R_omp_max = R_omp_sol + 20.0*(lqCN + lqCF)
        R_omp = np.linspace(R_omp_min, R_omp_max, 1000)
        Z_omp = np.zeros(R_omp.shape)


        # Evaluate B at outboard midplane
        Bp_omp = PFC.ep.BpFunc.ev(R_omp,Z_omp)
        Bt_omp = PFC.ep.BtFunc.ev(R_omp,Z_omp)
        B_omp = np.sqrt(Bp_omp**2 + Bt_omp**2)

        #Find coordinate transformation vector at midplane
        psiaxis = PFC.ep.g['psiAxis']
        psiedge = PFC.ep.g['psiSep']
        deltaPsi = np.abs(psiedge - psiaxis)
        gradPsi = Bp_omp*R_omp
        xfm = gradPsi / deltaPsi
        # transform hf width into flux space
        lqPN_hat = lqPN*xfm
        lqPF_hat = lqPF*xfm
        lqCN_hat = lqCN*xfm
        lqCF_hat = lqCF*xfm

        #Calculate flux at midplane using gfile
        psiN = PFC.ep.psiFunc.ev(R_omp,Z_omp)
        psi = PFC.ep.psiFunc_noN.ev(R_omp,Z_omp)
        PFC.psiMinLCFS = PFC.ep.psiFunc.ev(R_omp_sol,0.0)
        s_hat = psiN - PFC.psiMinLCFS

        #find locations in Private vs Common flux regions
        useP = np.where(s_hat < 0.0)[0]
        useC = np.where(s_hat >= 0.0)[0]

        #integral for SOL in flux space
        if len(useP)>0:
            qPN_hat = np.exp( s_hat[useP] / lqPN_hat[useP])
            qPF_hat = np.exp( s_hat[useP] / lqPF_hat[useP])
            intPN = simps(qPN_hat / B_omp[useP], psi[useP])
            intPF = simps(qPF_hat / B_omp[useP], psi[useP])
        else:
            qPN_hat = 0.0
            qPF_hat = 0.0
            intPN = 0.0
            intPF = 0.0

        if len(useC)>0:
            qCN_hat = np.exp(-s_hat[useC] / lqCN_hat[useC])
            qCF_hat = np.exp(-s_hat[useC] / lqCF_hat[useC])
            intCN = simps(qCN_hat / B_omp[useC], psi[useC])
            intCF = simps(qCF_hat / B_omp[useC], psi[useC])
        else:
            qCN_hat = 0.0
            qCF_hat = 0.0
            intCN = 0.0
            intCF = 0.0

        q0 = (self.Psol/(2*np.pi)) / (intCN*self.fracCN + intCF*self.fracCF +
                                      intPN*self.fracPN + intPF*self.fracPF)

        return q0



    def HFincidentAngle(self,PFC, MHD):
        """
        Calculates b_hat dot n_hat
        dot product between B field and pfc normal
        """
        xyz = PFC.centers
        r,z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        BNorms = MHD.Bfield_pointcloud(PFC.ep, r, z, phi, normal=True)
        PFC.bdotn = np.multiply(PFC.norms, BNorms).sum(1)
        return

    def q_div(self, PFC, MHD, q):
        """
        Calculate divertor heat flux, incorporating flux expansion and
        incident angle.  This takes an already calculated vector, q||, and
        applies it to the divertor tile.
        """
        psi = PFC.psimin
        xyz = PFC.centers

        R_div,Z_div,phi_div = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])

        R_omp = self.map_R_psi(psi,PFC)
        Z_omp = np.zeros(R_omp.shape)
        # Dot product between surface normal and B field creates PFC.bdotn (angle of incidence)
        self.HFincidentAngle(PFC, MHD)
        # Calculate Magnitude of B at Divertor
        Bp_div = MHD.ep.BpFunc.ev(R_div,Z_div)
        Bt_div = MHD.ep.BtFunc.ev(R_div,Z_div)
        B_div = np.sqrt(Bp_div**2 + Bt_div**2)
        # Evaluate B at outboard midplane
        Bp_omp = MHD.ep.BpFunc.ev(R_omp,Z_omp)
        Bt_omp = MHD.ep.BtFunc.ev(R_omp,Z_omp)
        B_omp = np.sqrt(Bp_omp**2 + Bt_omp**2)

#        Bt_omp = MHD.ep.BtFunc.ev(R_omp,Z_omp)
#        BR_omp = MHD.ep.BRFunc.ev(R_omp,Z_omp)
#        BZ_omp = MHD.ep.BZFunc.ev(R_omp,Z_omp)
#        B_omp = np.sqrt(Bt_omp**2 + BR_omp**2 + BZ_omp**2)
#
#        Bt_div = MHD.ep.BtFunc.ev(R_div,Z_div)
#        BR_div = MHD.ep.BRFunc.ev(R_div,Z_div)
#        BZ_div = MHD.ep.BZFunc.ev(R_div,Z_div)
#        B_div = np.sqrt(Bt_div**2 + BR_div**2 + BZ_div**2)


        #For Debugging, plot Bfield Ratio
        #import matplotlib.pyplot as plt
        #testB_div = B_div.reshape(self.grid['Nphi'],self.grid['Nswall']).T
        #testB_omp = B_omp.reshape(self.grid['Nphi'],self.grid['Nswall']).T
        #B_ratio = testB_div/testB_omp
        #CS = plt.contourf(self.grid['phi'], self.grid['Swall'],B_ratio,levels=30,cmap=plt.cm.cool)
        #plt.colorbar(CS, label=r'$B Ratio$')
        #plt.show()
        #Divertor heat flux
        q_div = np.zeros((len(xyz)))
        use = np.where(PFC.shadowed_mask == 0)[0]

        q_div[use] = q[use] * B_div[use]/B_omp * PFC.bdotn[use]

        #for i in range(len(q_div)):
        #	if q_div[i] > 8.0: q_div[i] = 0.0
        #Plot q|| and qdiv
        #import matplotlib.pyplot as plt
        #plt.scatter(self.grid['Swall'][:,0], q_div[0:self.grid['Nswall']], label='qdiv')
        #plt.scatter(self.grid['Swall'][:,0], q[0:self.grid['Nswall']], label='q||')
        #plt.legend()
        #plt.show()
        return np.abs(q_div)

    def write_shadow_pointcloud(self,centers,scalar,dataPath,tag=None):
        print("Creating Shadow Point Cloud")
        log.info("Creating Shadow Point Cloud")

        if tag is None:
            pcfile = dataPath + 'ShadowPointCloud.csv'
        else:
            pcfile = dataPath + 'ShadowPointCloud_'+tag+'.csv'

        #print("Shadow point cloud filename: "+pcfile)
        #log.info("Shadow point cloud filename: "+pcfile)

        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = scalar
        head = "X,Y,Z,ShadowMask"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

        #Now save a vtk file for paraviewweb
        if tag is None:
            tools.createVTKOutput(pcfile, 'points', 'ShadowMask')
        else:
            name = 'ShadowMask_'+tag
            tools.createVTKOutput(pcfile, 'points', name)
        return

    def write_heatflux_pointcloud(self,centers,hf,dataPath,name=None):
        print("Creating Heat Flux Point Cloud")
        log.info("Creating Heat Flux Point Cloud")
        pcfile = dataPath + 'HeatfluxPointCloud.csv'
        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = hf
        head = "X,Y,Z,HeatFlux"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

        #Now save a vtk file for paraviewweb
        tools.createVTKOutput(pcfile, 'points', 'HeatFlux')
        return


    def write_bdotn_pointcloud(self,centers,bdotn,dataPath,name=None):
        print("Creating bdotn Point Cloud")
        log.info("Creating bdotn Point Cloud")
        pcfile = dataPath + 'bdotnPointCloud.csv'
        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = bdotn
        head = "X,Y,Z,bdotn"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

        #Now save a vtk file for paraviewweb
        tools.createVTKOutput(pcfile, 'points', 'bdotn')
        return


    def write_psiN_pointcloud(self,centers,psiN,dataPath,name=None):
        print("Creating Psi Point Cloud")
        log.info("Creating Psi Point Cloud")
        pcfile = dataPath + 'psiPointCloud.csv'
        pc = np.zeros((len(centers), 4))
        pc[:,0] = centers[:,0]*1000.0
        pc[:,1] = centers[:,1]*1000.0
        pc[:,2] = centers[:,2]*1000.0
        pc[:,3] = psiN
        head = "X,Y,Z,psi"
        np.savetxt(pcfile, pc, delimiter=',',fmt='%.10f', header=head)

        #Now save a vtk file for paraviewweb
        tools.createVTKOutput(pcfile, 'points', 'psiN')
        return


    def write_openFOAM_boundary(self, centers, hf, openFoamDir):
        """
        Writes 2 files into <openFoamDir>/constant/boundaryData/
        1) points
        2) 0/T
        These files are then interpolated to the tile surface using the
        openFOAM timeVaryingMappedFixedValue boundary method
        hf should come in in [MW] (I convert to watts here)
        """
        print("Creating Heat Flux Boundary for OpenFoam")
        log.info("Creating Heat Flux Boundary for OpenFoam")

        centers *= 1000.0
        hf *= 1000000.0
        #openFoamDir = '/u/tlooby/OpenFOAM/tlooby-7/run/heatTestLaplace'
        pointFile = openFoamDir + '/constant/boundaryData/STLpatch/points'
        hfFile = openFoamDir + '/constant/boundaryData/STLpatch/0/HF'

        with open(pointFile, 'w') as f:
            f.write('{:d}\n'.format(len(centers[:,0])))
            f.write('(\n')
            for i in range(len(centers[:,0])):
                f.write('({:f} {:f} {:f})\n'.format(centers[i,0], centers[i,1], centers[i,2]))
            f.write(')\n')

        with open(hfFile, 'w') as f:
            f.write('{:d}\n'.format(len(hf)))
            f.write('(\n')
            for i in range(len(hf)):
                f.write('{:f}\n'.format(hf[i]))
            f.write(')\n')

        print("Wrote Heat Flux Boundary for OpenFoam")
        log.info("Wrote Heat Flux Boundary for OpenFoam")



    def PointCloudfromStructOutput(self,file):
        """
        Makes a point cloud from MAFOT structure output initial points
        """
        print("Creating Structure Point Cloud")
        xyz = self.readStructOutput(file)
        pc = np.zeros((int(len(xyz)/2.0),3))
        pc[:,0] = xyz[::2,0]*1000
        pc[:,1] = xyz[::2,1]*1000
        pc[:,2] = xyz[::2,2]*1000
        head = """X,Y,Z"""
        np.savetxt(file, pc, delimiter=',',fmt='%.10f', header=head)



    def findShadows_structure(self,MHD,PFC,targetMeshes, verbose=False, shadowMaskClouds=False):
        """
        Find shadowed faces for a given heatFluxPart object using MAFOT structure.
        Traces field lines from PFC surface, looking for intersections with
        triangles from mesh
        """
        use = np.where(PFC.shadowed_mask == 0)[0]

        print("\nFinding intersections for {:d} faces".format(len(PFC.centers[use])))
        log.info("\nFinding intersections for {:d} faces".format(len(PFC.centers[use])))
        controlfile = '_struct_CTL.dat'
        controlfilePath = MHD.dataPath + '/' + '{:06d}/'.format(PFC.t)
        structOutfile = MHD.dataPath + '/' + '{:06d}/struct.dat'.format(PFC.t)
        gridfile = MHD.dataPath + '/' + '{:06d}/struct_grid.dat'.format(PFC.t)

        numTargetFaces = 0
        targetPoints = []
        targetNorms = []
        print('Number of target parts: {:f}'.format(len(targetMeshes)))
        log.info('Number of target parts: {:f}'.format(len(targetMeshes)))
        for target in targetMeshes:
            numTargetFaces += target.CountFacets
            for face in target.Facets:
                targetPoints.append(face.Points)
                targetNorms.append(face.Normal)
        targetPoints = np.asarray(targetPoints)/1000.0 #scale to m
        targetNorms = np.asarray(targetNorms)

        #for debugging, save a shadowmask at each step up fieldline
        if shadowMaskClouds == True:
            self.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,controlfilePath,tag='original')

        #MAFOT always returns ordered points in CCW (from top) direction,
        #so we need to check which direction we are running so we read the output
        #correctly
        if MHD.MapDirectionStruct == -1:
            print('Tracing with reversed Map Direction')
            log.info('Tracing with reversed Map Direction')
            startIdx = 0
        else:
            print('Tracing with forward Map Direction')
            log.info('Tracing with forward Map Direction')
            startIdx = 1

        #===INTERSECTION TEST 1 (tricky frontface culling / first step up field line)
        dphi = 1.0
        MHD.ittStruct = 1.0
        numSteps = MHD.nTrace #actual trace is (numSteps + 1)*dphi degrees
        #If numSteps = 0, dont do intersection checking
        if numSteps > 0:
            MHD.writeControlFile(controlfile, PFC.t, mode='struct')

            #Perform first integration step
            MHD.writeMAFOTpointfile(PFC.centers[use],gridfile)
            MHD.getMultipleFieldPaths(dphi, gridfile, controlfilePath, controlfile)
            structData = self.readStructOutput(structOutfile)
            os.remove(structOutfile) #clean up

            #this is for printing information about a specific mesh element
#            ptIdx = np.where(use==16886)[0]
            ptIdx = None

            #First we do basic intersection checking
            intersect_mask = self.intersectTestBasic(structData,PFC.norms[use],
                                                    targetPoints,targetNorms,
                                                    MHD,ptIdx)

            PFC.shadowed_mask[use] = intersect_mask

            #for debugging, save a shadowmask at each step up fieldline
            if shadowMaskClouds == True:
                self.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,controlfilePath,tag='test0')

        #===INTERSECTION TEST 2 (multiple steps up field line)
        #Starts at second step up field line
        if numSteps > 1:
            MHD.writeControlFile(controlfile, PFC.t, mode='struct')
            use = np.where(PFC.shadowed_mask == 0)[0]
            intersect_mask2 = np.zeros((len(use)))

            #Perform fist integration step but dont use for finding intersections
            MHD.writeMAFOTpointfile(PFC.centers[use],gridfile)
            MHD.getMultipleFieldPaths(dphi, gridfile, controlfilePath, controlfile)
            structData = self.readStructOutput(structOutfile)
            os.remove(structOutfile) #clean up

            #Perform subsequent integration steps.  Use the point we left off at in
            #last loop iteration as the point we launch from in next loop iteration
            #This amounts to 'walking' up the field line looking for intersections,
            #which is important when field line curvature makes intersections happen
            #farther than 1-2 degrees from PFC surface.
            use2 = np.where(intersect_mask2 == 0)[0]
            for i in range(numSteps):
                print("\nIntersect Trace #2 Step {:d}".format(i))
                log.info("\nIntersect Trace #2 Step {:d}".format(i))
                useOld = use2
                use2 = np.where(intersect_mask2 == 0)[0]
                indexes = np.where([x==useOld for x in use2])[1] #map current steps's index back to intersect_mask2 index

                StartPoints = structData[startIdx::2,:][indexes] #odd indexes are second trace point
                MHD.writeMAFOTpointfile(StartPoints,gridfile)
                MHD.getMultipleFieldPaths(dphi, gridfile, controlfilePath, controlfile)
                structData = self.readStructOutput(structOutfile)
                os.remove(structOutfile) #clean up
                intersect_mask2[use2] = self.intersectTest2(structData,targetPoints,MHD.MapDirection)

                #for debugging, save a shadowmask at each step up fieldline
                if shadowMaskClouds == True:
                    PFC.shadowed_mask[use] = intersect_mask2
                    self.write_shadow_pointcloud(PFC.centers,PFC.shadowed_mask,controlfilePath,tag='test{:d}'.format(i+1))

            #Now revise shadowed_mask taking intersections into account
            PFC.shadowed_mask[use] = intersect_mask2
        return

    def longRangeIntersectCheck(self,PFC,qDiv,targetMeshes,thresh,distPhi,MHD):
        """
        Runs a long range intersection check.

        Because some PFC components may have STL faces that are in locations
        with very high Bt/(Br+Bz) ratio, it might take longer than 10-20 degrees
        for them to intersect with anything.

        This function requires a threshold power.  Any PFC face with power
        higher than the threshold is checked.

        The flagged faces are then checked for intersections over distPhi degrees
        """
        print("\nLong Range Trace Begins")
        log.info("\nLong Range Trace Begins")

        controlfile = '_struct_CTL.dat'
        controlfilePath = MHD.dataPath + '/' + '{:06d}/'.format(PFC.t)
        structOutfile = MHD.dataPath + '/' + '{:06d}/struct.dat'.format(PFC.t)
        gridfile = MHD.dataPath + '/' + '{:06d}/struct_grid.dat'.format(PFC.t)

        numTargetFaces = 0
        targetPoints = []
        print('Number of target parts: {:f}'.format(len(targetMeshes)))
        for target in targetMeshes:
            numTargetFaces += target.CountFacets
            for face in target.Facets:
                targetPoints.append(face.Points)
        targetPoints = np.asarray(targetPoints)/1000.0 #scale to m

        use = np.where(np.abs(qDiv) > thresh)[0]
        found = 0
        print('Number of points above threshold: {:f}'.format(len(use)))
        log.info('Number of points above threshold: {:f}'.format(len(use)))
        for idx in use:
            startPoint = PFC.centers[idx]
#            print("Start Point")
#            print(startPoint)
            MHD.writeMAFOTpointfile(startPoint,gridfile)
            MHD.getMultipleFieldPaths(distPhi, gridfile, controlfilePath, controlfile)
            structData = self.readStructOutput(structOutfile)
            os.remove(structOutfile) #clean up
            #create source format that intersectTest2 reads (odd start points, even end points)
            x = np.insert(structData[:,0], np.arange(len(structData[:,0])), structData[:,0])
            y = np.insert(structData[:,1], np.arange(len(structData[:,1])), structData[:,1])
            z = np.insert(structData[:,2], np.arange(len(structData[:,2])), structData[:,2])
            sources = np.array([x,y,z]).T[1:-1]
            intersects = self.intersectTest2(sources,targetPoints,MHD.MapDirectionStruct)
            if np.sum(intersects) > 0:
#                print("Found long range intersection")
#                log.info("Found long range intersection")
                qDiv[idx] = 0.0
                found +=1
        print("Long range trace found {:f} intersections".format(found))
        log.info("Long range trace found {:f} intersections".format(found))
        return qDiv

    def intersectTestBasic(self,sources,sourceNorms,
                        targets,targetNorms,MHD,ptIdx=None):
        """
        checks if any of the lines (field line traces) generated by MAFOT
        struct program intersect any of the target mesh faces.

        source represents the face where we want to find the heat flux on.
        here it is a 'source' for the field line launching that happens in
        MAFOT.  target represents the potential intersection faces, where a
        field line launched from the source face may intersect.

        Returns a boolean matrix (bitmask) of shape (N) that is true where
        sourceface i intersects with any targetface

        This function first checks for intersection with all of the faces that are
        labeled as 'shadowed' by the backface culling algorithm.  This will
        capture the majority of the intersections, but will miss intricate
        details like when the poloidal field reverses in a gap.

        The function then checks for intersections for the sources that the
        aforementioned check deemed to be heat loaded (no intersection), against
        all of the target faces.  It eliminates self intersections during this
        step by making sure that if an intersection occurs, it is with a face
        that is far away in any direction
        """
        q1 = sources[::2,:]  #even indexes are first trace point
        q2 = sources[1::2,:]  #odd indexes are second trace point
        p1 = targets[:,0,:]  #point 1 of mesh triangle
        p2 = targets[:,1,:]  #point 2 of mesh triangle
        p3 = targets[:,2,:]  #point 3 of mesh triangle

        N = len(q2)
        Nt = len(p1)
        print('{:d} Source Faces and {:d} Target Faces'.format(N,Nt))


        #Cull target front face (PFC front faces will only intersect with
        # other PFC shadowed faces)
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
        tools.targetCtrs = tools.faceCenters(x,y,z)
        targetBackfaceMask = self.backfaceCulling(tools.targetCtrs,targetNorms,MHD)

        #Potential intersection faces
        use = np.where(targetBackfaceMask == 1)[0]
        Nt_use = len(use)

        tools.q1 = q1
        tools.q2 = q2
        tools.p1 = p1[use]
        tools.p2 = p2[use]
        tools.p3 = p3[use]
        tools.Nt = Nt_use

        print('Entering basic intersection test for {:d} potential intersection faces'.format(Nt_use))
        log.info('Entering basic intersection test for {:d} potential intersection faces'.format(Nt_use))
        t0 = time.time()

        #Prepare intersectionTest across multiple cores
        Ncores = multiprocessing.cpu_count() -2 #reserve 2 cores for overhead
        print('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        log.info('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        #each worker receives a single start and end point (q1 and q2),
        #corresponding to one trace from the MAFOT structure output.
        #The worker accesses the many potential intersection triangle vertices
        #through tools class variables (p1,p2,p3).  Making p1,p2,p3 tools class
        #variables eliminates the overhead of transmitting these matrices to
        #each worker, which yields about an order of magnitude speedup.
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
        pool = multiprocessing.Pool(Ncores)
        mask = np.asarray(pool.map(tools.intersectTestParallel, np.arange(N)))
        pool.close()
        print('Found {:f} shadowed faces'.format(np.sum(mask)))
        log.info('Found {:f} shadowed faces'.format(np.sum(mask)))
        print('Time elapsed: {:f}'.format(time.time() - t0))
        log.info('Time elapsed: {:f}'.format(time.time() - t0))



        #Now remove trickier to identify intersections.  Test 2
        #Remove intersections that are outside of 5mm from us.
        #This is run with 'use' representing all faces that were identified
        #as not shadowed (mask = 0) above.
#        use = np.where( mask == 0 )[0]
#        tools.q1 = q1
#        tools.q2 = q2
#        tools.p1 = p1
#        tools.p2 = p2
#        tools.p3 = p3
#        tools.Nt = Nt
#        print('Entering intersection Test #2 for {:d} potential intersection faces'.format(Nt))
#        log.info('Entering intersection Test #2 for {:d} potential intersection faces'.format(Nt))
#        t0 = time.time()
#        print('Spawning tasks to workers')
#        log.info('Spawning tasks to workers')
#        pool = multiprocessing.Pool(Ncores)
#        mask[use] = np.asarray(pool.map(tools.intersectTestParallel_selfCheck, use))
#        pool.close()
#        print('Found {:f} shadowed faces'.format(np.sum(mask[use])))
#        log.info('Found {:f} shadowed faces'.format(np.sum(mask[use])))
#        print('Time elapsed: {:f}'.format(time.time() - t0))
#        log.info('Time elapsed: {:f}'.format(time.time() - t0))

        return mask




    def intersectTest2(self,sources,targets,mapDirection):
        """
        Run an intersection test against all possible source faces
        sources is endpoints of line
        targets is points that make up triangle faces we check for intersections
           against

        This test is called repeatedly from findShadows_structure, as we integrate
        up a field line.  Each time it is called, it checks if the sources line
        intersects a targets face
        """
        q1 = sources[::2,:] #even indexes are first trace point
        q2 = sources[1::2,:] #odd indexes are second trace point
        p1 = targets[:,0,:] #point 1 of mesh triangle
        p2 = targets[:,1,:] #point 2 of mesh triangle
        p3 = targets[:,2,:] #point 3 of mesh triangle

        tools.q1 = q1
        tools.q2 = q2
        tools.p1 = p1
        tools.p2 = p2
        tools.p3 = p3
        N = len(q2)
        Nt = len(p1)
        tools.Nt = Nt

        print('Entering intersection Test #1 for {:d} potential intersection faces'.format(Nt))
        log.info('Entering intersection Test #1 for {:d} potential intersection faces'.format(Nt))
        t0 = time.time()

        #Prepare intersectionTest across multiple cores
        Ncores = multiprocessing.cpu_count() -2 #reserve 2 cores for overhead
        print('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        log.info('Initializing parallel intersection check across {:d} cores'.format(Ncores))
        #each worker receives a single start and end point (q1 and q2),
        #corresponding to one trace from the MAFOT structure output.
        #The worker accesses the many potential intersection triangle vertices
        #through tools class variables (p1,p2,p3).  Making p1,p2,p3 tools class
        #variables eliminates the overhead of transmitting these matrices to
        #each worker, which yields about an order of magnitude speedup.
        print('Spawning tasks to workers')
        log.info('Spawning tasks to workers')
        pool = multiprocessing.Pool(Ncores)
        mask = np.asarray(pool.map(tools.intersectTestParallel, np.arange(N)))
        pool.close()

        print('Found {:f} shadowed faces'.format(np.sum(mask)))
        log.info('Found {:f} shadowed faces'.format(np.sum(mask)))
        print('Time elapsed: {:f}'.format(time.time() - t0))
        log.info('Time elapsed: {:f}'.format(time.time() - t0))

        return mask

    def power_sum_mesh(self, PFC):
        """
        Calculate power by summing over each mesh element.
        Scale to fraction of machine we are analyzing: deltaPhi/2pi
        """
        xyz = PFC.centers
        R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
        deltaPhi = phi.max() - phi.min()
        print('phiMin = {:f}'.format(phi.min()))
        print('phiMax = {:f}'.format(phi.max()))
        log.info('phiMin = {:f}'.format(phi.min()))
        log.info('phiMax = {:f}'.format(phi.max()))
        return np.sum(PFC.qDiv * PFC.areas ) * 2 * np.pi / deltaPhi

    class PFC:
        def __init__(self, MHD, mesh, centers, norms, areas, name=None):
            """
            Heat Flux object for a specific mesh object.
            """
            self.name = name
            self.mesh = mesh
            self.centers = centers / (1000.0)
            self.areas = areas / (1000.0**2)
            self.norms = norms
            self.Nfaces = self.mesh.CountFacets
            self.q = np.zeros((self.Nfaces))
            self.shadowed_mask = heatFlux.backfaceCulling(heatFlux,self.centers,self.norms,MHD)
#            self.shadowed_mask = np.zeros((self.Nfaces))
            #For now just take 1 timestep.  Needs to be adapted later
            try:
                self.t = MHD.timesteps[0]
            except:
                self.t = MHD.timesteps
            self.qBg = 0.0 #Background HF
            self.ep = MHD.ep

            return

#plasma3DClass.py
#Description:   Enable perturbed or M3DC1 3D plasmas in HEAT
#Engineer:      A. Wingen
#Date:          20230227
import sys
import pandas as pd
import numpy as np
import scipy.interpolate as scinter
import os, glob
import shutil
import logging
import subprocess
import toolsClass
tools = toolsClass.tools()

#==========================================================================================================================
#   plasma3D class
#==========================================================================================================================
class plasma3D:
	"""
	This class gets the unshadowed points, writes a points file, launches laminar_mpi,
	reads the output file and uses the penetration depth to calculate a 3D heat flux
	Example call:
	   plasma3D = plasma3DClass.plasma3D(inputFile = 'HEATinput.csv')
	   plasma3D.updatePoints(R, phi, Z)
	   plasma3D.launchLaminar(4, 'testrun')
	"""

	def __init__(self):
		self.R = None		# in m
		self.phi = None		# right-handed angle in degrees
		self.Z = None		# in m
		self.psimin = None
		self.Lc = None		# in km
		
		# Boundary Box limits
		self.bbRmin = None	
		self.bbRmax = None
		self.bbZmin = None
		self.bbZmax = None
		
		self.allowed_vars = ['plasma3Dmask','shot','time','tmax','gFile','itt','response','selectField','useIcoil','sigma','charge','Ekin','Lambda','Mass']
	
	
	def initializePlasma3D(self, shot, time, gFile = None, inputFile = None, cwd = None, inputDir = None):
		"""
		Set up basic input vars
		gfile should include the full path and file name
		inputFile is the main .csv file with input variables
		cwd is the HEAT data folder for this shot and pfc, typically ~/HEAT/data/<machine>_<shot>_<tag>/<time>/<pfcName>
		inputDir is the folder with input files, typically /root/terminal/<machine>
		"""
		self.shot = tools.makeInt(shot)
		self.time = tools.makeInt(time)
		if cwd is None: self.cwd = os.getcwd()
		else: self.cwd = cwd
		if inputDir is None: self.inputDir = self.cwd
		else: self.inputDir = inputDir
		if gFile is None: gFile = self.cwd + '/g' + format(int(self.shot),'06d') + '.' + format(int(self.time),'05d')
		self.gFile = gFile
		if inputFile is not None: self.read_input_file(inputFile)
		else: self.setMAFOTctl()	# just defaults
		self.readM3DC1supFile()


	def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
		"""
		sets up pythonic string number formats for shot and timesteps
		"""
		self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
		self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
		return
		
		
	def print_settings(self):
		"""
		Print all inputs
		"""
		print('#=============================================================')
		print('#                Equilibrium Variables')
		print('#=============================================================')
		print('shot = ' + str(self.shot))
		print('time = ' + str(self.time))
		print('gFile = ' + str(self.gFile))
		print('cwd = ' + str(self.cwd))
		print('#=============================================================')
		print('#                3D Plasma Variables')
		print('#=============================================================')
		print('plasma3Dmask = ' + str(self.plasma3Dmask))
		print('itt = ' + str(self.itt))
		print('useIcoil = ' + str(self.useIcoil))
		print('sigma = ' + str(self.sigma))
		print('charge = ' + str(self.charge))
		print('Ekin = ' + str(self.Ekin))
		print('Lambda = ' + str(self.Lambda))
		print('Mass = ' + str(self.Mass))
		print('#=============================================================')
		print('#                Boundary Box Variables')
		print('#=============================================================')
		print('Rmin = ' + str(self.bbRmin))
		print('Rmax = ' + str(self.bbRmax))
		print('Zmin = ' + str(self.bbZmin))
		print('Zmax = ' + str(self.bbZmax))
		print('#=============================================================')
		print('#                M3D-C1 Variables')
		print('#=============================================================')
		print('response = ' + str(self.response))
		print('selectField = ' + str(self.selectField))
		for i in range(len(self.C1Files)):
			print('File ' + str(i+1) + ' = ' + self.C1Files[i])
			print('   Scale = ' + str(self.C1scales[i]))
			print('   Phase = ' + str(self.C1phases[i]))
		

	def setMAFOTctl(self, itt = 300, response = 0, selectField = -1, useIcoil = 0, 
				sigma = 0, charge = -1, Ekin = 10, Lambda = 0.1, Mass = 2):
		"""
		Set the MAFOT specific class variables
		"""
		self.itt = tools.makeInt(itt) 						# toroidal iterations
		self.response = tools.makeInt(response) 			# M3D-C1 Plasma Response (0=no,>1=yes)
		self.selectField = tools.makeInt(selectField) 		# MHD fields to use (-3=VMEC,-2=SIESTA,-1=gfile,M3DC1:0=Eq,1=I-coil,2=both)
		self.useIcoil = tools.makeInt(useIcoil) 			# 0=no, 1=yes
		self.sigma = tools.makeInt(sigma) 					# Particle Direction (1=co-pass,-1=ctr-pass,0=field-lines)
		self.charge = tools.makeInt(charge) 				# Partile Charge (-1=electrons,>=1=ions)
		self.Ekin = tools.makeFloat(Ekin) 					# Particle Kinetic Energy in keV
		self.Lambda = tools.makeFloat(Lambda) 				# Ratio of perpendicular to parallel velocity
		self.Mass = tools.makeInt(Mass) 					# Particle Ion Mass (H=1, D=2, He=4)


	def setBoundaryBox(self, MHD, CAD):
		self.bbRmin = np.min([CAD.Rmin, MHD.ep[0].g['wall'][:,0].min()])
		self.bbRmax = np.max([CAD.Rmax, MHD.ep[0].g['wall'][:,0].max()])
		self.bbZmin = np.min([CAD.Zmin, MHD.ep[0].g['wall'][:,1].min()])
		self.bbZmax = np.max([CAD.Zmax, MHD.ep[0].g['wall'][:,1].max()])


	def setM3DC1input(self, C1Files = ['./C1.h5'], scales = [1], phases = None):
		"""
		Set the M3D-C1 specific class variables
		"""
		self.C1Files = C1Files 
		self.C1scales = scales
		if phases is None: self.C1phases = np.zeros(len(C1Files))
		else: self.C1phases = phases
	
	
	def readM3DC1supFile(self):
		"""
		Read M3D-C1 supplemental input file, if it already exists
		"""
		C1Files = []
		scales = []
		phases = []

		if not os.path.isfile(self.inputDir + '/' + 'm3dc1sup.in'): 
			print('m3dc1sup.in file not found!')
			self.setM3DC1input()
			return
		
		with open(self.inputDir + '/' + 'm3dc1sup.in') as f:
			lines = f.readlines()
		
		for line in lines:
			line = line.strip()
			if len(line) < 1: continue
			if line[0] == '#': continue
			words = line.split()
			c1file = words[0]
			if ('./' in c1file): c1file = c1file.replace('./', self.inputDir + '/')
			C1Files.append(c1file)
			scales.append(tools.makeFloat(words[1]))
			if len(words) > 2: phases.append(tools.makeFloat(words[2]))
			else: phases.append(0)
		
		if len(C1Files) < 1: 
			print('Error reading m3dc1sup.in')
			self.setM3DC1input()
			return
		else:
			self.setM3DC1input(C1Files, scales, phases)
			print('M3D-C1: ' + self.inputDir + '/' + 'm3dc1sup.in read successfully')
			return
		
	
	def read_input_file(self, file):
		"""
		Reads the 3D plasma csv input file
		Format for input file is comma delimited, # are comments.
		Example:
		#Important Comment
		variable_name, value
		"""
		if os.path.isfile(file): 
			tools.read_input_file(self, file)
			self.setTypes()
			print('Input file: ' + file + ' read successfully')
		else: 
			print('Input file: ' + file + ' not found!')
			self.setMAFOTctl()	# just defaults


	def setTypes(self):
		"""
		Set variable types for the stuff that isnt a string from the input file
		"""
		integers = ['plasma3Dmask','shot','time','tmax','itt','response','selectField','useIcoil','sigma','charge','Mass']
		floats = ['Ekin','Lambda']
		setAllTypes(self, integers, floats)     # this is not a typo, but the correct syntax for this call


	def updatePointsFromCenters(self, xyz):
		"""
		Converts xyz of centers into R,phi,Z, update class variables and write the points file
		"""
		if len(xyz.shape) > 1:
			R,Z,phi = tools.xyz2cyl(xyz[:,0],xyz[:,1],xyz[:,2])
		else:
			R,Z,phi = tools.xyz2cyl(xyz[0],xyz[1],xyz[2])

		phi = np.degrees(phi)
		self.updatePoints(R, phi, Z)
		
	
	def updatePoints(self, R, phi, Z):
		"""
		Get arrays of R,phi,Z, update class variables and write the points file
		"""
		self.R = R
		self.phi = phi 
		self.Z = Z
		self.writePoints()

	
	def writePoints(self, filename = 'points3DHF.dat'):
		"""
		Write the points file in CWD from the class variables
		"""
		R = self.R.flatten()
		phi = self.phi.flatten()
		Z = self.Z.flatten()
		N = len(R)
		
		with open(self.cwd + '/' + filename,'w') as f:
			f.write("# Number of points = " + str(N) + "\n")
			for i in range(N): 
				f.write(str(R[i]) + "\t" + str(phi[i]) + "\t" + str(Z[i]) + "\n")
				
				
	def launchLaminar(self, nproc, tag = None, MapDirection = 0):
		"""
		Write all inout files and launch MAFOT
		Read the output file when finished
		"""
		if tag is None: tag = ''
		self.writeControlFile(MapDirection)
		self.writeM3DC1supFile()
		self.writeCoilsupFile()
		
		if nproc > 20: nproc = 20
		self.nproc = nproc
		self.tag = tag
		print('-'*80)
		print('Launching 3D plasma field line tracing')
		
		bbLimits = str(self.bbRmin) + ',' + str(self.bbRmax) + ',' + str(self.bbZmin) + ',' + str(self.bbZmax)
		args = ['mpirun','-n',str(nproc),'heatlaminar_mpi','-P','points3DHF.dat','-B',bbLimits,'_lamCTL.dat',tag]
		current_env = os.environ.copy()        #Copy the current environment (important when in appImage mode)
		subprocess.run(args, env=current_env, cwd=self.cwd)
		#print('mpirun -n ' + str(nproc) + ' heatlaminar_mpi' + ' -P points3DHF.dat' + ' _lamCTL.dat' + ' ' + tag)
		
		#self.wait2finish(nproc, tag)
		self.readLaminar(tag)
		print('3D plasma field line tracing complete')
		print('-'*80)
		
	
	def readLaminar(self, tag = None):
		"""
		Read the MAFOT outputfile and set psimin and Lc class variables
		"""
		if tag is None: tag = ''    # this tag has len(tag) = 0
		
		file = self.cwd + '/' + 'lam_' + tag + '.dat'
		if os.path.isfile(file): 
			lamdata = np.genfromtxt(file,comments='#')
			self.Lc = lamdata[:,3]
			self.psimin = lamdata[:,4]
		else:
			print('MAFOT output file: ' + file + ' not found!')
			
			
	def cleanUp(self, tag = None):
		logs = self.cwd + '/' + 'log*'
		for f in glob.glob(logs): os.remove(f)


	def writeControlFile(self, MapDirection):
		"""
		Write MAFOT control file
		"""
		with open(self.cwd + '/' + '_lamCTL.dat', 'w') as f:
			f.write('# Parameterfile for HEAT Programs\n')
			f.write('# Shot: ' + format(int(self.shot),'06d') + '\tTime: ' + format(int(self.time),'04d') + 'ms\n')
			f.write('# Path: ' + self.gFile + '\n')
			f.write('NZ=\t10\n')
			f.write('itt=\t' + str(self.itt) + '\n')
			f.write('Rmin=\t1\n')
			f.write('Rmax=\t2\n')
			f.write('Zmin=\t-1\n')
			f.write('Zmax=\t1\n')
			f.write('NR=\t10\n')
			f.write('phistart(deg)=\t0\n')
			f.write('MapDirection=\t' + str(MapDirection) + '\n')
			f.write('PlasmaResponse(0=no,>1=yes)=\t' + str(self.response) + '\n')
			f.write('Field(-3=VMEC,-2=SIESTA,-1=gfile,M3DC1:0=Eq,1=I-coil,2=both)=\t' + str(self.selectField) + '\n')
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')
			f.write('createPoints(0=setR,3=setpsi)=\t0\n')    # This must be entry index 12
			f.write('unused=\t0\n')
			f.write('unused=\t0\n')
			f.write('unused=\t0\n')
			f.write('ParticleDirection(1=co-pass,-1=ctr-pass,0=field-lines)=\t' + str(self.sigma) + '\n')   # This must be entry index 16
			f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge) + '\n')
			f.write('Ekin[keV]=\t' + str(self.Ekin) + '\n')
			f.write('lambda=\t' + str(self.Lambda) + '\n')
			f.write('Mass=\t' + str(self.Mass) + '\n')    # This must be entry index 20
			f.write('unused=\t0\n')
			f.write('unused=\t0\n')
			f.write('dpinit=\t1.0\n')   # This must be entry index 23
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')


	def writeM3DC1supFile(self):
		"""
		Write M3D-C1 supplemental input file
		Overwrites any existing one.
		"""
		with open(self.cwd + '/' + 'm3dc1sup.in', 'w') as f:
			for i in range(len(self.C1Files)):
				f.write(self.C1Files[i] + '\t' + str(self.C1scales[i]) + '\t' + str(self.C1phases[i]) + '\n')


	def writeCoilsupFile(self, machine = None):
		"""
		This would be machine specific and needs updating in the future
		"""
		# with open(self.cwd + '/' + 'heatsup.in', 'w') as f:
		#	pass

		return
		
	
	def wait2finish(self, nproc, tag):
		import time
		print ('Waiting for job to finish...', end='')
		time.sleep(5)	# wait 5 seconds
		while(self.isProcessRunning()):
			time.sleep(60)		# wait 1 minute; no CPU usage
		print('done')
			
		if not self.isComplete():
			print('MAFOT run ended prematurely. Attempt restart...')
			subprocess.call(['mpirun','-n',str(nproc),'heatlaminar_mpi','-P','points.dat','_lamCTL.dat',tag])
			self.wait2finish(nproc, tag)
		else: return
	
	
	def isComplete(self, logsPath = None):
		if logsPath is None: 
			logsPath = self.cwd
		if not logsPath[-1] == '/': logsPath += '/'

		allFiles = os.listdir(logsPath)
		fileList = [f for f in allFiles if '_Master.dat' in f]
		for file in fileList:
			lines = subprocess.check_output(['tail', file])
			lines = lines.decode('UTF-8').split('\n')
			for line in lines: 
				if 'Program terminates normally' in line: 
					return True
			else: 
				return False
	
	
	def isProcessRunning(self):
		import getpass
		user = getpass.getuser()
		lines = subprocess.check_output('ps aux | grep heatlaminar_mpi', shell = True)	# for some reason only works with shell = True
		lines = lines.decode('UTF-8').split('\n')
		for line in lines:
			if ('mpirun' in line) & (user in line):
				self.pid = int(line.strip().split()[1])
				return True
		self.pid = -1
		return False
		

#==========================================================================================================================
#   heatflux3D class
#==========================================================================================================================
class heatflux3D:
	"""
	This class gets the penetration depth and connection length of unshadowed points
	and calculates the parallel heat flux
	still needs normalization to power balance
	still needs incident angle
	Example call:
	"""

	def __init__(self):
		self.psimin = None
		self.Lc = None		# in km
		self.N = 1
		self.q = np.zeros(self.N)
		self.ep = None	# equilParams_class instance for EFIT equilibrium
		self.HFS = None	# True: use high field side SOL, False: use low field side SOL
		self.allowed_vars = ['Lcmin', 'lcfs', 'lqEich', 'S', 'Psol', 'qBG']


	def intializeHF3D(self, Lc, psimin, ep, inputFile = None, T = None, HFS = False):
		"""
		Set up basic input vars
		"""
		self.psimin = psimin
		self.Lc = Lc		# in km
		self.N = len(self.Lc)
		self.q = np.zeros(self.N)
		self.ep = ep	# equilParams_class instance for EFIT equilibrium
		self.HFS = HFS	# True: use high field side SOL, False: use low field side SOL

		if inputFile is not None: self.read_input_file(inputFile)
		else: self.setHFctl()	# just defaults
		
		self.good = self.isGoodPoint()
		self.pfr = self.isPFR()

		# set T profile
		if T is None: T = 2							# temperature at top of pedestal in keV
		if isinstance(T, str):						# file name for T profile data
			TData = np.loadtxt(T)
			psiT = TData[:,0]
			T = TData[:,1]
			self.fT = scinter.UnivariateSpline(psiT, T, s = 0)
		elif isinstance(T, np.ndarray):				# array of T data assuming psi = [0, 1.1]
			psiT = np.linspace(0, 1.1, len(T))
			self.fT = scinter.UnivariateSpline(psiT, T, s = 0)
		else:										# any other option
			try:
				T = float(T)						# temperature at top of pedestal in keV
				self.fT = lambda x: Tprofile(x, T)	# generic temperature profile
			except:
				raise RuntimeError('Invalid T profile data')


	def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
		"""
		sets up pythonic string number formats for shot and timesteps
		"""
		self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
		self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
		return

	
	def print_settings(self):
		"""
		Print all inputs
		"""
		print('#=============================================================')
		print('#                Optical HF Variables')
		print('#=============================================================')
		print('lqEich = ' + str(self.lqEich))
		print('S = ' + str(self.S))
		print('Psol = ' + str(self.Psol))
		print('qBG = ' + str(self.qBG))
		print('#=============================================================')
		print('#                3D Plasma Variables')
		print('#=============================================================')
		print('Lcmin = ' + str(self.Lcmin))
		print('lcfs = ' + str(self.lcfs))
		

	def setHFctl(self, Lcmin = 0.075, lcfs = 0.97, lqEich = 5, S = 2, Psol = 10, qBG = 0):
		"""
		Set the specific class variables
		"""
		self.Lcmin = tools.makeFloat(Lcmin) 		# minimum connection length in SOL to separateout the PFR, in km
		self.lcfs = tools.makeFloat(lcfs) 			# psi of the Last Closed Flux Surface inside the stochastic layer
		self.lqEich = tools.makeFloat(lqEich) 		# heat flux layer width for Eich profile, in mm
		self.S = tools.makeFloat(S) 				# heat flux layer extension width in PFR, in mm
		self.Psol = tools.makeFloat(Psol) 			# total power into SOL, in MW
		self.qBG = tools.makeFloat(qBG) 			# background heat flux in MW/m^2
	
	
	def read_input_file(self, file):
		"""
		Reads the 3D plasma csv input file
		Format for input file is comma delimited, # are comments.
		Example:
		#Important Comment
		variable_name, value
		"""
		if os.path.isfile(file): 
			tools.read_input_file(self, file)
			self.time = self.tmax
			self.setTypes()
			print('Input file: ' + file + ' read successfully')
		else: 
			print('Input file: ' + file + ' not found!')
			self.setHFctl()	# just defaults


	def setTypes(self):
		"""
		Set variable types for the stuff that isnt a string from the input file
		"""
		integers = []
		floats = ['Lcmin', 'lcfs', 'lqEich', 'S', 'Psol', 'qBG']
		setAllTypes(self, integers, floats)
		
		
	def isGoodPoint(self):
		"""
		Returns a boolean mask for points to use or not
		True: good point
		False: point failed during laminar run
		"""
		mask = np.ones(self.N, dtype = bool)
		fails = np.where(self.psimin == 10)
		mask[fails] = False
		return mask
		
		
	def isPFR(self):
		"""
		Returns a boolean mask for point in PFR or not
		True: point in PFR
		False: point in SOL or lobes
		"""
		mask = np.zeros(self.N, dtype = bool)
		pfr = np.where((self.psimin < 1) & (self.Lc < Lcmin))
		mask[pfr] = True
		return mask


	def getq_conduct(self, kappa = 2000, T0 = 0, L = 1, limit = True, Tpfr = 0.01):
		"""
		Input:
		  kappa = electron heat conductivity in W/m/eV^3.5
		  T0 = electron temperature at sheath entrance near target in keV
		  L = conduction distance between target and LCFS in km
		Output:
		  updates self.q		
		"""
		T = self.fT(psimin)			# this is now temperature in keV
		
		if limit: 
			T[self.psimin < self.lcfs] = self.fT(self.lcfs)
		T[self.pfr] = Tpfr	# 10 eV in private flux region
		
		self.q = 2.0/7.0 * kappa/L * (T**3.5 - T0**3.5) * (1e+3)**3.5/1e+6   # in MW/m^2


	def getq_eich(self):
		"""
		Computes heat flux based on an Eich profile with lobes
		updates self.q		
		"""
		self.q, q0 = self.set_layer(self.psimin, self.lqEich, self.S, lcfs = self.lcfs, lobes = True)
		self.q[self.pfr],_ = self.set_layer(self.psimin[self.pfr], self.lqEich, self.S, q0 = q0)

	
	def set_layer(self, psi, lq, S, lcfs = 1.0, q0 = 1, lobes = False):
		"""
		psi is flat array of normalized flux
		lq is heat flux width at midplane in mm
		S is the private flux region spreading in mm
		returns flat array of heat flux based on Eich profile
		"""
		x = self.map_R_psi(psi)
		xsep = self.map_R_psi(1.0)

		if lobes:
			s0 = self.map_R_psi(lcfs)
			s = self.map_R_psi(np.linspace(lcfs-0.05,lcfs+0.1,10000))
		else:
			s0 = self.map_R_psi(1.0)
			s = self.map_R_psi(np.linspace(0.95,1.1,10000))
			
		qref = eich_profile(s, lq, S, s0, q0 = 1, qBG = self.qBG, fx = 1)
		idx = qref.argmax()
		smax = s[idx]
		qmax = qref[idx]
	
		if self.HFS:
			x *= -1
			s0 *= -1
			smax *= -1
			xsep *= -1
			x0 = smax
		else:
			x0 = s0 - (smax-s0)	# now the peak amplitude is at psi = lcfs; qlcfs = qmax too
					
		q = eich_profile(x, lq, S, x0, q0 = 1, qBG = self.qBG, fx = 1)
		qsep = eich_profile(xsep, lq, S, x0, q0 = 1, qBG = self.qBG, fx = 1)
		
		if lobes:
			q[psi < psilcfs] = qmax
		
		return q/q.max()*q0, qsep/q.max()*q0


	def map_R_psi(self, psi):
		"""
		Map normalized poloidal flux psi to R at midplane (Z = 0)
		psi is flat array
		return R(psi)
		"""
		if self.HFS:
			R = np.linspace(self.ep.g['RmAxis'], self.ep.g['R1'], 100)
		else:
			R = np.linspace(self.ep.g['RmAxis'], self.ep.g['R1'] + self.ep.g['Xdim'], 100)
			
		Z = self.ep.g['ZmAxis']*np.ones(len(R))
		p = self.ep.psiFunc.ev(R,Z)
		
		f = scinter.UnivariateSpline(p, R, s = 0, ext = 'const')	# psi outside of spline domain return the boundary value
		return f(psi)
	


#==========================================================================================================================
#   general functions
#==========================================================================================================================

def Tprofile(psi, Tped, deriv = False):
	"""
	A default T profile, used only if no T profile data is provided
	Input:
	  psi = normalized poloidal flux
	  Tped = temperature at top of pedestal in keV
	  deriv = bool, return derivative of profile, default is false
	Return:
	  T (,dT) = T(psi) profile (, derivative of profile)
	"""
	xs = 0.975		# Symmetry point in Pedestal
	dw = 0.04		# half of Pedestal width

	f = lambda x: 0.5*np.tanh(2*(xs - x)/dw) + 2*np.exp(-x*2)
	T0 = Tped/(f(xs-dw) - f(1.2))
	T = T0*(f(psi) - f(1.2))
	if deriv:
		dT = -T0/dw*(1 - np.tanh(2*(xs - psi)/dw)**2) - T0*4*np.exp(-psi*2)
		return T, dT
	else: return T


def setAllTypes(obj, integers, floats):
	"""
	Set data types for variales in obj
	"""
	for var in integers:
		if (getattr(obj, var) is not None) and (~np.isnan(float(getattr(obj, var)))):
			try:
				setattr(obj, var, tools.makeInt(getattr(obj, var)))
			except:
				print("Error with input file var "+var+".  Perhaps you have invalid input values?")
				#log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")
	for var in floats:
		if var is not None:
			if (getattr(obj, var) is not None) and (~np.isnan(float(getattr(obj, var)))):
				try:
					setattr(obj, var, tools.makeFloat(getattr(obj, var)))
				except:
					print("Error with input file var "+var+".  Perhaps you have invalid input values?")
					#log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")


def eich_profile(x, lq, S, s0, q0, qBG = 0, fx = 1):
	"""
	Based on the paper: T.Eich et al.,PRL 107, 215001 (2011)
	lq is heat flux width at midplane in mm
	S is the private flux region spreading in mm
	s0 is the separatrix location at Z = 0 in m
	q0 is the amplitude
	qBG is the background heat flux
	fx is the flux expansion between outer midplane and target plate
	x is in m
	return function q(x)
	
	in Eich paper: s (here x) and s0 are distances along target, mapped from midplane using flux expansion,
	so: s = s_midplane * fx; same for s0, with s0 the position of strikeline on target
	Here, use s_midplane directly, so set fx = 1 and identify s = s_midplane = R and s0 = Rsep
	"""
	from scipy.special import erfc
	lq *= 1e-3		# in m now
	S *= 1e-3		# in m now
	a = lq*fx
	b = 0.5*S/lq
	c = S*fx
	q = 0.5 * q0 * np.exp(b**2 - (x-s0)/a) * erfc(b - (x-s0)/c) + qBG
	return q



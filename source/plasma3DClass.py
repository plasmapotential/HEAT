#plasma3DClass.py
#Description:   Enable perturbed or M3DC1 3D plasmas in HEAT
#Engineer:      A. Wingen
#Date:          20230227
import sys
import pandas as pd
import numpy as np
import scipy.interpolate as scinter
import scipy.integrate as integ
import os, glob
import shutil
import logging
import subprocess
import toolsClass
try:
    from subprocess import DEVNULL  # Python 3.
except ImportError:
    DEVNULL = open(os.devnull, 'wb')
    
tools = toolsClass.tools()
log = logging.getLogger(__name__)

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
		self.useVertices = False
		
		# Boundary Box limits
		self.bbRmin = None	
		self.bbRmax = None
		self.bbZmin = None
		self.bbZmax = None
		
		# Default inputs
		self.plasma3Dmask = False
		self.NCPUs = 10
		self.loadHF = False
		self.loadBasePath = None
		
		
	def allowed_class_vars(self):
		"""
		.. Writes a list of recognized class variables to HEAT object
		.. Used for error checking input files and for initialization
		.. These variables are read in from the input file. 
		.. The call is in engine_class.loadInputs

		plasma3D Variables:
		----------------------------

		:plasma3Dmask:  boolean that determines if 3D plasmas should be turned on
		  If True, uses M3DC1 equilibrium
		:itt: integer number of toroidal iterations for the MAFOT 'heatlaminar_mpi' run
		:response: integer to select the M3D-C1 time_xxx.h5 file to use. 
		  For linear M3D-C1, time_001.h5 includes the plasma response and 
		  time_000.h5 is the vacuum field only
		:selectField: integer to switch between g-file (-1) or M3D-C1 (2) in MAFOT run
		:useIcoil: integer to turn on vacuum field perturbation coils, default is 0
		:sigma: integer to switch between tracing field lines (0), co-passing (+1) or 
		  counter-passing (-1) particles in MAFOT
		:charge: integer, charge of particles, -1 for electrons, >0 for ions
		  ignored for sigma = 0
		:Ekin: float, inetic energy of particles in keV; ignored for sigma = 0
		:Lambda: float, ratio of perpendicular to parallel particle velocity
		  ignored for sigma = 0
		:Mass: integer, mass number of particles, 1 for electrons and H, 2 for D, 4 for He
		  ignored for sigma = 0
		:loadHF: boolean, True means load previous heat flux results instead of 
		  running MAFOT, False means run MAFOT
		:loadBasePath: string, Path for find previous results if loadHF is True
		:NCPUs: integer, number of CPUs to use in MAFOT, default is 10
		"""
		self.allowed_vars = ['plasma3Dmask','itt','response',
				'selectField','useIcoil','sigma','charge','Ekin','Lambda','Mass','loadHF',
				'loadBasePath','NCPUs']
	
	
	def setTypes(self):
		"""
		Set variable types for the stuff that isnt a string from the input file
		"""
		integers = ['itt','response','selectField','useIcoil','sigma','charge','Mass','NCPUs']
		floats = ['Ekin','Lambda']
		bools = ['plasma3Dmask','loadHF']
		setAllTypes(self, integers, floats, bools)     # this is not a typo, but the correct syntax for this call


	def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
		"""
		sets up pythonic string number formats for shot and timesteps
		"""
		self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
		self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
		
		
	def initializePlasma3D(self, shot, time, gFile = None, inputDir = None):
		"""
		Set up basic input vars
		gfile should include the full path and file name
		inputFile is the main .csv file with input variables
		cwd is the HEAT data folder on the host machine for this shot and pfc, 
		  typically ~/HEAT/data/<machine>_<shot>_<tag>/<time>/<pfcName>
		inputDir is the folder in the docker container with input files, 
		  typically $HEAT_HOME/terminal/<machine>
		"""
		self.shot = tools.makeInt(shot)
		self.time = tools.makeInt(time)
		if inputDir is None: inputDir = os.getcwd()
		self.inputDir = inputDir
		if gFile is None: gFile = self.inputDir + '/g' + format(int(self.shot),'06d') + '.' + format(int(self.time),'05d')
		self.gFile = gFile		# this is not used for ep, but just as a string in the MAFOT control file
		self.readM3DC1supFile()


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

		log.info('#=============================================================')
		log.info('#                Equilibrium Variables')
		log.info('#=============================================================')
		log.info('shot = ' + str(self.shot))
		log.info('time = ' + str(self.time))
		log.info('gFile = ' + str(self.gFile))
		log.info('#=============================================================')
		log.info('#                3D Plasma Variables')
		log.info('#=============================================================')
		log.info('plasma3Dmask = ' + str(self.plasma3Dmask))
		log.info('itt = ' + str(self.itt))
		log.info('useIcoil = ' + str(self.useIcoil))
		log.info('sigma = ' + str(self.sigma))
		log.info('charge = ' + str(self.charge))
		log.info('Ekin = ' + str(self.Ekin))
		log.info('Lambda = ' + str(self.Lambda))
		log.info('Mass = ' + str(self.Mass))
		log.info('#=============================================================')
		log.info('#                Boundary Box Variables')
		log.info('#=============================================================')
		log.info('Rmin = ' + str(self.bbRmin))
		log.info('Rmax = ' + str(self.bbRmax))
		log.info('Zmin = ' + str(self.bbZmin))
		log.info('Zmax = ' + str(self.bbZmax))
		log.info('#=============================================================')
		log.info('#                M3D-C1 Variables')
		log.info('#=============================================================')
		log.info('response = ' + str(self.response))
		log.info('selectField = ' + str(self.selectField))
		for i in range(len(self.C1Files)):
			log.info('File ' + str(i+1) + ' = ' + self.C1Files[i])
			log.info('   Scale = ' + str(self.C1scales[i]))
			log.info('   Phase = ' + str(self.C1phases[i]))


	def updatePFCdata(self, cwd):
		"""
		Update class variables that are specific for each PFC
		"""
		self.cwd = cwd
		print('Plasma3D current working directory set to = ' + str(self.cwd))
		log.info('Plasma3D current working directory set to = ' + str(self.cwd))


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
			log.info('m3dc1sup.in file not found!')
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
			log.info('Error reading m3dc1sup.in')
			self.setM3DC1input()
			return
		else:
			self.setM3DC1input(C1Files, scales, phases)
			print('M3D-C1: ' + self.inputDir + '/' + 'm3dc1sup.in read successfully')
			log.info('M3D-C1: ' + self.inputDir + '/' + 'm3dc1sup.in read successfully')
			return
		

	def updatePointsFromVertices(self, xvertices, yvertices, zvertices, centers):
		"""
		Converts xyz of vertices and centers into R,phi,Z, update class variables and write the points file
		"""
		self.useVertices = True
		N = len(centers[:,0])
		R,Z,phi = np.zeros((N,4)),np.zeros((N,4)),np.zeros((N,4))
		R[:,0],Z[:,0],phi[:,0] = tools.xyz2cyl(centers[:,0],centers[:,1],centers[:,2])
		R[:,1],Z[:,1],phi[:,1] = tools.xyz2cyl(xvertices[:,0],yvertices[:,0],zvertices[:,0])
		R[:,2],Z[:,2],phi[:,2] = tools.xyz2cyl(xvertices[:,1],yvertices[:,1],zvertices[:,1])
		R[:,3],Z[:,3],phi[:,3] = tools.xyz2cyl(xvertices[:,2],yvertices[:,2],zvertices[:,2])
		R,Z,phi = R.flatten(),Z.flatten(),phi.flatten()
		phi = np.degrees(phi)
		self.updatePoints(R, phi, Z)


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
				
				
	def launchLaminar(self, NCPUs = None, tag = None, MapDirection = 0, verbose=False):
		"""
		Write all input files and launch MAFOT
		Read the output file when finished
		"""
		if tag is None: tag = ''
		self.writeControlFile(MapDirection)
		self.writeM3DC1supFile()
		self.writeCoilsupFile()
		
		self.tag = tag
		print('Launching 3D plasma field line tracing on ' + str(self.NCPUs) + ' cores')
		log.info('Launching 3D plasma field line tracing ' + str(self.NCPUs) + ' cores')
		
		bbLimits = str(self.bbRmin) + ',' + str(self.bbRmax) + ',' + str(self.bbZmin) + ',' + str(self.bbZmax)
		args = ['mpirun','-n',str(self.NCPUs),'heatlaminar_mpi','-P','points3DHF.dat','-B',bbLimits,'_lamCTL.dat',tag]
		current_env = os.environ.copy()        #Copy the current environment (important when in appImage mode)
		
		try:
			if verbose:
				proc = subprocess.Popen(args, env=current_env, cwd=self.cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,)
				for line in proc.stdout:
					sys.stdout.write(line)
				returncode = proc.wait()
				if returncode != 0:
					raise subprocess.CalledProcessError(returncode, args)
			else:
				result = subprocess.run(args, env=current_env, cwd=self.cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False,)
				if result.returncode != 0:
					log.error("Laminar run failed with code %s", result.returncode)
					log.error("Laminar output:\n%s", result.stdout)
					print.error("Laminar run failed with code %s", result.returncode)
					print.error("Laminar output:\n%s", result.stdout)
					raise subprocess.CalledProcessError(
                    	result.returncode, args, output=result.stdout
                		)
		except Exception as e:
			print("Error: laminar run failed")
			print(e.output)
			log.info("Error: laminar run failed")
			log.info(e.output)

		
		#if verbose == False:
		#	subprocess.run(args, env=current_env, cwd=self.cwd, stderr=DEVNULL)
		#else:
		#	subprocess.run(args, env=current_env, cwd=self.cwd) #dont suppress error messages
		#print('mpirun -n ' + str(self.NCPUs) + ' heatlaminar_mpi' + ' -P points3DHF.dat' + ' _lamCTL.dat' + ' ' + tag)
		
		#self.wait2finish(self.NCPUs, tag)
		self.readLaminar(tag)
		print('3D plasma field line tracing complete')
		log.info('3D plasma field line tracing complete')
		
	
	def readLaminar(self, tag = None, path = None):
		"""
		Read the MAFOT outputfile and set psimin and Lc class variables
		"""
		if tag is None: tag = ''    # this tag has len(tag) = 0
		if path is None: path = self.cwd
		
		file = path + '/' + 'lam_' + tag + '.dat'
		if os.path.isfile(file): 
			lamdata = np.genfromtxt(file,comments='#')
			if self.useVertices:
				Lc = lamdata[:,3]
				psimin = lamdata[:,4]
				N = int(len(Lc)/4)
				Lc = Lc.reshape(N,4)
				psimin = psimin.reshape(N,4)
				self.Lc = Lc.mean(1)
				self.psimin = psimin.mean(1)
			else:
				self.Lc = lamdata[:,3]
				self.psimin = lamdata[:,4]
		else:
			print('MAFOT output file: ' + file + ' not found!')
			log.info('MAFOT output file: ' + file + ' not found!')
		return


	def copyAndRead(self, path, tag = None):
		"""
		Copy all input files and MAFOT laminar data from path into self.cwd
		Read the output file when finished
		"""
		if tag is None: tag = ''    # this tag has len(tag) = 0
		
		#self.writeControlFile(MapDirection)
		src = path + '/' + '_lamCTL.dat'
		dst = self.cwd + '/' + '_lamCTL.dat'
		if os.path.isfile(src): 
			shutil.copy(src, dst)

		#self.writeM3DC1supFile()
		src = path + '/' + 'm3dc1sup.in'
		dst = self.cwd + '/' + 'm3dc1sup.in'
		if os.path.isfile(src): 
			shutil.copy(src, dst)

		#self.writeCoilsupFile()

		# normalization profile data
		for tag in ['LI','LO','UI','UO']:
			src = path + '/../' + 'lam_' + tag + '.dat'
			dst = self.cwd + '/../' + 'lam_' + tag + '.dat'
			if (not os.path.isfile(dst)) & os.path.isfile(src): 
				shutil.copy(src, dst)
		
		# main data
		src = path + '/' + 'lam_' + tag + '.dat'
		dst = self.cwd + '/' + 'lam_' + tag + '.dat'
		if os.path.isfile(src): 
			shutil.copy(src, dst)
		else:
			print('MAFOT output file: ' + src + ' not found!')
			log.info('MAFOT output file: ' + src + ' not found!')
		
		print('Copy and load 3D MAFOT Laminar data from file: ' + src)
		log.info('Copy and load 3D MAFOT Laminar data from file: ' + src)
		self.readLaminar(tag)
		return
		
		
	def checkValidOutput(self):
		""" 
		Check for invalid points in the laminar run: psimin > 2
		"""
		idx = np.where(self.psimin > 2.0)[0]
		invalid = np.zeros(len(self.psimin), dtype=bool)
		invalid[idx] = True
		print('Number of points for which Laminar run could not compute psimin:', np.sum(invalid))
		log.info('Number of points for which Laminar run could not compute psimin: ' + str(np.sum(invalid)))
		return invalid
			
			
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
		
	
	def wait2finish(self, NCPUs, tag):
		import time
		print ('Waiting for job to finish...', end='')
		time.sleep(5)	# wait 5 seconds
		while(self.isProcessRunning()):
			time.sleep(60)		# wait 1 minute; no CPU usage
		print('done')
			
		if not self.isComplete():
			print('MAFOT run ended prematurely. Attempt restart...')
			subprocess.call(['mpirun','-n',str(NCPUs),'heatlaminar_mpi','-P','points.dat','_lamCTL.dat',tag])
			self.wait2finish(NCPUs, tag)
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
		self.q0 = None
		self.ep = None	# equilParams_class instance for EFIT equilibrium
		self.HFS = None	# True: use high field side SOL, False: use low field side SOL
		
		#Default inputs
		self.NCPUs = 100
		self.teProfileData = None
		self.neProfileData = None
		
    
	def allowed_class_vars(self):
		"""
		.. Writes a list of recognized class variables to HEAT object
		.. Used for error checking input files and for initialization
		.. These variables are read in from the input file. 
		.. The call is in engine_class.loadInputs
		
		heatflux3D Variables:
		----------------------------

		:Lcmin:  float, maximum scrape-off-layer connection length
		:lcfs: float, <= 1, normalized poloidal flux of last closed flux surface location
		:lqCN: float, heat flux layer width lambda_q in mm for Eich profile
		  same as in 2D, see heatfluxClass.py
		:S: float, spread width in mm into private flux region for Eich profile
		  same as in 2D, see heatfluxClass.py
		:P: float, total power in MW into SOL; same as in 2D, see heatfluxClass.py
		:radFrac: float, 0<=x<=1, fraction of radiative power
		  same as in 2D, see heatfluxClass.py
		:qBG: float, background heat flux in MW/m^2; same as in 2D, see heatfluxClass.py
		:teProfileData: None (defaults to 2.0) 
		  or string (Name of Te data file) 
		  or float (scaler for generic Te profile) 
		  or comma-separated array (Te data vs psi=linspace(0,1.1,len(Te)))
		:neProfileData: None (defaults to 0.5) 
		  or string (Name of ne data file) 
		  or float (scaler for generic ne profile) 
		  or comma-separated array (ne data vs psi=linspace(0,1.1,len(ne)))
		:kappa: float, electron heat conductivity
		:model: string in [layer, conductive, convective] to select heat flux model 
		:NCPUs: integer, number of CPUs to use in MAFOT, default is 10
		  same as in plasma3D class
		"""
		self.allowed_vars = ['Lcmin', 'lcfs', 'lqCN', 'S', 'P', 'radFrac', 'qBG', 
				'teProfileData', 'neProfileData', 'kappa', 'model','NCPUs']


	def setTypes(self):
		"""
		Set variable types for the stuff that isnt a string from the input file
		"""
		integers = ['NCPUs']
		floats = ['Lcmin', 'lcfs', 'lqCN', 'S', 'P', 'radFrac', 'qBG', 'kappa']
		bools = []
		setAllTypes(self, integers, floats, bools)
		
		# data is an array or list
		if self.teProfileData is not None:
			if '[' in self.teProfileData:
				from ast import literal_eval
				self.teProfileData = self.teProfileData.replace(' ',',')
				self.teProfileData = np.array(literal_eval(self.teProfileData))
		if self.neProfileData is not None:
			if '[' in self.neProfileData:
				from ast import literal_eval
				self.neProfileData = self.neProfileData.replace(' ',',')
				self.neProfileData = np.array(literal_eval(self.neProfileData))
		
		# check if data is just a float
		try: self.teProfileData = float(self.teProfileData)
		except: pass	# data is a file name and remains a string or None
		try: self.neProfileData = float(self.neProfileData)
		except: pass	# data is a file name and remains a string or None
			
		
	def setupNumberFormats(self, tsSigFigs=6, shotSigFigs=6):
		"""
		sets up pythonic string number formats for shot and timesteps
		"""
		self.tsFmt = "{:."+"{:d}".format(tsSigFigs)+"f}"
		self.shotFmt = "{:0"+"{:d}".format(shotSigFigs)+"d}"
		return

	
	def initializeHF3D(self, inputDir = None):
		"""
		Set up basic input vars
		"""		
		if inputDir is None: inputDir = os.getcwd()
		self.inputDir = inputDir

		self.Psol = (1 - self.radFrac) * self.P
			
		T = self.teProfileData
		ne = self.neProfileData

		# set T profile
		if T is None: T = 2							# temperature at top of pedestal in keV
		if isinstance(T, str):						# file name for T profile data
			if ('./' in T) | ('/' not in T): path = self.inputDir + '/'
			else: path = ''
			if not os.path.isfile(path + T): 
				raise RuntimeError(path + T + ' file not found!')
			print('Loading T profile data from: ' + path + T)
			log.info('Loading T profile data from: ' + path + T)
			TData = np.loadtxt(path + T)
			psiT = TData[:,0]
			T = TData[:,1]
			self.fT = scinter.UnivariateSpline(psiT, T, s = 0, ext = 'const')
		elif isinstance(T, np.ndarray):				# array of T data assuming psi = [0, 1.1]
			psiT = np.linspace(0, 1.1, len(T))
			self.fT = scinter.UnivariateSpline(psiT, T, s = 0, ext = 'const')
		else:										# any other option
			try:
				T = float(T)						# temperature at top of pedestal in keV
				self.fT = lambda x: Tprofile(x, T)	# generic temperature profile
			except:
				raise RuntimeError('Invalid T profile data')

		# set density profile
		if ne is None: ne = 0.5							# electron density at top of pedestal in 1e-20/m^3
		if isinstance(ne, str):							# file name for density profile data
			if ('./' in ne) | ('/' not in ne): path = self.inputDir + '/'
			else: path = ''
			if not os.path.isfile(path + ne): 
				raise RuntimeError(path + ne + ' file not found!')
			print('Loading ne profile data from: ' + path + ne)
			log.info('Loading ne profile data from: ' + path + ne)
			neData = np.loadtxt(path + ne)
			nePsi = neData[:,0]
			ne = neData[:,1]
			self.fn = scinter.UnivariateSpline(nePsi, ne, s = 0, ext = 'const')
		elif isinstance(ne, np.ndarray):				# array of density data assuming psi = [0, 1.1]
			nePsi = np.linspace(0, 1.1, len(ne))
			self.fn = scinter.UnivariateSpline(nePsi, ne, s = 0, ext = 'const')
		else:											# any other option
			try:
				ne = float(ne)							# density at top of pedestal in 1e-20/m^3
				self.fn = lambda x: ne + x*0					# generic density profile
			except:
				raise RuntimeError('Invalid density profile data')


	def print_settings(self):
		"""
		Print all inputs
		"""
		print('#=============================================================')
		print('#                Optical HF Variables')
		print('#=============================================================')
		print('lqCN = ' + str(self.lqCN))
		print('S = ' + str(self.S))
		print('P = ' + str(self.P))
		print('radFrac = ' + str(self.radFrac))
		print('qBG = ' + str(self.qBG))
		print('kappa = ' + str(self.kappa))
		print('#=============================================================')
		print('#                3D Plasma Variables')
		print('#=============================================================')
		print('Lcmin = ' + str(self.Lcmin))
		print('lcfs = ' + str(self.lcfs))
		print('teProfileData = ' + str(self.teProfileData))
		print('neProfileData = ' + str(self.neProfileData))
		print('model = ' + str(self.model))

		log.info('#=============================================================')
		log.info('#                Optical HF Variables')
		log.info('#=============================================================')
		log.info('lqCN = ' + str(self.lqCN))
		log.info('S = ' + str(self.S))
		log.info('P = ' + str(self.P))
		log.info('radFrac = ' + str(self.radFrac))
		log.info('qBG = ' + str(self.qBG))
		log.info('kappa = ' + str(self.kappa))
		log.info('#=============================================================')
		log.info('#                3D Plasma Variables')
		log.info('#=============================================================')
		log.info('Lcmin = ' + str(self.Lcmin))
		log.info('lcfs = ' + str(self.lcfs))
		log.info('teProfileData = ' + str(self.teProfileData))
		log.info('neProfileData = ' + str(self.neProfileData))
		log.info('model = ' + str(self.model))

	
	def updatePFCdata(self, ep, cwd):
		"""
		Update class variables that are specific for each PFC
		"""
		self.ep = ep	# equilParams_class instance for EFIT equilibrium
		self.cwd = cwd
		

	def updateLaminarData(self, psimin, Lc):
		"""
		updates member variables for psimin, connection length, 
		checks for invalid points and finds private flux region points
		"""
		self.psimin = psimin
		self.Lc = Lc
		self.N = len(self.psimin)
		self.q = np.zeros(self.N)
		self.good = self.isGoodPoint()
		self.pfr = self.isPFR()
		
		
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
		pfr = np.where((self.psimin < 1) & (self.Lc < self.Lcmin))
		mask[pfr] = True
		return mask
		
	
	def heatflux(self, DivCode, powerFrac):
		"""
		computes self.q for the chosen model
		zeroes out invalid points
		updates self.q
		sets self.q0 this first time its called

		TO DO: make this function compatible with PFC divFracs as a mesh quantity
		"""
		print('3D Heat flux model type: ' + self.model)
		log.info('3D Heat flux model type: ' + self.model)
		if self.model in ['Layer', 'layer', 'eich', 'Eich', 'heuristic']:
			if 'O' in DivCode: HFS = False		# an Outer divertor is on low-field-side
			elif 'I' in DivCode: HFS = True		# an Inner divertor is on High-Field-Side
			else: raise ValueError('PFC Divertor Code cannot be identified. Check your PFC input file')
			self.HFS = HFS	# True: use high field side SOL, False: use low field side SOL
			print('Layer width lq =', self.lqCN)
			print('PFR spread S =', self.S)
			print('LCFS at', self.lcfs)
			print('Is on HFS:', self.HFS)
			log.info('Layer width lq = ' + str(self.lqCN))
			log.info('PFR spread S = ' + str(self.S))
			log.info('LCFS at ' + str(self.lcfs))
			log.info('Is on HFS: ' + str(self.HFS))
			q = self.getq_layer()	# normalized to qmax = 1
			if self.q0 is None: self.q0 = self.scale_layer(self.lqCN, self.S, self.Psol*powerFrac, DivCode)
		elif self.model in ['conduct', 'conductive']:
			L = np.mean(self.Lc[self.psimin > self.lcfs])*1e3	# average connection length in open field line area in m
			ratio = self.lqCN/self.S
			print('Conduction length L =', format(L,'.3f'), 'm')
			print('Ratio of SOL/PFR spread:', format(ratio,'.1f'))
			print('LCFS at', self.lcfs)
			log.info('Conduction length L = ' + format(L,'.3f') + 'm')
			log.info('Ratio of SOL/PFR spread: ' + format(ratio,'.1f'))
			log.info('LCFS at ' + str(self.lcfs))
			q = self.getq_conduct(self.psimin, kappa = self.kappa, L = L, pfr = self.pfr, ratio = ratio)
			if self.q0 is None: self.q0 = self.scale_conduct(self.Psol*powerFrac, self.kappa, L, ratio)
		else:
			raise ValueError('No valid model selected')
		
		q *= self.q0
		print('Scaling Factor q0 =', self.q0)	
		print('Background qBG =', self.qBG)
		log.info('Scaling Factor q0 = ' + str(self.q0))
		log.info('Background qBG = ' + str(self.qBG))
		self.q[self.good] = q[self.good]
		self.q += self.qBG


	def getq_conduct(self, psi, kappa = 2000, T0 = 0, L = 1, limit = True, pfr = None, ratio = 3):
		"""
		Input:
		  kappa = electron heat conductivity in W/m/eV^3.5
		  T0 = electron temperature at sheath entrance near target in keV
		  L = conduction distance between target and LCFS in m
		  scale: estimate of SOL/PFR spreading, like lq/S for Eich profiles
		Output:
		  updates self.q		
		"""
		T = self.fT(psi)			# this is now temperature in keV
		
		if limit: 
			T[psi < self.lcfs] = self.fT(self.lcfs)
		if pfr is not None: T[pfr] = self.fT(1 + ratio*(1-psi[pfr]))	# treat T in PFR as if in SOL: map psi<1 to psi>1 with ratio * dpsi
		
		q = 2.0/7.0 * kappa/L * (T**3.5 - T0**3.5) * (1e+3)**3.5/1e+6   # in MW/m^2
		return q
		
		
	def scale_conduct(self, P, kappa, L, ratio, T0 = 0):
		"""
		Get scale factor q||0 (q0) for heat flux via power balance:
		(input MW = output MW)
		Ignores wall psi and just creates a profile at OMP
		Creates a dense (1000pts) grid at the midplane to get higher resolution
		integral.  Integrates q_hat with respect to psi.
		q||0 = P_div / ( 2*pi* integral(q_hat dPsi ))
		return q0		
		"""
		psiN = np.linspace(0.85, 1.2, 1000)	# this is normalized
		T = self.fT(psiN)			# this is now temperature in keV
		
		pfr = psiN < 1.0
		T[pfr] = self.fT(1.0 + ratio*(1.0-psiN[pfr]))	# treat T in PFR as if in SOL: map psi<1 to psi>1 with ratio * dpsi
		
		q_hat = 2.0/7.0 * kappa/L * (T**3.5 - T0**3.5) * (1e+3)**3.5/1e+6   # in MW/m^2
		
		psi = psiN * (self.ep.g['psiSep']-self.ep.g['psiAxis']) + self.ep.g['psiAxis']	# this is flux
		P0 = 2*np.pi * integ.simpson(q_hat, psi)
		#account for nonphysical power
		if P0 < 0: P0 = -P0
		#Scale to input power
		q0 = P/P0
		return q0 #,q_hat,psiN


	def scale_conduct2(self, P, kappa, L, lq, S, ratio, T0 = 0, pfr = 1.0, verbose=False):
		"""
		"""		
		if pfr is None: pfr = self.lcfs
		runLaminar = True
		# Get a psi range that fully covers the profile for integration. Peak location does not matter, so use s0 from psi = 1.0
		Rlcfs = self.map_R_psi(self.lcfs)
		if self.HFS:
			Rmin = Rlcfs - 20.0*lq*(1e-3)		#in m
			if Rmin < min(self.ep.g['R']): Rmin = min(self.ep.g['R'])	#if Rmin outside EFIT grid, cap at minimum R of grid
			Rmax = Rlcfs + 20.0*S*(1e-3)		#in m
			if Rmax > self.ep.g['RmAxis']: Rmax = self.ep.g['RmAxis']	#if Rmax is outside the magnetic axis, psi would increase again, so cap at axis
		else:
			Rmin = Rlcfs - 20.0*S*(1e-3)		#in m
			if Rmin < self.ep.g['RmAxis']: Rmin = self.ep.g['RmAxis']	#if Rmin is inside the magnetic axis, psi would increase again, so cap at axis
			Rmax = Rlcfs + 20.0*lq*(1e-3)		#in m
			if Rmax > max(self.ep.g['R']): Rmax = max(self.ep.g['R'])	#if Rmax is outside EFIT grid, cap at maximum R of grid

		R = np.linspace(Rmin,Rmax,1000)
		Z = self.ep.g['ZmAxis']*np.ones(R.shape)
		
		# get q_hat from laminar		
		if self.HFS: tag = 'hfs_mp'
		else: tag = 'lfs_mp'
		file = self.cwd + '/../' + 'lam_' + tag + '.dat'
		if os.path.isfile(file): runLaminar = False

		if runLaminar:
			with open(self.cwd + '/' + 'points_' + tag + '.dat','w') as f:
				for i in range(len(R)):
					f.write(str(R[i]) + '\t' + str(0.0) + '\t' + str(Z[i]) + '\n')
					
			#nproc = 10
			args = ['mpirun','-n',str(self.NCPUs),'heatlaminar_mpi','-P','points_' + tag + '.dat','_lamCTL.dat',tag]
			current_env = os.environ.copy()        #Copy the current environment (important when in appImage mode)
			if verbose == False:
				subprocess.run(args, env=current_env, cwd=self.cwd, stderr=DEVNULL)
			else:
				subprocess.run(args, env=current_env, cwd=self.cwd) #dont suppress error messages
			for f in glob.glob(self.cwd + '/' + 'log*'): os.remove(f)		#cleanup
			# move one folder down
			src = self.cwd + '/' + 'lam_' + tag + '.dat'
			dst = self.cwd + '/../' + 'lam_' + tag + '.dat'
			if os.path.isfile(src): 
				shutil.move(src, dst)
		
		if os.path.isfile(file): 
			lamdata = np.genfromtxt(file,comments='#')
			psimin = lamdata[:,4]
		else:
			print('File', file, 'not found') 
			log.info('File ' + file + ' not found') 

		idx = np.abs(psimin - pfr).argmin()
		mask = np.zeros(len(R), dtype = bool)
		if self.HFS: pfr = np.where(R > R[idx])[0]
		else: pfr = np.where(R < R[idx])[0]
		mask[pfr] = True
		
		T = self.fT(psimin)			# this is now temperature in keV		
		T[psimin < self.lcfs] = self.fT(self.lcfs)
		T[mask] = self.fT(1.0 + ratio*(1.0-psimin[mask]))	# treat T in PFR as if in SOL: map psi<1 to psi>1 with ratio * dpsi
		
		q_hat = 2.0/7.0 * kappa/L * (T**3.5 - T0**3.5) * (1e+3)**3.5/1e+6   # in MW/m^2
		
		#Menard's method
		psiN = self.ep.psiFunc.ev(R,Z)	# this is normalized
		psi = psiN * (self.ep.g['psiSep']-self.ep.g['psiAxis']) + self.ep.g['psiAxis']	# this is flux
		P0 = 2*np.pi * integ.simpson(q_hat, psi)
		#account for nonphysical power
		if P0 < 0: P0 = -P0
		#Scale to input power
		q0 = P/P0
		return q0 #,q_hat,R,psiN,psi


	def getq_layer(self):
		"""
		Computes heat flux based on the flux layer profile with lobes
		updates self.q		
		"""
		q, q0 = self.set_layer(self.psimin, self.lqCN, self.S, lcfs = self.lcfs, lobes = True)
		if np.sum(self.pfr > 0): q[self.pfr],_ = self.set_layer(self.psimin[self.pfr], self.lqCN, self.S, q0 = q0)
		return q

	
	def set_layer(self, psi, lq, S, lcfs = 1.0, q0 = 1, lobes = False):
		"""
		psi is flat array of normalized flux
		lq is heat flux width at midplane in mm
		S is the private flux region spreading in mm
		returns flat array of heat flux based on Eich profile
		"""
		x = self.map_R_psi(psi)
		xsep = self.map_R_psi(1.0)

		# this only needs to resolve the peak well, no need to cover the entire profile, in case lq and S are large
		s0 = self.map_R_psi(lcfs)
		s = self.map_R_psi(np.linspace(lcfs-0.05,lcfs+0.1,10000))
			
		qref = eich_profile(s, lq, S, s0, q0 = 1, qBG = 0, fx = 1)
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
					
		q = eich_profile(x, lq, S, x0, q0 = 1, qBG = 0, fx = 1)
		qsep = eich_profile(xsep, lq, S, x0, q0 = 1, qBG = 0, fx = 1)
		
		if lobes:
			q[psi < lcfs] = qmax
		
		return q*q0/qmax, qsep*q0/qmax


	def map_R_psi(self, psi, HFS = None):
		"""
		Map normalized poloidal flux psi to R at midplane (Z = 0)
		psi is flat array
		return R(psi)
		"""
		if HFS is None: HFS = self.HFS
		if HFS:
			R = np.linspace(self.ep.g['RmAxis'], self.ep.g['R1'], 100)
		else:
			R = np.linspace(self.ep.g['RmAxis'], self.ep.g['R1'] + self.ep.g['Xdim'], 100)
			
		Z = self.ep.g['ZmAxis']*np.ones(len(R))
		p = self.ep.psiFunc.ev(R,Z)
		
		f = scinter.UnivariateSpline(p, R, s = 0, ext = 'const')	# psi outside of spline domain return the boundary value
		return f(psi)


	def scale_layer(self, lq, S, P, DivCode, verfyScaling = True, verbose=False):
		"""
		scales HF using a part of the limiter outline in the g-file 
		q-profile is obtained using laminar and apply the heat flux layer to psimin
		Get scale factor q||0 (q0) for heat flux via power balance:
		(input MW = output MW)
		Finds strike point on surface and sets a dense grid around it
		Integrates q_perp along surface and assumes axisymmetry.
		q||0 = P_div / ( 2*pi integral(R(s) * q_perp * ds))
		return q0
		"""
		# Parameter
		srange = 0.3
		ds = 0.0001
		Nphi = 36
		
		# strike lines
		d = self.ep.strikeLines()
		if d is None: 		
			# this means inner wall limited
			s0 = 0	# inner wall at Z = Zaxis 						!!!!!!! this needs to be computed properly !!!!!!!!!
			swall = np.arange(s0-srange, s0+srange, ds)
		else:				
			# find strike point for DivCode (this is to double check and prevent missmatches)
			if 'Rin2' in d: N = 4
			else: N = 2
			Rstr = np.zeros(N)
			Zstr = np.zeros(N)
			Sstr = np.zeros(N)
			keys = ['in','out','in2','out2']
			for i in range(N):
				Rstr[i] = d['R' + keys[i]]
				Zstr[i] = d['Z' + keys[i]]
				Sstr[i] = d['swall' + keys[i]]
		
			if 'L' in DivCode:					# Lower divertor, always 2 strike points
				Rtmp = Rstr[Zstr < 0]
				Ztmp = Zstr[Zstr < 0]
				Stmp = Sstr[Zstr < 0]
			elif 'U' in DivCode:				# Upper divertor
				Rtmp = Rstr[Zstr > 0]
				Ztmp = Zstr[Zstr > 0]
				Stmp = Sstr[Zstr > 0]
		
			if len(Rtmp) < 2: 
				raise RuntimeError('No strike points found for divertor ' + DivCode)
	
			# s0 is the strike point and s1 is the "other" strike point we don't want
			if 'I' in DivCode:
				if Rtmp[0] < Rtmp[1]: s0,s1 = Stmp[0],Stmp[1]
				else: s0,s1 = Stmp[1],Stmp[0]
			elif 'O' in DivCode:
				if Rtmp[0] < Rtmp[1]: s0,s1 = Stmp[1],Stmp[0]
				else: s0,s1 = Stmp[0],Stmp[1]

			# set swall range
			if s0 < s1:		# swall goes ccw for inner and cw for outer
				dpfr = s1 - s0
				swall = np.arange(s0 - srange, s0 + 0.4*dpfr, ds)
			else:
				dpfr = s0 - s1
				swall = np.arange(s0 - 0.4*dpfr, s0 + srange, ds)

		# get R,Z and write points file
		R,Z,nR,nZ = self.ep.all_points_along_wall(swall, get_normal = True)

		# Use MAFOT to get psimin
		runLaminar = True
		tag = DivCode
		file = self.cwd + '/../' + 'lam_' + tag + '.dat'
		if os.path.isfile(file): runLaminar = False		# MAFOT data already available

		if runLaminar:
			# write points file
			with open(self.cwd + '/' + 'points_' + tag + '.dat','w') as f:
				for j in range(Nphi):
					phi = j*(360.0/Nphi)
					for i in range(len(R)):
						f.write(str(R[i]) + '\t' + str(phi) + '\t' + str(Z[i]) + '\n')
			
			# set bounding box
			bbRmin = min([R.min()-0.1, self.ep.g['wall'][:,0].min()-0.1])
			bbRmax = max([R.max()+0.1, self.ep.g['wall'][:,0].max()+0.1])
			bbZmin = min([Z.min()-0.1, self.ep.g['wall'][:,1].min()-0.1])
			bbZmax = max([Z.max()+0.1, self.ep.g['wall'][:,1].max()+0.1])
			bbLimits = str(bbRmin) + ',' + str(bbRmax) + ',' + str(bbZmin) + ',' + str(bbZmax)
			
			# call MAFOT
			args = ['mpirun','-n',str(self.NCPUs),'heatlaminar_mpi','-P','points_' + tag + '.dat','-B',bbLimits,'_lamCTL.dat',tag]
			current_env = os.environ.copy()        #Copy the current environment (important when in appImage mode)
			if verbose == False:
				subprocess.run(args, env=current_env, cwd=self.cwd, stderr=DEVNULL)
			else:
				subprocess.run(args, env=current_env, cwd=self.cwd) #dont suppress error messages
			for f in glob.glob(self.cwd + '/' + 'log*'): os.remove(f)		#cleanup
			
			# move one folder down
			src = self.cwd + '/' + 'lam_' + tag + '.dat'
			dst = self.cwd + '/../' + 'lam_' + tag + '.dat'
			if os.path.isfile(src): 
				shutil.move(src, dst)
		
		# Read MAFOT data
		if os.path.isfile(file): 
			lamdata = np.genfromtxt(file,comments='#')
			Lc = lamdata[:,3]
			psimin = lamdata[:,4]
			#BR = lamdata[:,6]
			#BZ = lamdata[:,7]
			#Bt = lamdata[:,8]
		else:
			print('File', file, 'not found') 
			log.info('File ' + file + ' not found') 

		# Find PFR
		mask = np.zeros(len(psimin), dtype = bool)
		pfr = np.where((psimin < 1) & (Lc < self.Lcmin))
		mask[pfr] = True

		# get parallel heat flux
		qpar, q0tmp = self.set_layer(psimin, lq, S, lcfs = self.lcfs, lobes = True)
		if np.sum(mask > 0): qpar[mask],_ = self.set_layer(psimin[mask], lq, S, q0 = q0tmp)
		qpar = qpar.reshape(Nphi,len(R))

		# filter for outliers
		threshold = qpar.max(1).mean() + qpar.max(1).std()	# get max at each angle. outliers are > than average maximum + one standard deviation of all maxima
		idx = np.where(qpar > threshold)

		if verfyScaling:
			with open(self.cwd + '/../' + 'qpar_' + tag + '.dat','w') as f:
				f.write('# Parallel heat flux along g-file limiter for this divertor at multiple toroidal angles\n')
				f.write('# The field line tracing is in file: ' + file + '\n')			
				f.write('# Nphi = ' + str(Nphi) + '\n')
				f.write('# lq = ' + str(lq) + '\n')
				f.write('# S = ' + str(S) + '\n')
				f.write('# lcfs = ' + str(self.lcfs) + '\n')
				f.write('# Lcmin = ' + str(self.Lcmin) + '\n')
				f.write('# Number of outliers = ' + str(len(idx[0])) + '\n')
				f.write('# Filter threshold = ' + str(threshold) + '  Use: idx = np.where(qpar > threshold)' + '\n')
				if len(idx[0]) > 0: f.write('# Average outlier value = ' + str(qpar[idx].mean()) + '\n')
				else: f.write('# Average outlier value = None\n')
				f.write('# Each column of qpar is at another angle with phi[i] = i * 2pi/Nphi\n')
				f.write('#\n')
				f.write('# R[m]  Z[m]  swall[m]  qpar[phi,swall]\n')
				f.write('#\n')
				for i in range(len(R)):
					f.write(str(R[i]) + '\t' + str(Z[i]) + '\t' + str(swall[i]))
					for j in range(Nphi):
						f.write('\t' + str(qpar[j,i]))
					f.write('\n')

		if len(idx[0]) > 0:
			print('Filtering outliers: number = ' + str(len(idx[0])) + ', threshold = ' + str(threshold) + ', <outlier value> = ' + str(qpar[idx].mean()))
			log.info('Filtering outliers: number = ' + str(len(idx[0])) + ', threshold = ' + str(threshold) + ', <outlier value> = ' + str(qpar[idx].mean()))
			qpar[idx] = 0
		else:
			print('Filtering outliers: None found')
			log.info('Filtering outliers: None found')
			
		# average over the toroidal angles
		qparm = qpar.mean(0)

		# get incident angle
		BR = self.ep.BRFunc.ev(R,Z)
		Bt = self.ep.BtFunc.ev(R,Z)
		BZ = self.ep.BZFunc.ev(R,Z)
		
		#BR = BR.reshape(Nphi,len(R)); BR = BR.mean(0)
		#BZ = BZ.reshape(Nphi,len(R)); BZ = BZ.mean(0)
		#Bt = Bt.reshape(Nphi,len(R)); Bt = Bt.mean(0)

		B = np.sqrt(BR**2 + Bt**2 + BZ**2)
		nB = np.abs(nR*BR + nZ*BZ)/B
		
		# perpendicular heat flux
		q = qparm*nB
		
		# Integrate along line and along toroidal angle (axisymm) to get total power
		P0 = 2*np.pi * integ.simpson(R*q, swall)
		#account for nonphysical power
		if P0 < 0: P0 = -P0
		#Scale to input power
		q0 = P/P0
		return q0	#, q,mask,nB,qpar


	def scale_layer_circle(self, lq, S, P, DivCode):
		"""
		DEPRECATED
		This gives inconsistent results depending on where surface is placed
		scales HF using a circular surface 5cm below/above lower/upper x-point 
		q-profile is obtained using laminar and apply the heat flux layer to psimin
		Get scale factor q||0 (q0) for heat flux via power balance:
		(input MW = output MW)
		Finds strike point on circular poloidal line (radius r from magnetic axis) and sets a dense grid around it
		Integrates q_perp along line and assumes axisymmetry.
		q||0 = P_div / ( 2*pi*r integral(R(theta) * q_perp * dtheta))
		return q0
		"""
		# Define circular line
		from scipy.optimize import bisect
		
		if 'L' in DivCode:					# Lower divertor
			idx_Xpt = self.ep.g['lcfs'][:,1].argmin()
		elif 'U' in DivCode:				# Upper divertor
			idx_Xpt = self.ep.g['lcfs'][:,1].argmax()
			
		Rxpt = self.ep.g['lcfs'][idx_Xpt,0]
		Zxpt = self.ep.g['lcfs'][idx_Xpt,1]
		thetaxpt = self.ep.__get_theta__(Rxpt, Zxpt)
		radius = 0.05 + np.sqrt((Rxpt - self.ep.g['RmAxis'])**2 + (Zxpt - self.ep.g['ZmAxis'])**2)
		dth = 1e-4/radius
	
		f = lambda x: np.float64(self.ep.psiFunc.ev(radius*np.cos(x) + self.ep.g['RmAxis'],radius*np.sin(x) + self.ep.g['ZmAxis'])) - 1
			
		if 'L' in DivCode:
			if 'I' in DivCode:
				x0 = bisect(f,np.pi,thetaxpt)
				theta = np.arange(x0 - 0.2/radius, thetaxpt, dth)	# about 20 cm away from strike point to center of pfr
			elif 'O' in DivCode:
				x0 = bisect(f,thetaxpt,2*np.pi)
				theta = np.arange(thetaxpt, x0 + 0.2/radius, dth)
		elif 'U' in DivCode:
			if 'I' in DivCode:
				x0 = bisect(f,thetaxpt,np.pi)
				theta = np.arange(thetaxpt, x0 + 0.2/radius, dth)
			elif 'O' in DivCode:
				x0 = bisect(f,0,thetaxpt)
				theta = np.arange(x0 - 0.2/radius, thetaxpt, dth)

		R = radius*np.cos(theta) + self.ep.g['RmAxis']
		Z = radius*np.sin(theta) + self.ep.g['ZmAxis']
		
		# Use MAFOT to get psimin
		runLaminar = True
		tag = DivCode
		file = self.cwd + '/../' + 'lam_' + tag + '.dat'
		if os.path.isfile(file): runLaminar = False		# MAFOT data already available

		if runLaminar:
			# write points file
			with open(self.cwd + '/' + 'points_' + tag + '.dat','w') as f:
				for j in range(5):
					phi = j*(360/5.0)
					for i in range(len(R)):
						f.write(str(R[i]) + '\t' + str(phi) + '\t' + str(Z[i]) + '\n')
			
			# set bounding box
			bbRmin = min([R.min()-0.1, self.ep.g['wall'][:,0].min()-0.1])
			bbRmax = max([R.max()+0.1, self.ep.g['wall'][:,0].max()+0.1])
			bbZmin = min([Z.min()-0.1, self.ep.g['wall'][:,1].min()-0.1])
			bbZmax = max([Z.max()+0.1, self.ep.g['wall'][:,1].max()+0.1])
			bbLimits = str(bbRmin) + ',' + str(bbRmax) + ',' + str(bbZmin) + ',' + str(bbZmax)
			
			# call MAFOT
			args = ['mpirun','-n',str(self.NCPUs),'heatlaminar_mpi','-P','points_' + tag + '.dat','-B',bbLimits,'_lamCTL.dat',tag]
			current_env = os.environ.copy()        #Copy the current environment (important when in appImage mode)
			subprocess.run(args, env=current_env, cwd=self.cwd, stderr=DEVNULL)
			for f in glob.glob(self.cwd + '/' + 'log*'): os.remove(f)		#cleanup
			
			# move one folder down
			src = self.cwd + '/' + 'lam_' + tag + '.dat'
			dst = self.cwd + '/../' + 'lam_' + tag + '.dat'
			if os.path.isfile(src): 
				shutil.move(src, dst)
		
		# Read MAFOT data
		if os.path.isfile(file): 
			lamdata = np.genfromtxt(file,comments='#')
			Lc = lamdata[:,3]
			psimin = lamdata[:,4]
		else:
			print('File', file, 'not found') 
			log.info('File ' + file + ' not found') 

		# Find PFR
		mask = np.zeros(len(psimin), dtype = bool)
		pfr = np.where((psimin < 1) & (Lc < self.Lcmin))
		mask[pfr] = True

		# get parallel heat flux
		qpar, q0tmp = self.set_layer(psimin, lq, S, lcfs = self.lcfs, lobes = True)
		if np.sum(mask > 0): qpar[mask],_ = self.set_layer(psimin[mask], lq, S, q0 = q0tmp)

		# average over the toroidal angles
		qpar = qpar.reshape(5,len(R))
		qpar = qpar.mean(0)
		
		# get incident angle
		BR = self.ep.BRFunc.ev(R,Z)
		Bt = self.ep.BtFunc.ev(R,Z)
		BZ = self.ep.BZFunc.ev(R,Z)
		B = np.sqrt(BR**2 + Bt**2 + BZ**2)
		
		nR = self.ep.g['RmAxis'] - R
		nZ = self.ep.g['ZmAxis'] - Z
		norm = np.sqrt(nR**2 + nZ**2)
		nR = nR/norm
		nZ = nZ/norm
		
		nB = np.abs(nR*BR + nZ*BZ)/B
		
		# perpendicular heat flux
		q = qpar*nB
		
		# Integrate along line and along toroidal angle (axisymm) to get total power
		P0 = 2*np.pi * integ.simpson(radius*R*q, theta)
		#account for nonphysical power
		if P0 < 0: P0 = -P0
		#Scale to input power
		q0 = P/P0
		return q0	#, q,mask,nB,qpar


	def scale_layer_mpVar(self, lq, S, P):
		"""
		DEPRECATED
		scales HF using a R-profile along the midplane at phi = 0
		q-profile is obtained using laminar and apply the heat flux layer to psimin
		Get scale factor q||0 (q0) for heat flux via power balance:
		(input MW = output MW)
		Creates a dense (1000pts) R-grid at the midplane (Z = Zaxis) to get higher resolution
		integral.  Integrates q_hat with respect to psi.
		q||0 = P_div / ( 2*pi* integral(q_hat dPsi ))
		return q0
		"""		
		# Parameter
		dR = 0.0001		#(20*lq + 20*S)*(1e-6)		# 1000 points over the range of Rlcfs-20*lq <-> Rlcfs+20*S
		Nphi = 5

		# Get a psi range that fully covers the profile for integration. Peak location does not matter, so use s0 from psi = 1.0
		Rlcfs = self.map_R_psi(self.lcfs)
		if self.HFS:
			Rmin = min(self.ep.g['R']) + 0.01
			#Rmin = Rlcfs - 20.0*lq*(1e-3)		#in m
			#if Rmin < min(self.ep.g['R']): Rmin = min(self.ep.g['R'])	#if Rmin outside EFIT grid, cap at minimum R of grid
			Rmax = Rlcfs + 20.0*S*(1e-3)		#in m
			if Rmax > self.ep.g['RmAxis']: Rmax = self.ep.g['RmAxis']	#if Rmax is outside the magnetic axis, psi would increase again, so cap at axis
		else:
			Rmin = Rlcfs - 20.0*S*(1e-3)		#in m
			if Rmin < self.ep.g['RmAxis']: Rmin = self.ep.g['RmAxis']	#if Rmin is inside the magnetic axis, psi would increase again, so cap at axis
			#Rmax = Rlcfs + 20.0*lq*(1e-3)		#in m
			#if Rmax > max(self.ep.g['R']): Rmax = max(self.ep.g['R'])	#if Rmax is outside EFIT grid, cap at maximum R of grid
			Rmax = max(self.ep.g['R']) - 0.01

		R = np.arange(Rmin,Rmax,dR)
		Z = self.ep.g['ZmAxis']*np.ones(R.shape)
		
		# get q_hat from laminar		
		runLaminar = True
		if self.HFS: tag = 'hfs_mp'
		else: tag = 'lfs_mp'
		file = self.cwd + '/../' + 'lam_' + tag + '.dat'
		if os.path.isfile(file): runLaminar = False

		if runLaminar:
			with open(self.cwd + '/' + 'points_' + tag + '.dat','w') as f:
				for j in range(Nphi):
					phi = j*(360.0/Nphi)
					for i in range(len(R)):
						f.write(str(R[i]) + '\t' + str(phi) + '\t' + str(Z[i]) + '\n')
			
			#nproc = 10
			args = ['mpirun','-n',str(self.NCPUs),'heatlaminar_mpi','-P','points_' + tag + '.dat','_lamCTL.dat',tag]
			current_env = os.environ.copy()        #Copy the current environment (important when in appImage mode)
			subprocess.run(args, env=current_env, cwd=self.cwd, stderr=DEVNULL)
			for f in glob.glob(self.cwd + '/' + 'log*'): os.remove(f)		#cleanup
			# move one folder down
			src = self.cwd + '/' + 'lam_' + tag + '.dat'
			dst = self.cwd + '/../' + 'lam_' + tag + '.dat'
			if os.path.isfile(src): 
				shutil.move(src, dst)
		
		if os.path.isfile(file): 
			lamdata = np.genfromtxt(file,comments='#')
			psimin = lamdata[:,4]
		else:
			print('File', file, 'not found') 
			log.info('File ' + file + ' not found') 

		
		qpar,_ = self.set_layer(psimin, lq, S, lcfs = self.lcfs)

		# average over the toroidal angles
		qpar = qpar.reshape(Nphi,len(R))
		qpar = qpar.mean(0)
		
		#Menard's method
		psiN = self.ep.psiFunc.ev(R,Z)	# this is normalized
		psi = psiN * (self.ep.g['psiSep']-self.ep.g['psiAxis']) + self.ep.g['psiAxis']	# this is flux
		P0 = 2*np.pi * integ.simpson(qpar, psi)
		#account for nonphysical power
		if P0 < 0: P0 = -P0
		#Scale to input power
		q0 = P/P0
		return q0	#, q_hat,R,psiN,psi

	
	def scale_layer_mp(self, lq, S, P, pfr = 1.0):
		"""
		DEPRECATED
		scales HF using a R-profile along the midplane at phi = 0
		q-profile is obtained using laminar and apply the heat flux layer to psimin
		Get scale factor q||0 (q0) for heat flux via power balance:
		(input MW = output MW)
		Creates a dense (1000pts) R-grid at the midplane (Z = Zaxis) to get higher resolution
		integral.  Integrates q_hat with respect to psi.
		q||0 = P_div / ( 2*pi* integral(q_hat dPsi ))
		return q0
		"""		
		if pfr is None: pfr = self.lcfs
		runLaminar = True
		# Get a psi range that fully covers the profile for integration. Peak location does not matter, so use s0 from psi = 1.0
		Rlcfs = self.map_R_psi(self.lcfs)
		dR = 0.0001		#(20*lq + 20*S)*(1e-6)		# 1000 points over the range of Rlcfs-20*lq <-> Rlcfs+20*S
		if self.HFS:
			Rmin = min(self.ep.g['R']) + 0.01
			#Rmin = Rlcfs - 20.0*lq*(1e-3)		#in m
			#if Rmin < min(self.ep.g['R']): Rmin = min(self.ep.g['R'])	#if Rmin outside EFIT grid, cap at minimum R of grid
			Rmax = Rlcfs + 20.0*S*(1e-3)		#in m
			if Rmax > self.ep.g['RmAxis']: Rmax = self.ep.g['RmAxis']	#if Rmax is outside the magnetic axis, psi would increase again, so cap at axis
		else:
			Rmin = Rlcfs - 20.0*S*(1e-3)		#in m
			if Rmin < self.ep.g['RmAxis']: Rmin = self.ep.g['RmAxis']	#if Rmin is inside the magnetic axis, psi would increase again, so cap at axis
			#Rmax = Rlcfs + 20.0*lq*(1e-3)		#in m
			#if Rmax > max(self.ep.g['R']): Rmax = max(self.ep.g['R'])	#if Rmax is outside EFIT grid, cap at maximum R of grid
			Rmax = max(self.ep.g['R']) - 0.01

		R = np.arange(Rmin,Rmax,dR)
		Z = self.ep.g['ZmAxis']*np.ones(R.shape)
		
		# get q_hat from laminar		
		if self.HFS: tag = 'hfs_mp'
		else: tag = 'lfs_mp'
		file = self.cwd + '/../' + 'lam_' + tag + '.dat'
		if os.path.isfile(file): runLaminar = False

		if runLaminar:
			with open(self.cwd + '/' + 'points_' + tag + '.dat','w') as f:
				for i in range(len(R)):
					f.write(str(R[i]) + '\t' + str(0.0) + '\t' + str(Z[i]) + '\n')
			
			#nproc = 10
			args = ['mpirun','-n',str(self.NCPUs),'heatlaminar_mpi','-P','points_' + tag + '.dat','_lamCTL.dat',tag]
			current_env = os.environ.copy()        #Copy the current environment (important when in appImage mode)
			subprocess.run(args, env=current_env, cwd=self.cwd, stderr=DEVNULL)
			for f in glob.glob(self.cwd + '/' + 'log*'): os.remove(f)		#cleanup
			# move one folder down
			src = self.cwd + '/' + 'lam_' + tag + '.dat'
			dst = self.cwd + '/../' + 'lam_' + tag + '.dat'
			if os.path.isfile(src): 
				shutil.move(src, dst)
		
		if os.path.isfile(file): 
			lamdata = np.genfromtxt(file,comments='#')
			psimin = lamdata[:,4]
		else:
			print('File', file, 'not found') 
			log.info('File ' + file + ' not found') 

		idx = np.abs(psimin - pfr).argmin()
		mask = np.zeros(len(R), dtype = bool)
		if self.HFS: pfr = np.where(R > R[idx])[0]
		else: pfr = np.where(R < R[idx])[0]
		mask[pfr] = True
		
		q_hat, q0tmp = self.set_layer(psimin, lq, S, lcfs = self.lcfs, lobes = True)
		if np.sum(mask > 0): q_hat[mask],_ = self.set_layer(psimin[mask], lq, S, q0 = q0tmp)
		
		#Menard's method
		psiN = self.ep.psiFunc.ev(R,Z)	# this is normalized
		psi = psiN * (self.ep.g['psiSep']-self.ep.g['psiAxis']) + self.ep.g['psiAxis']	# this is flux
		P0 = 2*np.pi * integ.simpson(q_hat, psi)
		#account for nonphysical power
		if P0 < 0: P0 = -P0
		#Scale to input power
		q0 = P/P0
		return q0	#, q_hat,R,psiN,psi


	def fluxConversion(self, R):
		"""
		Returns the transformation factor xfm between the midplane distance s and the divertor flux psi.
		This also accounts for the flux expansion.
		"""
		if isinstance(R, np.ndarray): Z = self.ep.g['ZmAxis']*np.ones(R.shape)
		else: Z = self.ep.g['ZmAxis']
		Bp = self.ep.BpFunc.ev(R,Z)
		deltaPsi = np.abs(self.ep.g['psiSep'] - self.ep.g['psiAxis'])
		gradPsi = Bp*R
		xfm = gradPsi / deltaPsi
		return xfm
		
	
	def scale_layer2D(self, lq, S, P, HFS = False):
		"""
		DEPRECATED
		scales HF using a 2D profile, as if lcfs = 1.0
		Get scale factor q||0 (q0) for heat flux via power balance:
		(input MW = output MW)
		Ignores wall psi and just creates a profile at OMP
		Creates a dense (1000pts) grid at the midplane to get higher resolution
		integral.  Integrates q_hat with respect to psi.
		q||0 = P_div / ( 2*pi* integral(q_hat dPsi ))
		return q0
		"""		
		# Get a psi range that fully covers the profile for integration. Peak location does not matter, so use s0 from psi = 1.0
		if HFS:
			Rsep = self.ep.g['lcfs'][:,0].min()
			Rmin = Rsep - 20.0*lq*(1e-3)		#in m
			if Rmin < min(self.ep.g['R']): Rmin = min(self.ep.g['R'])	#if Rmin outside EFIT grid, cap at minimum R of grid
			Rmax = Rsep + 20.0*S*(1e-3)		#in m
			if Rmax > self.ep.g['RmAxis']: Rmax = self.ep.g['RmAxis']	#if Rmax is outside the magnetic axis, psi would increase again, so cap at axis
		else:
			Rsep = self.ep.g['lcfs'][:,0].max()
			Rmin = Rsep - 20.0*S*(1e-3)		#in m
			if Rmin < self.ep.g['RmAxis']: Rmin = self.ep.g['RmAxis']	#if Rmin is inside the magnetic axis, psi would increase again, so cap at axis
			Rmax = Rsep + 20.0*lq*(1e-3)		#in m
			if Rmax > max(self.ep.g['R']): Rmax = max(self.ep.g['R'])	#if Rmax is outside EFIT grid, cap at maximum R of grid

		R = np.linspace(Rmin,Rmax,1000)
		Z = self.ep.g['ZmAxis']*np.ones(R.shape)
		psiN = self.ep.psiFunc.ev(R,Z)	# this is normalized
		
		xfm = self.fluxConversion(R)
		q_hat = eich_profile(psiN, lq, S, 1.0, q0 = 1, qBG = 0, fx = xfm)	# becomes profile of psi by using xfm factor
		
		#Menard's method
		psi = psiN * (self.ep.g['psiSep']-self.ep.g['psiAxis']) + self.ep.g['psiAxis']	# this is flux
		P0 = 2*np.pi * integ.simpson(q_hat, psi)
		#account for nonphysical power
		if P0 < 0: P0 = -P0
		#Scale to input power
		q0 = P/P0
		return q0


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


def setAllTypes(obj, integers, floats, bools):
	"""
	Set data types for variales in obj
	"""
	for var in integers:
		if (getattr(obj, var) is not None) and (~np.isnan(float(getattr(obj, var)))):
			try:
				setattr(obj, var, tools.makeInt(getattr(obj, var)))
			except:
				print("Error with input file var "+var+".  Perhaps you have invalid input values?")
				log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")
	for var in floats:
		if var is not None:
			if (getattr(obj, var) is not None) and (~np.isnan(float(getattr(obj, var)))):
				try:
					setattr(obj, var, tools.makeFloat(getattr(obj, var)))
				except:
					print("Error with input file var "+var+".  Perhaps you have invalid input values?")
					log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")
	for var in bools:
		try:
			setattr(obj, var, tools.makeBool(getattr(obj, var)))
		except:
			print("Error with input file var "+var+".  Perhaps you have invalid input values?")
			log.info("Error with input file var "+var+".  Perhaps you have invalid input values?")


def eich_profile(s, lq, S, s0, q0, qBG = 0, fx = 1):
	"""
	Based on the paper: T.Eich et al.,PRL 107, 215001 (2011)
	lq is heat flux width at midplane in mm
	S is the private flux region spreading in mm
	s0 is the separatrix location at Z = 0 in m
	q0 is the amplitude
	qBG is the background heat flux
	fx is the flux expansion between outer midplane and target plate
	s is in m
	return function q(s)
	
	in Eich paper: s and s0 are distances along target, mapped from midplane using flux expansion fx,
	so: s = s_midplane * fx; same for s0, with s0 the position of strikeline on target
	Here, use s_midplane directly, so set fx = 1 and identify s = s_midplane = R and s0 = Rsep
	"""
	from scipy.special import erfc
	lq *= 1e-3		# in m now
	S *= 1e-3		# in m now
	a = lq*fx
	b = 0.5*S/lq
	c = S*fx
	q = 0.5 * q0 * np.exp(b**2 - (s-s0)/a) * erfc(b - (s-s0)/c) + qBG
	return q


def readShadowFile(f, PFC):
	"""
	read shadowMask.csv file
	"""
	#f = base + self.tsFmt.format(t) + '/' + PFC.name + '/shadowMask.csv
	try:
		df = pd.read_csv(f, names=['X','Y','Z','shadowMask'], skiprows=[0])
		if len(df['shadowMask'].values) != len(PFC.centers):
			print('shadowMask file mesh is not same length as STL file mesh.')
			print('Will not assign shadowMask to mismatched mesh')
			print("File length: {:d}".format(len(df['shadowMask'].values)))
			print("PFC STL mesh length: {:d}".format(len(PFC.centers)))
			val = -1
		else:
			PFC.shadowed_mask = df['shadowMask'].values
			print("Loaded Shadow Mask from file: "+f)
			val = 0
	except:
		print("COULD NOT READ FILE: "+f)
		print("Please point HEAT to a valid file,")
		print("which should be a .csv file with (X,Y,Z,shadowMask)")
		val = -1

	return val


def gridMaker(ep, DivCode, ds = 0.0001, srange = 0.3, write = False, cwd = None, Nphi = 1, plotme = True):
	"""
	Generates the grid on the g-file limiter for field line tracing to be used in scaling the heat flux models.
	"""
	# strike lines
	d = ep.strikeLines()
	if d is None: 		
		# this means inner wall limited
		s0 = 0	# inner wall at Z = Zaxis 						!!!!!!! this needs to be computed properly !!!!!!!!!
		swall = np.arange(s0-srange, s0+srange, ds)
	else:				
		# find strike point for DivCode (this is to double check and prevent missmatches)
		if 'Rin2' in d: N = 4
		else: N = 2
		Rstr = np.zeros(N)
		Zstr = np.zeros(N)
		Sstr = np.zeros(N)
		keys = ['in','out','in2','out2']
		for i in range(N):
			Rstr[i] = d['R' + keys[i]]
			Zstr[i] = d['Z' + keys[i]]
			Sstr[i] = d['swall' + keys[i]]
		
		if 'L' in DivCode:					# Lower divertor, always 2 strike points
			Rtmp = Rstr[Zstr < 0]
			Ztmp = Zstr[Zstr < 0]
			Stmp = Sstr[Zstr < 0]
		elif 'U' in DivCode:				# Upper divertor
			Rtmp = Rstr[Zstr > 0]
			Ztmp = Zstr[Zstr > 0]
			Stmp = Sstr[Zstr > 0]
		
		if len(Rtmp) < 2: 
			raise RuntimeError('No strike points found for divertor ' + DivCode)
	
		# s0 is the strike point and s1 is the "other" strike point we don't want
		if 'I' in DivCode:
			if Rtmp[0] < Rtmp[1]: s0,s1 = Stmp[0],Stmp[1]
			else: s0,s1 = Stmp[1],Stmp[0]
		elif 'O' in DivCode:
			if Rtmp[0] < Rtmp[1]: s0,s1 = Stmp[1],Stmp[0]
			else: s0,s1 = Stmp[0],Stmp[1]

		# set swall range
		if s0 < s1:		# swall goes ccw for inner and cw for outer
			dpfr = s1 - s0
			swall = np.arange(s0 - srange, s0 + 0.4*dpfr, ds)
		else:
			dpfr = s0 - s1
			swall = np.arange(s0 - 0.4*dpfr, s0 + srange, ds)

	# get R,Z and write points file
	R,Z,nR,nZ = ep.all_points_along_wall(swall, get_normal = True)

	if write:
		if cwd is None: cwd = os.getcwd()
		with open(cwd + '/' + 'points_' + DivCode + '.dat','w') as f:
			for j in range(Nphi):
				phi = j*(360.0/Nphi)
				for i in range(len(R)):
					f.write(str(R[i]) + '\t' + str(phi) + '\t' + str(Z[i]) + '\n')

	# plot stuff
	if plotme:
		import matplotlib.pyplot as plt
		ep.plot()
		plt.plot(R,Z,'r-')
	
		bbRmin = np.round(min([R.min()-0.1, ep.g['wall'][:,0].min()-0.1]),3)
		bbRmax = np.round(max([R.max()+0.1, ep.g['wall'][:,0].max()+0.1]),3)
		bbZmin = np.round(min([Z.min()-0.1, ep.g['wall'][:,1].min()-0.1]),3)
		bbZmax = np.round(max([Z.max()+0.1, ep.g['wall'][:,1].max()+0.1]),3)
		bbLimits = str(bbRmin) + ',' + str(bbRmax) + ',' + str(bbZmin) + ',' + str(bbZmax)
		print(bbLimits)
	
		plt.plot([bbRmin,bbRmax],[bbZmin,bbZmin],'g--')
		plt.plot([bbRmax,bbRmax],[bbZmin,bbZmax],'g--')
		plt.plot([bbRmin,bbRmax],[bbZmax,bbZmax],'g--')
		plt.plot([bbRmin,bbRmin],[bbZmin,bbZmax],'g--')
	
		plt.xlim(bbRmin-0.05,bbRmax+0.05)
		plt.ylim(bbZmin-0.05,bbZmax+0.05)

	return R,Z,nR,nZ,swall


def checkScaling(file, plotme =True):
	"""
	Reads the output file for heat flux scaling and plots the results to allow for visual verification
	"""
	data = np.genfromtxt(file,comments='#')
	R = data[:,0] 
	Z = data[:,1] 
	swall = data[:,2] 
	qpar = data[:,3::].T	# transpose to be consistent with array structure inside the class 
	Ns, Nphi = qpar.shape
	qparm = qpar.mean(0)
	
	threshold = qpar.max(1).mean() + qpar.max(1).std()	# get max at each angle. outliers are > than average maximum + one standard deviation of all maxima
	idx = np.where(qpar > threshold)	
	if len(idx[0]) > 0:
		print('Filtering outliers: number = ' + str(len(idx[0])) + ', threshold = ' + str(threshold) + ', <outlier value> = ' + str(qpar[idx].mean()))
		qpar[idx] = 0
	else:
		print('Filtering outliers: None found')
	
	# plot stuff
	if plotme:
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(swall,qparm,'k-', label = 'original')
		if len(idx[0]) > 0: plt.plot(swall,qpar.mean(0),'r-', label = 'filtered')
		plt.legend(fontsize = 14)
		plt.xlabel('S$_{wall}$ [m]')
		plt.ylabel('q$_{||}$ [a.u.]')
		
	return swall, qpar
	
	
	
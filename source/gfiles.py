#gfiles.py
#Description:   Grabs multiple gfiles from MDS+
#Engineer:      T Looby
#Date:          20190806

import numpy as np
import os
import shutil
import time
import sys


# ---write_mds_gfile -------------------------------------------
"""
Adapted from A. Wingen's EFIT.equilParams class for NSTX-U
reads g-file from MDS+ and writes g-file
Note: MDS+ data is only single precision!
specify shot, time (both as int) and ...
  tree (string)       ->  EFIT tree name, default = 'EFIT01'
further keywords:
  exact (bool)        ->  True: time must match time in EFIT tree, otherwise abort
                       False: EFIT time closest to time is used (default)
  Server (string)     ->  MDS+ server name or IP, default = 'atlas.gat.com' (for DIII-D)
  gpath (string)      ->  path where to save g-file, default = current working dir
"""
def write_gfile(shot, time, tree='EFIT01', exact=False, Server='skylark.pppl.gov',
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
#    RLIM = MDS.get(base + 'RLIM').data()[:, 0]
#    ZLIM = MDS.get(base + 'ZLIM').data()[:, 1]
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

#    dsign = np.sign(g['Ip'])
#    gsign = np.sign( (g['psiAxis'] - g['psiSep']) )
#    qsign = np.sign(g['Fpol'][-1]) #F_edge sign
#    g['FFprime'] *= dsign*gsign
#    g['Pprime'] *= dsign*gsign
#    g['psiRZ'] *= dsign*gsign
#    g['qpsi'] *= -qsign*dsign
#    g['Ip'] *= dsign
#    g['psiSep'] *= dsign*gsign
#    g['psiAxis'] *= dsign*gsign

#    print('dsign = {:f}'.format(dsign))
#    print('gsign = {:f}'.format(gsign))
#    print('qsign = {:f}'.format(qsign))
#    g['FFprime'] *= dsign*qsign*gsign
#    g['Pprime'] *= dsign*qsign*gsign
#    g['psiRZ'] *= dsign*qsign*gsign
#    g['qpsi'] *= -qsign*dsign
#    g['psiSep'] *= dsign*qsign*gsign
#    g['psiAxis'] *= dsign*qsign*gsign
#    g['Fpol'] *= dsign*qsign

    # Now, write to file using same style as J. Menard script (listed above)
    # Using function in WRITE_GFILE for reference
    with open(gpath + 'g' + format(shot, '06d') + '.' + format(time,'05d'), 'w') as f:
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
    return time
# --- _write_array -----------------------
# write numpy array in format used in g-file:
# 5 columns, 9 digit float with exponents and no spaces in front of negative numbers
def write_array(x, f):
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


def write_multiple_gfiles(machine, shot, tree='efit01', tmin=None, tmax=None,
                          rootDir=None, clobberwait=True):
    """
    Function for grabbing multiple gfiles associated with a single shot
    Mainly for time varying equilibrium stuff
    machine: string. either 'nstx' or 'd3d'
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

    try:
        # Create target Directory
        os.mkdir(dirName)
    except FileExistsError:
        #Directory Clobber checking enabled if clobberwait = True
        clobberlist = ['y','Y','yes','YES','Yes']
        if clobberwait:
            print("\nDirectory " , dirName ,  " already exists!")
            clobberflag = input('Should I overwrite this directory?  (DANGER!!!) [y/n]  ')
        else: clobberflag = 'y'
        #To clobber or not to clobber...
        if clobberflag in clobberlist:
            try: shutil.rmtree(dirName)
            except OSError as e:
                print ("Error: %s - %s." % (e.filename, e.strerror))
                sys.exit()

                os.mkdir(dirName)
                print("Directory " , dirName ,  " Created ")
            os.mkdir(dirName)
            preserveFlag = False
        else:
            ts = np.array(os.listdir(rootDir), dtype=int)
            preserveFlag = True
            return ts, preserveFlag

    #Pull gfiles and write into dirName/subdirectory where subdirectory
    #is exclusively for that timestep
    print('Pulling all gfiles between {:5d} and {:5d} [ms]'.format(int(tmin),int(tmax)))
    for t in ts:
        t_path = dirName + '/{:06d}'.format(t)
        os.mkdir(t_path)
        if machine == 'nstx':
            write_gfile(shot,t,tree=tree,Server=Server,gpath=t_path)
    preserveFlag = False
    return ts, preserveFlag

def loadgfile(machine,gfile,rootDir=None,clobberwait=True):
    """
    Copies existing local gfile to HEAT file tree
    clobber checking for existing gfile
    Also returns timesteps and preserveFlag
    """
    path, name = os.path.split(gfile)
    shot,time = name[:13].split('.')
    shot = int(shot.split('g')[1])
    ts = int(time)
    preserveFlag = True

    #Directory where we will save these gfiles if one is not declared
    if rootDir==None: dirName = machine + '_{:06d}'.format(shot)
    else: dirName = rootDir
    try:
        os.mkdir(dirName)
    except:
        pass

    newPath = dirName+'/{:06d}/'.format(ts)
    newgfile = newPath + name
    newgfile_noSuffix = newPath + name[:13]
    try:
        # Create target Directory
        os.mkdir(newPath)
        shutil.copyfile(gfile, newgfile)
        shutil.copyfile(gfile, newgfile_noSuffix) #MAFOT fails if gfile has suffix
    #except FileExistsError:
    except:
        #Directory Clobber checking enabled if clobberwait = True
        clobberlist = ['y','Y','yes','YES','Yes']
        if clobberwait:
            print("\nDirectory " , dirName ,  " already exists!")
            clobberflag = input('Should I overwrite this directory?  (DANGER!!!) [y/n]  ')
        else: clobberflag = 'y'
        #To clobber or not to clobber...
        if clobberflag in clobberlist:
            try: shutil.rmtree(newPath)
            except OSError as e:
                print ("Error: %s - %s." % (e.filename, e.strerror))
                sys.exit()

                os.mkdir(newPath)

            os.mkdir(newPath)
            preserveFlag = False
            shutil.copyfile(gfile, newgfile)
            shutil.copyfile(gfile, newgfile_noSuffix) #MAFOT fails if gfile has suffix
            print("Directory " , newPath ,  " Created ")
        else:
            preserveFlag = True
    return np.array([ts]), preserveFlag

if __name__ == '__main__':
    print('\nWARNING: If you call this module from CLI you need to edit')
    print('input arguments manually in the python file.\n\n')
    time.sleep(2)
    machine = 'nstx'
    tree = 'efit02'
    shot = 204118
    tmin = 950
    tmax = 965
    write_multiple_gfiles(machine,shot,tree,tmin,tmax)
    print("Pulled requested timesteps from MDS+.")
    print("Wrote gfiles.")

#launchHEAT.py
#Description:   entry point for HEAT
#Engineer:      T Looby
#Date:          20220131
"""
provides entry point to the Heat flux Engineering Analysis Toolkit (HEAT)

user runs this command from terminal to launch HEAT's graphical user
interface (GUI) or terminal user interface (TUI)

Sets up environment, depending upon user runMode: local, appImage, docker
"""
#========= VARIABLES THAT ARE SYSTEM DEPENDENT =================================
import os
import sys
import subprocess
import argparse


def loadEnviron():
    """
    loads the HEAT environment
    """
    try:
        runMode = os.environ["runMode"]
    except:
        runMode = 'local'
        os.environ["runMode"] = runMode

    #default HEAT output directory
    try:
        homeDir = os.path.expanduser("~")
    except:
        print("HOME env var not set.  Set before running HEAT!")
        print("Example:  export HOME=/home/tom")
        sys.exit()

    dataPath = homeDir + '/HEAT/data'
    OFversion = 'v2112'

    #=== Set up paths and environment vars
    ### appImage
    if runMode == 'appImage':
        print("Running in appImage mode")

        ### USER ROOT HEATDIR
        AppDir = os.environ["APPDIR"]
        #Root HEAT source code directory
        rootDir = AppDir + '/usr/src'

        ### PARAVIEW
        #Include the location of the paraview binaries
        #Specifically we need the python libs and pvpython
        PVPath = os.environ["PVPath"]
        pvpythonCMD = os.environ["pvpythonCMD"]

        ### OPENFOAM
        #openFOAM bashrc location v1912
        #OFbashrc = AppDir + '/usr/opt/openfoam/openfoam1912/etc/bashrc'
        #OFdir = AppDir+'/usr/opt/openfoam/openfoam1912'
        #openFOAM bashrc location v2106
        OFbashrc = AppDir + '/opt/openfoam/etc/bashrc'
        OFdir = AppDir+'/opt/openfoam'
        #python site packages where PyFoam resides
        pyFoamPath = AppDir + '/lib/python3.8/site-packages'

        ### FREECAD
        #default freecad path
        #FreeCADPath = AppDir + '/opt/freecad/squashfs-root/usr/lib'
        FreeCADPath = AppDir + '/usr/lib/freecad-python3/lib'

        ### ORNL EFIT MODULE
        #default source code location (EFIT class should be here)
        EFITPath = AppDir + '/usr/src'

        ### Open3D
        O3Dpath = AppDir + '/usr/lib/open3d/lib/python_package/open3d'


        ###  Docker container
    elif runMode == 'docker':
        print("Running in Docker mode")

        ### USER ROOT HEATDIR
        #Root HEAT source code directory
        rootDir = homeDir + '/source/HEAT/source'
        #default AppDir for when running in docker mode
        AppDir = os.environ["APPDIR"]

        ### PARAVIEW
        #Include the location of the paraview binaries.
        #Specifically we need the python libs and pvpython
        PVPath = homeDir + '/lib/python3.8/site-packages'
        pvpythonCMD = homeDir + '/opt/paraview/bin/pvpython'

        ### FREECAD
        #docker ubuntu repo freecad path
        FreeCADPath = '/usr/lib/freecad-python3/lib'
        #FreeCADPath = '/usr/lib/freecad-daily/lib'

        ### ORNL EFIT CLASS
        #default source code location (EFIT class should be here)
        EFITPath = homeDir + '/source'

        ### OPENFOAM
        #default openFOAM source path v1912
        #OFbashrc = '/opt/openfoam/openfoam-OpenFOAM-v1912/etc/bashrc'
        #default openFOAM source path v2106
        OFbashrc = AppDir + '/opt/openfoam/etc/bashrc'
        #python site packages where PyFoam resides
        pyFoamPath = homeDir + '/.local/lib/python3.8/site-packages'
        #pyFoam python scripts
        pyFoamPath = '/'

        ### Open3D
        O3Dpath = AppDir + '/opt/open3d/lib/python_package'

        #ENV VARS
        #create necessary environment variables when outside appImage
        os.environ["PVPath"] = PVPath
        os.environ["pvpythonCMD"] = pvpythonCMD

        #local development mode
    else:
        ###  If developing you will need to edit these manually!
        print("Running in local developer mode")
        print("You will need a manually compiled environment")
        ### USER ROOT HEATDIR
        #Root HEAT source code directory
        rootDir = homeDir + '/source/HEAT/github/source'

        ### PARAVIEW
        #Include the location of the paraview binaries.
        #Specifically we need the python libs and pvpython
        PVPath = '/opt/paraview/ParaView-5.10.1-MPI-Linux-Python3.9-x86_64/lib/python3.8/site-packages'
        pvpythonCMD = '/opt/paraview/ParaView-5.10.1-MPI-Linux-Python3.9-x86_64/bin/pvpython'

        ### FREECAD
        # daily build binary freecad path
        FreeCADPath = '/usr/lib/freecad-daily/lib'
        # downloaded appImage freecad path
        #FreeCADPath = '/opt/freecad/squashfs-root/usr/lib'
        # for ubuntu repo build
        #FreeCADPath = '/usr/lib/freecad-python3/lib'
        #FreeCADPath = '/usr/lib/freecad/lib'
        # for daily builds
        #FreeCADPath = '/usr/lib/freecad-daily-python3/lib'

        ### ORNL EFIT CLASS
        #default source code location (EFIT class should be here)
        EFITPath = homeDir + '/source'

        ### OPENFOAM
        #default openFOAM source path v1912
        #OFbashrc = '/opt/openfoam/openfoam-OpenFOAM-v1912/etc/bashrc'
        #default openFOAM source path v2112
        OFbashrc = '/opt/openfoam/OpenFOAM-v2112/etc/bashrc'
        #python site packages where PyFoam resides
        pyFoamPath = homeDir + '/.local/lib/python3.8/site-packages'
        #pyFoam python scripts
        pyFoamPath = '/'

        ### Open3D
        O3Dpath = '/opt/open3d/Open3D/build/lib/python_package/open3d'

        #ENV VARS
        #default AppDir for when running in dev mode
        AppDir = 'not in appImage mode'
        #create necessary environment variables when outside appImage
        os.environ["PVPath"] = PVPath
        os.environ["pvpythonCMD"] = pvpythonCMD

    #default logfile location
    logFile = dataPath + '/HEATlog.txt'

    #Now set the relevant environment variables
    os.environ["logFile"] = logFile
    os.environ["rootDir"] = rootDir
    os.environ["dataPath"] = dataPath
    os.environ["OFbashrc"] = OFbashrc
    os.environ["FreeCADPath"] = FreeCADPath
    os.environ["HEATchmod"] = '0o774' #chmod in base 8 (octal)
    os.environ["WM_PROJECT_VERSION"] = OFversion

    #clear uname mask for docker saving
    os.umask(0)

    #===========================================================================

    #=======UPDATE PATHS========================================================
    #orca installation location (for saving EQ plots)
    #pio.orca.config.executable='/usr/bin/orca'
    #append EFIT to python path
    sys.path.append(EFITPath)
    #append FreeCAD to python path
    sys.path.append(FreeCADPath)
    #append paraview to python path
    sys.path.append(PVPath)
    #append pyFoam site-packages location to python path
    sys.path.append(pyFoamPath)
    #append pvpython to binary path
    oldEnv = os.environ["PATH"]
    #os.environ["PATH"] = oldEnv + ':' + pvpythonCMD
    #append Open3D to python path
    sys.path.append(O3Dpath)
    #===============================================================================

    #Create dataPath
    if not os.path.exists(dataPath):
        os.makedirs(dataPath)
    return

def getLogger():
    """
    starts logFile logging
    """
    logFile = os.environ["logFile"]
    from pathlib import Path
    Path(logFile).touch()
    import logging
    logFlask = logging.getLogger('werkzeug')
    logFlask.disabled = True
    logging.basicConfig(filename=logFile, filemode="w", level=logging.INFO, format='%(message)s')
    log = logging.getLogger(__name__)
    return log

#=======HEAT Launch Point=======================================================
def launchHEAT(args):
    mode = vars(args)['m']
    address = vars(args)['a']
    port = vars(args)['p']
    batchFile = vars(args)['f']
    batchTemplate = vars(args)['sB']

    runMode = os.environ["runMode"]

    #list of tokamak flags that are options in HEAT (if adding new tokamak add flag to list)
    machineList = ['d3d','nstx','st40','step','sparc','west','kstar']

    #run HEAT in terminal mode
    if mode=='t':
        print('\nRunning HEAT via Terminal User Interface (TUI)...\n')
        log.info('\nRunning HEAT via Terminal User Interface (TUI)...\n')
        import terminalUI

        tui = terminalUI.TUI()
        tui.UImode = 't' #terminal mode
        tui.machineList = machineList

        if batchFile == None:
            if batchTemplate is not None:
                tui.saveBatchFile(batchTemplate)
                print("File saved...")
                sys.exit()
            else:
                print("No batchFile.dat found in arguments!  Required for terminal mode.  Aborting.")
                print("To save a batchFile template, run with --sB <batchPath> switch, ")
                print("where <batchPath> is where you want the batchFile saved")
                sys.exit()

        tui.simulationSchedule(batchFile)
        tui.runSimulations()

    #run HEAT in graphical mode
    else:
        print('\nRunning HEAT via Graphical User Interface (GUI)...\n')
        log.info('\nRunning HEAT via Graphical User Interface (GUI)...\n')
        import dashGUI as dashGUI
        dashGUI.generateLayout()
        dashGUI.machineList = machineList
        #use default IPv4 address and port unless user provided one
        if address == None:
            address = '127.0.0.1' #default
        if port == None:
            port = 8050 #default

        if runMode == 'local':
            dashGUI.app.run_server(
                            debug=True,
                            dev_tools_ui=True,
                            port=port,
                            host=address,
                            use_reloader=True, #this can be used in local developer mode only
                            dev_tools_hot_reload = True, #this can be used in local developer mode only
                            )
        else:
            dashGUI.app.run_server(
                            debug=True,
                            dev_tools_ui=True,
                            port=port,
                            host=address,
                            use_reloader=False, #this can be used in local developer mode only
                            dev_tools_hot_reload = False, #this can be used in local developer mode only
                            )
    return


if __name__ == '__main__':

    #parse command line arguments
    parser = argparse.ArgumentParser(description=""" HEAT...Use this command to launch HEAT.
                                                    You can run HEAT in terminal (t) or graphical (g) modes,\n
                                                     by using the --mode switch.  Default is graphical mode. """)
    parser.add_argument('--m', type=str, help='HEAT run mode', choices=['t', 'g'], required=False, default='g')
    parser.add_argument('--a', type=str, help='HEAT GUI IP address ', required=False)
    parser.add_argument('--p', type=str, help='HEAT GUI port # ', required=False)
    parser.add_argument('--f', type=str, help='Batch file path', required=False)
    parser.add_argument('--sB', type=str, help='Save Batch File', required=False)

    args = parser.parse_args()

    #initialize environment
    loadEnviron()

    #initialize logs
    log = getLogger()

    #launch HEAT
    launchHEAT(args)

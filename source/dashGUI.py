#dashGUI.py
#Description:   DASH python html gui
#Engineer:      T Looby
#Date:          20200615 (ish)
"""
This is the python - html interface, and the launch point for the HEAT code.
Running this script launches a web interface that can be used to run HEAT
The web interface is html and can be accessed via any web browser.

DASH is the python library that creates the html - python binding
Under the hood, DASH is running a flask server with associated proxy mapping,
  html, javascipt, etc.  We use decorators (@) to serve as javascript callbacks

There are css files that go with this program.  They are located in the
   ./assets directory

If you want to use this in a production environment serving multiple sessions
simultaneously, you will need to run a gunicorn server or something to isolate
different users' class variables from each other.

dashGUI.py should be called with no arguments to run on 127.0.0.1:8050 (default)
to run on address:port, use command line:
dashGUI.py <address> <port>

You will need to set a few variables below, based upon your system paths
rootDir, PVPath
"""
import os
import sys
import shutil
import base64
import io
import json
import numpy as np
import pandas as pd
#unsafe, removed 20220406
#import parser
import time
import copy
import dash
#import dash_html_components as html
from dash import html
#import dash_core_components as dcc
from dash import dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import visdcc
from flask import Flask, send_from_directory
import plotly.io as pio
#import dash_table
from dash import dash_table
import EFIT.equilParams_class as EP
import toolsClass
tools = toolsClass.tools()
from dash_extensions import Download
from dash_extensions.snippets import send_file
import ipaddress
import logging
log = logging.getLogger(__name__)

#get relevant environment variables
logFile = os.environ["logFile"]
rootDir = os.environ["rootDir"]
dataPath = os.environ["dataPath"]
OFbashrc = os.environ["OFbashrc"]
FreeCADPath = os.environ["FreeCADPath"]
PVPath = os.environ["PVPath"]
pvpythonCMD = os.environ["pvpythonCMD"]
try:
    AppDir = os.environ["APPDIR"]
except:
    AppDir = 'Not in appImage'

try:
    chmod = int(os.environ["HEATchmod"], 8) #convert from base 8
except:
    chmod = int(0o774, 8)

try:
    GID = int(os.environ["dockerGID"]) #group ID
except:
    GID = -1
try:
    UID = int(os.environ["dockerUID"]) #user ID
except:
    UID = -1

#Import HEAT engine
from engineClass import engineObj

server = Flask(__name__)
#app = dash.Dash(__name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}])
#Create our own server for downloading files

app = dash.Dash(server=server, meta_tags=[{"name": "viewport", "content": "width=device-width"}],
                prevent_initial_callbacks=False)

#Eventually need to fix this so that we are not using a global variable
#dash can acces Flask Cache so we should cache data by userID or something
#for R&D this works
gui = engineObj(logFile, rootDir, dataPath, OFbashrc, chmod, UID, GID)
gui.UImode = 'g' #graphical mode


"""
==============================================================================
Banner and Tabs
"""

def build_banner():
    return html.Div(
        id="banner",
        className="banner",
        children=[
            html.Div(
                id="banner-text",
                children=[
                    html.H4("Heat flux Engineering Analysis Toolkit (HEAT)"),
                ],
            ),
        ],
    )

# Dash doesnt let CSS overwrite tab settings so do it here
tabs_style = {
#    'width': '15%'
}
tab_style = {
    'fontWeight': 'bold',
    'color': '#ffffff',
    'backgroundColor': '#252626',
    'border': 'none',

}
tab_selected_style = {
    'backgroundColor': '#119DFF',
    'color': 'blue',
    'border': 'none',
}

def build_tabs():
    """
    returns tab bar
    """
    return html.Div(
        id="tabs",
#        className="tab-container",
        className="tabcontent",
        children=[
            dcc.Tabs(
                id="app-tabs",
                value="tab1",
                vertical=True,
                style=tabs_style,
                children=[
                    dcc.Tab(
                        id="input-tab",
                        label="Inputs",
                        value="tab1",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.Div(
#                            className='tabcontent',
                            className="innerTabContent",
                            children=[
                                buildDefaultPaths(),
                                buildButtonRibbon(),
                                buildMHDbox(),
                                buildCADbox(),
                                buildPFCbox(),
                                buildHFbox(),
                                buildGYRObox(),
                                buildOFbox(),
                                ]
                                )
                        ]
                    ),
                    dcc.Tab(
                        id="run-tab",
                        label="Run HEAT",
                        value="tab2",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.Div(
#                            className='tabcontent',
                            className="innerTabContent",
                            children=[
                                buildRunTab()
                                ]
                                )
                        ]
                    ),
                    dcc.Tab(
                        id="gfileCleaner-tab",
                        label="GEQDSK Tools",
                        value="tab3",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.Div(
#                            className='tabcontent',
                            className="innerTabContent",
                            children=[
                                buildGfileCleanerTab()
                                ]
                            )
                        ]
                    ),

                    dcc.Tab(
                        id="output-tab",
                        label="Output",
                        value="tab4",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.Div(
#                            className='tabcontent',
                            className="innerTabContentColumn",
                            children=[
                                buildOutputTab()
                                ]
                                )
                        ],
                    ),
                    dcc.Tab(
                        id="log-tab",
                        label="LogFile",
                        value="tab5",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.Div(
#                            className='tabcontent',
                            className="logTabContent",
                            children=[
                                buildLogTab()
                                ]
                                )
                        ],
                    ),


                ],
            )
        ],
    )

"""
==============================================================================
Tab Contents: inputs tab
"""
def buildButtonRibbon():
    """
    returns machine selector and load default buttons
    """
    return html.Div(
        id="buttonRib",
        className="buttonRibbon",
        children=[
            html.H5("Machine Selection and Input Files"),
            html.Div(
                children=[
                    buildMachineSelector(),
                    buildInputButtons(),
                    ],
                className="rowBox"
            )
        ],
    )

def buildDefaultPaths():
    """
    contains text boxes for HEAT relevent paths
    PVPath is path for paraview binaries and pvpython
    FreeCAD is location of freecad installation
    dataDir is location where HEAT output will be saved

    rootDir is location of HEAT source code and is not included in GUI
    because it would be impossible to run GUI if this was not already set.
    rootDir (and some other defaults) are hardcoded at the top of this file

    if className is Hidden, then the input exists in html but is hidden from user
    """
    return html.Div(
        id="defaultPaths",
        children=[
            html.Label("ParaVIEW Path", className="textInputHidden"),
            dcc.Input(id="PVPath", className="textInputHidden", value=PVPath),
            html.Label("FreeCAD Path", className="textInputHidden"),
            dcc.Input(id="FreeCADPath", className="textInputHidden", value=FreeCADPath),
            html.Label("Data Directory"),
            dcc.Input(id="dataPath", className="textInput", value=dataPath),
            html.Label("OpenFOAM bashrc file", className="textInputHidden"),
            dcc.Input(id="OFbashrc", className="textInputHidden", value=OFbashrc),
            html.Label("Image Mount Directory: "+AppDir)
        ],
        className="colBox"
    )

def buildMachineSelector():
    """
    returns machine selector dropdown
    """
    return html.Div(
            className="machineSelectBox",
            children=[
                html.Label(id="machLabel", children="Select a Tokamak"),
                dcc.Dropdown(
                    id='MachFlag',
                    className="machineSelect",
                    style={'backgroundColor': 'transparent', 'color':'transparent',
                            'align-items':'center'},
                    #style=dropdown_style,
                    options=[
                        {'label': 'NSTX-U', 'value': 'nstx'},
                        {'label': 'DIII-D', 'value': 'd3d'},
                        {'label': 'ST40', 'value': 'st40'},
                        {'label': 'STEP', 'value': 'step'},
                        {'label': 'SPARC', 'value': 'sparc'},
                        {'label': 'WEST', 'value': 'west'},
                        {'label': 'K-STAR', 'value': 'kstar'}
                        ],
                    value=None
                    ),
                html.Div(id="hiddenDivMachFlag"),
                dcc.Upload(
                    className="inputUpload",
                    id='input-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select input file')
                    ]),
                    style={
                        'width': '100%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed',
                        'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                        'align-items':'center'
                        },
                    multiple=True,
                    ),
            ],
        )

def buildInputButtons():
    """
    returns Load Defaults drag and drop and Upload Input buttons
    """
    return html.Div(
            id="buttonInputs",
            className="defaultButtonBox",
            children=[
                html.Div(id="hiddenDivInput"),
                html.Button("Load Defaults (optional)", id="loadDefaults", n_clicks=0, className="defaultButtons"),
                html.Button("Save Settings\nInto Input File", id="saveInputs",
                            n_clicks=0, className="defaultButtons"),
                Download(id="downloadInputs"),
                html.Div(id="hiddenDivSaveInput"),
            ],
        )

@app.callback(Output('hiddenDivMachFlag', 'children'),
              [Input('MachFlag', 'value')])
def machineSelector(MachFlag):
    """
    callback to handle machine selector drop down
    """
    if MachFlag == None:
        machFlagChosen = "Select a machine"
        return [html.Label(machFlagChosen, style={'color':'#fc0313'})]
    gui.machineSelect(MachFlag, machineList)
    machFlagChosen = "Selected "+MachFlag
    return [html.Label(machFlagChosen, style={'color':'#f5d142'})]


@app.callback([Output('hiddenDivInput', 'children'),
               Output('userInputFileData', 'data')],
             [Input('input-upload','filename')],
             [State('input-upload','contents'),
              State('MachFlag','value'),])
def inputDragDrop(file, contents, MachFlag):
    """
    callback to handle user input file drag and drop
    """
    if MachFlag is None:
        print("Select a machine before uploading input file")
        log.info("Select a machine before uploading input file")
        raise PreventUpdate
    if file is None:
        raise PreventUpdate
    else:
        outputDiv = html.Label("Loaded Input File", style={'color':'#f5d142'})
        newFile = gui.tmpDir + file[0]
        decoded = base64.b64decode(contents[0].split(',')[1])
        #Save user loaded file into tmp directory
        with open(newFile, "w") as f:
            f.write(decoded.decode('utf-8'))
        data = gui.loadDefaults(inFile=newFile)

    return [outputDiv, data]

@app.callback([Output('hiddenDivSaveInput','children'),
               Output('downloadInputs', 'data')],
              [Input('saveInputs','n_clicks')],
              [State('MachFlag', 'value'),
               State('shot', 'value'),
               State('tmin', 'value'),
               State('tmax', 'value'),
               State('nTrace', 'value'),
               State('gridRes', 'value'),
               State('lqEich', 'value'),
               State('S', 'value'),
               State('lqCN', 'value'),
               State('lqCF', 'value'),
               State('lqPN', 'value'),
               State('lqPF', 'value'),
               State('fracCN', 'value'),
               State('fracCF', 'value'),
               State('fracPN', 'value'),
               State('fracPF', 'value'),
               State('Pinj', 'value'),
               State('coreRadFrac', 'value'),
               State('fracUI', 'value'),
               State('fracUO', 'value'),
               State('fracLI', 'value'),
               State('fracLO', 'value'),
               State('qBG', 'value'),
               State('OFstartTime', 'value'),
               State('OFstopTime', 'value'),
               State('OFminMeshLev', 'value'),
               State('OFmaxMeshLev', 'value'),
               State('OFSTLscale', 'value'),
               State('OFdeltaT', 'value'),
               State('OFwriteDeltaT', 'value'),
               State('PVPath', 'value'),
               State('FreeCADPath', 'value'),
               State('dataPath', 'value'),
               State('OFbashrc', 'value'),
               State('N_gyroSteps','value'),
               State('gyroDeg','value'),
               State('gyroT_eV','value'),
               State('N_vSlice','value'),
               State('N_vPhase','value'),
               State('N_gyroPhase','value'),
               State('ionMassAMU','value'),
               State('vMode','value'),
               State('ionFrac','value')
               ]
               )
def saveGUIinputs(  n_clicks,
                    MachFlag,
                    shot,
                    tmin,
                    tmax,
                    nTrace,
                    gridRes,
                    lqEich,
                    S,
                    lqCN,
                    lqCF,
                    lqPN,
                    lqPF,
                    fracCN,
                    fracCF,
                    fracPN,
                    fracPF,
                    Pinj,
                    coreRadFrac,
                    fracUI,
                    fracUO,
                    fracLI,
                    fracLO,
                    qBG,
                    OFstartTime,
                    OFstopTime,
                    OFminMeshLev,
                    OFmaxMeshLev,
                    OFSTLscale,
                    OFdeltaT,
                    OFwriteDeltaT,
                    PVLoc,
                    FreeCADLoc,
                    dataLoc,
                    OFbashrcLoc,
                    N_gyroSteps,
                    gyroDeg,
                    gyroT_eV,
                    N_vSlice,
                    N_vPhase,
                    N_gyroPhase,
                    ionMassAMU,
                    vMode,
                    ionFrac
                ):
    """
    Saves GUI text boxes into an input file in the HEAT format
    first file is saved to the gui.tmpDir, then to client machine
    """
    if n_clicks < 1:
        raise PreventUpdate

    data = {}
    data['MachFlag'] = MachFlag
    data['shot'] = shot
    data['tmin'] = tmin
    data['tmax'] = tmax
    data['nTrace'] = nTrace
    data['gridRes'] = gridRes
    data['lqEich'] = lqEich
    data['S'] = S
    data['lqCN'] = lqCN
    data['lqCF'] = lqCF
    data['lqPN'] = lqPN
    data['lqPF'] = lqPF
    data['fracCN'] = fracCN
    data['fracCF'] = fracCF
    data['fracPN'] = fracPN
    data['fracPF'] = fracPF
    data['Pinj'] = Pinj
    data['coreRadFrac'] = coreRadFrac
    data['fracUI'] = fracUI
    data['fracUO'] = fracUO
    data['fracLI'] = fracLI
    data['fracLO'] = fracLO
    data['qBG'] = qBG
    data['OFtMin'] = OFstartTime
    data['OFtMax'] = OFstopTime
    data['meshMinLev'] = OFminMeshLev
    data['meshMaxLev'] = OFmaxMeshLev
    data['STLscale'] = OFSTLscale
    data['deltaT'] = OFdeltaT
    data['writeDeltaT'] = OFwriteDeltaT
    data['FreeCADPath'] = FreeCADLoc
    data['PVPath'] = PVLoc
    data['dataPath'] = dataLoc
    data['OFbashrc'] = OFbashrcLoc
    data['N_gyroSteps'] = N_gyroSteps
    data['gyroDeg'] = gyroDeg
    data['gyroT_eV'] = gyroT_eV
    data['N_vSlice'] = N_vSlice
    data['N_vPhase'] = N_vPhase
    data['N_gyroPhase'] = N_gyroPhase
    data['ionMassAMU'] = ionMassAMU
    data['vMode'] = vMode
    data['ionFrac'] = ionFrac
    tools.saveInputFile(data, gui.tmpDir, gui.rootDir, gui.dataPath)

    outputDiv = html.Label("Saved File", style={'color':'#f5d142'})
    return [outputDiv, send_file(gui.tmpDir + "HEATinput.csv")]



#==========MHD==========
def buildMHDbox():
    """
    MHD input parameters
    """
    return html.Div(
            id="MHDbox",
            children=[
                html.H6("MHD Settings"),
                html.Label(id="shotLabel", children="Shot Number  "),
                dcc.Input(id="shot", className="textInput"),
                html.Label(id="tMinLabel", children="Minimum Timestep [ms]"),
                dcc.Input(id="tmin", className="textInput"),
                html.Label(id="tMaxLabel", children="Maximum Timestep [ms]"),
                dcc.Input(id="tmax", className="textInput"),
                html.Label(id="nTraceLabel", children="Trace Distance (degrees)"),
                dcc.Input(id="nTrace", className="textInput"),
                dcc.Upload(
                    className="PFCupload",
                    id='gfiletable-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select gFiles')
                    ]),
                    style={
                        'width': '60%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed',
                        'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                        },
                    multiple=True,
                    ),
                html.Div(id="hiddenDivGfileUpload"),
                html.Br(),
                dcc.RadioItems(
                                id="plasma3Dmask",
                                options=[
                                    {'label': '2D Plasmas', 'value': 'plasma2D'},
                                    {'label': '3D Plasmas', 'value': 'plasma3D'},
                                        ],
                                        value='plasma2D'
                                ),
                html.Button("Load MHD", id="loadMHD", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
                html.Br(),
                #save EQ plots as png buttons / forms
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Label("Width (px)"),
                                dcc.Input(id="EQpixelX", className="textInput"),
                            ],
                            className="colBox"
                        ),
                        html.Div(
                            children=[
                                html.Label("Height (px)"),
                                dcc.Input(id="EQpixelY", className="textInput"),
                            ],
                            className="colBox"
                        ),
                        html.Button("Save EQs to .png", id="saveEQbutton", n_clicks=0, style={'margin':'30px 10px 10px 10px'}),
                    ],
                    className="rowBox"
                ),
                Download(id="downloadEQplots"),
                html.Div(id="hiddenDivSaveEQ")

            ],
            className="box",
        )


#==========MHD Callbacks

@app.callback([Output('timeSlider', 'min'),
               Output('timeSlider', 'max'),
               Output('timeSlider', 'marks'),
               Output('timeSlider', 'value'),
               Output('gFileTable2', 'data'),
               Output('MHDDataStorage', 'data')],
              [Input('loadMHD', 'n_clicks')],
              [State('shot', 'value'),
              State('tmin', 'value'),
              State('tmax', 'value'),
              State('nTrace', 'value'),
              State('gfiletable-upload', 'filename'),
              State('gfiletable-upload', 'contents'),
              State('plasma3Dmask', 'value'),
              State('dataPath', 'value'),
              State('MachFlag', 'value')]
              )
def loadMHD(n_clicks,shot,tmin,tmax,nTrace,gFileList,gFileData,plasma3Dmask,dataPath,MachFlag):
    """
    Load MHD
    """
    if MachFlag == None:
        print("Select a machine before loading MHD parameters")
        raise PreventUpdate

    #if data directory doesn't exist, create it
    tools.makeDir(dataPath, clobberFlag=False, mode=chmod, UID=UID, GID=GID)

    if plasma3Dmask == 'plasma3D': plasma3Dmask=1
    else:  plasma3Dmask=0

    if (shot is None) and (gFileList is None):
        raise PreventUpdate

    if shot is not None: shot = int(shot)
    if tmin is not None: tmin = int(tmin)
    if tmax is not None: tmax = int(tmax)
    if nTrace is not None: nTrace = int(nTrace)
    if gFileList is not None:
        if type(gFileList) is not list:
            gFileList = [gFileList]
    if gFileData is not None:
        if type(gFileData) is not list:
            gFileData = [gFileData]

    gui.getMHDInputs(shot=shot,
                     tmin=tmin,
                     tmax=tmax,
                     nTrace=nTrace,
                     gFileList=gFileList,
                     gFileData=gFileData,
                     plasma3Dmask=plasma3Dmask,
                    )

    ts = gui.MHD.timesteps
    tminMHD = ts.min()
    tmaxMHD = ts.max()
    if gFileList is None:
        tAll = np.linspace(int(tmin), int(tmax), (int(tmax)-int(tmin)+1))
        data = [dict([{'filename':'', 'timestep':''}][0])]
    else:
        tAll = ts
        keys=["filename"]
        interpData = pd.DataFrame(gFileList, columns=keys)
        interpData["timestep[ms]"] = ""
        data = interpData.to_dict('records')
        #interpData = [dict([{'filename':g, 'timestep':''} for g in gFileList])]
    marks = {}
    for t in ts:
        if t in tAll:
            marks.update({int(t):'{}'.format(t)})

    value = ts[0]


    MHDdata = {
        'Shot Number':shot,
        'Minimum Timestep [ms]':tmin,
        'Maximum Timesteo [ms]':tmax,
        'Trace Distance [degrees]':nTrace,
        '3D Plasma? [0=False]':plasma3Dmask,
        }

    return tminMHD, tmaxMHD, marks, value, data, MHDdata


@app.callback([Output('hiddenDivSaveEQ', 'children'),
               Output('downloadEQplots', 'data')],
              [Input('saveEQbutton', 'n_clicks')],
              [State('EQpixelX', 'value'),
               State('EQpixelY', 'value')])
def saveEQplots(n_clicks, x, y):
    if n_clicks < 1:
        raise PreventUpdate
    try: MachFlag = gui.MachFlag
    except:
        print("You didn't select a machine")
        raise PreventUpdate
    try: ts = gui.MHD.timesteps
    except:
        print('Please load MHD before saving EQ plots')
        return [html.Label("Load MHD First", style={'color':'#f51b60'})]

    #get png resolution
    if x == None or y == None:
        #default PV window size on toms computer,
        #more or less preserving NSTX aspect ratio
        x = 526
        y = 760
    else:
        x = float(x)
        y = float(y)

    shot = gui.MHD.shot
    #this import needs to be here (not at top of file) to prevent FreeCAD qt5
    #shared object conflict
    import GUIscripts.plot2DEQ as eqPlot
    allFiles = []
    for t in ts:
        print("Saving timestep {:5d}ms to PNG".format(t))
        log.info("Saving timestep {:5d}ms to PNG".format(t))
        fileName = gui.MHD.tmpDir + 'EQplot_{:05d}.png'.format(t)
        allFiles.append(fileName)
        idx = np.where(t==ts)[0][0]
        ep = gui.MHD.ep[idx]
        #write EQ plot using matplotlib
        plt = eqPlot.EQ2Dplot(ep,shot,t,MachFlag,height=y)
        plt.savefig(fileName, facecolor="#1b1f22")

        #Write Equilibrium plot using plotly
        #plot = plotly2DEQ.makePlotlyEQDiv(shot, t, MachFlag, ep)
        #plot.update_layout(
        #    title="{:05d}ms".format(t),
        #    xaxis_title="R [m]",
        #    yaxis_title="Z [m]",
        #    autosize=True,
        #    paper_bgcolor='#1b1f22',
        #    plot_bgcolor='#1b1f22',
        #    showlegend=False,
        #    font=dict(
    #   #         family="Courier New",
        #        size=22,
        #        color="#dcdce3"
        #    )
        #    )
        #plot.write_image(fileName, width=x, height=y)

    #now zip all these plots into a single file that the user may
    #download from GUI
    from zipfile import ZipFile
    from os.path import basename
    zipFile = gui.tmpDir + 'EQplots.zip'
    zipObj = ZipFile(zipFile, 'w')
    for f in allFiles:
        zipObj.write(f, basename(f))
    zipObj.close()


    return [html.Label("Saved EQs to file", style={'color':'#f5d142'}),
            send_file(zipFile)]

#Load gfile
@app.callback([Output('hiddenDivGfileUpload', 'children')],
              [Input('gfiletable-upload', 'filename')],
              [State('MachFlag', 'value')])
def gfileUpload(gFile, MachFlag):
    if MachFlag is None:
        raise PreventUpdate
    else:
        return [html.Label("Loaded gFile: "+gFile[0], style={'color':'#f5d142'})]

#==========CAD==========
def buildCADbox():
    """
    CAD input parameters
    """
    return html.Div(
            id="CADbox",
            draggable='yes',
            children=[
                html.H6("CAD Settings"),
                #html.Label(id="ROIgridResLabel", children="Heat Flux Resolution [mm]  "),
                #dcc.Input(id="ROIGridRes", className="textInput"),
                html.Label(id="gridResLabel", children="Intersect Resolution [mm]"),
                dcc.Input(id="gridRes", className="textInput"),
                html.Button("Load Res Settings", id="loadRes", style={'margin':'10px 10px 10px 10px'}),
                html.Div(id="hiddenDivCAD1"),
                html.Label(id="STPdropLabel", children="STP File Direct Upload:"),
                dcc.Upload(
                    className="PFCupload",
                    id='CAD-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select STP file')
                    ]),
                    style={
                        'width': '60%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed',
                        'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                        },
                    multiple=True,
                    ),
                html.Div(id="hiddenDivCAD2"),
            ],
            className="box",
        )

#Load res button connect
@app.callback([Output('hiddenDivCAD1', 'children'),
               Output('CADDataStorage', 'data')],
              [Input('loadRes', 'n_clicks')],
              [State('gridRes', 'value'),
               State('MachFlag', 'value')])
def loadRes(n_clicks, gridRes, MachFlag):
    if n_clicks is None:
        raise PreventUpdate
    if MachFlag is None:
        return [html.Label("Select a machine", style={'color':'#fc0313'})]
    gui.getCADResInputs(gridRes)

    CADdata = {
        'Intersect Max Edge Length [mm]':gridRes
        }

    return [html.Label("Loaded Resolution Settings", style={'color':'#f5d142'}), CADdata]

#Load CAD button connect
@app.callback([Output('hiddenDivCAD2', 'children')],
              [Input('CAD-upload', 'filename')],
              [State('CAD-upload', 'last_modified'),
               State('CAD-upload', 'contents'),
               State('MachFlag', 'value')])
def loadCAD(STPfile, ts, STPcontents, MachFlag):
    if STPfile is None:
        raise PreventUpdate
    if MachFlag is None:
        return [html.Label("Select a machine first", style={'color':'#fc0313'})]
    else:
        contents = STPcontents[0]
        content_type, content_string = contents.split(',')
        STPdata= base64.b64decode(content_string)
        gui.getCAD(STPfile=STPfile[0],STPdata=STPdata, ts=ts[0])
    return [html.Label("Loaded CAD: "+STPfile[0], style={'color':'#f5d142'})]

#==========HF==========
def buildHFbox():
    """
    Heat Flux input parameters
    """
    return html.Div(
            id="HFbox",
            children=[
                html.H6("HF Settings"),
                html.Label(id="hfModeLabel", children="Select a Heat Flux Profile"),
                dcc.Dropdown(
                id='hfMode',
                className="wideSelect",
                style={'backgroundColor': 'transparent', 'color':'transparent'},
                #style=dropdown_style,
                options=[
                    {'label': 'Gaussian Spreading', 'value': 'eich'},
                    {'label': 'Multi-Exponential', 'value': 'multiExp'},
                    {'label': 'Limiter', 'value': 'limiter'},
                    {'label': 'From qFiles', 'value': 'qFile'}
                    ],
                value='eich',
                ),
                html.Div(id='hfParameters',
                         children=[
                                    loadHFSettings(mode=None,hidden=True),
                                   ],
                        ),
                html.Br(),
                html.Button("Load HF", id="loadHF", style={'margin':'10px 10px 10px 10px'}),
                html.Div(id="hiddenDivHF"),

            ],
            className="HFbox",
        )


#Heat Flux Callbacks and Conditionals
@app.callback(Output('LRthreshDiv', 'children'),
              [Input('LRmask', 'value')])
def LRselector(mask):
    return [ loadLRsettings(mask, hidden=False) ]


def loadLRsettings(mask, hidden=False):
    if hidden==True or mask=='no':
        style={"display":"none"}
    else:
        style={}
    return html.Div(
                    children=
                        [
                            html.Label("Long Range Checking Power Threshold [MW]", style=style),
                            dcc.Input(id="LRthresh", className="textInput", style=style),
                        ],
                        )

@app.callback([Output('hfParameters', 'children')],
              [Input('hfMode', 'value')])
def hfParameters(mode):
    div = [loadHFSettings(mode=mode, hidden=False)]
    return [div]

def PsolInput(hidden=False):
    #do this hidden business so that we can always load defaults into these id's
    if hidden==True:
        className="hiddenBox"
    else:
        className="hfInput"


    row2 = html.Div(
                children=[
                    html.Div(
                        children=[
                            html.Label("Upper Inner Power Fraction", className="hfLabel"),
                            dcc.Input(id="fracUI", className="hfInput2"),
                            ],
                        className="colBox"
                    ),
                    html.Div(
                        children=[
                            html.Label("Upper Outer Power Fraction", className="hfLabel"),
                            dcc.Input(id="fracUO", className="hfInput2"),
                            ],
                        className="colBox"
                    ),
                    ],
                className="rowBox",
            )

    row3 = html.Div(
                children=[
                    html.Div(
                        children=[
                            html.Label("Lower Inner Power Fraction", className="hfLabel"),
                            dcc.Input(id="fracLI", className="hfInput2"),
                            ],
                        className="colBox"
                    ),
                    html.Div(
                        children=[
                            html.Label("Lower Outer Power Fraction", className="hfLabel"),
                            dcc.Input(id="fracLO", className="hfInput2"),
                            ],
                        className="colBox"
                    ),
                    ],
                className="rowBox",
            )


    return html.Div(
             className=className,
             children=[
                    html.Label("Power Injected [MW]"),
                    dcc.Input(id="Pinj", className="textInput"),
                    html.Label("Radiated Fraction of Injected Power"),
                    dcc.Input(id="coreRadFrac", className="textInput"),
                    row2,
                    row3,
                    html.Label("Long Range Intersection Checking?"),
                    dcc.RadioItems(
                                id="LRmask",
                                options=[
                                    {'label': 'Yes', 'value': 'yes'},
                                    {'label': 'No', 'value': 'no'},
                                        ],
                                        value='no'
                                ),

                    html.Div(id="LRthreshDiv",
                            children=[ loadLRsettings(mask='no', hidden=True), ]
                            ),
                        ],
                    )

def loadHFSettings(mode=None, hidden=False):
    #do this hidden business so that we can always load defaults into these id's
    #hideMask corresponds to the following parameters:
    #[eichProfile, commonRegion, privateRegion]
    hideMask = ['hiddenBox','hiddenBox','hiddenBox','hiddenBox']
    if mode=='eich':
        hideMask = ['hfInput','hiddenBox','hiddenBox','hiddenBox']
    elif mode=='limiter':
        hideMask = ['hiddenBox','hiddenBox','hfInput','hiddenBox'] #common flux region
    elif mode=='multiExp':
        hideMask = ['hiddenBox','hfInput','hiddenBox','hiddenBox'] #common + private flux region
    elif mode=='qFile':
        hideMask = ['hiddenBox','hiddenBox','hiddenBox','hfInput'] #common + private flux region
    if hidden==True or mode==None:
        hideMask=['hiddenBox','hiddenBox','hiddenBox','hiddenBox']


    #if we are using qFile, hide Psol inputs
    if mode == 'qFile':
        hidden=True

    return html.Div(
            children=[
                #gaussian spreading / eich
                html.Div(
                    #className=hideMask[0],
                    children=[
                        eichParameters(hideMask[0]),
                        ]
                ),
                #multiple exponentials
                html.Div(
                    #className=hideMask[1],
                    children=[
                        multiExpParameters(hideMask[1]),
                        ]
                ),
                #limiters
                html.Div(
                    #className=hideMask[2],
                    children=[
                        limiterParameters(hideMask[2]),
                        ]
                ),
                #from file
                html.Div(
                    #className=hideMask[2],
                    children=[
                        qFileParameters(hideMask[3]),
                        ]
                ),
                PsolInput(hidden),
                    ],
                    )


def eichParameters(className):
    row1 = html.Div(
        className='rowBox',
        children=[
            html.Div(
            className="colBox",
            children=[
                html.Label("Select Heat Flux Width source:", className="hfLabel"),
                dcc.Dropdown(
                id='eichlqCNMode',
                className="SelectorBoxInput",
                style={'backgroundColor': 'transparent', 'color':'transparent'},
                options=[
                    {'label': 'Eich #15', 'value': 'eich'},
                    {'label': 'User Defined', 'value': 'user'}
                    ],
                value='eich',
                ),
                ],
            ),
            html.Div(
            className="colBox",
            children=[
                html.Label("Select Gaussian Spreading source:", className="hfLabel"),
                dcc.Dropdown(
                id='eichSMode',
                className="SelectorBoxInput",
                style={'backgroundColor': 'transparent', 'color':'transparent'},
                #style=dropdown_style,
                options=[
                    {'label': 'From Makowski Scaling', 'value': 'makowski'},
                    {'label': 'User Defined', 'value': 'user'}
                    ],
                value='makowski',
                ),
                ],
            ),
            ])

    row2 = html.Div(
        className='rowBox',
        children=[
            html.Div(
            className="colBox",
            children=[
                html.Label("User Defined Heat Flux Width [mm]:", className="hfLabel"),
                dcc.Input(id="lqEich", className="hfInput2"),

                ],
            ),
            html.Div(
            className="colBox",
            children=[
                html.Label("User Defined Gaussian Spreading [mm]:", className="hfLabel"),
                dcc.Input(id="S", className="hfInput2"),
                ],
            ),
            ])

    row3 = html.Div(
        className="rowBox",
        children=[
            html.Div(
                className="colBox",
                children=[
                    html.Label("Background Heat Flux [MW/m^2]", className="hfLabel"),
                    dcc.Input(id="qBG", className="hfInput2"),
                ]),
            html.Div(
                className="colBox",
                children=[
                    html.Label("Greenwald Density Fraction", className="hfLabel"),
                    dcc.Input(id="fG", className="hfInput2", value=0.6),
                ]),
        ])


    div = html.Div(
        className=className,
        children=[
            row1,
            row2,
            row3
        ]
    )

    return div


def multiExpParameters(className):
    row1 = html.Div(
        className="rowBox",
        children = [
            html.Div(
                className="colBox",
                children=[
                    html.Label("Select Common Near Heat Flux Width source:"),
                    dcc.Dropdown(
                    id='multiExplqCNMode',
                    className="SelectorBoxInput",
                    style={'backgroundColor': 'transparent', 'color':'transparent'},
                    options=[
                        {'label': 'Eich #15', 'value': 'eich'},
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=None,
                    )
                ]),
            html.Div(
                className="colBox",
                children=[
                    html.Label("Select Common Far Heat Flux Width source:"),
                    dcc.Dropdown(
                    id='multiExplqCFMode',
                    className="SelectorBoxInput",
                    style={'backgroundColor': 'transparent', 'color':'transparent'},
                    options=[
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=None,
                    )
                ]),
            ]
    )
    row2 = html.Div(
        className="rowBox",
        children = [
            html.Div(
                className="colBox",
                children=[
                    html.Label("Select Private Near Heat Flux Width source:"),
                    dcc.Dropdown(
                    id='multiExplqPNMode',
                    className="SelectorBoxInput",
                    style={'backgroundColor': 'transparent', 'color':'transparent'},
                    options=[
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=None,
                    )
                ]),
            html.Div(
                className="colBox",
                children=[
                    html.Label("Select Private Far Heat Flux Width source:"),
                    dcc.Dropdown(
                    id='multiExplqPFMode',
                    className="SelectorBoxInput",
                    style={'backgroundColor': 'transparent', 'color':'transparent'},
                    options=[
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=None,
                    )
                ]),
            ]
    )
    row3 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Near Heat Flux Width [mm]"),
                        dcc.Input(id="lqCN", className="hfInput2"),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Near Power Fraction"),
                        dcc.Input(id="fracCN", className="hfInput2"),
                    ]),
            ])

    row4 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Far Heat Flux Width [mm]"),
                        dcc.Input(id="lqCF", className="hfInput2"),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Far Power Fraction"),
                        dcc.Input(id="fracCF", className="hfInput2"),
                    ]),
            ])

    row5 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Private Near Heat Flux Width [mm]"),
                        dcc.Input(id="lqPN", className="hfInput2"),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Private Near Power Fraction"),
                        dcc.Input(id="fracPN", className="hfInput2"),
                    ]),
            ])

    row6 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Private Far Heat Flux Width [mm]"),
                        dcc.Input(id="lqPF", className="hfInput2"),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Private Far Power Fraction"),
                        dcc.Input(id="fracPF", className="hfInput2"),
                    ]),
            ])



    div = html.Div(
        className=className,
        children=[
            row1, #commented for now because only user defined is allowed (no regression)
            row2,
            row3,
            row4,
            row5,
            row6
        ]
    )

    return div


def limiterParameters(className):
    #return if this div is supposed to be hidden to prevent duplicate IDs
    #if className== 'hiddenBox':
    #    return

    row1 = html.Div(
        className="rowBox",
        children = [
            html.Div(
                className="colBox",
                children=[
                    html.Label("Select Common Near Heat Flux Width source:"),
                    dcc.Dropdown(
                    id='limiterlqCNMode',
                    className="SelectorBoxInput",
                    style={'backgroundColor': 'transparent', 'color':'transparent'},
                    options=[
                        {'label': 'Eich #15', 'value': 'eich'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value='eich',
                    )
                ]),
            html.Div(
                className="colBox",
                children=[
                    html.Label("Select Common Far Heat Flux Width source:"),
                    dcc.Dropdown(
                    id='limiterlqCFMode',
                    className="SelectorBoxInput",
                    style={'backgroundColor': 'transparent', 'color':'transparent'},
                    options=[
                        {'label': 'Horacek Fig. 6a', 'value': 'horacek'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value='user',
                    )
                ]),

            ]
            )

    row2 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Near Heat Flux Width [mm]"),
                        dcc.Input(id="limlqCN", className="hfInput2"),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Far Heat Flux Width [mm]"),
                        dcc.Input(id="limlqCF", className="hfInput2"),
                    ]),
            ])

    row3 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Near Power Fraction"),
                        dcc.Input(id="limfracCN", className="hfInput2"),
                    ]),

                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Far Power Fraction"),
                        dcc.Input(id="limfracCF", className="hfInput2"),
                    ]),
            ])

    div = html.Div(
        className=className,
        children=[
            row1,
            row2,
            row3
            ]
        )
    return div


def  qFileParameters(className):
    #return if this div is supposed to be hidden to prevent duplicate IDs
    #if className== 'hiddenBox':
    #    return

    row1 = html.Div(
        className="rowBox",
        children = [
            html.Div(
                className="colBox",
                children=[
                    html.Label("Input qFile shot path:"),
                    dcc.Input(id="qFilePath", className="hfInput2"),
                ]),
            html.Div(
                className="colBox",
                children=[
                    html.Label("Input qFile tag (ie HF_optical.csv):"),
                    dcc.Input(id="qFileTag", className="hfInput2"),
                ]),

            ]
            )

    div = html.Div(
        className=className,
        children=[
            row1,
            ]
        )
    return div






def commonRegionParameters():
    """
    near and far heat flux widths and power sharing fractions
    """
    row1 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Near Heat Flux Width [mm]"),
                        dcc.Input(id="lqCN", className="hfInput2"),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Near Power Fraction"),
                        dcc.Input(id="fracCN", className="hfInput2"),
                    ]),
            ])

    row2 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Far Heat Flux Width [mm]"),
                        dcc.Input(id="lqCF", className="hfInput2"),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        html.Label("Common Far Power Fraction"),
                        dcc.Input(id="fracCF", className="hfInput2"),
                    ]),
            ])

    return row1,row2


#Load HF button connect
@app.callback([Output('hiddenDivHF', 'children'),
               Output('HFDataStorage', 'data')],
              [Input('loadHF', 'n_clicks')],
              [State('hfMode', 'value'),
               State('MachFlag', 'value'),
               State('lqEich', 'value'),
               State('S', 'value'),
               State('lqCN', 'value'),
               State('lqCF', 'value'),
               State('lqPN', 'value'),
               State('lqPF', 'value'),
               State('fracCN', 'value'),
               State('fracCF', 'value'),
               State('fracPN', 'value'),
               State('fracPF', 'value'),
               State('fracUI', 'value'),
               State('fracUO', 'value'),
               State('fracLI', 'value'),
               State('fracLO', 'value'),
               State('LRmask', 'value'),
               State('LRthresh', 'value'),
               State('eichlqCNMode', 'value'),
               State('eichSMode', 'value'),
               State('multiExplqCNMode', 'value'),
               State('multiExplqCFMode', 'value'),
               State('multiExplqPNMode', 'value'),
               State('multiExplqPFMode', 'value'),
               State('limiterlqCNMode', 'value'),
               State('limiterlqCFMode', 'value'),
               State('limlqCN', 'value'),
               State('limlqCF', 'value'),
               State('limfracCN', 'value'),
               State('limfracCF', 'value'),
               State('qBG', 'value'),
               State('Pinj', 'value'),
               State('coreRadFrac', 'value'),
               State('fG', 'value'),
               State('qFilePath', 'value'),
               State('qFileTag', 'value'),
               ])
def loadHF(n_clicks,hfMode,MachFlag,
            lqEich,S,lqCN,lqCF,lqPN,lqPF,
            fracCN,fracCF,fracPN,fracPF,
            fracUI,fracUO,fracLI,fracLO,
            LRmask,LRthresh,
            eichlqCNMode,SMode,
            multiExplqCNMode,multiExplqCFMode,multiExplqPNMode,multiExplqPFMode,
            limiterlqCNMode,limiterlqCFMode,limlqCN,limlqCF,limfracCN,limfracCF,
            qBG,Pinj,coreRadFrac,fG,
            qFilePath, qFileTag):
    if MachFlag is None:
        raise PreventUpdate
    else:
        #set up the heat flux configuration (which scalings to use)
        if hfMode == 'limiter':
            lqCNmode = limiterlqCNMode
            lqCFmode = limiterlqCFMode
            lqPNmode = None
            lqPFmode = None
            SMode = None
            lqCN = limlqCN
            lqCF = limlqCF
            lqPN = 0.0
            lqPF = 0.0
            fracCN = limfracCN
            fracCF = limfracCF
            fracPN = 0.0
            fracPF = 0.0
            qFileTag = None
            qFilePath = None
        elif hfMode == 'multiExp':
            #add this back in after brunner scaling is in HEAT:
            lqCNmode = multiExplqCNMode
            lqCFmode = multiExplqCFMode
            lqPNmode = multiExplqPNMode
            lqPFmode = multiExplqPFMode
            lqCN = lqCN #unnecessary, but here for clarity
            lqCF = lqCF #unnecessary, but here for clarity
            lqPN = lqPN #unnecessary, but here for clarity
            lqPF = lqCF #unnecessary, but here for clarity
            SMode = None
            qFileTag = None
            qFilePath = None
        elif hfMode == 'qFile':
            lqCNmode = eichlqCNMode
            lqCFmode = None
            lqPNmode = None
            lqPFmode = None
            lqCN = lqEich
            lqCF = 0.0
            lqPN = 0.0
            lqPF = 0.0
            SMode = SMode
        else: #eich mode is default
            lqCNmode = eichlqCNMode
            lqCFmode = None
            lqPNmode = None
            lqPFmode = None
            lqCN = lqEich
            lqCF = 0.0
            lqPN = 0.0
            lqPF = 0.0
            SMode = SMode
            qFileTag = None
            qFilePath = None
        #could add private flux scalings here if they ever exist


        gui.getHFInputs(hfMode,LRmask,LRthresh,
                        lqCN,lqCF,lqPN,lqPF,S,
                        fracCN,fracCF,fracPN,fracPF,
                        fracUI,fracUO,fracLI,fracLO,
                        lqCNmode,lqCFmode,lqPNmode,lqPFmode,SMode,
                        qBG,Pinj,coreRadFrac,fG,
                        qFilePath,qFileTag)


        #Update output tab table
        HFdata = gui.HF.HFdataDict
        #hfDict = [{'Parameter':i, 'Value':dataDict[i]} for i in inputTableData]


    return [html.Label("Loaded HF Settings", style={'color':'#f5d142'}), HFdata]



#==========PFC==========
def buildPFCbox():
    return html.Div(
            id="PFCbox",
            children=[
                html.H6("PFC Settings"),
                html.Div(id="pfcUploadDiv"),
                dcc.Upload(
                    className="PFCupload",
                    id='pfctable-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select Files')
                    ]),
                    style={
                        'width': '60%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed',
                        'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                        },
                    ),
                html.Div(children=loadPFCtable(), className="PFCtable"),
                html.Br(),
                html.Button("Load PFC Settings", id="loadPFC", n_clicks=0, style={'margin':'0 10px 10px 0'}),
#                html.Button('Add Row', id='add-rows-button', n_clicks=0, style={'margin':'0 10px 10px 0'}),
                html.Button("Download Default PFC file", id="downloadPFCbtn", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                Download(id="downloadPFC"),
                html.Div(id="hiddenDivDownloadPFC", style={"display": "none"}),
                html.Div(id="hiddenDivPFC", style={"display": "none"}),
                dcc.Input(id="hiddenDivPFC2", style={"display": "none"}),
                html.Div(id="hiddenDivPFC3")

            ],
            className="PFCbox",
            )

def loadPFCtable():
    params = ['timesteps','PFCname','resolution','DivCode','intersectName','excludeName']
    cols = [{'id': p, 'name': p} for p in params]
    data = [{}]
    return dash_table.DataTable(
        id='pfcTable',
        columns = cols,
        data = data,
        style_header={'backgroundColor': 'rgb(30, 30, 30)'},
        style_cell={
            'textAlign': 'left',
            'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'white'
        },
        editable=True,
        export_format='csv',
        row_deletable=True,
        )

#upload file area
def parse_contents(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    if 'csv' in filename:
        # Assume that the user uploaded a CSV file
        return pd.read_csv(io.StringIO(decoded.decode('utf-8')), comment='#')
    elif 'xls' in filename:
        # Assume that the user uploaded an excel file
        return pd.read_excel(io.BytesIO(decoded))


#Uploading files and button presses update table and HTML storage object
@app.callback([Output('PFCdataStorage', 'data'),
               Output('pfcTable', 'data'),
               Output('pfcTable', 'columns'),
               Output('hiddenDivPFC3', 'children'),
               Output('gyroSource', 'options')],
              [Input('loadPFC', 'n_clicks'),
               Input('pfctable-upload', 'filename')],
              [State('PFCdataStorage', 'data'),
               State('PFCdataStorage', 'modified_timestamp'),
               State('pfctable-upload', 'contents'),
               State('pfcTable', 'data'),
               State('pfcTable', 'columns'),
               State('MachFlag', 'value')])
def PFCtable(n_clicks, filename, dataStore, ts, uploadContents,
             tableData, tableColumns, MachFlag):
    tD = [{}]
    if dataStore==None:
        print('Initializing HTML PFC data storage object')
        dataStore = {}
        dataStore.update({'PFC_n_clicks':0})
        dataStore.update({'PFCfilename':''})

    #user has to load MHD and CAD for a specific machine before PFCs
    if MachFlag == None:
        hiddenDiv = [html.Label("Select a Tokamak and Load MHD/CAD First", style={'color':'#db362a'})]

    #Button Clicks
    elif n_clicks > 0 and n_clicks>dataStore.get('PFC_n_clicks'):
        dataStore['PFC_n_clicks'] = n_clicks
        #read default machine PFC file
        if (tableData == [{}]) and (uploadContents is None):
            gui.getPFCinputs(defaultMask=True)
            tableData = gui.timestepMap.to_dict('records')
            tableColumns = [{"name": i, "id": i} for i in gui.timestepMap.columns]
            tD = tableData.copy()

        #load PFC data from page and send it to HEAT
        else:
            gui.getPFCdataFromGUI(tableData)
            gui.getPFCinputs(defaultMask=False)
            tD = gui.timestepMap.to_dict('records')

        hiddenDiv = [html.Label("Loaded PFC Data into HEAT", style={'color':'#f5d142'})]

    #file dropper
    elif filename != dataStore['PFCfilename']:
        df = parse_contents(uploadContents, filename)
        tableData = df.to_dict('records')
        tableColumns = [{"name": i, "id": i} for i in df.columns]
        hiddenDiv = [html.Label("Loaded file: "+filename, style={'color':'#34b3ed'})]
        tD = tableData.copy()

    else:
        hiddenDiv = []

    #build dictionary for gyro orbit source dropdown
    ROIdict = []
    #ROIdict.append({'label': "", 'value': ''})
    if tD != [{}]:
        tD = pd.DataFrame.from_dict(tableData)[list (tableData[0].keys())]
        names = tD.rename(columns=lambda x: x.strip())['PFCname'].values
        for name in names:
            n = name.replace(' ', '')
            ROIdict.append({'label': n, 'value': n})

    return dataStore, tableData, tableColumns, hiddenDiv, ROIdict

#Download PFC Default file button connect
@app.callback([Output('downloadPFC', 'data'),
               Output('hiddenDivDownloadPFC', 'children')],
              [Input('downloadPFCbtn', 'n_clicks')])
def downloadPFCfile(n_clicks):
    if n_clicks < 1:
        raise PreventUpdate
    gui.savePFCfile()
    return [send_file(gui.tmpDir + "PFCinput.csv"),
            html.Label("Saved PFC Default File", style={'color':'#f5d142'})]






#==========gyro orbits==========
def buildGYRObox():
    return html.Div(
        id="gyrobox",
        children=[
            html.H6("Gyro Orbit Settings"),
            gyroInputBoxes(),
            html.Label("Velocity Method"),
            dcc.Dropdown(
                id='vMode',
                className="wideSelect",
                style={'backgroundColor': 'transparent', 'color':'transparent'},
                options=[
                    {'label': 'Single Value', 'value': 'single'},
                    {'label': 'From 3D File', 'value': 'file'}, # to be added at future date
                    ],
                value='single',
                ),
            html.Div(id='gyroVparams',
                     children=[
                                loadGyroSettings(mode=None,hidden=True),
                               ],
                    ),
            html.Label("Gyro Orbit Power Source"),
            dcc.Dropdown(
                id='sourceMode',
                className="wideSelect",
                style={'backgroundColor': 'transparent', 'color':'transparent'},
                options=[
                    {'label': 'All ROI PFCs', 'value': 'allROI'},
                    {'label': 'User Selected PFCs', 'value': 'custom'},
                    ],
                value=None,
                ),
            html.Div(id='gyroSourceParams',
                     children=[
                                loadSourceSettings(mode=None,hidden=True),
                               ],
                    ),
            html.Label("(Note: Sources are OMITTED from Temp Calculation)"),
            html.Br(),
            html.Button("Load Gyro Settings", id="loadGYRO", n_clicks=0, style={'margin':'0 10px 10px 0'}),
            html.Div(id="hiddenDivGyro")
        ],
        className="HFbox",
    )


def gyroInputBoxes():
    return html.Div(
            children=[
            html.Div(
                children=[
                    html.Label("# steps per helix"),
                    dcc.Input(id="N_gyroSteps", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label("# gyroPhase Angles"),
                    dcc.Input(id="N_gyroPhase", className="textInput"),
                ],
                className="OFInput"
            ),

            html.Div(
                children=[
                    html.Label("Gyro Trace Length [deg]"),
                    dcc.Input(id="gyroDeg", className="textInput"),
                ],
                className="OFInput"
            ),
            html.Div(
                children=[
                    html.Label("Ion Power Fraction [0-1]"),
                    dcc.Input(id="ionFrac", className="textInput"),
                ],
                className="OFInput"
            ),
            html.Div(
                children=[
                    html.Label("Ion Mass [AMU]"),
                    dcc.Input(id="ionMassAMU", className="textInput"),
                ],
                className="OFInput"
            ),
            ],
            className="wideBoxNoColor",
            )

@app.callback([Output('gyroVparams', 'children')],
              [Input('vMode', 'value')])
def hfParameters(mode):
    #select velocity mode
    div = [loadGyroSettings(mode=mode, hidden=False)]
    return [div]

def loadGyroSettings(mode=None, hidden=False):
    #do this hidden business so that we can always load defaults into these id's
    #hideMask corresponds to the following parameters:
    #[eichProfile, commonRegion, privateRegion]
    hideMask = ['hiddenBox','hiddenBox']
    if mode=='single':
        hideMask = ['hfInput','hiddenBox']
    elif mode=='file':
        hideMask = ['hiddenBox','hfInput'] #common flux region
    if hidden==True or mode==None:
        hideMask=['hiddenBox','hiddenBox']

    return html.Div(
            children=[
                #use single temperature to define entire divertor
                html.Div(
                    #className=hideMask[0],
                    children=[
                        singleVelocity(hideMask[0]),
                        ]
                ),
                #read temperature pointcloud from file and interpolate to PFC surface
                html.Div(
                    #className=hideMask[1],
                    children=[
                        velocityFromFile(hideMask[1]),
                        ]
                ),
                    ],
                    )


def singleVelocity(className):
    return html.Div(
            children = [
                html.Div(
                    children=[
                    html.Div(
                        children=[
                            html.Label("Average Temperature [eV]"),
                            dcc.Input(id="gyroT_eV", className="textInput"),
                            ],
                        className="OFInput"
                        ),

                    html.Div(
                        children=[
                            html.Label("# Velocity Slices"),
                            dcc.Input(id="N_vSlice", className="textInput"),
                            ],
                        className="OFInput"
                        ),

                    html.Div(
                        children=[
                            html.Label("# Velocity Phases"),
                            dcc.Input(id="N_vPhase", className="textInput"),
                            ],
                        className="OFInput"
                        ),


                        ],
                    className="wideBoxNoColor"
                    )
                ],
            className=className,
            )

def velocityFromFile(className):
    return html.Div(
        children=[
            html.Label("3D plasma temperature interpolation not yet available", style={'color':'#f5d142'})
        ],
        className = className
    )





@app.callback([Output('gyroSourceParams', 'children')],
              [Input('sourceMode', 'value')],
              [State('gyroSource', 'options'),
               State('loadPFC', 'n_clicks')]
              )
def gyroParameters(mode, options, n_clicks):
    if n_clicks == 0:
        raise PreventUpdate
    #select gyroSource mode
    div = [loadSourceSettings(mode=mode, hidden=False, options=options)]
    return [div]

def loadSourceSettings(mode=None, hidden=False, options=None):
    #do this hidden business so that we can always load defaults into these id's
    #hideMask corresponds to the following parameters:
    #[eichProfile, commonRegion, privateRegion]
    hideMask = ['hiddenBox','hiddenBox']
    if mode=='allROI':
        hideMask = ['hfInput','hiddenBox']
    elif mode=='custom':
        hideMask = ['hiddenBox','hfInput']
    if hidden==True or mode==None:
        hideMask=['hiddenBox','hiddenBox']

    return html.Div(
            children=[
                #use single temperature to define entire divertor
                html.Div(
                    #className=hideMask[0],
                    children=[
                        allROISource(hideMask[0]),
                        ]
                ),
                #read temperature pointcloud from file and interpolate to PFC surface
                html.Div(
                    #className=hideMask[1],
                    children=[
                        customSource(hideMask[1], options),
                        ]
                ),
                    ],
                    )


def customSource(className, options):
    if options == None: options = [{'label':'', 'value':''}]
    div = html.Div(
            children = [
                html.Div(
                    children=[
                        dcc.Checklist(
                            id='gyroSource',
                            className="checkListBox",
                            options=options,
                            value=['allROI'],
                        ),
                    ],
                    className="wideBoxNoColor"
                    )
                ],
            className=className,
            )
    return div

def allROISource(className):
    return html.Div(
        children=[
            html.Label("Using all ROI PFCs as gyro power sources", style={'color':'#f5d142'})
        ],
        className = className
    )





#Load GYRO button connect
@app.callback([Output('hiddenDivGyro', 'children'),
               Output('gyroPhasePlot', 'children'),
               Output('vPhasePlot', 'children'),
               Output('vSlicePlot', 'children'),
               Output('cdfSlicePlot', 'children'),
               Output('GYRODataStorage', 'data')
               ],
              [Input('loadGYRO', 'n_clicks')],
              [State('N_gyroSteps', 'value'),
               State('N_gyroPhase', 'value'),
               State('gyroDeg', 'value'),
               State('ionMassAMU', 'value'),
               State('vMode','value'),
               State('gyroT_eV', 'value'),
               State('N_vPhase', 'value'),
               State('N_vSlice', 'value'),
               State('ionFrac', 'value'),
               State('gyroSource', 'value')
              ])
def loadGYRO(n_clicks,N_gyroSteps,N_gyroPhase,gyroDeg,ionMassAMU,vMode,gyroT_eV,
             N_vPhase, N_vSlice, ionFrac, gyroSource):
    """
    sets up GYRO module
    """
    if n_clicks == 0:
        raise PreventUpdate

    #handle gyro orbit dropdown (user can only select all or custom PFCs)
    if ('allROI' in gyroSource) and (len(gyroSource)>1):
        gyroSource.remove('allROI')

    gui.getGyroInputs(N_gyroSteps,N_gyroPhase,gyroDeg,ionMassAMU,vMode,gyroT_eV,
                      N_vPhase, N_vSlice, ionFrac, gyroSource)
    gyroPhaseFig = gyroPhasePlots(update=True)
    vPhaseFig = vPhasePlots(update=True)
    vSliceFig = vSlicePlots(update=True)
    cdfSliceFig = cdfSlicePlots(update=True)



    GYROdata = {
            'Number of steps per helix period':N_gyroSteps,
            'Number of samples in gyro phase space':N_gyroPhase,
            'Gyro trace length [degrees] (goes both directions)': gyroDeg,
            'Ion effective mass [AMU]': ionMassAMU,
            'Velocity / Temperature Mode':vMode,
            'Ion temperature at PFC surface [eV]':gyroT_eV,
            'Number of samples in velocity phase space':N_vPhase,
            'Number of samples from Maxwellian velocity distribution':N_vSlice,
            'Fraction of power carried by ions':ionFrac,
            'Source for gyro orbit power':':'.join(gyroSource),
            }

    return [html.Label("Loaded Gyro Orbit Data into HEAT", style={'color':'#f5d142'}),
            gyroPhaseFig,
            vPhaseFig,
            vSliceFig,
            cdfSliceFig,
            GYROdata]



#==========openFOAM==========
def buildOFbox():
    return html.Div(
        id="OFbox",
        children=[
            html.H6("openFOAM Settings"),
            OFinputBoxes(),
            html.Br(),
            html.Button("Load OF Settings", id="loadOF", n_clicks=0, style={'margin':'0 10px 10px 0'}),
            html.Div(id="hiddenDivOF")
        ],
        className="HFbox",
    )

def OFinputBoxes():
    return html.Div(
            children=[
            html.Div(
                children=[
                    html.Label("Start Time [ms]"),
                    dcc.Input(id="OFstartTime", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label("Stop Time [ms]"),
                    dcc.Input(id="OFstopTime", className="textInput"),
                ],
                className="OFInput"
            ),
            html.Div(
                children=[
                    html.Label("Minimum Mesh Refinement Level"),
                    dcc.Input(id="OFminMeshLev", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label("Maximum Mesh Refinement Level"),
                    dcc.Input(id="OFmaxMeshLev", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label("STL scaling"),
                    dcc.Input(id="OFSTLscale", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label("deltaT [s]"),
                    dcc.Input(id="OFdeltaT", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label(" Write deltaT [s]"),
                    dcc.Input(id="OFwriteDeltaT", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label("Material Selection"),

                    dcc.Dropdown(
                        id='materialSelect',
                        className="machineSelect",
                        style={'backgroundColor': 'transparent', 'color':'transparent',
                                'align-items':'center'},
                                #style=dropdown_style,
                                options=[
                                {'label': 'SGLR6510 Graphite', 'value': 'SGL'},
                                {'label': 'ATJ Graphite', 'value': 'ATJ'},
                                {'label': 'Molybdenum', 'value': 'MOLY'},
                                {'label': 'Tungsten', 'value': 'TUNG'},

                                ],
                                value='SGL'
                                ),
                ],
                className="OFInput",
            ),
            ],
            className="wideBoxNoColor",
            )

#Load OF button connect
@app.callback([Output('hiddenDivOF', 'children'),
               Output('OFDataStorage', 'data')],
              [Input('loadOF', 'n_clicks')],
              [State('OFstartTime', 'value'),
               State('OFstopTime', 'value'),
               State('OFminMeshLev', 'value'),
               State('OFmaxMeshLev', 'value'),
               State('OFSTLscale', 'value'),
               State('OFdeltaT', 'value'),
               State('OFwriteDeltaT', 'value'),
               State('OFbashrc', 'value'),
               State('materialSelect', 'value')
              ])
def loadOF(n_clicks,OFstartTime,OFstopTime,
            OFminMeshLev,OFmaxMeshLev,OFSTLscale,OFdeltaT,OFwriteDeltaT,
            OFbashrcLoc,materialSelect):
    """
    sets up openFOAM for an analysis
    """
    if n_clicks < 1:
        raise PreventUpdate
    else:
        gui.loadOF(OFstartTime,OFstopTime,
                   OFminMeshLev,OFmaxMeshLev,
                   OFSTLscale,OFbashrcLoc,OFdeltaT,OFwriteDeltaT,materialSelect)
        OFdata = {
            'OpenFOAM start time [ms]':OFstartTime,
            'OpenFOAM stop time [ms]':OFstopTime,
            'Minimum FVM mesh refinement level':OFminMeshLev,
            'Maximum FVM mesh refinement level':OFmaxMeshLev,
            'Scale multiplier for FVM meshes':OFSTLscale,
            'OpenFOAM simulation timestep size [s]': OFdeltaT,
            'OpenFOAM write output timestep size [s]':OFwriteDeltaT,
            'OpenFOAM material selection': materialSelect
            }
        return [html.Label("Loaded OF Data into HEAT", style={'color':'#f5d142'}), OFdata]





"""
==============================================================================
Tab Contents: run tab
"""
def buildRunTab():
    return html.Div(
            id="runTab",
            children=[runChildren()],
            className="runContentBox",
            )

def runChildren():
    return html.Div(
            children = [
                html.H4("HEAT Run Settings", className="buttonRibbon"),
                html.Br(),
                html.H6("Acceleration, Filters, Settings: "),
                buildRunSettings(),
                html.H6("Point Clouds at Tile Surface:"),
                runTabChecklist(),
                html.H6("Traces from Points:"),
                runTabTraces(),
                html.Button("Run HEAT", id="runHEAT", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                html.Div(id="hiddenDivRun")
                ],
            className="wideBoxDark",
                )


#==========HEAT Run Settings==========
def buildRunSettings():
    """
    returns user options for the HEAT run
    """
    return html.Div(
            id="runSettings",
            className="buttonRibbon",
            children=[
                html.Label(children="Acceleration Filters:  "),
                dcc.Checklist(
                    options=[
                        {'label': 'Toroidal Filter', 'value': 'torFilt'},
                        {'label': 'Psi Filter', 'value': 'psiFilt'},
                    ],
                    value=['torFilt'],
                    id='accFilters',
                    className="PCbox",
                )
            ],
        )



def runTabChecklist():
    return html.Div(
            id="runTabChecklist",
            children=[
                dcc.Checklist(
                    options=[
                        {'label': 'B-field point cloud ', 'value': 'B'},
                        {'label': 'Normal vector point cloud', 'value': 'norm'},
                        {'label': 'powerDir point cloud', 'value': 'pwrDir'},
                        {'label': 'psiN point cloud', 'value': 'psiN'},
                        {'label': 'bdotn point cloud', 'value': 'bdotn'},
                        {'label': 'Heat flux point cloud', 'value': 'hfOpt'},
                        {'label': 'Gyro Orbit heat flux point cloud', 'value': 'hfGyro'},
                        {'label': 'openFOAM thermal analysis', 'value': 'T'},
                        ],
                        value=['hfOpt'],
                        id='checklistPC',
                        ),
                    ],
                className="PCbox",
                )
def runTabTraces():
    return html.Div(
            id="runTabTraces",
            children=[
                dcc.Checklist(
                    options=[{'label': 'B-field Trace ', 'value': 'Btrace'}],
                    value=[''],
                    id="Btrace",
                        ),
                html.Div(id="bFieldTracePoints", children=[loadBfieldTrace(False)]),
                dcc.Checklist(
                    options=[{'label': 'OpenFOAM Temp Probe ', 'value': 'OFtrace'}],
                    value=[''],
                    id="OFtrace",
                        ),
                html.Div(id="OFTracePoints", children=[loadOFTrace(False)]),
                dcc.Checklist(
                    options=[{'label': 'Gyro Orbit Trace ', 'value': 'gyrotrace'}],
                    value=[''],
                    id="gyrotrace",
                        ),
                html.Div(id="gyroTracePoints", children=[loadGyroTrace(False)]),
                    ],
                className="PCbox",
                )

@app.callback([Output('bFieldTracePoints', 'children')],
              [Input('Btrace', 'value')])
def bfieldTracePoint(value):
    if 'Btrace' in value:
        val = True
    else:
        val = False
    return [loadBfieldTrace(val)]


#this function enables the inputs to be rendered on page load but hidden
def loadBfieldTrace(display):
    if (display == True):
        style={}
    else:
        style={"display":"none"}
    return html.Div(
        children=[
            dcc.Upload(
                className="PFCupload",
                id='Btrace-upload',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select File')
                    ]),
                    style={
                        'width': '100%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed',
                        'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                        },
                multiple=False,
                ),
            html.Div(children=loadBtraceTable(), className="PFCtable"),
            html.Div(id="hiddenDivBtrace"),
        ],
        className="xyzBoxVert",
        style=style,
    )

def loadBtraceTable():
    params = ['x[mm]','y[mm]','z[mm]','traceDirection','Length[deg]']
    cols = [{'id': p, 'name': p} for p in params]
    data = [{}]
    return dash_table.DataTable(
        id='BtraceTable',
        columns = cols,
        data = data,
        style_header={'backgroundColor': 'rgb(30, 30, 30)'},
        style_cell={
            'textAlign': 'left',
            'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'white'
        },
        editable=True,
        export_format='csv',
        row_deletable=True,
        persistence=True,
        )

#Uploading files and button presses update table and HTML storage object
@app.callback([Output('BtraceDataStorage', 'data'),
               Output('BtraceTable', 'data'),
               Output('BtraceTable', 'columns'),
               Output('hiddenDivBtrace', 'children')],
              [Input('Btrace-upload', 'filename')],
              [State('BtraceDataStorage', 'data'),
               State('Btrace-upload', 'contents'),
               State('BtraceTable', 'data'),
               State('BtraceTable', 'columns'),
               State('MachFlag', 'value'),
               State('Btrace', 'value')])
def BtraceTable(filename, dataStore, uploadContents,
             tableData, tableColumns, MachFlag, BtraceList):
    if 'Btrace' in BtraceList:
        trace = True
    else:
        trace = False

    if dataStore==None:
        print('Initializing HTML PFC data storage object')
        dataStore = {}
        dataStore.update({'BtraceFileName':None})
        dataStore.update({'BtraceContents':None})
        dataStore.update({'Btrace':trace})

    #user has to load MHD and CAD for a specific machine before PFCs
    if MachFlag == None:
        hiddenDiv = [html.Label("Select a Tokamak and Load MHD First", style={'color':'#db362a'})]
    elif filename is None:
        hiddenDiv = []
    #file dropper
    elif (trace == True) and (uploadContents!= dataStore['BtraceContents']):
            df = parse_contents(uploadContents, filename)
            tableData = df.to_dict('records')
            #tableColumns = [{"name": i, "id": i} for i in df.columns]
            params = ['x[mm]','y[mm]','z[mm]','traceDirection','Length[deg]']
            tableColumns = [{'id': p, 'name': p} for p in params]
            hiddenDiv = [html.Label("Loaded file: "+filename, style={'color':'#34b3ed'})]
            dataStore.update({'BtraceFileName':filename})
            dataStore.update({'BtraceContents':uploadContents})
    else:
        hiddenDiv = []

    return dataStore, tableData, tableColumns, hiddenDiv



@app.callback(Output('OFTracePoints', 'children'),
              [Input('OFtrace', 'value')])
def OFTracePoint(value):
    if 'OFtrace' in value:
        val = True
    else:
        val = False
    return [loadOFTrace(val)]


#this function enables the inputs to be rendered on page load but hidden
def loadOFTrace(display):
    if (display == True):
        style={}
    else:
        style={"display":"none"}
    return html.Div(
                children=[
                    html.Div(
                        children = [
                            html.Label("x [mm]"),
                            dcc.Input(id="xOFtrace", className="xyzBoxInput"),
                            html.Label("y [mm]"),
                            dcc.Input(id="yOFtrace", className="xyzBoxInput"),
                            html.Label("z [mm]"),
                            dcc.Input(id="zOFtrace", className="xyzBoxInput"),
                            #html.Label("t [ms]"),
                            #dcc.Input(id="tOFtrace", className="xyzBoxInput"),
                        ],
                        className="xyzBox"
                        )
                    ],
                style=style,
                className="xyzBoxVert",
                    )

@app.callback(Output('gyroTracePoints', 'children'),
              [Input('gyrotrace', 'value')])
def gyroTracePoints(value):
    if 'gyrotrace' in value:
        val = True
    else:
        val = False
    return [loadGyroTrace(val)]

#this function enables the inputs to be rendered on page load but hidden
def loadGyroTrace(display):
    if (display == True):
        style={}
    else:
        style={"display":"none"}
    return html.Div(
                children=[
                    html.Div(
                        children=[
                            html.Label("x [mm]"),
                            dcc.Input(id="xGyroTrace", className="xyzBoxInput"),
                            html.Label("y [mm]"),
                            dcc.Input(id="yGyroTrace", className="xyzBoxInput"),
                            html.Label("z [mm]"),
                            dcc.Input(id="zGyroTrace", className="xyzBoxInput"),
                        ],
                        className="xyzBox"
                    ),
                    html.Div(
                        children=[
                            html.Label("Temperature [eV]"),
                            dcc.Input(id="gyroT_eV_trace", className="xyzBoxInput"),
                            html.Label("Gyro Phase [deg]"),
                            dcc.Input(id="gyroPhase_trace", className="xyzBoxInput"),
                            html.Label("Trace Length [deg]"),
                            dcc.Input(id="gyroDeg_trace", className="xyzBoxInput"),
                            html.Label("# Helix Discretizations"),
                            dcc.Input(id="N_gyroSteps_trace", className="xyzBoxInput"),
                            html.Label("Trace Direction"),
                            dcc.Input(id="gyroDir_trace", className="xyzBoxInput"),
                        ],
                        className="xyzBox"
                    )
                    ],
                style=style,
                className="xyzBoxVert",
                    )



@app.callback([Output('hiddenDivRun', 'children'),
               Output('qDivDist', 'children'),
               Output('OFmaxTplot', 'children'),
               Output('OFTprobePlot', 'children')],
              [Input('runHEAT', 'n_clicks')],
              [State('checklistPC','value'),
               State('Btrace','value'),
               State('OFtrace','value'),
               State('gyrotrace','value'),
               State('BtraceTable', 'data'),
               State('xOFtrace','value'),
               State('yOFtrace','value'),
               State('zOFtrace','value'),
               State('xGyroTrace','value'),
               State('yGyroTrace','value'),
               State('zGyroTrace','value'),
               State('timeSlider','value'),
               State('gyroT_eV_trace','value'),
               State('gyroDeg_trace','value'),
               State('N_gyroSteps_trace','value'),
               State('gyroDir_trace','value'),
               State('gyroPhase_trace','value'),
               State('accFilters','value'),
               ])
def runHEAT(n_clicks,runList,Btrace,OFtrace,gyrotrace,
            BtraceTableData,
            xOFtrace,yOFtrace,zOFtrace,
            xGyroTrace,yGyroTrace,zGyroTrace,
            t,
            gyroT_eV_trace,gyroDeg_trace,N_gyroSteps_trace,
            gyroDir_trace,gyroPhase_trace, accFilters):
    if n_clicks == 0:
        raise PreventUpdate

    #acceleration structure filters
    gui.loadAccFilters(accFilters)

    #Bfield trace
    if 'Btrace' in Btrace:
        gui.BtraceMultiple(BtraceTableData, t)

    #gyro orbit trace
    if 'gyrotrace' in gyrotrace:
        gui.gyroTrace(xGyroTrace,yGyroTrace,zGyroTrace,t,gyroPhase_trace,
                      gyroDeg_trace,N_gyroSteps_trace,gyroDir_trace,gyroT_eV_trace)

    gui.runHEAT(runList)

    if 'hfOpt' in runList:
        #load HF distribution plots on output page
        qDistFig = hfDistPlots(update=True)
    else:
        qDistFig = hfDistPlots(update=False)

    if 'T' in runList:
        gui.runOpenFOAM()
        OFminmaxFig = OFmaxTPlots(update=True)
    else:
        OFminmaxFig = OFmaxTPlots(update=False)

    if 'OFtrace' in OFtrace:
        OFTprobeFig = OFTprobePlots(update=True,x=xOFtrace,y=yOFtrace,z=zOFtrace)
        #do these again here so that the plots on output tab dont go blank
        #if the hfOpt and T boxes arent checked
        OFminmaxFig = OFmaxTPlots(update=True)
        qDistFig = hfDistPlots(update=True)
    else:
        OFTprobeFig = OFTprobePlots(update=False)

    return ([html.Label("HEAT Run Complete", style={'color':'#f5d142'})],
            qDistFig,
            OFminmaxFig,
            OFTprobeFig)



"""
==============================================================================
Tab Contents: gfile cleaner tab
"""
def buildGfileCleanerTab():
        return html.Div(
            id="gfileCleanerTab",
            children=[gfileChildren()],
            className="runContentBox",
            )

def gfileChildren():
    return html.Div(
        children = [
            html.H4("GEQDSK Tools", style={"text-align":"center", "width":"100%"}),
            #html.Div( children=MHDplot(), className="MHDplotBox2" ),
            html.H6("Loaded gFile Parameters:"),
            html.Div( children=buildGfileTable(), className="gfileTable" ),
            html.Div( children=buildGfilePlots(), className="gfileTable" ),
            html.Div( children=buildBFieldPlots(), className="gfileTable" ),
            html.H6("gFile Multipliers:"),
            gFileMultipiers(),
            html.H6("Re-define LCFS:"),
            gFileNewSep(),
            html.H6("Save New gFile:"),
            saveNewGfile(),
            html.H6("Interpolate gFile:"),
            interpolateGfile(),
            ],
        className="wideBoxDark",
        )

def buildGfilePlots():
    """
    gFile plot for Fpol, psi, etc
    """
    return html.Div(
            id="gFilePlots",
            children=[dcc.Graph(id="gFilePlot1", className=""),],
            className="gfileBox"
            )

def buildBFieldPlots():
    """
    gFile plot for Fpol, psi, etc
    """
    return html.Div(
            id="bFieldPlots",
            children=[dcc.Graph(id="bFieldPlot", className=""),],
            className="gfileBox"
            )

def buildGfileTable():
    cols = getGfileColumns()
    data = [{}]
    return dash_table.DataTable(
        id='gFileTable',
        columns = cols,
        data = data,
        style_header={'backgroundColor': 'rgb(30, 30, 30)'},
        style_cell={
            'textAlign': 'left',
            'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'white'
        },
        editable=False,
        row_deletable=False,
        )


#generic row data
def getGfileColumns():
    params = ['Parameter','Value']
    return [{'id': p, 'name': p} for p in params]

#load data
def getGfileData(t=None):
    if t is None:
        return [{}]
    idx = np.where(t==gui.MHD.timesteps)[0][0]
    g = gui.MHD.ep[idx].g
    #parameters from gfile we want to display in table
    keepValues = ('psiSep',
                  'psiAxis',
                  'shot',
                  'time',
                  'NR',
                  'NZ',
                  'RmAxis',
                  'ZmAxis',
                  'Ip')
    dict ={k: g[k] for k in keepValues}
    psiMax = g['psiRZ'].max()
    psiMin = g['psiRZ'].min()
    FpolMax = g['Fpol'].max()
    FpolMin = g['Fpol'].min()
    dict.update({'psiRZmin':psiMin})
    dict.update({'psiRZmax':psiMax})
    dict.update({'FpolMin':FpolMin})
    dict.update({'FpolMax':FpolMax})
    dict = [{'Parameter':i, 'Value':dict[i]} for i in dict]
    return dict

def gFileMultipiers():
    return html.Div(
        id="gfileMult",
        children=[
            html.Label("psiRZ Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiRZMult", className="gfileBoxInput", value="1.0"),
            html.Label("psiRZ Addition", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiRZAdd", className="gfileBoxInput", value="0.0"),
            html.Label("psiSep Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiSepMult", className="gfileBoxInput", value="1.0"),
            html.Label("psiSep Addition", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiSepAdd", className="gfileBoxInput", value="0.0"),
            html.Label("psiAxis Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiAxisMult", className="gfileBoxInput", value="1.0"),
            html.Label("psiAxis Addition", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiAxisAdd", className="gfileBoxInput", value="0.0"),
            html.Label("Fpol Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="FpolMult", className="gfileBoxInput", value="1.0"),
            html.Label("Fpol Addition", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="FpolAdd", className="gfileBoxInput", value="0.0"),
            html.Label("Bt0 Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="Bt0Mult", className="gfileBoxInput", value="1.0"),
            html.Label("Bt0 Addition", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="Bt0Add", className="gfileBoxInput", value="0.0"),
            html.Label("Ip Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="IpMult", className="gfileBoxInput", value="1.0"),
            html.Label("Ip Addition", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="IpAdd", className="gfileBoxInput", value="0.0"),
            html.Button("Apply Corrections", id="applyMult", n_clicks=0, style={'margin':'0 10px 10px 0'}),
            html.Div(id="hiddenDivMult"),
            html.Label("*Sign of Bt0 and Ip checked for helicity in traces (MAFOT)", style={'margin':'0 10px 0 10px'}),
            html.Label("**Sign of psiRZ, psiSep, psiAxis, checked for helicity in point clouds (HEAT)", style={'margin':'0 10px 0 10px'}),
        ],
        className="gfileBox",
    )

def gFileNewSep():
    return html.Div(
        id="gfileNewSep",
        children=[
            html.Label("New LCFS R Value [m]", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="newLCFSr", className="gfileBoxInput", value="NA"),
            html.Label("New LCFS Z Value [m]", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="newLCFSz", className="gfileBoxInput", value="NA"),
            html.Label("(Value of NA in Z will choose minimum psi)"),
            html.Br(),
            dcc.RadioItems(
                            id="newLCFSradio",
                            options=[
                                {'label': 'Apply to All Timesteps', 'value': 'all'},
                                {'label': 'Apply to Currently Selected Timestep', 'value': 'single'},
                                    ],
                                    value='all'
                            ),
            html.Button("Re-Define LCFS", id="newLCFSbutton", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
            html.Div(id="hiddenDivSep"),
            html.Button("Find LCFS From PFCs", id="findLCFSbutton", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
            html.Div(id="hiddenDivSep2")
        ],
        className="gfileBox",
    )

def saveNewGfile():
    return html.Div(
        id="saveGfile",
        children=[
            html.Label("New gFile Name", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="newGfileName", className="gfileBoxInput"),
            html.Button("Save New gFile", id="saveGfileButton", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
            Download(id="downloadNewGfile"),
            html.Div(id="hiddenDivSaveGfile")
        ],
        className="gfileBox",
    )

def interpolateGfile():
    params = ['filename', 'timestep[ms]']
    data = [dict([{'filename':'', 'timestep':''}][0])]
    return html.Div(
            id="interpGfile",
            children=[
                html.Label("Interpolation by Timestep", style={'margin':'0 10px 0 10px'}),
                dcc.Input(id="interpTime", className="gfileBoxInput"),
                html.Button("Interpolate this Timestep", id="interpButton", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                html.Div(id="hiddenDivInterp1"),
                html.Br(),
                dash_table.DataTable(
                    id='gFileTable2',
                    columns = ([{'id': p, 'name': p} for p in params]),
                    data = data,
                    style_header={'backgroundColor': 'rgb(30, 30, 30)'},
                    style_cell={
                        'textAlign': 'left',
                        'backgroundColor': 'rgb(50, 50, 50)',
                        'color': 'white'
                                },
                    editable=True,
                    row_deletable=False,

                ),
                html.Br(),

                html.Label("Interpolate N steps between gFiles", style={'margin':'0 10px 0 10px'}),
                html.Div(id='interpTable', className="gfileTable"), #updated from MHD button callback
                dcc.Input(id="interpN", className="gfileBoxInput", placeholder='Enter N steps'),
                html.Button("Interpolate these Timesteps", id="interpButton2", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                Download(id="download1InterpGfile"),
                Download(id="downloadInterpGfiles"),
                html.Div(id="hiddenDivInterp2")
            ],
            className="gfileBox",
            )


@app.callback([Output('hiddenDivInterp1', 'children'),
               Output('download1InterpGfile', 'data')],
              [Input('interpButton', 'n_clicks')],
              [State('interpTime','value'),
              ])
def interpolate(n_clicks, t):
    """
    interpolate gfile at user defined timestep
    """
    if n_clicks < 1:
        raise PreventUpdate
    name = gui.interpolateGfile(t)
    return [html.Label("Gfile Interpolated", style={'color':'#f5d142'}),
            send_file(name)]


@app.callback([Output('hiddenDivInterp2', 'children'),
               Output('downloadInterpGfiles', 'data')],
              [Input('interpButton2', 'n_clicks')],
              [State('interpN','value'),
               State('gFileTable2','data'),
              ])
def interpolateNsteps(n_clicks, N, data):
    """
    interpolate gfile at user defined steps between two gfiles
    """
    if n_clicks < 1:
        raise PreventUpdate
    #load interpolation table data
    df = pd.DataFrame(data)
    df = df.sort_values('timestep[ms]')
    print(df)
    print(df['timestep[ms]'].values)
    #interpolate N steps between each point
    gui.interpolateNsteps(df['filename'].values, pd.to_numeric(df['timestep[ms]']).values,int(N))
    zipFile = gui.tmpDir + 'InterpolatedGfiles.zip'
    return [html.Label("gFiles Interpolated", style={'color':'#f5d142'}),
            send_file(zipFile)]





@app.callback(Output('hiddenDivMult', 'children'),
              [Input('applyMult', 'n_clicks')],
              [State('psiRZMult','value'),
               State('psiSepMult','value'),
               State('psiAxisMult','value'),
               State('FpolMult','value'),
               State('psiRZAdd','value'),
               State('psiSepAdd','value'),
               State('psiAxisAdd','value'),
               State('FpolAdd','value'),
               State('Bt0Mult','value'),
               State('Bt0Add','value'),
               State('IpMult','value'),
               State('IpAdd','value'),
               State('timeSlider', 'value')])
def applyMult(n_clicks, psiRZMult, psiSepMult, psiAxisMult, FpolMult,
              psiRZAdd,psiSepAdd,psiAxisAdd,FpolAdd,
              Bt0Mult,Bt0Add,IpMult,IpAdd,t):
    """
    apply multiplier to psiRZ, psiSep, psiAxis, Fpol for currently
    selected equilibrium timestep
    """
    if n_clicks < 1:
        raise PreventUpdate
    #parse user formulas and convert then to number via python compiler
    pi = np.pi

    #THIS METHOD NOT SAFE BECAUSE ANY CODE CAN BE INSERTED
#    psiRZMult = eval(parser.expr(psiRZMult).compile())
#    psiSepMult = eval(parser.expr(psiSepMult).compile())
#    psiAxisMult = eval(parser.expr(psiAxisMult).compile())
#    FpolMult = eval(parser.expr(FpolMult).compile())
#    Bt0Mult = eval(parser.expr(Bt0Mult).compile())
#    IpMult = eval(parser.expr(IpMult).compile())

#    psiRZAdd = eval(parser.expr(psiRZAdd).compile())
#    psiSepAdd = eval(parser.expr(psiSepAdd).compile())
#    psiAxisAdd = eval(parser.expr(psiAxisAdd).compile())
#    FpolAdd = eval(parser.expr(FpolAdd).compile())
#    Bt0Add = eval(parser.expr(Bt0Add).compile())
#    IpAdd = eval(parser.expr(IpAdd).compile())


    try:
        psiRZMult = float(psiRZMult)
        psiSepMult = float(psiSepMult)
        psiAxisMult = float(psiAxisMult)
        FpolMult = float(FpolMult)
        Bt0Mult = float(Bt0Mult)
        IpMult = float(IpMult)

        psiRZAdd = float(psiRZAdd)
        psiSepAdd = float(psiSepAdd)
        psiAxisAdd = float(psiAxisAdd)
        FpolAdd = float(FpolAdd)
        Bt0Add = float(Bt0Add)
        IpAdd = float(IpAdd)
        returnDiv = [html.Label("Corrections Applied", style={'color':'#f5d142'})]

    except:
        print("Could not evaluate multipliers or adders")
        print("Please ensure values are of type float")
        print("and retry")
        returnDiv = [html.Label("Error with inputs", style={'color':'#f70c0c'})]


    gui.gfileClean(psiRZMult,psiSepMult,psiAxisMult,FpolMult,
                   psiRZAdd,psiSepAdd,psiAxisAdd,FpolAdd,
                   Bt0Mult,Bt0Add,IpMult,IpAdd,t)
    return returnDiv

@app.callback(Output('hiddenDivSep', 'children'),
              [Input('newLCFSbutton', 'n_clicks')],
              [State('newLCFSr','value'),
               State('newLCFSz','value'),
               State('timeSlider', 'value'),
               State('newLCFSradio', 'value')])
def newLCFSbutton(n_clicks, rNew, zNew, t, radio):
    if n_clicks < 1:
        raise PreventUpdate
    if radio=='all':
        gui.newLCFSallTimesteps(rNew, zNew)
    else:
        gui.newLCFS(t, rNew, zNew)
    return [html.Label("Redefined LCFS", style={'color':'#f5d142'})]

@app.callback(Output('hiddenDivSep2', 'children'),
              [Input('findLCFSbutton', 'n_clicks')],
              [State('timeSlider', 'value'),
               State('loadPFC', 'n_clicks'),
               State('newLCFSr','value')])
def findLCFSbutton(n_clicks, t, PFC_n_clicks, r):
    if n_clicks < 1:
        raise PreventUpdate
    if PFC_n_clicks < 1:
        print("You must load a PFC before running this function")
        raise PreventUpdate
    gui.findPsiSepfromPFCs(t, r)
#    gui.findPsiSepfromEQ(t)
    return [html.Label("Iteratively Found LCFS", style={'color':'#f5d142'})]



@app.callback([Output('hiddenDivSaveGfile', 'children'),
                Output('downloadNewGfile', 'data')],
              [Input('saveGfileButton', 'n_clicks')],
              [State('newGfileName','value'),
               State('timeSlider', 'value'),
               State('shot', 'value')])
def saveG(n_clicks, filename, t, shot):
    if n_clicks < 1:
        raise PreventUpdate
    gui.writeGfile(filename, shot, t)
    return [html.Label("Saved gFile", style={'color':'#f5d142'}),
            send_file(gui.tmpDir + filename)]

"""
==============================================================================
Tab Contents: Output Tab
"""
def buildOutputTab():
        return html.Div(
            id="outputTab",
            children=[outputChildren()],
            )

def outputChildren():
    return html.Div(
        children = [
            html.H4("HEAT outputs", style={"text-align":"center", "width":"100%"}),
            html.Button("Download HEAT Results", id="downloadResults", n_clicks=0, style={'margin':'0 10px 10px 0', "width":"95%" }),
            html.Br(),
            html.Div(id="hiddenDivDownloadResults", style={"width":"100%"}),
            Download(id="downloadResultsDir"),
            html.H6("Input Parameters:"),
            html.Div( children=buildInputsTable(), className="gfileTable" ),
            #qDiv plot
            html.Div(
                children=[
                    html.H6("HEAT qDiv Distributions (run HEAT for plot):"),
                    html.Div(children=hfDistPlots(update=False), id="qDivDist"),
                    ],
                className="wideBoxDarkColumn",
            ),
            #Temp(t) plot for max(T) in OF run
            html.Div(
                children=[
                    html.H6("openFOAM max(T) Evolutions (run openFOAM for plot):"),
                    html.Div(children=OFmaxTPlots(update=False), id="OFmaxTplot"),
                    ],
                className="wideBoxDarkColumn",
            ),
            #Temp(t) for probes from GUI
            html.Div(
                children=[
                    html.H6("openFOAM Tprobe Evolution (run openFOAM Tprobe for plot):"),
                    html.Div(children=OFTprobePlots(update=False), id="OFTprobePlot"),
                    ],
                className="wideBoxDarkColumn",
            ),
            #gyroPhase plot
            html.Div(
                children=[
                    html.H6("HEAT gyroPhase angles (load gyro settings for plot):"),
                    html.Div(children=gyroPhasePlots(update=False), id="gyroPhasePlot"),
                    ],
                className="wideBoxDarkColumn",
            ),
            #vPhase plot
            html.Div(
                children=[
                    html.H6("HEAT vPhase angles (run HEAT for plot):"),
                    html.Div(children=vPhasePlots(update=False), id="vPhasePlot"),
                    ],
                className="wideBoxDarkColumn",
            ),
            #vSlice plot
            html.Div(
                children=[
                    html.H6("HEAT vSlices (run HEAT for plot):"),
                    html.Div(children=vSlicePlots(update=False), id="vSlicePlot"),
                    ],
                className="wideBoxDarkColumn",
            ),
            #cdfSlice plot
            html.Div(
                children=[
                    html.H6("HEAT PDF vs. CDF Slices (run HEAT for plot):"),
                    html.Div(children=cdfSlicePlots(update=False), id="cdfSlicePlot"),
                    ],
                className="wideBoxDarkColumn",
            ),
            ],
        className="wideBoxDark",
        )

@app.callback([Output('hiddenDivDownloadResults', 'children'),
               Output('downloadResultsDir', 'data')],
              [Input('downloadResults', 'n_clicks')],
              [State('MachFlag','value'),
               State('shot', 'value')])
def saveResults(n_clicks, MachFlag, shot):
    if n_clicks is None:
        raise PreventUpdate
    if MachFlag is None:
        raise PreventUpdate
    #create zip archive of data folder for this shot #
    file = gui.tmpDir + 'HEATresults'
    print("Creating HEAT results zip.  This may take a while for large OF runs")
    log.info("Creating HEAT results zip.  This may take a while for large OF runs")
    shutil.make_archive(file, 'zip', gui.MHD.shotPath)
    print("Zipped results")
    log.info("Zipped results")
    return [html.Label("Saved HEAT output", style={'color':'#f5d142'}),
            send_file(file+'.zip')]



def buildInputsTable(data=None):
    cols = getOutputColumns()
    if data==None:
        data = [{}]
    return dash_table.DataTable(
        id='inputsTable',
        columns = cols,
        data = data,
        style_header={'backgroundColor': 'rgb(30, 30, 30)'},
        style_cell={
            'textAlign': 'left',
            'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'white'
        },
        editable=False,
        row_deletable=False,
        export_format="csv",
        )


#generic row data (table data is created in loadHF button connect)
def getOutputColumns():
    params = ['Parameter','Value']
    return [{'id': p, 'name': p} for p in params]

@app.callback([Output('inputsTable', 'data')],
              [Input('MHDDataStorage', 'data'),
               Input('CADDataStorage', 'data'),
               Input('HFDataStorage', 'data'),
               Input('GYRODataStorage', 'data'),
               Input('OFDataStorage', 'data')])
def inputsTable(MHDdata,CADdata,HFdata,GYROdata,OFdata):
    newData = {}
    if MHDdata != None:
        for key,val in MHDdata.items():
            newData[key] = val
    if CADdata != None:
        for key,val in CADdata.items():
            newData[key] = val
    if HFdata != None:
        for key,val in HFdata.items():
            newData[key] = val
    if GYROdata != None:
        for key,val in GYROdata.items():
            newData[key] = val
    if OFdata != None:
        for key,val in OFdata.items():
            newData[key] = val
    newData = [{'Parameter': p, 'Value': newData[p]} for p in newData]
    return [newData]




def hfDistPlots(update=False):
    """
    div for heat flux distribution plots

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.getHFdistPlots()

        return html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        return html.Div(
            children=[
                html.Label("Run HEAT to get qDiv plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )


def OFmaxTPlots(update=False):
    """
    div for maximum openFOAM temperature plots

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.getOFMinMaxPlots()

        return html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        return html.Div(
            children=[
                html.Label("Run OpenFOAM to get max(T) plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )


def OFTprobePlots(update=False,x=None,y=None,z=None):
    """
    div for openFOAM probe plots

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    x,y,z are coordinate locations [mm] of T probe

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.TprobeOF(float(x),float(y),float(z))

        return html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        return html.Div(
            children=[
                html.Label("Run OpenFOAM Temp Probe to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )


def gyroPhasePlots(update=False):
    """
    div for gyro orbit phase angle plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.gyroPhasePlot()

        return html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        return html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )

def vPhasePlots(update=False):
    """
    div for gyro orbit velocity phase angle plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.vPhasePlot()

        return html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        return html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )

def vSlicePlots(update=False):
    """
    div for gyro orbit velocity slice / distribution plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.vSlicePlot()

        return html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        return html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )

def cdfSlicePlots(update=False):
    """
    div for gyro orbit CDF slice / distribution plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.cdfSlicePlot()

        return html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        return html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )




"""
==============================================================================
Tab Contents: logfile tab
"""
def buildLogTab():
        return html.Div(
            id="logTab",
            children=[
                        html.H4("HEAT Log File Updated Every 5 Seconds"),
                        dcc.Textarea(id="logData", value='', className="logBox",
                            readOnly=True),
                    ],

            className="wideBoxDark"
            )

#callback for updating log tab every 5 seconds
@app.callback([Output('logData', 'value'),
               #Output('javascriptLog', 'run')
               ],
              [Input('intervalLog', 'n_intervals')])
def updateLogFile(n):
    if n < 1:
        raise PreventUpdate
    with open(logFile, "r") as f:
        content = f.read()
    #for visdcc to scroll to bottom every 5 sec
    logCMD = '''
             var textarea = document.getElementById('logData');
             textarea.scrollTop = textarea.scrollHeight;

             '''
    return [content]#, logCMD


"""
==============================================================================
Graphs, plots, etc.
"""
def build_graphs():
    return html.Div(
        id="graph-container",
        children=[
            MHDplot(),
        ],
    )

def MHDplot():

    return html.Div(
        className="MHDplotBox",
        children=[
            dcc.Graph(id="2DEQ", className="EQplot"),
            html.Br(),
            dcc.Slider(
                id='timeSlider',
                min=0,
                max=100,
                step=None,
                value=10,
                marks={50: 'Load MHD to see timesteps'},
            ),
        ],
    )

#this callback is fired any time user moves time slider under EQ plot,
#or when user loads mhd (which changes time slider value)
@app.callback([Output('2DEQ', 'figure'),
               Output('gFileTable', 'data'),
               Output('gFilePlot1', 'figure'),
               Output('bFieldPlot', 'figure')],
              [Input('timeSlider', 'value'),
               Input('hiddenDivMult', 'children'),
               Input('hiddenDivSep', 'children'),])
def slideEQplot(value, dummy1, tom):
    try: MachFlag = gui.MachFlag
    except:
        print("You didn't select a machine")
        raise PreventUpdate
    try: ts = gui.MHD.timesteps
    except:
        print('Please load MHD')
        raise PreventUpdate
    idx = np.where(value==gui.MHD.timesteps)[0][0]
    ep = gui.MHD.ep[idx]
    shot = gui.MHD.shot
    t = gui.MHD.timesteps[idx]
    #Update Equilibrium plot
    #this import needs to be here (not at top of file) to prevent FreeCAD qt5
    #shared object conflict
    import GUIscripts.plotly2DEQ as plotly2DEQ
    plot = plotly2DEQ.makePlotlyEQDiv(shot, t, MachFlag, ep)
    data = getGfileData(t)

    #update plot on gfile cleaner tab with Fpol, psi. etc
    plot2 = plotly2DEQ.makePlotlyGfilePlot(ep)

    #update plot on gfile cleaner tab with Bp, Bt
    plot3= plotly2DEQ.makePlotlyBpBt(ep, MachFlag)

    return plot, data, plot2, plot3








"""
==============================================================================
Main Layout
"""
def build_simulator():
    return html.Div(
        id="simulator-container",
        children=[
            build_tabs(),
            build_graphs(),
        ],
    )

logJScmd = ""
app.layout = html.Div(
        id="big-app-container",
        children=[
            #store the session data until browser tab is closed
            dcc.Store(id='session', storage_type='memory'),
            dcc.Store(id='userInputFileData', storage_type='memory'),
            dcc.Store(id='PFCdataStorage', storage_type='memory'),
            dcc.Store(id='gFileListStorage', storage_type='memory'),
            dcc.Store(id='BtraceDataStorage', storage_type='memory'),
            dcc.Store(id='MHDDataStorage', storage_type='memory'),
            dcc.Store(id='CADDataStorage', storage_type='memory'),
            dcc.Store(id='HFDataStorage', storage_type='memory'),
            dcc.Store(id='GYRODataStorage', storage_type='memory'),
            dcc.Store(id='OFDataStorage', storage_type='memory'),

            #interval for updating logFile tab every 5 seconds
            dcc.Interval(
                id='intervalLog',
                interval=5*1000, # in milliseconds
                n_intervals=0
                ),
            #visdcc to run js to scroll log box to bottom
            #visdcc.Run_js(id = 'javascriptLog', run = ""),
            build_banner(),
            build_simulator(),
            html.Div(id="hiddenDiv", style={'display': 'none'}),
        ],
)

app.title = 'HEAT'



"""
==============================================================================
Session storage callbacks and functions
"""
# Load Default callback and input file drag and drop
@app.callback([Output('shot', 'value'),
               Output('tmin', 'value'),
               Output('tmax', 'value'),
               Output('nTrace', 'value'),
               Output('gridRes', 'value'),
               Output('lqEich', 'value'),
               Output('S', 'value'),
               Output('lqCN', 'value'),
               Output('lqCF', 'value'),
               Output('lqPN', 'value'),
               Output('lqPF', 'value'),
               Output('limlqCN', 'value'),
               Output('limlqCF', 'value'),
               Output('limfracCN', 'value'),
               Output('limfracCF', 'value'),
               Output('fracCN', 'value'),
               Output('fracCF', 'value'),
               Output('fracPN', 'value'),
               Output('fracPF', 'value'),
               Output('Pinj', 'value'),
               Output('coreRadFrac', 'value'),
               Output('fracUI', 'value'),
               Output('fracUO', 'value'),
               Output('fracLI', 'value'),
               Output('fracLO', 'value'),
               Output('qBG', 'value'),
               Output('fG', 'value'),
               Output('qFilePath', 'value'),
               Output('qFileTag', 'value'),
               Output('OFstartTime', 'value'),
               Output('OFstopTime', 'value'),
               Output('OFminMeshLev', 'value'),
               Output('OFmaxMeshLev', 'value'),
               Output('OFSTLscale', 'value'),
               Output('OFdeltaT', 'value'),
               Output('OFwriteDeltaT', 'value'),
               Output('N_gyroSteps','value'),
               Output('gyroDeg','value'),
               Output('gyroT_eV','value'),
               Output('N_vSlice','value'),
               Output('N_vPhase','value'),
               Output('N_gyroPhase','value'),
               Output('ionMassAMU','value'),
               #Output('vMode','value'),
               Output('ionFrac','value'),
               Output('session', 'data'),
#               Output('hiddenDivMachFlag', 'children')
               ],
               [Input('loadDefaults', 'n_clicks'),
                Input('userInputFileData', 'modified_timestamp')],
               [State('session', 'modified_timestamp'),
                State('MachFlag', 'value'),
                State('session', 'data'),
                State('userInputFileData', 'data')])
def session_data(n_clicks, inputTs, ts, MachFlag, data, inputFileData):
    #default case
    if ts is None or MachFlag not in machineList:
        print('Initializing Data Store')
        if MachFlag not in machineList:
            print("Machine not in machine list")
            log.info("Machine not in machine list")
        data = gui.getDefaultDict()
        data.update({'default_n_clicks':n_clicks})

    #Let the user know if this worked or if we still need a MachFlag
    if MachFlag not in machineList:
        outputDiv = html.Label("Select Machine First", style={'color':'#f51b60'})
    else:
        outputDiv = html.Label("Loaded Input File", style={'color':'#f5d142'})

    #load defaults
    if n_clicks > 0 and n_clicks>data.get('default_n_clicks') and MachFlag is not None:
        print("Loading Default Input File")
        log.info("Loading Default Input File")
        data = gui.loadDefaults()
        data['default_n_clicks'] = n_clicks
    elif inputTs == None or inputTs == -1:
        pass
    #use data we saved into storage object that we got from user input file
    else:
        print(inputTs)
        print("Loading User Input File")
        log.info("Loading User Input File")
        inputFileData['default_n_clicks'] = n_clicks
        data = inputFileData

    print("Updating Session Store")
    return [data.get('shot', ''),
            data.get('tmin', ''),
            data.get('tmax', ''),
            data.get('nTrace', ''),
            data.get('gridRes', ''),
            data.get('lqEich', ''),
            data.get('S', ''),
            data.get('lqCN', ''),
            data.get('lqCF', ''),
            data.get('lqPN', ''),
            data.get('lqPF', ''),
            data.get('limlqCN', ''),
            data.get('limlqCF', ''),
            data.get('limfracCN', ''),
            data.get('limfracCF', ''),
            data.get('fracCN', ''),
            data.get('fracCF', ''),
            data.get('fracPN', ''),
            data.get('fracPF', ''),
            data.get('Pinj', ''),
            data.get('coreRadFrac', ''),
            data.get('fracUI', ''),
            data.get('fracUO', ''),
            data.get('fracLI', ''),
            data.get('fracLO', ''),
            data.get('qBG', ''),
            data.get('fG', ''),
            data.get('qFilePath', ''),
            data.get('qFileTag', ''),
            data.get('OFtMin', ''),
            data.get('OFtMax', ''),
            data.get('meshMinLev', ''),
            data.get('meshMaxLev', ''),
            data.get('STLscale', ''),
            data.get('deltaT', ''),
            data.get('writeDeltaT', ''),
            data.get('N_gyroSteps',''),
            data.get('gyroDeg',''),
            data.get('gyroT_eV',''),
            data.get('N_vSlice',''),
            data.get('N_vPhase',''),
            data.get('N_gyroPhase',''),
            data.get('ionMassAMU',''),
            #data.get('vMode',''),
            data.get('ionFrac',''),
            data,
#            outputDiv
            ]



if __name__ == '__main__':
# dashGUI.py should be called with no arguments to run on 127.0.0.1:8050 (default)
# to run on address:port, use command line:
# dashGUI.py <address> <port>
    print("Running as MAIN")

    try:
        ipaddress.ip_address(sys.argv[1]) #validate it is an ip address
        address = sys.argv[1]
    except:
        address = '127.0.0.1' #default
    try:
        port = int(sys.argv[2])
    except:
        port = 8050 # default


    app.run_server(
                    debug=True,
                    dev_tools_ui=True,
                    port=port,
                    host=address
                   )

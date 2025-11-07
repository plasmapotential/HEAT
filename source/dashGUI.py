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
rootDir
"""
import os
import sys
import shutil
import base64
import io
import json
import numpy as np
import pandas as pd
import time
import copy
from dash import Dash, dcc, html, dash_table, ctx
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template, template_from_url
from dash_bootstrap_templates import ThemeSwitchAIO, ThemeChangerAIO
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import visdcc
from flask import Flask, send_from_directory
import plotly.io as pio
import EFIT.equilParams_class as EP
import toolsClass
tools = toolsClass.tools()
import ipaddress
import logging
log = logging.getLogger(__name__)

#get relevant environment variables
logFile = os.environ["logFile"]
rootDir = os.environ["rootDir"]
dataPath = os.environ["dataPath"]
OFbashrc = os.environ["OFbashrc"]
FreeCADPath = os.environ["FreeCADPath"]

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

#app = dash.Dash(server=server, meta_tags=[{"name": "viewport", "content": "width=device-width"}],
#                prevent_initial_callbacks=False)

app = Dash(server=server, meta_tags=[{"name": "viewport", "content": "width=device-width"}],
                prevent_initial_callbacks=False, update_title='Processing...')


#load the HEAT flame icon into browser tab
iconFile = "flame.ico" #located in assets folder
if os.path.exists(rootDir + "/assets/flame.ico"):
    app._favicon = ("flame.ico")
else:
    print("Cannot find icon file:")
    print(rootDir + "/assets/flame.ico\n")


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
                className="rowBox",
                style={"margin": "10px 10px 10px 10px"},
                children=[
                    html.Div(style={'width':'10%'}),
                    html.H4("Heat flux Engineering Analysis Toolkit (HEAT)", style={'width':'80%'}),
                    html.Div(themeToggle(), style={'width':'10%', 'justify-content':'flex-end'})
                ],
            ),
        ],
    )

def build_tabs():
    """
    returns tab bar
    """
    return html.Div(
        id="tabs",
        className="tabWindow",
        children=[
            dbc.Tabs(
                id="app-tabs",
                children=[
                    dbc.Tab(
                        id="input-tab",
                        label="Inputs",
                        children=[
                            html.Div(
                                [build_accordian()],
                                className="tabContent"

                                )
                                ]
                            ),
                    dbc.Tab(
                        id="run-tab",
                        label="Run HEAT",
                        children=[
                            html.Div(
                                [buildRunTab()],
                                className="tabContent",
                                )
                                ]
                            ),
                    dbc.Tab(
                        id="gfileCleaner-tab",
                        label="GEQDSK Tools",
                        children=[
                            html.Div(
                                [buildGfileCleanerTab()],
                                className="tabContent",
                                )
                                ]
                                ),

                    dbc.Tab(
                        id="output-tab",
                        label="Output",
                        children=[
                            html.Div(
                                [buildOutputTab()],
                                className="tabContent",
                                )
                            ],
                        ),

                ],
            )
        ],
    )



def build_accordian():
    """
    returns accordian
    """
    return dbc.Accordion(
        id="accordian",
        className="accordion",
        start_collapsed=True,
        flush=True,
        children=[
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody([buildButtonRibbon()])
                    )
                ],
                title="Machine Selection and Input Files",

            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody([buildMHDbox()])
                    )
                ],
                title="MHD Settings",
            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody(buildCADbox())
                    )
                ],
                title="CAD Settings"
            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody(buildPFCbox(),)
                    )
                ],
                title="PFC Settings"
            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody(buildHFbox())
                    )
                ],
                title="Optical HF Settings"
            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody(buildGYRObox(),)
                    )
                ],
                title="Gyro HF Settings"
            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody(buildRADbox(),)
                    )
                ],
                title="Photon HF Settings"
            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody(buildOFbox(),),
                    )
                ],
                title="OpenFOAM Settings"
            ),
            dbc.AccordionItem(
                [
                dbc.Card(
                    dbc.CardBody([buildDefaultPaths()]),
                    #color='secondary',
                    )
                ],
                title="Default Paths (read only)",
            ),
            ],
        )


"""
==============================================================================
Tab Contents: inputs tab
"""
def themeToggle():
    """
    returns a div with theme toggling
    """

    themes_list = [
        dbc.themes.COSMO,
        dbc.themes.SLATE
        ]

    themeToggle = html.Div(
        [
            #for toggling between two themes
            #ThemeChangerAIO(aio_id="theme", themes=themes_list, )
            #for switching between multiple themes
            dbc.Row(ThemeChangerAIO(aio_id="theme",
                                    radio_props={"value":dbc.themes.SLATE},
                                    offcanvas_props={'placement':'end'}
                                    ),
                    ),
            html.Div(id="hiddenDivTheme", style={'display':'None'}),
        ]
        )

    return themeToggle


@app.callback(Output("themeStorage", 'data'),
              Input(ThemeChangerAIO.ids.radio("theme"), "value"))
def switchStyleSheet(theme):
    """
    callback to switch figure theme
    """
    data = {"theme":theme}
    return data




def buildButtonRibbon():
    """
    returns machine selector and load default buttons
    """
    return html.Div(
        id="buttonRib",
        children=[
            #html.H5("Machine Selection and Input Files"),
            html.Div(
                children=[
                    buildMachineSelector(),
                    buildInputButtons(),
                    ],
            )
        ],
    )

def buildDefaultPaths():
    """
    contains text boxes for HEAT relevent paths
    FreeCAD is location of freecad installation
    dataDir is location where HEAT output will be saved

    rootDir is location of HEAT source code and is not included in GUI
    because it would be impossible to run GUI if this was not already set.
    rootDir (and some other defaults) are hardcoded at the top of this file

    if className is Hidden, then the input exists in html but is hidden from user
    """
    return html.Div(
        #id="defaultPaths",
        className="column1",
        children=[
            dbc.Label("ParaVIEW Path"),
            dbc.Label("FreeCAD Path"),
            dbc.Input(id="FreeCADPath", value=FreeCADPath, readonly=True),
            dbc.Label("Data Directory"),
            dbc.Input(id="dataPath", value=dataPath, readonly=True),
            dbc.Label("OpenFOAM bashrc file"),
            dbc.Input(id="OFbashrc", value=OFbashrc, readonly=True),
        ],
    )

def buildMachineSelector():
    """
    returns machine selector dropdown
    """
    return html.Div(
            className="column1",
            children=[
                dbc.Label(id="machLabel", children="Select a Tokamak"),
                dcc.Dropdown(
                    id='MachFlag',
                    className="SelectorBoxInput",
                    options=[
                        {'label': 'SPARC', 'value': 'sparc'},
                        {'label': 'ARC', 'value': 'arc'},
                        {'label': 'CMOD', 'value': 'cmod'},
                        {'label': 'NSTXU', 'value': 'nstx'},
                        {'label': 'DIIID', 'value': 'd3d'},
                        {'label': 'ST40', 'value': 'st40'},
                        {'label': 'STEP', 'value': 'step'},
                        {'label': 'WEST', 'value': 'west'},
                        {'label': 'KSTAR', 'value': 'kstar'},
                        {'label': 'AUG', 'value': 'aug'},
                        {'label': 'TCV', 'value': 'tcv'},
                        {'label': 'OTHER', 'value': 'other'}
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
            className="column2",
            children=[
                html.Div(id="hiddenDivInput"),
                dbc.Button("Load Defaults (optional)", id="loadDefaults", n_clicks=0, color="primary"),
                html.Div(id="hiddenDivLoadDefaults"),
                html.Br(),
                dbc.Button("Save Settings\nInto Input File", id="saveInputs",
                            n_clicks=0, color="secondary"),
                dcc.Download(id="downloadInputs"),
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
        return [html.Label(machFlagChosen, className="text-warning")]
    gui.machineSelect(MachFlag, machineList)
    machFlagChosen = "Selected "+MachFlag
    return [html.Label(machFlagChosen, className="text-success")]


@app.callback([Output('hiddenDivInput', 'children'),
               Output('userInputFileData', 'data'),],
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
        outputDiv = dbc.Label("Loaded Input File", className="text-success")
        newFile = gui.tmpDir + file[0]
        decoded = base64.b64decode(contents[0].split(',')[1])
        #Save user loaded file into tmp directory
        with open(newFile, "w") as f:
            f.write(decoded.decode('utf-8'))
        data = gui.loadInputs(inFile=newFile)

    return [outputDiv, data]



@app.callback([Output('hiddenDivSaveInput','children'),
               Output('downloadInputs', 'data')],
              [Input('saveInputs','n_clicks')],
              [State('shot', 'value'),
               State('traceLength', 'value'),
               State('dpinit', 'value'),
               State('psiMult1', 'value'),
               State('BtMult1', 'value'),
               State('IpMult1', 'value'),
               State('gridRes', 'value'),
               State('hfMode', 'value'),
               State('eichlqCNmode', 'value'),
               State('eichSMode', 'value'),
               State('multiExplqCNmode', 'value'),
               State('multiExplqCFmode', 'value'),
               State('multiExplqPNmode', 'value'),
               State('multiExplqPFmode', 'value'),
               State('limiterlqCNmode', 'value'),
               State('limiterlqCFmode', 'value'),
               State('lqEich', 'value'),
               State('S', 'value'),
               State('lqCN', 'value'),
               State('lqCF', 'value'),
               State('lqPN', 'value'),
               State('lqPF', 'value'),
               State('lqTopHat', 'value'),
               State('fracCN', 'value'),
               State('fracCF', 'value'),
               State('fracPN', 'value'),
               State('fracPF', 'value'),
               State('limlqCN', 'value'),
               State('limlqCF', 'value'),
               State('limfracCN', 'value'),
               State('limfracCF', 'value'),
               State('P', 'value'),
               State('radFrac', 'value'),
               State('fracUI', 'value'),
               State('fracUO', 'value'),
               State('fracLI', 'value'),
               State('fracLO', 'value'),
               State('qBG', 'value'),
               State('fG', 'value'),
               State('qFilePath', 'value'),
               State('qFileTag', 'value'),
               State('OFstartTime', 'value'),
               State('OFstopTime', 'value'),
               State('OFminMeshLev', 'value'),
               State('OFmaxMeshLev', 'value'),
               State('OFSTLscale', 'value'),
               State('OFdeltaT', 'value'),
               State('OFwriteDeltaT', 'value'),
               State('dataPath', 'value'),
               State('N_gyroSteps','value'),
               State('gyroTraceLength','value'),
               State('gyroT_eV','value'),
               State('N_vSlice','value'),
               State('N_vPhase','value'),
               State('N_gyroPhase','value'),
               State('ionMassAMU','value'),
               State('vMode','value'),
               State('ionFrac','value'),
               State('phiMin', 'value'),
               State('phiMax','value'),
               State('Ntor','value'),
               State('Nref','value'),
               ]
               )
def saveGUIinputs(  n_clicks,
                    shot,
                    traceLength,
                    dpinit,
                    psiMult,
                    BtMult,
                    IpMult,
                    gridRes,
                    hfMode,
                    eichlqCNmode,
                    eichSMode,
                    multiExplqCNmode,
                    multiExplqCFmode,
                    multiExplqPNmode,
                    multiExplqPFmode,
                    limiterlqCNmode,
                    limiterlqCFmode,
                    lqEich,
                    S,
                    lqCN,
                    lqCF,
                    lqPN,
                    lqPF,
                    lqTopHat,
                    fracCN,
                    fracCF,
                    fracPN,
                    fracPF,
                    limlqCN,
                    limlqCF,
                    limFracCN,
                    limFracCF,
                    P,
                    radFrac,
                    fracUI,
                    fracUO,
                    fracLI,
                    fracLO,
                    qBG,
                    fG,
                    qFilePath,
                    qFileTag,
                    OFstartTime,
                    OFstopTime,
                    OFminMeshLev,
                    OFmaxMeshLev,
                    OFSTLscale,
                    OFdeltaT,
                    OFwriteDeltaT,
                    dataLoc,
                    N_gyroSteps,
                    gyroTraceLength,
                    gyroT_eV,
                    N_vSlice,
                    N_vPhase,
                    N_gyroPhase,
                    ionMassAMU,
                    vMode,
                    ionFrac,
                    phiMin,
                    phiMax,
                    Ntor,
                    Nref,
                ):
    """
    Saves GUI text boxes into an input file in the HEAT format
    first file is saved to the gui.tmpDir, then to client machine
    """
    if n_clicks < 1:
        raise PreventUpdate


    #get an empty dictionary with all input file parameters as None
    data = gui.getDefaultDict()
    #cad variables
    data['gridRes'] = gridRes
    #mhd variables
    data['shot'] = shot
    data['traceLength'] = traceLength
    data['dpinit'] = dpinit
    data['dataPath'] = dataLoc
    data['psiMult'] = psiMult
    data['BtMult'] = BtMult
    data['IpMult'] = IpMult

    #hf variables
    data['hfMode'] = hfMode
    if hfMode == 'multiExp':
        data['lqCN'] = lqCN
        data['lqCF'] = lqCF
        data['lqPN'] = lqPN
        data['lqPF'] = lqPF
        data['lqCNmode'] = multiExplqCNmode
        data['lqCFmode'] = multiExplqCFmode
        data['lqPNmode'] = multiExplqPNmode
        data['lqPFmode'] = multiExplqPFmode
        data['fracCN'] = fracCN
        data['fracCF'] = fracCF
        data['fracPN'] = fracPN
        data['fracPF'] = fracPF

    elif hfMode == 'limiter':
        data['lqCN'] = limlqCN
        data['lqCF'] = limlqCF
        data['lqCNmode'] = limiterlqCNmode
        data['lqCFmode'] = limiterlqCFmode
        data['fracCN'] = limFracCN
        data['fracCF'] = limFracCF

    elif hfMode == 'qFile':
        data['qFilePath'] = qFilePath
        data['qFileTag'] = qFileTag

    elif hfMode == 'eich': #gaussian spreading
        data['lqCN'] = lqEich
        data['lqCNmode'] = eichlqCNmode
        data['S'] = S
        data['SMode'] = eichSMode

    elif hfMode == 'tophat':
        data['lqCN'] = lqTopHat
        data['lqCNmode'] = 'user'
    
    # elif hfMode == 'rzqprofile':
    #     data['rzqFile'] = rzqFile

    data['fracUI'] = fracUI
    data['fracUO'] = fracUO
    data['fracLI'] = fracLI
    data['fracLO'] = fracLO
    data['P'] = P
    data['radFrac'] = radFrac
    data['qBG'] = qBG
    data['fG'] = fG
    #openfoam variables
    data['OFtMin'] = OFstartTime
    data['OFtMax'] = OFstopTime
    data['meshMinLevel'] = OFminMeshLev
    data['meshMaxLevel'] = OFmaxMeshLev
    data['STLscale'] = OFSTLscale
    data['deltaT'] = OFdeltaT
    data['writeDeltaT'] = OFwriteDeltaT
    #gyro variables
    data['N_gyroSteps'] = N_gyroSteps
    data['gyroTraceLength'] = gyroTraceLength
    data['gyroT_eV'] = gyroT_eV
    data['N_vSlice'] = N_vSlice
    data['N_vPhase'] = N_vPhase
    data['N_gyroPhase'] = N_gyroPhase
    data['ionMassAMU'] = ionMassAMU
    data['vMode'] = vMode
    data['ionFrac'] = ionFrac
    # data['radFile'] = radFile
    data['phiMin'] = phiMin
    data['phiMax'] = phiMax
    data['Ntor'] = Ntor
    data['Nref'] = Nref

    tools.saveInputFile(data, gui.tmpDir, gui.rootDir, gui.dataPath)

    outputDiv = html.Label("Saved File", className="text-success")
    return [outputDiv, dcc.send_file(gui.tmpDir + "HEATinput.csv")]




#==========MHD==========
def buildMHDbox():
    """
    MHD input parameters
    """
    return html.Div(
            #id="MHDbox",
            className="column1",
            children=[
                dbc.Label(id="shotLabel", children="Shot Number  "),
                dbc.Input(id="shot"),
                dbc.Label("Trace Distance [degrees]"),
                dbc.Input(id="traceLength"),
                dbc.Label("Toroidal Step Size [degrees]"),
                dbc.Input(id="dpinit"),
                dbc.Label("Multiplier for Psi (RZ, sep, axis) Values"),
                dbc.Input(id="psiMult1"),
                dbc.Label("Multiplier for Bt (Bt0, Fpol) Values"),
                dbc.Input(id="BtMult1"),
                dbc.Label("Multiplier for Ip Values"),
                dbc.Input(id="IpMult1"),
                html.Br(),
                dbc.Label("Load MHD Equilibria:"),
                dcc.Upload(
                    className="PFCupload",
                    id='gfiletable-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select GEQDSKs')
                    ]),
                    style={
                        'width': '60%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed',
                        'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                        'justify-content':'center',
                        },
                    multiple=True,
                    ),
                html.Div(id="hiddenDivGfileUpload"),
                html.Br(),
                dbc.Label("Select Plasma Type:"),
                dbc.RadioItems(
                                id="plasma3Dmask",
                                options=[
                                    {'label': '2D Plasmas', 'value': 'plasma2D'},
                                    {'label': '3D Plasmas', 'value': 'plasma3D'},
                                        ],
                                        value='plasma2D'
                                ),
                html.Br(),
                html.Br(),
                dbc.Button("Load MHD", id="loadMHD", n_clicks=0, color="secondary"),
                html.Br(),
            ],
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
              State('traceLength', 'value'),
              State('dpinit', 'value'),
              State('gfiletable-upload', 'filename'),
              State('gfiletable-upload', 'contents'),
              State('plasma3Dmask', 'value'),
              State('dataPath', 'value'),
              State('MachFlag', 'value'),
              State('psiMult1', 'value'),
              State('BtMult1', 'value'),
              State('IpMult1', 'value')]
              )
def loadMHD(n_clicks,shot,traceLength,dpinit,eqList,eqData,plasma3Dmask,dataPath,MachFlag,
            psiMult, BtMult, IpMult):
    """
    Load MHD
    """
    if MachFlag == None:
        print("Select a machine before loading MHD parameters")
        raise PreventUpdate

    #if data directory doesn't exist, create it
    tools.makeDir(dataPath, clobberFlag=False, mode=chmod, UID=UID, GID=GID)

    if plasma3Dmask == 'plasma3D': plasma3Dmask = True
    else:  plasma3Dmask = False

    if (shot is None) and (eqList is None):
        raise PreventUpdate

    if shot is not None: shot = int(shot)
    if traceLength is not None:
        traceLength = float(traceLength)
    if dpinit is not None:
        dpinit = float(dpinit)
    else:
        dpinit = 1.0
    if eqList is not None:
        if type(eqList) is not list:
            eqList = [eqList]
    if eqData is not None:
        if type(eqData) is not list:
            eqData = [eqData]

    gui.getMHDInputsForGUI(shot=shot,
                     traceLength=traceLength,
                     dpinit=dpinit,
                     eqList=eqList,
                     eqData=eqData,
                     plasma3Dmask=plasma3Dmask,
                     psiMult=psiMult,
                     BtMult=BtMult,
                     IpMult=IpMult,
                    )

    ts = gui.MHD.timesteps
    tminMHD = ts.min()
    tmaxMHD = ts.max()
    if eqList is None:
        tAll = np.linspace(int(tminMHD), int(tmaxMHD), (int(tmaxMHD)-int(tminMHD)+1))
        data = [dict([{'filename':'', 'timestep':''}][0])]
    else:
        tAll = ts
        keys=["filename"]
        interpData = pd.DataFrame(eqList, columns=keys)
        interpData["timestep[ms]"] = ""
        data = interpData.to_dict('records')
        #interpData = [dict([{'filename':g, 'timestep':''} for g in eqList])]
    marks = {}
    tick_interval = int(len(ts)/10) #max 10 tick marks
    for tIdx,t in enumerate(ts):
        if t in tAll:
            #display all tick marks
            #marks.update({t:'{}'.format(t)})
            #hide all tick marks
            #marks.update({t:'{}'.format("")})
            #only display some tick marks
            if len(ts) > 10:
                marks.update({t: str(t) if tIdx % tick_interval == 0 else ""})
            else:
                marks.update({t:'{}'.format(t)})
    value = ts[0]

    MHDdata = {
        'Shot Number':shot,
        'Trace Distance [degrees]':traceLength,
        'Toroidal Step Size [degrees]':dpinit,
        '3D Plasma? [0=False]':plasma3Dmask,
        'Psi Multiplier':psiMult,
        'Bt Multiplier':BtMult,
        'Ip Multiplier':IpMult,
        }

    return tminMHD, tmaxMHD, marks, value, data, MHDdata

#Load gfile
@app.callback([Output('hiddenDivGfileUpload', 'children')],
              [Input('gfiletable-upload', 'filename')],
              [State('MachFlag', 'value')])
def gfileUpload(gFile, MachFlag):
    if MachFlag is None:
        raise PreventUpdate
    else:
        if len(gFile) > 1:
            div = dbc.Label("Loaded multiple gfiles", className="text-success")
        else:
            div = dbc.Label("Loaded gFile: "+gFile[0], className="text-success")
        return [div]

#==========CAD==========
def buildCADbox():
    """
    CAD input parameters
    """
    return html.Div(
            id="CADbox",
            draggable='yes',
            children=[
                dbc.Label(id="gridResLabel", children="Intersect Resolution [mm]"),
                dbc.Input(id="gridRes", className="textInput"),

                dbc.Label("Global Translation in X direction [mm]"),
                dbc.Input(id="gTx", className="textInput"),
                dbc.Label("Global Translation in Y direction [mm]"),
                dbc.Input(id="gTy", className="textInput"),
                dbc.Label("Global Translation in Z direction [mm]"),
                dbc.Input(id="gTz", className="textInput"),

                dbc.Button("Load CAD Settings", id="loadRes", style={'margin':'10px 10px 10px 10px'}),
                html.Div(id="hiddenDivCAD1"),
                dbc.Label(id="STPdropLabel", children="STP File Direct Upload:"),
                dcc.Upload(
                    className="PFCupload",
                    id='CAD-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select CAD file')
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
               State('gTx', 'value'),
               State('gTy', 'value'),
               State('gTz', 'value'),
               State('MachFlag', 'value')])
def loadRes(n_clicks, gridRes, gTx, gTy, gTz, MachFlag):
    if n_clicks is None:
        raise PreventUpdate
    if MachFlag is None:
        return [dbc.Label("Select a machine", style={'color':'#fc0313'})]
    gui.getCADinputs(gridRes, gTx, gTy, gTz, mode='gui')

    CADdata = {
        'Intersect Max Edge Length [mm]':gridRes
        }

    return [dbc.Label("Loaded Resolution Settings", className="text-success"), CADdata]

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
        return [dbc.Label("Select a machine first", style={'color':'#fc0313'})]
    else:
        contents = STPcontents[0]
        content_type, content_string = contents.split(',')
        STPdata= base64.b64decode(content_string)
        gui.getCAD(STPfile=STPfile[0],STPdata=STPdata, ts=ts[0])
    return [dbc.Label("Loaded CAD: "+STPfile[0], className="text-success")]

#==========HF==========
def buildHFbox():
    """
    Heat Flux input parameters
    """
    return html.Div(
            id="HFbox",
            children=[
                dbc.Label(id="hfModeLabel", children="Select a Heat Flux Profile"),
                dcc.Dropdown(
                id='hfMode',
                options=[
                    {'label': 'Gaussian Spreading', 'value': 'eich'},
                    {'label': 'Multi-Exponential', 'value': 'multiExp'},
                    {'label': 'Limiter', 'value': 'limiter'},
                    {'label': 'Top Hat', 'value': 'tophat'},
                    {'label': 'From qFiles', 'value': 'qFile'},
                    {'label': 'From rzq profile', 'value': 'rzqprofile'}
                    ],
                value=None,
                ),
                html.Div(id='hfParameters',
                         children=[
                                    loadHFSettings(mode=None,hidden=True),
                                   ],
                        ),
                html.Br(),
                dbc.Button("Load HF", id="loadHF", style={'margin':'10px 10px 10px 10px'}),
                html.Div(id="hiddenDivHF"),

            ],
        )


#Heat Flux Callbacks and Conditionals
@app.callback([Output('hfParameters', 'children')],
              [Input('hfMode', 'value')],
              [State('session', 'data')]
              )
def hfParameters(mode, sessionData):
    div = [loadHFSettings(mode=mode, hidden=False, sessionData=sessionData)]
    return [div]

def PsolInput(hidden=False, data=None):
    #do this hidden business so that we can always load defaults into these id's
    if hidden==True:
        className="hiddenBox"
    else:
        className="visibleBox"


    row2 = html.Div(
                className="rowBox",
                children=[
                    html.Div(
                        className="colBox",
                        children=[
                            dbc.Label("Upper Inner Power Fraction", className="hfLabel"),
                            dbc.Input(id="fracUI", value=data['fracUI']),
                            ],
                    ),
                    html.Div(
                        className="colBox",
                        children=[
                            dbc.Label("Upper Outer Power Fraction", className="hfLabel"),
                            dbc.Input(id="fracUO", value=data['fracUO']),
                            ],
                    ),
                    ],
            )

    row3 = html.Div(
                className="rowBox",
                children=[
                    html.Div(
                        className="colBox",
                        children=[
                            dbc.Label("Lower Inner Power Fraction", className="hfLabel"),
                            dbc.Input(id="fracLI", value=data['fracLI']),
                            ],
                    ),
                    html.Div(
                        className="colBox",
                        children=[
                            dbc.Label("Lower Outer Power Fraction", className="hfLabel"),
                            dbc.Input(id="fracLO", value=data['fracLO']),
                            ],
                    ),
                    ],
            )


    return html.Div(
             className=className,
             children=[
                    dbc.Label("Source Power (PSOL or Psep) [MW]", className="psolInput"),
                    dbc.Input(id="P", value=data['P'], className="psolInput"),
                    dbc.Label("Fraction of Source Power Radiated by Photons", className="psolInput"),
                    dbc.Input(id="radFrac", value=data['radFrac'], className="psolInput"),
                    row2,
                    row3,
                    ],
                    )

def loadHFSettings(mode=None, hidden=False, sessionData=None):
    #do this hidden business so that we can always load defaults into these id's
    #hideMask corresponds to the following parameters:
    #[eichProfile, commonRegion, privateRegion]
    hideMask = ['hiddenBox','hiddenBox','hiddenBox','hiddenBox', 'hiddenBox', 'hiddenBox']
    if mode=='eich':
        hideMask = ['hfInput','hiddenBox','hiddenBox','hiddenBox','hiddenBox', 'hiddenBox']
    elif mode=='limiter':
        hideMask = ['hiddenBox','hiddenBox','hfInput','hiddenBox','hiddenBox', 'hiddenBox'] #common flux region
    elif mode=='multiExp':
        hideMask = ['hiddenBox','hfInput','hiddenBox','hiddenBox','hiddenBox','hiddenBox'] #common + private flux region
    elif mode=='qFile':
        hideMask = ['hiddenBox','hiddenBox','hiddenBox','hfInput','hiddenBox','hiddenBox'] #common + private flux region
    elif mode=='tophat':
        hideMask = ['hiddenBox','hiddenBox','hiddenBox', 'hiddenBox','hfInput','hiddenBox'] #common + private flux region
    elif mode=='rzqprofile': 
        hideMask = ['hiddenBox','hiddenBox','hiddenBox','hiddenBox','hiddenBox', 'hfInput'] #common + private flux region
    if hidden==True or mode==None:
        hideMask=['hiddenBox','hiddenBox','hiddenBox','hiddenBox','hiddenBox','hiddenBox']


    #if we are using qFile, hide Psol inputs
    if mode == 'qFile':
        hidden=True

    #load session data that was provided in input file
    if sessionData == None:
        sessionData = {}
        sessionData['hfMode'] = None
        sessionData['lqCN'] = None
        sessionData['lqCF'] = None
        sessionData['lqPN'] = None
        sessionData['lqPF'] = None
        sessionData['lqCNmode'] = None
        sessionData['lqCFmode'] = None
        sessionData['lqPNmode'] = None
        sessionData['lqPFmode'] = None
        sessionData['S'] = None
        sessionData['SMode'] = None
        sessionData['fracCN'] = None
        sessionData['fracCF'] = None
        sessionData['fracPN'] = None
        sessionData['fracPF'] = None
        sessionData['fracUI'] = None
        sessionData['fracUO'] = None
        sessionData['fracLI'] = None
        sessionData['fracLO'] = None
        sessionData['P'] = None
        sessionData['radFrac'] = None
        sessionData['qBG'] = None
        sessionData['fG'] = None
        sessionData['qFilePath'] = None
        sessionData['qFileTag'] = None
        sessionData['rzqFile'] = None


    return html.Div(
            children=[
                #gaussian spreading / eich
                html.Div(
                    #className=hideMask[0],
                    children=[
                        eichParameters(hideMask[0], sessionData),
                        ]
                ),
                #multiple exponentials
                html.Div(
                    #className=hideMask[1],
                    children=[
                        multiExpParameters(hideMask[1], sessionData),
                        ]
                ),
                #limiters
                html.Div(
                    #className=hideMask[2],
                    children=[
                        limiterParameters(hideMask[2], sessionData),
                        ]
                ),
                #from file
                html.Div(
                    #className=hideMask[2],
                    children=[
                        qFileParameters(hideMask[3], sessionData),
                        ]
                ),
                #tophat
                html.Div(
                    #className=hideMask[0],
                    children=[
                        tophatParameters(hideMask[4], sessionData),
                        ]
                ),
                # rzqFile
                html.Div(
                    children=[
                        rzqprofileParameters(hideMask[5]),
                    ]
                ),
                PsolInput(hidden, sessionData),
                    ],
                    )


def eichParameters(className, data):
    row1 = html.Div(
    className='rowBox',
        children=[
            html.Div(
            className='colBox',
            children=[
                dbc.Label("Select Heat Flux Width source:", className="hfLabel"),
                dcc.Dropdown(
                id='eichlqCNmode',
                className="SelectorBoxInput",
                options=[
                    {'label': 'Eich #15', 'value': 'eich'},
                    {'label': 'User Defined', 'value': 'user'}
                    ],
                value=data['lqCNmode'],
                ),
                ],
            ),
            html.Div(
            className='colBox',
            children=[
                dbc.Label("Select Gaussian Spreading source:", className="hfLabel"),
                dcc.Dropdown(
                id='eichSMode',
                className="SelectorBoxInput",
                options=[
                    {'label': 'From Makowski Scaling', 'value': 'makowski'},
                    {'label': 'User Defined', 'value': 'user'}
                    ],
                value=data['SMode'],
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
                dbc.Label("User Defined Heat Flux Width [mm]:", className="hfLabel"),
                dbc.Input(id="lqEich", value=data['lqCN'], className="hfInput2"),

                ],
            ),
            html.Div(
            className="colBox",
            children=[
                dbc.Label("User Defined Gaussian Spreading [mm]:", className="hfLabel"),
                dbc.Input(id="S", value=data['S'], className="hfInput2"),
                ],
            ),
            ])

    row3 = html.Div(
        className='rowBox',
        children=[
            html.Div(
                className="colBox",
                children=[
                    dbc.Label("Background Heat Flux [MW/m^2]", className="hfLabel"),
                    dbc.Input(id="qBG", value=data['qBG'], className="hfInput2"),
                ]),
            html.Div(
                className="colBox",
                children=[
                    dbc.Label("Greenwald Density Fraction", className="hfLabel"),
                    dbc.Input(id="fG", value=data['fG'], className="hfInput2"),
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


def multiExpParameters(className, data):
    row1 = html.Div(
        className='rowBox',
        children = [
            html.Div(
                className='colBox',
                children=[
                    dbc.Label("Select Common Near Heat Flux Width source:"),
                    dcc.Dropdown(
                    className="SelectorBoxInput",
                    id='multiExplqCNmode',
                    options=[
                        {'label': 'Eich #15', 'value': 'eich'},
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=data['lqCNmode'],
                    )
                ]),
            html.Div(
                className='colBox',
                children=[
                    dbc.Label("Select Common Far Heat Flux Width source:"),
                    dcc.Dropdown(
                    className="SelectorBoxInput",
                    id='multiExplqCFmode',
                    options=[
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=data['lqCFmode'],
                    )
                ]),
            ]
    )
    row2 = html.Div(
        className='rowBox',
        children = [
            html.Div(
                className='colBox',
                children=[
                    dbc.Label("Select Private Near Heat Flux Width source:"),
                    dcc.Dropdown(
                    className="SelectorBoxInput",
                    id='multiExplqPNmode',
                    options=[
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=data['lqPNmode'],
                    )
                ]),
            html.Div(
                className='colBox',
                children=[
                    dbc.Label("Select Private Far Heat Flux Width source:"),
                    dcc.Dropdown(
                    className="SelectorBoxInput",
                    id='multiExplqPFmode',
                    options=[
                        #{'label': 'From Brunner Scaling', 'value': 'brunner'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=data['lqPFmode'],
                    )
                ]),
            ]
    )
    row3 = html.Div(
        className='rowBox',
            children=[
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Common Near Heat Flux Width [mm]"),
                        dbc.Input(id="lqCN", value=data['lqCN']),
                    ]),
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Common Near Power Fraction"),
                        dbc.Input(id="fracCN", value = data['fracCN']),
                    ]),
            ])

    row4 = html.Div(
            className='rowBox',
            children=[
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Common Far Heat Flux Width [mm]"),
                        dbc.Input(id="lqCF", value=data['lqCF']),
                    ]),
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Common Far Power Fraction"),
                        dbc.Input(id="fracCF", value=data['fracCF']),
                    ]),
            ])

    row5 = html.Div(
            className='rowBox',
            children=[
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Private Near Heat Flux Width [mm]"),
                        dbc.Input(id="lqPN", value=data['lqPN']),
                    ]),
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Private Near Power Fraction"),
                        dbc.Input(id="fracPN", value=data['fracPN']),
                    ]),
            ])

    row6 = html.Div(
            className='rowBox',
            children=[
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Private Far Heat Flux Width [mm]"),
                        dbc.Input(id="lqPF", value=data['lqPF']),
                    ]),
                html.Div(
                    className='colBox',
                    children=[
                        dbc.Label("Private Far Power Fraction"),
                        dbc.Input(id="fracPF", value=data['fracPF']),
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


def limiterParameters(className, data):
    #return if this div is supposed to be hidden to prevent duplicate IDs
    #if className== 'hiddenBox':
    #    return

    row1 = html.Div(
        className="rowBox",
        children = [
            html.Div(
                className="colBox",
                children=[
                    dbc.Label("Select Common Near Heat Flux Width source:"),
                    dcc.Dropdown(
                    className="SelectorBoxInput",
                    id='limiterlqCNmode',
                    options=[
                        {'label': 'Eich #15', 'value': 'eich'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=data['lqCNmode'],
                    )
                ]),
            html.Div(
                className="colBox",
                children=[
                    dbc.Label("Select Common Far Heat Flux Width source:"),
                    dcc.Dropdown(
                    className="SelectorBoxInput",
                    id='limiterlqCFmode',
                    options=[
                        {'label': 'Horacek Fig. 6a', 'value': 'horacek'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=data['lqCFmode'],
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
                        dbc.Label("Common Near Heat Flux Width [mm]"),
                        dbc.Input(id="limlqCN", value=data['lqCN']),
                    ]),
                html.Div(
                    className="colBox",
                    children=[
                        dbc.Label("Common Far Heat Flux Width [mm]"),
                        dbc.Input(id="limlqCF", value=data['lqCF']),
                    ]),
            ])

    row3 = html.Div(
            className="rowBox",
            children=[
                html.Div(
                    className="colBox",
                    children=[
                        dbc.Label("Common Near Power Fraction"),
                        dbc.Input(id="limfracCN", value=data['fracCN']),
                    ]),

                html.Div(
                    className="colBox",
                    children=[
                        dbc.Label("Common Far Power Fraction"),
                        dbc.Input(id="limfracCF", value=data['fracCF']),
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


def  qFileParameters(className, data):
    #return if this div is supposed to be hidden to prevent duplicate IDs
    #if className== 'hiddenBox':
    #    return

    row1 = html.Div(
        className="rowBox",
        children = [
            html.Div(
            className="colBox",
                children=[
                    dbc.Label("Input qFile shot path:"),
                    dbc.Input(id="qFilePath", value=data['qFilePath']),
                ]),
            html.Div(
                className="colBox",
                children=[
                    dbc.Label("Input qFile tag (ie HF_optical.csv):"),
                    dbc.Input(id="qFileTag", value=data['qFileTag']),
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


def  tophatParameters(className, data):
    #return if this div is supposed to be hidden to prevent duplicate IDs
    #if className== 'hiddenBox':
    #    return

    row1 = html.Div(
        className="rowBox",
        children = [
            html.Div(
            className="colBox",
                children=[
                    dbc.Label("Heat Flux Width [mm]:"),
                    dbc.Input(id="lqTopHat", value=data['lqCN']),
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

def rzqprofileParameters(className):
    row1 = html.Div(
            className="colBox",
                id="rzqprofilebox",
                children=[
                    dcc.Upload(
                        className="rowBox",
                        id='rzqprofile-upload',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select rzq profile Files')
                        ]),
                        style={
                            'width': '100%', 'height': '60px', 'lineHeight': '60px',
                            'borderWidth': '1px', 'borderStyle': 'dashed',
                            'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                            'align-items':'left'
                            },
                        multiple=True,
                        ),
                html.Br(),
                html.Div(id="hiddenDivrzqprofileUpload"),
                html.Br(),

        ],
        )
    div = html.Div(
        className=className,
        children=[
            row1,
            ]
        )
    return div

@app.callback([Output('hiddenDivrzqprofileUpload', 'children')],
              [Input('rzqprofile-upload', 'filename')],
              [State('MachFlag', 'value')])

def rzqprofileUpload(rzqFile, MachFlag):
    if MachFlag is None:
        raise PreventUpdate
    if rzqFile is None:
        raise PreventUpdate
    else:
        return[html.Label("Loaded rzq profile: "+ rzqFile[0], className="text-success")]


#def commonRegionParameters():
#    """
#    near and far heat flux widths and power sharing fractions
#    """
#    row1 = html.Div(
#            className="rowBox",
#            children=[
#                html.Div(
#                    className="colBox",
#                    children=[
#                        html.Label("Common Near Heat Flux Width [mm]"),
#                        dcc.Input(id="lqCN", className="hfInput2"),
#                    ]),
#                html.Div(
#                    className="colBox",
#                    children=[
#                        html.Label("Common Near Power Fraction"),
#                        dcc.Input(id="fracCN", className="hfInput2"),
#                    ]),
#            ])
#
#    row2 = html.Div(
#            className="rowBox",
#            children=[
#                html.Div(
#                    className="colBox",
#                    children=[
#                        html.Label("Common Far Heat Flux Width [mm]"),
#                        dcc.Input(id="lqCF", className="hfInput2"),
#                    ]),
#                html.Div(
#                    className="colBox",
#                    children=[
#                        html.Label("Common Far Power Fraction"),
#                        dcc.Input(id="fracCF", className="hfInput2"),
#                    ]),
#            ])
#
#    return row1,row2


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
               State('lqTopHat', 'value'),
               State('fracCN', 'value'),
               State('fracCF', 'value'),
               State('fracPN', 'value'),
               State('fracPF', 'value'),
               State('fracUI', 'value'),
               State('fracUO', 'value'),
               State('fracLI', 'value'),
               State('fracLO', 'value'),
               State('eichlqCNmode', 'value'),
               State('eichSMode', 'value'),
               State('multiExplqCNmode', 'value'),
               State('multiExplqCFmode', 'value'),
               State('multiExplqPNmode', 'value'),
               State('multiExplqPFmode', 'value'),
               State('limiterlqCNmode', 'value'),
               State('limiterlqCFmode', 'value'),
               State('limlqCN', 'value'),
               State('limlqCF', 'value'),
               State('limfracCN', 'value'),
               State('limfracCF', 'value'),
               State('qBG', 'value'),
               State('P', 'value'),
               State('radFrac', 'value'),
               State('fG', 'value'),
               State('qFilePath', 'value'),
               State('qFileTag', 'value'),
               State('rzqprofile-upload', 'filename'),
               State('rzqprofile-upload', 'contents')
               ])
def loadHF(n_clicks,hfMode,MachFlag,
            lqEich,S,lqCN,lqCF,lqPN,lqPF,lqTopHat,
            fracCN,fracCF,fracPN,fracPF,
            fracUI,fracUO,fracLI,fracLO,
            eichlqCNmode,SMode,
            multiExplqCNmode,multiExplqCFmode,multiExplqPNmode,multiExplqPFmode,
            limiterlqCNmode,limiterlqCFmode,limlqCN,limlqCF,limfracCN,limfracCF,
            qBG,P,radFrac,fG,
            qFilePath, qFileTag,
            rzqFile, rzqFiledata):
    if MachFlag is None:
        raise PreventUpdate
    else:
        #set up the heat flux configuration (which scalings to use)
        if hfMode == 'limiter':
            lqCNmode = limiterlqCNmode
            lqCFmode = limiterlqCFmode
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
            lqCNmode = multiExplqCNmode
            lqCFmode = multiExplqCFmode
            lqPNmode = multiExplqPNmode
            lqPFmode = multiExplqPFmode
            lqCN = lqCN #unnecessary, but here for clarity
            lqCF = lqCF #unnecessary, but here for clarity
            lqPN = lqPN #unnecessary, but here for clarity
            lqPF = lqPF #unnecessary, but here for clarity
            SMode = None
            qFileTag = None
            qFilePath = None
        elif hfMode == 'qFile':
            lqCNmode = eichlqCNmode
            lqCFmode = None
            lqPNmode = None
            lqPFmode = None
            lqCN = lqEich
            lqCF = 0.0
            lqPN = 0.0
            lqPF = 0.0
            SMode = SMode
        elif hfMode == 'tophat':
            lqCNmode = 'user'
            lqCFmode = None
            lqPNmode = None
            lqPFmode = None
            lqCN = lqTopHat
            lqCF = 0.0
            lqPN = 0.0
            lqPF = 0.0          
            Smode = None 
            qFileTag = None
            qFilePath = None
        elif hfMode == 'rzqprofile': 
            lqCNmode = None
            lqCFmode = None
            lqPNmode = None
            lqPFmode = None
            lqCN = 0.0
            lqCF = 0.0
            lqPN = 0.0
            lqPF = 0.0
            SMode = None
            qFileTag = None
            qFilePath = None
            #only support one rzqFile
            rzqFile = rzqFile[0]
            rzqFiledata = rzqFiledata[0]
        else: #eich mode is default
            lqCNmode = eichlqCNmode
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


        gui.getHFInputs(hfMode,
                        lqCN,lqCF,lqPN,lqPF,S,
                        fracCN,fracCF,fracPN,fracPF,
                        fracUI,fracUO,fracLI,fracLO,
                        lqCNmode,lqCFmode,lqPNmode,lqPFmode,SMode,
                        qBG,P,radFrac,fG,
                        qFilePath,qFileTag,
                        rzqFile, rzqFiledata)


        #Update output tab table
        HFdata = gui.HF.HFdataDict
        #hfDict = [{'Parameter':i, 'Value':dataDict[i]} for i in inputTableData]

    return [html.Label("Loaded HF Settings", className="text-success"), HFdata]



#==========PFC==========
def buildPFCbox():
    return html.Div(
            id="PFCbox",
            children=[
                dbc.Label("Load a PFC file:"),
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
                html.Br(),
                html.Div(children=loadPFCtable(), className="PFCtable"),
                html.Br(),
                dbc.Button("Load PFC Settings", id="loadPFC", n_clicks=0, style={'margin':'0 10px 10px 0'}),
#                html.Button('Add Row', id='add-rows-button', n_clicks=0, style={'margin':'0 10px 10px 0'}),
                dbc.Button("Download Default PFC file", id="downloadPFCbtn", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                dcc.Download(id="downloadPFC"),
                html.Div(id="hiddenDivDownloadPFC", style={"display": "none"}),
                html.Div(id="hiddenDivPFC", style={"display": "none"}),
                dbc.Input(id="hiddenDivPFC2", style={"display": "none"}),
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
        style_header={'backgroundColor': 'rgb(30, 30, 30)', 'color':'white'},
        style_cell={
            'textAlign': 'left',
            #'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'black' #text color
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
        hiddenDiv = [html.Label("Select a Tokamak and Load MHD/CAD First", className="text-warning")]

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

        hiddenDiv = [html.Label("Loaded PFC Data into HEAT.  Meshing complete.", className="text-success")]

    #file dropper
    elif filename != dataStore['PFCfilename']:
        df = parse_contents(uploadContents, filename)
        tableData = df.to_dict('records')
        tableColumns = [{"name": i, "id": i} for i in df.columns]
        hiddenDiv = [html.Label("Loaded file: "+filename+", creating meshes now...", className="text-info")]
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
    return [dcc.send_file(gui.tmpDir + "PFCinput.csv"),
            html.Label("Saved PFC Default File", className="text-success")]






#==========gyro orbits==========
def buildGYRObox():
    return html.Div(
        id="gyrobox",
        children=[
            gyroInputBoxes(),
            dbc.Label("Velocity Method"),
            dcc.Dropdown(
                id='vMode',
                className="wideSelect",
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
            dbc.Label("Gyro Orbit Power Source"),
            dcc.Dropdown(
                id='sourceMode',
                className="wideSelect",
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
            dbc.Label("(Note: Sources are OMITTED from Temp Calculation)"),
            html.Br(),
            dbc.Button("Load Gyro Settings", id="loadGYRO", n_clicks=0, style={'margin':'0 10px 10px 0'}),
            html.Div(id="hiddenDivGyro")
        ],
        className="HFbox",
    )


def gyroInputBoxes():
    return html.Div(
            children=[
            html.Div(
                children=[
                    dbc.Label("# steps per helix"),
                    dbc.Input(id="N_gyroSteps", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    dbc.Label("# gyroPhase Angles"),
                    dbc.Input(id="N_gyroPhase", className="textInput"),
                ],
                className="OFInput"
            ),

            html.Div(
                children=[
                    dbc.Label("Gyro Trace Length [deg]"),
                    dbc.Input(id="gyroTraceLength", className="textInput"),
                ],
                className="OFInput"
            ),
            html.Div(
                children=[
                    dbc.Label("Ion Power Fraction [0-1]"),
                    dbc.Input(id="ionFrac", className="textInput"),
                ],
                className="OFInput"
            ),
            html.Div(
                children=[
                    dbc.Label("Ion Mass [AMU]"),
                    dbc.Input(id="ionMassAMU", className="textInput"),
                ],
                className="OFInput"
            ),
            ],
            className="wideBoxNoColor",
            )

@app.callback([Output('gyroVparams', 'children')],
              [Input('vMode', 'value')])
def gyroParameters(mode):
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
                            dbc.Label("Average Temperature [eV]"),
                            dbc.Input(id="gyroT_eV", className="textInput"),
                            ],
                        className="OFInput"
                        ),

                    html.Div(
                        children=[
                            dbc.Label("# Velocity Slices"),
                            dbc.Input(id="N_vSlice", className="textInput"),
                            ],
                        className="OFInput"
                        ),

                    html.Div(
                        children=[
                            dbc.Label("# Velocity Phases"),
                            dbc.Input(id="N_vPhase", className="textInput"),
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
            dbc.Label("3D plasma temperature interpolation not yet available", className="text-warning")
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
                        dbc.Checklist(
                            id='gyroSource',
#                            className="checkListBox",
                            options=options,
                            value=['allROI'],
                        ),
                    ],
#                    className="wideBoxNoColor"
                    )
                ],
            className=className,
            )
    return div

def allROISource(className):
    return html.Div(
        children=[
            dbc.Label("Using all ROI PFCs as gyro power sources", className="text-success")
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
               State('gyroTraceLength', 'value'),
               State('ionMassAMU', 'value'),
               State('vMode','value'),
               State('gyroT_eV', 'value'),
               State('N_vPhase', 'value'),
               State('N_vSlice', 'value'),
               State('ionFrac', 'value'),
               State('gyroSource', 'value')
              ])
def loadGYRO(n_clicks,N_gyroSteps,N_gyroPhase,gyroTraceLength,ionMassAMU,vMode,gyroT_eV,
             N_vPhase, N_vSlice, ionFrac, gyroSource):
    """
    sets up GYRO module
    """
    if n_clicks == 0:
        raise PreventUpdate

    #handle gyro orbit dropdown (user can only select all or custom PFCs)
    if ('allROI' in gyroSource) and (len(gyroSource)>1):
        gyroSource.remove('allROI')

    gui.getGyroInputs(N_gyroSteps,N_gyroPhase,gyroTraceLength,ionMassAMU,vMode,gyroT_eV,
                      N_vPhase, N_vSlice, ionFrac, gyroSource)
    gyroPhaseFig = gyroPhasePlots(update=True)
    vPhaseFig = vPhasePlots(update=True)
    vSliceFig = vSlicePlots(update=True)
    cdfSliceFig = cdfSlicePlots(update=True)



    GYROdata = {
            'Number of steps per helix period':N_gyroSteps,
            'Number of samples in gyro phase space':N_gyroPhase,
            'Gyro trace length [degrees] (goes both directions)': gyroTraceLength,
            'Ion effective mass [AMU]': ionMassAMU,
            'Velocity / Temperature Mode':vMode,
            'Ion temperature at PFC surface [eV]':gyroT_eV,
            'Number of samples in velocity phase space':N_vPhase,
            'Number of samples from Maxwellian velocity distribution':N_vSlice,
            'Fraction of power carried by ions':ionFrac,
            'Source for gyro orbit power':':'.join(gyroSource),
            }

    return [dbc.Label("Loaded Gyro Orbit Data into HEAT", className="text-success"),
            gyroPhaseFig,
            vPhaseFig,
            vSliceFig,
            cdfSliceFig,
            GYROdata]




#==========RAD==========
def buildRADbox():
    """
    MHD input parameters
    """
    return html.Div(
            id="RADbox",
            children=[
                dcc.Upload(
                    className="PFCupload",
                    id='RAD-upload',
                    children=html.Div([
                        'Drag and Drop or ',
                        html.A('Select radFile')
                    ]),
                    style={
                        'width': '60%', 'height': '60px', 'lineHeight': '60px',
                        'borderWidth': '1px', 'borderStyle': 'dashed',
                        'borderRadius': '5px', 'textAlign': 'center', 'margin': '10px',
                        },
                    multiple=True,
                    ),
                html.Div(id="hiddenDivRadUpload"),
                html.Div(
                    children=[
                    html.Div(
                        children=[
                            dbc.Label("Minimum phi of source [deg] "),
                            dbc.Input(id="phiMin", className="textInput"),
                        ],
                        className="OFInput",
                    ),
                    html.Div(
                        children=[
                            dbc.Label("Maximum phi of source [deg] "),
                            dbc.Input(id="phiMax", className="textInput"),
                        ],
                        className="OFInput",
                    ),
                    html.Div(
                        children=[
                            dbc.Label("Number of toroidal repetitions  "),
                            dbc.Input(id="Ntor", className="textInput"),
                        ],
                    className="OFInput",
                    ),
                    html.Div(
                        children=[
                            dbc.Label("Number of reflections (currently disabled)"),
                            dbc.Input(id="Nref", className="textInput", disabled=True),
                        ],
                        className="OFInput",
                    ),
                    ],
                    className="wideBoxNoColor"
                    ),

                html.Br(),
                dbc.Button("Load RAD Settings", id="loadRAD", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
                html.Div(id="hiddenDivRAD"),
            ],
            className="HFbox",
        )




#Load RAD button connect
@app.callback([Output('hiddenDivRAD', 'children'),
               Output('RADDataStorage', 'data')
               ],
              [Input('loadRAD', 'n_clicks')],
              [State('phiMin', 'value'),
               State('phiMax', 'value'),
               State('Ntor', 'value'),
               State('Nref', 'value'),
               State('RAD-upload', 'filename'),
               State('RAD-upload', 'contents'),
              ])
def loadRAD(n_clicks,phiMin,phiMax,Ntor,Nref,radFile,radData):
    """
    sets up RAD module
    """
    if n_clicks == 0:
        raise PreventUpdate

    if radFile == None:
        outDiv = dbc.Label("Load radFile before loading settings...", style={'color':'#f50707'})
        RADdata = {}
    else:
        outDiv = dbc.Label("Loaded RAD Data into HEAT", className="text-success")
        Ntor = int(Ntor)
        Nref = int(Nref)
        phiMin = float(phiMin)
        phiMax = float(phiMax)
        #these values are hard coded for GUI:
        Prad_mult=1.0
        rayTracer='open3d'
        saveRadFrac=False

        gui.getRADInputs(radFile[0], Ntor, Nref, phiMin, phiMax, rayTracer, Prad_mult, saveRadFrac, radData[0])


        RADdata = {
            'Minimum phi of radiation source [deg]':phiMin,
            'Maximum phi of radiation source [deg]':phiMax,
            'Number of toroidal repetitions of radiation source':Ntor,
            'Number of photon reflections to trace':Nref,
            'Raytracer to use for photon calc':rayTracer,
            'Multiplier for Prad in R,Z,P file':Prad_mult,
            'Should we save map between emission source and targets?':saveRadFrac,
            }
    return [outDiv, RADdata]

#Load radFile
@app.callback([Output('hiddenDivRadUpload', 'children')],
              [Input('RAD-upload', 'filename')],
              [State('MachFlag', 'value')])
def gfileUpload(radFile, MachFlag):
    if MachFlag is None:
        raise PreventUpdate
    else:
        return [dbc.Label("Loaded radFile: "+radFile[0], className="text-success")]



#==========openFOAM==========
def buildOFbox():
    return html.Div(
        id="OFbox",
        children=[
            OFinputBoxes(),
            html.Br(),
            dbc.Button("Load OF Settings", id="loadOF", n_clicks=0, style={'margin':'0 10px 10px 0'}),
            html.Div(id="hiddenDivOF")
        ],
#        className="HFbox",
    )

def OFinputBoxes():
    return html.Div(
            children=[
            html.Div(
                children=[
                    dbc.Label("Start Time [s]"),
                    dbc.Input(id="OFstartTime", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    dbc.Label("Stop Time [s]"),
                    dbc.Input(id="OFstopTime", className="textInput"),
                ],
                className="OFInput"
            ),
            html.Div(
                children=[
                    dbc.Label("Minimum Mesh Refinement Level"),
                    dbc.Input(id="OFminMeshLev", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    dbc.Label("Maximum Mesh Refinement Level"),
                    dbc.Input(id="OFmaxMeshLev", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    dbc.Label("STL scaling"),
                    dbc.Input(id="OFSTLscale", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    dbc.Label("Timestep Size [s]"),
                    dbc.Input(id="OFdeltaT", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    dbc.Label(" Write Timestep [s] (disabled)"),
                    dbc.Input(id="OFwriteDeltaT", className="textInput", value="NA", disabled=True),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    dbc.Label("Material Selection"),
                    dcc.Dropdown(
                        id='materialSelect',
                        className="SelectorBoxInput",
                        options=[
                            {'label': 'SGLR6510 Graphite', 'value': 'SGL'},
                            {'label': 'ATJ Graphite', 'value': 'ATJ'},
                            {'label': 'Molybdenum', 'value': 'MOLY'},
                            {'label': 'Tungsten', 'value': 'TUNG'},
                            {'label': 'Tungsten (SPARC)', 'value': 'TUNG_SPARC'},
                            {'label': 'User Defined', 'value': 'USER'},
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
        #we choose the same delta t for writing and the solver, to avoid
        #mismatched timesteps for solver and writer.  Advanced users
        #can uncomment the OFwriteDeltaT input div above and adapt this function
        #if desired
        OFwriteDeltaT = OFdeltaT
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
        return [dbc.Label("Loaded OF Data into HEAT", className="text-success"), OFdata]





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
                html.H4("HEAT Run Settings", style={"width":"100%"}),
                html.Br(),
                html.H5("Point Clouds at Tile Surface:"),
                runTabChecklist(),
                html.H5("Traces from Points:"),
                runTabTraces(),
                html.H5("Output File Settings:"),
                outputSettings(),
                dbc.Button("Run HEAT", id="runHEAT", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                html.Div(id="hiddenDivRun")
                ],
            className="wideBoxNoColor",
                )


#==========HEAT Run Settings==========
def runTabChecklist():
    return html.Div(
            id="runTabChecklist",
            children=[
                dbc.Card(
                    [dbc.CardBody(
                        dbc.Checklist(
                        options=[
                            {'label': 'B-field point cloud ', 'value': 'B'},
                            {'label': 'Normal vector point cloud', 'value': 'norm'},
                            {'label': 'powerDir point cloud', 'value': 'pwrDir'},
                            {'label': 'psiN point cloud', 'value': 'psiN'},
                            {'label': 'bdotn point cloud', 'value': 'bdotn'},
                            {'label': 'Heat flux point cloud', 'value': 'hfOpt'},
                            {'label': 'Ion gyro-orbit heat flux point cloud', 'value': 'hfGyro'},
                            {'label': 'Photon radiation heat flux point cloud', 'value': 'hfRad'},
                            {'label': 'openFOAM thermal analysis', 'value': 'T'},
                        ],
                        value=['hfOpt'],
                        id='checklistPC',
                        switch=True,
                        ),
                    )],
                    className="card100",
                )
            ],
                className="PCbox",
            )


def runTabTraces():
    return html.Div(
            id="runTabTraces",
            children=[
                dbc.Card(
                    [dbc.CardBody(
                        [
                        dbc.Checklist(
                            options=[{'label': 'B-field Trace ', 'value': 'Btrace'}],
                            value=[''],
                            id="Btrace",
                            switch=True,
                            ),
                        html.Div(id="bFieldTracePoints", children=[loadBfieldTrace(False)]),
                        ],
                        className="card100",
                    )
                    ],
                    className="card100",
                ),
                dbc.Card(
                    [dbc.CardBody(
                        [
                        dbc.Checklist(
                            options=[{'label': 'OpenFOAM Temp Probe ', 'value': 'OFtrace'}],
                            value=[''],
                            id="OFtrace",
                            switch=True,
                            ),
                        html.Div(id="OFTracePoints", children=[loadOFTrace(False)]),
                        ],
                        )
                    ],
                    className="card100",
                ),

                dbc.Card(
                    [dbc.CardBody(
                        [
                        dbc.Checklist(
                            options=[{'label': 'Gyro Orbit Trace ', 'value': 'gyrotrace'}],
                            value=[''],
                            id="gyrotrace",
                            switch=True,
                        ),
                        html.Div(id="gyroTracePoints", children=[loadGyroTrace(False)]),
                        ],
                        className="card100",
                    )
                    ],
                    className="card100",
                ),



            ],
            className="PCbox",
            )

def outputSettings():
    return html.Div(
            id="outputSettings",
            children=[
                dbc.Card(
                    [dbc.CardBody(
                        dbc.Checklist(
                        options=[
                            {'label': 'VTP Point Cloud / Glyphs ', 'value': 'vtpPC'},
                            {'label': 'VTP Mesh File ', 'value': 'vtpMesh'},
                            {'label': 'CSV File ', 'value': 'csv'},
                        ],
                        value=['csv', 'vtpMesh', 'vtpPC'],
                        id='checklistOutput',
                        switch=True,
                        ),
                    )],
                    className="card100",
                )
            ],
                className="PCbox",
            )

#=== magnetic field line traces
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
    params = ['x[mm]','y[mm]','z[mm]','traceDirection','Length[deg]','stepSize[deg]']
    cols = [{'id': p, 'name': p} for p in params]
    data = [{}]
    return dash_table.DataTable(
        id='BtraceTable',
        columns = cols,
        data = data,
#        style_header={'backgroundColor': 'rgb(30, 30, 30)'},
        style_cell={
            'textAlign': 'left',
            'color': 'black'
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
            params = ['x[mm]','y[mm]','z[mm]','traceDirection','Length[deg]','stepSize[deg]']
            tableColumns = [{'id': p, 'name': p} for p in params]
            hiddenDiv = [html.Label("Loaded file: "+filename, style={'color':'#34b3ed'})]
            dataStore.update({'BtraceFileName':filename})
            dataStore.update({'BtraceContents':uploadContents})
    else:
        hiddenDiv = []

    return dataStore, tableData, tableColumns, hiddenDiv


#=== openfoam temperature trace
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
                            dbc.Label("x [mm]"),
                            dbc.Input(id="xOFtrace", className="xyzBoxInput"),
                            dbc.Label("y [mm]"),
                            dbc.Input(id="yOFtrace", className="xyzBoxInput"),
                            dbc.Label("z [mm]"),
                            dbc.Input(id="zOFtrace", className="xyzBoxInput"),
                            #html.Label("t [ms]"),
                            #dcc.Input(id="tOFtrace", className="xyzBoxInput"),
                        ],
                        className="xyzBox"
                        )
                    ],
                style=style,
                className="xyzBoxVert",
                    )

#=== Gyro orbit traces
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
            dcc.Upload(
                className="PFCupload",
                id='gyroTrace-upload',
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
            html.Div(children=loadGyroTable(), className="PFCtable"),
            html.Div(id="hiddenDivGyroTrace"),
        ],
        className="xyzBoxVert",
        style=style,
    )

def loadGyroTable():
    params = ['x[mm]','y[mm]','z[mm]','T[eV]','gPhase[deg]', 'vPhase[deg]','N_helix','traceDirection','Length[deg]','stepSize[deg]']
    cols = [{'id': p, 'name': p} for p in params]
    data = [{}]
    return dash_table.DataTable(
        id='gyroTraceTable',
        columns = cols,
        data = data,
#        style_header={'backgroundColor': 'rgb(30, 30, 30)'},
        style_cell={
            'textAlign': 'left',
            'color': 'black',
        },
        editable=True,
        export_format='csv',
        row_deletable=True,
        persistence=True,
        )



#Uploading files and button presses update table and HTML storage object
@app.callback([Output('gyroTraceDataStorage', 'data'),
               Output('gyroTraceTable', 'data'),
               Output('gyroTraceTable', 'columns'),
               Output('hiddenDivGyroTrace', 'children')],
              [Input('gyroTrace-upload', 'filename')],
              [State('gyroTraceDataStorage', 'data'),
               State('gyroTrace-upload', 'contents'),
               State('gyroTraceTable', 'data'),
               State('gyroTraceTable', 'columns'),
               State('MachFlag', 'value'),
               State('gyrotrace', 'value')])
def gyroTraceTable(filename, dataStore, uploadContents,
             tableData, tableColumns, MachFlag, gyroTraceList):
    if 'gyrotrace' in gyroTraceList:
        trace = True
    else:
        trace = False

    if dataStore==None:
        print('Initializing HTML PFC data storage object')
        dataStore = {}
        dataStore.update({'gyroTraceFileName':None})
        dataStore.update({'gyroTraceContents':None})
        dataStore.update({'gyroTrace':trace})

    #user has to load MHD and CAD for a specific machine before PFCs
    if MachFlag == None:
        hiddenDiv = [html.Label("Select a Tokamak and Load MHD First", style={'color':'#db362a'})]
    elif filename is None:
        hiddenDiv = []
    #file dropper
    elif (trace == True) and (uploadContents!= dataStore['gyroTraceContents']):
            df = parse_contents(uploadContents, filename)
            tableData = df.to_dict('records')
            #tableColumns = [{"name": i, "id": i} for i in df.columns]
            params = ['x[mm]','y[mm]','z[mm]','T[eV]','gPhase[deg]', 'vPhase[deg]','N_helix','traceDirection','Length[deg]','stepSize[deg]']
            tableColumns = [{'id': p, 'name': p} for p in params]
            hiddenDiv = [html.Label("Loaded file: "+filename, style={'color':'#34b3ed'})]
            dataStore.update({'gyroTraceFileName':filename})
            dataStore.update({'gyroTraceContents':uploadContents})
    else:
        hiddenDiv = []

    return dataStore, tableData, tableColumns, hiddenDiv



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
               State('gyroTraceTable','data'),
               State('timeSlider','value'),
               State('inputsTable', 'data'),
               State('checklistOutput', 'value'),
               ])
def runHEAT(n_clicks,runList,Btrace,OFtrace,gyrotrace,
            BtraceTableData,
            xOFtrace,yOFtrace,zOFtrace,
            gyroTraceTableData, t, inputData,
            outputList):
    if n_clicks == 0:
        raise PreventUpdate

    #Bfield trace
    if 'Btrace' in Btrace:
        print("Tracing magnetic field lines...")
        log.info("Tracing magnetic field lines...")
        gui.BtraceMultiple(t, BtraceTableData)

    #gyro orbit trace
    if 'gyrotrace' in gyrotrace:
        gui.gyroTraceMultiple(gyroTraceTableData, t)

    #set up output file stream
    gui.getIOInputs(outputList)

    #run HEAT
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


    gui.writeInputTable(inputData)

    #set tree permissions
    tools.recursivePermissions(gui.MHD.shotPath, gui.UID, gui.GID, gui.chmod)
    print("\nReturned to GUI.  HEAT run complete.\n")
    log.info("\nReturned to GUI.  HEAT run complete.\n")

    return ([html.Label("HEAT Run Complete", className="text-success")],
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
            html.Br(),
            html.H6("Loaded gFile Parameters:", style={"text-align":"left", "width":"100%"}),
            html.Div( children=buildGfileTable(), className="gfileTable" ),
            dbc.Accordion(
                [
                dbc.AccordionItem(
                    [
                    dbc.Card(
                        dbc.CardBody(
                            [
                            html.Div( children=buildGfilePlots(), className="gfileTable" ),
                            html.Div( children=buildBFieldPlots(), className="gfileTable" ),
                            ]
                            )
                    )
                    ],
                    title="Plots",
                ),
                dbc.AccordionItem(
                    [
                        dbc.Card(
                            dbc.CardBody(
                            gFileMultipiers(),
                            )
                        )
                    ],
                    title="Multiply and Add"
                ),
                dbc.AccordionItem(
                    [
                        dbc.Card(
                            dbc.CardBody(
                                gFileNewSep(),
                            )
                        )
                    ],
                    title="Re-define LCFS"
                ),
                dbc.AccordionItem(
                    [
                        dbc.Card(
                            dbc.CardBody(
                                saveNewGfile(),
                            )
                        )
                    ],
                    title="Save New gFile"
                ),
                dbc.AccordionItem(
                    [
                        dbc.Card(
                            dbc.CardBody(
                                interpolateGfile(),
                            )
                        )
                    ],
                    title="Interpolate gFile"
                ),
                ],
                start_collapsed=True
            )
            ],
        className="wideBoxNoColor",
        )

def buildGfilePlots():
    """
    gFile plot for Fpol, psi, etc
    """
    return html.Div(
            id="gFilePlots",
            children=[dcc.Graph(id="gFilePlot1", className="card100"),],
            className="gfileBox"
            )

def buildBFieldPlots():
    """
    gFile plot for Fpol, psi, etc
    """
    return html.Div(
            id="bFieldPlots",
            children=[dcc.Graph(id="bFieldPlot", className="card100"),],
            className="gfileBox"
            )

def buildGfileTable():
    cols = getGfileColumns()
    data = [{}]
    return dash_table.DataTable(
        id='gFileTable',
        columns = cols,
        data = data,
        style_header={'backgroundColor': 'rgb(60, 60, 60)', 'color': 'white'},
        style_cell={
            'textAlign': 'left',
#            'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'black'
        },
        editable=False,
        row_deletable=False,
        )


#generic row data
def getGfileColumns():
    params = ['Parameter','Value']
    return [{'id': p, 'name': p} for p in params]

#load data
def geteqData(t=None):
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
            dbc.Label("psiRZ Multiplier", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="psiRZMult", className="gfileBoxInput", value="1.0"),
            dbc.Label("psiRZ Addition", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="psiRZAdd", className="gfileBoxInput", value="0.0"),
            dbc.Label("psiSep Multiplier", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="psiSepMult", className="gfileBoxInput", value="1.0"),
            dbc.Label("psiSep Addition", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="psiSepAdd", className="gfileBoxInput", value="0.0"),
            dbc.Label("psiAxis Multiplier", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="psiAxisMult", className="gfileBoxInput", value="1.0"),
            dbc.Label("psiAxis Addition", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="psiAxisAdd", className="gfileBoxInput", value="0.0"),
            dbc.Label("Fpol Multiplier", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="FpolMult", className="gfileBoxInput", value="1.0"),
            dbc.Label("Fpol Addition", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="FpolAdd", className="gfileBoxInput", value="0.0"),
            dbc.Label("Bt0 Multiplier", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="Bt0Mult", className="gfileBoxInput", value="1.0"),
            dbc.Label("Bt0 Addition", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="Bt0Add", className="gfileBoxInput", value="0.0"),
            dbc.Label("Ip Multiplier", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="IpMult", className="gfileBoxInput", value="1.0"),
            dbc.Label("Ip Addition", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="IpAdd", className="gfileBoxInput", value="0.0"),
            html.Br(),
            dbc.Checklist(options=[{'label':'Apply to all timesteps?', 'value':'all'}], id='correctAllts'),
            html.Br(),
            dbc.Button("Apply Corrections", id="applyMult", n_clicks=0, style={'margin':'0 10px 10px 0'}),
            html.Div(id="hiddenDivMult"),
            dbc.Label("*Sign of Bt0 and Ip checked for helicity in traces (MAFOT)", style={'margin':'0 10px 0 10px'}),
            dbc.Label("**Sign of psiRZ, psiSep, psiAxis, checked for helicity in point clouds (HEAT)", style={'margin':'0 10px 0 10px'}),
        ],
        className="gfileBox",
    )

def gFileNewSep():
    return html.Div(
        id="gfileNewSep",
        children=[
            dbc.Label("New LCFS R Value [m]", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="newLCFSr", className="gfileBoxInput", value="NA"),
            dbc.Label("New LCFS Z Value [m]", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="newLCFSz", className="gfileBoxInput", value="NA"),
            dbc.Label("(Value of NA in Z will choose minimum psi)"),
            html.Br(),
            dbc.RadioItems(
                            id="newLCFSradio",
                            options=[
                                {'label': 'Apply to All Timesteps', 'value': 'all'},
                                {'label': 'Apply to Currently Selected Timestep', 'value': 'single'},
                                    ],
                                    value='all'
                            ),
            dbc.Button("Re-Define LCFS", id="newLCFSbutton", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
            html.Div(id="hiddenDivSep"),
            dbc.Button("Find LCFS From PFCs", id="findLCFSbutton", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
            html.Div(id="hiddenDivSep2")
        ],
        className="gfileBox",
    )

def saveNewGfile():
    return html.Div(
        id="saveGfile",
        children=[
            dbc.Label("New gFile Name (suffix for multiple)", style={'margin':'0 10px 0 10px'}),
            dbc.Input(id="newGfileName", className="gfileBoxInput"),
            dbc.Button("Save New gFile", id="saveGfileButton", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),
            html.Br(),
            dbc.Checklist(options=[{'label':'Save all timesteps?', 'value':'all'}], id='saveAllGs'),
            html.Br(),
            dcc.Download(id="downloadNewGfile"),
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
                dbc.Label("Interpolation by Timestep", style={'margin':'0 10px 0 10px'}),
                dbc.Input(id="interpTime", className="gfileBoxInput"),
                dbc.Button("Interpolate this Timestep", id="interpButton", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                html.Div(id="hiddenDivInterp1"),
                html.Br(),
                dash_table.DataTable(
                    id='gFileTable2',
                    columns = ([{'id': p, 'name': p} for p in params]),
                    data = data,
                    style_header={'backgroundColor': 'rgb(60, 60, 60)', 'color': 'white'},
                    style_cell={
                        'textAlign': 'left',
#                        'backgroundColor': 'rgb(50, 50, 50)',
#                        'color': 'white'
                                },
                    editable=True,
                    row_deletable=False,

                ),
                html.Br(),

                dbc.Label("Interpolate N steps between gFiles", style={'margin':'0 10px 0 10px'}),
                html.Div(id='interpTable', className="gfileTable"), #updated from MHD button callback
                dbc.Input(id="interpN", className="gfileBoxInput", placeholder='Enter N steps'),
                dbc.Button("Interpolate these Timesteps", id="interpButton2", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                dcc.Download(id="download1InterpGfile"),
                dcc.Download(id="downloadInterpGfiles"),
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
    return [dbc.Label("Gfile Interpolated", className="text-success"),
            dcc.send_file(name)]


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
    df['timestep[ms]'] = df['timestep[ms]'].astype(int)
    df = df.sort_values('timestep[ms]', ascending=True)
    print("GEQDSK Dataframe:")
    print(df)
    print("New timesteps:")
    print(df['timestep[ms]'].values)
    #interpolate N steps between each point
    gui.interpolateNsteps(df['filename'].values, pd.to_numeric(df['timestep[ms]']).values,int(N))
    zipFile = gui.tmpDir + 'InterpolatedGfiles.zip'
    return [html.Label("gFiles Interpolated", className="text-success"),
            dcc.send_file(zipFile)]





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
               State('timeSlider', 'value'),
               State('correctAllts', 'value')])
def applyMult(n_clicks, psiRZMult, psiSepMult, psiAxisMult, FpolMult,
              psiRZAdd,psiSepAdd,psiAxisAdd,FpolAdd,
              Bt0Mult,Bt0Add,IpMult,IpAdd,t,correctAllts):
    """
    apply multiplier to psiRZ, psiSep, psiAxis, Fpol for currently
    selected equilibrium timestep
    """
    if n_clicks < 1:
        raise PreventUpdate
    #parse user formulas and convert then to number via python compiler
    pi = np.pi

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
        returnDiv = [html.Label("Corrections Applied", className="text-success")]

    except:
        print("Could not evaluate multipliers or adders")
        print("Please ensure values are of type float")
        print("and retry")
        returnDiv = [html.Label("Error with inputs", style={'color':'#f70c0c'})]


    gui.gfileClean(psiRZMult,psiSepMult,psiAxisMult,FpolMult,
                   psiRZAdd,psiSepAdd,psiAxisAdd,FpolAdd,
                   Bt0Mult,Bt0Add,IpMult,IpAdd,t,correctAllts)
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
    return [html.Label("Redefined LCFS", className="text-success")]

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
    return [html.Label("Iteratively Found LCFS", className="text-success")]



@app.callback([Output('hiddenDivSaveGfile', 'children'),
                Output('downloadNewGfile', 'data')],
              [Input('saveGfileButton', 'n_clicks')],
              [State('newGfileName','value'),
               State('timeSlider', 'value'),
               State('shot', 'value'),
               State('saveAllGs','value')])
def saveG(n_clicks, name, t, shot, saveAllGs):
    if n_clicks < 1:
        raise PreventUpdate

    if saveAllGs is None:
        allMask = False
        gui.writeGfile(name, shot, t)
        sendBack = gui.tmpDir + name
        outDiv = html.Label("Saved gFile", className="text-success")
    elif 'all' in saveAllGs:
        allMask = True
        sendBack = gui.createGfileZip()
        outDiv = html.Label("Saved gFiles", className="text-success")


    return [outDiv,
            dcc.send_file(sendBack)]

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
            dbc.Button("Download HEAT Results", id="downloadResults", n_clicks=0, style={'margin':'0 10px 10px 0', "width":"95%" }),
            html.Br(),
            html.Div(id="hiddenDivDownloadResults", style={"width":"100%"}),
            dcc.Download(id="downloadResultsDir"),
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
    """
    saves all HEAT results in the shotDir
    """
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
    return [html.Label("Saved HEAT output", className="text-success"),
            dcc.send_file(file+'.zip')]



def buildInputsTable(data=None):
    cols = getOutputColumns()
    if data==None:
        data = [{}]
    return dash_table.DataTable(
        id='inputsTable',
        columns = cols,
        data = data,
        style_header={'backgroundColor': 'rgb(60, 60, 60)', 'color':'white'},
        style_cell={
            'textAlign': 'left',
#            'backgroundColor': 'rgb(50, 50, 50)',
            'color': 'black'
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
               Input('RADDataStorage', 'data'),
               Input('OFDataStorage', 'data')])
def inputsTable(MHDdata,CADdata,HFdata,GYROdata,RADdata,OFdata):
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
    if RADdata != None:
        for key,val in RADdata.items():
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

        div = html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        div = html.Div(
            children=[
                html.Label("Run HEAT to get qDiv plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )

    return dbc.Card(dbc.CardBody(div))

def OFmaxTPlots(update=False):
    """
    div for maximum openFOAM temperature plots

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.getOFMinMaxPlots()

        div = html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        div = html.Div(
            children=[
                html.Label("Run OpenFOAM to get max(T) plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )

    return dbc.Card(dbc.CardBody(div))


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

        div =  html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        div = html.Div(
            children=[
                html.Label("Run OpenFOAM Temp Probe to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )
    return dbc.Card(dbc.CardBody(div))


def gyroPhasePlots(update=False):
    """
    div for gyro orbit phase angle plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.gyroPhasePlot()

        div = html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        div =  html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )

    return dbc.Card(dbc.CardBody(div))

def vPhasePlots(update=False):
    """
    div for gyro orbit velocity phase angle plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.vPhasePlot()

        div = html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        div = html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )
    return dbc.Card(dbc.CardBody(div))

def vSlicePlots(update=False):
    """
    div for gyro orbit velocity slice / distribution plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.vSlicePlot()

        div = html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        div = html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )
    return dbc.Card(dbc.CardBody(div))

def cdfSlicePlots(update=False):
    """
    div for gyro orbit CDF slice / distribution plot

    if update is False, just return an empty div, so that nothing happens on
    page load.  If update=True, get qDivs and update the plot.

    This is called at the end of a HEAT run (runHEAT callback) button click
    """
    if update==True:
        fig = gui.cdfSlicePlot()

        div = html.Div(
            className="plotBox",
            children=[
                dcc.Graph(id="", figure=fig),
                ],
                )
    else:

        div = html.Div(
            children=[
                html.Label("Run gyro orbit calculation to get plot", style={'color':'#52caeb'})
                ],
            className="gfileBox"
            )
    return dbc.Card(dbc.CardBody(div))




"""
==============================================================================
Tab Contents: logfile tab

DEPRECATED
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

"""
==============================================================================
Graphs, plots, etc.
"""
def build_graphs():
    return html.Div(
        id="graph-container",
        className="graphWindow",
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
            #the dcc.Slider values below are overwritten in the loadMHD function
            #these are for initialization only (min, max, value, marks)
            dcc.Slider(
                id='timeSlider',
                min=0, 
                max=100,
                step=None,
                value=None,
                marks={20: 'Load MHD to see timesteps'},
                tooltip={"placement": "bottom", "always_visible": True},
                updatemode='drag',
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
               Input('hiddenDivSep', 'children'),
               Input('themeStorage', 'data')])
def slideEQplot(value, dummy1, tom, themeData):
    try: MachFlag = gui.MachFlag
    except:
        print("You didn't select a machine")
        raise PreventUpdate
    try: ts = gui.MHD.timesteps
    except:
        print('Please load MHD')
        raise PreventUpdate
    try: idx = np.where(value==gui.MHD.timesteps)[0][0]
    except:
        print(value)
        print(gui.MHD.timesteps)
        print("Trouble loading MHD EQ")
        raise PreventUpdate
    ep = gui.MHD.ep[idx]
    shot = gui.MHD.shot
    t = gui.MHD.timesteps[idx]
    try:
        gName = gui.MHD.eqList[idx]
    except:
        gName = gui.MHD.eqList[0] #this is for when all EQ are in a single file (ie a JSON)
    #Update Equilibrium plot
    #this import needs to be here (not at top of file) to prevent FreeCAD qt5
    #shared object conflict
    import GUIscripts.plotly2DEQ as plotly2DEQ
    load_figure_template(template_from_url(themeData['theme']))
    plot = plotly2DEQ.makePlotlyEQDiv(shot, t, MachFlag, ep, gName=gName)
    data = geteqData(t)

    #update plot on gfile cleaner tab with Fpol, psi. etc
    plot2 = plotly2DEQ.makePlotlyGfilePlot(ep)

    #update plot on gfile cleaner tab with Bp, Bt
    plot3= plotly2DEQ.makePlotlyBpBt(ep, MachFlag)

    #this allows zoom to be preserved across button clicks
    plot.update_layout(uirevision='neverUpdate')
    plot2.update_layout(uirevision='neverUpdate')
    plot3.update_layout(uirevision='neverUpdate')

    return plot, data, plot2, plot3




"""
==============================================================================
Main Layout
"""
def generateLayout():
    app.layout = html.Div(
        id="big-app-container",
        className="entireWindow",
        children=[
            #store the session data until browser tab is closed
            dcc.Store(id='session', storage_type='memory'),
            dcc.Store(id='userInputFileData', storage_type='memory'),
            dcc.Store(id='PFCdataStorage', storage_type='memory'),
            dcc.Store(id='eqListStorage', storage_type='memory'),
            dcc.Store(id='BtraceDataStorage', storage_type='memory'),
            dcc.Store(id='gyroTraceDataStorage', storage_type='memory'),
            dcc.Store(id='MHDDataStorage', storage_type='memory'),
            dcc.Store(id='CADDataStorage', storage_type='memory'),
            dcc.Store(id='HFDataStorage', storage_type='memory'),
            dcc.Store(id='GYRODataStorage', storage_type='memory'),
            dcc.Store(id='RADDataStorage', storage_type='memory'),
            dcc.Store(id='OFDataStorage', storage_type='memory'),
            dcc.Store(id='themeStorage', storage_type='memory'),

            build_banner(),
            build_simulator(),
            html.Div(id="hiddenDiv", style={'display': 'none'}),
        ],
    )

    app.title = 'HEAT'

    return


def build_simulator():
    return html.Div(
        id="simulator",
        className="appWindow",
        children=[
            build_tabs(),
            #testDiv(), #use this to test out CSS classes
            build_graphs(),
        ],
    )

def testDiv():
    div = html.Div(
        className="tabWindow",
        children=[
            dbc.Card(
                [dbc.CardBody(dbc.Label("TEST"))],
                className="testCard",
                ),
            ]
    )
    return div

"""
==============================================================================
Session storage callbacks and functions
"""
# Load Default callback and input file drag and drop
@app.callback([Output('shot', 'value'),
               Output('traceLength', 'value'),
               Output('dpinit', 'value'),
               Output('psiMult1', 'value'),
               Output('BtMult1', 'value'),
               Output('IpMult1', 'value'),
               Output('gridRes', 'value'),
               Output('hfMode', 'value'), #this causes undefined vars
               Output('eichlqCNmode', 'value'), #this causes undefined vars
               Output('multiExplqCNmode', 'value'),
               Output('multiExplqCFmode', 'value'),
               Output('multiExplqPNmode', 'value'),
               Output('multiExplqPFmode', 'value'),
               Output('limiterlqCNmode', 'value'),
               Output('limiterlqCFmode', 'value'),
               Output('lqEich', 'value'),
               Output('S', 'value'),
               Output('eichSMode', 'value'),
               Output('lqCN', 'value'),
               Output('lqCF', 'value'),
               Output('lqPN', 'value'),
               Output('lqPF', 'value'),
               Output('lqTopHat', 'value'),
               Output('limlqCN', 'value'),
               Output('limlqCF', 'value'),
               Output('limfracCN', 'value'),
               Output('limfracCF', 'value'),
               Output('fracCN', 'value'),
               Output('fracCF', 'value'),
               Output('fracPN', 'value'),
               Output('fracPF', 'value'),
               Output('P', 'value'),
               Output('radFrac', 'value'),
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
               Output('gyroTraceLength','value'),
               Output('gyroT_eV','value'),
               Output('N_vSlice','value'),
               Output('N_vPhase','value'),
               Output('N_gyroPhase','value'),
               Output('ionMassAMU','value'),
               #Output('vMode','value'), #this causes undefined vars
               Output('ionFrac','value'),
               Output('phiMin','value'),
               Output('phiMax','value'),
               Output('Ntor','value'),
               Output('Nref','value'),
               Output('session', 'data'),
               Output('hiddenDivLoadDefaults', 'children')
               ],
               [Input('loadDefaults', 'n_clicks'),
                Input('userInputFileData', 'modified_timestamp')],
               [State('session', 'modified_timestamp'),
                State('MachFlag', 'value'),
                State('session', 'data'),
                State('userInputFileData', 'data')])
def session_data(n_clicks, inputTs, ts, MachFlag, data, inputFileData):
    #for some reason triggers are happening from weird things.  If nothing was
    #triggered, don't fire callback
    if ctx.triggered[0]['value']==None:
        raise PreventUpdate

    #default case
    if ts is None or MachFlag not in machineList:
        print('Initializing Data Store')
        if MachFlag not in machineList:
            print("Machine not in machine list")
            log.info("Machine not in machine list")
        data = gui.getDefaultDict()
        data.update({'default_n_clicks':0})

    #Let the user know if this worked or if we still need a MachFlag
    if MachFlag not in machineList:
        outputDiv = html.Label("Select Machine First", className="text-danger")
    else:
        outputDiv = html.Label("Loaded Input File", className="text-success")

    btnTest = n_clicks > 0 and n_clicks>data.get('default_n_clicks')
    #load defaults
    if btnTest and (MachFlag is not None):
        print("Loading Default Input File")
        log.info("Loading Default Input File")
        data = gui.loadInputs()
        data['default_n_clicks'] = n_clicks
        defaultDiv = [dbc.Label("Loaded Defaults", className="text-success")]
    elif inputTs == None or inputTs == -1:
        pass
    #use data we saved into storage object that we got from user input file
    else:
        print("Loading User Input File")
        log.info("Loading User Input File")
        inputFileData['default_n_clicks'] = n_clicks
        data = inputFileData
        defaultDiv = [dbc.Label("Loaded From File", className="text-success")]

    print("Updating Session Store")
    return [data.get('shot', ''),
            data.get('traceLength', ''),
            data.get('dpinit', ''),
            data.get('psiMult', ''),
            data.get('BtMult', ''),
            data.get('IpMult', ''),
            data.get('gridRes', ''),
            data.get('hfMode', ''), #these dropdowns cause undefined text boxes
            data.get('lqCNmode', ''),
            data.get('lqCNmode', ''),
            data.get('lqCFmode', ''),
            data.get('lqPNmode', ''),
            data.get('lqPFmode', ''),
            data.get('lqCNmode', ''),
            data.get('lqCFmode', ''),
            data.get('lqCN', ''), #eich
            data.get('S', ''),
            data.get('SMode', ''),
            data.get('lqCN', ''),
            data.get('lqCF', ''),
            data.get('lqPN', ''),
            data.get('lqPF', ''),
            data.get('lqCN', ''), #tophat
            data.get('lqCN', ''), #lim
            data.get('lqCF', ''), #lim
            data.get('fracCN', ''), #lim
            data.get('fracCF', ''), #lim
            data.get('fracCN', ''),
            data.get('fracCF', ''),
            data.get('fracPN', ''),
            data.get('fracPF', ''),
            data.get('P', ''),
            data.get('radFrac', ''),
            data.get('fracUI', ''),
            data.get('fracUO', ''),
            data.get('fracLI', ''),
            data.get('fracLO', ''),
            data.get('qBG', ''),
            data.get('fG', ''),
            data.get('qFilePath', ''),
            data.get('qFileTag', ''),
            # data.get('radFile',''),
            data.get('OFtMin', ''),
            data.get('OFtMax', ''),
            data.get('meshMinLevel', ''),
            data.get('meshMaxLevel', ''),
            data.get('STLscale', ''),
            data.get('deltaT', ''),
            data.get('writeDeltaT', ''),
            data.get('N_gyroSteps',''),
            data.get('gyroTraceLength',''),
            data.get('gyroT_eV',''),
            data.get('N_vSlice',''),
            data.get('N_vPhase',''),
            data.get('N_gyroPhase',''),
            data.get('ionMassAMU',''),
            #data.get('vMode',''),
            data.get('ionFrac',''),
            # data.get('radFile',''),
            data.get('phiMin',''),
            data.get('phiMax',''),
            data.get('Ntor',''),
            data.get('Nref',''),
            data,
            defaultDiv
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


    app.run(
            debug=True,
            dev_tools_ui=True,
            port=port,
            host=address
            )

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
different users' class variables from each other.  For now this is just set up to
serve a single session @ 127.0.0.1:8050

You will need to set a few variables below, based upon your system paths
rootDir, PVPath
"""

import base64
import io
import sys
import json
import numpy as np
import pandas as pd
import parser
import time
import copy
import dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_table
import EFIT.equilParams_class as EP

#========= VARIABLES THAT ARE SYSTEM DEPENDENT =================================
#Include the location of the paraview binaries.  Specifically we need 'pvpython'
PVPath = '/opt/paraview/ParaView-5.7.0-MPI-Linux-Python3.7-64bit/lib/python3.7/site-packages'
sys.path.append(PVPath)
#Root HEAT directory
rootDir = '/u/tlooby/source/HEAT/rev9/source/'
#openFOAM bashrc location
OFbashrc = '/opt/OpenFOAM/OpenFOAM-v1912/etc/bashrc'
#===============================================================================


#Create log files that will be displayed in the HTML GUI
import logging
logFile=rootDir + 'HEATlog.txt'
logFlask = logging.getLogger('werkzeug')
logFlask.disabled = True
logging.basicConfig(filename=logFile, filemode="w", level=logging.INFO, format='%(message)s')
log = logging.getLogger(__name__)

#Make sure all our python scripts are in the path
from GUIclass import GUIobj
app = dash.Dash(__name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}])

#Eventually need to fix this so that we are not using a global variable
#dash can acces Flask Cache so we should cache data by userID or something
#for R&D this works
gui = GUIobj(logFile, rootDir)

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
                                buildButtonRibbon(),
                                buildMHDbox(),
                                buildCADbox(),
                                buildPFCbox(),
                                buildHFbox(),
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
                        label="gFile Tools",
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
                        id="log-tab",
                        label="LogFile Output",
                        value="tab4",
                        style=tab_style,
                        selected_style=tab_selected_style,
                        children=[
                            html.Div(
#                            className='tabcontent',
                            className="innerTabContent",
                            children=[
                                html.Div([html.H3('logfile')], className="box")
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
            buildMachineSelector(),
            buildInputButtons(),
        ],
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
                style={'backgroundColor': 'transparent', 'color':'transparent'},
                #style=dropdown_style,
                options=[
                    {'label': 'NSTX-U', 'value': 'nstx'},
                    {'label': 'DIII-D', 'value': 'd3d'},
                    {'label': 'ST40', 'value': 'st40'}
        ],
        value=None
    ),
            ],
        )

def buildInputButtons():
    """
    returns Load Defaults and Upload Inputs buttons
    """
    return html.Div(
            id="buttonInputs",
            className="defaultButtonBox",
            children=[
                html.Button("Load Defaults", id="loadDefaults", n_clicks=0, className="defaultButtons"),
                html.Button("Upload Input File", id="uploadInputs", n_clicks=0, className="defaultButtons"),
            ],
        )

@app.callback(Output('hiddenDiv', 'children'),
              [Input('MachFlag', 'value')])
def machineSelector(MachFlag):
    gui.machineSelect(MachFlag)
    return




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
                html.Label(id="nTraceLabel", children="Number of Trace Steps (degrees)"),
                dcc.Input(id="nTrace", className="textInput"),
#                html.Label(id="ionDirLabel", children="Ion Direction"),
#                dcc.Input(id="ionDir", className="textInput"),
#                html.Label(id="gfileLabel", children="Path to gFile"),
#                dcc.Input(id="gfilePath", type='text', className="textInput"),
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

                dcc.RadioItems(
                                id="plasma3Dmask",
                                options=[
                                    {'label': '2D Plasmas', 'value': 'plasma2D'},
                                    {'label': '3D Plasmas', 'value': 'plasma3D'},
                                        ],
                                        value='plasma2D'
                                ),
                html.Br(),
                html.Button("Load MHD", id="loadMHD", n_clicks=0, style={'margin':'10px 10px 10px 10px'}),

            ],
            className="box",
        )


#==========MHD Callbacks
@app.callback([Output('timeSlider', 'min'),
               Output('timeSlider', 'max'),
               Output('timeSlider', 'marks'),
               Output('timeSlider', 'value'),
               Output('gFileTable2', 'data')],
              [Input('loadMHD', 'n_clicks')],
              [State('shot', 'value'),
              State('tmin', 'value'),
              State('tmax', 'value'),
              State('nTrace', 'value'),
              #State('ionDir', 'value'),
              #State('gfilePath', 'value'),
              State('gfiletable-upload', 'filename'),
              State('gfiletable-upload', 'contents'),
              State('plasma3Dmask', 'value')]
              )
def loadMHD(n_clicks,shot,tmin,tmax,nTrace,gFileList,gFileData,plasma3Dmask):
    """
    Load MHD
    """
    try: MachFlag = gui.MachFlag
    except:
        print("You didn't select a machine")
        raise PreventUpdate

    if plasma3Dmask == 'plasma3D': plasma3Dmask=1
    else:  plasma3Dmask=0

    if (shot is None) and (gFileList is None):
        raise PreventUpdate

    if shot is not None: shot = int(shot)
    if tmin is not None: tmin = int(tmin)
    if tmax is not None: tmax = int(tmax)
    if nTrace is not None: nTrace = int(nTrace)
#    if ionDir is not None: ionDir = int(ionDir)
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
                     plasma3Dmask=plasma3Dmask
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
    return tminMHD, tmaxMHD, marks, value, data



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
                html.Label(id="ROIgridResLabel", children="Heat Flux Resolution [mm]  "),
                dcc.Input(id="ROIGridRes", className="textInput"),
                html.Label(id="gridResLabel", children="Intersect Resolution [mm]"),
                dcc.Input(id="gridRes", className="textInput"),
                html.Label(id="STPfileLabel", children="STP File Path"),
                dcc.Input(id="STPfile", className="textInput"),
                html.Br(),
                html.Button("Load CAD", id="loadCAD", style={'margin':'10px 10px 10px 10px'}),
                html.Div(id="hiddenDivCAD")
            ],
            className="box",
        )

#Load CAD button connect
@app.callback([Output('hiddenDivCAD', 'children')],
              [Input('loadCAD', 'n_clicks')],
              [State('ROIGridRes', 'value'),
               State('gridRes', 'value'),
               State('STPfile', 'value'),
               State('MachFlag', 'value')])
def loadCAD(n_clicks, ROIGridRes, gridRes, STPfile, MachFlag):
    if MachFlag is None:
        raise PreventUpdate
    else:
        gui.getCADInputs(ROIGridRes,gridRes,STPfile)
    return [html.Label("Loaded CAD File", style={'color':'#f5d142'})]



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
                className="hfSelect",
                style={'backgroundColor': 'transparent', 'color':'transparent'},
                #style=dropdown_style,
                options=[
                    {'label': 'Gaussian Spreading', 'value': 'eich'},
                    {'label': 'Multi-Exponential', 'value': 'multiExp'},
                    {'label': 'Limiter', 'value': 'limiter'}
                    ],
                ),
                html.Div(id='hfParameters',
                         children=[
                                    loadHFSettings(mode=None,hidden=True),
                                   ],
                        ),
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
    return html.Div(
             className=className,
             children=[
                    html.Label("Power Crossing Separatrix [MW]"),
                    dcc.Input(id="Psol", className="textInput"),
                    ],
                    )

def loadHFSettings(mode=None, hidden=False):
    #do this hidden business so that we can always load defaults into these id's
    #hideMask corresponds to the following parameters:
    #[eichProfile, commonRegion, privateRegion]
    hideMask = ['hiddenBox','hiddenBox','hiddenBox']
    if mode=='eich':
        hideMask = ['hfInput','hiddenBox','hiddenBox']
    elif mode=='limiter':
        hideMask = ['hiddenBox','hiddenBox','hfInput'] #common flux region
    elif mode=='multiExp':
        hideMask = ['hiddenBox','hfInput','hiddenBox'] #common + private flux region
    if hidden==True or mode==None:
        hideMask=['hiddenBox','hiddenBox','hiddenBox']

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
                    {'label': 'From Eich Scaling', 'value': 'eich'},
                    {'label': 'User Defined', 'value': 'user'}
                    ],
                value=None,
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
                value=None,
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

    row3 = html.Div(
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

    row4 = html.Div(
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

    row5 = html.Div(
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
            #row1, #commented for now because only user defined is allowed (no regression)
            row2,
            row3,
            row4,
            row5
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
                        {'label': 'From Eich Scaling', 'value': 'eich'},
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
                    id='limiterlqCFMode',
                    className="SelectorBoxInput",
                    style={'backgroundColor': 'transparent', 'color':'transparent'},
                    options=[
                        {'label': 'From Horaceck Scaling', 'value': 'horaceck'},
                        {'label': 'User Defined', 'value': 'user'}
                        ],
                    value=None,
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
@app.callback([Output('hiddenDivHF', 'children')],
              [Input('loadHF', 'n_clicks')],
              [State('hfMode', 'value'),
               State('lqEich', 'value'),
               State('S', 'value'),
               State('qBG', 'value'),
               State('lqCN', 'value'),
               State('lqCF', 'value'),
               State('lqPN', 'value'),
               State('lqPF', 'value'),
               State('fracCN', 'value'),
               State('fracCF', 'value'),
               State('fracPN', 'value'),
               State('fracPF', 'value'),
               State('Psol', 'value'),
               State('LRmask', 'value'),
               State('LRthresh', 'value'),
               State('MachFlag', 'value'),
               State('eichlqCNMode', 'value'),
               #State('multiExplqCNMode', 'value'),
               State('limiterlqCNMode', 'value'),
               State('limiterlqCFMode', 'value'),
               State('eichSMode', 'value'),
               State('fG', 'value'),
               State('limlqCN', 'value'),
               State('limlqCF', 'value'),
               State('limfracCN', 'value'),
               State('limfracCF', 'value'),
               ])
def loadHF(n_clicks,hfMode,lqEich,S,qBG,lqCN,lqCF,lqPN,lqPF,
            fracCN,fracCF,fracPN,fracPF,Psol,LRmask,LRthresh,MachFlag,
            eichlqCNMode,limiterlqCNMode,limiterlqCFMode,SMode,fG,
            limlqCN,limlqCF,limfracCN,limfracCF):
    if MachFlag is None:
        raise PreventUpdate
    else:
        #set up the heat flux configuration (which scalings to use)
        if hfMode == 'limiter':
            lqCNmode = limiterlqCNMode
            lqCFmode = limiterlqCFMode
            lqCN = limlqCN
            lqCF = limlqCF
            fracCN = limfracCN
            fracCF = limfracCF
        elif hfMode == 'multiExp':
            lqCNmode = multiExplqCNMode
        else: #eich mode is default
            lqCNmode = eichlqCNMode
            lqCFmode = None
            SMode = SMode
        #could add private flux scalings here if they ever exist

        #set up HF object in HEAT
        gui.getHFInputs(lqEich,S,Psol,qBG,lqPN,lqPF,lqCN,lqCF,
                        fracPN,fracPF,fracCN,fracCF,hfMode,LRmask,LRthresh,
                        lqCNmode,lqCFmode,SMode,fG)

    return [html.Label("Loaded HF Settings", style={'color':'#f5d142'})]





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
                html.Div(id="hiddenDivPFC", style={"display": "none"}),
                dcc.Input(id="hiddenDivPFC2", style={"display": "none"}),
                html.Div(id="hiddenDivPFC3")

            ],
            className="PFCbox",
            )

def loadPFCtable():
    params = ['timesteps','PFCname','MapDirection','intersectName']
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
               Output('hiddenDivPFC3', 'children')],
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
        #load PFC data from page and send it to HEAT
        else:
            gui.getPFCdataFromGUI(tableData)
            gui.getPFCinputs(defaultMask=False)

        hiddenDiv = [html.Label("Loaded PFC Data into HEAT", style={'color':'#f5d142'})]

    #file dropper
    elif filename != dataStore['PFCfilename']:
        df = parse_contents(uploadContents, filename)
        tableData = df.to_dict('records')
        tableColumns = [{"name": i, "id": i} for i in df.columns]
        hiddenDiv = [html.Label("Loaded file: "+filename, style={'color':'#34b3ed'})]
    else:
        hiddenDiv = []


    return dataStore, tableData, tableColumns, hiddenDiv


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
                    html.Label("Minimum Resh Refinement Level"),
                    dcc.Input(id="OFminMeshLev", className="textInput"),
                ],
                className="OFInput",
            ),
            html.Div(
                children=[
                    html.Label("Maximum Resh Refinement Level"),
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
                    html.Label("deltaT [ms]"),
                    dcc.Input(id="OFdeltaT", className="textInput"),
                ],
                className="OFInput",
            ),
            ],
            className="wideBoxNoColor",
            )

#Load OF button connect
@app.callback([Output('hiddenDivOF', 'children')],
              [Input('loadOF', 'n_clicks')],
              [State('OFstartTime', 'value'),
               State('OFstopTime', 'value'),
               State('OFminMeshLev', 'value'),
               State('OFmaxMeshLev', 'value'),
               State('OFSTLscale', 'value'),
               State('OFdeltaT', 'value')
              ])
def loadOF(n_clicks,OFstartTime,OFstopTime,
            OFminMeshLev,OFmaxMeshLev,OFSTLscale,OFdeltaT):
    """
    sets up openFOAM for an analysis
    """
    if n_clicks == 0:
        raise PreventUpdate
    gui.loadOF(OFstartTime,OFstopTime,OFminMeshLev,OFmaxMeshLev,OFSTLscale,OFbashrc,OFdeltaT)
    return [html.Label("Loaded OF Data into HEAT", style={'color':'#f5d142'})]


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
                html.H6("Point Clouds at Tile Surface:"),
                runTabChecklist(),
                html.H6("Traces from Points:"),
                runTabTraces(),
                html.Button("Run HEAT", id="runHEAT", n_clicks=0, style={'margin':'0 10px 10px 0'}),
                html.Div(id="hiddenDivRun")
                ],
            className="wideBoxDark",
                )


def runTabChecklist():
    return html.Div(
            id="runTabChecklist",
            children=[
                dcc.Checklist(
                    options=[
                        {'label': 'B-field point cloud ', 'value': 'Bpc'},
                        {'label': 'Normal vector point cloud', 'value': 'NormPC'},
                        {'label': 'ShadowMask point cloud', 'value': 'shadowPC'},
                        {'label': 'psiN point cloud', 'value': 'psiPC'},
                        {'label': 'bdotn point cloud', 'value': 'bdotnPC'},
                        {'label': 'Heat flux point cloud', 'value': 'HFpc'},
                        {'label': 'openFOAM thermal analysis', 'value': 'OFpc'}
                        ],
                        value=['HFpc'],
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
                html.Div(id="bFieldTracePoints", children=[loadBfieldTrace(hidden=True)]),
                dcc.Checklist(
                    options=[{'label': 'OpenFOAM Temp Probe ', 'value': 'OFtrace'}],
                    value=[''],
                    id="OFtrace",
                        ),
                html.Div(id="OFTracePoints", children=[loadOFTrace(hidden=True)]),

                    ],
                className="PCbox",
                )

@app.callback(Output('bFieldTracePoints', 'children'),
              [Input('Btrace', 'value')],)
def bfieldTracePoint(value):
    return [loadBfieldTrace(Btrace=value)]

#this function enables the inputs to be rendered on page load but hidden
def loadBfieldTrace(Btrace=None, hidden=False):
    if (Btrace is 'Btrace') or (hidden is True):
        style={}
    else:
        style={"display":"hidden"}
    return html.Div(
                children=[
                    html.Label("x [mm]"),
                    dcc.Input(id="xBtrace", className="xyzBoxInput"),
                    html.Label("y [mm]"),
                    dcc.Input(id="yBtrace", className="xyzBoxInput"),
                    html.Label("z [mm]"),
                    dcc.Input(id="zBtrace", className="xyzBoxInput"),
                    html.Label("ionDirection"),
                    dcc.Input(id="ionDir", className="xyzBoxInput"),
                    ],
                style=style,
                className="xyzBox",
                    )

@app.callback(Output('OFTracePoints', 'children'),
              [Input('OFtrace', 'value')])
def OFTracePoint(value):
    return [loadOFTrace(OFtrace=value)]

#this function enables the inputs to be rendered on page load but hidden
def loadOFTrace(OFtrace=None, hidden=False):
    if (OFtrace is 'OFtrace') or (hidden is True):
        style={}
    else:
        style={"display":"hidden"}
    return html.Div(
                children=[
                    html.Label("x [mm]"),
                    dcc.Input(id="xOFtrace", className="xyzBoxInput"),
                    html.Label("y [mm]"),
                    dcc.Input(id="yOFtrace", className="xyzBoxInput"),
                    html.Label("z [mm]"),
                    dcc.Input(id="zOFtrace", className="xyzBoxInput"),
                    #html.Label("t [ms]"),
                    #dcc.Input(id="tOFtrace", className="xyzBoxInput"),
                    ],
                style=style,
                className="xyzBox",
                    )


@app.callback(Output('hiddenDivRun', 'children'),
              [Input('runHEAT', 'n_clicks')],
              [State('checklistPC','value'),
               State('Btrace','value'),
               State('OFtrace','value'),
               State('xBtrace','value'),
               State('yBtrace','value'),
               State('zBtrace','value'),
               State('xOFtrace','value'),
               State('yOFtrace','value'),
               State('zOFtrace','value'),
               State('timeSlider', 'value'),
               State('ionDir', 'value'),
               ])
def runHEAT(n_clicks,runList,Btrace,OFtrace,
            xBtrace,yBtrace,zBtrace,
            xOFtrace,yOFtrace,zOFtrace,t,ionDir):
    if n_clicks == 0:
        raise PreventUpdate

    if 'Btrace' in Btrace:
        gui.Btrace(xBtrace,yBtrace,zBtrace,t,ionDir)

    gui.runHEAT(runList)

    if 'OFpc' in runList:
        gui.runOpenFOAM()

#NEED TO ADD THIS BACK IN WHEN OF MODULE BUILT INTO DASH
#    if thermal_flag in checklist:
#        parts, intersects = gui.loadPFCDefaults()
#        part = parts[0]+'_'+gui.CAD.ROIGridRes+'mm'
#        gui.runOpenFOAM(part)
#    if Tprobe_flag == 'true':
#        gui.TprobeOF(xProbe,yProbe,zProbe)
    return [html.Label("HEAT Run Complete", style={'color':'#f5d142'})]



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
            html.H4("gFile Tools", style={"text-align":"center", "width":"100%"}),
            html.H6("Loaded gFile Parameters:"),
            html.Div( children=buildGfileTable(), className="gfileTable" ),
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
            html.Label("psiSep Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiSepMult", className="gfileBoxInput", value="1.0"),
            html.Label("psiAxis Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="psiAxisMult", className="gfileBoxInput", value="1.0"),
            html.Label("Fpol Multiplier", style={'margin':'0 10px 0 10px'}),
            dcc.Input(id="FpolMult", className="gfileBoxInput", value="1.0"),
            html.Button("Apply Multipliers", id="applyMult", n_clicks=0, style={'margin':'0 10px 10px 0'}),
            html.Div(id="hiddenDivMult")
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
                html.Div(id="hiddenDivInterp2")
            ],
            className="gfileBox",
            )


@app.callback(Output('hiddenDivInterp1', 'children'),
              [Input('interpButton', 'n_clicks')],
              [State('interpTime','value'),
              ])
def interpolate(n_clicks, t):
    """
    interpolate gfile at user defined timestep
    """
    if n_clicks < 1:
        raise PreventUpdate
    gui.interpolateGfile(t)
    return [html.Label("Gfile Interpolated", style={'color':'#f5d142'})]


@app.callback(Output('hiddenDivInterp2', 'children'),
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
    #interpolate N steps between each point
    gui.interpolateNsteps(df['filename'].values, pd.to_numeric(df['timestep[ms]']).values,int(N))
    return [html.Label("Gfile Interpolated", style={'color':'#f5d142'})]





@app.callback(Output('hiddenDivMult', 'children'),
              [Input('applyMult', 'n_clicks')],
              [State('psiRZMult','value'),
               State('psiSepMult','value'),
               State('psiAxisMult','value'),
               State('FpolMult','value'),
               State('timeSlider', 'value')])
def applyMult(n_clicks, psiRZMult, psiSepMult, psiAxisMult, FpolMult,t):
    """
    apply multiplier to psiRZ, psiSep, psiAxis, Fpol for currently
    selected equilibrium timestep
    """
    if n_clicks < 1:
        raise PreventUpdate
    #parse user formulas and convert then to number via python compiler
    pi = np.pi
    psiRZMult = eval(parser.expr(psiRZMult).compile())
    psiSepMult = eval(parser.expr(psiSepMult).compile())
    psiAxisMult = eval(parser.expr(psiAxisMult).compile())
    FpolMult = eval(parser.expr(FpolMult).compile())

    gui.gfileClean(psiRZMult,psiSepMult,psiAxisMult,FpolMult, t)
    return [html.Label("Multipliers Applied", style={'color':'#f5d142'})]

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



@app.callback(Output('hiddenDivSaveGfile', 'children'),
              [Input('saveGfileButton', 'n_clicks')],
              [State('newGfileName','value'),
               State('timeSlider', 'value'),
               State('shot', 'value')])
def saveG(n_clicks, filename, t, shot):
    if n_clicks < 1:
        raise PreventUpdate
    gui.writeGfile(filename, shot, t)
    return [html.Label("Saved gFile", style={'color':'#f5d142'})]

"""
==============================================================================
Tab Contents: logfile tab
"""







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

@app.callback([Output('2DEQ', 'figure'),
               Output('gFileTable', 'data')],
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
    import GUIscripts.plotly2DEQ as plotly2DEQ
    plot = plotly2DEQ.makePlotlyEQDiv(shot, t, MachFlag, ep)
    data = getGfileData(t)
    return plot, data





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


app.layout = html.Div(
        id="big-app-container",
        children=[
            #store the session data until browser tab is closed
            dcc.Store(id='session', storage_type='session'),
            dcc.Store(id='PFCdataStorage', storage_type='session'),
            dcc.Store(id='gFileListStorage', storage_type='session'),
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
# output the stored data in the GUI
@app.callback([Output('shot', 'value'),
               Output('tmin', 'value'),
               Output('tmax', 'value'),
               Output('nTrace', 'value'),
               Output('ionDir', 'value'),
               Output('ROIGridRes', 'value'),
               Output('gridRes', 'value'),
               Output('STPfile', 'value'),
               Output('lqEich', 'value'),
               Output('S', 'value'),
               Output('lqCN', 'value'),
               Output('lqCF', 'value'),
               Output('lqPN', 'value'),
               Output('lqPF', 'value'),
               Output('fracCN', 'value'),
               Output('fracCF', 'value'),
               Output('fracPN', 'value'),
               Output('fracPF', 'value'),
               Output('Psol', 'value'),
               Output('qBG', 'value'),
               Output('OFstartTime', 'value'),
               Output('OFstopTime', 'value'),
               Output('OFminMeshLev', 'value'),
               Output('OFmaxMeshLev', 'value'),
               Output('OFSTLscale', 'value'),
               Output('OFdeltaT', 'value'),
               Output('session', 'data'),
               ],
               [Input('loadDefaults', 'n_clicks')],
               [State('session', 'modified_timestamp'),
                State('MachFlag', 'value'),
                State('session', 'data')])
def session_data(n_clicks, ts, MachFlag, data):
    if ts is None or MachFlag not in ['nstx', 'd3d', 'st40']:
        print('Initializing Data Store')
        data = gui.getDefaultDict()
        data.update({'default_n_clicks':n_clicks})

    #load defaults
    if n_clicks > 0 and n_clicks>data.get('default_n_clicks') and MachFlag is not None:
        data = gui.loadDefaults()
        data['default_n_clicks'] = n_clicks
    else:
        pass
    print("Updating Session Store")
    return [data.get('shot', ''),
            data.get('tmin', ''),
            data.get('tmax', ''),
            data.get('nTrace', ''),
            data.get('ionDirection', ''),
            data.get('ROIGridRes', ''),
            data.get('gridRes', ''),
            data.get('STPfile', ''),
            data.get('lqEich', ''),
            data.get('S', ''),
            data.get('lqCN', ''),
            data.get('lqCF', ''),
            data.get('lqPN', ''),
            data.get('lqPF', ''),
            data.get('fracCN', ''),
            data.get('fracCF', ''),
            data.get('fracPN', ''),
            data.get('fracPF', ''),
            data.get('Psol', ''),
            data.get('qBG', ''),
            data.get('tMin', ''),
            data.get('tMax', ''),
            data.get('meshMinLevel', ''),
            data.get('meshMaxLevel', ''),
            data.get('STLscale', ''),
            data.get('deltaT', ''),
            data
            ]


if __name__ == '__main__':
    app.run_server(debug=True, dev_tools_ui=True)

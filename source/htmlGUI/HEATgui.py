#HEATgui.py
#Description:   GUI Module
#Engineer:      T Looby
#Date:          20200115
"""
This is the FLASK interface between the HTML/JS/CSS scripts and the python
scripts.  FLASK uses jinja to build the html on the fly.
We call the functions below from the client side using a javascript ajax request,
where the url is the function name.  To make such a call, we must wrap the function
in the @app.route decorator.
In the if __name__ == '__main__' section we define the flask server's ip address.
The HTML/JS/CSS scripts that the flask is binding to is located in this directory.
The directory structure matters, so don't mess with it unless you know what you
are doing!
"""


from flask import Flask, render_template, url_for, request, Response, redirect
from flask import jsonify, make_response, send_file
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import random
import io
import base64
import sys
import parser
import numpy as np
import os
import shutil

#Make sure all our python scripts are in the path
PYTHONPATH = '../'
sys.path.append(PYTHONPATH)
from GUIclass import create_app, GUIobj

#Include the location of the paraview binaries.  Specifically we need 'pvpython'
PVPath = '/opt/paraview/ParaView-5.7.0-MPI-Linux-Python3.7-64bit/lib/python3.7/site-packages'
sys.path.append(PVPath)

#Root HEAT directory
rootDir = '/u/tlooby/source/HEAT/rev6/source/'

#Create log files that will be displayed in the HTML GUI
import logging
logFile=rootDir + 'HEATlog.txt'
logFlask = logging.getLogger('werkzeug')
logFlask.disabled = True
logging.basicConfig(filename=logFile, filemode="w", level=logging.INFO, format='%(message)s')
log = logging.getLogger(__name__)



#2D Equilibrium plot outfile (html file for interactive plots)
defaultEQfile = rootDir + 'htmlGUI/static/html/defaultEQ.html'
#html2DEQfile = rootDir + 'htmlGUI/static/html/EQ2D.html'
html2DEQfile = rootDir + 'htmlGUI/templates/EQ2D.html'

#Initialize GUI and FLASK objects
gui = GUIobj(logFile, rootDir)
app = create_app(gui)

@app.route('/')
def index():
    """
    Load index.html page
    """
    return render_template('index.html')

@app.route('/run')
def run():
    """
    Load Run Page
    """
    return render_template('run.html')

@app.route('/wiki')
def wiki():
    """
    Load wiki page
    """
    return render_template('wiki.html')

@app.route('/equil')
def equil():
    """
    Load 2D MHD equilibrium plot (made via plotly)
    """
    print('Re-render equilibrium plot')
    return render_template('EQ2D.html')

@app.after_request
def add_header(response):
    """
    No caching at all for API endpoints.
    """
    # response.cache_control.no_store = True
    if 'Cache-Control' not in response.headers:
        response.headers['Cache-Control'] = 'no-store'
    return response

@app.route('/plot.png')
def plot_png():
    """
    Generate a plot example random figure
    """
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure():
    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)
    xs = range(100)
    ys = [random.randint(1, 50) for x in xs]
    axis.plot(xs, ys)
    return fig


@app.route('/loadMHD', methods=['POST'])
def loadMHD():
    """
    Load MHD from gfile or MDS+ based upon MHD tab in GUI
    """
    #Erase contents of log file:
    open(logFile, 'w').close()

    print("Loading MHD")
    log.info("Loading MHD")
    #Update MHD variables
    #vars = request.get_json()
    if request.form.get('shot'): shot = int(request.form.get('shot'))
    else: shot = None
    if request.form.get('tmin'): tmin = int(request.form.get('tmin'))
    else: tmin = None
    if request.form.get('tmax'): tmax = int(request.form.get('tmax'))
    else: tmax = None
    if request.form.get('nTrace'): nTrace = int(request.form.get('nTrace'))
    else: nTrace = None
    if request.form.get('gfilePath'): gfilePath = request.form.get('gfilePath')
    else:gfilePath = None

    gui.getMHDInputs(shot,tmin,tmax,nTrace,gfilePath)

    #Update Equilibrium plot on client
    #Now using plotly interactive plots
    import GUIscripts.plotly2DEQ as plotly2DEQ
    ep = app.config['GUI'].MHD.ep
    shot = app.config['GUI'].MHD.shot
    t = app.config['GUI'].t
    if os.path.exists(html2DEQfile):
        os.remove(html2DEQfile)
    plotly2DEQ.writePlotlyEQ(shot, t, html2DEQfile, ep=ep, gfile=None)

#    height = float(request.form.get('height')) #available height in html page

#    import GUIscripts.plot2DEQ as EQplot
#    ep = app.config['GUI'].MHD.ep
#    shot = app.config['GUI'].MHD.shot
##    t = app.config['GUI'].MHD.timesteps[0]
#    t = app.config['GUI'].t
#    plt = EQplot.EQ2Dplot(ep,shot,t,gui.MachFlag,height)
#    #Create plt figure in memory buffer
#    output = io.BytesIO()
#    plt.savefig(output, format='png', transparent=True)
#    #Now Create JSON response for ajax running on client side.
#    #We change output.getvalue() to base64 encoded variable so that ajax
#    #script can read it directly on client after it is returned
#    response = make_response(base64.b64encode(output.getvalue()))
#    response.headers.set('Content-Type', 'image/png')
    print("MHD Load Complete")
    log.info("MHD Load Complete")
#    response = {'outFile': html2DEQfile}
#    return jsonify(response)
#    return render_template('EQ2D.html')
    return ('', 204)

@app.route('/gfileClean', methods=['POST'])
def gfileClean():
    """
    Apply multiplication to a few variables in the MHD object
    """
    print("Cleaning Gfile")
    log.info("Cleaning Gfile")
    #parse user formulas and convert then to number via python compiler
    pi = np.pi

    psiRZMult = eval(parser.expr(request.form.get('psiRZ')).compile())
    psiSepMult = eval(parser.expr(request.form.get('psiSep')).compile())
    psiAxisMult = eval(parser.expr(request.form.get('psiAxis')).compile())
    FpolMult = eval(parser.expr(request.form.get('Fpol')).compile())

    gui.gfileClean(psiRZMult,psiSepMult,psiAxisMult,FpolMult)

    height = float(request.form.get('height')) #available height in html page
    #Update Equilibrium plot on client
    import GUIscripts.plot2DEQ as EQplot
    ep = app.config['GUI'].MHD.ep
    print(ep.g['psiRZ'])
    shot = app.config['GUI'].MHD.shot
#    t = app.config['GUI'].MHD.timesteps[0]
    t = app.config['GUI'].t
    plt = EQplot.EQ2Dplot(ep,shot,t,gui.MachFlag,height)
    #Create plt figure in memory buffer
    output = io.BytesIO()
    plt.savefig(output, format='png', transparent=True)
    #Now Create JSON response for ajax running on client side.
    #We change output.getvalue() to base64 encoded variable so that ajax
    #script can read it directly on client after it is returned
    response = make_response(base64.b64encode(output.getvalue()))
    response.headers.set('Content-Type', 'image/png')
    log.info("MHD Load Complete")
    return response


@app.route('/writeGfile', methods=['POST'])
def writeGfile():
    """
    write gfile.  used in gfile cleaning page after user applies multipliers
    with the gfileClean function
    """
    print("Writing gFile")
    log.info("Writing gFile")
    newShot = request.form.get('newShot')
    newTime = request.form.get('newTime')
    newGfile = request.form.get('newGfile')

    gui.writeGfile(newShot, newTime, newGfile)

    print("gFile write complete")
    log.info("gFile write complete")
    return ('', 204)

@app.route('/loadCAD', methods=['POST'])
def loadCAD():
    """
    Load CAD STP file and generate meshes based upon CAD tab in GUI
    """
    #Erase contents of log file:
    open(logFile, 'w').close()

    print("Loading CAD")
    log.info("Loading CAD")
#    vars = request.get_json()
    ROIGridRes = request.form.get('ROIGridRes')
    gridRes = request.form.get('gridRes')
    STPfile = request.form.get('STPfile')
    gui.getCADInputs(ROIGridRes, gridRes, STPfile)
    log.info("CAD Load Complete")
    return ('', 204)

@app.route('/loadHF', methods=['POST'])
def loadHF():
    """
    Load HF based upon HF tab in GUI
    """
    #Erase contents of log file:
    open(logFile, 'w').close()

    print("Loading HF Variables")
    log.info("Loading HF Variables")
    lq = float(request.form.get('lq'))
    S = float(request.form.get('s'))
    P = float(request.form.get('p'))
    qBG = float(request.form.get('qBG'))
    gui.getHFInputs(lq,S,P,qBG)
    log.info("HF Load Complete")
    return ('', 204)

@app.route('/runHEAT', methods=['POST'])
def runHEAT():
    """
    Run HEAT based upon desired outputs from RUN tab in GUI
    """
    #Erase contents of log file:
    open(logFile, 'w').close()

    print('Heat Run initiated...')
    log.info('Heat Run initiated...')
    Bpc_flag = request.form.get('Bpc')
    Btrace_flag = request.form.get('Btrace')
    Normpc_flag = request.form.get('Normpc')
    HFpc_flag = request.form.get('HFpc')
    psiPC_flag = request.form.get('psiPC')
    shadow_flag = request.form.get('shadow')
    xCoord = float(request.form.get('XCoord'))
    yCoord = float(request.form.get('YCoord'))
    zCoord = float(request.form.get('ZCoord'))

    thermal_flag = False
    Tprobe_flag = False
    thermal_flag = request.form.get('thermal')
    Tprobe_flag = request.form.get('Tprobe')
    xProbe = float(request.form.get('xProbe'))*1000.0
    yProbe = float(request.form.get('yProbe'))*1000.0
    zProbe = float(request.form.get('zProbe'))*1000.0


    if Bpc_flag == 'true':
        gui.bfieldAtSurface()
    if Btrace_flag == 'true':
        gui.Btrace(xCoord,yCoord,zCoord)
    if Normpc_flag == 'true':
        gui.NormPC()
    if HFpc_flag == 'true':
        gui.HFPC()
    if shadow_flag == 'true':
        gui.shadowPC(HFpc_flag)
    if psiPC_flag == 'true':
        gui.psiPC()
    if thermal_flag == 'true':
        parts, intersects = gui.loadPFCDefaults()
        part = parts[0]+'_'+gui.CAD.ROIGridRes+'mm'
        gui.runOpenFOAM(part)
    if Tprobe_flag == 'true':
        gui.TprobeOF(xProbe,yProbe,zProbe)



    log.info("HEAT Run Complete")
    return ('', 204)


@app.route('/setupOpenFOAM', methods=['POST'])
def setupOpenFOAM():
    """
    Run OpenFOAM based upon desired outputs from FEM tab in GUI
    """
    print('Preparing OpenFOAM Simulation')
    log.info('Preparing OpenFOAM Simulation')

    xMin = request.form.get('xMinOF')
    xMax = request.form.get('xMaxOF')
    yMin = request.form.get('yMinOF')
    yMax = request.form.get('yMaxOF')
    zMin = request.form.get('zMinOF')
    zMax = request.form.get('zMaxOF')
    tMin = request.form.get('tMinOF')
    tMax = request.form.get('tMaxOF')
    deltaT = request.form.get('deltaT')
    writeDeltaT = request.form.get('writeDeltaT')
    STLscale = request.form.get('STLscale')
    meshMinLevel = request.form.get('meshMinLevel')
    meshMaxLevel = request.form.get('meshMaxLevel')
    xMid = request.form.get('xMidOF')
    yMid = request.form.get('yMidOF')
    zMid = request.form.get('zMidOF')

    dict = {'xMin':xMin,'xMax':xMax,'yMin':yMin,'yMax':yMax,
            'zMin':zMin,'zMax':zMax,'tMin':tMin,'tMax':tMax,
            'deltaT':deltaT,'writeDeltaT':writeDeltaT,'STLscale':STLscale,
            'meshMinLevel':meshMinLevel,'meshMaxLevel':meshMaxLevel,
            'xMid':xMid,'yMid':yMid,'zMid':zMid}

    #initialize OF
    gui.getOpenFOAMinputs(dict)
    parts, intersects = gui.loadPFCDefaults()
    part = parts[0]+'_'+gui.CAD.ROIGridRes+'mm'
    gui.setUpOpenFOAM(part)

    log.info("OpenFOAM Ready for Execution")
    return ('', 204)


@app.route('/loadDefaults', methods=['POST'])
def loadDefaults():
    """
    load the infile and return the values that need to be displayed in the GUI
    when user selects Load Defaults button
    """
    #Erase contents of log file:
    open(logFile, 'w').close()

    print("Loading Default Variables")
    log.info("Loading Default Variables")
    defaults = gui.loadDefaults()
    log.info("Defaults Loaded")
    return jsonify(defaults)


@app.route('/refreshLog', methods=['POST'])
def refreshLog():
    """
    parse the HEATlog.txt file and push its contacts back to the GUI for
    display in the integrated terminal window thing
    """
    try:
        with open(logFile, "r") as f:
            content = f.read()
        return Response(content, mimetype='text/plain')
    except:
        return ('', 204)


@app.route('/uploadFile', methods=['POST'])
def uploadFile():
    """
    Upload a client provided infile to server
    """
    log.info('Uploading Input File...')
    data = request.form.get('result')
    #name = request.form.get('fileName')
    name = 'inputs/uploaded_input.csv'
    with open(rootDir + name,'w') as file:
        file.write(data)
    log.info('Upload Complete.')
    defaults = gui.loadDefaults(rootDir+name)
    log.info('Loaded parameters from file')
    return jsonify(defaults)


#Load potential PFC part numbers from input file
@app.route('/loadPFCParts', methods=['POST'])
def loadPFCParts():
    """
    Import a list of PFC parts from the PartsFile that we use to populate
    GUI checkboxes when the html page is initially loaded
    """
    log.info('Importing PFC Parts...')
    parts = gui.loadPFCParts()
    log.info('Successfully Loaded Parts...')
    return jsonify(parts)

@app.route('/loadPFCDefaults', methods=['POST'])
def loadPFCDefaults():
    """
    returns a dictionary with two arrays
    1) PFC HF calculation parts
    2) Intersection calculation parts
    these are loaded from default Parts and Intersect files, respectively
    """
    log.info('Updating PFC Parts...')
    parts, intersects = gui.loadPFCDefaults()
    keys = ['parts','intersects']
    lists = [parts,intersects]
    dictionary = dict(zip(keys,lists))
    log.info('Successfully Updated Parts...')
    return jsonify(dictionary)

@app.route('/writePFCs', methods=['POST'])
def writePFCs():
    """
    write a new PartsFile and IntersectFile, based upon user selections in the GUI
    """
    log.info('Writing PFC Parts...')
    parts = request.form.get('parts')
    intersects = request.form.get('intersects')
    gui.writePFCs(parts, intersects)
    log.info('Successfully Wrote Parts...')
    return ('', 204)

@app.route('/machineSelect', methods=['POST'])
def machineSelect():
    """
    Launches HEAT with the default settings associated with a particular machine
    """
    machine = request.form.get('machine')
    gui.machineSelect(machine)
    print('Launching HEAT for '+gui.MachFlag)
    run()

    return ('', 204)


if __name__ == '__main__':
    #Setup default plots before we load MHD
    shutil.copy(defaultEQfile, html2DEQfile)
#Uncomment this to run as a server on network
#    app.run(host='tlooby-cent', debug = False)
    log.info("Initializing Log File...")
#Uncomment this to run as locahost
    app.run(debug=True)
    #Logfile for troubleshooting and for passing back to html page

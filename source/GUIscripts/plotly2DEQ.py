#visualizationClass.py
#Description:   Base HEAT Viz module
#Engineer:      T Looby
#Date:          20200513

import sys
import os

import numpy as np
import MDSplus
import EFIT.equilParams_class as EP
from scipy import interpolate
from scipy.interpolate import interp1d
import json
import matplotlib.pyplot as plt
import logging


def makePlotlyEQDiv(shot, time, MachFlag, ep, gfile=None, logFile=False):
    """
    returns a DASH object for use directly in dash app
    """
    if logFile is True:
        log = logging.getLogger(__name__)

    #Use Equilparamsclass to create EQ object
    if ep is None:
        print('Note:  no EP object')
        if gfile is None:
            print("Error generating EQ plot: no EQ object or gfile")
            sys.exit()
        else:
            ep = EP.equilParams(gfile)

    r = ep.g['R']
    z = ep.g['Z']
    psi = ep.g['psiRZn']
    R,Z = np.meshgrid(r, z)


    rbdry = ep.g['lcfs'][:,0]
    zbdry = ep.g['lcfs'][:,1]

    if MachFlag == 'nstx':
        rlim, zlim = nstxu_wall(oldwall=False) #FOR NSTXU
    else:
        rlim = ep.g['wall'][:,0]
        zlim = ep.g['wall'][:,1]


    levels = sorted(np.append([0.0,0.05,0.1,0.25,0.5,0.75,1.0], np.linspace(1.01,psi.max(),10)))

#    aspect = (z.max()-z.min()) / (r.max()-r.min())
#    width = (1.0/aspect)*height

    import plotly
    import plotly.graph_objects as go
    import plotly.express as px
    #psi data
    fig = go.Figure(data =
        go.Contour(
            z=psi,
            x=r, # horizontal axis
            y=z, # vertical axis
            colorscale='cividis',
            contours_coloring='heatmap',
            name='psiN',
            showscale=False,
        ))

    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    )
            )
            )

    #Seperatrix in red.  Sometimes this fails if psi is negative
    #so we try and except.
    #if try fails, just plot using rbdry,zbdry from gfile
    try:
        CS = plt.contourf(R,Z,psi,levels,cmap=plt.cm.cividis)
        lcfsCS = plt.contour(CS, levels = [1.0])
        for i in range(len(lcfsCS.allsegs[0])):
            rlcfs = lcfsCS.allsegs[0][i][:,0]
            zlcfs = lcfsCS.allsegs[0][i][:,1]

            fig.add_trace(
                go.Scatter(
                    x=rlcfs,
                    y=zlcfs,
                    mode="lines",
                    name="LCFS",
                    line=dict(
                        color="red",
                        width=4,

                            )
                    )
                    )
    except:
        print("Could not create contour plot.  Psi levels must be increasing.")
        print("Try flipping psi sign and replotting.")
        print("plotting rbdry, zbdry from gfile (not contour)")
        if logFile is True:
            log.info("Could not create contour plot.  Psi levels must be increasing.")
            log.info("Try flipping psi sign and replotting.")
            log.info("plotting rbdry, zbdry from gfile (not contour)")

        fig.add_trace(
            go.Scatter(
                x=rbdry,
                y=zbdry,
                mode="lines",
                name="LCFS",
                line=dict(
                    color="red",
                    width=4,

                        )
                )
                )


    fig.update_layout(
        title="{:06d}@{:05d}ms".format(shot,time),
        xaxis_title="R [m]",
        yaxis_title="Z [m]",
        autosize=True,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False,
        font=dict(
#            family="Courier New",
            size=18,
            color="#dcdce3"
        )
        )
    return fig







def writePlotlyEQ(shot, time, outFile, MachFlag, ep=None, gfile=None, logFile=False):
    """
    saves a plotly webpage to a file that can be imported to html via iframe
    """
    if logFile is True:
        log = logging.getLogger(__name__)

    #Use Equilparamsclass to create EQ object
    if ep is None:
        if gfile is None:
            print("Error generating EQ plot: no EQ object or gfile")
            sys.exit()
        else:
            ep = EP.equilParams(gfile)

    r = ep.g['R']
    z = ep.g['Z']
    psi = ep.g['psiRZn']
    R,Z = np.meshgrid(r, z)


    rbdry = ep.g['lcfs'][:,0]
    zbdry = ep.g['lcfs'][:,1]

    if MachFlag == 'nstx':
        rlim, zlim = nstxu_wall(oldwall=False) #FOR NSTXU
    else:
        rlim = ep.g['wall'][:,0]
        zlim = ep.g['wall'][:,1]


    levels = sorted(np.append([0.0,0.05,0.1,0.25,0.5,0.75,1.0], np.linspace(1.01,psi.max(),10)))

    import plotly
    import plotly.graph_objects as go
    import plotly.express as px

    #psi data
    fig = go.Figure(data =
        go.Contour(
            z=psi,
            x=r, # horizontal axis
            y=z, # vertical axis
            colorscale='cividis',
            contours_coloring='heatmap',
            name='psiN',
            showscale=False,
        ))

    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    )
            )
            )

    #Seperatrix in red.  Sometimes this fails if psi is negative
    #so we try and except.
    #if try fails, just plot using rbdry,zbdry from gfile
    try:
        CS = plt.contourf(R,Z,psi,levels,cmap=plt.cm.cividis)
        lcfsCS = plt.contour(CS, levels = [1.0])
        for i in range(len(lcfsCS.allsegs[0])):
            rlcfs = lcfsCS.allsegs[0][i][:,0]
            zlcfs = lcfsCS.allsegs[0][i][:,1]

            fig.add_trace(
                go.Scatter(
                    x=rlcfs,
                    y=zlcfs,
                    mode="lines",
                    name="LCFS",
                    line=dict(
                        color="red",
                        width=4,

                            )
                    )
                    )
    except:
        print("Could not create contour plot.  Psi levels must be increasing.")
        print("Try flipping psi sign and replotting.")
        print("plotting rbdry, zbdry from gfile (not contour)")
        if logFile is True:
            log.info("Could not create contour plot.  Psi levels must be increasing.")
            log.info("Try flipping psi sign and replotting.")
            log.info("plotting rbdry, zbdry from gfile (not contour)")

        fig.add_trace(
            go.Scatter(
                x=rbdry,
                y=zbdry,
                mode="lines",
                name="LCFS",
                line=dict(
                    color="red",
                    width=4,

                        )
                )
                )


    fig.update_layout(
        title="{:06d}@{:05d}ms".format(shot,time),
        xaxis_title="R [m]",
        yaxis_title="Z [m]",
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=False,
        font=dict(
            family="Anurati",
            size=18,
            color="#dcdce3"
        )
        )

    print("writing EQ output to: "+outFile)
    fig.write_html(outFile)


def nstxu_wall(oldwall=False):
    """
    returns simplified wall.  Uses two different wall versions
    """
    if oldwall:
        R = np.array([0.1851, 0.1851, 0.2794, 0.2794, 0.2979, 0.5712,
                    1.0433, 1.3192, 1.3358,
                    1.4851, 1.4791, 1.5174, 1.5313, 1.5464, 1.5608,
                    1.567, 1.5657, 1.5543, 1.5341, 1.5181, 1.4818,
                    1.4851, 1.3358, 1.3192, 1.0433,
                    0.5712, 0.2979, 0.2794, 0.2794, 0.1851, 0.1851])
        Z = np.array([0.0, 1.0081, 1.1714, 1.578, 1.6034, 1.6034,
                    1.43, 1.0397, 0.9976,
                    0.545, 0.4995, 0.306, 0.2355, 0.1586, 0.0801,
                    0.0, -0.0177, -0.1123, -0.221, -0.3026, -0.486,
                    -0.545, -0.9976, -1.0397, -1.43,
                    -1.6034, -1.6034, -1.578, -1.1714, -1.0081, 0])
    else:
      R = np.array([ 0.3147568,  0.3147568,  0.4441952,  0.4441952,  0.443484 ,
           0.443484 ,  0.6000496,  0.7672832,  0.8499856,  1.203452,  1.3192,  1.3358,  1.4851,  1.489 ,
           1.5638,  1.57  ,  1.5737,  1.575 ,  1.5737,  1.57  ,  1.5638,
           1.489 ,  1.4851,  1.3358,  1.3192,  1.203452 ,  0.8499856,  0.7672832,  0.6000496,  0.443484 ,
           0.443484 ,  0.4441952,  0.4441952,  0.3147568,  0.3147568 ])
      Z = np.array([ 0.       ,  1.0499344,  1.2899136,  1.5104872,  1.5104872,
            1.6028416,  1.6028416,  1.5367   ,  1.5367   ,  1.397508,  1.0397,  0.9976,  0.545 ,  0.49  ,
            0.1141,  0.0764,  0.0383,  0.    , -0.0383, -0.0764, -0.1141,
            -0.49  , -0.545 , -0.9976, -1.0397, -1.397508 , -1.5367   , -1.5367   , -1.6028416, -1.6028416,
            -1.5104872, -1.5104872, -1.2899136, -1.0499344,  0.])
    return R,Z



if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('You must call this function with 4 input arguments:')
        print('-gfile path (absolute)')
        print('-shot number')
        print('-timestep in milliseconds')
        print('-output file name\n')
        print('Example Command:  python3 plotlyEQ.py ~/test/g0000001.000002 1 2 output.html\n')
        sys.exit()

    else:
        gfile = str(sys.argv[1])
        shot = int(sys.argv[2])
        time = int(sys.argv[3])
        outFile = str(sys.argv[4])

    writePlotlyEQ(gfile, shot, time, outFile)

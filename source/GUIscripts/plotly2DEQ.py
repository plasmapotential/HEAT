#plotly2DEQ.py
#Description:   Base HEAT equilibrium plot
#Engineer:      T Looby
#Date:          20200513

import sys
import os

import numpy as np
import EFIT.equilParams_class as EP
from scipy import interpolate
from scipy.interpolate import interp1d
import json
import matplotlib.pyplot as plt
import logging
from dash_bootstrap_templates import template_from_url


def makePlotlyEQDiv(shot, time, MachFlag, ep, height=None, gfile=None,
                    logFile=False, bg = None, xRange=None, yRange=None,
                    tsFmt="{:.6f}", shotFmt="{:06d}"):
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

# for aspect ratio
#    if height==None:
#        height=1000
    if height is not None:
        aspect = (z.max()-z.min()) / (r.max()-r.min())
        width = (1.0/aspect)*height



    levels = sorted(np.append([0.0,0.05,0.1,0.25,0.5,0.75,0.95, 1.0], np.linspace(0.99,psi.max(),15)))

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
            ncontours=20,
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
                    ),
            )
            )

    #white lines around separatrix
    levelsAtLCFS = np.linspace(0.95,1.05,15)
    CS = plt.contourf(R,Z,psi,levelsAtLCFS,cmap=plt.cm.cividis)
    for i in range(len(levelsAtLCFS)):
        levelsCS = plt.contour(R,Z,psi,levels=[levelsAtLCFS[i]])
        for j in range(len(levelsCS.allsegs[0])):
            r = levelsCS.allsegs[0][j][:,0]
            z = levelsCS.allsegs[0][j][:,1]
            fig.add_trace(
                go.Scatter(
                    x=r,
                    y=z,
                    mode="lines",
                    line=dict(
                        color="white",
                        width=1,
                        dash='dot',
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

    #set bkgrd color if not default (clear)
#    if bg is None:
#        bg = 'rgba(0,0,0,0)'


    fig.update_layout(
        title=shotFmt.format(shot)+"@"+tsFmt.format(time)+"s",
        xaxis_title="R [m]",
        yaxis_title="Z [m]",
        autosize=True,
        #for aspect ratio
        #autosize=False,
        #width=width*1.1,
        #height=height,
#        paper_bgcolor=bg,
#        plot_bgcolor=bg,
        showlegend=False,
        font=dict(
#            family="Courier New",
            size=18,
#            color="primary"
        ),
        margin=dict(
            l=10,
            r=10,
            b=10,
            t=100,
            pad=4
        )
        )

    fig.update_yaxes(scaleanchor = "x",scaleratio = 1,)

    #set height if we use this for a png image
    if height is not None:
        fig.update_layout(height=height, width=width)

    #set ranges if not none
    if xRange is not None:
        fig.update_layout(xaxis_range=xRange)
    if yRange is not None:
        fig.update_layout(yaxis_range=yRange)

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


    levels = sorted(np.append([0.0,0.05,0.1,0.25,0.5,0.75,1.0], np.linspace(1.01,psi.max(),20)))

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


def makePlotlyGfilePlot(ep):
    """
    plots Fpol, psiRZ (along midplane), psiSep, psiAxis, all on one plot

    ep is equilibrium object (from equilParams_class)
    """
    #assuming plasma is centered in machine here
    zMin = ep.g['ZmAxis'] - 0.25
    zMax = ep.g['ZmAxis'] + 0.25
    zWall = np.linspace(zMin, zMax, 1000)
    zLCFS = ep.g['lcfs'][:,1]
    #this prevents us from getting locations not at midplane
    idx = np.where(np.logical_and(zLCFS>zMin,zLCFS<zMax))
    RmaxLCFS = ep.g['lcfs'][:,0][idx].max()
    RminLCFS = ep.g['lcfs'][:,0][idx].min()
    rAxis = ep.g['RmAxis']
    rMax = max(ep.g['R'])
    zMid = ep.g['ZmAxis']
    Fpol = ep.g['Fpol']
    N = len(Fpol)
    rN = np.linspace(rAxis, RmaxLCFS, N)
    psiAxis = ep.g['psiAxis']
    psiSep = ep.g['psiSep']
    zMidplane = np.ones(len(ep.g['R']))*zMid

    #psiRZN = ep.psiFunc.ev(ep.g['R'],zMidplane)
    psiRZ = ep.g['psiRZ'][int(ep.g['NZ']/2)][:]
    psiRZN = (psiRZ - psiAxis)/(psiSep - psiAxis)

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(go.Scatter(x=rN, y=Fpol,
                    mode='lines',
                    name='Fpol'),
                    secondary_y=True)
    fig.add_trace(go.Scatter(x=ep.g['R'], y=psiRZN,
                    mode='lines',
                    name='psiN'))
    fig.add_trace(go.Scatter(x=ep.g['R'], y=psiRZ,
                    mode='lines',
                    name='psi'))
    fig.add_trace(go.Scatter(x=[rAxis], y=[psiAxis],
                    mode='markers',
                    name='psiAxis'))
    fig.add_trace(go.Scatter(x=[RmaxLCFS], y=[psiSep],
                    mode='markers',
                    name='psiSep'))

    fig.update_layout(
    title="GEQDSK Parameters at Midplane",
    xaxis_title="R [m]",
    font=dict(
        family="Arial",
        size=16,
        )
    )
    fig.update_yaxes(title_text="psi", secondary_y=False)
    fig.update_yaxes(title_text="Fpol", secondary_y=True)
    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01
    ))
    return fig





def makePlotlyBpBt(ep, MachFlag, logFile=False):
    """
    returns a DASH object for use directly in dash app
    """

    if logFile is True:
        log = logging.getLogger(__name__)

    r = ep.g['R']
    z = ep.g['Z']
    R,Z = np.meshgrid(r, z)

    #Bp = ep.BpFunc.ev(R,Z)
    #Bt = ep.BtFunc.ev(R,Z)
    Bp = ep.Bp_2D
    Bt = ep.Bt_2D
    Br = ep.B_R
    Bz = ep.B_Z

    if MachFlag == 'nstx':
        rlim, zlim = nstxu_wall(oldwall=False) #FOR NSTXU
    else:
        rlim = ep.g['wall'][:,0]
        zlim = ep.g['wall'][:,1]

    #dont waste contour space on xfmr coil field if its way higher than Bt0
    BtMax = 3*ep.g['Bt0']
    if ep.g['Bt0'] < 0:
        BtMin = np.max(Bt)
    else:
        BtMin = np.min(Bt)

    import plotly
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.05, vertical_spacing=0.05, shared_yaxes=True,
                    subplot_titles=("Bt [T]", "Bp [T]", "Br [T]", "Bz [T]"))

    fig.add_trace(
        go.Contour(
            z=Bt,
            x=r, # horizontal axis
            y=z, # vertical axis
            #colorscale='cividis',
            contours_coloring='heatmap',
            name='Bt',
            showscale=False,
            ncontours=30,
            contours=dict(
                start=BtMin,
                end=BtMax,
                ),
            ),
        row=1,
        col=1
        )
    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    ),
            ),
        row=1,
        col=1
        )

    fig.add_trace(
        go.Contour(
            z=Bp,
            x=r, # horizontal axis
            y=z, # vertical axis
            colorscale='viridis',
            contours_coloring='heatmap',
            name='Bp',
            showscale=False,
            ncontours=40,
            ),
        row=1,
        col=2
        )

    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    ),
            ),
        row=1,
        col=2
        )



    fig.add_trace(
        go.Contour(
            z=Br,
            x=r, # horizontal axis
            y=z, # vertical axis
            colorscale='Purples',
            contours_coloring='heatmap',
            name='Br',
            showscale=False,
            ncontours=30,
            ),
        row=2,
        col=1
        )
    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    ),
            ),
        row=2,
        col=1
        )

    fig.add_trace(
        go.Contour(
            z=Bz,
            x=r, # horizontal axis
            y=z, # vertical axis
            colorscale='Plotly3',
            contours_coloring='heatmap',
            name='Bz',
            showscale=False,
            ncontours=30,
            ),
        row=2,
        col=2
        )
    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    ),
            ),
        row=2,
        col=2
        )






    fig.update_layout(showlegend=False,
        margin=dict(
        l=5,
        r=5,
        b=30,
        t=50,
        pad=2,
        ),
        height=800,
    )
    return fig





def makePlotlyBrBz(ep, MachFlag, logFile=False):
    """
    returns a DASH object for use directly in dash app
    """

    if logFile is True:
        log = logging.getLogger(__name__)

    r = ep.g['R']
    z = ep.g['Z']
    R,Z = np.meshgrid(r, z)

    #Bp = ep.BpFunc.ev(R,Z)
    #Bt = ep.BtFunc.ev(R,Z)
    Br = ep.B_R
    Bz = ep.B_Z

    if MachFlag == 'nstx':
        rlim, zlim = nstxu_wall(oldwall=False) #FOR NSTXU
    else:
        rlim = ep.g['wall'][:,0]
        zlim = ep.g['wall'][:,1]

    import plotly
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=2, cols=1, horizontal_spacing=0.05, shared_yaxes=True,
                    subplot_titles=("Br [T]", "Bz [T]"))

    fig.add_trace(
        go.Contour(
            z=Br,
            x=r, # horizontal axis
            y=z, # vertical axis
            #colorscale='cividis',
            contours_coloring='heatmap',
            name='Br',
            showscale=False,
            ncontours=30,
            ),
        row=1,
        col=1
        )
    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    ),
            ),
        row=1,
        col=1
        )

    fig.add_trace(
        go.Contour(
            z=Bz,
            x=r, # horizontal axis
            y=z, # vertical axis
            #colorscale='cividis',
            contours_coloring='heatmap',
            name='Bz',
            showscale=False,
            ncontours=30,
            ),
        row=1,
        col=2
        )
    #Wall in green
    fig.add_trace(
        go.Scatter(
            x=rlim,
            y=zlim,
            mode="markers+lines",
            name="Wall",
            line=dict(
                color="#19fa1d"
                    ),
            ),
        row=1,
        col=2
        )






    fig.update_layout(showlegend=False,
        margin=dict(
        l=5,
        r=5,
        b=30,
        t=50,
        pad=2,
        ),
    )
    return fig











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

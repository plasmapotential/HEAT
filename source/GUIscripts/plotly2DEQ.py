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

def writePlotlyEQ(shot, time, outFile, ep=None, gfile=None):

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

    rlim = ep.g['wall'][:,0]
    zlim = ep.g['wall'][:,1]


    levels = sorted(np.append([0.0,0.05,0.1,0.25,0.5,0.75,1.0], np.linspace(1.01,psi.max(),10)))
    CS = plt.contourf(R,Z,psi,levels,cmap=plt.cm.cividis)
    lcfsCS = plt.contour(CS, levels = [1.0])

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

    #Seperatrix in red
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

#plotlyGUIplots.py
#Description:   Base HEAT GUI plots that are outputs
#Engineer:      T Looby
#Date:          20200727
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
import plotly.express as px

def plotlyOpenFOAMplot(data,names):
    """
    returns a DASH object for use directly in dash app

    data is a list of dataframes, where each list item corresponds to a different
    PFC, whose names are in names

    here we plot the output from min/max T against the min/max HF
    """
    if type(data) is not list:
        data = [data]
    if type(names) is not list:
        names = [names]

    #if we want multiple fields
#    fields = ['T', 'qFVM', 'qINV']
#    units = [' [K]', ' [W/m^2]', ' [W/m^2]']
#    y2 = [False, True, True]
    fields = ['T', 'qFVM']
    units = [' [K]', ' [W/m^2]']
    y2 = [False, True, True]

    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])

    for pfcIdx,name in enumerate(names):
        df = data[pfcIdx]
        for i,field in enumerate(fields):
            mask = df['field'] == field
            t = df[mask].sort_values('# Time')['# Time'].values
            varMax = df[mask].sort_values('# Time')['max'].values
            lineName = field+" "+name
            fig.add_trace(go.Scatter(x=t, y=varMax, name=lineName),
                            secondary_y=y2[i],
                          )


    fig.update_layout(title='Time History of openFOAM FVM Variables',
                xaxis_title='<b>Time [s]</b>',
                yaxis_title='',
                font=dict(
                    family='Arial, sans-serif',
                    size=14,
                    )
                    )
    # Set y-axes titles
    fig.update_yaxes(title_text="<b>Maximum PFC Heat Flux [W/m^2]</b>", secondary_y=True)
    fig.update_yaxes(title_text="<b>Maximum PFC Temperature [K]</b>", secondary_y=False)

    return fig


def plotlyqDivPlot(heatFlux, labels, logPlot=False):
    """
    returns a DASH object for use directly in dash app

    heatFlux is an arrays of heat fluxes on PFC surface

    heatFlux will be filtered so no zero points are in array

    labels is a list of same length as heatFlux

    if logPlot is true plot a log plot


    if heatFlux and labels are lists, then we add a plot for each item in list to figure
    this is useful when running multiple PFCs
    """
    fig = go.Figure()
    for i,hf in enumerate(heatFlux):
        use = np.where(hf>0)
        hf = hf[use]
        label = labels[i]

        if logPlot is True:
            hflog = np.log10(hf)

#        if i==0:
#            fig = px.histogram(x=hflog,
#                               marginal="box",
#                               nbins=20,
#                               opacity=0.7,
#                               labels=label,
#                               )
#        else:
#            pass
        fig.add_trace(go.Histogram(
                            x=hflog,
                            #marginal="box",
                            nbinsx=20,
                            opacity=0.7,
                            name=label,
                            )
                        )

    fig.update_layout(
        title="qDiv Distribution Across PFC Surface",
        xaxis_title="log10(qDiv) [MW/m^2]",
        yaxis_title="# of Points (Count)",
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
            ),
        font=dict(
            family="Arial",
            size=18,
            color="Black"
            ),
        )

    return fig



def plotlyTprobes(t,T,names):
    """
    returns a DASH object for use directly in dash app

    t is list of timesteps

    T is list of temperature at each corresponding t item

    names is name of PFC for that T (element-wise)

    both t and T are numpy arrays of the same length at each list index
    """
    if type(t) is not list:
        t = [t]
    if type(T) is not list:
        T = [T]
    if type(names) is not list:
        names = [names]


    fig = go.Figure()
    for i,T in enumerate(T):
        name = 'T{:d} '.format(i) + names[i]
        fig.add_trace(go.Scatter(x=t[i], y=T, name=name))

    xMin = min([min(arr) for arr in t])
    xMax = max([max(arr) for arr in t])
    fig.add_trace(go.Scatter(
        x=[xMin, xMax],
        y=[1873, 1873],
        mode="lines+markers+text",
        name="Sublimation T",
        text=["Limit", "Limit"],
        textposition="top center",
        line=dict(color='firebrick', width=3, dash='dot'),
        textfont=dict(family="Arial", size=16, color="firebrick"),

    ))

    fig.update_layout(
        title="Temperature Probe Time Evolution",
        xaxis_title="Time [s]",
        yaxis_title="Temperature [K]",
        font=dict(
        family="Arial",
        size=18,
        color="Black"
    ),
    )



    return fig


def plotlyGyroPhasePlot(gyroPhases):
    """
    returns a DASH object for use directly in dash app

    gyroPhase is an array of gyroPhases

    """
    gyroTicks = [0, 45, 90, 135, 180, 225, 270, 315]
    fig = go.Figure()
    for i in range(len(gyroPhases)):
        fig.add_trace(go.Scatterpolar(
                r = np.array([0.0,1.0]),
                theta = np.array([gyroPhases[i],gyroPhases[i]]),
                mode = 'lines+markers',
                marker={'symbol':"circle", 'size':16},
                name = '{:0.1f}'.format(gyroPhases[i]),
                text=['{:0.1f}'.format(gyroPhases[i])],
                textposition="bottom center",
                #line_color= px.colors.qualitative.Dark2[0],
                line={'width':7},
                ))

    fig.update_layout(
        title={'text': "Gyro Phase Angles",'y':0.94,'x':0.5,'xanchor': 'center','yanchor': 'top'},
        showlegend = False,
        polar = dict(radialaxis=dict(visible=False), angularaxis=dict(tickvals=gyroTicks)),
        )

    return fig

def plotlyVPhasePlot(vPhases):
    """
    returns a DASH object for use directly in dash app
    vPhase is an array of velocity phases
    """

    fig = go.Figure()
    for i in range(len(vPhases)):
        fig.add_trace(go.Scatterpolar(
                r = np.array([0.0,1.0]),
                theta = np.array([vPhases[i],vPhases[i]]),
                mode = 'lines+markers',
                marker={'symbol':"square", 'size':16},
                name = '{:0.1f}'.format(vPhases[i]),
                text=['{:0.1f}'.format(vPhases[i])],
                textposition="bottom center",
                #line_color= px.colors.qualitative.Dark2[0],
                line={'width':7},
                ))



    fig.add_annotation(x=0.2, y=0.5,
                text="V||",
                font=dict(size=16),
                showarrow=False,
                )

    fig.update_annotations(font_size=16)

    fig.update_layout(
        title={'text': "Velocity Phase Angles",'y':0.94,'x':0.5,'xanchor': 'center','yanchor': 'top'},
        showlegend = False,
        polar = dict(
            sector = [0,90],
            radialaxis=dict(title=dict(text="V\u22A5",font=dict(size=16))),
            #angularaxis=dict(title=dict(text='$V_{||}$',font=dict(size=24)))
            ),
            )
    return fig


def plotlyVSlicePlot(m,c,T0,vSlices,v):
    """
    returns a DASH object for use directly in dash app
    m is mass [eV]
    c is speed of light [m/s]
    T0 is temperature [eV]
    vSlices is slice velocities [m/s]
    v is all velocities in distribution scan (ie vScan) [m/s]
    """

    #generate the (here maxwellian) PDF
    pdf = lambda x: (m/c**2) / (T0) * np.exp(-(m/c**2 * x**2) / (2*T0) )
    v_pdf = v * pdf(v)
    vSlice_pdf = vSlices * pdf(vSlices)



    fig = go.Figure()
    fig.add_trace(go.Scatter(x=v, y=v_pdf,
                        mode='lines',
                        line={'width':6},
                        name='Maxwellian PDF'))

    fig.add_trace(go.Scatter(x=vSlices, y=vSlice_pdf,
                        mode='markers',
                        marker={'symbol':"square", 'size':16},
                        name='Slices'))


    fig.update_layout(
        title="Velocity Distribution (vSlices)",
        yaxis= dict(showticklabels=False),
        xaxis_title="Velocity",
        font=dict(size=18),
        )

    return fig

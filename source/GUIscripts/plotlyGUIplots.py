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


def plotlyVSlicePlot(m,c,T0,vSlices,vBounds,v):
    """
    returns a DASH object for use directly in dash app
    m is mass [eV]
    c is speed of light [m/s]
    T0 is temperature [eV]
    vSlices is slice velocities [m/s]
    v is all velocities in distribution scan (ie vScan) [m/s]
    """

    #generate the (here maxwellian) PDF
    pdf = lambda x: ( (m/c**2) / (2 * np.pi * T0) )**(3.0/2.0) * np.exp(-(m/c**2 * x**2) / (2*T0) )
    v_pdf = 4*np.pi * v**2 * pdf(v)
    vSlice_pdf = 4*np.pi * vSlices**2 * pdf(vSlices)
    vThermal = np.sqrt(2.0*T0/(m/c**2))
    vMax = 5.0 * vThermal

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=v, y=v_pdf,
                        mode='lines',
                        line={'width':6},
                        name='Maxwellian PDF'))

    fig.add_trace(go.Scatter(x=vSlices, y=vSlice_pdf,
                        mode='markers',
                        marker={'symbol':"square", 'size':16},
                        name='Slices'))

    #add bin boundaries
    #for i in range(len(vBounds)):
    #    fig.add_vline(x=vBounds[i], line_width=3, line_dash="dash", line_color="green")


    fig.update_layout(
        title="Velocity Distribution (vSlices)",
        yaxis= dict(showticklabels=False),
        yaxis_range=[0,max(v_pdf)*1.05],
        xaxis_range=[0,vMax],
        xaxis_title="Velocity [m/s]",
        font=dict(size=18),
        )

    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.95,
        xanchor="right",
        x=0.95
        ),
        )

    return fig



def plotlycdfSlicePlot(m,c,T0,vSlices,vBounds,v,N_vSlice):
    """
    returns a DASH object for use directly in dash app
    m is mass [eV]
    c is speed of light [m/s]
    T0 is temperature [eV]
    vSlices is slice velocities [m/s]
    v is all velocities in distribution scan (ie vScan) [m/s]
    """

    #generate the (here maxwellian) PDF
    pdf = lambda x: ( (m/c**2) / (2 * np.pi * T0) )**(3.0/2.0) * np.exp(-(m/c**2 * x**2) / (2*T0) )
    v_pdf = 4*np.pi * v**2 * pdf(v)
    vThermal = np.sqrt(2.0*T0/(m/c**2))
    vMax = 5.0 * vThermal

    #generate the CDF
    v_cdf = np.cumsum(v_pdf[1:])*np.diff(v)
    v_cdf = np.insert(v_cdf, 0, 0)
    #CDF location of vSlices and bin boundaries
    cdfBounds = np.linspace(0,v_cdf[-1],N_vSlice+1)
    cdfSlices = np.diff(cdfBounds)/2.0 + cdfBounds[:-1]


    fig = go.Figure()
    fig.add_trace(go.Scatter(x=v, y=v_cdf,
                        mode='lines',
                        line={'width':6, 'color':'purple'},
                        name='Maxwellian CDF'))

    fig.add_trace(go.Scatter(x=vSlices, y=cdfSlices,
                        mode='markers',
                        marker={'symbol':"square", 'size':16},
                        name='Slices'))

    for i in range(len(cdfBounds)):
        fig.add_vline(x=vBounds[i], line_width=3, line_dash="dash", line_color="green")


    fig.update_layout(
        title="Cumulative Distribution Function",
        yaxis= dict(showticklabels=True),
        yaxis_range=[0,1.05],
        xaxis_range=[0,vMax],
        xaxis_title="Velocity [m/s]",
        font=dict(size=18),
        )

    fig.update_layout(legend=dict(
        yanchor="bottom",
        y=0.05,
        xanchor="right",
        x=0.95
        ),
        )

    return fig


def plotlyallSlicePlot(m,c,T0,vSlices,vBounds,v,N_vSlice):
    """
    returns a DASH object for use directly in dash app
    m is mass [eV]
    c is speed of light [m/s]
    T0 is temperature [eV]
    vSlices is slice velocities [m/s]
    v is all velocities in distribution scan (ie vScan) [m/s]
    """

    #generate the (here maxwellian) PDF
    pdf = lambda x: ( (m/c**2) / (2 * np.pi * T0) )**(3.0/2.0) * np.exp(-(m/c**2 * x**2) / (2*T0) )
    v_pdf = 4*np.pi * v**2 * pdf(v)
    vThermal = np.sqrt(2.0*T0/(m/c**2))
    vMax = 5.0 * vThermal
    vSlice_pdf = 4*np.pi * vSlices**2 * pdf(vSlices)
    #generate the CDF
    v_cdf = np.cumsum(v_pdf[1:])*np.diff(v)
    v_cdf = np.insert(v_cdf, 0, 0)
    #CDF location of vSlices and bin boundaries
    cdfBounds = np.linspace(0,v_cdf[-1],N_vSlice+1)
    cdfSlices = np.diff(cdfBounds)/2.0 + cdfBounds[:-1]




    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.05, x_title="Velocity [m/s]",)

    #PDF Traces
    fig.add_trace(go.Scatter(x=v, y=v_pdf,
                        mode='lines',
                        line={'width':6},
                        name='Maxwellian PDF',
                        showlegend=True),
                        row=1, col=1)

    fig.add_trace(go.Scatter(x=vSlices, y=vSlice_pdf,
                        mode='markers',
                        marker={'symbol':"circle", 'size':16, "color":"black"},
                        name='PDF Slices'),
                        row=1, col=1)

    #CDF Traces
    fig.add_trace(go.Scatter(x=v, y=v_cdf,
                        mode='lines',
                        line_dash="dashdot",
                        line={'width':6, 'color':'purple'},
                        name='Maxwellian CDF',
                        showlegend=True),
                        row=2, col=1)

    fig.add_trace(go.Scatter(x=vSlices, y=cdfSlices,
                        mode='markers',
                        marker={'symbol':"cross", 'size':16, "color":"darkgreen"},
                        name='CDF Slices'),
                        row=2, col=1)

    #Bin boundaries (vertical for PDF, horizontal for CDF)
    for i in range(len(cdfBounds)):
        #fig.add_vline(dict(x=vBounds[i], line_width=3, line_dash="dash", line_color="green"), row=2, col=1)
        fig.add_shape(dict(type="line", x0=vBounds[i], x1=vBounds[i], y0=0, y1=max(v_pdf)*1.05, line_color="firebrick", line_width=3, line_dash="dot"), row=1, col=1)
        fig.add_shape(dict(type="line", x0=0, x1=vMax, y0=cdfBounds[i], y1=cdfBounds[i], line_color="firebrick", line_width=3, line_dash="dot"), row=2, col=1)


    fig.update_layout(
        title="Velocity PDF vs. CDF",
#        yaxis= dict(showticklabels=True),
#        yaxis_range=[0,1.05*max(v_pdf)],
#        xaxis_range=[0,vMax],
        font=dict(size=18),
        )

    fig.update_yaxes(range=[-1e-6,max(v_pdf)*1.05], showticklabels=False, row=1,col=1)
    fig.update_yaxes(range=[-0.05,1.05], showticklabels=True, row=2,col=1)

    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.98,
        xanchor="right",
        x=0.95,
        font=dict(size=14),
        ),
        )
    #fig.update_layout(showlegend=False)

    return fig

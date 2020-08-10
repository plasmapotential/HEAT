#plotlyGUIplots.py
#Description:   Base HEAT GUI plots that are outputs
#Engineer:      T Looby
#Date:          20200727
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff

def plotlyOpenFOAMplot(data,names):
    """
    returns a DASH object for use directly in dash app

    data is a list of dataframes, where each list item corresponds to a different
    PFC, whose names are in names

    here we plot the output from min/max T against the min/max HF
    """
    if type(data) is not list:
        data = [data]



    fields = ['T', 'qFVM', 'qINV']
    units = [' [k]', ' [W/m^2]', ' [W/m^2]']
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


    fig.update_layout(title='Time History of openFOAM FVM Temperature',
                xaxis_title='<b>Time [s]</b>',
                yaxis_title='',
                font=dict(
                    family='Arial, sans-serif',
                    size=24,
                    )
                    )
    # Set y-axes titles
    #fig.update_yaxes(title_text="<b>Maximum PFC Heat Flux [W/m^2]</b>", secondary_y=True)
    fig.update_yaxes(title_text="<b>Maximum PFC Temperature [K]</b>", secondary_y=False)

    return fig


def plotlyqDivPlot(hf, label):
    """
    returns a DASH object for use directly in dash app

    hf is a list of arrays of heat fluxes.  returns
    a histogram for each list item.  This can be used to pass multiple
    PFC qDivs and have them all shown on the same plot.

    hf should not contain any zero points!

    labels is a list of same length as hf
    """
    if type(hf) is not list:
        hf = [hf]
    if type(label) is not list:
        label = [label]


    fig = ff.create_distplot(hf, label, curve_type='normal', show)
    fig.update_layout(xaxis_type="log", yaxis_type="log")
    fig = px.histogram(df, y="HeatFlux",marginal="box")
    return fig

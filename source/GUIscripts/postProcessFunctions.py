#helper script to do plots and other postprocessing functions
#usually I just edit as needed, then copy and paste into a python terminal


#===============================================================================
#                       Plot Tprobe Data
#===============================================================================

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
#file="T"
#data = np.genfromtxt(file,comments="#", autostrip=True)
#fig = px.line(x=data[:,0],y=data[:,1])

root = '/home/tom/results/WGPFC_Validation/staticScans/memo010_1.1/nstx_116313_LwrOut/openFoam/heatFoam/'
file1 = root+'SOLID844/postProcessing/probes/0/T'
file2 = root+'E-ED1408-284/postProcessing/probes/0/T'

data1 = np.genfromtxt(file1,comments="#", autostrip=True)
data2 = np.genfromtxt(file2,comments="#", autostrip=True)
fig = go.Figure()
fig.add_trace(go.Scatter(x=data1[:,0], y=data1[:,1], name="IBDH", line=dict(color='royalblue', width=6, dash='dot')))
fig.add_trace(go.Scatter(x=data2[:,0],y=data2[:,1],name="OBD", line=dict(color='magenta', width=6)))
fig.add_trace(go.Scatter(
    x=[0, 6],
    y=[1873, 1873],
    mode="lines+markers+text",
    name="Limit",
    text=["Limit", "Limit"],
    textposition="top center",
    line=dict(color='firebrick', width=3, dash='dot'),
    textfont=dict(family="Arial", size=18, color="firebrick"),

))

fig.update_layout(
#title="Temperature Probe Time Evolution",
xaxis_title="Time [s]",
yaxis_title="Temperature [K]",
font=dict(
    family="Arial",
    size=24,
    color="Black"
),
margin=dict(
    l=5,
    r=5,
    b=5,
    t=5,
    pad=2
),
)

fig.update_layout(legend=dict(
    yanchor="top",
    y=0.99,
    xanchor="left",
    x=0.01
))

fig.add_annotation(x=3.0, y=1000,
            text="IBDH T(5s) = 2533 K",
            showarrow=False,
            yshift=20,
            font=dict(
            family="Arial",
            size=24,
            color="royalblue",
            ),
            )

fig.add_annotation(x=3.0, y=850,
            text=" OBD T(5s) = 1935 K",
            showarrow=False,
            yshift=20,
            font=dict(
            family="Arial",
            size=24,
            color="magenta",
            ),
            )



fig.show()

fig.write_image("/home/tom/phd/papers/2020Fall_FST_HEATintro/HEATintro/figures/hf_case1_2.svg")


import pandas as pd
data1 = pd.read_csv('qDiv_E-ED1408-284.csv')
data2 = pd.read_csv('qDiv_E-ED1408-573.csv')
data3 = pd.read_csv('qDiv_SOLID844.csv')
hfs = [data1['HeatFlux'].values, data2['HeatFlux'].values, data3['HeatFlux'].values]
nombres = ['E-ED1408-284','E-ED1408-573','SOLID844']
import plotlyGUIplots as pgp
fig = pgp.plotlyqDivPlot(hfs,nombres,logPlot=True)



#===============================================================================
#                       Plot maxT, maxHF, etc.
#===============================================================================
#read data file
import plotlyGUIplots as pgp
import pandas as pd
nombres = ['T002','T101A']
names = ['Tprobe_T002','Tprobe_T101A']
root = '/home/tom/HEAT/data/d3d_000002/openFoam/heatFoam/'
data = []
for name in nombres:
    outfile = root+name+'/postProcessing/fieldMinMax1/0/minMaxTnoTab.dat'
    tmp = pd.read_csv(outfile, header=1)
    tmp.columns = tmp.columns.str.strip()
    tmp = tmp.sort_values('field')
    tmp['field'] = tmp['field'].str.strip()
    data.append(tmp)




###FINISH THIS SECTION FOR PLOTTING!  there is a bug
if type(data) is not list:
    data = [data]
if type(names) is not list:
    names = [names]
#if we want multiple fields
fields = ['T']
units = [' [K]']
y2 = [False]
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

fig.update_yaxes(title_text="<b>Maximum PFC Temperature [K]</b>", secondary_y=False)

fig.show()

#===============================================================================
#                       Plot maxT for different sweeps / OF results
#===============================================================================
#read data file
import plotlyGUIplots as pgp
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
root = '/home/tom/HEAT/data/d3d_000002/openFoam/heatFoam/'
file1 = root+'T101A/postProcessing/fieldMinMax1/0/minMaxTnoTab.dat'
file2 = root+'T002/postProcessing/fieldMinMax1/0/minMaxTnoTab.dat'
file3 = root+'T005/postProcessing/fieldMinMax1/0/minMaxTnoTab.dat'
files = [file1, file2, file3]
nombres = ['T101A', 'T002', 'T005']

dashes = ['longdash', 'dashdot', 'solid'] #dashdot
#colors = ['#0dd5db', '#ab0ddb', '#eb9215', '#1e0ddb']
fields = ['T']
data = []
for file in files:
    tmp = pd.read_csv(file, header=1)
    tmp.columns = tmp.columns.str.strip()
    tmp = tmp.sort_values('field')
    tmp['field'] = tmp['field'].str.strip()
    data.append(tmp)

fig = go.Figure()
for idx,name in enumerate(nombres):
    df = data[idx]
    for i,field in enumerate(fields):
        mask = df['field'] == field
        t = df[mask].sort_values('# Time')['# Time'].values
        varMax = df[mask].sort_values('# Time')['max'].values
        fig.add_trace(go.Scatter(x=t, y=varMax, name=name, line=dict(width=6, dash=dashes[idx],
                                 #color=colors[idx]
                                 ))
                                 )

fig.add_trace(go.Scatter(
    x=[0.05, t[-1]],
    y=[1873, 1873],
    mode="lines+markers+text",
    name="Sublimation T",
    text=["Limit", "Limit"],
    textposition="top center",
    line=dict(color='firebrick', width=3, dash='dot'),
    textfont=dict(family="Arial", size=16, color="firebrick"),

))


#fig.update_layout(title='Peak Temp for Various Strike Point Frequencies',
fig.update_layout(title='Peak Temp for negD Tiles',
            xaxis_title='<b>Time [s]</b>',
            yaxis_title='<b>Temp [K]</b>',
            font=dict(
                family='Arial, sans-serif',
                size=18,
                )
                )
fig.show()


#===============================================================================
#                       generate sweep gfiles
#===============================================================================
from shutil import copyfile
import numpy as np
root = '/u/tlooby/results/WGPFC_Validation/dynamicScans/memo010_case2scan4/'
sweepDir = root + '5s_10.0Hz/'
shot = 135111
newFreq = 10.0 #Hz
duration = 5.0 #s
maxOrigt = 10 #N of original gFiles per period

origFiles = ['g135111.00001',
             'g135111.00002',
             'g135111.00003',
             'g135111.00004',
             'g135111.00005',
             'g135111.00006',
             'g135111.00007',
             'g135111.00008',
             'g135111.00009',
             'g135111.00010',]

#origFiles = ['g116313.00001',
#             'g116313.00002',
#             'g116313.00003',
#             'g116313.00004',
#             'g116313.00005',
#             'g116313.00006',
#             'g116313.00007',
#             'g116313.00008']


Norig = len(origFiles)
N = duration * newFreq * Norig
ts = np.linspace(0,duration,N+1)*1000 #[ms]

origCount = 1
for t in ts:
    infile = root + origFiles[origCount-1]
    outfile = sweepDir + 'g{:06d}.{:05d}'.format(shot,int(t))
    copyfile(infile, outfile)
    if origCount == maxOrigt:
        origCount = 1
    else:
        origCount += 1

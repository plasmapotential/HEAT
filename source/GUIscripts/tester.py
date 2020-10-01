import numpy as np
import plotly.express as px
import plotly.graph_objects as go
#file="T"
#data = np.genfromtxt(file,comments="#", autostrip=True)
#fig = px.line(x=data[:,0],y=data[:,1])

root='/u/tlooby/results/WGPFC_Validation/dynamicScans/memo010_case2scan4/5s_0.1Hz/HEAToutput/LwrOuterDiv/openFoam/heatFoam/'
file1 = root+'E-ED1408-573/postProcessing/probes/0/T'
file2 = root+'SOLID844/postProcessing/probes/0/T'

data1 = np.genfromtxt(file1,comments="#", autostrip=True)
data2 = np.genfromtxt(file2,comments="#", autostrip=True)
fig = go.Figure()
fig.add_trace(go.Scatter(x=data1[:,0], y=data1[:,1], name="OBD", line=dict(color='royalblue', width=6, dash='dashdot')))
fig.add_trace(go.Scatter(x=data2[:,0],y=data2[:,1],name="IBDH", line=dict(color='magenta', width=6)))
fig.add_trace(go.Scatter(
    x=[0, 6],
    y=[1873, 1873],
    mode="lines+markers+text",
    name="SGLR6510 Sublimation T",
    text=["Limit", "Limit"],
    textposition="top center",
    line=dict(color='firebrick', width=3, dash='dot'),
    textfont=dict(family="Arial", size=16, color="firebrick"),

))

fig.update_layout(
title="Temperature Probe Time Evolution: Case 1.21",
xaxis_title="Time [s]",
yaxis_title="Temperature [K]",
font=dict(
family="Arial",
size=24,
color="Black"
),
)
fig.show()



import pandas as pd
data1 = pd.read_csv('qDiv_E-ED1408-284.csv')
data2 = pd.read_csv('qDiv_E-ED1408-573.csv')
data3 = pd.read_csv('qDiv_SOLID844.csv')
hfs = [data1['HeatFlux'].values, data2['HeatFlux'].values, data3['HeatFlux'].values]
nombres = ['E-ED1408-284','E-ED1408-573','SOLID844']
import plotlyGUIplots as pgp
fig = pgp.plotlyqDivPlot(hfs,nombres,logPlot=True)



#plot maxT
#read data file
nombres = ['E-ED1408-284','E-ED1408-573','SOLID844']
names = ['OBD1','OBD2','SOLID844']
root = '/u/tlooby/results/WGPFC_Validation/dynamicScans/memo010_case2scan4/5s_0.1Hz/HEAToutput/LwrOuterDiv/openFoam/heatFoam/'
dataSet = []
for name in nombres:
    outfile = root+name+'/postProcessing/fieldMinMax1/0/minMaxTnoTab.dat'
    data = pd.read_csv(outfile, header=1)
    data.columns = data.columns.str.strip()
    data = data.sort_values('field')
    data['field'] = data['field'].str.strip()
    dataSet.append(data)
import plotlyGUIplots as pgp
fig = pgp.plotlyOpenFOAMplot(dataSet,names)
fig.show()

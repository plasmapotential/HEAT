import openFOAMclass
import GUIscripts.plotlyGUIplots as pg
OF = openFOAMclass.OpenFOAM()
file = './GUIscripts/fieldMinMax.dat'
data = OF.getMinMaxData(file)
fig = pg.plotlyOpenFOAMplot(data, ['SOLID844'])
fig.show()




import pandas as pd
df = pd.read_csv('qDiv_allLwrOuter.csv')
label=['OBDH']
mask = df['HeatFlux'].values > 0
hf = [df['HeatFlux'].values[mask]]

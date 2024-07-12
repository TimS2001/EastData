from Read import GetMatlabData, Get2DMatlabData
from Read import LaCl3_ID, BC501_ID, CLYC_ID, CLYC_OUTSIDE_ID
import matplotlib.pyplot as plt

from MyPlot import plotAndFind, plotWithout

str = '143550'
Hist = GetMatlabData(str, 'Energy', BC501_ID)

#plt.plot(x,y)

#plt.xlim(0,300)

#plotAndFind(Hist)
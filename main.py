import numpy as np
from scipy.signal import savgol_filter as sgf
from scipy.optimize import curve_fit as cf
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.stats import linregress
from matplotlib.colors import LogNorm


def plotMatlabData(filename, channel=1):
    data = loadmat(filename,struct_as_record=False,simplify_cells=True)
    x = data['hgS_070000']['children'][channel]['children'][0]['properties']['XData']
    y = data['hgS_070000']['children'][channel]['children'][0]['properties']['YData']
    return x, y


x,y = plotMatlabData('11-07-2024-lacl3\\143604-CH1-19088\\143604_Energy.fig',channel=8)
plt.plot(x,y)

plt.xlim(0,300)

plt.show()

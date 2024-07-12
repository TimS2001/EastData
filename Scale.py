from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
#

from MyPlot import plotAndFind, plotWithout

def GetData(Energy):
    fileName = 'Read_scale_data\\exp_data\\BC501A_N' + np.str_(Energy) + '.mat'
    data = loadmat(fileName,struct_as_record=False,simplify_cells=True)

    data = data['draw_Recoil_Proton']
    N = len(data)
    Hist = [[],[]]
    for i in range(0, N):
        Hist[0].append(data[i][0])
        Hist[1].append(data[i][1])
    return np.array(Hist)

Energy = 5510
Hist = GetData(Energy)

plotWithout(Hist, np.str_(Energy), koef = 32)
#E = plotAndFind(Hist, np.str_(Energy), [720, 1000], [0, 3e3], [0, 10.5e2], modSave = 0, koef = 16)

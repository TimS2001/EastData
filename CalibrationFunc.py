import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
from Analysis import RealDeriv, find 
from scipy.io import loadmat


def linCon(data, par = 2):
    win = par
    filt = np.ones(win)/win
    mov = win//2
    res = np.convolve(data[1], filt, mode='same')
    Hist1 = [(data[0][mov:-mov]),(res[mov:-mov])]
    return Hist1  


#Get data from .mat files

def GetData(Energy, Detector = 'BC501A'):
    fileName = 'Read_scale_data\\exp_data\\' + Detector + '_N' + np.str_(Energy) + '.mat'
    data = loadmat(fileName,struct_as_record=False,simplify_cells=True)

    data = data['draw_Recoil_Proton']
    N = len(data)
    Hist = [[],[]]
    for i in range(0, N):
        Hist[0].append(data[i][0])
        Hist[1].append(data[i][1])
    return np.array(Hist)


def plotWithout(Hist, str1, koef = 2):
    Hist = linCon(Hist, koef)
    #DHist = linCon(RealDeriv(Hist), 3)
    DHist = RealDeriv(Hist)



    fig, ax = plt.subplots(figsize=(5, 5), layout='constrained')
    k = np.max(Hist[1]) / np.max(DHist[1])

    plt.plot(Hist[0], Hist[1], '.', color ='skyblue', label = 'данные')
    plt.plot(DHist[0], DHist[1] * k, color = 'navy', label = 'производная')



    ax.set_ylim(0, np.max(Hist[1]))
    ax.set_xlim(0, np.max(Hist[0]))

    ax.grid(which='major')
    ax.set_ylabel('Отсчеты', fontsize = 14)
    ax.set_xlabel('световыход усл.ед.', fontsize=14)
    ax.set_title('E = ' + str1 + ' МэВ', fontsize=14)
    matplotlib.rc('font', size=14)
    plt.legend(loc='upper right')

    plt.show()
 


#levels - where is a peak, limits - sizes of graph
def plotAndFind(Hist, str1, levels, limitsY, limitsX, modSave = 0, koef = 2, Detector = 'BC501A'): 
    Hist = linCon(Hist, koef)
    #Hist = np.array(Hist)
    DHist = RealDeriv(Hist)



    fig, ax = plt.subplots(figsize=(5, 5), layout='constrained')
    


    #max_ = np.max(Hist[0])


    #'''
    aprox, E, S, err = find(DHist, levels, 0.9)
    #print(E, ' ', S)
    
    R = 235 * S / E
    #'''


    DHist = linCon(DHist, 6)

    k = 0.5 * limitsY[1] / np.max(aprox[1])
    plt.plot(Hist[0], Hist[1], '.', color ='skyblue', label = 'data')
    plt.plot(DHist[0], DHist[1] * k, color = 'navy', label = 'restored')
    plt.plot(aprox[0], aprox[1] * k, label = 'aproximation'+ '\nL = ' + "{:.0f}".format(round(E))  + ' ch\n' + r'$\sigma$' + '/L = ' + "{:.1f}".format(R / 2.35) + '%', color = 'red')
    

    ax.set_ylim(limitsY[0], 1.2 * limitsY[1])
    ax.set_xlim(limitsX[0], limitsX[1])

    ax.grid(which='major')
    ax.set_ylabel('Amount', fontsize = 14)
    ax.set_xlabel('channel', fontsize=14)
    ax.set_title('E = ' + str1 + ' keV', fontsize=14)
    matplotlib.rc('font', size=14)
    plt.legend(loc='upper right')

    plt.show()
    if(modSave == 1):
        fig.savefig('CalibrationPictures/' + Detector + '_' + str1 + 'keV' + '.png', bbox_inches='tight', pad_inches=0, dpi=600)
    
    return E
        



from scipy.optimize import minimize

def F(E, v):
    L = v[0] * E - v[1] * (1.0 - math.exp(-v[2] * math.pow(E, v[3])))
    return L

def Error(Energy, Light):
    def err(v):
        df = 0
        for i in range(0, len(Energy)):
            E = Energy[i]
            L = v[0] * E - v[1] * (1.0 - math.exp(-v[2] * math.pow(E, v[3])))
            df += (L - Light[i]) * (L - Light[i])
        return df
    return err


def GetKoefs(Energies, Channels, modSave = 0, Detector = 'BC501A'): #E = MeV, Ch ~ 1

    Light = Channels
    Energy = Energies

    Light = np.array(Light)

    v0 = np.array([1.0, 1.0, 1.0, 1.0])

    res = minimize(Error(Energy, Light), v0)
    a0 = res.x[0]
    a1 = res.x[1]
    a2 = res.x[2]
    a3 = res.x[3]

    v = [a0, a1, a2, a3]

    X1 = []
    Y1 = []
    dx = 0.5
    x = 1.
    
    while(x < 8):
        y1 = F(x, v)
        Y1.append(y1)
        X1.append(x)
        x += dx

    figure = plt.figure(figsize=(5, 5))
    ax = figure.add_subplot()

    ax.plot(Energy, Light, 'o', linewidth=2, label = 'cal. points')

    str = 'data'
    ax.plot(X1, Y1, linewidth=2.0, label = 'approx')

    ax.set_ylabel('Channel', fontsize = 12)
    ax.set_xlabel('Energy, keV', fontsize=12)
    plt.legend(loc='upper right')
    ax.grid(which='major')

    plt.show()
    if(modSave == 1):
        figure.savefig('CalibrationPictures/' + Detector + '_Calibration' + '.png', bbox_inches='tight', pad_inches=0, dpi=600)

    #v = [3.285, 103.4, 0.03017, 1.00428]
    v[0] = round(v[0] * 1e3) / 1e3
    v[1] = round(v[1] * 1e1) / 1e1
    v[2] = round(v[2] * 1e5) / 1e5
    v[3] = round(v[3] * 1e5) / 1e5

    return v

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from Analysis import RealDeriv, find 

def linCon(data, par = 2):
    win = par
    filt = np.ones(win)/win
    mov = win//2
    res = np.convolve(data[1], filt, mode='same')
    Hist1 = [(data[0][mov:-mov]),(res[mov:-mov])]
    return Hist1  


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
 



def plotAndFind(Hist, str1, levels, limitsY, limitsX, modSave = 0, koef = 2):#levels - where is a peak, limits - sizes of graph
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
        fig.savefig('CalibrationPictures/' + 'Resp_' + str1 + 'keV' + '.png', bbox_inches='tight', pad_inches=0, dpi=600)
    
    return E
        



import numpy as np
import math
from scipy.optimize import minimize
from scipy.integrate import quad

def Gauss(Ex, sigma):
    k = 1. / (sigma * math.sqrt(2*math.pi))
    s = -0.5 / (sigma * sigma)
    def f(x):
        return k * math.exp(s * (x - Ex)*(x - Ex))
    return f


from Constants import BIN_E, MAX_E, MIN_CH, EPS, ENERGY_ID, AMOUNT_ID, LEN, KOEF, EfHist


#Функция возращающая поправку на эффективность детектора для энергиия в МэВ
def GetEfEnergy(X, h):
    N = len(EfHist['L'])
    i = 1
    while((EfHist['L'][i] <= X)and(i < N - 1)):
        i += 1
    k = (EfHist['S'][i] - EfHist['S'][i - 1])/(EfHist['L'][i] - EfHist['L'][i - 1])
    b = EfHist['S'][i] - k * EfHist['L'][i]
    sig_C = k * X + b
    
    sig_H = 0.
    if(X <= 3.8):
        sig_H = 53.86 * math.pow(X, -0.366)
    elif(X <= 4.0):
        sig_H = -0.00986 * math.pow(X, 0.631) + 4.9162
    else:
        sig_H = 42.69 * math.pow(X, -0.316) - 0.587

    k = sig_H / (sig_H + sig_C) * (1. - math.exp(- (sig_H + sig_C)* h))
    return  1./ k 

###############

#Функция для центарльной производной
def deriv(Hist):
    N = len(Hist[ENERGY_ID])
    Amount = [0]
    
    for i in range(1, N - 1):
        k = 0.
        if(Hist[ENERGY_ID][i] > MIN_CH):
            dx = Hist[ENERGY_ID][i + 1] - Hist[ENERGY_ID][i - 1]
            dy = Hist[AMOUNT_ID][i + 1] - Hist[AMOUNT_ID][i - 1]
            k = 1. * dy / dx
        Amount.append(-k)
    Amount.append(0.)
    Amount = np.array(Amount)

    #Normalize by origin Hist
    K1 = np.max(Hist[AMOUNT_ID])
    K2 = np.max(Amount)
    Amount *= K1 / K2
    return np.array([Hist[ENERGY_ID],Amount])

#Функция - решение уравнения ГГ 2 рода
def RealDeriv(Hist):
    DHist = deriv(Hist)
    N = len(DHist[ENERGY_ID])
    for i in range(0, N):
        DHist[AMOUNT_ID][i] *= GetEfEnergy(DHist[ENERGY_ID][i], LEN)
        if(DHist[AMOUNT_ID][i] < 0):
            DHist[AMOUNT_ID][i] = 0
    #DHist[AMOUNT_ID] *= -1.
    return DHist


#Функиця возращающая апроксимирующую гистограмму 
#Также возращает квадрат ошибки аппроксимации
def SquareErr(E, S, M, Hist):
    N = len(Hist[ENERGY_ID])
    bin_ = Hist[ENERGY_ID][1] - Hist[ENERGY_ID][0]
    AprAmount = []

    for m in range(0, N):
        x0 = Hist[ENERGY_ID][0] + (m - 0.5) * bin_
        x1 = Hist[ENERGY_ID][0] + (m + 0.5) * bin_
        tmp = quad(Gauss(E, S), x0, x1)
        AprAmount.append(tmp[0])
    AprHist = np.array([Hist[ENERGY_ID], AprAmount])
    IdealSum = np.sum(Hist[AMOUNT_ID])
    Sum = np.sum(AprHist[AMOUNT_ID])
    AprHist[AMOUNT_ID] *= M * IdealSum / Sum

    err = ((AprHist[AMOUNT_ID] - Hist[AMOUNT_ID]) ** 2)

    return AprHist, np.sum(err) / IdealSum / IdealSum

#Функция для МНК 
#Возвращает функция от трех переменных (e, s, M) f
#f - квадрат разности аппроксимации в точке (e, s, M) и данных
def Error(Hist):
    N = len(Hist[ENERGY_ID])
    bin_ = Hist[ENERGY_ID][1] - Hist[ENERGY_ID][0]
    
    def f(v): #(e, s, M)
        df = 0.
        AprAmount = []

        if(v[0] <= 0)or(v[1] <= 0)or(v[2] <= 0):
            return 1.e10

        for m in range(0, N):
            x0 = Hist[ENERGY_ID][0] + (m - 0.5) * bin_
            x1 = Hist[ENERGY_ID][0] + (m + 0.5) * bin_
            tmp = quad(Gauss(v[0], v[1]), x0, x1)
            AprAmount.append(tmp[0])
        AprHist = np.array([Hist[ENERGY_ID], AprAmount])

        IdealSum = np.sum(Hist[AMOUNT_ID])
        Sum = np.sum(AprHist[AMOUNT_ID])
        k = (v[2] * IdealSum / Sum)
        AprHist[AMOUNT_ID] *= k

        
        for m in range(0, N):
            df += ((AprHist[AMOUNT_ID][m] - Hist[AMOUNT_ID][m]) ** 2)
        
        return df
    return f
        
#Функция по поиску пиков для данных 
def findPeaks(Hist):
    N = len(Hist[ENERGY_ID])

    DHist = deriv(Hist)
    
    peaks = []
    
    M = len(DHist[ENERGY_ID])
    for i in range(1, M - 1):
        if(DHist[AMOUNT_ID][i - 1] < 0)and(DHist[AMOUNT_ID][i + 1] > 0):
            peaks.append([DHist[ENERGY_ID][i], 'r'])#rise
        
        elif(DHist[AMOUNT_ID][i - 1] > 0)and(DHist[AMOUNT_ID][i + 1] < 0):
            peaks.append([DHist[ENERGY_ID][i], 'f'])#fall
    peaks.append([Hist[ENERGY_ID][N - 1],'n'])
    
    
    tmpHist = [[],[]]
    Peaks = []
    k = 0
    for i in range(0, N):
        while(peaks[k][1] != 'r')and(k < len(peaks) - 1):
            k += 1

        if(Hist[ENERGY_ID][i] > peaks[k][0]):
            Peaks.append(tmpHist)
            tmpHist = [[],[]]
            #print(peaks[k][0])
            k += 1
        else:
            tmpHist[ENERGY_ID].append(Hist[ENERGY_ID][i])
            tmpHist[AMOUNT_ID].append(Hist[AMOUNT_ID][i])
    Peaks.append(tmpHist)


    return Peaks 

#Функция возвращает мат ожидание, дисперсию и ... для пика на основе статистики
def GetKoefStat(Hist):
    N = len(Hist[ENERGY_ID])
    
    Ex0 = np.sum(Hist[AMOUNT_ID])
    if(Ex0 == 0):
        return 0, 0
    Ex1 = 0
    Ex2 = 0
    
    for i in range(0, N):
        Ex1 += Hist[1][i] * Hist[0][i]
        Ex2 += Hist[1][i] * (Hist[0][i] ** 2)

    #Normilize
    Ex1 /= Ex0
    Ex2 /= Ex0
        
    
    Dx = (Ex2 - Ex1 * Ex1)

    if(Dx <= 0):
        return 0, 0

    S = math.sqrt(N / (N - 1) * Dx) #square
   
    return Ex1, S
#####

#Функция возвращает мат ожидание, дисперсию и ошибку для пика на основе аппроксимации
def GetKoefLSP(Hist):
    M = 1.
    E, S = GetKoefStat(Hist)
    v0 = np.array([E, S, M])
    
    res = minimize(Error(Hist), v0)
    E = res.x[0]
    S = res.x[1]
    M = res.x[2]
    return E, S, M

#Функция по поиску главного пика
def Resp(Hist, levels):
    min_ = levels[0]
    max_ = levels[1]

    ################
    Peaks = findPeaks(Hist)
    numb = []
    for i in range(0, len(Peaks)):
        Peak = Peaks[i]
        E1, S1 = GetKoefStat(Peak)
        if(min_ < E1 < max_):#and(0. < S1 * S1 < 3.):
            numb.append(i)
    
    MainPeak = [[],[]]
    if(len(numb) > 0):
        Energy_min = Peaks[numb[0]][ENERGY_ID][0]
        N = len(numb) - 1
        M = len(Peaks[numb[N]][ENERGY_ID]) - 1
        Energy_max = Peaks[numb[N]][ENERGY_ID][M]
        
        i = 0
        while(Hist[ENERGY_ID][i] < Energy_min):
            i += 1
        while(Hist[ENERGY_ID][i] + BIN_E < Energy_max):
            MainPeak[ENERGY_ID].append(Hist[ENERGY_ID][i])
            MainPeak[AMOUNT_ID].append(Hist[AMOUNT_ID][i])
            i += 1
    else:
        print('NO PEAKS')
        return [[0],[0]]
    return MainPeak

#Функция производит анализ гистограммы 
#Также поиск пика, температуры, ошибки и характера реакции

from scipy.stats import t
from scipy.stats import chi2
#from scipy.stats import norm

def find(Hist, levels, eps):
    Hist = Resp(Hist, levels)

    E, S, M = GetKoefLSP(Hist)
    AprHist, errGrad = SquareErr(E, S, M, Hist)
    err = errGrad
    '''
    N = len(Hist[AMOUNT_ID])
    X = np.sum(Hist) / N
    err = 1. / math.sqrt(2 * N * X)
    
    E1, S1 = GetKoefStat(Hist)
    
    q = t.ppf((1 - eps) * 0.5, N - 1)
    e = S1 / math.sqrt(N - 1)
    E_arr = [round(X + q * e), round(X - q * e)]

    q1 = chi2.ppf(eps * 0.5, N - 1)
    q2 = chi2.ppf(1 - eps * 0.5, N - 1) 
    S_arr = [round(math.sqrt(N / q2) * S1), round(math.sqrt(N / q1) * S1)]

    #new = chi2 + norm
    q1 = round(2350 * S_arr[0] / E_arr[0])
    q2 = round(2350 * S_arr[1] / E_arr[1])
    q3 = round(2350 * S_arr[0] / E_arr[1])
    q4 = round(2350 * S_arr[1] / E_arr[0])
    q = np.array([q1/10, q2/10, q3/10, q4/10])
    R_arr = [np.min(q), np.max(q)]

    print('E = ', E_arr)
    print('S = ', S_arr)
    print('R = ', R_arr)
    '''
    return AprHist, E, S, err
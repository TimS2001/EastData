import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from scipy.optimize import minimize


#Cyclotrone data
#Energy = [1.88, 2.48, 3.13, 4.11, 5.8, 6.57, 6.8]
#Light = [2.71, 4.20, 5.99, 9.72, 16.46, 19.79, 20.75]

#EAST data from article
#Light =	[0.26, 0.338, 0.405, 0.52, 0.628, 0.719, 0.973, 1.51, 2.308, 2.936]
#Energy = [1.23, 1.5, 1.7, 2.05, 2.34, 2.61, 3.243, 4.355, 5.79, 6.65]

#EAST data 12.07.2024
Light = [0.247, 0.343, 0.674, 0.813, 1.070]
Energy = [2.040, 2.520, 3.940, 4.530, 5.510]

Light = np.array(Light)

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



v0 = np.array([1.0, 100.0, 0.01, 1.0])

res = minimize(Error(Energy, Light), v0)
a0 = res.x[0]
a1 = res.x[1]
a2 = res.x[2]
a3 = res.x[3]

v = [a0, a1, a2, a3]
print(v)
#EAST
#v = [3.285, 103.4, 0.03017, 1.00428]

#Cyclotrone
#v = [3.758, 5.145, 3.8608, 3.504]


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

ax.set_ylabel('Channel 10^3', fontsize = 12)
ax.set_xlabel('Energy, MeV', fontsize=12)
plt.legend(loc='upper right')
ax.grid(which='major')

plt.show()



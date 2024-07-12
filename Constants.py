#const
#Energy_bins
BIN_E = 0.05 #MeV
MAX_E = 18.0 #MeV
MIN_E = 0.0 #MeV
MIN_CH = 50

EPS = 'f8'
ENERGY_ID = 0
AMOUNT_ID = 1

LEN = 2. #cm

KOEF = 0.0 #const of Gauss filter

import numpy as np
EfHist = np.genfromtxt('sys_data/Carbon_cross_section.txt', dtype=[('L', '<f8'), ('S', '<f8')])

import numpy as np
import matplotlib.pyplot as plt
#

from CalibrationFunc import GetData, plotAndFind, plotWithout, GetKoefs

Detector = 'EJ301'

Energy = 5510

Hist = GetData(Energy, Detector)

plotWithout(Hist, 'data', koef = 8)

#Channel = plotAndFind(Hist, np.str_(Energy), [190, 300], [0, 4e3], [0, 350], modSave = 0, koef = 2)
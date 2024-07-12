from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
#

from CalibrationFunc import plotAndFind, plotWithout, GetKoefs, GetData


E = [2040, 2520, 3940, 4530, 5510]
Ch = [246.90068039812957, 343.21905023358335, 673.7484029249783, 813.019744529321, 1070.0095842600263]

GetKoefs(E, Ch, 0)
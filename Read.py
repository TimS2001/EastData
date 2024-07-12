from scipy.io import loadmat
import numpy as np
#time in sec

### Read Matlab Dat
'''
Channels: 
1) LaCl3
2) CLYC
3) 
4)
5) CLYC (outside)
6)
7)
8) BC501 (outside)
'''

LaCl3_ID = 1
CLYC_ID = 2
CLYC_OUTSIDE_ID = 5
BC501_ID = 8



def GetMatlabData(fileName,  mod = 'Energy', channel=1):
    fileName = 'WorkData\\' + fileName + '\\' + fileName + '_' + mod + '.fig'
    data = loadmat(fileName,struct_as_record=False,simplify_cells=True)
    x = data['hgS_070000']['children'][channel]['children'][0]['properties']['XData']
    y = data['hgS_070000']['children'][channel]['children'][0]['properties']['YData']
    return [x, y]

def Get2DMatlabData(filename, channel=1):
    mod = 'PSD'
    fileName = 'WorkData\\' + fileName + '\\' + fileName + '_' + mod + '.fig'
    data = loadmat(filename,struct_as_record=False,simplify_cells=True)
    x = data['hgS_070000']['children'][channel]['children'][0]['properties']['XData']
    y = data['hgS_070000']['children'][channel]['children'][0]['properties']['YData']
    z = data['hgS_070000']['children'][channel]['children'][0]['properties']['CData']
    z = np.pow(z, 10)
    return [x, y, z]


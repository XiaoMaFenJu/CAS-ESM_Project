import numpy as np
import netCDF4 as nc

filename1 = "t2m.nc"
data1 = nc.Dataset(filename1)
print(data1['t2m'])

filename = "D:/data/B20TRC5X_IAPDGVM-colm-2001-01.nc"
data = nc.Dataset(filename)
print(data['lat'])
import numpy as np
import netCDF4 as nc

f_t2m = nc.Dataset('./t2m.nc','w',format = 'NETCDF4')
f_prc = nc.Dataset('./prc.nc','w',format = 'NETCDF4')
f_prl = nc.Dataset('./prl.nc','w',format = 'NETCDF4')

filename = "D:/data/B20TRC5X_IAPDGVM-colm-2001-01.nc"
data = nc.Dataset(filename)
lat = data['lat'][:]
lon = data['lon'][:]
lat.long_name = "latitude"
lat.units = "degrees_north"
lon.long_name = "longitude"
lon.units = "degrees_east"

time = range(780)

f_t2m.createDimension('lat',128)
f_t2m.createDimension('lon',256)
f_t2m.createDimension('time',780)
f_t2m.createVariable('lat', np.float32, ('lat'))
f_t2m.createVariable('lon', np.float32, ('lon'))
f_t2m.createVariable('time', np.int32, ('time'))
f_t2m.variables['lat'][:] = lat
f_t2m.variables['lon'][:] = lon
f_t2m.variables['time'][:] = time
f_t2m.variables['lat'] = data.variables['lat']
f_t2m.variables['lon'] = data.variables['lon']
f_t2m.variables['time'] = data.variables['time']

# f_t2m.close()

f_prc.createDimension('lat',128)
f_prc.createDimension('lon',256)
f_prc.createDimension('time',780)
f_prc.createVariable('lat', np.float32, ('lat'))
f_prc.createVariable('lon', np.float32, ('lon'))
f_prc.createVariable('time', np.int32, ('time'))
f_prc.variables['lat'][:] = lat
f_prc.variables['lon'][:] = lon
f_prc.variables['time'][:] = time
# f_prc.close()

f_prl.createDimension('lat',128)
f_prl.createDimension('lon',256)
f_prl.createDimension('time',780)
f_prl.createVariable('lat', np.float32, ('lat'))
f_prl.createVariable('lon', np.float32, ('lon'))
f_prl.createVariable('time', np.int32, ('time'))
f_prl.variables['lat'][:] = lat
f_prl.variables['lon'][:] = lon
f_prl.variables['time'][:] = time
# f_prl.close()
data.close()

t_t = np.zeros((780,128,256))
prc_t = np.zeros((780,128,256))
prl_t = np.zeros((780,128,256))

for i in range(2014-1950+1):
    i = i + 1950
    for j in range(12):
        j = j + 1
        if j < 10:
            filename = "D:/data/B20TRC5X_IAPDGVM-colm-"+str(i)+"-0"+str(j)+".nc"
        else:
            filename = "D:/data/B20TRC5X_IAPDGVM-colm-" + str(i) + "-" + str(j) + ".nc"
        print(filename)

        # filename = "./B20TRC5X_IAPDGVM-colm-2001-01.nc"

        data = nc.Dataset(filename)
        t_name = 'tref'
        prc_name = 'prc'
        prl_name = 'prl'

        t = data[t_name][:]
        prc = data[prc_name][:]
        prl = data[prl_name][:]
        # lat = data['lat'][:]
        # lon = data['lon'][:]

        t_t[i+j-1951,:,:] = np.array(t[0,:,:])
        prc_t[i+j-1951,:,:] = np.array(prc[0,:,:])
        prl_t[i+j-1951,:,:] = np.array(prl[0,:,:])

        data.close()

# f_t2m = nc.Dataset('t2m.nc', 'a')
f_t2m.createVariable('t2m', np.float32, ('time','lat', 'lon'), fill_value=-9999.0)
f_t2m.variables['t2m'][:] = t_t[:, :, :]
f_t2m.close()

# f_prc = nc.Dataset('prc.nc', 'a')
f_prc.createVariable('prc', np.float32, ('time','lat', 'lon'), fill_value=-9999.0)
f_prc.variables['prc'][:] = prc_t[:, :, :]
f_prc.close()

# f_prl = nc.Dataset('prl.nc', 'a')
f_prl.createVariable('prl', np.float32, ('time','lat', 'lon'), fill_value=-9999.0)
f_prl.variables['prl'][:] = prl_t[:, :, :]
f_prl.close()

        #输出
        # outfile_t2m = "mq/t2m-"+str(i)+filename[23:]
        # outfile_prc = "mq/prc-" + filename[23:]
        # outfile_prl = "mq/prl-" + filename[23:]
        #
        # f_t2m = nc.Dataset(outfile_t2m, 'w', format='NETCDF4')
        # f_prc = nc.Dataset(outfile_prc, 'w', format='NETCDF4')
        # f_prl = nc.Dataset(outfile_prl, 'w', format='NETCDF4')

        # f_t2m.createDimension('lat', 128)
        # f_t2m.createDimension('lon', 256)
        # f_t2m.createVariable('lat', np.float32, ('lat'))
        # f_t2m.createVariable('lon', np.float32, ('lon'))
        # f_t2m.variables['lat'][:] = lat
        # f_t2m.variables['lon'][:] = lon
        #
        # f_prc.createDimension('lat', 128)
        # f_prc.createDimension('lon', 256)
        # f_prc.createVariable('lat', np.float32, ('lat'))
        # f_prc.createVariable('lon', np.float32, ('lon'))
        # f_prc.variables['lat'][:] = lat
        # f_prc.variables['lon'][:] = lon
        #
        # f_prl.createDimension('lat', 128)
        # f_prl.createDimension('lon', 256)
        # f_prl.createVariable('lat', np.float32, ('lat'))
        # f_prl.createVariable('lon', np.float32, ('lon'))
        # f_prl.variables['lat'][:] = lat
        # f_prl.variables['lon'][:] = lon

        # f_t2m = nc.Dataset('t2m.nc', 'a')
        # f_t2m.createVariable('t2m', np.float32, ('lat', 'lon'))
        # f_t2m.variables['t2m'][:] = t[0, :, :]
        # f_t2m.close()
        #
        # f_prc = nc.Dataset('prc.nc', 'a')
        # f_prc.createVariable('prc', np.float32, ('lat', 'lon'))
        # f_prc.variables['prc'][:] = prc[0, :, :]
        # f_prc.close()
        #
        # f_prl = nc.Dataset('prl.nc', 'a')
        # f_prl.createVariable('prl', np.float32, ('lat', 'lon'))
        # f_prl.variables['prl'][:] = prl[0, :, :]
        # f_prl.close()


















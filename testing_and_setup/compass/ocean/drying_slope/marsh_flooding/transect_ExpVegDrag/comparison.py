#!/usr/bin/env python
# coding: utf-8
# Authors: Zhendong Cao, Phillip Wolfram
# Date: May 2020

# import libs
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# plot tidal amplitude at open boundary

fname1 = 'NoVeg'
fname2 = 'VegHeight0.1'
fname3 = 'VegHeight0.5'

# mpas-o
ds1 = xr.open_dataset(fname1+'.nc')
ds2 = xr.open_dataset(fname2+'.nc')
ds3 = xr.open_dataset(fname3+'.nc')

x = np.linspace(0, 1.0, len(ds1.xtime))*48.0
ympas1 = ds1.ssh.where(ds1.tidalInputMask).mean('nCells').values
ympas2 = ds2.ssh.where(ds2.tidalInputMask).mean('nCells').values
ympas3 = ds3.ssh.where(ds3.tidalInputMask).mean('nCells').values

plt.plot(x, ympas1, marker='o', label=fname1)
plt.plot(x, ympas2, marker='^', label=fname2)
plt.plot(x, ympas3, marker='*', label=fname3)

# analytical
x = np.linspace(0,48.0,100)
y = 1.5*np.sin(x*np.pi/12.0) - 3.0
plt.plot(x, y, 'k-', lw=3, label='analytical')

plt.legend(frameon=False, loc='upper right')
plt.ylabel('Tidal amplitude (m)')
plt.xlabel('Time (hrs)')
plt.suptitle('Tidal amplitude forcing (right side) for test cases and analytical')
plt.savefig('tidalcomparison.png',dpi=300)

# plot water elevation at selected time slices
times = np.linspace(0.0,2.0,21)
for atime in times:
    plottime = int((float(atime)/0.2)*24.0)
    timestr = ds1.isel(Time=plottime).xtime.values
    ymean = ds1.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
    x = ymean.yCell.values/1000.0
    y = ymean.ssh.values
    plt.plot(x, y,'-k',label=fname1)

    ymean = ds2.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
    x = ymean.yCell.values/1000.0
    y = ymean.ssh.values
    plt.plot(x, y,'-b',label=fname2)

    ymean = ds3.isel(Time=plottime).groupby('yCell').mean(dim=xr.ALL_DIMS)
    x = ymean.yCell.values/1000.0
    y = ymean.ssh.values
    plt.plot(x, y,'-g',label=fname3)

    plt.fill_between(x=[0, 2.27, 6.25, 7.1],y1=[0,-2.94,-2.94,-4.6], y2=[-4.6,-4.6,-4.6,-4.6],color='grey',zorder=9)
    plt.xlabel('Cross shore distance (km)')
    plt.ylabel('Bathymetry (m)')
    plt.title('time='+'%.2f'%atime+'days')
    plt.legend(frameon=False, loc='upper right')
    atime='%.2f'%atime
    plt.savefig('WaterLevel@time='+atime+'days.png', dpi=300)
    plt.close()
    print('Time '+atime+ ' done!')


import xarray as xr
import numpy as np
import os

import matplotlib as mpl
#mpl.use('pdf')
    
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import datetime as dt
import calendar

#from utilities import read_mooring, get_composite_dates, MonClim, get_reanalysis
import utilities as util

def plot_panel(da, nrow, ncol, index, norm=None, cmap=None, levels=None, title=''):

    ax = plt.subplot(nrow, ncol, index, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180.,180.,45.,90.], ccrs.PlateCarree())

    da.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, add_colorbar=False)

    ax.coastlines()
    ax.gridlines()
    ax.set_title(title)
    
    return ax

def main(highest=True):

    cmap = mpl.cm.get_cmap('RdYlBu_r')
    levels = np.linspace(-20,20,11)

    transport = util.read_mooring(type='TRANSPORT', column='MeanCorr')

    slp = util.get_reanalysis('SLP')

    #def anom(x):
    #    return x - x.sel(time=slice('1981','2010')).mean(dim='time')
    
    #slpAnom = slp.groupby('time.month').apply(anom)
    slpAnom = util.monAnom(slp)
    
    theseDates = util.get_composite_dates(transport, N=14, highest=highest)

    nrow=4
    ncol=4
    ax = []

    width = 10
    height = 10
    fig = plt.figure(figsize=(width, height))
    for i, d in enumerate(theseDates):
        ax.append( plot_panel( slpAnom.sel(time=d, method='nearest'),
                               nrow, ncol, i+1, levels=levels, cmap=cmap,
                               title=d.strftime('%Y%m%d')) )
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)
    cbar_ax = fig.add_axes([0.25, 0.05, 0.6, 0.02])
    cb = mpl.colorbar.ColorbarBase(cbar_ax, norm=norm, cmap=cmap,
                                   orientation='horizontal', extend='both')
    cb.set_label('hPa')

    plt.show()

    return

if __name__ == "__main__":
    main(highest=False)
    



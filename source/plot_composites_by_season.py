#----------------------------------------------------------------------
# Plots composites and composite differences of SLP for five highest
# and lowest heat fluxes for May, June, July and August.
#
# 2018-07-17 A.P.Barrett
#----------------------------------------------------------------------

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

from utilities import read_mooring, get_composite_dates, MonClim
from constants import REANALYSIS_DIRPATH

def get_reanalysis(variable='SLP'):
    """
    Returns time series of monthly reanalysis grids for a given variable
    (SLP, U10M, V10M).
    """
    fileGlob = os.path.join(REANALYSIS_DIRPATH, variable, 
                            'MERRA2_???.tavgM_2d_slv_Nx.{:s}.??????.nc4'.format(variable))
    ds = xr.open_mfdataset(fileGlob)
    da = ds[variable].load()
    da = da * 1e-2 # Convert to hPa
    da.attrs['units'] = "hPa"
    return da

def plot_panel(da, nrow, ncol, index, norm=None, cmap=None, levels=None):

    ax = plt.subplot(nrow, ncol, index, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180.,180.,45.,90.], ccrs.PlateCarree())

    da.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, add_colorbar=False)

    ax.coastlines()
    ax.gridlines()
    ax.set_title('')
    
    return ax

def lagDate(dates, lag):
    """
    Gets dates for month lag: e.g. if Month is 6 and lag=-1, returns month 5
    """
    if lag > 12:
        print ("Doesn't handle lags greater than 12 months")
        return
    
    if lag == 0: return dates

    newdates = []
    for y, m in zip(dates.year, dates.month):
        if (m == 1) & (lag < 0):
            m = 13+lag
            y = y-1
        else:
            m = m - 1
        newdates.append(dt.datetime(y,m,1))
    return newdates
    
def make_one_composite(index, field, month=None, highest=True, N=5, lag=0):

    if month:
        index = index[index.index.month == month]

    theseDates = get_composite_dates(index, highest=highest)
    if (lag != 0): theseDates = lagDate(theseDates, lag)
    
    return field.sel(time=theseDates, method='nearest').mean(dim='time')

def month_composite(index, field, highest=True, N=5, lag=0):
    """
    Calculates composites for each month
    """
    result = xr.concat([make_one_composite(index, field, month=mo, highest=highest, N=N, lag=lag)
                        for mo in np.arange(1,13)],
                       dim='month')
    result['month'] = np.arange(1,13)
    return result

def add_construction_grid(fig, hori=None, vert=None, novert=False, nohori=False):

    if not hori: hori=np.linspace(0,1,11)
    if not vert: vert=np.linspace(0,1,11)

    print (vert)
    
    ll = []
    if not nohori:
        for h in hori:
            ll.append( mpl.lines.Line2D([0,1],[h,h], transform=fig.transFigure, figure=fig, ls='--', lw=1) )

    if not novert:
        for v in vert:
            ll.append( mpl.lines.Line2D([v,v],[0,1], transform=fig.transFigure, figure=fig, ls='--', lw=1) )

    fig.lines.extend(ll)
    
    return
    
def plot_composites(season, variable='SLP', lag=0):

    monid = {'DJF': [12,1,2],
             'MAM': [3,4,5],
             'JJA': [6,7,8],
             'SON': [9,10,11],
            }
    
    moName = calendar.month_name[:][1:]
    if lag != 0:
        moName = np.roll(moName,lag*-1)
    moName = [moName[i-1] for i in monid[season]]
    
    diff_cmap = mpl.cm.get_cmap('RdYlBu_r')
    diff_levels = np.linspace(-8,8,17)
    diff_norm = mpl.colors.BoundaryNorm(diff_levels, diff_cmap.N)

    # Index to composite on
    index = read_mooring(type='TRANSPORT', column='MeanCorr') # 

    # Variable to composite
    var = get_reanalysis(variable=variable)
    varClm = var.sel(time=slice('1981','2010')).groupby('time.month').mean(dim='time')

    hiVar = month_composite(index, var, highest=True, lag=lag)
    loVar = month_composite(index, var, highest=False, lag=lag)

    hiVarAnom = hiVar - varClm
    loVarAnom = loVar - varClm
    
    diffVar = hiVar - loVar

    width = 10.
    height = 11.
    
    fig = plt.figure(figsize=(width,height))

    nrow = 3
    ncol = 3
    ax = []

    for i, mo in enumerate(monid[season]):
        ax.append( plot_panel( hiVarAnom.sel(month=mo), nrow, ncol, 1+(i*ncol), levels=diff_levels, cmap=diff_cmap ) )
        ax.append( plot_panel( loVarAnom.sel(month=mo), nrow, ncol, 2+(i*ncol), levels=diff_levels, cmap=diff_cmap ) )
        ax.append( plot_panel( diffVar.sel(month=mo), nrow, ncol, 3+(i*ncol), levels=diff_levels, cmap=diff_cmap ) )
    
    fig.subplots_adjust(top=0.90, bottom=0.1, left=0.20, right=0.95, hspace=0.05,
                        wspace=0.1)

    #add_construction_grid(fig, vert=[0.2, 0.325, 0.45, 0.575, 0.7, 0.825, 0.95], nohori=True)
    
    norm = mpl.colors.BoundaryNorm(diff_levels, diff_cmap.N)
    cbar_ax = fig.add_axes([0.25, 0.05, 0.6, 0.02])
    cb = mpl.colorbar.ColorbarBase(cbar_ax, norm=norm, cmap=diff_cmap,
                                   orientation='horizontal', extend='both')
    cb.set_label('hPa')

    plt.figtext(0.325, 0.91, 'High', fontsize=20, rotation='horizontal', ha='center')
    plt.figtext(0.575, 0.91, 'Low', fontsize=20, rotation='horizontal', ha='center')
    plt.figtext(0.825, 0.91, 'Difference', fontsize=20, rotation='horizontal', ha='center')

    dy = (0.8/3.)
    y0 = 0.9

    plt.figtext(0.15, y0-(0.5*dy), moName[0], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.15, y0-(1.5*dy), moName[1], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.15, y0-(2.5*dy), moName[2], fontsize=20, rotation='vertical', va='center')

#    fig.set_size_inches(width, height)
#    fig.savefig('plot.pdf')
    
    #fig.savefig('june_to_sept_composites_lag{:d}.eps'.format(lag), dpi=600)
    
    plt.show()
    
    return

if __name__ == "__main__":

    #variable = 'SLP'
    #lag = 0
    #season = 'DJF'

    import argparse

    parser = argparse.ArgumentParser(description = 'Generates composites of variable for a season')
    parser.add_argument('season', type=str, help='Three letter code for season: e.g. DJF, MAM,...')
    parser.add_argument('--variable', type=str, default='SLP',
                          help='Name of variable to composite')
    parser.add_argument('--lag', '-l', type=int, default=0,
                          help='Lag in months')
    args = parser.parse_args()
    
    plot_composites(args.season, variable=args.variable, lag=args.lag)
    

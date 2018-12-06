#----------------------------------------------------------------------
# Plots composites and composite differences of SLP for five highest
# and lowest heat fluxes for May, June, July and August.
#
# 2018-07-17 A.P.Barrett
#----------------------------------------------------------------------

import xarray as xr
import numpy as np
import os
import random

#import matplotlib as mpl
#mpl.use('agg')
#mpl.use('pdf')
    
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import datetime as dt
import calendar

import utilities as util

#from utilities import read_mooring, get_composite_dates, MonClim, get_reanalysis
monid = {'DJF': [12,1,2],
         'MAM': [3,4,5],
         'JJA': [6,7,8],
         'SON': [9,10,11],}

def plot_panel(da, nrow, ncol, index, winds=None, bss=False, norm=None, cmap=None, levels=None):

    ax = plt.subplot(nrow, ncol, index, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180.,180.,45.,90.], ccrs.PlateCarree())

    da.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, add_colorbar=False)

    ax.coastlines()
    ax.gridlines()
    if winds:
        if bss:
            ax.set_title('V(BS)={:4.1f}, U(ESS)={:4.1f}, U(BSS)={:4.1f}'.format(*winds), fontsize=8)
        else:
            ax.set_title('V(BS)={:4.1f}, U(ESS)={:4.1f}, U(GOA)={:4.1f}'.format(*winds), fontsize=8)
    else:
        ax.set_title('')
    
    return ax

def make_one_composite(index, field, season=None, highest=True, N=5, lag=0):

    monid = {'DJF': [12,1,2],
             'MAM': [3,4,5],
             'JJA': [6,7,8],
             'SON': [9,10,11],}
    
    theseDates = util.get_composite_dates( index[ index.index.month.isin(monid[season]) ], highest=highest)
    if (lag != 0): theseDates = util.lagDate(theseDates, lag)
    
    return field.sel(time=theseDates, method='nearest').mean(dim='time')

def get_mcsample(field, season, nmonth=5, nsamples=100):
    """generates a monte carlo sample of months"""
    random.seed(1234)
    index = field.time[field.time.dt.month.isin(monid[season])].values
    sample = xr.concat([field.sel(time=random.choices(index, k=nsamples)).mean(dim='time')
                      for i in range(0,nsamples)], dim='sample')
#    sample = xr.concat([field.sel(time=random.sample(index, k=nsamples)).mean(dim='time')
#                      for i in range(0,nsamples)], dim='sample')
    return sample
    
def get_quantile(field, season, nmonth=5, nsamples=100, q=[0.05,0.95]):
    """Gets upper and lower quantiles of MC sampling of composites"""
    """monid = {'DJF': [12,1,2],
             'MAM': [3,4,5],
             'JJA': [6,7,8],
             'SON': [9,10,11],}
    
    index = field.time[field.time.dt.month.isin(monid[season])].values
     
    mc = xr.concat([field.sel(time=random.choices(index, k=nsamples)).mean(dim='time') for i in range(0,nsamples)], dim='sample')"""
    mcsample = get_mcsample(field, season, nmonth=nmonth, nsamples=nsamples)
    return mcsample.quantile(q, dim='sample', interpolation='linear')

def get_confidence_limits(field, nmonth=15, nsamples=100, q=[0.05,0.95]):
    """Returns confidence limits for composites"""
    season = ['DJF','MAM','JJA','SON']
    result = xr.concat([get_quantile(field, ssn, nmonth=nmonth, nsamples=nsamples, q=q) for ssn in season], dim='season')
    result['season'] = season
    result['quantile'] = ['lower','upper']
    return result

def season_composite(index, field, highest=True, N=15, lag=0):                                     
    """
    Calculates composites for each month
    """
    season = ['DJF','MAM','JJA','SON']
    result = xr.concat([make_one_composite(index, field, season=ssn, highest=highest, N=N, lag=lag)
                        for ssn in season],
                       dim='season')
    result['season'] = season
    return result

def plot_composites(variable='SLP', lag=0, bss=False):
                                     
    if bss:
        # Dictionary of regional wind speeds: V(BS), U(ESS), U(BSS)
        windd = {'DJF': {'hi': [-1.5, -0.4, -2.2], 'lo': [-5.5,  1.4, -1.3]},
                 'MAM': {'hi': [-0.8, -1.4, -1.4], 'lo': [-5.2, -0.2, -1.3]},
                 'JJA': {'hi': [-0.2, -1.4,  0.0], 'lo': [-0.7,  1.6,  0.8]},
                 'SON': {'hi': [-1.6, -1.8,  0.0], 'lo': [-4.1,  1.3,  0.0]}}
    else:
        windd = {'DJF': {'hi': [-1.5, -0.4,  1.2], 'lo': [-5.5,  1.4,  2.3]},
                 'MAM': {'hi': [-0.8, -1.4,  0.6], 'lo': [-5.2, -0.2,  2.3]},
                 'JJA': {'hi': [-0.2, -1.4,  1.7], 'lo': [-0.7,  1.6,  1.8]},
                 'SON': {'hi': [-1.6, -1.8,  2.8], 'lo': [-4.1,  1.3,  3.5]}}

    print (windd['DJF']['hi'])
    
    #### GET DATA AND GENERATE COMPOSITES ####
                                     
    # Index to composite on
    index = util.read_mooring(type='TRANSPORT', column='MeanCorr') # 

    # Variable to composite
    var = util.get_reanalysis(variable=variable)
    varAnom = util.monAnom(var) 

    # Get high and low composites
    hiVarAnom = season_composite(index, varAnom, highest=True, lag=lag)
    loVarAnom = season_composite(index, varAnom, highest=False, lag=lag)

    # Get confidence limits using Monte Carlo sampling of data
    season_confidence = get_confidence_limits(varAnom, nsamples=1000)

    # Set values within upper (0.95) and lower (0.05) bounds to missing
    hiVarAnom = hiVarAnom.where( (hiVarAnom > season_confidence.sel(quantile='upper')) |
                                 (hiVarAnom < season_confidence.sel(quantile='lower')) )
    loVarAnom = loVarAnom.where( (loVarAnom > season_confidence.sel(quantile='upper')) |
                                 (loVarAnom < season_confidence.sel(quantile='lower')) )

    # Get difference
    diffVarAnom = hiVarAnom - loVarAnom
        

    #### MAKE PLOT #####
                                     
    width = 10.
    height = 11.
    
#    cmap = mpl.cm.get_cmap('RdYlBu_r')
    cmap = mpl.cm.get_cmap('coolwarm')
    levels = np.linspace(-10,10,11)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    fig = plt.figure(figsize=(width,height))

    nrow = 4
    ncol = 3
    ax = []

    season = ['DJF','MAM','JJA','SON']
    for i, ssn in enumerate(season):
        ax.append( plot_panel( hiVarAnom.sel(season=ssn), nrow, ncol, 1+(i*ncol),
                               levels=levels, cmap=cmap, winds=windd[ssn]['hi'], bss=bss ) )
        ax.append( plot_panel( loVarAnom.sel(season=ssn), nrow, ncol, 2+(i*ncol),
                               levels=levels, cmap=cmap, winds=windd[ssn]['lo'], bss=bss ) )
        ax.append( plot_panel( diffVarAnom.sel(season=ssn), nrow, ncol, 3+(i*ncol),
                               levels=levels, cmap=cmap ) )
    
    fig.subplots_adjust(top=0.875, bottom=0.1, left=0.20, right=0.95, hspace=0.15,
                        wspace=0.1)

    #add_construction_grid(fig, vert=[0.2, 0.325, 0.45, 0.575, 0.7, 0.825, 0.95], nohori=True)
    
    cbar_ax = fig.add_axes([0.25, 0.05, 0.6, 0.02])
    cb = mpl.colorbar.ColorbarBase(cbar_ax, norm=norm, cmap=cmap,
                                   orientation='horizontal', extend='both')
    cb.set_label('hPa')

    plt.figtext(0.325, 0.91, 'High', fontsize=20, rotation='horizontal', ha='center')
    plt.figtext(0.575, 0.91, 'Low', fontsize=20, rotation='horizontal', ha='center')
    plt.figtext(0.825, 0.91, 'Difference', fontsize=20, rotation='horizontal', ha='center')

    dy = (0.8/4.)
    y0 = 0.9

    plt.figtext(0.15, y0-(0.5*dy), season[0], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.15, y0-(1.5*dy), season[1], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.15, y0-(2.5*dy), season[2], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.15, y0-(3.5*dy), season[3], fontsize=20, rotation='vertical', va='center')

#    fig.set_size_inches(width, height)
#    fig.savefig('plot.pdf')
    
    fig.savefig('june_to_sept_composites_lag{:d}_with_cf.png'.format(lag), dpi=600)
    
    #plt.show()
    
    return

if __name__ == "__main__":

    #variable = 'SLP'
    #lag = 0
    #season = 'DJF'

    import argparse

    parser = argparse.ArgumentParser(description = 'Generates composites of variable for a season')
    parser.add_argument('--variable', type=str, default='SLP',
                          help='Name of variable to composite')
    parser.add_argument('--lag', '-l', type=int, default=0,
                          help='Lag in months')
    parser.add_argument('--bss', type=bool, default=False,
                        help='Use Bering Sea Shelf winds rather than Gulf of Alaska')
    args = parser.parse_args()
    
    plot_composites(variable=args.variable, lag=args.lag, bss=args.bss)
    

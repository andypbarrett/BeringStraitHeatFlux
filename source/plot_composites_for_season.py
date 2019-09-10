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

import matplotlib as mpl
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

def make_title(windspd, windsig, fmtstr):
    """
    Generates title for composite plots showing
    significant composite windspeeds in bold
    """
    strarr = []
    for spd, sig, fmt in zip(windspd,windsig,fmtstr):
        if sig:
            strarr.append(r"$\bf{"+fmt.format(spd)+"}$")
        else:
            strarr.append(fmt.format(spd))
    return ', '.join(strarr)

def plot_panel(da, nrow, ncol, index, windspd=None, windsig=None, bss=False, norm=None, cmap=None, levels=None):

    ax = plt.subplot(nrow, ncol, index, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180.,180.,45.,90.], ccrs.PlateCarree())

    da.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, add_colorbar=False)

    ax.coastlines()
    ax.gridlines()
    
    if bss:
        fmtstr = ['V(BS)={:4.1f}', 'U(ESS)={:4.1f}', 'U(BSS)={:4.1f}']
    else:
        fmtstr = ['V(BS)={:4.1f}', 'U(ESS)={:4.1f}', 'U(GOA)={:4.1f}']
        
    if winds:
        ax.set_title( make_title(windspd, windsig, fmtstr), fontsize=8 ) 
        #'V(BS)={:4.1f}, U(ESS)={:4.1f}, U(BSS)={:4.1f}'.format(*winds)
    else:
        ax.set_title('')
    
    return ax

def plot_composites(variable='SLP', lag=0, bss=False,
                    infile='slp_composites_for_season.nc4',
                    outfile='slp_composites_for_season', verbose=False):
                                     
    if bss:
        # Dictionary of regional wind speeds: V(BS), U(ESS), U(BSS)
        windd = {'DJF': {'hi': [-1.5, -0.4, -0.3], 'lo': [-5.5,  1.4, -4.1],
                         'hisig': [True,True,True], 'losig': [True,True,True]},
                 'MAM': {'hi': [-0.8, -1.4,  0.1], 'lo': [-5.2, -0.2, -3.3],
                         'hisig': [True,True,True], 'losig': [True,False,True]},
                 'JJA': {'hi': [-0.2, -1.4,  0.6], 'lo': [-0.7,  1.6, -0.2],
                         'hisig': [False,True,False], 'losig': [False,True,True]},
                 'SON': {'hi': [-1.6, -1.8, -0.3], 'lo': [-4.1,  1.3, -2.5],
                         'hisig': [True,True,True], 'losig': [True,True,True]}}
    else:
        windd = {'DJF': {'hi': [-1.5, -0.4,  1.2], 'lo': [-5.5,  1.4,  2.3],
                         'hisig': [True,True,False], 'losig': [True,True,True]},
                 'MAM': {'hi': [-0.8, -1.4,  0.6], 'lo': [-5.2, -0.2,  2.3],
                         'hisig': [True,True,False], 'losig': [True,False,True]},
                 'JJA': {'hi': [-0.2, -1.4,  1.7], 'lo': [-0.7,  1.6,  1.8],
                         'hisig': [False,True,False], 'losig': [False,True,False]},
                 'SON': {'hi': [-1.6, -1.8,  2.8], 'lo': [-4.1,  1.3,  3.5],
                         'hisig': [True,True,False], 'losig': [True,True,False]}}

    # Get data
    ds = xr.open_dataset(infile)

    if verbose: print ('Making plot')
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
        ax.append( plot_panel( ds.hiVarAnom.sel(season=ssn), nrow, ncol, 1+(i*ncol),
                               levels=levels, cmap=cmap,
                               windspd=windd[ssn]['hi'], windsig=windd[ssn]['hisig'],
                               bss=bss ) )
        ax.append( plot_panel( ds.loVarAnom.sel(season=ssn), nrow, ncol, 2+(i*ncol),
                               levels=levels, cmap=cmap,
                               windspd=windd[ssn]['lo'], windsig=windd[ssn]['losig'],
                               bss=bss ) )
        ax.append( plot_panel( ds.diffVarAnom.sel(season=ssn), nrow, ncol, 3+(i*ncol),
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
    
    #fig.savefig('june_to_sept_composites_lag{:d}_with_cf_95percent.png'.format(lag), dpi=600)
    fig.savefig(outfile+'_lag{:d}_with_cf_95percent.png'.format(lag), dpi=600)
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
    parser.add_argument('--outfile', type=str, default='slp_composite',
                          help='Name of outfile')
    parser.add_argument('--bss', action='store_true',
                        help='Use Bering Sea Shelf winds rather than Gulf of Alaska')
    parser.add_argument('--verbose', '-v', action='store_true')
    
    args = parser.parse_args()
    
    plot_composites(variable=args.variable, lag=args.lag, bss=args.bss, verbose=args.verbose)
    

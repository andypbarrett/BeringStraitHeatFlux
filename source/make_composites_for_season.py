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

def make_composites(variable='SLP', lag=0, bss=False, outfile='bering_strait_slp_composites.nc4', verbose=False):
                                     
    # Index to composite on
    if verbose: print ('Getting mooring data')
    index = util.read_mooring(type='TRANSPORT', column='MeanCorr') # 

    # Variable to composite
    if verbose: print ('Getting reanalysis data') 
    var = util.get_reanalysis(variable=variable)
    varAnom = util.monAnom(var) 

    # Get high and low composites
    if verbose: print ('Generating composites')
    hiVarAnom = season_composite(index, varAnom, highest=True, lag=lag)
    loVarAnom = season_composite(index, varAnom, highest=False, lag=lag)

    # Get confidence limits using Monte Carlo sampling of data
    if verbose: print ('Estimating confidence limits')
    season_confidence = get_confidence_limits(varAnom, nsamples=1000, q=[0.025,0.975])

    # Set values within upper (0.95) and lower (0.05) bounds to missing
    hiVarAnom = hiVarAnom.where( (hiVarAnom > season_confidence.sel(quantile='upper')) |
                                 (hiVarAnom < season_confidence.sel(quantile='lower')) )
    loVarAnom = loVarAnom.where( (loVarAnom > season_confidence.sel(quantile='upper')) |
                                 (loVarAnom < season_confidence.sel(quantile='lower')) )

    # Get difference
    diffVarAnom = hiVarAnom - loVarAnom
        
    # Save data to file
    dsOut = xr.DataSet({'hiVarAnomaly': hiVarAnom,
                        'loVarAnomaly': loVarAnom,
                        'diffVarAnomaly': diffVarAnom})
    dsOut.to_netcdf(outfile)
    
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
    parser.add_argument('--outfile', type=str, default='slp_composite_for_season',
                          help='Name of outfile')
    parser.add_argument('--bss', action='store_true',
                        help='Use Bering Sea Shelf winds rather than Gulf of Alaska')
    parser.add_argument('--verbose', '-v', action='store_true')
    
    args = parser.parse_args()
    
    plot_composites(variable=args.variable, lag=args.lag, bss=args.bss, verbose=args.verbose)
    

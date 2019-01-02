# A script to test required sample size for MC sampling of SLP

#import matplotlib
#matplotlib.use('agg')

import matplotlib.pyplot as plt
import xarray as xr

import utilities as util
from plot_composites_for_season import get_quantile, get_mcsample, season_composite   

def main():

    season = 'DJF'
    lat = 70.
    lon = 0.
    
    lag = 0
    
    print ('Getting data...')
    da = util.get_reanalysis()
    test = da.sel(lat=lat, lon=lon, method='nearest')
    testAnom = util.monAnom(test)

    # Index to composite on
    index = util.read_mooring(type='TRANSPORT', column='MeanCorr') # 

    # Get high and low composites
    hiTestAnom = season_composite(index, testAnom, highest=True, lag=lag)
    loTestAnom = season_composite(index, testAnom, highest=False, lag=lag)
#    print (hiTestAnom)
    
    sample = get_mcsample(testAnom, season, nmonth=5, nsamples=1000)
    quantile = sample.quantile([0.05,0.95], interpolation='linear')
    
    fig, ax = plt.subplots(figsize=(10,6))

    ax.hist(sample, bins=20, label='Random Samples')
    ax.axvline(quantile[0], c='k')
    ax.axvline(quantile[1], c='k')
    ax.plot(hiTestAnom.sel(season=season), 5, 'o', color='red', label='High Transport Composite')
    ax.plot(loTestAnom.sel(season=season), 5, 'o', color='cyan', label='Low Transport Composite')
    
    ax.set_title('15 month {:s} composites: lat={:4.1f}N, lon={:5.1f}E'.format(season,lat,lon))

    plt.legend()
    
    plt.show()
    fig.savefig('composite_percentile_test.{:s}.{:4.1f}N_{:5.1f}E.png'.format(season,lat,lon))

if __name__ == "__main__":
    main()
    
    

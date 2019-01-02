# A script to test required sample size for MC sampling of SLP

import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import xarray as xr

import utilities as util
from plot_composites_for_season import get_quantile   

def main():

    nsamples = range(0,1501,100)[1:]

    print ('Getting data...')
    da = util.get_reanalysis()
    test = da.sel(lat=70., lon=0., method='nearest')
    
    upper_mean = []
    upper_min = []
    upper_max = []
    lower_mean = []
    lower_min = []
    lower_max = []
    
    for n in nsamples:

        print ('Getting estimate for {:d} samples'.format(n))

        cf_estimate = xr.concat([get_quantile(test, season='DJF', nsamples=n)
                                 for i in range(100)], dim='ensemble')
        upper_mean.append(cf_estimate.sel(quantile=0.95).mean(dim='ensemble'))
        upper_min.append(cf_estimate.sel(quantile=0.95).min(dim='ensemble'))
        upper_max.append(cf_estimate.sel(quantile=0.95).max(dim='ensemble'))
        lower_mean.append(cf_estimate.sel(quantile=0.05).mean(dim='ensemble'))
        lower_min.append(cf_estimate.sel(quantile=0.05).min(dim='ensemble'))
        lower_max.append(cf_estimate.sel(quantile=0.05).max(dim='ensemble'))

    fig, ax = plt.subplots(figsize=(10,6))

    ax.fill_between(nsamples, upper_min, upper_max, color='lightblue', alpha=0.5)
    ax.plot(nsamples, upper_mean, ls='-', marker='o', color='b', zorder=3, label='Upper CF')
    ax.fill_between(nsamples, lower_min, lower_max, color='darksalmon', alpha=0.5)
    ax.plot(nsamples, lower_mean, ls='-', marker='o', color='red', zorder=4, label='Lower CF')
    plt.legend()

    fig.savefig('sample_size.png')

if __name__ == "__main__":
    main()
    
    

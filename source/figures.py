import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from constants import MAP_PROJ, MAP_EXTENT


def plot_composite_panel(slp, u, v, nrow, ncol, index,
               norm=None, cmap=None, levels=None,
               title=None):
    """Generates a single map of SLP and wind anomalies for figures 04 and 05

    slp: sea level pressure, xarray dataset like (nlon, nlat)
    u: zonal wind, xarray dataset (nlon, nlat)
    v: meridional wind, xarray dataset (nlon, nlat)
    nrow: number of rows in plot grid
    ncol: number of rows in plot grid
    index: position index of subplot in nrow x ncol grid - index starts at 1
    """

    ax = plt.subplot(nrow, ncol, index, projection=ccrs.NorthPolarStereo())
    ax.set_extent([-180., 180., 45., 90.], ccrs.PlateCarree())

    slp.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels, cmap=cmap, add_colorbar=False)
    #ax.quiver(u.lon, u.lat, u.values, v.values, transform=ccrs.PlateCarree(), zorder=3)

    ax.coastlines()
    ax.gridlines()

    if title:
        ax.set_title(title)
    else:
        ax.set_title('')

    return ax


def plot_composite_figure(slp, u, v, verbose=False):
    """Generates figures 4 and 5"""

    if verbose: print('Making plot')
    width = 10.
    height = 11.

    cmap = mpl.cm.get_cmap('coolwarm')
    levels = np.linspace(-10, 10, 11)
    norm = mpl.colors.BoundaryNorm(levels, cmap.N)

    fig = plt.figure(figsize=(width, height))

    nrow = 4
    ncol = 3
    ax = []

    season = ['DJF', 'MAM', 'JJA', 'SON']
    for i, ssn in enumerate(season):
        ax.append(plot_composite_panel(slp.hiVarAnom.sel(season=ssn),
                                       u.hiVarAnom.sel(season=ssn),
                                       v.hiVarAnom.sel(season=ssn),
                                       nrow, ncol, 1 + (i * ncol),
                                       levels=levels, cmap=cmap))
        ax.append(plot_composite_panel(slp.loVarAnom.sel(season=ssn),
                                       u.loVarAnom.sel(season=ssn),
                                       v.loVarAnom.sel(season=ssn),
                                       nrow, ncol, 2 + (i * ncol),
                                       levels=levels, cmap=cmap))
        ax.append(plot_composite_panel(slp.diffVarAnom.sel(season=ssn),
                                       u.diffVarAnom.sel(season=ssn),
                                       v.diffVarAnom.sel(season=ssn),
                                       nrow, ncol, 3 + (i * ncol),
                                       levels=levels, cmap=cmap))

    fig.subplots_adjust(top=0.875, bottom=0.1, left=0.20, right=0.95, hspace=0.15,
                        wspace=0.1)

    # add_construction_grid(fig, vert=[0.2, 0.325, 0.45, 0.575, 0.7, 0.825, 0.95], nohori=True)

    cbar_ax = fig.add_axes([0.25, 0.05, 0.6, 0.02])
    cb = mpl.colorbar.ColorbarBase(cbar_ax, norm=norm, cmap=cmap,
                                   orientation='horizontal', extend='both')
    cb.set_label('hPa')

    tx0 = np.mean([ax[0].get_position().x0, ax[0].get_position().x1])
    tx1 = np.mean([ax[1].get_position().x0, ax[1].get_position().x1])
    tx2 = np.mean([ax[2].get_position().x0, ax[2].get_position().x1])

    plt.figtext(tx0, 0.89, 'High', fontsize=20, rotation='horizontal', ha='center')
    plt.figtext(tx1, 0.89, 'Low', fontsize=20, rotation='horizontal', ha='center')
    plt.figtext(tx2, 0.89, 'Difference', fontsize=20, rotation='horizontal', ha='center')

    dy = (0.8 / 4.)
    y0 = 0.9

    plt.figtext(0.17, y0 - (0.5 * dy), season[0], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.17, y0 - (1.5 * dy), season[1], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.17, y0 - (2.5 * dy), season[2], fontsize=20, rotation='vertical', va='center')
    plt.figtext(0.17, y0 - (3.5 * dy), season[3], fontsize=20, rotation='vertical', va='center')

#    fig.savefig(outfile + '_lag{:d}_with_cf_95percent.png'.format(lag), dpi=600)
    # plt.show()

    return fig

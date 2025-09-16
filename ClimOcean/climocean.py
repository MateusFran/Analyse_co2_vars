import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

class ClimOcean:
    def __init__(self, ds):
        """
        ds: xarray.Dataset
        """
        self.ds = ds

    def plot_map(self, data, lat, lon, title, cmap, cbar_label,
                 vmin=None, vmax=None, extent=None, savepath=None,
                 add_contours=True, contour_levels=20, contour_color='k',
                 contour_linewidth=0.5, contour_labels=False, contour_fmt="%.1f",
                 aspect='auto', levels=20):
        
        fig = plt.figure(figsize=(10,6))
        ax = plt.subplot(111, projection=ccrs.PlateCarree())
        if extent is not None:
            ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.set_aspect(aspect)

        ax.coastlines(resolution="110m", linewidth=1)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, facecolor="lightgrey")
        gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlocator = plt.MultipleLocator(10)
        gl.ylocator = plt.MultipleLocator(5)
        im = data.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            levels=levels,
            vmin=vmin, vmax=vmax,
            cbar_kwargs={"label": cbar_label, "shrink": 1, "aspect": 20, "pad": 0.05}
        )
        if add_contours:
            contours = ax.contour(
                lon, lat, data,
                levels=contour_levels,
                colors=contour_color,
                linewidths=contour_linewidth,
                transform=ccrs.PlateCarree()
            )
            if contour_labels:
                ax.clabel(contours, inline=True, fontsize=8, fmt=contour_fmt)
        plt.title(title, fontsize=14)
        if savepath is not None:
            plt.savefig(savepath, dpi=300, bbox_inches="tight")
        plt.show()

    def get_stats(self, var, mean_data, monthly_climatology, anomalies):
        return {
            "mean": {
                "min": float(mean_data[var].min().values),
                "max": float(mean_data[var].max().values),
            },
            "climatology": {
                "min": float(monthly_climatology[var].min().values),
                "max": float(monthly_climatology[var].max().values),
            },
            "anomaly": {
                "min": float(anomalies[var].min().values),
                "max": float(anomalies[var].max().values),
            },
        }

    def get_mean(self, var):
        return self.ds[var].mean(dim="time")

    def get_monthly_climatology(self, var):
        return self.ds[var].groupby("time.month").mean("time")


    def get_anomaly(self, var):
        monthly_climatology = self.get_monthly_climatology(var)
        return self.ds[var].groupby("time.month") - monthly_climatology

    def get_seasonal_climatology(self, var):
        """
        Calculates seasonal climatology (3-month seasons).
        Returns a dictionary with means for each season.
        Southern Hemisphere seasons:
            - Summer: Dec, Jan, Feb
            - Autumn: Mar, Apr, May
            - Winter: Jun, Jul, Aug
            - Spring: Sep, Oct, Nov
        """
        months_seasons = {
            'Summer':   [12, 1, 2],
            'Autumn':  [3, 4, 5],
            'Winter': [6, 7, 8],
            'Spring': [9, 10, 11]
        }
        result = {}
        for season, months in months_seasons.items():
            result[season] = self.ds[var].sel(time=self.ds['time.month'].isin(months)).mean(dim="time")
        return result

    @staticmethod
    def robust_sym_limits(arr, pct=98):
        a = np.asarray(arr)
        vmax = np.nanpercentile(np.abs(a), pct)
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        return -vmax, vmax

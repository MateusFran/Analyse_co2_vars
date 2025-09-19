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
                 contour_linewidth=0.2, contour_labels=False, contour_fmt="%.1f",
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

    def plot_seasonal_combined(self, var, lat, lon, cmap, cbar_label,
                              vmin=None, vmax=None, extent=None, savepath=None,
                              add_contours=False, levels=20, aspect='equal'):
        """
        Plot all 4 seasons in a single figure (2x2 grid)
        """
        seasonal_data = self.get_seasonal_climatology(var)
        
        # Season order and labels
        seasons = ['Summer', 'Autumn', 'Winter', 'Spring']
        season_labels = ['DJF', 'MAM', 'JJA', 'SON']
        
        fig = plt.figure(figsize=(14, 10))
        fig.suptitle('Chlorophyll Concentration', fontsize=16, fontweight='bold', y=0.95)
        
        for i, (season, label) in enumerate(zip(seasons, season_labels)):
            ax = plt.subplot(2, 2, i+1, projection=ccrs.PlateCarree())
            
            if extent is not None:
                ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.set_aspect(aspect)
            
            # Map features
            ax.coastlines(resolution="110m", linewidth=1)
            ax.add_feature(cfeature.BORDERS, linewidth=0.5)
            ax.add_feature(cfeature.LAND, facecolor="lightgrey")
            
            # Gridlines - show labels only on left and bottom edges
            gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlocator = plt.MultipleLocator(10)
            gl.ylocator = plt.MultipleLocator(5)
            
            # Show lat/lon labels only on appropriate edges
            if i in [0, 2]:  # Left column
                gl.left_labels = True
            else:
                gl.left_labels = False
                
            if i in [2, 3]:  # Bottom row
                gl.bottom_labels = True
            else:
                gl.bottom_labels = False
            
            # Plot data
            im = seasonal_data[season].plot(
                ax=ax,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                levels=levels,
                vmin=vmin, vmax=vmax,
                add_colorbar=False
            )
            
            # Contours if requested
            if add_contours:
                contours = ax.contour(
                    lon, lat, seasonal_data[season],
                    levels=levels,
                    colors='k',
                    linewidths=0.2,
                    transform=ccrs.PlateCarree()
                )
            
            ax.set_title(f'{label}', fontsize=14, fontweight='bold')
        
        # Add colorbar with adjusted spacing
        plt.subplots_adjust(bottom=0.15, right=0.85, wspace=0, hspace=0.15, top=0.9)
        cbar_ax = fig.add_axes([0.87, 0.2, 0.02, 0.6])
        cbar = plt.colorbar(im, cax=cbar_ax)
        cbar.set_label(cbar_label, fontsize=12)
        
        if savepath is not None:
            plt.savefig(savepath, dpi=300, bbox_inches="tight")
        plt.show()

    def get_trend(self, var):
        """
        Calculate linear trend over time for a variable.
        Returns the slope (trend per year) as xarray DataArray.
        """
        # Create time coordinate as numeric (years since start)
        data_copy = self.ds[var].copy()
        
        # Assign numeric time coordinates (0, 1, 2, ... representing months)
        time_months = np.arange(len(data_copy.time))
        data_copy = data_copy.assign_coords(time=time_months)
        
        # Calculate linear trend using polyfit
        trend = data_copy.polyfit(dim='time', deg=1, skipna=True)
        slope_per_month = trend.polyfit_coefficients.sel(degree=1)
        
        # Convert from per-month to per-year (12 months = 1 year)
        slope_per_year = slope_per_month * 12
        
        return slope_per_year

    def plot_trend(self, var, lat, lon, cmap, cbar_label, title=None,
                   vmin=None, vmax=None, extent=None, savepath=None,
                   add_contours=False, levels=20, aspect='equal'):
        """
        Plot linear trend of a variable over time.
        """
        trend_data = self.get_trend(var)
        
        fig = plt.figure(figsize=(10, 6))
        ax = plt.subplot(111, projection=ccrs.PlateCarree())
        
        if extent is not None:
            ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.set_aspect(aspect)
        
        # Map features
        ax.coastlines(resolution="110m", linewidth=1)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)
        ax.add_feature(cfeature.LAND, facecolor="lightgrey")
        
        # Gridlines
        gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlocator = plt.MultipleLocator(10)
        gl.ylocator = plt.MultipleLocator(5)
        
        # Plot trend data
        im = trend_data.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            levels=levels,
            vmin=vmin, vmax=vmax,
            cbar_kwargs={"label": cbar_label, "shrink": 1, "aspect": 20, "pad": 0.05}
        )
        
        # Contours if requested
        if add_contours:
            contours = ax.contour(
                lon, lat, trend_data,
                levels=levels,
                colors='k',
                linewidths=0.2,
                transform=ccrs.PlateCarree()
            )
        
        # Set title
        plot_title = title if title is not None else f'{var} Linear Trend'
        plt.title(plot_title, fontsize=14, fontweight='bold')
        
        if savepath is not None:
            plt.savefig(savepath, dpi=300, bbox_inches="tight")
        plt.show()

    @staticmethod
    def robust_sym_limits(arr, pct=98):
        a = np.asarray(arr)
        vmax = np.nanpercentile(np.abs(a), pct)
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        return -vmax, vmax

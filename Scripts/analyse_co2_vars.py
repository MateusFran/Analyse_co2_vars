# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 17:03:08 2025

@author: Mateus Francisco
"""

#%% Imports

import os
import pandas as pd
import xarray as xr
import numpy as np
import netCDF4 as nc
import cmocean.cm as cmo
import matplotlib.pyplot as plt
from scipy.stats import linregress
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#%% Configurações iniciais

# Selecionar variável
VARIAVEL = 'SST'

os.getcwd()
os.chdir(r'D:\Users\Mateus Francisco\OneDrive\Documentos\UFPE\Artigos\artigo_sdena_co2\Data')

savepath='D:\\Users\\Mateus Francisco\\OneDrive\\Documentos\\UFPE\\Artigos\\artigo_sdena_co2\\Results\\'

fileName = 'cmems_obs-mob_glo_bgc-car_my_irr-i_1755720033553.nc'
ds = xr.open_dataset(fileName)

time = ds['time']
lat = ds['latitude']
lon = ds['longitude']

extent = [-70, 0, -60, -30]

#%% Vars setup
# shape(time, lat, lon) - Variáveis
var_config = {
    "ph": {
        "title": "Surface pH",
        "units": "pH",
        "cmap": "plasma",
        "limits": {
            "mean":       {"vmin": None, "vmax": None},
            "climatology":{"vmin": 7.98, "vmax": 8.2},
            "anomaly":    {"vmin": None, "vmax": None},
        }
    },
    "spco2": {
        "title": "Surface pCO₂",
        "units": "µatm",
        "cmap": cmo.ice,
        "limits": {
            "mean":       {"vmin": None, "vmax": None},
            "climatology":{"vmin": 275, "vmax": 455},
            "anomaly":    {"vmin": None, "vmax": None},
        }
    },
    "talk": {
        "title": "Total Alkalinity",
        "units": "µmol/kg",
        "cmap": cmo.matter,
        "limits": {
            "mean":       {"vmin": None, "vmax": None},
            "climatology":{"vmin": 1985, "vmax": 2390},
            "anomaly":    {"vmin": None, "vmax": None},
        }
    },
    "tco2": {
        "title": "Dissolved Inorganic Carbon (DIC)",
        "units": "µmol/kg",
        "cmap": "viridis",
        "limits": {
            "mean":       {"vmin": None, "vmax": None},
            "climatology":{"vmin": 1770, "vmax": 2220},
            "anomaly":    {"vmin": None, "vmax": None},
        }
    },
    "fgco2": {
        "title": "Air–Sea CO₂ Flux",
        "units": "mmol m⁻² day⁻¹",
        "cmap": cmo.balance,
        "limits": {
            "mean":       {"vmin": -5, "vmax": 8},
            "climatology":{"vmin": -3, "vmax": 6},
            "anomaly":    {"vmin": None, "vmax": None},
        }
    },
}

month_names = {
    1:"Jan", 2:"Fev", 3:"Mar", 4:"Abr", 5:"Mai", 6:"Jun",
    7:"Jul", 8:"Ago", 9:"Set",10:"Out",11:"Nov",12:"Dez"
}

#%% Médias
means = {}

for var in var_config.keys():
    means[var] = ds[var].mean(dim="time")


#%% Climatologia mensal por variável
clim_monthly = {
    var: ds[var].groupby("time.month").mean("time")
    for var in var_config.keys()
}

#%% Anomalia = dado - climatologia mensal correspondente
anomalies = {
    var: ds[var].groupby("time.month") - clim_monthly[var]
    for var in var_config.keys()
}

#%%
stats = {}

for var in var_config.keys():
    stats[var] = {
        "mean": {
            "min": float(means[var].min().values),
            "max": float(means[var].max().values),
        },
        "climatology": {
            "min": float(clim_monthly[var].min().values),
            "max": float(clim_monthly[var].max().values),
        },
        "anomaly": {
            "min": float(anomalies[var].min().values),
            "max": float(anomalies[var].max().values),
        },
    }

#%% Limites de anomalia
def robust_sym_limits(arr, pct=98):
    # usa percentil robusto e impõe simetria ao redor de zero
    a = np.asarray(arr)
    vmax = np.nanpercentile(np.abs(a), pct)
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    return -vmax, vmax

#%% Tendências linear reestruturada

# tco2_m = tco2.copy()
# tco2_m = tco2_m.assign_coords(time=("time", np.arange(tco2_m.sizes["time"])))  # 0..467 (meses)

# fit = tco2.polyfit(dim="time", deg=1)
# slope_m = fit.polyfit_coefficients.sel(degree=1)   # var / mês
# annual_slope = slope_m * 12                              # var / ano

# r_time = pd.DataFrame(tco2_m['time'])

# print("Valor mínimo:", float(annual_slope.min().values))
# print("Valor máximo:", float(annual_slope.max().values))

# slope_min = float(annual_slope.min().values)
# slope_max = float(annual_slope.max().values)

#%% Função de plotagem

def plot_map(
    data, lat, lon, title, cmap, cbar_label,
    vmin=None, vmax=None, extent=None, savepath=None,
    add_contours=True,
    contour_levels=20,            # pode ser int (auto) ou lista/np.array de valores
    contour_color='k',
    contour_linewidth=0.5,
    contour_labels=False,
    contour_fmt="%.1f"
):
    
    fig = plt.figure(figsize=(10,6))
    ax = plt.subplot(111, projection=ccrs.PlateCarree())

    # Limites (lon_min, lon_max, lat_min, lat_max)
    if extent is not None:
        ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.set_aspect('auto')

    # Base do mapa
    ax.coastlines(resolution="110m", linewidth=1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor="lightgrey")

    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    
    gl.xlocator = plt.MultipleLocator(10)   # longitude a cada 5°
    gl.ylocator = plt.MultipleLocator(5)   # latitude a cada 5°
    
    levels = np.arange(-2, 2.5, 0.5)
    
    # Campo preenchido
    im = data.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        cmap=cmap,
        levels=20,
        vmin=vmin, vmax=vmax,
        cbar_kwargs={"label": cbar_label, "shrink": 1, "aspect": 20, "pad": 0.05}
    )

    # Isolinhas (opcional)
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

#%%

for var, cfg in var_config.items():
    plot_map(
        means[var],
        lat, lon,
        f"{cfg['title']} Mean (1985–2023)",
        cmap=cfg["cmap"],
        cbar_label=cfg["units"],
        extent=extent,
        add_contours=True,
        vmin=cfg["limits"]["mean"]["vmin"],
        vmax=cfg["limits"]["mean"]["vmax"],
        savepath=f"{savepath}_20levels-{var}_climatology.png",
    )
    
#%% Plotagem de climatologia

for var, cfg in var_config.items():
    da_clim = clim_monthly[var]  # (month, lat, lon)
    for m in range(1, 13):
        field = da_clim.sel(month=m)

        plot_map(
            field,
            lat, lon,
            f"{cfg['title']} Climatology – {month_names[m]} (1985–2023)",
            cmap=cfg["cmap"],
            cbar_label=cfg["units"],
            extent=extent,
            vmin=cfg["limits"]["climatologia"]["vmin"],
            vmax=cfg["limits"]["climatologia"]["vmax"],
            savepath=f"{savepath}/Climatology/{var}_{m:02d}_climatology.png",
        )

#%% Plotagem de anomalia

target_time = pd.to_datetime(ds["time"].values[-1])  # último tempo disponível

for var, cfg in var_config.items():
    da_anom = anomalies[var].sel(time=target_time)  # anomalia no tempo escolhido
    vmin, vmax = robust_sym_limits(da_anom.values, pct=98)

    plot_map(
        da_anom,
        lat, lon,
        f"{cfg['title']} Anomaly ({target_time.strftime('%Y-%m')})",
        cmap=cmo.balance,          # divergente, centrado em 0
        cbar_label=cfg["units"],   # mesmas unidades da variável
        extent=extent,
        vmin=vmin, vmax=vmax,      # NÃO usar limites da climatologia
        # savepath=f"{savepath}/{var}_anomaly_{target_time.strftime('%Y%m')}.png",
    )

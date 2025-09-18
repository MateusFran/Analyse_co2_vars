import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import xarray as xr
from ClimOcean.climocean import ClimOcean
import cmocean.cm as cmo

# Load dataset
file_path = r'D:\Users\Mateus Francisco\OneDrive\Documentos\UFPE\Artigos\artigo_sdena_co2\Data\chlorophyll_cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_1758131488437.nc'
ds = xr.open_dataset(file_path)

# Instantiate library
climo = ClimOcean(ds)

# Parameters
variable = 'CHL'
latitude = ds['latitude']
longitude = ds['longitude']
extent = [-70, -40, -60, -30]
colormap = cmo.speed
colorbar_label = 'mg m$^{-3}$'
vmin = 0
vmax = 5

# Output directory for seasonal means
output_dir = r'D:\Users\Mateus Francisco\OneDrive\Documentos\UFPE\Artigos\artigo_sdena_co2\Results\Seasonal_mean'
os.makedirs(output_dir, exist_ok=True)

# Plot combined seasonal figure
combined_savepath = os.path.join(output_dir, f'{variable}_seasonal_combined.png')
climo.plot_seasonal_combined(
    variable,
    latitude, longitude,
    colormap,
    colorbar_label,
    vmin=vmin,
    vmax=vmax,
    extent=extent,
    savepath=combined_savepath,
    add_contours=False,
    levels=np.arange(vmin, vmax+0.5, 0.5),
    aspect='equal'
)
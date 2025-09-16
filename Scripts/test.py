import numpy as np
import xarray as xr
from ClimOcean.climocean import ClimOcean
import cmocean.cm as cmo
import os

# Load dataset
file_path = r'D:\Users\Mateus Francisco\OneDrive\Documentos\UFPE\Artigos\artigo_sdena_co2\Data\DIATO_cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M_1757962746741.nc'
ds = xr.open_dataset(file_path)

# Instantiate library
climo = ClimOcean(ds)

# Parameters
variable = 'DIATO'
latitude = ds['latitude']
longitude = ds['longitude']
extent = [-70, -40, -60, -30]
colormap = 'turbo'
colorbar_label = 'mg m$^{-3}$'
vmin = 0
vmax = 5

# Output directory for seasonal means
output_dir = r'D:\Users\Mateus Francisco\OneDrive\Documentos\UFPE\Artigos\artigo_sdena_co2\Results\Seasonal_mean'
os.makedirs(output_dir, exist_ok=True)

# Seasonal climatology
seasonal_data = climo.get_seasonal_climatology(variable)

# Mapping months per season
months_seasons = {
    'Summer':   [12, 1, 2],
    'Autumn':  [3, 4, 5],
    'Winter': [6, 7, 8],
    'Spring': [9, 10, 11]
}
month_names = {1:'Jan', 2:'Feb', 3:'Mar', 4:'Apr', 5:'May', 6:'Jun', 7:'Jul', 8:'Aug', 9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec'}

# Plot maps for each season using month initials and save to folder
def get_month_initials(months):
    return ''.join([month_names[m][0] for m in months])

for season, data in seasonal_data.items():
    months = months_seasons[season]
    months_str = get_month_initials(months)
    title = f'Diatoms Concentration - {months_str}'
    filename = f'{variable}_mean_{months_str}.png'
    savepath = os.path.join(output_dir, filename)
    climo.plot_map(
        data,
        latitude, longitude,
        title,
        colormap,
        colorbar_label,
        vmin=vmin,
        vmax=vmax,
        extent=extent,
        add_contours=False,
        savepath=savepath,
        aspect='equal',
        levels=np.arange(vmin, vmax+0.5, 0.25)
    )

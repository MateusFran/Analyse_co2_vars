# Surface ocean carbon fields - Analyse

This repository contains analyses of surface ocean carbon fields and diatom concentrations in the South Atlantic Ocean, covering the period 1985-2023.

## Current variables analyzed:

### Carbon system variables:
- Surface partial pressure of carbon dioxide in sea water (pCO₂)
- Total alkalinity in sea water (Talk)
- Dissolved inorganic carbon in sea water (DIC)
- Sea water pH reported on total scale
- Air-sea CO₂ flux

### Biological variables:
- Diatom concentration (DIATO)

## Analysis approach:

### Spatial analysis:
- Calculation of temporal means for each variable (1985-2023)
- Monthly climatology computation
- Anomaly analysis (data minus monthly climatology)
- Statistical summary (min/max values for means, climatology, and anomalies)

### Visualization:
- Contour maps with coastlines and geographic features
- Histograms of variable distributions (linear and logarithmic scales)
- Custom color maps optimized for each variable type

## Study region:
- **Longitude**: 70°W to 0°
- **Latitude**: 60°S to 30°S
- **Focus**: South Atlantic Ocean

## Scripts:
- `analyse_co2_vars.py`: Main analysis script for carbon system variables
- `analyse_diatom_vars.ipynb`: Jupyter notebook for diatom concentration analysis

## Data sources:
- Carbon system data: CMEMS ocean biogeochemistry products
- Diatom data: CMEMS ocean color biogeochemistry products

## Key findings:
- Comprehensive spatial patterns of carbon system variables in the South Atlantic
- Seasonal variability through monthly climatologies
- Statistical characterization of variable ranges and distributions

## Author:
**Mateus Francisco**  
Universidade Federal de Pernambuco (UFPE)  
Department of Oceanography  

*Research focus: Ocean biogeochemistry, carbon cycle dynamics, and marine ecosystems in the South Atlantic Ocean.*
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 17:03:08 2025

@author: Mateus Francisco
"""

import os
import numpy as np
import pandas as pd
import netCDF4 as nc

# Selecionar variável
VARIAVEL = 'SST'

os.getcwd()
os.chdir(r'D:\Users\Mateus Francisco\OneDrive\Documentos\UFPE\Artigos\artigo_sdena_co2\Data')

fileName = 'cmems_obs-mob_glo_bgc-car_my_irr-i_1752867988386.nc'
ds = nc.Dataset(fileName)

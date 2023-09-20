#!/usr/bin/env python3
import os
import xarray as xr

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from data_parse import timestep
from plot_gen import plt_lat_lon, plt_lev_var, plt_lev_lon,plt_lev_lat,plt_lev_time,plt_lat_time



dataset = '/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/tiegcm_res5.0_decsol_smax_sech_005.nc'
timestep_array=timestep('/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/', 'sech')
print(timestep_array[119][1])
dir = '/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/'
file = 'tiegcm_res5.0_decsol_smax_sech_005.nc'
file2 = 'tiegcm_res5.0_decsol_smax_sech_004.nc'
file3 = 'tiegcm_res5.0_decsol_smax_sech_003.nc'
file4 = 'tiegcm_res5.0_decsol_smax_sech_002.nc'
file5 = 'tiegcm_res5.0_decsol_smax_sech_001.nc'
dataset = xr.open_dataset(os.path.join(dir,file))


dataset2 = xr.open_dataset(os.path.join(dir,file2))
dataset3 = xr.open_dataset(os.path.join(dir,file3))
dataset4 = xr.open_dataset(os.path.join(dir,file4))
dataset5 = xr.open_dataset(os.path.join(dir,file5))

datasets = [dataset5,dataset4,dataset3,dataset2, dataset]
print("-------------------------Plot Generation--------------------------")
fig1 = plt_lat_lon(dataset, "TN", timestep_array[119][1], -4.0)
fig2 = plt_lat_lon(dataset, "O2", timestep_array[119][1], 0.0)
fig3 = plt_lat_lon(dataset, "VN", timestep_array[119][1], 2.0)
fig4 = plt_lat_lon(dataset, "ZG", timestep_array[119][1], 4.0)
fig5 = plt_lev_var(dataset, "UN", timestep_array[119][1], -62.50, 0.0)
fig6 = plt_lev_var(dataset, "ZG", timestep_array[119][1], -62.50, -180.0)
fig7 =plt_lev_lon(dataset, "ZG", timestep_array[119][1], -62.50)
fig8 =plt_lev_lon(dataset, "UN", timestep_array[119][1], -62.50)
fig9 =plt_lev_lat(dataset, "UN", timestep_array[119][1], 0.0)
fig10 =plt_lev_lat(dataset, "UN", timestep_array[119][1], -180.0)
fig11 = plt_lev_time(datasets, "TN", -62.50, -180.0)
fig12 = plt_lat_time(datasets, "TN", -4.00, 0.0)
#fig8 = test(dataset, "TN", timestep_array[119][1], -4.25)
pdf_multi_path = "/glade/u/home/nikhilr/tiegcm_func/tiegcm2.0/benchmarks/p3postproc/test.pdf"
with PdfPages(pdf_multi_path) as pdf:
    pdf.savefig(fig1, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig2, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig3, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig4, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig5, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig6, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig7, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig8, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig9, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig10, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig11, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig12, bbox_inches='tight', pad_inches=0.5)
#print(lat_lon_lev(dataset, "TN", "2002-12-21T00:00:00", 2.25))
#print(lat_lon_ilev(dataset, "NE", "2002-12-21T00:00:00", -6.5))

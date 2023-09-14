#!/usr/bin/env python3
import os
import xarray as xr

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from data_parse import timestep
from plot_gen import plt_lat_lon, plt_lev_var, plt_lev_lon



dataset = '/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/tiegcm_res5.0_decsol_smax_sech_005.nc'
timestep_array=timestep('/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/', 'sech')
print(timestep_array[119][1])






print("-------------------------Plot Generation--------------------------")
#fig1 = plt_lat_lon(dataset, "TN", timestep_array[119][1], -4.0)
#fig2 = plt_lat_lon(dataset, "O2", timestep_array[119][1], 0.0)
#fig3 = plt_lat_lon(dataset, "VN", timestep_array[119][1], 2.0)
#fig4 = plt_lat_lon(dataset, "TN", timestep_array[119][1], 4.0)
#fig5 = plt_lev_var(dataset, "UN", timestep_array[119][1], -62.50, 0.0)
#fig6 = plt_lev_var(dataset, "UN", timestep_array[119][1], -62.50, -180.0)
fig7 =plt_lev_lon(dataset, "TN", timestep_array[119][1], -62.50)
pdf_multi_path = "/glade/u/home/nikhilr/tiegcm_func/tiegcm2.0/benchmarks/p3postproc/test.pdf"
with PdfPages(pdf_multi_path) as pdf:
    #pdf.savefig(fig1, bbox_inches='tight', pad_inches=0.5)
    #pdf.savefig(fig2, bbox_inches='tight', pad_inches=0.5)
    #pdf.savefig(fig3, bbox_inches='tight', pad_inches=0.5)
    #pdf.savefig(fig4, bbox_inches='tight', pad_inches=0.5)
    #pdf.savefig(fig5, bbox_inches='tight', pad_inches=0.5)
    #pdf.savefig(fig6, bbox_inches='tight', pad_inches=0.5)
    pdf.savefig(fig7, bbox_inches='tight', pad_inches=0.5)
#print(lat_lon_lev(dataset, "TN", "2002-12-21T00:00:00", 2.25))
#print(lat_lon_ilev(dataset, "NE", "2002-12-21T00:00:00", -6.5))

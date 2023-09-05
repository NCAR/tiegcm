#!/usr/bin/env python3
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from data_parse import lat_lon_lev, lat_lon_ilev
from plot_gen import plt_var_time_lev

dataset = '/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/tiegcm_res5.0_decsol_smax_prim_001.nc'


print("-------------------------Plot Generation--------------------------")
fig1 = plt_var_time_lev(dataset, "TE", "2002-12-21T00:00:00", -4.0)
fig2 = plt_var_time_lev(dataset, "Z", "2002-12-21T00:00:00", -2.0)
fig3 = plt_var_time_lev(dataset, "Z", "2002-12-21T00:00:00", 0.0)

pdf_multi_path = "/glade/u/home/nikhilr/tiegcm_func/tiegcm2.0/benchmarks/p3postproc/test.pdf"
with PdfPages(pdf_multi_path) as pdf:
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    pdf.savefig(fig3)
#print(lat_lon_lev(dataset, "TN", "2002-12-21T00:00:00", 2.25))
#print(lat_lon_ilev(dataset, "NE", "2002-12-21T00:00:00", -6.5))
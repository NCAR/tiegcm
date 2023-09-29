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
print(type(timestep_array[119][1]))
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

str = '2002-12-26T00:00:00'
npstr = np.datetime64(str)
print(npstr)
#files = [os.path.join(dir,file5),os.path.join(dir,file4)]
#datasetsss = xr.open_mfdataset(files, combine='by_coords')


datasets = [[dataset5,'tiegcm_res5.0_decsol_smax_sech_001.nc'],[dataset4,'tiegcm_res5.0_decsol_smax_sech_002.nc'],[dataset3,'tiegcm_res5.0_decsol_smax_sech_003.nc'],[dataset2,'tiegcm_res5.0_decsol_smax_sech_004.nc'], [dataset,'tiegcm_res5.0_decsol_smax_sech_005.nc']]
print("-------------------------Plot Generation--------------------------")
fig1 = plt_lat_lon(datasets, "TN", -4.0, None, mtime=[360,0,0])

'''
fig2 = plt_lat_lon(datasets, "O2", 0.0, time=timestep_array[119][1])
fig3 = plt_lat_lon(datasets, "VN", 2.0, timestep_array[119][1], )
fig4 = plt_lat_lon(datasets, "ZG", 4.0, timestep_array[119][1], )
fig5 = plt_lev_var(datasets, "UN", -62.50, 0.0, timestep_array[119][1], )
fig6 = plt_lev_var(datasets, "ZG", -62.50, -180.0, timestep_array[119][1], )
fig7 =plt_lev_lon(datasets, "ZG", -62.50, timestep_array[119][1], )
fig8 =plt_lev_lon(datasets, "UN", -62.50, timestep_array[119][1] ,)
fig9 =plt_lev_lat(datasets, "UN", 0.0, timestep_array[119][1],)
fig10 =plt_lev_lat(datasets, "UN", -180.0, timestep_array[119][1],)
fig11 = plt_lev_time(datasets, "TN", -42.50, -180.0)
fig12 = plt_lat_time(datasets, "TN", 0.0, -180.0)
'''

pdf_multi_path = "/glade/u/home/nikhilr/tiegcm_func/tiegcm2.0/benchmarks/p3postproc/test.pdf"
with PdfPages(pdf_multi_path) as pdf:
    pdf.savefig(fig1, bbox_inches='tight', pad_inches=0.5)
    '''
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
    '''
#print(lat_lon_lev(dataset, "TN", "2002-12-21T00:00:00", 2.25))
#print(lat_lon_ilev(dataset, "NE", "2002-12-21T00:00:00", -6.5))

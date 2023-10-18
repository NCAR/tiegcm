#!/usr/bin/env python3
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from src.data_parse import timestep
from src.plot_gen import plt_lat_lon, plt_lev_var, plt_lev_lon,plt_lev_lat,plt_lev_time,plt_lat_time



#dataset = '/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/tiegcm_res5.0_decsol_smax_sech_005.nc'
#timestep_array=timestep('/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/', 'sech')


#dir = '/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/'


dir = '/glade/work/nikhilr/work/'


datasets = []

dataset_filter = 'sech'
files = sorted(os.listdir(dir))
for file in files:
    if file.endswith('.nc') and (dataset_filter is None or dataset_filter in file):
        file_path = os.path.join(dir, file)
        datasets.append([xr.open_dataset(file_path), file])
# tiepy -dir /glade/work/nikhilr/work/ -dsf sech -plt lat_lon -var TN -zp mean -mtime 235 14 0 0
#print(datasets)
#datasets = [[dataset5,'tiegcm_res5.0_decsol_smax_sech_001.nc'],[dataset4,'tiegcm_res5.0_decsol_smax_sech_002.nc'],[dataset3,'tiegcm_res5.0_decsol_smax_sech_003.nc'],[dataset2,'tiegcm_res5.0_decsol_smax_sech_004.nc'], [dataset,'tiegcm_res5.0_decsol_smax_sech_005.nc']]
print("-------------------------Plot Generation--------------------------")

#
# tiegcm 3.0 plots
#
fig1 = plt_lat_lon(datasets, "TN", level = 'mean', mtime = [235 , 14, 0 , 0])#time=np.datetime64('2005-08-23T13:00:00.000000000'))
fig2 = plt_lat_lon(datasets, "O2", level = 0.0,mtime = [235 , 14, 0 , 0])
fig3 = plt_lat_lon(datasets, "VN", level = 2.0, mtime = [235 , 14, 0 , 0] )
fig4 = plt_lat_lon(datasets, "ZG", level = 4.0, mtime = [235 , 14, 0 , 0] )
fig5 = plt_lev_var(datasets, "UN",latitude = -62.50, localtime = 0.0, mtime = [235 , 14, 0 , 0])
fig6 = plt_lev_var(datasets, "ZG", latitude = -62.50, localtime = 'mean', mtime = [235 , 14, 0 , 0])

fig7 = plt_lev_lon(datasets, "ZG", latitude = -62.50 , mtime = [235 , 14, 0 , 0], contour_intervals = 20,variable_unit= 'km')
fig8 =plt_lev_lon(datasets, "UN", latitude = 'mean', mtime = [235 , 14, 0 , 0], contour_value= 1.00e+03)

fig9 =plt_lev_lat(datasets, "TN", longitude = 'mean', mtime = [235 , 14, 0 , 0])
fig10 =plt_lev_lat(datasets, "UN", longitude = -180.0, mtime = [235 , 14, 0 , 0])
'''
fig11 = plt_lev_time(datasets, "TN", latitude = -42.50, longitude = -180.0)

fig12 = plt_lat_time(datasets, "TN", level = 0.0, longitude = -180.0)
'''
'''

#
# tiegcm 2.0 plots
#
fig1 = plt_lat_lon(datasets, "TN", level = 'mean', mtime=[360,0,0])
fig2 = plt_lat_lon(datasets, "O2", level = 0.0,mtime = [360, 0, 0])
fig3 = plt_lat_lon(datasets, "VN", level = 2.0, mtime = [360, 0, 0] )
fig4 = plt_lat_lon(datasets, "ZG", level = 4.0, mtime = [360, 0, 0] )
fig5 = plt_lev_var(datasets, "UN",latitude = -62.50, localtime = 0.0, mtime = [360, 0, 0])
fig6 = plt_lev_var(datasets, "ZG", latitude = -62.50, localtime = 'mean', mtime = [360, 0, 0])

fig7 = plt_lev_lon(datasets, "ZG", latitude = -62.50 , mtime = [360, 0, 0], contour_intervals = 20,variable_unit= 'km')
fig8 =plt_lev_lon(datasets, "UN", latitude = 'mean', mtime = [360, 0, 0], contour_value= 1.00e+03)

fig9 =plt_lev_lat(datasets, "TN", longitude = 'mean', mtime = [360, 0, 0])
fig10 =plt_lev_lat(datasets, "UN", longitude = -180.0, mtime = [360, 0, 0])
fig11 = plt_lev_time(datasets, "TN", latitude = -42.50, longitude = -180.0)
fig12 = plt_lat_time(datasets, "TN", level = 0.0, longitude = -180.0)
'''






pdf_multi_path = "/glade/u/home/nikhilr/tiegcm_func/tiegcm2.0/benchmarks/tgcm3postproc/test.pdf"
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
    '''
    pdf.savefig(fig11, bbox_inches='tight', pad_inches=0.5)
    
    pdf.savefig(fig12, bbox_inches='tight', pad_inches=0.5)
    '''
#print(lat_lon_lev(dataset, "TN", "2002-12-21T00:00:00", 2.25))
#print(lat_lon_ilev(dataset, "NE", "2002-12-21T00:00:00", -6.5))
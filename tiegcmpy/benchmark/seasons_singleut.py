import subprocess
from tiegcmpy import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import time
bench_names=['decsol_smax', 'decsol_smin', 'junsol_smax', 'junsol_smin', 'mareqx_smax', 'mareqx_smin', 'sepeqx_smax', 'sepeqx_smin']
dirs=[]
outs=[]
dir="/glade/work/nikhilr/tiegcm3.0/benchmarks/1.25/seasons/"
dataset_filter = "sech"

for bench_name in bench_names:
    dirs.append(dir+'/'+bench_name+'/hist')
    outs.append(dir+'/'+bench_name+'/hist'+'/'+bench_name+'_singleut.pdf')


all_mtimes=[
     [[360, 0, 0, 0]],
     [[360, 0, 0, 0]],
     [[177, 0, 0, 0]],
     [[177, 0, 0, 0]],
     [[85, 0, 0, 0]],
     [[85, 0, 0, 0]],
     [[269, 0, 0, 0]],
     [[269, 0, 0, 0]],  
]
fullTime= time.time()
for x in range(len(bench_names)):
    if bench_names[x] in ["mareqx_smax","sepeqx_smax"]:
        continue
    directory = dirs[x]
    output = outs[x]
    mtime_set = all_mtimes[x]

    datasets = load_datasets(directory,dataset_filter)

    # for select UT, [lat_lon,lev_lat,lev_lon,lev_var,]
    command = ["tiegcmpy", "-rec", 
            "-dir", directory, 
            "--dataset_filter", dataset_filter, 
            "--multiple_output", output]
    inputs = []
    initTime= time.time()
    with PdfPages(output) as pdf:
        
        #
        # Plot: plt_lat_lon
        #
        # Define the list of variable_names and levels to be sent to tiepy when it requests input
        variable_names = ['TN', 'UN', 'VN','O2','O1','NO','HE','NE','TE','TI','OP','POTEN','UI_ExB','VI_ExB','WI_ExB','HMF2','NMF2', 'Z']
        gen_levels = [-7.00, -4.00, 2.00, 6.00]
        mtimes = mtime_set#[[360, 0, 0, 0]]
        plot_types = ['lat_lon']

        startTime = time.time()

        for mtime in mtimes:
            for var in variable_names:
                if var == 'HMF2' or var == 'NMF2':
                    startTimeplot = time.time()
                    plot = plt_lat_lon(datasets, var, mtime = mtime)
                    pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                    plt.close(plot)
                    print("Took %f seconds" %(time.time() - startTimeplot))
                else:
                    if var == 'POTEN':
                        levels = [2.00]
                    else:
                        levels = gen_levels
                    for lev in levels:
                        startTimeplot = time.time()
                        plot = plt_lat_lon(datasets, var,level=lev, mtime = mtime)
                        pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                        plt.close(plot)
                        print("Took %f seconds" %(time.time() - startTimeplot))
        print("plt_lat_lon took %f seconds" %(time.time() - startTime))

        #
        # Plot: plt_lev_lat
        #
        # Define the list of variable_names and levels to be sent to tiepy when it requests input
        variable_names = ['TN', 'UN', 'VN','WN','HE','TE','TI', 'Z']
        localtimes = [0.0, 12.0]
        mtimes =  mtime_set#[[360, 0, 0, 0]]
        startTime = time.time()
        for mtime in mtimes:
            for var in variable_names:
                for localtime in localtimes:
                        startTimeplot = time.time()
                        plot = plt_lev_lat(datasets, var,mtime=mtime, localtime = localtime)
                        pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                        plt.close(plot)
                        print("Took %f seconds" %(time.time() - startTimeplot))
        print("plt_lev_lat took %f seconds" %(time.time() - startTime))

        #
        # Plot: plt_lev_lon
        #
        # Define the list of variable_names and levels to be sent to tiepy when it requests input
        variable_names = ['TN', 'UN', 'VN','WN','HE','NE','TE','TI','POTEN', 'Z']
        latitudes = ['mean', -62.50, -42.50, -2.50, 37.50, 57.50]
        mtimes =  mtime_set#[[360, 0, 0, 0]]
        startTime = time.time()
        for mtime in mtimes:
            for var in variable_names:
                for latitude in latitudes:
                        startTimeplot = time.time()
                        plot = plt_lev_lon(datasets, var,mtime=mtime, latitude = latitude)
                        pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                        plt.close(plot)
                        print("Took %f seconds" %(time.time() - startTimeplot))
        print("plt_lev_lon took %f seconds" %(time.time() - startTime))
            
        #
        # Plot: plt_lev_var  #########issue with lev var
        #
        # Define the list of variable_names and levels to be sent to tiepy when it requests input
        variable_names = ['TN', 'UN', 'VN','WN','NE']
        gen_latitudes = [-62.50, -42.50 , 37.50, 57.50]
        localtimes = [0.0, 12.0]
        mtimes =  mtime_set#[[360, 0, 0, 0]]

        startTime = time.time()

        for var in variable_names:
            for mtime in mtimes:
                if var == 'NE':
                        latitudes = [-62.50, -42.50,-22.50 , 0, 17.50, 37.50, 57.50]
                else:
                    latitudes = gen_latitudes
                for latitude in latitudes:
                    for localtime in localtimes:
                        startTimeplot = time.time()
                        plot = plt_lev_var(datasets, var,mtime=mtime, latitude = latitude,localtime = localtime)
                        pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                        plt.close(plot)
                        print("Took %f seconds" %(time.time() - startTimeplot))
        print("plt_lev_var took %f seconds" %(time.time() - startTime))
    print(bench_names[x]+" took %f seconds" %(time.time() - initTime))    
print("Totally took %f seconds" %(time.time() - fullTime))    
    
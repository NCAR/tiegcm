import subprocess
from tiegcmpy import *
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import time


bench_names=['dec2006_heelis_gpi', 'dec2006_weimer_imf', 'jul2000_heelis_gpi', 'jul2000_weimer_imf', 'nov2003_heelis_gpi', 'nov2003_weimer_imf', 'whi2008_heelis_gpi', 'whi2008_weimer_imf']
dirs=[]
outs=[]
dir="/glade/work/nikhilr/tiegcm3.0/benchmarks/2.5/storms"
dataset_filter = "sech"
for bench_name in bench_names:
    dirs.append(dir+'/'+bench_name+'/hist')
    outs.append(dir+'/'+bench_name+'/hist'+'/'+bench_name+'_utlat.pdf')


fullTime= time.time()
for x in range(len(bench_names)):
    if bench_names[x] in ["nov2003_heelis_gpi"]:
        continue
    directory = dirs[x]
    output = outs[x]
    print(directory)
    datasets = load_datasets(directory,dataset_filter)
    initTime= time.time()
    with PdfPages(output) as pdf:
        #
        # Plot: plt_lat_time  
        #
        # Define the list of variable_names and levels to be sent to tiepy when it requests input
        variable_names = ['UN', 'VN','WN','HE','NE','TEC','OP','POTEN','UI_ExB','VI_ExB','WI_ExB','HMF2','NMF2', 'Z']
        gen_levels = [-4.00, 0.00, 4.00]
        gen_localtimes = [0.0, 12.0]
        startTime = time.time()
        for var in variable_names:
            if var in ['UN', 'VN','WN']:
                levels = gen_levels
                localtimes = [0.0,12.0,'mean']
            elif var in ['HMF2','NMF2','TEC']:
                if var == 'TEC':
                    localtimes = [0.0,12.0,15.0]
                else:
                    localtimes = gen_localtimes
                for localtime in localtimes:
                    startTimeplot = time.time()
                    plot = plt_lat_time(datasets,var, localtime = localtime)
                    pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                    plt.close(plot)
                    print("Took %f seconds" %(time.time() - startTimeplot))
                continue
            elif var == 'OP':
                levels = [0.0,4.0]
                localtimes = gen_localtimes
            elif var in ['POTEN','UI_ExB','VI_ExB','WI_ExB']:
                levels = [4.0]
                localtimes = gen_localtimes
            else:
                levels = gen_levels
                localtimes = gen_localtimes
            for lvl in levels:
                for localtime in localtimes:
                    startTimeplot = time.time()
                    plot = plt_lat_time(datasets, var, level = lvl,localtime = localtime)
                    pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                    plt.close(plot)
                    print("Took %f seconds" %(time.time() - startTimeplot))
        print("plt_lat_time took %f seconds" %(time.time() - startTime))
    print(bench_names[x]+" took %f seconds" %(time.time() - initTime))    
print("Totally took %f seconds" %(time.time() - fullTime))    
     
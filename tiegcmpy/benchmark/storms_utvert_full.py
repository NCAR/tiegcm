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
    outs.append(dir+'/'+bench_name+'/hist'+'/'+bench_name+'_utvert.pdf')

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
        # Plot: plt_lev_time  
        #
        # Define the list of variable_names and levels to be sent to tiepy when it requests input
        variable_names = ['UN', 'VN','WN','HE','NE','OP','POTEN','UI_ExB','VI_ExB','WI_ExB','Z']
        latitudes = [-62.50, -42.50, -2.50, 37.50, 57.50]
        gen_localtimes = [0.0, 12.0]
        startTime = time.time()
        for var in variable_names:
            localtimes = gen_localtimes
            for lat in latitudes:
                for localtime in localtimes:
                    startTimeplot = time.time()
                    plot = plt_lev_time(datasets, var, latitude = lat,localtime = localtime)
                    pdf.savefig(plot, bbox_inches='tight', pad_inches=0.5)
                    plt.close(plot)
                    print("Took %f seconds" %(time.time() - startTimeplot))
        print("plt_lat_time took %f seconds" %(time.time() - startTime))
    print(bench_names[x]+" took %f seconds" %(time.time() - initTime))    
print("Totally took %f seconds" %(time.time() - fullTime))    

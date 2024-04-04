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
    outs.append(dir+'/'+bench_name+'/hist'+'/'+bench_name+'_utvert.pdf')


fullTime= time.time()
for x in range(len(bench_names)):
    if bench_names[x] in ["mareqx_smax","sepeqx_smax"]:
        continue
    directory = dirs[x]
    output = outs[x]

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
        print("plt_lev_time took %f seconds" %(time.time() - startTime))
    print(bench_names[x]+" took %f seconds" %(time.time() - initTime))    
print("Totally took %f seconds" %(time.time() - fullTime))    
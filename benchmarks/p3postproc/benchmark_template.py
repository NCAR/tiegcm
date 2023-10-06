import subprocess


directory = '/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/'
dataset_filter = "sech"
output_filename= 'test1223'

# Define the command to run tiepy
command = ["tiepy", "-rec", 
           "-dir", directory, 
           "--dataset_filter", dataset_filter, 
           "--multiple_output", output_filename]
inputs = []
'''
#
# Plot: plt_lat_lon
#
# Define the list of variable_names and levels to be sent to tiepy when it requests input
variable_names = ['TN', 'UN', 'VN','WN','HE','NE','TE','TI','POTEN', 'Z']
levels = [-4.00, 0.00, 2.00, 4.00]
mtimes = [[360, 0, 0]]
plot_types = ['lat_lon']

for plot_type in plot_types:
    for mtime in mtimes:
        for var in variable_names:
            for lev in levels:
                    inputs.append(f"--plot_type {plot_type}  -var {var} -mtime {' '.join(map(str, mtime))} -lev {lev: .2f}")

#
# Plot: plt_lev_lat
#
# Define the list of variable_names and levels to be sent to tiepy when it requests input
variable_names = ['TN', 'UN', 'VN','WN','HE','NE','TE','TI','POTEN', 'Z']
localtimes = [0.0, 12.0]
mtimes = [[360, 0, 0]]
plot_types = ['lev_lat']
for plot_type in plot_types:
    for mtime in mtimes:
        for var in variable_names:
            for localtime in localtimes:
                    inputs.append(f"--plot_type {plot_type}  -var {var} -mtime {' '.join(map(str, mtime))} -ut {localtime: .2f}")
#
# Plot: plt_lev_lon
#
# Define the list of variable_names and levels to be sent to tiepy when it requests input
variable_names = ['TN', 'UN', 'VN','WN','HE','NE','TE','TI','POTEN', 'Z']
latitudes = [-62.50, -42.50, 37.50, 57.50]
mtimes = [[360, 0, 0]]
plot_types = ['lev_lon']
for plot_type in plot_types:
    for mtime in mtimes:
        for var in variable_names:
            for latitude in latitudes:
                    inputs.append(f"--plot_type {plot_type}  -var {var} -mtime {' '.join(map(str, mtime))} -lat {latitude: .2f}")
'''
#
# Plot: plt_lev_var  #########issue with lev var
#
# Define the list of variable_names and levels to be sent to tiepy when it requests input
variable_names = ['TN', 'UN', 'VN','WN','HE','NE','TE','TI','POTEN', 'Z']
latitudes = [-62.50, -42.50,-2.50 , 37.50, 57.50]
localtimes = [0.0, 12.0, 'mean']
mtimes = [[360, 0, 0]]
plot_types = ['lev_var']
for plot_type in plot_types:
    for mtime in mtimes:
        for var in variable_names:
            for latitude in latitudes:
                for localtime in localtimes:
                    localtime_str = f"{localtime: .2f}" if isinstance(localtime, (int, float)) else localtime
                    inputs.append(f"--plot_type {plot_type}  -var {var} -mtime {' '.join(map(str, mtime))} -lat {latitude: .2f} -ut {localtime_str}")


                
# Start the tiepy script using subprocess.Popen
process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)


for input_str in inputs:
    while True:
        line = process.stdout.readline()
        if not line:
            break 

        if "Enter command or 'exit' to terminate: " in line:

            process.stdin.write(input_str + '\n')
            process.stdin.flush()
            process.stdout.flush()
            break  
        print(line.strip())  # Print tiepy's output to console

process.stdin.write('exit\n')

process.stdin.flush()

# Wait for the tiepy script to finish and print the remaining output
for line in process.stdout:
    if "Enter command or 'exit' to terminate: " in line:
        break
    print(line.strip())

# Close the stdin and stdout pipes
process.stdin.close()
process.stdout.close()

# Wait for the tiepy script to finish
process.wait()
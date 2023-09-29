import subprocess

# Define the command to run tiepy
command = ["./tiepy", "-rec", "-dir", "/glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/", "--dataset_filter", "sech", "--multiple_output", "test"]

# Define the list of inputs to be sent to tiepy when it requests input
inputs = [
    "--plot_type lat_lon  -var TN -mtime 360 0 0 -lev -4.00",
    "--plot_type lat_lon  -var TN -mtime 360 0 0 -lev -3.00",
    "--plot_type lat_lon  -var TN -mtime 360 0 0 -lev -2.00",
    "--plot_type lat_lon  -var TN -mtime 360 0 0 -lev -1.00",
    "--plot_type lat_lon  -var TN -mtime 360 0 0 -lev  0.00",
]

# Start the tiepy script using subprocess.Popen
process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)

# Read tiepy's stdout line by line and send inputs when it's ready to receive
for input_str in inputs:
    while True:
        line = process.stdout.readline()
        if not line:
            break  # EOF
        #print(line.strip())  # Print tiepy's output to console
        if "Enter command or 'exit' to terminate: " in line:
            #print(line.strip())
            process.stdin.write(input_str + '\n')
            process.stdin.flush()
            process.stdout.flush()
            break  # Break from the inner loop to send the next input
        print(line.strip())  # Print tiepy's output to console
# Send 'exit' command to terminate tiepy script after sending all inputs

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
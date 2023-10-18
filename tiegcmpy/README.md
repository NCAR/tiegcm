# TIEGCMPY

tiegcmpy is a tool used for TIE-GCM post processing and plot generation.

## Installation

```bash
pip install tiegcmpy
```

## Requrements 
- Python >= 3.8.0
- Python pip >= 23.2.1
  - cartopy >= 0.21.1
  - numpy >= 1.24.4
  - matplotlib >= 3.7.2
  - xarray >= 2023.1.0


## Usage

tiegcmpy can be run in two modes:

1) Mode: Single Plot

   Load database for generation of a single plot to a single file
2) Mode: Multiple Plot

   Load database for generation of multiple plots to a single or multiple files

### Mode: Single Plot 

Example:
```bash
tiegcmpy --plot_type {plot_type} -dir {directory of datasets} --dataset_filter {primary or secondary files} --output_format {format of output plot} --[Other optional arguments for specific plots]
```


### Mode: Multiple Plot
#### Option 1: Initiate Interactive mode to generate multiple plots from different datasets
```bash
tiegcmpy -rec
```
Wait for the command input request. 
```bash
Entering Interactive Mode
Enter command or 'exit' to terminate:
```
Type the arguments to request a plot with the dataset and output file information
```bash
--plot_type {plot_type} -dir {directory of datasets} --dataset_filter {primary or secondary files} --output_format {format of output plot} --[Other optional arguments for specific plots]
```
Wait for input request. Type another command or 'exit'

#### Option 2: Load dataset for multiple plot generation to multiple files
```bash
tiegcmpy -rec -dir {directory of datasets} --dataset_filter {primary or secondary files} 
```
Loads datasets into memory and requests for input.
```bash
Entering Interactive Mode
Loading datasets globally.
Enter command or 'exit' to terminate: 
```
Type the arguments to request a plot with the output file information
```bash
--plot_type {plot_type} --output_format {format of output plot} --[Other optional arguments for specific plots]
```
Wait for input request. Type another command or 'exit'
#### Option 3: Load dataset for multiple plot generation to a single PDF file 
```bash
tiegcmpy -rec -dir {directory of datasets} --dataset_filter {primary or secondary files} --multiple_output {Output PDF file name}
```
Loads datasets into memory, generates the PDF file and requests for input.
```bash
Entering Interactive Mode
Loading datasets globally.
Enter command or 'exit' to terminate: 
```
Type the arguments to request a plot.
```bash
--plot_type {plot_type} --[Other optional arguments for specific plots]
```
Wait for input request. Type another command or 'exit'

## License 
"""License Information"""
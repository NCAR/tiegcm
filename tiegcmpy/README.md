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

1) Mode: API

   For use in custom python scripts or juypter notebooks

2) Mode: Command Line Interface

   For a command line interface with arguments

   a) Single Plot

   Load database for generation of a single plot to a single file

   b) Multiple Plot

   Load database for generation of multiple plots to a single or multiple files


### Mode: API
#### Importing tiegcmpy
```python
import tiegcmpy as ty
```
#### Load datasets

a) Loading a single dataset
```python
ty.load_dataset(directory, dataset_filter)
```
b) Loading multiple datasets
```python
ty.load_datasets(directory, dataset_filter)
```
#### Plot generation
The following plots can be made:

a)Latitude vs Longitude plots: plt_lat_lon

b)Pressure level vs Variable Value plots: plt_lev_var

c)Pressure level vs Longitude pltos: plt_lev_lon

d)Pressure level vs Latitude plots: plt_lev_lat

e)Pressure level vs Time plots: plt_lev_time

f)Latitude vs Time plots: plt_lat_time

Example:
```python
ty.plt_lat_lon(datasets, variable_name, level, localtime)
```

Look at functionality section for list of all plot types with required and optional arguments.

### Mode: CLI (Command Line Interface)
#### Single Plot 

Example:
```bash
tiegcmpy --plot_type {plot_type} -dir {directory of datasets} --dataset_filter {primary or secondary files} --output_format {format of output plot} --[Other optional arguments for specific plots]
```


#### Multiple Plots
##### Option 1: Initiate Interactive mode to generate multiple plots from different datasets
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

##### Option 2: Load dataset for multiple plot generation to multiple files
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
##### Option 3: Load dataset for multiple plot generation to a single PDF file 
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
## Functionality
### Latitude vs Longitude Contour Plots

This function generates a contour plot of a variable against latitude and longitude.

ty.plt_lat_lon (datasets, variable_name, time= None, mtime=None, level = None,  variable_unit = None, contour_intervals = None, contour_value = None, coastlines=False, latitude_minimum = None, latitude_maximum = None, longitude_minimum = None, longitude_maximum = None, localtime_minimum = None, localtime_maximum = None )

   **Parameters:** 
* **datasets** (xarray): The loaded dataset/s using xarray.
* **variable_name** (str): The name of the variable with latitude, longitude, ilev dimensions.
* **time** (np.datetime64, optional): The selected time e.g., '2022-01-01T12:00:00'.
* **mtime** (array, optional): The selected time as a list e.g., [1, 12, 0] for 1st day, 12 hours, 0 mins.
* **level** (float, optional): The selected lev/ilev value.
* **variable_unit** (str, optional): The desired unit of the variable.
* **contour_intervals** (int, optional): The number of contour intervals. Defaults to 20.
* **contour_value** (int, optional): The value between each contour interval.
* **coastlines** (bool, optional): Shows coastlines on the plot. Defaults to False.
* **latitude_minimum** (float, optional): Minimum latitude to slice plots. Defaults to -87.5.
* **latitude_maximum** (float, optional): Maximum latitude to slice plots. Defaults to 87.5.
* **longitude_minimum** (float, optional): Minimum longitude to slice plots. Defaults to -180.
* **longitude_maximum** (float, optional): Maximum longitude to slice plots. Defaults to 175.
* **localtime_minimum** (float, optional): Minimum localtime to slice plots.
* **localtime_maximum** (float, optional): Maximum localtime to slice plots.

**Example:** 

1. Load datasets:
```python
datasets = ty.load_datasets(directory, dataset_filter)
```
2. Set variable values:
```python
variable_name = 'TN'
value_of_mtime = [360, 0, 0, 0]
pressure_level = 4.0
unit_of_variable = 'K'
intervals = 20
```
3. Generate Latitude vs Longitude contour plot:
```python
plot = ty.plt_lat_lon (
    datasets, 
    variable_name, 
    mtime=value_of_mtime, 
    level = pressure_level,  
    variable_unit = unit_of_variable,       
    contour_intervals = intervals
    )
```

### Pressure Level vs Variable Line Plot

This function generates a line plot of a variable at a specific latitude and optional longitude, time, and local time.

ty.plt_lev_var(datasets, variable_name, latitude, time= None, mtime=None, longitude = None, localtime = None, variable_unit = None, level_minimum = None, level_maximum = None)

#### Parameters:

- **datasets** (xarray): The loaded dataset(s) using xarray.
- **variable_name** (str): The name of the variable with latitude, longitude, and ilev dimensions.
- **latitude** (float): The specific latitude value for the plot.
- **time** (np.datetime64, optional): The selected time, e.g., '2022-01-01T12:00:00'.
- **mtime** (array, optional): The selected time as a list, e.g., [1, 12, 0] for the 1st day, 12 hours, 0 minutes.
- **longitude** (float, optional): The specific longitude value for the plot.
- **localtime** (float, optional): The specific local time value for the plot.
- **variable_unit** (str, optional): The desired unit of the variable.
- **level_minimum** (float, optional): Minimum level value for the plot. Defaults to -8.
- **level_maximum** (float, optional): Maximum level value for the plot. Defaults to 8.


**Example:** 

1. Load datasets:
```python
datasets = ty.load_datasets(directory, dataset_filter)
```
2. Set variable values:
```python
variable_name = 'TN'
latitude = 30.0
time_value = '2022-01-01T12:00:00'
longitude_value = 45.0
unit_of_variable = 'K'
```
3. Generate a Level vs Variable line plot:
```python
plot = ty.plt_lev_var(
    datasets,
    variable_name,
    latitude,
    time=time_value,
    longitude=longitude_value,
    variable_unit=unit_of_variable,
    )
```


### Pressure level vs Longitude Contour Plot

This function generates a contour plot of a variable at a specific latitude against longitude, with optional time and local time.

ty.plt_lev_lon(datasets, variable_name, latitude, time= None, mtime=None, variable_unit = None, contour_intervals = 20, contour_value = None,  level_minimum = None, level_maximum = None, longitude_minimum = None, longitude_maximum = None, localtime_minimum = None, localtime_maximum = None)

#### Parameters:

- **datasets** (xarray): The loaded dataset(s) using xarray.
- **variable_name** (str): The name of the variable with latitude, longitude, and ilev dimensions.
- **latitude** (float): The specific latitude value for the plot.
- **time** (np.datetime64, optional): The selected time, e.g., '2022-01-01T12:00:00'.
- **mtime** (array, optional): The selected time as a list, e.g., [1, 12, 0] for the 1st day, 12 hours, 0 minutes.
- **variable_unit** (str, optional): The desired unit of the variable.
- **contour_intervals** (int, optional): The number of contour intervals. Defaults to 20.
- **contour_value** (int, optional): The value between each contour interval.
- **level_minimum** (float, optional): Minimum level value for the plot. Defaults to -6.75.
- **level_maximum** (float, optional): Maximum level value for the plot. Defaults to 6.75.
- **longitude_minimum** (float, optional): Minimum longitude value for the plot. Defaults to -180.
- **longitude_maximum** (float, optional): Maximum longitude value for the plot. Defaults to 175.
- **localtime_minimum** (float, optional): Minimum localtime value for the plot.
- **localtime_maximum** (float, optional): Maximum localtime value for the plot.

#### Example:
1. Load datasets:

```python
datasets = ty.load_datasets(directory, dataset_filter)
```
2. Set variable values:
```python
variable_name = 'TN'
latitude = 30.0
time_value = '2022-01-01T12:00:00'
unit_of_variable = 'K'
contour_intervals = 20
```
3. Generate a Level vs Longitude contour plot:
```python
plot = ty.plt_lev_lon(
    datasets,
    variable_name,
    latitude,
    time=time_value,
    variable_unit=unit_of_variable,
    contour_intervals=contour_intervals,
)
```
### Pressure Level vs Latitude Contour Plot

This function generates a contour plot of a variable at a specified level against latitude, with optional time, longitude, and local time.

ty.plt_lev_lat(datasets, variable_name, time= None, mtime=None, longitude = None, localtime = None, variable_unit = None, contour_intervals = 20, contour_value = None, level_minimum = None, level_maximum = None, latitude_minimum = None,latitude_maximum = None):


#### Parameters:

- **datasets** (xarray): The loaded dataset(s) using xarray.
- **variable_name** (str): The name of the variable with latitude, longitude, and ilev dimensions.
- **time** (np.datetime64, optional): The selected time, e.g., '2022-01-01T12:00:00'.
- **mtime** (array, optional): The selected time as a list, e.g., [1, 12, 0] for the 1st day, 12 hours, 0 minutes.
- **longitude** (float, optional): The specific longitude value for the plot.
- **localtime** (float, optional): The specific local time value for the plot.
- **variable_unit** (str, optional): The desired unit of the variable.
- **contour_intervals** (int, optional): The number of contour intervals. Defaults to 20.
- **contour_value** (int, optional): The value between each contour interval.
- **level_minimum** (float, optional): Minimum level value for the plot. Defaults to -6.75.
- **level_maximum** (float, optional): Maximum level value for the plot. Defaults to 6.75.
- **latitude_minimum** (float, optional): Minimum latitude value for the plot. Defaults to -87.5.
- **latitude_maximum** (float, optional): Maximum latitude value for the plot. Defaults to 87.5.

#### Example:
1. Load datasets:

```python
datasets = ty.load_datasets(directory, dataset_filter)
```
2. Set variable values:
```python
variable_name = 'TN'
time_value = '2022-01-01T12:00:00'
unit_of_variable = 'K'
contour_intervals = 20
```
3. Generate a Level vs Latitude contour plot:
```python
plot = ty.plt_lev_lat(
    datasets,
    variable_name,
    time=time_value,
    variable_unit=unit_of_variable,
    contour_intervals=contour_intervals,
)
```

### Pressure Level vs Time Contour Plot

This function generates a contour plot of a variable at a specified level against time, with optional latitude, longitude, and local time.

ty.plt_lev_lat(datasets, variable_name, time= None, mtime=None, longitude = None, localtime = None, variable_unit = None, contour_intervals = 20, contour_value = None, level_minimum = None, level_maximum = None, latitude_minimum = None,latitude_maximum = None)

#### Parameters:

- **datasets** (xarray): The loaded dataset(s) using xarray.
- **variable_name** (str): The name of the variable with latitude, longitude, time, and ilev dimensions.
- **latitude** (float): The specific latitude value for the plot.
- **longitude** (float, optional): The specific longitude value for the plot.
- **localtime** (float, optional): The specific local time value for the plot.
- **variable_unit** (str, optional): The desired unit of the variable.
- **contour_intervals** (int, optional): The number of contour intervals. Defaults to 20.
- **contour_value** (int, optional): The value between each contour interval.
- **level_minimum** (float, optional): Minimum level value for the plot. Defaults to -6.75.
- **level_maximum** (float, optional): Maximum level value for the plot. Defaults to 6.75.


#### Example:
1. Load datasets:

```python
datasets = ty.load_datasets(directory, dataset_filter)
```
2. Set variable values:
```python
variable_name = 'TN'
latitude_value = 30.0
longitude_value = 45.0
unit_of_variable = 'K'
contour_intervals = 20
```
3. Generate a Level vs Time contour plot:
```python
plot = ty.plt_lev_time(
    datasets,
    variable_name,
    latitude=latitude_value,
    longitude=longitude_value,
    variable_unit=unit_of_variable,
    contour_intervals=contour_intervals,
)
```

### Latitude vs Time Contour Plot

This function generates a contour plot of a variable at a specified latitude against time, with optional level, longitude, and local time.

#### Parameters:

- **datasets** (xarray): The loaded dataset(s) using xarray.
- **variable_name** (str): The name of the variable with latitude, longitude, time, and ilev dimensions.
- **level** (float): The specific level value for the plot.
- **longitude** (float, optional): The specific longitude value for the plot.
- **localtime** (float, optional): The specific local time value for the plot.
- **variable_unit** (str, optional): The desired unit of the variable.
- **contour_intervals** (int, optional): The number of contour intervals. Defaults to 10.
- **contour_value** (int, optional): The value between each contour interval.
- **latitude_minimum** (float, optional): Minimum latitude value for the plot. Defaults to -87.5.
- **latitude_maximum** (float, optional): Maximum latitude value for the plot. Defaults to 87.5.


#### Example:
1. Load datasets:

```python
datasets = ty.load_datasets(directory, dataset_filter)
```
2. Set variable values:
```python
variable_name = 'TN'
level_value = 4.0
longitude_value = 45.0
unit_of_variable = 'K'
contour_intervals = 10

```
3. Generate a Latitude vs Time contour plot:
```python
plot = ty.plt_lat_time(
    datasets,
    variable_name,
    level=level_value,
    longitude=longitude_value,
    variable_unit=unit_of_variable,
    contour_intervals=contour_intervals,
)
```
## File Structure

    
    ├── src                         # Directory for all tiegcmpy source files
    │   ├── tiegcmpy          
    │       ├── __init__.py         # Initialize functions for API
    │       ├── convert_units.py    # Contains unit conversion functions
    │       ├── data_parse.py       # Contains data extraction and parsing functions
    │       ├── plot_gen.py         # Contains plot generation functions
    │       ├── io.py               # Contains Input Output functions for API
    │       ├── getoptions.py       # Contains argument parser for the Command Line Interface
    │       └── main.py             # Main python file to run
    ├── README.md                   # README   
    ├── benchmark_template.py       # Template for running benchmark routines     
    ├── p3postproc.py               # Testing file    
    ├── requirements.txt            # List of required libraries     
    └── setup.py                    # PIP package builder  

## License 
"""License Information"""
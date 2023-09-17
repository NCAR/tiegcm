import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr



def timestep(directory, type):
    """
    Reads all .nc files in the given directory and generates an array with id, timestamp, and filename.
    
    Args:
    - directory (str): Path to the directory containing the .nc files.
    - type (str): Type of datasets (prim or sech)
    
    Returns:
    - timestep_array (list): List of [id, timestamp, filename] for each .nc file in the directory.
    """
    
    # Get all .nc files with the type in the name from the directory
    nc_files = [f for f in os.listdir(directory) if f.endswith('.nc') and type in f]
    
    # Extract timestamps from each file
    timestamps = []
    for nc_file in nc_files:
        ds = xr.open_dataset(os.path.join(directory, nc_file))
        for timestamp in ds['time'].values:
            timestamps.append((timestamp, nc_file))
    
    # Sort files by timestamp
    sorted_files = sorted(timestamps, key=lambda x: x[0])
    
    # Generate the timestep array
    timestep_array = [[idx, ts_file[0], ts_file[1]] for idx, ts_file in enumerate(sorted_files)]
    
    #print(timestep_array[0][1])

    return timestep_array



def lev_lon (ds, variable_name, selected_time, selected_lat):
    """
    Extract data from the dataset based on the given variable name, timestamp, and lev value.
    
    Args:
    - variable_name (str): Name of the variable to extract.
        - valid variables: ['TN', 'UN', 'VN', 'O2', 'O1', 'N2', 'NO', 'N4S', 'HE', 'TE', 'TI', 'O2P', 'OP', 'QJOULE']    
    - selected_time (str): Timestamp to filter the data.
    - selected_lev (float): Level value to filter the data.
    
    Returns:
    - var_lat_lon (array): 2D array of [variable values, lat, lon] for the given timestamp and lev.
    - variable_unit (str): Unit of the variable.
    - variable_long_name (str): Long name of the variable.
    - selected_ut (float): UT value in hours for selected_time.
    """

    # Load the dataset using xarray
    #ds = xr.open_dataset(dataset)   

    # Extract variable attributes
    variable_unit = ds[variable_name].attrs.get('units', 'N/A')
    variable_long_name = ds[variable_name].attrs.get('long_name', 'N/A')
    selected_ut = ds['ut'].sel(time=selected_time).values.item() / (1e9 * 3600)
    selected_mtime = get_mtime(ds,selected_time)
    

    # Extract the data for the given selected_time and lat  
    #var_lev_lat=get_var_lev_lat (ds, variable_name, selected_time, selected_lat)

    data = ds[variable_name].sel(time=selected_time, lat=selected_lat, method='nearest')
    unique_lons = data.lon.values

    not_all_nan_indices = ~np.isnan(data.values).all(axis=1)
    var_data = data.values[not_all_nan_indices, :]
    unique_levs = data.lev.values[not_all_nan_indices]
    
    return(var_data,unique_lons,unique_levs, selected_lat, variable_unit, variable_long_name, selected_ut, selected_mtime)



def lat_lon_lev(ds, variable_name, selected_time, selected_lev):
    """
    Extract data from the dataset based on the given variable name, timestamp, and lev value.
    
    Args:
    - variable_name (str): Name of the variable to extract.
        - valid variables: ['TN', 'UN', 'VN', 'O2', 'O1', 'N2', 'NO', 'N4S', 'HE', 'TE', 'TI', 'O2P', 'OP', 'QJOULE']    
    - selected_time (str): Timestamp to filter the data.
    - selected_lev (float): Level value to filter the data.
    
    Returns:
    - var_lat_lon (array): 2D array of [variable values, lat, lon] for the given timestamp and lev.
    - variable_unit (str): Unit of the variable.
    - variable_long_name (str): Long name of the variable.
    - selected_ut (float): UT value in hours for selected_time.
    """

    # Load the dataset using xarray
    #ds = xr.open_dataset(dataset)   
    
    if 'lev' not in ds[variable_name].dims:
        raise ValueError("The variable "+variable_name+" doesn't use the dimensions 'lat', 'lon', 'lev'")
        return 0


    # Extract variable attributes
    variable_unit = ds[variable_name].attrs.get('units', 'N/A')
    variable_long_name = ds[variable_name].attrs.get('long_name', 'N/A')
    selected_ut = ds['ut'].sel(time=selected_time).values.item() / (1e9 * 3600)
    selected_mtime = get_mtime(ds,selected_time)
    


    lev_ilev = 'lev'
    # Extract the data for the given selected_time and lev
    if selected_lev in ds['lev'].values:
        data = ds[variable_name].sel(time=selected_time, lev=selected_lev).values
        unique_lons = data.lon.values
        unique_lats = data.lat.values
        var_lat_lon = data.values
    else:
        print(f"The {lev_ilev} {selected_lev} isn't in the listed valid values.")
        sorted_levs = sorted(ds['lev'].values, key=lambda x: abs(x - selected_lev))
        closest_lev1 = sorted_levs[0]
        closest_lev2 = sorted_levs[1]
        print(f"Averaging from the closest valid {lev_ilev}s: {closest_lev1} and {closest_lev2}")
        # Extract data for the two closest lev values using .sel()
        data1 = ds[variable_name].sel(time=selected_time, lev=closest_lev1)
        unique_lons = data1.lon.values
        unique_lats = data1.lat.values
        var_lat_lon_1 = data1.values

        data2 = ds[variable_name].sel(time=selected_time, lev=closest_lev2)
        var_lat_lon_2 = data2.values
        # Return the averaged data
        var_lat_lon = (var_lat_lon_1 + var_lat_lon_2) / 2

    return var_lat_lon, selected_lev, unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime


def lat_lon_ilev(ds, variable_name, selected_time, selected_ilev):
    """
    Extract data from the dataset based on the given variable name, timestamp, and ilev value.
    
    Args:
    - variable_name (str): Name of the variable to extract.
        - valid variables: ['WN', 'NE', 'POTEN', 'UI_ExB', 'VI_ExB', 'WI_ExB', 'DEN', 'Z', 'ZG']
    - selected_time (str): Timestamp to filter the data.
    - selected_ilev (float): ilevel value to filter the data.
    
    Returns:
    - var_lat_lon (array): 2D array of [variable values, lat, lon] for the given timestamp and ilev.
    - variable_unit (str): Unit of the variable.
    - variable_long_name (str): Long name of the variable.
    - selected_ut (float): UT value in hours for selected_time.
    """

    # Load the dataset using xarray
    #ds = xr.open_dataset(dataset)   
    
    if 'ilev' not in ds[variable_name].dims:
        raise ValueError("The variable "+variable_name+" doesn't use the dimensions 'lat', 'lon', 'ilev'")
        return 0
 
    

    # Extract variable attributes
    variable_unit = ds[variable_name].attrs.get('units', 'N/A')
    variable_long_name = ds[variable_name].attrs.get('long_name', 'N/A')
    selected_ut = ds['ut'].sel(time=selected_time).values.item() / (1e9 * 3600)
    selected_mtime=get_mtime(ds,selected_time)

    lev_ilev = 'lev'
    # Extract the data for the given selected_time and lev
    if selected_ilev in ds['ilev'].values:
        data = ds[variable_name].sel(time=selected_time, ilev=selected_ilev).values
        unique_lons = data.lon.values
        unique_lats = data.lat.values
        var_lat_lon = data.values
    else:
        print(f"The {lev_ilev} {selected_ilev} isn't in the listed valid values.")
        sorted_levs = sorted(ds['ilev'].values, key=lambda x: abs(x - selected_ilev))
        closest_lev1 = sorted_levs[0]
        closest_lev2 = sorted_levs[1]
        print(f"Averaging from the closest valid {lev_ilev}s: {closest_lev1} and {closest_lev2}")
        # Extract data for the two closest lev values using .sel()
        data1 = ds[variable_name].sel(time=selected_time, ilev=closest_lev1)
        unique_lons = data1.lon.values
        unique_lats = data1.lat.values
        var_lat_lon_1 = data1.values

        data2 = ds[variable_name].sel(time=selected_time, ilev=closest_lev2)
        var_lat_lon_2 = data2.values
        # Return the averaged data
        var_lat_lon = (var_lat_lon_1 + var_lat_lon_2) / 2

    return var_lat_lon, selected_ilev, unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime


def var_lev(ds, variable_name, selected_time, selected_lat, selected_lon): #make var_ilev next
    """
    Extracts data from the dataset for a given variable name, latitude, longitude, and time.
    
    Parameters:
    - variable_name (str): Name of the variable to retrieve.
    - selected_lat (float): Latitude value.
    - selected_lon (float): Longitude value.
    - selected_time (int): Index of the time dimension.
    
    Returns:
    - list: A list of [lev, var_val] pairs.
    """
    # Load the dataset using xarray
    #ds = xr.open_dataset(dataset) 

    # Extract the variable data for the specified time, latitude, and longitude
    data = ds[variable_name].sel(time=selected_time, lat=selected_lat, lon=selected_lon, method="nearest")
    
    variable_unit = ds[variable_name].attrs.get('units', 'N/A')
    variable_long_name = ds[variable_name].attrs.get('long_name', 'N/A')
    selected_ut = ds['ut'].sel(time=selected_time).values.item() / (1e9 * 3600)
    selected_mtime=get_mtime(ds,selected_time)
    # Create a list of [var_val, lev] pairs
    #var_lev_arr = [[data_val.item(), lev_val] for data_val, lev_val in zip(data.values, ds['lev'].values)]
    

    #valid_indices = ~np.isnan(data.values)
    var_values = data.values#[valid_indices]
    lev_values = ds['lev'].values#[valid_indices]

    return var_values , lev_values, variable_unit, variable_long_name, selected_ut, selected_mtime 





def valid_lev_ilev(ds, selected_lev_ilev, lev_ilev):
    """
    Checks if the given lev or ilev value is in the dataset and if not gets the two closest valid values.
    
    Parameters:
    - ds (xarray): The loaded dataset opened using xarray.
    - selected_lev_ilev (float): Desired lev or ilev value.
    - lev_ilev (str): Dimension name ('lev' or 'ilev') in the dataset.
    
    Returns:
    - selected_levs (list): A list containing one or two values of lev or ilev present in the dataset.
    """
    valid_levs_ilevs = ds[lev_ilev].values

    if selected_lev_ilev not in valid_levs_ilevs:
        print(f"The {lev_ilev} {selected_lev_ilev} isn't in the listed valid values.")
        valid_levs_ilevs = np.asarray(valid_levs_ilevs)
        idx = (np.abs(valid_levs_ilevs - selected_lev_ilev)).argsort()[:2]
        selected_levs_ilevs = valid_levs_ilevs[idx]
        print(f"Averaging from the closest valid {lev_ilev}s: {selected_levs_ilevs[0]} and {selected_levs_ilevs[1]}")
    else:
        selected_levs_ilevs = [selected_lev_ilev]

    return selected_levs_ilevs





def calc_avg_ht(ds, selected_time, lev):
    """
    Compute the average Z value for a given set of lat, lon, and lev from a dataset.
    
    Parameters:
    - data_array (list): A list of [varval, lat, lon].
    - lev (float): The level for which to retrieve data.
    - dataset: The dataset containing the Z variable.
    
    Returns:
    - float: The average Z value for the given conditions.
    """
    if lev in ds['ilev'].values:
        heights = ds['ZG'].sel(time=selected_time, ilev=lev).values
    else:
        sorted_levs = sorted(ds['lev'].values, key=lambda x: abs(x - lev))
        closest_lev1 = sorted_levs[0]
        closest_lev2 = sorted_levs[1]
        
        # Extract data for the two closest lev values using .sel()
        data1 = ds['ZG'].sel(time=selected_time, ilev=closest_lev1).values
        data2 = ds['ZG'].sel(time=selected_time, ilev=closest_lev2).values
        
        # Return the averaged data
        heights = (data1 + data2) / 2
    avg_ht= round(heights.mean()/ 100000, 2)
    return avg_ht

def min_max(data):
    """Find the minimum and maximum values of varval from the 2D array
    
    Parameters:
    - data_array (list): A list of [varval, lat, lon].
    
    Returns:
    - min_val (float): Minimum value of the variable in the array.
    - max_val (float): Maximum value of the variable in the array.
    """
    
    return np.nanmin(data), np.nanmax(data)



def get_mtime(ds, selected_time):
    """Find the mtime array for the corresponding selected time.
    
    Parameters:
    - ds (xarray): The loaded dataset opened using xarray.
    - selected_time (str): Timestamp to filter the data.
    
    Returns:
    - array: mtime array containing [Day,Hour,Min].
    """
    # Convert input string to numpy datetime64 format
    selected_time = np.datetime64(selected_time)

    # Extract time and mtime variables from the dataset
    time_variable = ds['time'].values
    mtime_variable = ds['mtime'].values

    # Find the index corresponding to the input time
    index = np.where(time_variable == selected_time)

    # Extract and return the corresponding mtime value
    if index[0].size > 0:
        return mtime_variable[index[0][0]]
    else:
        return None


def avg_var_lat_lon(var_lat_lon1, var_lat_lon2):
    """
    Computes the average variable value for each lat-lon pair from two var_lat_lon arrays.
    
    Args:
    - var_lat_lon1 (array): First list of [variable values, lat, lon].
    - var_lat_lon2 (array): Second list of [variable values, lat, lon].
    
    Returns:
    - avg_var_lat_lon (array): List of averaged [variable values, lat, lon] for each lat-lon pair.
    """
    
    # Create a dictionary with lat-lon pairs as keys and their values as lists
    lat_lon_dict = {}
    
    for var_val, lat, lon in var_lat_lon1:
        if (lat, lon) in lat_lon_dict:
            lat_lon_dict[(lat, lon)].append(var_val)
        else:
            lat_lon_dict[(lat, lon)] = [var_val]
    
    for var_val, lat, lon in var_lat_lon2:
        if (lat, lon) in lat_lon_dict:
            lat_lon_dict[(lat, lon)].append(var_val)
        else:
            lat_lon_dict[(lat, lon)] = [var_val]
            
    # Compute the average for each lat-lon pair
    avg_var_lat_lon = [[sum(lat_lon_dict[(lat, lon)]) / len(lat_lon_dict[(lat, lon)]), lat, lon] 
                       for lat, lon in lat_lon_dict]
    
    return avg_var_lat_lon



def get_avg_ht_arr(ds, time, lat, lon):
    """
    Extracts ZG values for a given time, latitude, and longitude.
    
    Args:
    - time (str): Desired time in the format 'YYYY-MM-DDTHH:MM:SS'
    - lat (float): Desired latitude value
    - lon (float): Desired longitude value
    - dataset (xarray.Dataset): The dataset containing the ZG values. Default is the loaded dataset.
    
    Returns:
    - list: A list of lists where each inner list contains [ilev_val, ZG_val]
    """

    #ds = xr.open_dataset(dataset)
    # Extract the ZG values for the specified time, latitude, and longitude
    selected_zg = ds['ZG'].sel(time=time, lat=lat, lon=lon)

    # Convert the values from cm to km
    selected_zg_km = selected_zg / 100000  # 1 km = 100000 cm
    
    # Extract the ilev values for the same selection
    ilev_values = selected_zg['ilev'].values

    # Combine the ilev and ZG values into a single list
    combined_values = list(zip(ilev_values, selected_zg_km.values))

    averaged_array = []
    
    # Iterate over the zg_array and calculate the average for consecutive values
    for i in range(len(combined_values) - 1):
        avg_lev = (combined_values[i][0] + combined_values[i+1][0]) / 2
        avg_zg = (combined_values[i][1] + combined_values[i+1][1]) / 2
        averaged_array.append((avg_lev, avg_zg))
    averaged_array.append((7.25, float('nan')))
    zg_values_array = [item[1] for item in averaged_array]

    return zg_values_array

    #return combined_values
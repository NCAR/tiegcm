import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def lat_lon_lev(dataset, variable_name, selected_time, selected_lev):
    """Extracts a 2D array of variable values, latitude, and longitude for a given variable name, time, and lev.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lon, lev dimensions.
            - Valid variables:['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM']
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_lev (float): The selected level value.
    
    Returns:
        - 2D array containing [variable values, latitude, longitude] for the selected time and level.
    """
    valid_variables = ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM']
    if variable_name not in valid_variables:
        raise ValueError("The variable "+variable_name+" doesn't use the dimensions 'lat', 'lon', 'lev'")
    valid_lev= [-6.75, -6.25, -5.75, -5.25, -4.75, -4.25, -3.75, -3.25, -2.75,
       -2.25, -1.75, -1.25, -0.75, -0.25,  0.25,  0.75,  1.25,  1.75,
        2.25,  2.75,  3.25,  3.75,  4.25,  4.75,  5.25,  5.75,  6.25,
        6.75,  7.25]
    if selected_lev not in valid_lev:
        print("The lev "+ str(selected_lev)+" isn't in the listed valid values")
        valid_lev = np.asarray(valid_lev)
        idx = (np.abs(valid_lev - selected_lev)).argmin()
        selected_lev = valid_lev[idx]
        print("Selecting the closest valid lev "+ str(selected_lev))   
    # Open the NetCDF file
    ds = nc.Dataset(dataset, 'r')

    variable_unit=ds.variables[variable_name].units
    variable_long_name=ds.variables[variable_name].long_name
    # Extract time, latitude, longitude, level and variable_unit variables
    variable_unit=ds.variables[variable_name].units
    time_var = ds.variables['time']
    lat_var = ds.variables['lat']
    lon_var = ds.variables['lon']
    lev_var = ds.variables['lev']
    
    # Convert time variable to string format for comparison
    time_values = nc.num2date(time_var[:], time_var.units)
    time_str_values = [t.isoformat() for t in time_values]
    
    # Find the index of the selected time and level
    time_idx = time_str_values.index(selected_time)
    lev_idx = list(lev_var[:]).index(selected_lev)
    
    # Extract the data for the given variable, time, and level
    data_var = ds.variables[variable_name]
    data_values = data_var[time_idx, lev_idx, :, :]
    
    # Create a 2D array of [variable values, lat, lon]
    result = []
    for i, lat in enumerate(lat_var[:]):
        for j, lon in enumerate(lon_var[:]):
            result.append([data_values[i][j], lat, lon])
    
    # Close the dataset
    ds.close()
    
    return result,selected_lev, variable_unit,variable_long_name



def lat_lon_ilev(dataset, variable_name, selected_time, selected_ilev):
    """Extracts a 2D array of variable values, latitude, and longitude for a given variable name, time, and ilev.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lon, ilev dimensions.
            - Valid variables:['NE', 'OMEGA', 'Z', 'POTEN']
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_ilev (float): The selected ilevel value.
    
    Returns:
        - 2D array containing [variable values, latitude, longitude] for the selected time and ilevel.
    """
    valid_variables = ['NE', 'OMEGA', 'Z', 'POTEN']
    if variable_name not in valid_variables:
        raise ValueError("The variable "+variable_name+" doesn't use the dimensions 'lat', 'lon', 'ilev'")

    valid_ilev= [-7. , -6.5, -6. , -5.5, -5. , -4.5, -4. , -3.5, -3. , -2.5, -2. ,
       -1.5, -1. , -0.5,  0. ,  0.5,  1. ,  1.5,  2. ,  2.5,  3. ,  3.5,
        4. ,  4.5,  5. ,  5.5,  6. ,  6.5,  7. ]
    if selected_ilev not in valid_ilev:
        print("The lev "+ str(selected_ilev)+" isn't in the listed valid values")
        valid_ilev = np.asarray(valid_ilev)
        idx = (np.abs(valid_ilev - selected_ilev)).argmin()
        selected_ilev = valid_ilev[idx]
        print("Selecting the closest valid lev "+ str(selected_ilev))  

    # Open the NetCDF file
    ds = nc.Dataset(dataset, 'r')
    
    variable_unit=ds.variables[variable_name].units
    variable_long_name=ds.variables[variable_name].long_name
    # Extract time, latitude, longitude, and ilevel variables
    time_var = ds.variables['time']
    lat_var = ds.variables['lat']
    lon_var = ds.variables['lon']
    ilev_var = ds.variables['ilev']
    
    # Convert time variable to string format for comparison
    time_values = nc.num2date(time_var[:], time_var.units)
    time_str_values = [t.isoformat() for t in time_values]
    
    # Find the index of the selected time and ilevel
    time_idx = time_str_values.index(selected_time)
    ilev_idx = list(ilev_var[:]).index(selected_ilev)
    
    # Extract the data for the given variable, time, and ilevel
    data_var = ds.variables[variable_name]
    data_values = data_var[time_idx, ilev_idx, :, :]
    
    # Create a 2D array of [variable values, lat, lon]
    result = []
    for i, lat in enumerate(lat_var[:]):
        for j, lon in enumerate(lon_var[:]):
            result.append([data_values[i][j], lat, lon])
    
    # Close the dataset
    ds.close()
    
    return result,selected_ilev, variable_unit,variable_long_name


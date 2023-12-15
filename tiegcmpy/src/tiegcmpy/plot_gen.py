import sys
import numpy as np
import matplotlib.pyplot as plt
from .data_parse import arr_lat_lon,arr_lev_var,arr_lev_lon, arr_lev_lat,arr_lev_time,arr_lat_time, calc_avg_ht, min_max, get_time
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
from datetime import datetime, timezone
import matplotlib.ticker as mticker


def longitude_to_local_time(longitude):
    """
    Convert longitude to local time.
    
    Parameters:
        - longitude (float): Longitude value.
    
    Returns:
        - local_time (float): Local time corresponding to the given longitude.
    """
    local_time = (longitude / 15) % 24
    return local_time

def local_time_to_longitude(local_time):
    """
    Convert local time to longitude.
    
    Parameters:
        - local_time (float): Local time value.
    
    Returns:
        - longitude (float): Longitude corresponding to the given local time.
    """
    if local_time == 'mean':
        longitude = 'mean'
    else:
        #
        # Each hour of local time corresponds to 15 degrees of longitude
        #
        longitude = (local_time * 15) % 360
        #
        # Adjusting the longitude to be between -180 and 180 degrees
        #
        if longitude > 180:
            longitude = longitude - 360

    return longitude

def color_scheme(variable_name):
    """
    Sets color scheme for plots.
    
    Parameters:
        - variable_name (str): The name of the variable with latitude, longitude, ilev dimensions.
    
    Returns:
        - cmap_color (str): Color scheme of the countour map.
        - line_color (str): Color scheme of conutour lines.
    """
    #
    # Setting type of variable 
    #
    density_type = ['NE', 'DEN', 'O2', 'O1', 'N2', 'NO', 'N4S', 'HE']
    temp_type = ['TN', 'TE', 'TI', 'QJOULE']
    wind_type = ['WN', 'UI_ExB', 'VI_ExB', 'WI_ExB', 'UN', 'VN']
    #
    # Color scheme for density type variables
    #
    if variable_name in density_type:
        cmap_color = 'viridis'
        line_color = 'white'
    #
    # Color scheme for temprature type variables
    #
    elif variable_name in temp_type:
        cmap_color = 'inferno'
        line_color = 'white'
    #
    # Color scheme for wind type variables
    #
    elif variable_name in wind_type:
        cmap_color = 'bwr'
        line_color = 'black'
    #
    # Color scheme for all other types of variables
    #
    else:
        cmap_color = 'viridis'
        line_color = 'white'
    return cmap_color, line_color


def test(datasets, variable_name, time, level=None, selected_unit=None):

    variable_values, level, lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime, filename = arr_lat_lon(
        datasets,
        variable_name,
        time,
        level,
        selected_unit,
        plot_mode=True
    )

    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, edgecolor='white', linewidth=3)

    # Contour plot
    contour = ax.contourf(unique_lons, lats, variable_values, transform=ccrs.PlateCarree())
    contour_lines = ax.contour(unique_lons, lats, variable_values, colors='white', transform=ccrs.PlateCarree())
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical')
    cbar.set_label(variable_unit)  # Label for the colorbar
    plt.title(f'{variable_long_name} ({variable_unit})\nTime: {time}, Level: {level}')

    # Setting labels and ticks
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xticks([value for value in unique_lons if value % 30 == 0], crs=ccrs.PlateCarree())  # Adjust the range and step as per your data
    ax.set_yticks(range(-90, 91, 30), crs=ccrs.PlateCarree())   # Adjust the range and step as per your data

    plt.show()

    return fig

def plt_lat_lon(datasets, variable_name, time= None, mtime=None, level = None,  variable_unit = None, contour_intervals = None, contour_value = None, cmap_color = None, line_color = 'white', coastlines=False, nightshade=False, gm_equator=False, latitude_minimum = None, latitude_maximum = None, longitude_minimum = None, longitude_maximum = None, localtime_minimum = None, localtime_maximum = None ):

    """
    Generates a Latitude vs Longitude contour plot for a variable.
    
    Parameters:
        datasets (xarray): The loaded dataset/s using xarray.
        variable_name (str): The name of the variable with latitude, longitude, ilev dimensions.
        time (np.datetime64, optional): The selected time e.g., '2022-01-01T12:00:00'.
        mtime (array, optional): The selected time as a list e.g., [1, 12, 0] for 1st day, 12 hours, 0 mins.
        level (float, optional): The selected lev/ilev value.
        variable_unit (str, optional): The desired unit of the variable.
        contour_intervals (int, optional): The number of contour intervals. Defaults to 20.
        contour_value (int, optional): The value between each contour interval.
        cmap_color (str, optional): The color map of the conutour. Defaults to 'viridis' for Density,'inferno' for Temp, 'bwr' for Wind, 'viridis' for undefined.
        line_color (str, optional): The color for all lines in the on the plot. Defaults to 'white'.
        coastlines (bool, optional): Shows coastlines on the plot. Defaults to False.
        nightshade (bool, optional): Shows nighshade on the plot. Defaults to False.
        gm_equator (bool, optional): Shows geomagmetic equator on the plot. Defaults to False.
        latitude_minimum (float, optional): Minimum latitude to slice plots. Defaults to -87.5.
        latitude_maximum (float, optional): Maximum latitude to slice plots. Defaults to 87.5.
        longitude_minimum (float, optional): Minimum longitude to slice plots. Defaults to -180.
        longitude_maximum (float, optional): Maximum longitude to slice plots. Defaults to 175.
        localtime_minimum (float, optional): Minimum localtime to slice plots.
        localtime_maximum (float, optional): Maximum localtime to slice plots.
    
    Returns:
        Contour plot.
    """
    # Printing Execution data
    if time == None:
        time = get_time(datasets, mtime)
    if localtime_minimum != None:
        longitude_minimum = local_time_to_longitude(localtime_minimum)
    if localtime_maximum != None:
        longitude_maximum = local_time_to_longitude(localtime_maximum)
    if contour_intervals == None:
        contour_intervals = 20
    print("---------------["+variable_name+"]---["+str(time)+"]---["+str(level)+"]---------------")
    # Generate 2D arrays, extract variable_unit
    '''
    if level != None:
        try:
            data, level,  unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime, filename =lat_lon_lev(datasets, variable_name, time, level, variable_unit)
        except ValueError:
            data, level,  unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime, filename =lat_lon_ilev(datasets, variable_name, time, level, variable_unit)
        if level != 'mean':
            avg_ht=calc_avg_ht(datasets, time,level)
    else:
        data, unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime, filename =lat_lon(datasets, variable_name, time)
    '''
    if isinstance(time, str):
        time = np.datetime64(time, 'ns')

    data, level,  unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime, filename =arr_lat_lon(datasets, variable_name, time, selected_lev_ilev = level, selected_unit = variable_unit, plot_mode = True)
    if level != 'mean' and level != None:
            avg_ht=calc_avg_ht(datasets, time,level)
    if latitude_minimum == None:
        latitude_minimum = np.nanmin(unique_lats)
    if latitude_maximum == None:
        latitude_maximum = np.nanmax(unique_lats)
    if longitude_minimum == None:
        longitude_minimum = np.nanmin(unique_lons)
    if longitude_maximum == None:   
        longitude_maximum = np.nanmax(unique_lons)

    min_val, max_val = min_max(data)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]

    if cmap_color == None:
        cmap_color, line_color = color_scheme(variable_name)
    # Extract values, latitudes, and longitudes from the array
    if contour_value is not None:
        contour_levels = np.arange(min_val, max_val + contour_value, contour_value)
    else:
        contour_levels = np.linspace(min_val, max_val, contour_intervals)
    # Generate contour plot
    
    interval_value = contour_value if contour_value else (max_val - min_val) / (contour_intervals - 1)

    # Generate contour plot
    plot = plt.figure(figsize=(20, 9))

    
    # Check if add_coastlines parameter is True
    if coastlines:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, edgecolor=line_color, linewidth=3)
    if nightshade:

        print(datetime.fromtimestamp(time.astype('O')/1e9, tz=timezone.utc))
        ax.add_feature(Nightshade(datetime.fromtimestamp(time.astype('O')/1e9, tz=timezone.utc), alpha=0.4))
    if gm_equator:
        ax.plot(unique_lons, [0]*len(unique_lons), color=line_color, linestyle='--', transform=ccrs.Geodetic())
    else:
        ax = plt.gca()
    
    contour_filled = plt.contourf(unique_lons, unique_lats, data, cmap=cmap_color, levels=contour_levels)
    contour_lines = plt.contour(unique_lons, unique_lats, data, colors=line_color, linewidths=0.5, levels=contour_levels)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=line_color)
    cbar = plt.colorbar(contour_filled, label=variable_name + " [" + variable_unit + "]",fraction=0.046, pad=0.04)
    cbar.set_label(variable_name + " [" + variable_unit + "]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    
    
    plt.xlabel('Longitude (Deg)', fontsize=28)
    plt.xticks([value for value in unique_lons if value % 30 == 0],fontsize=18)
    plt.ylabel('Latitude (Deg)', fontsize=28)
    plt.yticks(range(-90, 91, 30),fontsize=18)
    plt.xlim(longitude_minimum,longitude_maximum)
    plt.ylim(latitude_minimum,latitude_maximum)

    plt.tight_layout()
    
    if coastlines == False and nightshade == False and gm_equator == False:
        # Add Local Time secondary x-axis
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2_xticks = ax.get_xticks()
        ax2.set_xticks(ax2_xticks)
        ax2.set_xticklabels([str(int(longitude_to_local_time(longitude) % 24)) for longitude in ax2_xticks],fontsize=18)
        ax2.set_xlabel('Local Time (Hrs)', labelpad=15, fontsize=28)


    # Add plot title
    plt.title(variable_long_name + ' ' + variable_name + ' (' + variable_unit + ') ' + '\n\n', fontsize=36)
    # Add plot subtitle
    if level == 'mean':
        plt.text(0, 115, 'UT=' + str(selected_ut) + '  ZP=' + str(level), ha='center', va='center', fontsize=28)
    elif level != None:
        plt.text(0, 115, 'UT=' + str(selected_ut) + '  ZP=' + str(level)+' AVG HT=' + str(avg_ht), ha='center', va='center', fontsize=28)
    else:
        plt.text(0, 115, 'UT=' + str(selected_ut), ha='center', va='center', fontsize=28)
    

    # Add subtext to the plot
    plt.text(-90, -115, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28)
    plt.text(90, -115, "Contour Interval = "+str("{:.2e}".format(interval_value)), ha='center', va='center',fontsize=28)
    plt.text(0, -125, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28)
    plt.text(0, -135, str(filename), ha='center', va='center',fontsize=28)
    
    return(plot)



def plt_lev_var(datasets, variable_name, latitude, time= None, mtime=None, longitude = None, localtime = None, variable_unit = None, level_minimum = None, level_maximum = None):
    """
    Generates a Level vs Variable line plot for a given latitude.
    
    Parameters:
        datasets (xarray): The loaded dataset/s using xarray.
        variable_name (str): The name of the variable with latitude, longitude, and ilev dimensions.
        latitude (float): The specific latitude value for the plot.
        time (np.datetime64, optional): The selected time, e.g., '2022-01-01T12:00:00'.
        mtime (array, optional): The selected time as a list, e.g., [1, 12, 0] for 1st day, 12 hours, 0 mins.
        longitude (float, optional): The specific longitude value for the plot.
        localtime (float, optional): The specific local time value for the plot.
        variable_unit (str, optional): The desired unit of the variable.
        level_minimum (float, optional): Minimum level value for the plot. Defaults to -8.
        level_maximum (float, optional): Maximum level value for the plot. Defaults to 8.
    
    Returns:
        Line plot.
    """
    # Printing Execution data
    if time == None:
        time = get_time(datasets, mtime)
    if longitude == None:
        longitude = local_time_to_longitude(localtime)
    print("---------------["+variable_name+"]---["+str(time)+"]---["+str(latitude)+"]---["+str(longitude)+"]---------------")


    variable_values , levs_ilevs, variable_unit, variable_long_name, selected_ut, selected_mtime, filename = arr_lev_var(datasets, variable_name, time, latitude, longitude,  variable_unit, plot_mode = True)

    if level_minimum == None:
        level_minimum = np.nanmin(levs_ilevs)
    if level_maximum == None:
        level_maximum = np.nanmax(levs_ilevs)

    min_val, max_val = min_max(variable_values)
    #print(min_val, max_val)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]


    #print(len(zg_values))

    
    # Plotting
    
    plot = plt.figure(figsize=(22, 12))
    plt.plot(variable_values, levs_ilevs)
    plt.xlabel(variable_long_name, fontsize=28, labelpad=15)
    plt.ylabel('LN(P0/P) (INTERFACES)', fontsize=28)
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   

    plt.ylim(level_minimum, level_maximum)



    if longitude == 'mean' and latitude == 'mean':
        plt.text(0.5, 1.08,'UT='+str(selected_ut) +"  LAT= Mean SLT= Mean", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    elif longitude == 'mean':
        plt.text(0.5, 1.08,'UT='+str(selected_ut) +'  LAT='+str(latitude)+" SLT= Mean", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    elif latitude == 'mean':
        plt.text(0.5, 1.08,'UT='+str(selected_ut) +'  LAT= Mean'+" SLT="+str(longitude_to_local_time(longitude))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    else:
        plt.text(0.5, 1.08,'UT='+str(selected_ut) +'  LAT='+str(latitude)+" SLT="+str(longitude_to_local_time(longitude))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.25, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.5, -0.3, str(filename), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    return(plot)


def plt_lev_lon(datasets, variable_name, latitude, time= None, mtime=None, variable_unit = None, contour_intervals = 20, contour_value = None, cmap_color = None, line_color = 'white',  level_minimum = None, level_maximum = None, longitude_minimum = None, longitude_maximum = None, localtime_minimum = None, localtime_maximum = None):
    """
    Generates a Level vs Longitude contour plot for a given latitude.
    
    Parameters:
        datasets (xarray): The loaded dataset/s using xarray.
        variable_name (str): The name of the variable with latitude, longitude, and ilev dimensions.
        latitude (float): The specific latitude value for the plot.
        time (np.datetime64, optional): The selected time, e.g., '2022-01-01T12:00:00'.
        mtime (array, optional): The selected time as a list, e.g., [1, 12, 0] for 1st day, 12 hours, 0 mins.
        variable_unit (str, optional): The desired unit of the variable.
        contour_intervals (int, optional): The number of contour intervals. Defaults to 20.
        contour_value (int, optional): The value between each contour interval.
        cmap_color (str, optional): The color map of the conutour. Defaults to 'viridis' for Density,'inferno' for Temp, 'bwr' for Wind, 'viridis' for undefined.
        line_color (str, optional): The color for all lines in the on the plot. Defaults to 'white'.
        level_minimum (float, optional): Minimum level value for the plot. Defaults to -6.75.
        level_maximum (float, optional): Maximum level value for the plot. Defaults to 6.75.
        longitude_minimum (float, optional): Minimum longitude value for the plot. Defaults to -180.
        longitude_maximum (float, optional): Maximum longitude value for the plot. Defaults to 175.
        localtime_minimum (float, optional): Minimum localtime value for the plot.
        localtime_maximum (float, optional): Maximum localtime value for the plot.
    
    Returns:
        Contour plot.
    """
    # Printing Execution data
    if time == None:
        time = get_time(datasets, mtime)
    if localtime_minimum != None:
        longitude_minimum = local_time_to_longitude(localtime_minimum)
    if localtime_maximum != None:
        longitude_maximum = local_time_to_longitude(localtime_maximum)
    if contour_intervals == None:
        contour_intervals = 20    
    print("---------------["+variable_name+"]---["+str(time)+"]---["+str(latitude)+"]---------------")
    # Generate 2D arrays, extract variable_unit
    variable_values, unique_lons, unique_levs,latitude, variable_unit, variable_long_name, selected_ut, selected_mtime, filename = arr_lev_lon(datasets, variable_name, time, latitude, variable_unit, plot_mode = True)

    if level_minimum == None:
        level_minimum = np.nanmin(unique_levs)
    if level_maximum == None:
        level_maximum = np.nanmax(unique_levs)
    if longitude_minimum == None:
        longitude_minimum = np.nanmin(unique_lons)
    if longitude_maximum == None:   
        longitude_maximum = np.nanmax(unique_lons)

    min_val, max_val = min_max(variable_values)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]

    if cmap_color == None:
        cmap_color, line_color = color_scheme(variable_name)

    if contour_value is not None:
        contour_levels = np.arange(min_val, max_val + contour_value, contour_value)
    else:
        contour_levels = np.linspace(min_val, max_val, contour_intervals)
    # Generate contour plot
    
    interval_value = contour_value if contour_value else (max_val - min_val) / (contour_intervals - 1)


    plot=plt.figure(figsize=(24, 12))
    contour_filled = plt.contourf(unique_lons, unique_levs, variable_values, cmap= cmap_color, levels=contour_levels)
    contour_lines = plt.contour(unique_lons, unique_levs, variable_values, colors=line_color, linewidths=0.5, levels=contour_levels)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=line_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    if latitude == 'mean':
        plt.text(0.5, 1.18,'UT='+str(selected_ut) +' ZONAL MEANS', ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    else:
        plt.text(0.5, 1.18,'UT='+str(selected_ut) +'  LAT='+str(latitude), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    plt.xlabel('Longitude (Deg)',fontsize=28)
    plt.xticks([value for value in unique_lons if value % 30 == 0],fontsize=18)  
    plt.yticks(fontsize=18) 
    plt.xlim(longitude_minimum,longitude_maximum)
    plt.ylim(level_minimum, level_maximum)

    # Add Local Time secondary x-axis
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2_xticks = ax.get_xticks()
    ax2.set_xticks(ax2_xticks)
    ax2.set_xticklabels([str(int(longitude_to_local_time(longitude) % 24)) for longitude in ax2_xticks],fontsize=18)
    ax2.set_xlabel('Local Time (Hrs)', labelpad=15, fontsize=28)

    # Add subtext to the plot
    plt.text(0.25, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.75, -0.2, "Contour Interval = "+str("{:.2e}".format(interval_value)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.75, -0.25, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.25, -0.25, str(filename), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    
    #plot, ax = plt.subplots()

    
    
    return(plot)


def plt_lev_lat(datasets, variable_name, time= None, mtime=None, longitude = None, localtime = None, variable_unit = None, contour_intervals = 20, contour_value = None, cmap_color = None, line_color = 'white', level_minimum = None, level_maximum = None, latitude_minimum = None,latitude_maximum = None):
    """
    Generates a Level vs Latitude contour plot for a specified time and/or longitude.
    
    Parameters:
        datasets (xarray): The loaded dataset/s using xarray.
        variable_name (str): The name of the variable with latitude, longitude, and ilev dimensions.
        time (np.datetime64, optional): The selected time, e.g., '2022-01-01T12:00:00'.
        mtime (array, optional): The selected time as a list, e.g., [1, 12, 0] for 1st day, 12 hours, 0 mins.
        longitude (float, optional): The specific longitude value for the plot.
        localtime (float, optional): The specific local time value for the plot.
        variable_unit (str, optional): The desired unit of the variable.
        contour_intervals (int, optional): The number of contour intervals. Defaults to 20.
        contour_value (int, optional): The value between each contour interval.
        cmap_color (str, optional): The color map of the conutour. Defaults to 'viridis' for Density,'inferno' for Temp, 'bwr' for Wind, 'viridis' for undefined.
        line_color (str, optional): The color for all lines in the on the plot. Defaults to 'white'.
        level_minimum (float, optional): Minimum level value for the plot. Defaults to -6.75.
        level_maximum (float, optional): Maximum level value for the plot. Defaults to 6.75.
        latitude_minimum (float, optional): Minimum latitude value for the plot. Defaults to -87.5.
        latitude_maximum (float, optional): Maximum latitude value for the plot. Defaults to 87.5.
    
    Returns:
        Contour plot.
    """
    # Printing Execution data
    if time == None:
        time = get_time(datasets, mtime)
    if longitude == None:
        longitude = local_time_to_longitude(localtime)
    if contour_intervals == None:
        contour_intervals = 20
    print("---------------["+variable_name+"]---["+str(time)+"]---["+str(longitude)+"]---------------")
    # Generate 2D arrays, extract variable_unit
    variable_values, unique_lats, unique_levs,longitude, variable_unit, variable_long_name, selected_ut, selected_mtime, filename = arr_lev_lat(datasets, variable_name, time, longitude,  variable_unit, plot_mode = True)

    if level_minimum == None:
        level_minimum = np.nanmin(unique_levs)
    if level_maximum == None:
        level_maximum = np.nanmax(unique_levs)
    if latitude_minimum == None:
        latitude_minimum = np.nanmin(unique_lats)
    if latitude_maximum == None:
        latitude_maximum = np.nanmax(unique_lats)

    min_val, max_val = min_max(variable_values)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]

    if cmap_color == None:
        cmap_color, line_color = color_scheme(variable_name)

    if contour_value is not None:
        contour_levels = np.arange(min_val, max_val + contour_value, contour_value)
    else:
        contour_levels = np.linspace(min_val, max_val, contour_intervals)
    
    
    interval_value = contour_value if contour_value else (max_val - min_val) / (contour_intervals - 1)
    
    # Generate contour plot
    plot=plt.figure(figsize=(24, 12))
    contour_filled = plt.contourf(unique_lats, unique_levs, variable_values, cmap= cmap_color, levels=contour_levels)
    contour_lines = plt.contour(unique_lats, unique_levs, variable_values, colors=line_color, linewidths=0.5, levels=contour_levels)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=line_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    if longitude == 'mean':
          plt.text(0.5, 1.08,'UT='+str(selected_ut) +' ZONAL MEANS', ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    else:
        plt.text(0.5, 1.08,'UT='+str(selected_ut) +'  LON='+str(longitude)+" SLT="+str(longitude_to_local_time(longitude))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    plt.xlabel('Latitude (Deg)',fontsize=28)
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 
    plt.xlim(latitude_minimum,latitude_maximum)
    plt.ylim(level_minimum,level_maximum)
    
    # Add subtext to the plot
    plt.text(0.25, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.75, -0.2, "Contour Interval = "+str("{:.2e}".format(interval_value)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.50, -0.25, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.50, -0.3, str(filename), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    
    #plot, ax = plt.subplots()

    
    
    return(plot)




def plt_lev_time(datasets, variable_name, latitude, longitude = None, localtime = None, variable_unit = None, contour_intervals = 20, contour_value = None, cmap_color = None, line_color = 'white',  level_minimum = None, level_maximum = None):
    """
    Generates a Level vs Time contour plot for a specified latitude and/or longitude.
    
    Parameters:
        datasets (xarray): The loaded dataset/s using xarray.
        variable_name (str): The name of the variable with latitude, longitude, time, and ilev dimensions.
        latitude (float): The specific latitude value for the plot.
        longitude (float, optional): The specific longitude value for the plot.
        localtime (float, optional): The specific local time value for the plot.
        variable_unit (str, optional): The desired unit of the variable.
        contour_intervals (int, optional): The number of contour intervals. Defaults to 20.
        contour_value (int, optional): The value between each contour interval.
        cmap_color (str, optional): The color map of the conutour. Defaults to 'viridis' for Density,'inferno' for Temp, 'bwr' for Wind, 'viridis' for undefined.
        line_color (str, optional): The color for all lines in the on the plot. Defaults to 'white'.
        level_minimum (float, optional): Minimum level value for the plot. Defaults to -6.75.
        level_maximum (float, optional): Maximum level value for the plot. Defaults to 6.75.
    
    Returns:
        Contour plot.
    """
    if longitude is None:
        longitude = local_time_to_longitude(localtime)
    if contour_intervals == None:
        contour_intervals = 20
    #print(datasets)
    variable_values_all, levs_ilevs, mtime_values, longitude, variable_unit, variable_long_name = arr_lev_time(datasets, variable_name, latitude, longitude, variable_unit, plot_mode = True)
    
    if level_minimum == None:
        level_minimum = np.nanmin(levs_ilevs)
    if level_maximum == None:
        level_maximum = np.nanmax(levs_ilevs)

    print("---------------["+variable_name+"]---["+str(latitude)+"]---["+str(longitude)+"]---------------")
    
    min_val, max_val = np.nanmin(variable_values_all), np.nanmax(variable_values_all)

    if cmap_color == None:
        cmap_color, line_color = color_scheme(variable_name)

    if contour_value is not None:
        contour_levels = np.arange(min_val, max_val + contour_value, contour_value)
    else:
        contour_levels = np.linspace(min_val, max_val, contour_intervals)
    
    
    interval_value = contour_value if contour_value else (max_val - min_val) / (contour_intervals - 1)
    mtime_tuples = [tuple(entry) for entry in mtime_values]
    try:    # Modify this part to show both day and hour
        unique_times = sorted(list(set([(day, hour) for day, hour, _, _ in mtime_values])))
        time_indices = [i for i, (day, hour, _, _) in enumerate(mtime_tuples) if i == 0 or mtime_tuples[i-1][:2] != (day, hour)]
    except:
        unique_times = sorted(list(set([day for day, _, _ in mtime_values])))
        time_indices = [i for i, (day, _, _) in enumerate(mtime_values) if i == 0 or mtime_values[i-1][0] != day]

    plot=plt.figure(figsize=(20, 12))
    X, Y = np.meshgrid(range(len(mtime_values)), levs_ilevs)
    contour_filled = plt.contourf(X, Y, variable_values_all, cmap=cmap_color, levels=contour_levels)
    contour_lines = plt.contour(X, Y, variable_values_all, colors=line_color, linewidths=0.5, levels=contour_levels)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=line_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    try:
        plt.xticks(time_indices, ["{}-{:02d}h".format(day, hour) for day, hour in unique_times], rotation=45)
        plt.xlabel("Model Time (Day-Hour) from "+str(unique_times[0])+" to "+str(unique_times[-1]), fontsize=28) 
    except:
        plt.xticks(time_indices, unique_times, rotation=45)
        plt.xlabel("Model Time (Day) from "+str(np.nanmin(unique_times))+" to "+str(np.nanmax(unique_times)) ,fontsize=28)
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    
    plt.tight_layout()
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 
    plt.ylim(level_minimum,level_maximum)


    # Add subtext to the plot
    if longitude == 'mean' and latitude == 'mean':
        plt.text(0.5, 1.08,'  LAT= Mean SLT= Mean', ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    elif longitude == 'mean':
        plt.text(0.5, 1.08,'  LAT='+str(latitude)+" SLT= Mean", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    elif latitude == 'mean':
        plt.text(0.5, 1.08,'  LAT= Mean'+" SLT="+str(longitude_to_local_time(longitude))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    else:
        plt.text(0.5, 1.08,'  LAT='+str(latitude)+" SLT="+str(longitude_to_local_time(longitude))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.5, -0.25, "Contour Interval = "+str("{:.2e}".format(interval_value)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    return(plot)



def plt_lat_time(datasets, variable_name, level = None, longitude = None, localtime = None,  variable_unit = None, contour_intervals = 10, contour_value = None, cmap_color = None, line_color = 'white', latitude_minimum = None,latitude_maximum = None):
    """
    Generates a Latitude vs Time contour plot for a specified level and/or longitude.
    
    Parameters:
        datasets (xarray): The loaded dataset/s using xarray.
        variable_name (str): The name of the variable with latitude, longitude, time, and ilev dimensions.
        level (float): The specific level value for the plot.
        longitude (float, optional): The specific longitude value for the plot.
        localtime (float, optional): The specific local time value for the plot.
        variable_unit (str, optional): The desired unit of the variable.
        contour_intervals (int, optional): The number of contour intervals. Defaults to 20.
        contour_value (int, optional): The value between each contour interval.
        latitude_minimum (float, optional): Minimum latitude value for the plot. Defaults to -87.5.
        latitude_maximum (float, optional): Maximum latitude value for the plot. Defaults to 87.5.
    
    Returns:
        Contour plot.
    """
    if longitude is None:
        longitude = local_time_to_longitude(localtime)
    if contour_intervals == None:
        contour_intervals = 20
    print("---------------["+variable_name+"]---["+str(level)+"]---["+str(longitude)+"]---------------")
    '''
    if level != None: 
        try:
            variable_values_all, unique_lats, mtime_values, longitude, variable_unit, variable_long_name, filename = lat_time_lev(datasets, variable_name, level, longitude, variable_unit)
        except:
            variable_values_all, unique_lats, mtime_values, longitude, variable_unit, variable_long_name, filename = lat_time_ilev(datasets, variable_name, level, longitude, variable_unit)
    else:
        variable_values_all, unique_lats, mtime_values, longitude, variable_unit, variable_long_name, filename = lat_time(datasets, variable_name, longitude, variable_unit)
    '''
    variable_values_all, unique_lats, mtime_values, longitude, variable_unit, variable_long_name, filename = arr_lat_time(datasets, variable_name, longitude, level, variable_unit, plot_mode = True)
    # Assuming the levels are consistent across datasets, but using the minimum size for safety
    
    if latitude_minimum == None:
        latitude_minimum = np.nanmin(unique_lats)
    if latitude_maximum == None:
        latitude_maximum = np.nanmax(unique_lats)

    
    
    min_val, max_val = np.nanmin(variable_values_all), np.nanmax(variable_values_all)
    
    if cmap_color == None:
        cmap_color, line_color = color_scheme(variable_name)

    if contour_value is not None:
        contour_levels = np.arange(min_val, max_val + contour_value, contour_value)
    else:
        contour_levels = np.linspace(min_val, max_val, contour_intervals)
    

    interval_value = contour_value if contour_value else (max_val - min_val) / (contour_intervals - 1)

    mtime_tuples = [tuple(entry) for entry in mtime_values]
    try:    # Modify this part to show both day and hour
        unique_times = sorted(list(set([(day, hour) for day, hour, _, _ in mtime_values])))
        time_indices = [i for i, (day, hour, _, _) in enumerate(mtime_tuples) if i == 0 or mtime_tuples[i-1][:2] != (day, hour)]
        if len(time_indices) >=10:
            unique_times = sorted(list(set([day for day, _, _, _ in mtime_values])))
            time_indices = [i for i, (day, _, _, _) in enumerate(mtime_values) if i == 0 or mtime_values[i-1][0] != day]
    except:
        unique_times = sorted(list(set([day for day, _, _ in mtime_values])))
        time_indices = [i for i, (day, _, _) in enumerate(mtime_values) if i == 0 or mtime_values[i-1][0] != day]

    plot = plt.figure(figsize=(20, 12))
    X, Y = np.meshgrid(range(len(mtime_values)), unique_lats)
    contour_filled = plt.contourf(X, Y, variable_values_all, cmap=cmap_color, levels=contour_levels)
    contour_lines = plt.contour(X, Y, variable_values_all, colors=line_color, linewidths=0.5, levels=contour_levels)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=line_color)
    cbar = plt.colorbar(contour_filled, label=variable_name + " [" + variable_unit + "]")
    cbar.set_label(variable_name + " [" + variable_unit + "]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    try:
        plt.xticks(time_indices, ["{}-{:02d}h".format(day, hour) for day, hour in unique_times], rotation=45)
        plt.xlabel("Model Time (Day-Hour) from "+str(unique_times[0])+" to "+str(unique_times[-1]), fontsize=28) 
    except:
        plt.xticks(time_indices, unique_times, rotation=45)
        plt.xlabel("Model Time (Day) from "+str(np.nanmin(unique_times))+" to "+str(np.nanmax(unique_times)) ,fontsize=28)
    
    plt.ylabel('Latitude (Deg)',fontsize=28)
    
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    
    plt.tight_layout()
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylim(latitude_minimum, latitude_maximum)

    # Add subtext to the plot
    if level == 'mean' and longitude == 'mean':
        plt.text(0.5, 1.08, '  ZP= Mean SLT= Mean', ha='center', va='center', fontsize=28, transform=plt.gca().transAxes)
    elif longitude == 'mean':
        plt.text(0.5, 1.08, '  ZP=' + str(level) + " SLT= Mean", ha='center', va='center', fontsize=28, transform=plt.gca().transAxes)
    elif level == 'mean':
        plt.text(0.5, 1.08, '  ZP= Mean' + " SLT=" + str(longitude_to_local_time(longitude)) + "Hrs", ha='center', va='center', fontsize=28, transform=plt.gca().transAxes)
    else:
        plt.text(0.5, 1.08, '  ZP=' + str(level) + " SLT=" + str(longitude_to_local_time(longitude)) + "Hrs", ha='center', va='center', fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.5, -0.2, "Min, Max = " + str("{:.2e}".format(min_val)) + ", " + str("{:.2e}".format(max_val)), ha='center', va='center', fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.5, -0.25, "Contour Interval = " + str("{:.2e}".format(interval_value)), ha='center', va='center', fontsize=28, transform=plt.gca().transAxes)
    plt.close(plot)
    return plot

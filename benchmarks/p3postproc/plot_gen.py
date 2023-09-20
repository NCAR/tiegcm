
import numpy as np
import matplotlib.pyplot as plt
from data_parse import lat_lon_lev, lat_lon_ilev,calc_avg_ht, min_max, lev_ilev_var, get_avg_ht_arr, lev_ilev_lon, lev_ilev_lat,lev_ilev_time, lat_time_lev,lat_time_ilev

def longitude_to_local_time(longitude):
    """
    Convert longitude to local time.
    
    Parameters:
        - longitude (float): Longitude value.
    
    Returns:
        - local_time (float): Local time corresponding to the given longitude.
    """
    # Longitude is given in degrees, and every 15 degrees corresponds to 1 hour difference in local time
    # Longitude 0 corresponds to local time 0
    local_time = (longitude / 15) % 24
    return local_time

def color_scheme(variable_name):
    density_type = ['NE', 'DEN', 'O2', 'O1', 'N2', 'NO', 'N4S', 'HE']
    temp_type = ['TN', 'TE', 'TI', 'QJOULE']
    wind_type = ['WN', 'UI_ExB', 'VI_ExB', 'WI_ExB', 'UN', 'VN']
    if variable_name in density_type:
        cmap_color = 'viridis'
        contour_color = 'white'
    elif variable_name in temp_type:
        cmap_color = 'inferno'
        contour_color = 'white'
    elif variable_name in wind_type:
        cmap_color = 'bwr'
        contour_color = 'black'
    else:
        cmap_color = 'viridis'
        contour_color = 'white'
    return cmap_color, contour_color

def plt_lat_lon(dataset, variable_name, selected_time, selected_lev_ilev):
    """
    Generates a contour plot for the given 2D array of variable values, latitude, and longitude.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lon, ilev dimensions.
            - Valid variables:['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 
                                'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 
                                'Z', 'POTEN']
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_ilev (float): The selected ilevel value.
    
    Returns:
        - Contour plot.
    """
    # Printing Execution data
    
    print("---------------["+variable_name+"]---["+str(selected_time)+"]---["+str(selected_lev_ilev)+"]---------------")
    # Generate 2D arrays, extract variable_unit

    try:
        data, selected_lev_ilev,  unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime =lat_lon_lev(dataset, variable_name, selected_time, selected_lev_ilev)
    except ValueError:
        data, selected_lev_ilev,  unique_lats, unique_lons, variable_unit, variable_long_name, selected_ut, selected_mtime=lat_lon_ilev(dataset, variable_name, selected_time, selected_lev_ilev)

    
    avg_ht=calc_avg_ht(dataset, selected_time,selected_lev_ilev)
    min_val, max_val = min_max(data)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]

    cmap_color, contour_color = color_scheme(variable_name)
    # Extract values, latitudes, and longitudes from the array

    # Generate contour plot
    plot = plt.figure(figsize=(20, 12))
    contour_filled = plt.contourf(unique_lons, unique_lats, data, cmap=cmap_color, levels=20)
    contour_lines = plt.contour(unique_lons, unique_lats, data, colors=contour_color, linewidths=0.5, levels=20)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=contour_color)
    cbar = plt.colorbar(contour_filled, label=variable_name + " [" + variable_unit + "]")
    cbar.set_label(variable_name + " [" + variable_unit + "]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.title(variable_long_name + ' ' + variable_name + ' (' + variable_unit + ') ' + '\n\n', fontsize=36)
    plt.text(0, 110, 'UT=' + str(selected_ut) + '  ZP=' + str(selected_lev_ilev)+' AVG HT=' + str(avg_ht), ha='center', va='center', fontsize=28)
    plt.ylabel('Latitude', fontsize=28)
    plt.xlabel('Longitude (Deg)', fontsize=28)
    plt.xticks([value for value in unique_lons if value % 30 == 0],fontsize=18)
    plt.yticks(fontsize=18)
    plt.tight_layout()
    plt.show()

    # Add Local Time secondary x-axis
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2_xticks = ax.get_xticks()
    ax2.set_xticks(ax2_xticks)
    ax2.set_xticklabels([str(int(longitude_to_local_time(lon) % 24)) for lon in ax2_xticks],fontsize=18)
    ax2.set_xlabel('Local Time (Hrs)', labelpad=15, fontsize=28)

    # Add subtext to the plot
    plt.text(-90, -115, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28)
    plt.text(90, -115, "Contour Interval = "+str("{:.2e}".format((max_val-min_val)/20)), ha='center', va='center',fontsize=28)
    plt.text(90, -125, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28)
    #plt.text(-90, -125, str(dataset.split("/")[-1]), ha='center', va='center',fontsize=28)

    plt.show()
    #plot, ax = plt.subplots()
    #plot, ax = plt.subplots()

    
    
    return(plot)



def plt_lev_var(dataset, variable_name, selected_time, selected_lat, selected_lon):
    """
    Plots the given data as a line plot.
    
    Parameters:
    - data_array (list): A list of [lev, var_val] pairs.
    
    Returns:
    - None
    """
    # Printing Execution data
    
    print("---------------["+variable_name+"]---["+str(selected_time)+"]---["+str(selected_lat)+"]---["+str(selected_lon)+"]---------------")


    variable_values , levs_ilevs, variable_unit, variable_long_name, selected_ut, selected_mtime = lev_ilev_var(dataset, variable_name, selected_time, selected_lat, selected_lon)

    min_val, max_val = min_max(variable_values)
    #print(min_val, max_val)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]


    #zg_values = get_avg_ht_arr(dataset, selected_time, selected_lat, selected_lon)
    #print(len(zg_values))

    
    # Plotting
    
    plot = plt.figure(figsize=(22, 12))
    plt.plot(variable_values, levs_ilevs)
    plt.xlabel(variable_long_name, fontsize=28, labelpad=15)
    plt.ylabel('LN(P0/P) (INTERFACES)', fontsize=28)
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   

    plt.ylim(-8, 8)

    '''
    ax = plt.gca()
    ax2 = ax.twinx()
    ax2.plot(variable_values, zg_values, 'r--', alpha=0)  # Plot with alpha=0 to make it invisible
    ax2.set_ylabel('Height (in km)', fontsize=28, labelpad=15, color='black')
    ax2.tick_params(axis='y', labelcolor='black', labelsize=18)
    plt.show()
    '''

    plt.text(0.5, 1.08,'UT='+str(selected_ut) +'  LAT='+str(selected_lat)+" SLT="+str(longitude_to_local_time(selected_lon))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.25, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    #plt.text(0.5, -0.3, str(dataset.split("/")[-1]), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    return(plot)


def plt_lev_lon(dataset, variable_name, selected_time, selected_lat):
    """
    Generates a contour plot for the given 2D array of variable values, latitude, and longitude.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lon, ilev dimensions.
            - Valid variables:['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 
                                'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 
                                'Z', 'POTEN']le
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_ilev (float): The selected ilevel value.
    
    Returns:
        - Contour plot.
    """
    # Printing Execution data
    
    print("---------------["+variable_name+"]---["+str(selected_time)+"]---["+str(selected_lat)+"]---------------")
    # Generate 2D arrays, extract variable_unit
    variable_values, unique_lons, unique_levs,selected_lat, variable_unit, variable_long_name, selected_ut, selected_mtime = lev_ilev_lon(dataset, variable_name, selected_time, selected_lat)

    
    min_val, max_val = min_max(variable_values)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]

    cmap_color, contour_color = color_scheme(variable_name)

    
    # Generate contour plot
    plot=plt.figure(figsize=(24, 12))
    contour_filled = plt.contourf(unique_lons, unique_levs, variable_values, cmap= cmap_color, levels=20)
    contour_lines = plt.contour(unique_lons, unique_levs, variable_values, colors=contour_color, linewidths=0.5, levels=20)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=contour_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    plt.text(0, 9,'UT='+str(selected_ut) +'  LAT='+str(selected_lat), ha='center', va='center',fontsize=28) 
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    plt.xlabel('Longitude (Deg)',fontsize=28)
    plt.xticks([value for value in unique_lons if value % 30 == 0],fontsize=18)  
    plt.yticks(fontsize=18) 

    # Add Local Time secondary x-axis
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2_xticks = ax.get_xticks()
    ax2.set_xticks(ax2_xticks)
    ax2.set_xticklabels([str(int(longitude_to_local_time(lon) % 24)) for lon in ax2_xticks],fontsize=18)
    ax2.set_xlabel('Local Time (Hrs)', labelpad=15, fontsize=28)

    # Add subtext to the plot
    plt.text(-90, -9, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28)
    plt.text(90, -9, "Contour Interval = "+str("{:.2e}".format((max_val-min_val)/20)), ha='center', va='center',fontsize=28)
    plt.text(90, -10, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28)
    #plt.text(-90, -10, str(dataset.split("/")[-1]), ha='center', va='center',fontsize=28)

    plt.show()
    #plot, ax = plt.subplots()

    
    
    return(plot)


def plt_lev_lat(dataset, variable_name, selected_time, selected_lon):
    """
    Generates a contour plot for the given 2D array of variable values, latitude, and latgitude.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lat, ilev dimensions.
            - Valid variables:['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 
                                'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 
                                'Z', 'POTEN']
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_ilev (float): The selected ilevel value.
    
    Returns:
        - Contour plot.
    """
    # Printing Execution data
    
    print("---------------["+variable_name+"]---["+str(selected_time)+"]---["+str(selected_lon)+"]---------------")
    # Generate 2D arrays, extract variable_unit
    variable_values, unique_lats, unique_levs,selected_lon, variable_unit, variable_long_name, selected_ut, selected_mtime = lev_ilev_lat(dataset, variable_name, selected_time, selected_lon)

    
    min_val, max_val = min_max(variable_values)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]

    cmap_color, contour_color = color_scheme(variable_name)

    
    # Generate contour plot
    plot=plt.figure(figsize=(24, 12))
    contour_filled = plt.contourf(unique_lats, unique_levs, variable_values, cmap= cmap_color, levels=20)
    contour_lines = plt.contour(unique_lats, unique_levs, variable_values, colors=contour_color, linewidths=0.5, levels=20)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=contour_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    plt.text(0, 8,'UT='+str(selected_ut) +'  LON='+str(selected_lon)+" SLT="+str(longitude_to_local_time(selected_lon))+"Hrs", ha='center', va='center',fontsize=28) 
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    plt.xlabel('Latitude (Deg)',fontsize=28)
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 



    # Add subtext to the plot
    plt.text(-45, -9, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28)
    plt.text(45, -9, "Contour Interval = "+str("{:.2e}".format((max_val-min_val)/20)), ha='center', va='center',fontsize=28)
    plt.text(45, -10, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=28)
    #plt.text(-90, -10, str(dataset.split("/")[-1]), ha='center', va='center',fontsize=28)

    plt.show()
    #plot, ax = plt.subplots()

    
    
    return(plot)





def plt_lev_time(datasets, variable_name, selected_lat, selected_lon):
    """
    Generates a contour plot for the given 2D array of variable values, latitude, and latgitude.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lat, ilev dimensions.
            - Valid variables:['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 
                                'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 
                                'Z', 'POTEN']
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_ilev (float): The selected ilevel value.
    
    Returns:
        - Contour plot.
    """

    variable_values_all, levs_ilevs, mtime_values, selected_lon, variable_unit, variable_long_name = lev_ilev_time(datasets, variable_name, selected_lat, selected_lon)
    
    # Assuming the levels are consistent across datasets, but using the minimum size for safety
    

    print("---------------["+variable_name+"]---["+str(selected_lat)+"]---["+str(selected_lon)+"]---------------")
    
    min_val, max_val = np.nanmin(variable_values_all), np.nanmax(variable_values_all)

    cmap_color, contour_color = color_scheme(variable_name)

    unique_days = sorted(list(set([day for day, _, _ in mtime_values])))
    day_indices = [i for i, (day, _, _) in enumerate(mtime_values) if i == 0 or mtime_values[i-1][0] != day]

    plot=plt.figure(figsize=(24, 12))
    X, Y = np.meshgrid(range(len(mtime_values)), levs_ilevs)
    contour_filled = plt.contourf(X, Y, variable_values_all, cmap=cmap_color, levels=20)
    contour_lines = plt.contour(X, Y, variable_values_all, colors=contour_color, linewidths=0.5, levels=20)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=contour_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.xticks(day_indices, unique_days, rotation=45)
    plt.xlabel("Model Time (Day) from "+str(np.nanmin(unique_days))+" to "+str(np.nanmax(unique_days)) ,fontsize=28)
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    
    plt.tight_layout()
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 



    # Add subtext to the plot
    plt.text(0.5, 1.08,'  LAT='+str(selected_lat)+" SLT="+str(longitude_to_local_time(selected_lon))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.5, -0.25, "Contour Interval = "+str("{:.2e}".format((max_val-min_val)/20)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    return(plot)



def plt_lev_time(datasets, variable_name, selected_lat, selected_lon):
    """
    Generates a contour plot for the given 2D array of variable values, latitude, and latgitude.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lat, ilev dimensions.
            - Valid variables:['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 
                                'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 
                                'Z', 'POTEN']
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_ilev (float): The selected ilevel value.
    
    Returns:
        - Contour plot.
    """

    variable_values_all, levs_ilevs, mtime_values, selected_lon, variable_unit, variable_long_name = lev_ilev_time(datasets, variable_name, selected_lat, selected_lon)
    
    # Assuming the levels are consistent across datasets, but using the minimum size for safety
    

    print("---------------["+variable_name+"]---["+str(selected_lat)+"]---["+str(selected_lon)+"]---------------")
    
    min_val, max_val = np.nanmin(variable_values_all), np.nanmax(variable_values_all)

    cmap_color, contour_color = color_scheme(variable_name)

    unique_days = sorted(list(set([day for day, _, _ in mtime_values])))
    day_indices = [i for i, (day, _, _) in enumerate(mtime_values) if i == 0 or mtime_values[i-1][0] != day]

    plot=plt.figure(figsize=(20, 12))
    X, Y = np.meshgrid(range(len(mtime_values)), levs_ilevs)
    contour_filled = plt.contourf(X, Y, variable_values_all, cmap=cmap_color, levels=20)
    contour_lines = plt.contour(X, Y, variable_values_all, colors=contour_color, linewidths=0.5, levels=20)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=contour_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.xticks(day_indices, unique_days, rotation=45)
    plt.xlabel("Model Time (Day) from "+str(np.nanmin(unique_days))+" to "+str(np.nanmax(unique_days)) ,fontsize=28)
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    
    plt.tight_layout()
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 



    # Add subtext to the plot
    plt.text(0.5, 1.08,'  LAT='+str(selected_lat)+" SLT="+str(longitude_to_local_time(selected_lon))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.5, -0.25, "Contour Interval = "+str("{:.2e}".format((max_val-min_val)/20)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    return(plot)



def plt_lat_time(datasets, variable_name, selected_lev, selected_lon):
    """
    Generates a contour plot for the given 2D array of variable values, latitude, and latgitude.
    
    Parameters:
        - dataset (str): Path to the NetCDF file.
        - variable_name (str): The name of the variable with lat, lat, ilev dimensions.
            - Valid variables:['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D','TI', 'TE', 'O2P', 'TN_NM', 
                                'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 
                                'Z', 'POTEN']
        - selected_time (str): The selected datetime in the format 'YYYY-MM-DDTHH:MM:SS'.
        - selected_ilev (float): The selected ilevel value.
    
    Returns:
        - Contour plot.
    """
    print("---------------["+variable_name+"]---["+str(selected_lev)+"]---["+str(selected_lon)+"]---------------")

    try:
        variable_values_all, lats, mtime_values, selected_lon, variable_unit, variable_long_name = lat_time_lev(datasets, variable_name, selected_lev, selected_lon)
    except:
        variable_values_all, lats, mtime_values, selected_lon, variable_unit, variable_long_name = lat_time_ilev(datasets, variable_name, selected_lev, selected_lon)
    # Assuming the levels are consistent across datasets, but using the minimum size for safety
    

    
    
    min_val, max_val = np.nanmin(variable_values_all), np.nanmax(variable_values_all)

    cmap_color, contour_color = color_scheme(variable_name)

    unique_days = sorted(list(set([day for day, _, _ in mtime_values])))
    day_indices = [i for i, (day, _, _) in enumerate(mtime_values) if i == 0 or mtime_values[i-1][0] != day]

    plot=plt.figure(figsize=(20, 12))
    X, Y = np.meshgrid(range(len(mtime_values)), lats)
    contour_filled = plt.contourf(X, Y, variable_values_all, cmap=cmap_color, levels=20)
    contour_lines = plt.contour(X, Y, variable_values_all, colors=contour_color, linewidths=0.5, levels=20)
    plt.clabel(contour_lines, inline=True, fontsize=16, colors=contour_color)
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=28, labelpad=15)
    cbar.ax.tick_params(labelsize=18)
    plt.xticks(day_indices, unique_days, rotation=45)
    plt.xlabel("Model Time (Day) from "+str(np.nanmin(unique_days))+" to "+str(np.nanmax(unique_days)) ,fontsize=28)
    plt.ylabel('LN(P0/P) (INTERFACES)',fontsize=28)
    
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=36 )   
    
    plt.tight_layout()
    plt.xticks(fontsize=18)  
    plt.yticks(fontsize=18) 



    # Add subtext to the plot
    plt.text(0.5, 1.08,'  ZP='+str(selected_lev)+" SLT="+str(longitude_to_local_time(selected_lon))+"Hrs", ha='center', va='center',fontsize=28, transform=plt.gca().transAxes) 
    plt.text(0.5, -0.2, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)
    plt.text(0.5, -0.25, "Contour Interval = "+str("{:.2e}".format((max_val-min_val)/20)), ha='center', va='center',fontsize=28, transform=plt.gca().transAxes)

    return(plot)
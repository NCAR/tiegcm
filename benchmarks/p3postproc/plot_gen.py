
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from data_parse import lat_lon_lev, lat_lon_ilev,calc_avg_ht, min_max

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
        data_array, selected_lev_ilev, variable_unit, variable_long_name, selected_ut, selected_mtime =lat_lon_lev(dataset, variable_name, selected_time, selected_lev_ilev)
    except ValueError:
        data_array, selected_lev_ilev, variable_unit, variable_long_name, selected_ut, selected_mtime=lat_lon_ilev(dataset, variable_name, selected_time, selected_lev_ilev)

    #print(data_array)
    avg_ht=calc_avg_ht(data_array, selected_time,selected_lev_ilev, dataset)
    min_val, max_val = min_max(data_array)
    selected_day=selected_mtime[0]
    selected_hour=selected_mtime[1]
    selected_min=selected_mtime[2]
    # Extract values, latitudes, and longitudes from the array
    values = [row[0] for row in data_array]
    lats = [row[1] for row in data_array]
    lons = [row[2] for row in data_array]
    
    # Convert lists to 2D arrays for plotting
    unique_lats = sorted(list(set(lats)))
    unique_lons = sorted(list(set(lons)))
    values_2d = np.array(values).reshape(len(unique_lats), len(unique_lons))
    
    # Generate contour plot
    plot=plt.figure(figsize=(24, 12))
    contour_filled = plt.contourf(unique_lons, unique_lats, values_2d, cmap='viridis', levels=20)
    contour_lines = plt.contour(unique_lons, unique_lats, values_2d, colors='black', linewidths=0.5, levels=20)
    plt.clabel(contour_lines, inline=True, fontsize=12, colors='black')
    cbar = plt.colorbar(contour_filled, label=variable_name+" ["+variable_unit+"]")
    cbar.set_label(variable_name+" ["+variable_unit+"]", size=24)
    cbar.ax.tick_params(labelsize=12)
    plt.title(variable_long_name+' '+variable_name+' ('+variable_unit+') '+'\n\n',fontsize=32 )   
    plt.text(0, 110,'UT='+str(selected_ut) +'  ZP='+str(selected_lev_ilev)+" Avg HT="+str(avg_ht), ha='center', va='center',fontsize=24) #selected_lev_ilev is not the used lev_ilev fix this
    plt.xlabel('Longitude (Deg)',fontsize=24)
    plt.ylabel('Latitude (Deg)',fontsize=24)
    plt.xticks([value for value in unique_lons if value % 30 == 0],fontsize=12)  
    plt.yticks(fontsize=12) 

    # Add Local Time secondary x-axis
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2_xticks = ax.get_xticks()
    ax2.set_xticks(ax2_xticks)
    ax2.set_xticklabels([str(int(longitude_to_local_time(lon) % 24)) for lon in ax2_xticks],fontsize=12)
    ax2.set_xlabel('Local Time (Hrs)',fontsize=24)

    # Add subtext to the plot
    plt.text(-90, -110, "Min, Max = "+str("{:.2e}".format(min_val))+", "+str("{:.2e}".format(max_val)), ha='center', va='center',fontsize=24)
    plt.text(90, -110, "Contour Interval = "+str("{:.2e}".format((max_val-min_val)/20)), ha='center', va='center',fontsize=24)
    plt.text(90, -120, "Day, Hour, Min = "+str(selected_day)+","+str(selected_hour)+","+str(selected_min), ha='center', va='center',fontsize=24)
    plt.text(-90, -120, str(dataset.split("/")[-1]), ha='center', va='center',fontsize=24)

    plt.show()
    #plot, ax = plt.subplots()

    
    
    return(plot)


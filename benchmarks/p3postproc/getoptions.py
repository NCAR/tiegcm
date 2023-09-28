import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Generate different types of plots based on user input.')
    
    parser.add_argument('--plot_type', type=str, required=True, choices=['lat_lon', 'lev_var', 'lev_lon', 'lev_lat', 'lev_time', 'lat_time'], help='Type of the plot to be generated.')
    parser.add_argument('-var','--variable_name', type=str, required=True, help='Name of the variable to be plotted.')
    parser.add_argument('-time','--time', type=str, help='Selected time for the plot in YYYY-MM-DDTHH:MM:SS format.')
    parser.add_argument('-mtime','--mtime',  nargs=3, type=int, help='Selected time for the plot in [Day, Hour, Min] format.')
    parser.add_argument('-lev','--level', type=float, help='Selected lev/ilev for the plot.')
    parser.add_argument('-lat', '--latitude', type=float, help='Selected latitude for the plot.')
    parser.add_argument('-lon','--longitude', type=float, help='Selected longitude for the plot.')
    parser.add_argument('-ut','--localtime', type=float, help='Selected localtime / longitude for the plot.')

    parser.add_argument('-dir','--directory', type=str, help='Directory path containing the datasets.')
    parser.add_argument('--dataset_filter', type=str, help='Filter to load datasets.')
    parser.add_argument('--output_format', type=str, choices=['jpeg', 'pdf'], help='Format to save the plots.')
    parser.add_argument('--coastlines', type=bool, help='Add coast lines to the lat_lon plots')

    parser.add_argument('-zpmin','--level_minimum', type=float, help='Minimum level to plots [lev_var, lev_lon, lev_lat, lev_time]')
    parser.add_argument('-zpmax','--level_maximum', type=float, help='Maximum level to plots [lev_var, lev_lon, lev_lat, lev_time]')

    parser.add_argument('-latmin','--latitude_minimum', type=float, help='Minimum latitude to plots [lat_lon, lev_lat, lat_time]')
    parser.add_argument('-latmax','--latitude_maximum', type=float, help='Maximum latitude to plots [lat_lon, lev_lat, lat_time]')

    parser.add_argument('-lonmin','--longitude_minimum', type=float, help='Minimum longitude to plots [lat_lon, lev_lon]')
    parser.add_argument('-lonmax','--longitude_maximum', type=float, help='Maximum longitude to plots [lat_lon, lev_lon]')
    parser.add_argument('-utmin','--localtime_minimum', type=float, help='Minimum localtime to plots [lat_lon, lev_lon]')
    parser.add_argument('-utmax','--localtime_maximum', type=float, help='Maximum localtime to plots [lat_lon, lev_lon]')

    args = parser.parse_args()
    
    plot_requirements = {
        'lat_lon': {'required': ['variable_name', 'time', 'mtime', 'level'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_levels': [-7.  , -6.75, -6.5 , -6.25, -6.  , -5.75, -5.5 , -5.25, -5.  , -4.75, -4.5 , -4.25, -4.  , -3.75, -3.5 , -3.25, -3.  , -2.75, -2.5 , -2.25, -2.  , -1.75, -1.5 , -1.25, -1.  , -0.75, -0.5 , -0.25,  0.  ,  0.25,  0.5 ,  0.75,  1.  ,  1.25,  1.5 ,  1.75, 2.  ,  2.25,  2.5 ,  2.75,  3.  ,  3.25,  3.5 ,  3.75,  4.  , 4.25,  4.5 ,  4.75,  5.  ,  5.25,  5.5 ,  5.75,  6.  ,  6.25, 6.5 ,  6.75,  7.  ,  7.25]
                    },
        
        'lev_var': {'required': ['variable_name', 'time', 'mtime', 'latitude', 'longitude', 'localtime'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_latitudes': [-87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5,  42.5, 47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5,  82.5,  87.5],
                    'valid_longitudes': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.]
                    },
        'lev_lon': {'required': ['variable_name', 'time', 'mtime', 'latitude'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_latitudes': [-87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5,  42.5, 47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5,  82.5,  87.5]
                    },
        'lev_lat': {'required': ['variable_name', 'time', 'mtime', 'longitude', 'localtime'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_longitude': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.]
                    },
        'lev_time': {'required': ['variable_name', 'latitude', 'longitude', 'localtime'], 
                     'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                     'valid_latitudes': [-87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5,  42.5, 47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5,  82.5,  87.5],
                     'valid_longitudes': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.]
                     },
        'lat_time': {'required': ['variable_name', 'level', 'longitude', 'localtime'], 
                     'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                     'valid_longitudes': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.],
                     'valid_levels': [-7.  , -6.75, -6.5 , -6.25, -6.  , -5.75, -5.5 , -5.25, -5.  , -4.75, -4.5 , -4.25, -4.  , -3.75, -3.5 , -3.25, -3.  , -2.75, -2.5 , -2.25, -2.  , -1.75, -1.5 , -1.25, -1.  , -0.75, -0.5 , -0.25,  0.  ,  0.25,  0.5 ,  0.75,  1.  ,  1.25,  1.5 ,  1.75, 2.  ,  2.25,  2.5 ,  2.75,  3.  ,  3.25,  3.5 ,  3.75,  4.  , 4.25,  4.5 ,  4.75,  5.  ,  5.25,  5.5 ,  5.75,  6.  ,  6.25, 6.5 ,  6.75,  7.  ,  7.25]
                     },
    }
    

    if args.plot_type in plot_requirements:
        requirements = plot_requirements[args.plot_type]
        for arg in requirements['required']:
            if arg in ['time', 'mtime', 'longitude', 'localtime']:  # Skip checking for 'time', 'mtime', 'longitude', and 'localtime' here
                continue
            if getattr(args, arg) is None:
                parser.error(f"{args.plot_type} requires the argument --{arg}")
        
        # Check if 'longitude' or 'localtime' is in requirements and at least one of them is provided
        if 'longitude' in requirements['required'] or 'localtime' in requirements['required']:
            if getattr(args, 'longitude') is None and getattr(args, 'localtime') is None:
                parser.error(f"{args.plot_type} requires either the argument --longitude or --localtime")
        
        # Enforce that either 'time' or 'mtime' must be provided
        if 'time' in requirements['required'] or 'mtime' in requirements['required']:
            if getattr(args, 'time') is None and getattr(args, 'mtime') is None:
                parser.error(f"{args.plot_type} requires either the argument --time or --mtime")
        
        if args.variable_name not in requirements['valid_variables']:
            parser.error(f"Invalid variable_name for {args.plot_type}. Valid variables are {', '.join(requirements['valid_variables'])}.")

    level_bound_plots = ['lev_var', 'lev_lon', 'lev_lat', 'lev_time']
    if args.plot_type not in level_bound_plots:
        if args.level_minimum is not None or args.level_maximum is not None:
            parser.error(f"--level_minimum and --level_maximum can only be used with plot types: {', '.join(level_bound_plots)}")     

    longitude_bound_plots = ['lat_lon', 'lev_lon']
    if args.plot_type not in longitude_bound_plots:
        if args.longitude_minimum is not None or args.longitude_maximum is not None:
            parser.error(f"--longitude_minimum and --longitude_maximum can only be used with plot types: {', '.join(longitude_bound_plots)}") 

    latitude_bound_plots = ['lat_lon', 'lev_lat', 'lat_time']
    if args.plot_type not in latitude_bound_plots:
        if args.latitude_minimum is not None or args.latitude_maximum is not None:
            parser.error(f"--latitude_minimum and --latitude_maximum can only be used with plot types: {', '.join(latitude_bound_plots)}") 

    return args

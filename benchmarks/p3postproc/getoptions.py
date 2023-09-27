import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Generate different types of plots based on user input.')
    
    parser.add_argument('--plot_type', type=str, required=True, choices=['lat_lon', 'lev_var', 'lev_lon', 'lev_lat', 'lev_time', 'lat_time'], help='Type of the plot to be generated.')
    parser.add_argument('-var','--variable_name', type=str, required=True, help='Name of the variable to be plotted.')
    parser.add_argument('-time','--time', type=str, help='Selected time for the plot.')
    parser.add_argument('-mtime','--mtime',  nargs=3, type=int, help='Selected time for the plot.')
    parser.add_argument('-lev','--level', type=float, help='Selected lev/ilev for the plot.')
    parser.add_argument('-lat', '--latitude', type=float, help='Selected latitude for the plot.')
    parser.add_argument('-lon','--longitude', type=float, help='Selected longitude for the plot.')
    parser.add_argument('--dir', type=str, help='Directory path containing the datasets.')
    parser.add_argument('--dataset_filter', type=str, help='Filter to load datasets.')
    parser.add_argument('--output_format', type=str, choices=['jpeg', 'pdf'], help='Format to save the plots.')
    parser.add_argument('--coastlines', type=bool, help='Add coast lines to the plots.')
    parser.add_argument('-levmin','--level_minimum', type=float, help='Minimum level to plot')
    parser.add_argument('-levmax','--level_maximum', type=float, help='Maximum level to plot')
    args = parser.parse_args()
    
    plot_requirements = {
        'lat_lon': {'required': ['variable_name', 'time', 'mtime', 'level'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_levels': [-7.  , -6.75, -6.5 , -6.25, -6.  , -5.75, -5.5 , -5.25, -5.  , -4.75, -4.5 , -4.25, -4.  , -3.75, -3.5 , -3.25, -3.  , -2.75, -2.5 , -2.25, -2.  , -1.75, -1.5 , -1.25, -1.  , -0.75, -0.5 , -0.25,  0.  ,  0.25,  0.5 ,  0.75,  1.  ,  1.25,  1.5 ,  1.75, 2.  ,  2.25,  2.5 ,  2.75,  3.  ,  3.25,  3.5 ,  3.75,  4.  , 4.25,  4.5 ,  4.75,  5.  ,  5.25,  5.5 ,  5.75,  6.  ,  6.25, 6.5 ,  6.75,  7.  ,  7.25]
                    },
        
        'lev_var': {'required': ['variable_name', 'time', 'mtime', 'latitude', 'longitude'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_latitudes': [-87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5,  42.5, 47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5,  82.5,  87.5],
                    'valid_longitudes': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.]
                    },
        'lev_lon': {'required': ['variable_name', 'time', 'mtime', 'latitude'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_latitudes': [-87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5,  42.5, 47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5,  82.5,  87.5]
                    },
        'lev_lat': {'required': ['variable_name', 'time', 'mtime', 'longitude'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                    'valid_longitude': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.]
                    },
        'lev_time': {'required': ['variable_name', 'latitude', 'longitude'], 
                     'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                     'valid_latitudes': [-87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, -37.5, -32.5, -27.5, -22.5, -17.5, -12.5,  -7.5,  -2.5, 2.5, 7.5,  12.5,  17.5,  22.5,  27.5,  32.5,  37.5,  42.5, 47.5,  52.5,  57.5,  62.5,  67.5,  72.5,  77.5,  82.5,  87.5],
                     'valid_longitudes': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.]
                     },
        'lat_time': {'required': ['variable_name', 'level', 'longitude'], 
                     'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN'],
                     'valid_longitudes': [-180., -175., -170., -165., -160., -155., -150., -145., -140., -135., -130., -125., -120., -115., -110., -105., -100.,  -95., -90.,  -85.,  -80.,  -75.,  -70.,  -65.,  -60.,  -55.,  -50., -45.,  -40.,  -35.,  -30.,  -25.,  -20.,  -15.,  -10.,   -5.,  0.,    5.,   10.,   15.,   20.,   25.,   30.,   35.,   40., 45.,   50.,   55.,   60.,   65.,   70.,   75.,   80.,   85., 90.,   95.,  100.,  105.,  110.,  115.,  120.,  125.,  130., 135.,  140.,  145.,  150.,  155.,  160.,  165.,  170.,  175.],
                     'valid_levels': [-7.  , -6.75, -6.5 , -6.25, -6.  , -5.75, -5.5 , -5.25, -5.  , -4.75, -4.5 , -4.25, -4.  , -3.75, -3.5 , -3.25, -3.  , -2.75, -2.5 , -2.25, -2.  , -1.75, -1.5 , -1.25, -1.  , -0.75, -0.5 , -0.25,  0.  ,  0.25,  0.5 ,  0.75,  1.  ,  1.25,  1.5 ,  1.75, 2.  ,  2.25,  2.5 ,  2.75,  3.  ,  3.25,  3.5 ,  3.75,  4.  , 4.25,  4.5 ,  4.75,  5.  ,  5.25,  5.5 ,  5.75,  6.  ,  6.25, 6.5 ,  6.75,  7.  ,  7.25]
                     },
    }
    

    if args.plot_type in plot_requirements:
        requirements = plot_requirements[args.plot_type]
        for arg in requirements['required']:
            if arg == 'time' or arg == 'mtime':  # Skip checking for 'time' and 'mtime' here
                continue
            if getattr(args, arg) is None:
                parser.error(f"{args.plot_type} requires the argument --{arg}")
        
        # Enforce that either 'time' or 'mtime' must be provided
        if getattr(args, 'time') is None and getattr(args, 'mtime') is None:
            parser.error(f"{args.plot_type} requires either the argument --time or --mtime")
        
        if args.variable_name not in requirements['valid_variables']:
            parser.error(f"Invalid variable_name for {args.plot_type}. Valid variables are {', '.join(requirements['valid_variables'])}.")
            
    return args

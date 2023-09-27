import argparse

def get_options():
    parser = argparse.ArgumentParser(description='Generate different types of plots based on user input.')
    
    parser.add_argument('--plot_type', type=str, required=True, help='Type of the plot to be generated.')
    parser.add_argument('--dataset', type=str, help='Path to the dataset.')
    parser.add_argument('--variable_name', type=str, required=True, help='Name of the variable to be plotted.')
    parser.add_argument('--time', type=str, help='Selected time for the plot.')
    parser.add_argument('-lev','level', type=str, help='Selected lev/ilev for the plot.')
    parser.add_argument('-lat', '--latitude', type=str, help='Selected latitude for the plot.')
    parser.add_argument('-lon','--longitude', type=str, help='Selected longitude for the plot.')
    parser.add_argument('--dir', type=str, help='Directory path containing the datasets.')
    parser.add_argument('--dataset_filter', type=str, help='Filter to load datasets.')
    parser.add_argument('--output_format', type=str, choices=['jpeg', 'pdf'], help='Format to save the plots.')
    
    args = parser.parse_args()
    
    plot_requirements = {
        'lat_lon': {'required': ['variable_name', 'time', 'level'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN']},
        'lev_var': {'required': ['variable_name', 'time', 'latitude', 'longitude'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN']},
        'lev_lon': {'required': ['variable_name', 'time', 'latitude'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN']},
        'lev_lat': {'required': ['variable_name', 'time', 'longitude'], 
                    'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN']},
        'lev_time': {'required': ['variable_name', 'latitude', 'longitude'], 
                     'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN']},
        'lat_time': {'required': ['variable_name', 'level', 'longitude'], 
                     'valid_variables': ['TN', 'UN', 'VN', 'O2', 'O1', 'N4S', 'NO', 'HE', 'AR', 'OP', 'N2D', 'TI', 'TE', 'O2P', 'TN_NM', 'UN_NM', 'VN_NM', 'O2_NM', 'O1_NM', 'N4S_NM', 'NO_NM', 'OP_NM', 'HE_NM', 'AR_NM', 'NE', 'OMEGA', 'Z', 'POTEN']},
    }
    
    if args.plot_type in plot_requirements:
        requirements = plot_requirements[args.plot_type]
        for arg in requirements['required']:
            if getattr(args, arg) is None:
                parser.error(f"{args.plot_type} requires the argument --{arg}")
        
        if args.variable_name not in requirements['valid_variables']:
            parser.error(f"Invalid variable_name for {args.plot_type}. Valid variables are {', '.join(requirements['valid_variables'])}.")
    
    return args

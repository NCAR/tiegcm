conversion_units = {
    'cm-3': {'m-3': 1000000.0},
    'km': {'m': 1000.0},
    'cm/s': {'m/s': 0.01, 'km/s': 1e-5},
    'erg/g/s': {'J/kg/s': 1e-7},
    'ergs/cm2/s': {'J/m2/s': 0.001},
    'millibars': {'pascals': 100},
    'microbars': {'pascals': 0.1},
    'g/cm3': {'kg/m3': 1000.0},
    'degrees': {'radians': 0.0174533},
    'degrees_east': {'radians': 0.0174533},
    'degrees_north': {'radians': 0.0174533},
    'km/s': {'m/s': 1000.0, 'cm/s': 100000.0},
    'GW': {'MW': 1000.0, 'kW': 1000000.0, 'W': 1000000000.0},
    'keV': {'eV': 1000.0, 'J': 1.60218e-16},
    'nT': {'µT': 0.001, 'T': 1e-9},
    'cm': {'m': 0.01, 'km': 1e-5},
    'K': {
        'C': {'factor': 1, 'offset': -273.15},  # Celsius = Kelvin - 273.15
        'F': {'factor': 9/5, 'offset': -459.67}  # Fahrenheit = Kelvin * 9/5 - 459.67
    }
}

def convert_units(data, from_unit, to_unit):
    """
    Convert data from one unit to another based on predefined conversion factors.
    
    Parameters:
    - data: Numeric data to be converted
    - from_unit: The current unit of the data
    - to_unit: The desired unit to convert to
    
    Returns:
    - Converted data
    """
    
    if from_unit == to_unit:
        return data, from_unit
    
    # Check if conversion is possible
    if from_unit in conversion_units and to_unit in conversion_units[from_unit]:
        conversion = conversion_units[from_unit][to_unit]
        print(f"Converting from {from_unit} to {to_unit}")

        # Check if conversion is a dictionary or direct float
        if isinstance(conversion, dict):
            factor = conversion.get('factor', 1)
            offset = conversion.get('offset', 0)

            return data * factor + offset, to_unit
        else:
            return data * conversion, to_unit
    else:
        #raise ValueError(f"No conversion factor found for {from_unit} to {to_unit}")
        print(f"No conversion factor found for {from_unit} to {to_unit}")
        return data, from_unit

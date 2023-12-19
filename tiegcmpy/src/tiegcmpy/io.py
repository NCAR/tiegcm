import os
import sys  
import inspect
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np



def load_datasets(directory,dataset_filter = None):
    datasets=[]
    if os.path.isdir(directory):
        files = sorted(os.listdir(directory)) 
        print("Loading datasets globally.") 
        for file in files:
            if file.endswith('.nc') and (dataset_filter is None or dataset_filter in file):
                file_path = os.path.join(directory, file)
                datasets.append([xr.open_dataset(file_path), file])
    else:
        file = os.path.basename(directory)
        datasets.append([xr.open_dataset(directory), file])
    return(datasets)
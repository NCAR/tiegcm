import os
import sys  
import inspect
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np


def load_datasets(directory,dataset_filter):
    datasets=[]
    files = sorted(os.listdir(directory)) 
    print("Loading datasets globally.") 
    for file in files:
        if file.endswith('.nc') and (dataset_filter is None or dataset_filter in file):
            file_path = os.path.join(directory, file)
            datasets.append([xr.open_dataset(file_path), file])
    return(datasets)

def load_dataset(file_path):
    datasets=[]
    #directory = os.path.dirname(file_path)
    file = os.path.basename(file_path)
    datasets.append([xr.open_dataset(file_path), file])
    return(datasets)
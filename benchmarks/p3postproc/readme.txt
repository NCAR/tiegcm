For p3postproc here are the requirements:
    -module load conda
    -module load geos for cartopy (Not yet used)
    -conda create  --name {env_name} --clone /glade/work/nikhilr/conda-envs/p3postproc
    -conda activate {env_name}
    - chmod +x * on the directory
    - ./tiepy with the arguments  
        - Example: ./tiepy --plot_type lev_lat -dir /glade/u/home/nikhilr/tiegcm_func/postproc/tiegcm_res5.0_decsol_smax/hist/ --dataset_filter sech -var TN -mtime 360 0 0 -lon 0.0 --output_format jpeg


cond env contains 
    -Python 3.8
    -netcdf4
    -cartopy
    -numpy
    -xarray

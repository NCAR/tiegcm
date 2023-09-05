For p3postproc here are the requirements:
    -module load conda
    -module load geos for cartopy (Not yet used)
    -conda create  --name {env_name} --clone /glade/work/nikhilr/conda-envs/p3postproc
    -conda activate {env_name}


cond env contains 
    -Python 3.8
    -netcdf4
    -cartopy
    -numpy

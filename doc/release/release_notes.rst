
Release Notes for |modeluc| |version|
=====================================

:Release date: March, 2016
:Contact: discussion group email list `tgcmgroup@ucar.edu <http://mailman.ucar.edu/mailman/listinfo/tgcmgroup>`_
:Download page: http://www.hao.ucar.edu/modeling/tgcm/download.php

This is a summary of modifications made to the |model| since the release of 
|model|\1.95 (June 21, 2013)

New Features
------------

.. describe:: 2.5-degree resolution

   Dynamical filters have been tuned for 2.5-degree horizontal resolution 
   (with four grid-points per scale height in the vertical coordinate). 
   Users are encouraged to run the model at the 2.5-degree resolution, as 
   it can produce a more realistic ionosphere than the 5-degree model. 
   However, 2.5-degree resolution is not recommended for earlier versions 
   of the model. All benchmark runs and source history files are available 
   at both 5-degree and 2.5-degree resolutions.

.. describe:: Helium

   Helium is now calculated as a major species. If a source history without 
   helium is used, helium will be initialized to a global value (1.154e-6 mmr), 
   and will evolve from there. Helium is always saved on primary histories, 
   and can be saved on secondary histories as ‘HE’. The namelist input flag 
   CALC_HELIUM can toggle this calculation on or off, but it is on by default. 
   Thanks to Eric Sutton for leading this effort. See Sutton et al. (2015), 
   JGR, 120, 6884, doi:10.1002/2015JA021223 for details of the implementation 
   and results.

.. describe:: Argon

   Argon is now calculated as a minor species. Argon is always saved on primary 
   histories, and can be saved on secondary histories as ‘AR’.

.. describe:: Parallel Dynamo

   The electro-dynamo code was parallelized with pure MPI. Transformations between 
   the geographic and geomagnetic grids in the MPI environment are accomplished 
   with the Earth System Modeling Framework (ESMF library). This results in a 
   performance speed-up of about 25% at the 5-degree resolution, and almost 50% 
   at the 2.5-degree resolution when using 64-processors on the NCAR |ncarsuper| system.

.. describe:: IGRF12

   The geomagnetic field was updated to the International Geomagnetic Reference 
   Field version 12, and annual secular variation is included for the years 1900-2015. 
   Extrapolation of secular variation through the year 2020 is performed.

.. describe:: Lower Boundary Condition

   A zonal mean climatology of u, v, T, and z was implemented at the lower boundary 
   by Astrid Maute, following work by Mac Jones using a combination of MSIS and HWM 
   empirical models, and UARS data. See Jones et al. (2014), JGR, 119, 2197, 
   doi:10.1002/2013JA019744. The climatology is based on monthly means, and is applied 
   by specifying an input file (similar to GSWM and other tidal specifications). 
   If no input file is specified, a flat lower boundary (u=v=0, Tn=181 K, z=96.4 km) 
   is employed as in the past. Other zonal mean climatologies can be used by generating 
   and specifying a different file.

.. describe:: Non-Migrating Tides

   GSWM non-migrating tides are turned on in the default namelist inputs for 2.5-degree 
   resolution only. They can be used with either resolution but are only recommended at 
   2.5-degree. Migrating tides are on be default at both resolutions. All GSWM tides are 
   switchable on/off simply by specifying or not specifying an input file.

.. describe:: Additional Diagnostic Fields

   Many additional diagnostic fields are available as optional outputs on secondary files, 
   including RHO, N2, HE, AR, ZGMID, EFLUX, NFLUX, ALFA, CUSP, and DRIZZLE, TEC, 
   conductivities, ExB velocities, currents, cooling and heating rates. See the Namelist 
   Input File section of the User Guide, scripts/master.inp, and src/diags.F.

.. describe:: Pressure/Altitude Coordinates

   For altitude-registering of output fields, see
   :base_url:`Grid Structure and Resolution <userguide/html/grid.html>`
   section of the User Guide for an explanation of the pressure grids, reference pressure, 
   geometric v. geopotential height, and interfaces v. midpoints. 
   See also :base_url:`Altitude Coordinates in the TIEGCM <userguide/html/grid.html#altitude-coordinates-in-the-ncar-tie-gcm-and-time-gcm>`
   for a detailed explanation of the relationship between the Pressure and Altitude
   coordinates in the TIEGCM. 

Functional Changes
------------------

.. describe:: Namelist Input

   The comment character in namelist input files is now and exclamation point (!) 
   (formerly it was a semi-colon). This was changed to conform to the Fortran standard. 
   If you have many input files with the ‘;’ comment character, use the change_nlcomment 
   script in the scripts directory to replace the comment chars. However, we strongly 
   urge that for v. 2.0, users start over with new namelist input files using the examples 
   in scripts/*.inp as templates. Also use the new example job files in scripts/*.job. 
   Reading the input file was simplified in the source code (see src/input_read.F), and 
   in the job scripts, the input file is an argument to the program rather than redirected 
   as unit 5.

.. describe:: Time Step

   The recommended and default time steps are now 60 s at 5-degree resolution and 30 s 
   at 2.5-degree resolution. We strongly urge users to use the recommended time steps. 
   Changing the time step will change the model results. Sometimes it may be absolutely 
   necessary to shorten the time step, e.g., during a major storm, but in that case it 
   is advisable to do so for only as long as absolutely necessary to get through the peak 
   of the storm.

.. describe:: Shapiro Filter

   The Shapiro filter factor was reduced to 3.0e-3 for the default time step, and now 
   changes if the time step is changed. This minimizes but does not eliminate varying 
   results for different time steps. The filter factor = 3.0e-3 x (time step)/(default time step), 
   i.e., it is reduced for shorter steps (see "shapiro" in src/cons.F).

.. describe:: FFT Filters

   The polar FFT filters were re-tuned for the 2.5-degree resolution model 
   (no change at 5-degree).

.. describe:: O+ Floor

   A double-Gaussian shaped floor (in latitude and altitude) is applied to O+ at 
   low-to-mid latitudes in the F-region in order to keep the model stable when the 
   ionosphere gets very low in density. It can be turned off in the namelist inputs 
   but this is not recommended. 
   See :base_url:`Namelist parameter ENFORCE_OPFLOOR <userguide/html/namelist.html#enforce-opfloor>`

.. describe:: O+ Diffusion Limiter

   An optional O+ diffusion coefficient limit can be supplied as a namelist input. 
   This can improve model stability in the topside F-region, but it is only 
   recommended as a last resort since it will change model results.
   See :base_url:`Namelist parameter OPDIFFCAP <userguide/html/namelist.html#opdiffcap>`

Bug Fixes
---------

.. describe:: Reduced Temperature

   The reduced temperature (Ti+Tn)/2 was erroneously coded in the O+ transport 
   routine as (Ti+Te)/2 (which is actually the plasma temperature). This bug has 
   been there for a long time, probably dating back to the TIGCM in the 1980’s. 
   This was corrected, which makes small but significant changes in the ionosphere, 
   particularly in the vicinity of the F2 peak. NmF2 generally changes by <10%, 
   but this is nevertheless a significant correction. The model is also somewhat 
   less stable at 2.5-degrees due to this correction. Note that the plasma temperature 
   is also used in O+ transport but this was correctly defined in the code.

.. describe:: Other Miscellaneous Corrections  

   * The domain decomposition was adjusted so that the model can run with 4 processors
   * ZG, ZGmid and DEN extrapolation to the top level was corrected.
   * GLAT and ALT input/output to apex.F90 was updated
   * A problem using ESMF with a single processor was fixed.

Make/Build System
-----------------

.. describe:: Platforms

   The TIE-GCM can be run on the NCAR yellowstone supercomputer or on Linux desktop 
   computers. Some users have run it on PC and Mac, and desktops and laptops, but 
   the setup may require some re-configuration by the user. Compiling with MPI is 
   now required on any platform, but the model may still be run using a single processor. 
   We recommend 4 to 16 processors for the 5-degree resolution model and 8 to 64 
   processors for the 2.5-degree resolution model. 

.. describe:: Compilers

   The TIE-GCM can be run on three Fortran compilers for Linux desktop systems: 
   Intel, PGI, and GNU gfortran. On the NCAR yellowstone system, only Intel is 
   used. Intel is the default compiler, because it consistently out-performs the 
   other compilers.

.. describe:: ESMF Library
 
   The parallel dynamo now requires that the ESMF library be linked. 
   The NetCDF library is also still required. See 
   :base_url:`Grid Structure and Resolution <userguide/html/build.html>`
   for more information.

.. note::
  A warning for user's of previous revisions of TIEGCM: Do not use old namelist
  input files or job scripts from previous revisions. Copy the default files
  from the :term:`scripts/` directory, and modify them for your own runs.
  Also, for initial runs, do not build/execute the model in an old :term:`execdir`. 
  Instead, allow the job script to make a new execdir for you.

Benchmarking Utility
--------------------

.. describe:: Tgcmrun

   Tgcmrun is a Python code (in the tgcmrun/ directory) that can be used to make 
   selected benchmark runs of the model in a semi-automated fashion. Tgcmrun can 
   be executed interactively on the command line, or a series of runs can be submitted 
   from a simple shell script using command-line arguments. See example scripts ‘run_*’ 
   or type ‘tgcmrun -h’ in the tgcmrun directory a usage statement. The tgcmrun 
   directory takes the place of the former ‘tests’ 

License
-------

.. describe:: Academic Research License

   The open-source academic research license was updated. This license requires that 
   any derivatives of the model also be open-source, that the model can be used for 
   academic research purposes only, and that the code cannot be sold. The updates 
   clarify that the prohibition on sales includes derivative products as well as code, 
   and that operational use is not permitted. Any existing licenses for v. 1.9* are 
   replaced by this license for v. 2.0. See the file 
   :download:`tiegcmlicense.txt <../../src/tiegcmlicense.txt>` here and in the src 
   directory.

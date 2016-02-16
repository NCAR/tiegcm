
Release Notes for |modeluc| |version|
==================================

:Release date: [not yet released as of January, 2016]
:Contact: discussion group email list `tgcmgroup@ucar.edu <http://mailman.ucar.edu/mailman/listinfo/tgcmgroup>`_
:Download page: http://www.hao.ucar.edu/modeling/tgcm/download.php

This is a summary of modifications made to the |model| since the release of 
|model|\1.95 (June 21, 2013)

Primary Modifications to Source Code
-------------------------------------

  * **2.5-degree resolution**. 
    The model has now been tuned to our satisfaction for community use at the 2.5-degree
    horizontal resolution (with four grid-points per scale height in the vertical coordinate).
    User's are encouraged to run the model at the 2.5-degree resolution, as it can produce
    a more realistic ionosphere than the 5-degree model. All benchmark runs and source history 
    files are available at both 5-degree and 2.5-degree resolutions.

  * **Parallel dynamo**. The electro-dynamo code was parallelized with pure MPI. Transformations
    between the geographic and geomagnetic grids in the MPI environment are accomplished with 
    the Earth System Modeling Framework (ESMF library). This resulted in a performance speed-up
    of about 25% at the 5-degree resolution, and almost 50% at the 2.5-degree resolution
    when using 64-processors on the NCAR |ncarsuper| system.

  * **Helium**. Helium is now calculated as a major species at both resolutions.
    If a source history without helium is used, helium will be initialized to a 
    global value (0.1154e-5), and will evolve from there. Helium is always saved 
    on primary histories, and can be saved on secondary histories as 'HE'. Thanks 
    to Eric Sutton for leading this effort.

  * **Argon**. Argon is now calculated as a minor species. Argon is always saved on 
    primary histories, and can be saved on secondary histories as 'AR'.

  * **Additional Diagnostic Fields**. Several new diagnostic fields have been added, 
    including N2, ZGMID, CUSP, DRIZZLE, ALFA, NFLUX, and EFLUX (see src/diags.F). 
 
  * **Namelist input**. The comment character in namelist input files is now '!'
    (formerly it was ';'). This was changed to conform to the fortran standard.
    If you have many input files with the ';' comment character, use the ``change_nlcomment``
    script in the scripts directory to replace the comment chars. Reading the input 
    file was simplified in the source code (see src/input_read.F), and in the job
    scripts, the input file is an argument to the program rather than redirected as unit 5.

Make/Build System
-----------------

  * **Compilers**. We now support three compilers for Linux desktop systems: 
    Intel, PGI, and GNU gfortran. On the NCAR |ncarsuper| system, only Intel
    is supported.  Intel is the default compiler on both systems, because it
    consistently out-performs the other compilers.

  * **ESMF Library**. The parallel dynamo now requires that the ESMF library
    be linked. See the User's Guide for more information. 

Benchmarking Utility
------------------------

  Tgcmrun is a Python code (in the tgcmrun/ directory) that can be used to 
  make selected benchmark runs of the model in a semi-automated fashion. 
  Tgcmrun can be executed interactively on the command line, or a series
  of runs can be submitted from a simple shell script using command-line
  arguments. See example scripts 'run_*' or type 'tgcmrun -h' in the 
  tgcmrun directory for more information. The tgcmrun directory takes
  the place of the former 'tests' directory.

Performance
-----------
  Please see :download:`Performance Table <_static/perf.table>` for estimates
  of performance timings on Linux desktop systems, and the NCAR |ncarsuper|
  system, using a variety of processor counts.


Release Notes for tiegcm |version|
==================================

:Release date: [not yet released as of January, 2016]
:Contact: discussion group email list `tgcmgroup@ucar.edu <http://mailman.ucar.edu/mailman/listinfo/tgcmgroup>`_
:Download page: http://www.hao.ucar.edu/modeling/tgcm/download.php

Summary of modifications made to the source code since the release of tiegcm1.95 (June 21, 2013)
------------------------------------------------------------------------------------------------

**Primary Modifications:**

  * **2.5-degree resolution**. 
    The model has now been tuned to our satisfaction for community use at the 2.5-degree
    horizontal resolution (with four grid-points per scale height in the vertical coordinate).
    User's are encouraged to run the model at the 2.5-degree resolution. All benchmark runs
    and source history files are available at both 5-degree and 2.5-degree resolutions.

  * **Parallel dynamo**. The electro-dynamo code was parallelized with pure MPI. Transformations
    between the geographic and geomagnetic grids in the MPI environment are accomplished with 
    the Earth System Modeling Framework (ESMF library). This resulted in a performance speed-up
    of about 25% at the 5-degree resolution, and almost 50% at the 2.5-degree resolution
    when using 64-processors on the NCAR |ncarsuper| system.

  * **Helium**. Helium is now calculated as a major species at both resolutions.
    Thanks to Eric Sutton for leading this effort. If a source history without helium
    is used, helium will be initialized to a global value (0.1154e-5), and will evolve
    from there.

  * **Argon**. Argon is now calculated as a minor species.

  * **Additional diagnostics**. Several new diagnostic fields have been added (see diags.F).

**Make/Build System:**

  * **Compilers**. We now support three compilers for Linux desktop systems: 
    Intel, PGI, and GNU gfortran. On the NCAR |ncarsuper| system, only Intel
    is supported.  Intel is the default compiler on both systems, because it
    consistently out-performs the other compilers.

  * **ESMF Library**. The parallel dynamo now requires that the ESMF library
    be linked. See the User's Guide for more information. 

**Performance:**::
    Please see :download:`Performance Table <_static/perf.table>` for estimates
    of performance timings on Linux desktop systems, and the NCAR |ncarsuper|
    system, using a variety of processor counts.

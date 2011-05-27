
.. _tests:

Benchmark Test Runs
===================

A series of benchmark runs are made for each release of the TIEGCM
(since v1.93).  Job scripts and namelist input files for these runs are 
available in subdirectories in the :term:`tests/` directory.

Netcdf files with the first history of each benchmark run are available in
the :ref:`data download file <download>`.  These files can be used as start-up 
:ref:`SOURCE <SOURCE>` files to reproduce the runs.  The namelist input files and 
job scripts used to make the runs are provided in the :term:`tests/` directory
of the release.

For full history file output and "sanity check" plots for benchmark runs made
by version |version|, please see `Release Documentation 
<http://download.hao.ucar.edu/pub/tgcm/tiegcm1.94/release/_build/html/>`_

Following is a summary of benchmark runs made by TIEGCM version |version|:

#. **control**: 5-day control runs, started from steady-state histories at both equinoxes
   and both solstices, and at solar minimum and maximum conditions::

    mareqx: March Equinox (day 80) 
    junsol: June Solstice (day 172) 
    sepeqx: September Equinox (day 264) 
    decsol: December Solstice (day 355)

   * Solar Minimum::

      POWER   = 18.
      CTPOTEN = 30.
      F107    = 70.
      F107A   = 70.

   * Solar Maximum::

      POWER   = 39.
      CTPOTEN = 60.
      F107    = 200.
      F107A   = 200.

#. **climatology**: Full-year Climatology with constant solar forcing:

   * Heelis potential model with constant solar forcing::

      POWER   = 18.
      CTPOTEN = 30.
      F107    = 100.
      F107A   = 100.

#. **dec2006**: December, 2006 "AGU" storm case:

   * Heelis potential model with GPI (Kp) data
   * Weimer potential model with IMF data (F10.7 from GPI)

#. **nov2003**: November 19-24 (days 323-328), 2003 storm case:

   * Heelis potential model with GPI (Kp) data
   * Weimer potential model with IMF data (F10.7 from GPI)

#. **whi2008**: Whole Heliosphere interval (WHI) (March 21 to April 16, 2008)

   * Heelis potential model with GPI (Kp) data
   * Weimer potential model with IMF data (F10.7 from GPI)


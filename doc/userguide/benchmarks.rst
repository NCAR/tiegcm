
.. _benchmarks:

Benchmark Runs
==============

A series of benchmark runs are made for each major release of the |model|
These runs (for |tgcm_version|) were made using the python scripts in the 
:term:`tgcmrun/` directory.

Netcdf files with the first history of each benchmark run are available in
the :ref:`data download file <download>`.  These files can be used as start-up 
:ref:`SOURCE <SOURCE>` files to reproduce the runs.  

Following is a summary of benchmark runs made by |model| version |version|.
For full history file output, steady-state start-up files, and "sanity check" plots 
of the benchmark runs made by version |version|, please see 
:base_url:`Release Documentation <release/html>`

Seasonal
--------

**Seasonal**: 5-day control runs, started from steady-state histories at both equinoxes
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

    POWER   = 40.
    CTPOTEN = 60.
    F107    = 200.
    F107A   = 200.

Climatology
-----------

**climatology**: Full-year Climatology with constant solar forcing:

 * Heelis potential model with constant solar forcing::

    POWER   = 18.
    CTPOTEN = 30.
    F107    = 100.
    F107A   = 100.

December, 2006 "AGU" Storm Case
-------------------------------

**dec2006**: December, 2006 "AGU" storm case:

 * Heelis potential model with GPI (Kp) data
 * Weimer potential model with IMF data (F10.7 from GPI)

November, 2003 Storm Case
-------------------------

**nov2003**: November 19-24 (days 323-328), 2003 storm case:

 * Heelis potential model with GPI (Kp) data
 * Weimer potential model with IMF data (F10.7 from GPI)

Whole Heliosphere Interval
--------------------------

**whi2008**: Whole Heliosphere interval (WHI) (March 21 to April 16, 2008)

 * Heelis potential model with GPI (Kp) data
 * Weimer potential model with IMF data (F10.7 from GPI)

.. note::
   For more detailed information and access to history file output, and
   preliminary post-processing of these runs, 
   see :base_url:`Release Documentation <release/html>`

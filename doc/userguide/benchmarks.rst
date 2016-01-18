
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

Making Benchmark Runs
---------------------

The :term:`tgcmrun <tgcmrun/>` directory under the model root directory 
:term:`(modeldir) <modeldir>` contains Python code that semi-automates
submission of selected benchmark model runs on the NCAR supercomputer
system (|ncarsuper|). The tgcmrun command can be executed interactively
on the command line, or from a shell script. Type "tgcmrun -h" on the 
command line for a detailed usage message.  Typing "tgcmrun" on the
command line will cause the program to print the available benchmark 
runs and prompt the user as follows::

  The following runs are available:

  NUMBER	NAME		DESCRIPTION
  ------	----		-----------
  0 	default_run 	Default run
  1 	decsol_smax 	December Solstice Solar Maximum
  2 	decsol_smin 	December Solstice Solar Minimum
  3 	junsol_smax 	June Solstice Solar Maximum
  4 	junsol_smin 	June Solstice Solar Minimum
  5 	mareqx_smax 	March Equinox Solar Maximum
  6 	mareqx_smin 	March Equinox Solar Minimum
  7 	sepeqx_smax 	September Equinox Solar Maximum
  8 	sepeqx_smin 	September Equinox Solar Minimum
  9 	nov2003_heelis_gpi 	November 2003 "Halloween Storm", Heelis potential model, GPI data
  10 	nov2003_weimer_imf 	November 2003 "Halloween Storm", Weimer potential model, IMF, GPI data
  11 	dec2006_heelis_gpi 	December 2006 "AGU storm", Heelis potential model, GPI data
  12 	dec2006_weimer_imf 	December 2006 "AGU storm", Weimer potential model, IMF and GPI data
  13 	whi2008_heelis_gpi 	2008 "Whole Heliosphere Interval", Heelis potential model, GPI data
  14 	whi2008_weimer_imf 	2008 "Whole Heliosphere Interval", Weimer potential model, IMF, GPI data
  15 	climatology_smin 	Climatology run with constant solar minimum conditions (Jan 1-5)
  16 	climatology_smax 	Climatology run with constant solar maximum conditions (Jan 1-5)

  Enter number of desired run (0-16) ('q' to quit, 'p' to print list, default=0): 

At this point the user can enter an integer 0 to 16, specifying the desired run.
The user will then be prompted for a few additional parameters (tiegcm or timegcm model,
resolution, model root directory, etc).  However, it is easiest to set a few environment
variables before executing tgcmrun, to minimize the need to enter long file paths at the 
prompt:

Environment variables to set before using the tgcmrun utility:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **TGCMTEMP**: Path to a large temporary directory where the model can be built, 
  executed, and output stored.
* **TGCMDATA**: Path to a directory containing data files required by the model 
  (netcdf data and start-up history files)
* **TIEGCM_ROOT**: Path to the tiegcm model root directory containing source code, 
  scripts, tgcmrun, etc. (not necessary if making only TIMEGCM runs)
* **TIMEGCM_ROOT**: Path to the timegcm model root directory containing source code, 
  scripts, tgcmrun, etc. (not necessary if making only TIEGCM runs)
 
Source history files (start-up netcdf files with a single history) to start these
runs are provided in the :ref:`data download <data_download>` (there are separate
data downloads available for each model resolution). These source files should be 
located in the :term:`TGCMDATA` directory.

The tgcmrun program can also be executed from a shell script. There are several
example tcsh scripts in the tgcmrun directory that make a series of runs for
various purposes:

* **run_climatology**: Make two short climatology runs, one for each resolution.
* **run_compilers**: Make three runs, each with a different compiler (linux desktop systems only)
* **run_perf**: Make several runs using different processor (MPI task) counts (super systems only)
* **run_scriptsonly**: This only makes the namelist input and job scripts (does not submit the jobs)
* **run_seasons**: Make the 8 seasonal benchmark runs (at both resolutions).
* **run_storms**: Make storm case benchmark runs (at both resolutions).

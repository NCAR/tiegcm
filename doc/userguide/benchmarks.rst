
.. _benchmarks:

Benchmark Runs
==============

.. note::

  Benchmark results (pdf plot files) for version |tgcm_version| are available here: 
  :base_url:`Release Benchmarks Results <release/html/benchmarks.html>`

A series of benchmark runs are made for each major release of the |model|
These runs (for |tgcm_version|) were made using the python scripts in the 
:term:`tgcmrun/` directory.

Netcdf files with the first history of each benchmark run are available in
the :ref:`data download file <download>`.  These files can be used as start-up 
:ref:`SOURCE <SOURCE>` files to reproduce the runs.  

Benchmark runs are provided in three groups: full-year climatology, seasonal, and solar storm events.
Here is a table of runs at 2.5-degree model resolution (there is also a set at the 5-degree resolution)::
 
 Full-year climatologies (solar min,max):
   tiegcm_res2.5_climatology_smax	Full-year climatology at solar maximum conditions
   tiegcm_res2.5_climatology_smin	Full-year climatology at solar minimum conditions
 
 Seasons (Equinox, Solstice, solar min,max):
   tiegcm_res2.5_decsol_smax		December Solstice, solar maximum (days 355-360)
   tiegcm_res2.5_decsol_smin		December Solstice, solar minimum (days 355-360)
   tiegcm_res2.5_junsol_smax		June Solstice, solar maximum     (days 172-177)
   tiegcm_res2.5_junsol_smin		June Solstice, solar minimum     (days 172-177)
   tiegcm_res2.5_mareqx_smax		March Equinox, solar maximum     (days 80-85)
   tiegcm_res2.5_mareqx_smin		March Equinox, solar minimum     (days 80-85)
   tiegcm_res2.5_sepeqx_smax		September Equinox, solar maximum (days 264-269)
   tiegcm_res2.5_sepeqx_smin		September Equinox, solar minimum (days 264-269)
  
 Solar storm events (Heelis/GPI and Weimer/IMF,GPI):
   tiegcm_res2.5_dec2006_heelis_gpi	December 2006 storm, Heelis/GPI          (days 330-360)
   tiegcm_res2.5_dec2006_weimer_imf	December 2006 storm, Weimer/IMF,GPI      (days 330-360)
   tiegcm_res2.5_jul2000_heelis_gpi	July 2000 storm,     Heelis/GPI          (days 192-202)
   tiegcm_res2.5_jul2000_weimer_imf	July 2000 storm,     Weimer/IMF,GPI      (days 192-202)
   tiegcm_res2.5_nov2003_heelis_gpi	November 2003 storm, Heelis/GPI          (days 323-328)
   tiegcm_res2.5_nov2003_weimer_imf	November 2003 storm, Weimer/IMF,GPI      (days 323-328)
   tiegcm_res2.5_whi2008_heelis_gpi	Whole Helio Interval 2008 Heelis/GPI     (days 81-106)
   tiegcm_res2.5_whi2008_weimer_imf	Whole Helio Interval 2008 Weimer/IMF,GPI (days 81-106)

.. note::

  Seasonal runs and full-year Climatologies were run with constant solar forcing, as follows:

 * Solar Minimum:

  | POWER   = 18.
  | CTPOTEN = 30.
  | F107    = 70.
  | F107A   = 70.

 * Solar Maximum:

  | POWER   = 40.
  | CTPOTEN = 60.
  | F107    = 200.
  | F107A   = 200.

For comprehensive plots of all benchmark runs, please see the |tgcm_version_uc|
:base_url:`Release Documentation <release_html>`

Seasonal (equinoxes, solstices, solar min, max):

  * decsol: December Solstice (days 355-360)
  * junsol: June Solstice (days 172-177) 
  * mareqx: March Equinox (days 80-85) 
  * sepeqx: September Equinox (days 264-269) 

Climatologies (full-year runs with daily histories):

  * climatology_smax: Climatologies at Solar Maximum
  * climatology_smin: Climatologies at Solar Minimum

December, 2006 "AGU" Storm Case (days 330-360):

 * dec2006: Heelis potential model with GPI (Kp) data
 * dec2006: Weimer potential model with IMF data (F10.7 from GPI)

November 19-24, 2003 Storm Case (days 323-328)

 * nov2003: Heelis potential model with GPI (Kp) data
 * nov2003: Weimer potential model with IMF data (F10.7 from GPI)

July 11-21, 2000 "Bastille Day" Storm Case (days 192-202)

 * jul2000: Heelis potential model with GPI (Kp) data
 * jul2000: Weimer potential model with IMF data (F10.7 from GPI)

Whole Heliosphere Interval (WHI) (March 21 to April 16, 2008)

 * whi2008: Heelis potential model with GPI (Kp) data
 * whi2008: Weimer potential model with IMF data (F10.7 from GPI)

.. note::
   For more detailed information and access to history file output, and
   extensive post-processing of these runs, see the |tgcm_version_uc|
   :base_url:`Release Documentation <release_html>`

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
 9 	nov2003_heelis_gpi 	November 2003 storm case, Heelis potential model, GPI data
 10 	nov2003_weimer_imf 	November 2003 storm case, Weimer potential model, IMF, GPI data
 11 	dec2006_heelis_gpi 	December 2006 "AGU storm", Heelis potential model, GPI data
 12 	dec2006_weimer_imf 	December 2006 "AGU storm", Weimer potential model, IMF and GPI data
 13 	whi2008_heelis_gpi 	2008 "Whole Heliosphere Interval", Heelis potential model, GPI data
 14 	whi2008_weimer_imf 	2008 "Whole Heliosphere Interval", Weimer potential model, IMF, GPI data
 15 	jul2000_heelis_gpi 	July 2000 "Bastille Day" storm, Heelis potential model, GPI data
 16 	jul2000_weimer_imf 	July 2000 "Bastille Day" storm, Weimer potential model, IMF, GPI data
 17 	climatology_smin 	Climatology run with constant solar minimum conditions (Jan 1-5)
 18 	climatology_smax 	Climatology run with constant solar maximum conditions (Jan 1-5)
 
 Enter number of desired run (0-18) ('q' to quit, 'p' to print list, default=0): 

At this point the user can enter an integer 0 to 18, specifying the desired run.
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
located in the :term:`TGCMDATA` directory (or the path may be specified in the 
job script with the :term:`tgcmdata` shell variable).

The tgcmrun program can also be executed from a shell script. There are several
example tcsh scripts in the tgcmrun directory that make series of runs for
various purposes. The scripts optionally run at one or both model resolutions.
History files, stdout log files, and job scripts used, are stored in a directory
tree below the working directory. 

Standard 18 benchmark runs (as in the interactive tgcmrun command above):

* **run_climatology**: Start climatology runs (smin,smax). These can be extended to a full year by the user.
* **run_seasons**: Make seasonal benchmark runs (equinoxes, solstices, at smin, smax)
* **run_storms**: Make storm case benchmark runs (heelis_gpi and weimer_imf)

Additional runs for testing compilers, performance, etc.:

* **run_compilers**: Make three runs, each with a different compiler (linux desktop systems only)
* **run_perf**: Make several runs using different processor (MPI task) counts (super systems only)
* **run_scriptsonly**: This only makes the namelist input and job scripts (does not submit the jobs)

.. _benchmark_history_files:

Model Output History Files of the |tgcm_version| Benchmark runs
---------------------------------------------------------------

Model output history files are stored in CF-compliant netCDF format (see :ref:`historyoutput`).
Benchmark history files are available via `Globus research data sharing service <https://www.globus.org>`_.
The tiegcm benchmark history files are stored at the "NCAR Data Sharing Service" :term:`Globus` 
shared endpoint (for users with an NCAR/CISL login: this endpoint is /glade/u/datashare/tgcm).

See these CISL docs for information regarding the NCAR Data Sharing Service:

  * `NCAR Data Sharing Service <https://www2.cisl.ucar.edu/resources/storage-and-file-systems/using-the-ncar-data-sharing-service>`_

  * `Globus file transfers     <https://www2.cisl.ucar.edu/resources/storage-and-file-systems/globus-file-transfers>`_
    (see especially "Transferring files with the webh interface")

  * `Retreiving data from a shared endpoint <https://www2.cisl.ucar.edu/resources/storage-and-file-systems/using-the-ncar-data-sharing-service#retrieve>`_

Here is a summary procedure for accessing the |tgcm_version| benchmark data:

.. note::

  You do *NOT* have to have an NCAR user account or token to retrieve this data.


* You must have or create a `Globus <https://www.globus.org>`_ account. If your 
  institution/organization has a Globus data sharing endpoint, you can use your institutional 
  authorization to login to Globus. Otherwise, you can create a 
  `Globus personal account <https://www.globus.org/SignUp>`_ to transfer files
  to your personal laptop or desktop computer.

* Log  in to your Globus account, and click on "File Transfer"

* To reach the NCAR/TIEGCM source endpoint, click in the "Endpoint" text box on the left, 
  and type "TIEGCM v2.0". It should retrieve directory contents, and show a "benchmarks" folder.

* Next, establish your destination endpoint on the right. This is either your institutional
  endpoint, or the username of your personal Globus login.

* Select the locations/files you want to download from the left side, and the destination
  location on the right, then click the right arrow at the top to begin the transfer.

Here's a screen shot of a Globus file transfer from the TIEGCM v2.0 endpoint to my personal
Macbook Pro: :download:`Globus_screenshot.png <_static/globus_screenshot.png>`

In each of the 6 benchmark groups are folders for each run, with folders containing
the history files (hist), post-processing (proc), and scripts and log files (stdout). 
Individual files or whole directories can be downloaded.

.. note::
  Users wanting to use their NCAR authentication rather than personal GlobusID, apparently 
  need to have a login on the `NCAR RDA <http://rda.ucar.edu/>`_ (Research Data Archives) 
  to access the NCAR GLADE endpoint on Globus.


.. _benchmarks:

|modeluc| Version |version| Benchmark Runs
==========================================

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

 Climatology and Seasonal runs were made with constant solar conditions as follows:

 * Solar Minimum Conditions:
   
   | POWER   = 18.
   | CTPOTEN = 30.
   | F107    = 70.
   | F107A   = 70.
   
 * Solar Maximum Conditions:            
   
   | POWER   = 40.
   | CTPOTEN = 60.
   | F107    = 200.
   | F107A   = 200.
 
.. note::
  
  Explanation of terms for Storm Case Benchmarks:

    * :term:`Heelis`: The Heelis Potential Model was used for high-latitude ion drift velocities

    * :term:`Weimer`: The Weimer Potential Model was used for high-latitude ion drift velocities

    * :term:`GPI`: Runs used GeoPhysical Indices data with Kp, and solar flux f10.7d, f10.7a 

    * :term:`IMF`: Runs used IMF/OMNI data with Solar Wind density and velocity

.. _benchmark_results:

Selected Results of |tgcm_version_uc| Benchmark Runs 
----------------------------------------------------

Results of the benchmark runs are available in 2d plots of time series, grid slices, and maps.
These are multi-page PDF files, arranged as follows:

  * Files with "singleut" in their names contain instantaneous maps and slices
  * Files with "utvert" in their names contain ut vs pressure plots (daily or hourly histories)
  * Files with "utlat" in their names contain ut vs latitude plots (daily or hourly histories)

Instantaneous and time series plots for **all** fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**All** fields are: 
  | 'TN','UN','VN','WN','O2','O1','N2','NO','N4S','HE','NE','TE','TI','TEC','O2P','OP',
  | 'POTEN','UI_ExB','VI_ExB','WI_ExB','DEN','QJOULE','HMF2','NMF2','Z'

At 5-degree model resolution:

  * :benchmarks_url:`Climatology Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_climatology/PDFs_all_fields>`
  * :benchmarks_url:`Seasonal Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_seasons/PDFs_all_fields>`
  * :benchmarks_url:`Storm Case Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_storms/PDFs_all_fields>`

At 2.5-degree model resolution:

  * :benchmarks_url:`Climatology Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_climatology/PDFs_all_fields>`
  * :benchmarks_url:`Seasonal Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_seasons/PDFs_all_fields>`
  * :benchmarks_url:`Storm Case Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_storms/PDFs_all_fields>`

Instantaneous and time series plots for **selected** fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Selected** fields are:
  | 'TN','UN','VN','WN','HE','NE','TE','TI','HMF2','NMF2','POTEN','Z'

At 5-degree model resolution:

  * :benchmarks_url:`Climatology Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_climatology/PDFs_select_fields>`
  * :benchmarks_url:`Seasonal Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_seasons/PDFs_select_fields>`
  * :benchmarks_url:`Storm Case Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_storms/PDFs_select_fields>`

At 2.5-degree model resolution:

  * :benchmarks_url:`Climatology Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_climatology/PDFs_select_fields>`
  * :benchmarks_url:`Seasonal Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_seasons/PDFs_select_fields>`
  * :benchmarks_url:`Storm Case Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_storms/PDFs_select_fields>`

Namelist input files used for the 2.5-deg benchmark runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These are the namelist input files used for making the 2.5-degree resolution benchmark runs.
Note that a few of the 2.5-deg runs need to be run at a shorter timestep than the default of
30 seconds to maintain numerical stability. Also note that the 5-degree full-year climatology
run can be completed on |ncarsuper| in 12 hours WC, but the 2.5-deg model has to be restarted 
several times to reach a full year.

* Climatology Benchmark Runs (full-year):

  * :download:`Solar max climatology <_static/namelist_files/tiegcm2.0_res2.5_climatology_smax.inp>`
  * :download:`Solar min climatology <_static/namelist_files/tiegcm2.0_res2.5_climatology_smin.inp>`

* Seasonal Benchmark Runs (5-day runs):

  * :download:`December solstice solar max <_static/namelist_files/tiegcm2.0_res2.5_decsol_smax.inp>`
  * :download:`December solstice solar min <_static/namelist_files/tiegcm2.0_res2.5_decsol_smin.inp>`
  * :download:`June solstice solar max <_static/namelist_files/tiegcm2.0_res2.5_junsol_smax.inp>`
  * :download:`June solstice solar min <_static/namelist_files/tiegcm2.0_res2.5_junsol_smin.inp>`
  * :download:`March equinox solar max <_static/namelist_files/tiegcm2.0_res2.5_mareqx_smax.inp>`
  * :download:`March equinox solar min <_static/namelist_files/tiegcm2.0_res2.5_mareqx_smin.inp>`
  * :download:`September equinox solar max <_static/namelist_files/tiegcm2.0_res2.5_sepeqx_smax.inp>`
  * :download:`September equinox solar min <_static/namelist_files/tiegcm2.0_res2.5_sepeqx_smin.inp>`

* Storm Simulations: 

  * :download:`December, 2006 storm (Heelis/gpi)     <_static/namelist_files/tiegcm2.0_res2.5_dec2006_heelis_gpi.inp>`
  * :download:`December, 2006 storm (Weimer/gpi,imf) <_static/namelist_files/tiegcm2.0_res2.5_dec2006_weimer_imf.inp>`
  * :download:`July, 2000 storm (Heelis/gpi)     <_static/namelist_files/tiegcm2.0_res2.5_jul2000_heelis_gpi.inp>`
  * :download:`July, 2000 storm (Weimer/gpi,imf) <_static/namelist_files/tiegcm2.0_res2.5_jul2000_weimer_imf.inp>`
  * :download:`November, 2003 storm (Heelis/gpi)     <_static/namelist_files/tiegcm2.0_res2.5_nov2003_heelis_gpi.inp>`
  * :download:`November, 2003 storm (Weimer/gpi,imf) <_static/namelist_files/tiegcm2.0_res2.5_nov2003_weimer_imf.inp>`
  * :download:`2008 Whole Heliosphere Interval (Heelis/gpi)    <_static/namelist_files/tiegcm2.0_res2.5_whi2008_heelis_gpi.inp>`
  * :download:`2008 Whole Heliosphere Interval (Weimer/gpi,imf <_static/namelist_files/tiegcm2.0_res2.5_whi2008_weimer_imf.inp>`

.. _benchmark_history_files:

Model Output History Files of the |tgcm_version| Benchmark runs
---------------------------------------------------------------

Model output history files are stored in CF-compliant netCDF format.
Benchmark history files for |tgcm_version| are available via the 
`Globus Research Data Sharing Service <https://www.globus.org>`_.
The history files are stored at the "NCAR Data Sharing Service" :term:`Globus` shared endpoint.
(for users with an NCAR/CISL login: this endpoint is /glade/u/datashare/tgcm)

See these CISL docs for information regarding the NCAR Data Sharing Service:

  * `NCAR Data Sharing Service <https://www2.cisl.ucar.edu/resources/storage-and-file-systems/using-the-ncar-data-sharing-service>`_

  * `Globus file transfers     <https://www2.cisl.ucar.edu/resources/storage-and-file-systems/globus-file-transfers>`_
    (see especially "Transferring files with the webh interface")

  * `Retreiving data from a shared endpoint <https://www2.cisl.ucar.edu/resources/storage-and-file-systems/using-the-ncar-data-sharing-service#retrieve>`_

Here is a summary procedure for accessing the |tgcm_version| benchmark data:

.. note::

  You do *NOT* need to have an NCAR user account or token to retrieve this data.

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

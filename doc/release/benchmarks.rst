
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

Selected Results of Benchmark Runs
----------------------------------

Results of the benchmark runs are available in 2d plots of time series, grid slices, and maps.
These are multi-page PDF files (you may need to use the rotate tool in your browser or previewer):

Runs at 5.0-degree resolution:

  * :benchmarks_url:`Climatology Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_climatology/PDFs>`
  * :benchmarks_url:`Seasonal Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_seasons/PDFs>`
  * :benchmarks_url:`Storm Case Benchmark Runs at 5.0-deg resolution <tiegcm_res5.0_storms/PDFs>`

Runs at 2.5-degree resolution:

  * :benchmarks_url:`Climatology Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_climatology/PDFs>`
  * :benchmarks_url:`Seasonal Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_seasons/PDFs>`
  * :benchmarks_url:`Storm Case Benchmark Runs at 2.5-deg resolution <tiegcm_res2.5_storms/PDFs>`


Availability of Model Output History Files 
------------------------------------------

Model output history files are stored in CF-compliant netCDF format.
Benchmark history files are available from the Globus data-sharing service.

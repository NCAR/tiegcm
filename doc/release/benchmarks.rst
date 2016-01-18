
.. _benchmarks:

Benchmark Runs
==============

For seasonal start-up files, and control runs, the definition
of "Solar Minimum" and "Solar Maximum" conditions are as follows:

* Solar Minimum Conditions::

   POWER   = 18.
   CTPOTEN = 30.
   F107    = 70.
   F107A   = 70.

* Solar Maximum Conditions::

   POWER   = 40.
   CTPOTEN = 60.
   F107    = 200.
   F107A   = 200.

Steady-State Seasonal Start-up Files
------------------------------------

Single-history start-up files for the four seasons (at both model resolutions) are available
at the `tgcm download page <http://www.hao.ucar.edu/modeling/tgcm/download.php>`_ 
These files can be used as SOURCE files in the namelist input for initial runs.
Files listed below are for the 5.0-degree resolution. Files for benchmark runs
at the 2.5-degree resolution are the same, but have 'res2.5' instead of 'res5.0'.

Start-up files at Solar Minimum Conditions:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

=============================== ====================================
Date				Filenames (5.0 and 2.5-deg res)
=============================== ====================================
March Equinox (day 80)          tiegcm_res5.0_decsol_smin_prim.nc  
                                tiegcm_res2.5_decsol_smin_prim.nc
June Solstice (day 172)         tiegcm_res5.0_junsol_smin_prim.nc  
                                tiegcm_res2.5_junsol_smin_prim.nc
September Equinox (day 264)     tiegcm_res5.0_sepeqx_smin_prim.nc  
                                tiegcm_res2.5_sepeqx_smin_prim.nc
December Solstice (day 355)     tiegcm_res5.0_decsol_smin_prim.nc  
                                tiegcm_res2.5_decsol_smin_prim.nc
=============================== ====================================

Start-up files at Solar Maximum Conditions:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

=============================== ====================================
Date				Filename (5.0 and 2.5-deg res)
=============================== ====================================
March Equinox (day 80)          tiegcm_res5.0_decsol_smax_prim.nc
                                tiegcm_res2.5_decsol_smax_prim.nc
June Solstice (day 172)         tiegcm_res5.0_junsol_smax_prim.nc
                                tiegcm_res2.5_junsol_smax_prim.nc
September Equinox (day 264)     tiegcm_res5.0_sepeqx_smax_prim.nc
                                tiegcm_res2.5_sepeqx_smax_prim.nc
December Solstice (day 355)     tiegcm_res5.0_decsol_smax_prim.nc
                                tiegcm_res2.5_decsol_smax_prim.nc
=============================== ====================================


.. _seasonal_control:

Seasonal Control Runs
---------------------

**control**: 5-day control runs, started from the above steady-state histories at 
both equinoxes and solstices, and at solar minimum and maximum conditions.

* "Sanity check" plots for solar min control runs (global maps on last day of 5-day run): 

  * :ref:`December Solstice, Solar Minimum <control_decsol_smin>`
  * :ref:`June Solstice, Solar Minimum <control_junsol_smin>`
  * :ref:`March Equinox, Solar Minimum <control_mareqx_smin>`
  * :ref:`September Equinox, Solar Minimum <control_sepeqx_smin>`

* "Sanity check" plots for solar max control runs (global maps on last day of 5-day run): 

  * :ref:`December Solstice, Solar Maximum <control_decsol_smax>`
  * :ref:`June Solstice, Solar Maximum <control_junsol_smax>`
  * :ref:`March Equinox, Solar Maximum <control_mareqx_smax>`
  * :ref:`September Equinox, Solar Maximum <control_sepeqx_smax>`

* These history files are bundled with the tiegcm1.95 data download:

  * `<http://www.hao.ucar.edu/modeling/tgcm/download.php>`_

.. _climatology_info:

Climatologies
-------------

**climatology smax**: Full-year Climatology with constant solar maximum forcing,
  with the Heelis potential model::

   POWER   = 40.
   CTPOTEN = 60.
   F107    = 200.
   F107A   = 200.

* :ref:`"Sanity check" plots for Climatology at Solar Max <climatology_smax>`

**climatology smin**: Full-year Climatology with constant solar minimum forcing,
  with the Heelis potential model::

   POWER   = 18.
   CTPOTEN = 30.
   F107    = 100.
   F107A   = 100.

* :ref:`"Sanity check" plots for Climatology at Solar Min <climatology_smin>`

.. _dec2006:

December, 2006 "AGU Storm"
--------------------------

**dec2006**: December, 2006 "AGU" storm case:

* Heelis potential model with GPI (Kp) data
* Weimer potential model with IMF data (F10.7 from GPI)
* "Sanity check" plots for dec2006:

  * :ref:`Heelis/GPI (daily, days 330-360) <dec2006_heelis_daily>`
  * :ref:`Heelis/GPI (hourly, days 346-352) <dec2006_heelis_hourly>`
  * :ref:`Weimer/IMF (daily, days 330-360) <dec2006_weimer_daily>`
  * :ref:`Weimer/IMF (hourly, days 346-352) <dec2006_weimer_hourly>`

.. _nov2003:

November, 2003 Storm Case
-------------------------

**nov2003**: November 19-24 (days 323-328), 2003 storm case:

* Heelis potential model with GPI (Kp) data
* Weimer potential model with IMF data (F10.7 from GPI)
* "Sanity check" plots for nov2003:

  * :ref:`Heelis/GPI (hourly) <nov2003_heelis_hourly>`
  * :ref:`Weimer/IMF (hourly) <nov2003_weimer_hourly>`

.. _whi2008:

Whole Heliosphere Interval (WHI)
--------------------------------

**whi2008**: Whole Heliosphere interval (WHI) (March 21 to April 16, 2008)

* Heelis potential model with GPI (Kp) data
* Weimer potential model with IMF data (F10.7 from GPI)
* "Sanity check" plots for whi2008:

  * :ref:`Heelis/GPI (daily) <whi2008_heelis_daily>`
  * :ref:`Weimer/IMF (daily) <whi2008_weimer_daily>`

* Download history files at `TIEGCM Version 1.95 Whole Heliosphere Interval 2008 <http://cdp.ucar.edu/browse/browse.htm?uri=http%3a%2f%2fdataportal.ucar.edu%2fmetadata%2ftgcm%2fTIEGCM_Version_1_95%2fWhole_Heliosphere_Interval_2008%2fWhole_Heliosphere_Interval_2008.thredds.xml>`_

History files on the NCAR HPSS 
------------------------------

Seasonal start-up files and complete history files for the 
benchmark runs are on the `NCAR HPSS <http://www2.cisl.ucar.edu/docs/hpss>`_
in directory /home/tgcm/tiegcm\ |version|. Here is a 
:download:`complete catalog listing <_static/benchmarks.hpss>`
including "contents" annotations.


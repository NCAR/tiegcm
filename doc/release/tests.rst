
.. _tests:

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

   POWER   = 39.
   CTPOTEN = 60.
   F107    = 200.
   F107A   = 200.

Steady-State Seasonal Start-up Files
------------------------------------

Single-history start-up files for the four seasons are in the **data download file**, 
which is available from the `tgcm download page <http://www.hao.ucar.edu/modeling/tgcm/download.php>`_ 
These files can be used as SOURCE files in the namelist input for initial runs:

Start-up files at Solar Minimum Conditions:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

=============================== ====================================
Date				Filename
=============================== ====================================
March Equinox (day 80)          TGCM.tiegcm1.95.pcntr_mareqx_smin.nc
June Solstice (day 172)         TGCM.tiegcm1.95.pcntr_junsol_smin.nc
September Equinox (day 264)     TGCM.tiegcm1.95.pcntr_sepeqx_smin.nc
December Solstice (day 355)     TGCM.tiegcm1.95.pcntr_decsol_smin.nc
=============================== ====================================

Start-up files at Solar Maximum Conditions:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

=============================== ====================================
Date				Filename
=============================== ====================================
March Equinox (day 80)          TGCM.tiegcm1.95.pcntr_mareqx_smax.nc
June Solstice (day 172)         TGCM.tiegcm1.95.pcntr_junsol_smax.nc
September Equinox (day 264)     TGCM.tiegcm1.95.pcntr_sepeqx_smax.nc
December Solstice (day 355)     TGCM.tiegcm1.95.pcntr_decsol_smax.nc
=============================== ====================================

* Download these files at `TIEGCM Version 1.95 Steady State Source Histories <http://cdp.ucar.edu/browse/browse.htm?uri=http%3a%2f%2fdataportal.ucar.edu%2fmetadata%2ftgcm%2fTIEGCM_Version_1_95%2fSteady_State_Source_Histories%2fSteady_State_Source_Histories.thredds.xml>`_

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

Climatology
-----------

**climatology**: Full-year Climatology with constant solar forcing:

* Heelis potential model with constant solar forcing::

   POWER   = 18.
   CTPOTEN = 30.
   F107    = 100.
   F107A   = 100.

* :ref:`"Sanity check" plots for Climatology <climatology>`

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

.. 6/21/13 btf: Comment out v1.95 movies for now:

.. * Heelis/GPI Movies (~3.5M gif):
..
  * :ref:`Heelis/GPI hourly movies: TN at Zp -4 <dec2006_heelis_TN_zp-4_movie>`
  * :ref:`Heelis/GPI hourly movies: TN at Zp +2 <dec2006_heelis_TN_zp+2_movie>`
  * :ref:`Heelis/GPI hourly movies: NE at Zp +2 <dec2006_heelis_NE_zp+2_movie>`

.. * Weimer/IMF Movies (~3.5 gif):
..
  * :ref:`Weimer/IMF hourly movies: TN at Zp -4 <dec2006_weimer_TN_zp-4_movie>`
  * :ref:`Weimer/IMF hourly movies: TN at Zp +2 <dec2006_weimer_TN_zp+2_movie>`
  * :ref:`Weimer/IMF hourly movies: NE at Zp +2 <dec2006_weimer_NE_zp+2_movie>`

.. * :ref:`AVI Movies <dec2006_avi_movies>`

.. * Download history files at `TIEGCM Version 1.95 December, 2006 AGU Storm <http://cdp.ucar.edu/browse/browse.htm?uri=http%3a%2f%2fdataportal.ucar.edu%2fmetadata%2ftgcm%2fTIEGCM_Version_1_95%2fDecember_2006_AGU_Storm_Simulations%2fDecember_2006_AGU_Storm_Simulations.thredds.xml>`_

.. 6/21/13 btf: End commented movie section

.. _nov2003:

November, 2003 Storm Case
-------------------------

**nov2003**: November 19-24 (days 323-328), 2003 storm case:

* Heelis potential model with GPI (Kp) data
* Weimer potential model with IMF data (F10.7 from GPI)
* "Sanity check" plots for nov2003:

  * :ref:`Heelis/GPI (hourly) <nov2003_heelis_hourly>`
  * :ref:`Weimer/IMF (hourly) <nov2003_weimer_hourly>`

.. 6/21/13 btf: Comment out v1.95 movies for now:

.. * Heelis/GPI Movies (~3-5M gif):
..
  * :ref:`Heelis/GPI hourly movies: TN at Zp -4 <nov2003_heelis_TN_zp-4_movie>`
  * :ref:`Heelis/GPI hourly movies: TN at Zp +2 <nov2003_heelis_TN_zp+2_movie>`
  * :ref:`Heelis/GPI hourly movies: NE at Zp +2 <nov2003_heelis_NE_zp+2_movie>`

.. * Weimer/IMF Movies (~3-5M gif):
..
  * :ref:`Weimer/IMF hourly movies: TN at Zp -4 <nov2003_weimer_TN_zp-4_movie>`
  * :ref:`Weimer/IMF hourly movies: TN at Zp +2 <nov2003_weimer_TN_zp+2_movie>`
  * :ref:`Weimer/IMF hourly movies: NE at Zp +2 <nov2003_weimer_NE_zp+2_movie>`

.. * :ref:`AVI Movies <nov2003_avi_movies>`

.. * Download history files at `TIEGCM Version 1.95 November 2003 Storm Simulations <http://cdp.ucar.edu/browse/browse.htm?uri=http%3a%2f%2fdataportal.ucar.edu%2fmetadata%2ftgcm%2fTIEGCM_Version_1_95%2fNovember_2003_Storm_Simulations%2fNovember_2003_Storm_Simulations.thredds.xml>`_

.. 6/21/13 btf: End commented movie section

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


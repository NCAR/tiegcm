
.. _tests:

TIEGCM v\ |version| Benchmark Runs
==================================

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

Steady-State Seasonal Start-up Files:
-------------------------------------

Single-history start-up files are available at the four seasons, and may be
used as SOURCE files in the namelist input for initial runs (these files are
about 14.5 MB each):

* Start-up files at Solar Minimum:

 * `March Equinox (day 80) Solar Minimum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_mareqx_smin.nc>`_
 * `June Solstice (day 172) Solar Minimum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_junsol_smin.nc>`_
 * `September Equinox (day 264) Solar Minimum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_sepeqx_smin.nc>`_
 * `December Solstice (day 355) Solar Minimum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_decsol_smin.nc>`_

* Start-up files at Solar Maximum:

 * `March Equinox (day 80) Solar Maximum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_mareqx_smax.nc>`_
 * `June Solstice (day 172) Solar Maximum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_junsol_smax.nc>`_
 * `September Equinox (day 264) Solar Maximum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_sepeqx_smax.nc>`_
 * `December Solstice (day 355) Solar Maximum 
   <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/pcntr_decsol_smax.nc>`_

Seasonal Control Runs
---------------------

**control**: 5-day control runs, started from the above steady-state histories at 
both equinoxes and both solstices, and at solar minimum and maximum conditions.

* "Sanity check" plots for solar min control runs 
  (global maps on last day of 5-day run): 

 * :ref:`December Solstice, Solar Minimum <control_decsol_smin>`
 * :ref:`June Solstice, Solar Minimum <control_junsol_smin>`
 * :ref:`March Equinox, Solar Minimum <control_mareqx_smin>`
 * :ref:`September Equinox, Solar Minimum <control_sepeqx_smin>`

* "Sanity check" plots for solar max control runs 
  (global maps on last day of 5-day run): 

 * :ref:`December Solstice, Solar Maximum <control_decsol_smax>`
 * :ref:`June Solstice, Solar Maximum <control_junsol_smax>`
 * :ref:`March Equinox, Solar Maximum <control_mareqx_smax>`
 * :ref:`September Equinox, Solar Maximum <control_sepeqx_smax>`

* `History file output for control runs <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/control>`_.

Climatology
-----------

**climatology**: Full-year Climatology with constant solar forcing:

* Heelis potential model with constant solar forcing::

   POWER   = 18.
   CTPOTEN = 30.
   F107    = 100.
   F107A   = 100.

* :ref:`"Sanity check" plots for Climatology <climatology>`
* `History file output for climatology run 
  <http://download.hao.ucar.edu/pub/tgcm/data/tiegcm1.94/climatology>`_.

December, 2006 "AGU Storm"
--------------------------

**dec2006**: December, 2006 "AGU" storm case:

* Heelis potential model with GPI (Kp) data
* Weimer potential model with IMF data (F10.7 from GPI)
* "Sanity check" plots for dec2006:

 * :ref:`Heelis/GPI (daily) <dec2006_heelis_daily>`
 * :ref:`Heelis/GPI (hourly) <dec2006_heelis_hourly>`
 * :ref:`Weimer/IMF (daily) <dec2006_weimer_daily>`
 * :ref:`Weimer/IMF (hourly) <dec2006_weimer_hourly>`

* `History file output for dec2006 simulation
  <http://download.hao.ucar.edu/pub/tgcm/tiegcm1.94/dec2006>`_.

November, 2003 Storm Case
-------------------------

**nov2003**: November 19-24 (days 323-328), 2003 storm case:

* Heelis potential model with GPI (Kp) data
* Weimer potential model with IMF data (F10.7 from GPI)
* "Sanity check" plots for nov2003:

 * :ref:`Heelis/GPI (hourly) <nov2003_heelis_hourly>`
 * :ref:`Weimer/IMF (hourly) <nov2003_weimer_hourly>`

* `History file output for nov2003 simulation 
  <http://download.hao.ucar.edu/pub/tgcm/tiegcm1.94/nov2003>`_.

Whole Heliosphere Interval (WHI)
--------------------------------

**whi2008**: Whole Heliosphere interval (WHI) (March 21 to April 16, 2008)

* Heelis potential model with GPI (Kp) data
* Weimer potential model with IMF data (F10.7 from GPI)
* "Sanity check" plots for whi2008:

 * :ref:`Heelis/GPI (daily) <whi2008_heelis_daily>`
 * :ref:`Weimer/IMF (daily) <whi2008_weimer_daily>`

* `History file output for whi2008 simulation
  <http://download.hao.ucar.edu/pub/tgcm/tiegcm1.94/whi2008>`_.


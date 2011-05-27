
Release Notes for TIEGCM v\ |version|
=====================================

New validated features in the public release
--------------------------------------------

* Weimer 2005 electric potential model, optionally with IMF data.
* Dynamic critical cross-over latitudes for transition between
  dynamo and empirical models at high-latitude.
* For a more detailed description of the use of the Weimer 2005 model in
  the TIEGCM, and the dynamic cross-over latitudes, please see 
  :ref:`Notes on Weimer05 in TIEGCM <tiegcm_weimer05>`.
* New diagnostics module (diags.F).  22 "sanctioned" diagnostic fields
  are available, and can be saved on secondary histories
  (see diags.table in doc directory).
* HPSS archive script can be saved by the model, to be executed after
  a run to save history files to the NCAR HPSS with "contents" annotations.

New Build System (scripts/ directory)
-------------------------------------

* Make.machine files containing platform and machine-specific compiler
  and link options, library locations, etc. are available in scripts directory.
* New simplified Makefile without platform-dependent conditionals.
* Support for Intel ifort/mpif90 compiler on 64-bit HAO Linux systems
  (Intel-built code is significantly faster than PGI-built code on systems
  with Intel hardware)
* Write svnversion to history files.

Utilities (scripts/ directory)
------------------------------

* *tgcm_contents:* Print "contents" (annotations) of history files
* *tgcm_ncdump:*   Print ncdump of history files, with data for scalars and 1d vectors
* *tgcm_put:*      Archive history files on the NCAR HPSS storage system
* *mknamelist:*    Make a namelist input file
* *mkjob:*         Make a job script

Validation and Benchmark Runs (tests/ directory)
------------------------------------------------

* Namelist input files and job scripts for benchmark runs are available
  in the tests/ directory.
* Seasonal steady-state source histories are provided for both equinoxes and 
  solstices, at solar minimum and solar maximum conditions.

* The following benchmark runs were made with tiegcm v\ |version|:

  * **Control:** 5-day Control Runs are made, starting from the steady-state 
    seasonal source histories.
  * **Climatology:** Full-year Climatology with constant solar forcing.
  * **Dec2006:** December, 2006 "AGU" storm case.
  * **Nov2003:** November 19-24 (days 323-328), 2003 storm case.
  * **Whi2008:** Whole Heliosphere Interval (WHI), March 21 to April 16, 2008.

* For more detailed information and access to the history files and
  post-processing, please see :ref:`Benchmark/test runs <tests>`.
 
Changes to namelist inputs
--------------------------

* Weimer potential model with IMF data from IMF_NCFILE, and F10.7 data
  from GPI_NCFILE.
* HPSS_PATH specifies destination path to the NCAR HPSS, for history file
  archive script. (All references to the old NCAR MSS have been removed)
* Time-dependent namelist option for F107a and F107d
* The following namelist parameters are deprecated: SAVE, SECSAVE, DISPOSE

Changes in the documentation directory doc/
-------------------------------------------

* All new User's Guide (pdf, html) is now provided in the doc directory.
* Model Description (pdf) also provided in the doc directory.
* README.download: Information and instructions for downloading and
  executing a default run with the current version.
* Release_Notes: Release notes for the current version.

Bug Fixes
---------

* Change nmlon to nlonp1 where defining dynpot, line 292 of magfield.F
* Change op(k,i) to nop(k,i) in calculation of ion chemistry heating (qic).
  Line 119 of qjion.F

Other
-----

* Changed units of f107d, f107a, ctpoten, and hpower on the histories.
* Divide by avogadro's number instead of multiplying by 1.66e-24 in
  sub calczg (addiag.F).


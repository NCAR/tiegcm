
.. _source_section:

Model Source Code
=================

The source code is in the :term:`src/` subdirectory of the model root directory 
(:term:`modeldir`), which is provided in the model :ref:`download <download>` file.


The Academic License Agreement
------------------------------

The TIEGCM :download:`Open Source Academic Research License Agreement <_static/tiegcmlicense.txt>`
specifies the terms and restrictions under which the NCAR/UCAR grants permission to use the
model, including the source code, for research, academic, and non-profit purposes. 

Source Code Flow Diagram
------------------------

A detailed flow diagram and calling tree of the source code structure is available
in single and multi-page pdf files:

.. Warning::

  Some details of these flow charts may be out of date with respect to TIEGCM version |version|

* :base_url:`TIEGCM Code Structure (multi-page pdf) <tiegcm_codestruct.pdf>`

* :base_url:`TIEGCM Code Structure (single-page pdf) <tiegcm_code_poster.pdf>`

.. _resolution:

Grid Structure and Resolution
-----------------------------

The TIEGCM can be configured for two spatial/temporal resolutions (use the :ref:`modelres <modelres>`
shell variable in the :ref:`job script <jobscript>` to set the model resolution):

 * | 5 degrees lat x lon, 2 grid points per scale height 
   | (default time step = 120 secs)
 * | 2.5 degrees lat x lon, 4 grid points per scale height 
   | (default time step = 60 secs)
 * The vertical coordinate ``lev``, or ``Zp``, is a log-pressure scale ``ln(p0/p)``, where p 
   is pressure and p0 is a reference pressure. Fields are calculated at either "interface" levels 
   (``ilev``), or at "midpoint" levels (``lev``) (see lev and ilev coordinate definitions below).
 
   * At 5.0 degree horizontal, Zp at interfaces = -7 to +5 by 0.5
   * At 2.5 degree horizontal, Zp at interfaces = -7 to +5 by 0.25

   .. note::
      To interpolate model fields to constant height surfaces, you should use
      geometric height, which is available on the 3d model grid as "ZG" on secondary 
      histories.

The spatial coordinates at the 5-degree resolution are defined as follows::

  double lon(lon) ;
    lon:long_name = "geographic longitude (-west, +east)" ;
    lon:units = "degrees_east" ;

  lon = -180, -175, -170, -165, -160, -155, -150, -145, -140, -135, -130, 
    -125, -120, -115, -110, -105, -100, -95, -90, -85, -80, -75, -70, -65, 
    -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 
    20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 
    110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175 ;


  double lat(lat) ;
    lat:long_name = "geographic latitude (-south, +north)" ;
    lat:units = "degrees_north" ;

  lat = -87.5, -82.5, -77.5, -72.5, -67.5, -62.5, -57.5, -52.5, -47.5, -42.5, 
    -37.5, -32.5, -27.5, -22.5, -17.5, -12.5, -7.5, -2.5, 2.5, 7.5, 12.5, 
    17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 
    77.5, 82.5, 87.5 ;

  double lev(lev) ;
    lev:long_name = "midpoint levels" ;
    lev:short_name = "ln(p0/p)" ;
    lev:units = "" ;
    lev:positive = "up" ;
    lev:standard_name = "atmosphere_ln_pressure_coordinate" ;
    lev:formula_terms = "p0: p0 lev: lev" ;
    lev:formula = "p(k) = p0 * exp(-lev(k))" ;

  lev = -6.75, -6.25, -5.75, -5.25, -4.75, -4.25, -3.75, -3.25, -2.75, -2.25, 
    -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 
    3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25 ;

  double ilev(ilev) ;
    ilev:long_name = "interface levels" ;
    ilev:short_name = "ln(p0/p)" ;
    ilev:units = "" ;
    ilev:positive = "up" ;
    ilev:standard_name = "atmosphere_ln_pressure_coordinate" ;
    ilev:formula_terms = "p0: p0 lev: ilev" ;
    ilev:formula = "p(k) = p0 * exp(-ilev(k))" ;

  ilev = -7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 
    0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7 ;

.. note::
   The 2.5 degree configuration ("double-resolution"), is not fully tuned and validated
   in TIEGCM version |version|.

.. _modifying_source:

Modifying the Source Code
-------------------------

As a community user, student, research scientist or developer, you may need to modify the model
source code. It is best to do this after building and at least making a default execution 
of the model (see the :ref:`QuickStart <quickstart>` Section). To change one or more 
source files, simply go to the :term:`src/` subdirectory in the model root directory
:term:`modeldir`, and edit the files as necessary. Then return to the working directory 
:term:`workdir` and re-execute the job script. It will recompile the modified files, and 
any other source files that depend on the modified files, and re-execute the model. 
Alternatively, you can enter the execution directory :term:`execdir`, and recompile 
the code by typing "gmake" on the command line, then return to the working directory 
and re-execute the job script.


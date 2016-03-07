
.. _grid:
.. _resolution:

Grid Structure and Resolution
=============================

The TIEGCM can be configured for two spatial/temporal resolutions (use the :ref:`modelres <modelres>`
shell variable in the :ref:`job script <jobscript>` to set the model resolution):

 * | 5 degrees lat x lon, 2 grid points per scale height 
   | (default time step = 60 secs)
 * | 2.5 degrees lat x lon, 4 grid points per scale height 
   | (default time step = 30 secs)
 * The vertical coordinate ``lev``, or :term:`Zp`, is a log-pressure scale ``ln(p0/p)``, where p 
   is pressure and p0 is a reference pressure. Fields are calculated at either "interface" levels 
   (``ilev``), or at "midpoint" levels (``lev``) (see lev and ilev coordinate definitions below).
 
   * At 5.0 degree horizontal, Zp at interfaces = -7 to +7 by 0.5
   * At 2.5 degree horizontal, Zp at interfaces = -7 to +7 by 0.25

   .. note::
      To interpolate model fields to constant height surfaces, you should use
      geometric height, which is available on the 3d model grid as "ZG" on secondary 
      histories. See the section below on :ref:`Altitude Coordinates the NCAR TIEGCM <altcoords>` 
      for a detailed explanation of the relationship between Zp and Altitude.

.. _geocoords:

Geographic 3d spatial coordinates at 5-degree horizontal resolution 
-------------------------------------------------------------------

Following are spatial coordinates for the 5x5-degree latxlon horizontal
grid, with two grid points per scale height in the vertical (delta :term:`Zp` = 0.5)::

  dimensions:
        lon = 72 ;
        lat = 36 ;
        lev = 29 ;
        ilev = 29 ;

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

Geographic 3d spatial coordinates at 2.5-degree resolution 
----------------------------------------------------------

Following are spatial coordinates for the 2.5x2.5-degree latxlon horizontal
grid, with four grid points per scale height in the vertical (delta :term:`Zp` = 0.25)::

 dimensions:
	lon = 144 ;
	lat = 72 ;
	lev = 57 ;
	ilev = 57 ;

	double lon(lon) ;
		lon:long_name = "geographic longitude (-west, +east)" ;
		lon:units = "degrees_east" ;

 lon = -180, -177.5, -175, -172.5, -170, -167.5, -165, -162.5, -160, -157.5, 
    -155, -152.5, -150, -147.5, -145, -142.5, -140, -137.5, -135, -132.5, 
    -130, -127.5, -125, -122.5, -120, -117.5, -115, -112.5, -110, -107.5, 
    -105, -102.5, -100, -97.5, -95, -92.5, -90, -87.5, -85, -82.5, -80, 
    -77.5, -75, -72.5, -70, -67.5, -65, -62.5, -60, -57.5, -55, -52.5, -50, 
    -47.5, -45, -42.5, -40, -37.5, -35, -32.5, -30, -27.5, -25, -22.5, -20, 
    -17.5, -15, -12.5, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5, 15, 
    17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 
    52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70, 72.5, 75, 77.5, 80, 82.5, 85, 
    87.5, 90, 92.5, 95, 97.5, 100, 102.5, 105, 107.5, 110, 112.5, 115, 117.5, 
    120, 122.5, 125, 127.5, 130, 132.5, 135, 137.5, 140, 142.5, 145, 147.5, 
    150, 152.5, 155, 157.5, 160, 162.5, 165, 167.5, 170, 172.5, 175, 177.5 ;

	double lat(lat) ;
		lat:long_name = "geographic latitude (-south, +north)" ;
		lat:units = "degrees_north" ;

 lat = -88.75, -86.25, -83.75, -81.25, -78.75, -76.25, -73.75, -71.25, 
    -68.75, -66.25, -63.75, -61.25, -58.75, -56.25, -53.75, -51.25, -48.75, 
    -46.25, -43.75, -41.25, -38.75, -36.25, -33.75, -31.25, -28.75, -26.25, 
    -23.75, -21.25, -18.75, -16.25, -13.75, -11.25, -8.75, -6.25, -3.75, 
    -1.25, 1.25, 3.75, 6.25, 8.75, 11.25, 13.75, 16.25, 18.75, 21.25, 23.75, 
    26.25, 28.75, 31.25, 33.75, 36.25, 38.75, 41.25, 43.75, 46.25, 48.75, 
    51.25, 53.75, 56.25, 58.75, 61.25, 63.75, 66.25, 68.75, 71.25, 73.75, 
    76.25, 78.75, 81.25, 83.75, 86.25, 88.75 ;

	double lev(lev) ;
		lev:long_name = "midpoint levels" ;
		lev:short_name = "ln(p0/p)" ;
		lev:units = "" ;
		lev:positive = "up" ;
		lev:standard_name = "atmosphere_ln_pressure_coordinate" ;
		lev:formula_terms = "p0: p0 lev: lev" ;
		lev:formula = "p(k) = p0 * exp(-lev(k))" ;

 lev = -6.875, -6.625, -6.375, -6.125, -5.875, -5.625, -5.375, -5.125, 
    -4.875, -4.625, -4.375, -4.125, -3.875, -3.625, -3.375, -3.125, -2.875, 
    -2.625, -2.375, -2.125, -1.875, -1.625, -1.375, -1.125, -0.875, -0.625, 
    -0.375, -0.125, 0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 
    2.125, 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125, 4.375, 
    4.625, 4.875, 5.125, 5.375, 5.625, 5.875, 6.125, 6.375, 6.625, 6.875, 
    7.125 ;

	double ilev(ilev) ;
		ilev:long_name = "interface levels" ;
		ilev:short_name = "ln(p0/p)" ;
		ilev:units = "" ;
		ilev:positive = "up" ;
		ilev:standard_name = "atmosphere_ln_pressure_coordinate" ;
		ilev:formula_terms = "p0: p0 lev: ilev" ;
		ilev:formula = "p(k) = p0 * exp(-ilev(k))" ;
 ilev = -7, -6.75, -6.5, -6.25, -6, -5.75, -5.5, -5.25, -5, -4.75, -4.5, 
    -4.25, -4, -3.75, -3.5, -3.25, -3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, 
    -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 
    2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 
    5.75, 6, 6.25, 6.5, 6.75, 7 ;

.. _magcoords:

Geomagnetic 3d spatial coordinates
----------------------------------

The longitude geomagnetic coordinate is from -180 to +180 by 4.5 degrees.
The latitude coordinate is non-regular, with resolution increasing toward
the magnetic equator. The vertical :term:`Zp` (ln(p0/p)) interface coordinate 
is from -8.5 to 7 by 0.25::

 dimensions:
	mlon = 81 ;
	mlat = 97 ;
	mlev = 63 ;
	imlev = 63 ;

 mlon = -180, -175.5, -171, -166.5, -162, -157.5, -153, -148.5, -144, -139.5, 
    -135, -130.5, -126, -121.5, -117, -112.5, -108, -103.5, -99, -94.5, -90, 
    -85.5, -81, -76.5, -72, -67.5, -63, -58.5, -54, -49.5, -45, -40.5, -36, 
    -31.5, -27, -22.5, -18, -13.5, -9, -4.5, 0, 4.5, 9, 13.5, 18, 22.5, 27, 
    31.5, 36, 40.5, 45, 49.5, 54, 58.5, 63, 67.5, 72, 76.5, 81, 85.5, 90, 
    94.5, 99, 103.5, 108, 112.5, 117, 121.5, 126, 130.5, 135, 139.5, 144, 
    148.5, 153, 157.5, 162, 166.5, 171, 175.5, 180 ;

 mlat = -90, -88.1238292398491, -86.2386359278657, -84.3344382773342, 
    -82.4013318763435, -80.4295344892688, -78.4094552099168, 
    -76.331796630125, -74.1876988925388, -71.9689341802758, 
    -69.6681589022773, -67.2792279882741, -64.7975706790533, 
    -62.2206194320588, -59.5482728298363, -56.7833601290164, 
    -53.9320608459732, -51.0042204168578, -48.0134966005524, 
    -44.9772754602266, -41.916313892128, -38.8540980954293, 
    -35.8159497801506, -32.8279553674349, -29.9158266703621, 
    -27.1038148776609, -24.4137889090065, -21.8645574169981, 
    -19.4714697638694, -17.2462861630082, -15.1972697734841, 
    -13.3294282264571, -11.6448185129562, -10.142824406667, 
    -8.82031765103987, -7.67162666281269, -6.68827297583048, 
    -5.85851734698832, -5.16689314460211, -4.5940469432968, 
    -4.11722526306697, -3.71151170575937, -3.35148255039153, 
    -3.01257883277328, -2.67136426606314, -2.3036287214954, 
    -1.87754943767857, -1.32687203939232, -7.72840966450717e-08, 
    1.32687203939232, 1.87754943767857, 2.3036287214954, 2.67136426606314, 
    3.01257883277328, 3.35148255039153, 3.71151170575936, 4.11722526306697, 
    4.59404694329679, 5.16689314460211, 5.85851734698832, 6.68827297583048, 
    7.67162666281268, 8.82031765103987, 10.142824406667, 11.6448185129562, 
    13.3294282264571, 15.1972697734841, 17.2462861630082, 19.4714697638694, 
    21.8645574169981, 24.4137889090064, 27.1038148776609, 29.9158266703621, 
    32.8279553674348, 35.8159497801506, 38.8540980954293, 41.916313892128, 
    44.9772754602266, 48.0134966005524, 51.0042204168578, 53.9320608459731, 
    56.7833601290163, 59.5482728298363, 62.2206194320588, 64.7975706790533, 
    67.2792279882741, 69.6681589022773, 71.9689341802758, 74.1876988925387, 
    76.331796630125, 78.4094552099168, 80.4295344892687, 82.4013318763434, 
    84.3344382773342, 86.2386359278657, 88.123829239849, 90 ;

 mlev = -8.25, -8, -7.75, -7.5, -7.25, -7, -6.75, -6.5, -6.25, -6, -5.75, 
    -5.5, -5.25, -5, -4.75, -4.5, -4.25, -4, -3.75, -3.5, -3.25, -3, -2.75, 
    -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 0.25, 
    0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 
    4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25 ;

 imlev = -8.5, -8.25, -8, -7.75, -7.5, -7.25, -7, -6.75, -6.5, -6.25, -6, 
    -5.75, -5.5, -5.25, -5, -4.75, -4.5, -4.25, -4, -3.75, -3.5, -3.25, -3, 
    -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, -0.75, -0.5, -0.25, 0, 
    0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 
    3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7 ;


.. _altcoords:

Altitude Coordinates in the NCAR TIE-GCM and TIME-GCM
-----------------------------------------------------

Author: Stan Solomon 
Date:   April, 2016

The purpose of this document is to define the altitude coordinate systems 
used in the NCAR Thermosphere-Ionosphere-Electrodynamics General Circulation Model 
(TIE-GCM) and  Thermosphere-Ionosphere-Mesosphere-Electrodynamics General Circulation 
Model (TIME-GCM), especially to inform model users as to how to register model output 
in the vertical dimension.

The TIE-GCM and TIME-GCM use a log-pressure coordinate system, with each pressure 
level defined as ln(P0/P), where P0 = 5x10-4 dynes/cm2 = 5x10-5 Pascal = 5x10-7 
hPa = 5x10-7 mb.  (Native units in these models are cgs, i.e., dynes/cm2.) 
This pressure occurs at ~200 km altitude, depending on conditions. 

The TIE-GCM vertical coordinate extends from -7 to +7 (~97 km to ~600 km) and the 
TIME-GCM vertical coordinate extends from -17 to +7 (~30 km to ~600 km).  Each integer 
interval in pressure level is one scale height apart, so the low-resolution (5°x5°xH/2) 
versions are spaced at half-integer intervals and the high-resolution (2.5°x2.5°xH/4) 
versions of the models are spaced at quarter-integer intervals:

=================  ==========  =============  ============  =========  =======  =======
Model/Resolution   Num Levels  Level Spacing  Bottom Level  Top Level  Min Alt  Max Alt
=================  ==========  =============  ============  =========  =======  =======
Low-Res TIE-GCM    29          0.5            -7            +7         97 km    600 km
High-Res TIE-GCM   57          0.25           -7            +7         97 km    600 km
Low-Res TIME-GCM   49          0.5            -17           +7         30 km    600 km
High-Res TIME-GCM  97          0.25           -17           +7         30 km    600 km
=================  ==========  =============  ============  =========  =======  =======

The height of the pressure surface is defined at each grid point in arrays provided in 
output history files (in cm).  Unfortunately, there are four different possibilities for 
altitude definition, all slightly different.

First, we define the geopotential height z.  Geopotential height is the height that the 
pressure surface would be, assuming that the acceleration due to gravity g is constant 
at the value used in the model calculations (870 cm/s2 for the TIE-GCM and 950 cm/s2 for 
the TIME-GCM).  It is registered to the altitude of the model lower boundary, which can 
vary horizontally due to the tidal and climatological lower boundary specification.  
This is the native coordinate system for the models, and so z is included in all history 
files.  However, it is not the appropriate altitude coordinate for comparison with 
real-world data.  Also note that this definition of geopotential height is not the same 
as what is used in, e.g., tropospheric meteorology, because it is referenced to value of 
g that is different from the value of g at the surface (~980 cm/s2).

We can correct the geopotential height z to obtain geometric height zg.  This is performed 
inside the models by subroutine zgcalc, using an empirical formulation of the variation of 
g over the globe (including centripetal force), and vertical integration, to account for 
the variation with altitude.  It can also be done, using the same subroutine, in the 
Fortran model processers, and is also available in various IDL processing routines.  
Geometric height zg is now forced onto secondary histories (i.e., it is output whether you 
request it or not) but not on primary histories (because primary histories contain only 
what is necessary to re-start the model).  However, some older secondary history files 
may not include zg which necessitates that it be calculated in the post-processing if 
needed for data comparison.

Now we come to the final complication, which is the distinction between model interfaces 
and model mid-points.  The interfaces are the native coordinate system of the model grid, 
as defined in the table above, i.e., at -7.0, -6.5, -6.0, etc.;  z and zg are defined on 
these interfaces.  However, most model output quantities are actually reported at the 
midpoints, half-way between interfaces in pressure, i.e., at -6.75, -6.25, -5.75, etc.  
Each midpoint is a half-interval above the corresponding interface.  All temperatures, 
winds, neutral densities, etc., are defined at these midpoints.  However, electron 
density and electric potential are defined at the interfaces:

============= === === === === === === === === === === === === === === === === === ===
Field         Z   Zg  Zm  Tn  Un  Vn  O2  O   N2  NO  N   N2D He  Ne  Te  Ti  OM  Pot
============= === === === === === === === === === === === === === === === === === ===
Specified at  I   I   M   M   M   M   M   M   M   M   M   M   M   I   M   M   I   I
============= === === === === === === === === === === === === === === === === === ===

In order to register midpoint quantities in altitude, it is therefore necessary to 
interpolate from the midpoints to the interfaces.  Alternatively, it may be simpler 
to interpolate zg from the interfaces to the midpoints.  For TIE-GCM 2.0, a new output 
variable has been added, zm, which is geometric height that has been interpolated to 
the mid points.  However, older history files do not include zm.  As with zg, it is 
available on secondary histories but not on primary histories.

In output histories, quantities specified at interfaces are defined by the ilev 
coordinate variable and quantities specified at midpoints are defined by the lev 
coordinate variable.  These quantities are generally numerically identical, but their 
definitions in the files can serve as a reminder of what is defined where. 


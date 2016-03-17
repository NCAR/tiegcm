
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

Geographic 3d spatial coordinates at 5-degree resolution 
--------------------------------------------------------

Following are spatial coordinates for the 5x5-degree latxlon horizontal
grid, with two grid points per scale height in the vertical (delta :term:`Zp` = 0.5):

 * Dimensions:

   * lon = 72 ;
   * lat = 36 ;
   * lev = 29 ;
   * ilev = 29 ;

 * Coordinates:

   * lon = -180W to +180E by 5 degrees
   * lat = -87.5S to +87.5N by 5 degrees
   * lev = -6.75 to +7.25 by 0.5 
   * ilev = -7 to +7 by 0.5

Geographic 3d spatial coordinates at 2.5-degree resolution 
----------------------------------------------------------

Following are spatial coordinates for the 2.5x2.5-degree latxlon horizontal
grid, with four grid points per scale height in the vertical (delta :term:`Zp` = 0.25):

 * Dimensions:

   * lon = 144 ;
   * lat = 72 ;
   * lev = 57 ;
   * ilev = 57 ;

 * Coordinates:

   * lon = -180W to 177.5E by 2.5 degrees
   * lat = -88.75S to 88.75N by 2.5 degrees
   * lev = -6.875 to 7.125 by 0.25 
   * ilev = -7 to +7 by 0.25

.. _magcoords:

Geomagnetic 3d spatial coordinates
----------------------------------

The longitude geomagnetic coordinate is from -180 to +180 by 4.5 degrees.
The latitude coordinate is non-regular, with resolution increasing toward
the magnetic equator. The vertical :term:`Zp` (ln(p0/p)) interface coordinate 
is from -8.5 to 7 by 0.25:

 * Dimensions:

   * mlon = 81 ;
   * mlat = 97 ;
   * mlev = 63 ;
   * imlev = 63 ;

 * Coordinates:

   * mlon = -180W to 180E by 4.5 degrees
   * mlat = -90S to 90N: irregular, increasing resolution equatorward
   * mlev = -8.25 to 7.25 by 0.25
   * imlev = -8.5 to 7.0 by 0.25

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


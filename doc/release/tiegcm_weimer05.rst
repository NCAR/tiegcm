
.. _tiegcm_weimer05:

Notes on the use of the Weimer 2005 model in TIEGCM
===================================================

Description of the Weimer 2005 model in TIEGCM v1.94, by Barbara Emery and Dan Weimer (May 26, 2011)
----------------------------------------------------------------------------------------------------

The TIEGCM uses an imposed electric field model such as the Weimer (2005) or Heelis (1982) electric 
potential models at high magnetic latitudes.  In version 1.94 of the TIEGCM, the Weimer 2005 model 
can now be used with or without optional IMF files of 15-min averages where t=-20 to -5 minute values 
are used for time t=0.  Sometimes the filled values are not reasonable, especially for large gaps in 
the IMF (e.g. during the Halloween storm, 03301.6-03304.3, Oct 28-31, 2003).  The convection radius 
increases from 11 to 19 degrees from Kp 0 to 9 in the Heelis model, but can reach 25 or more degrees 
in the Weimer 2005 model, especially with large solar wind velocities.  There were a very small number 
of high velocity cases used in the construction of the Weimer 2001 and 2005 models.  Thus the electric 
potential patterns for large solar wind velocities (>900 km/s, but especially ~1100 km/s) can be 
unrealistic.   The data used in the construction of the models also had very few cases with IMF magnitudes 
greater than 15 nT.  Through the use of a saturation curve, the 2005 model appears to be realistic 
(ie, lower potential drops consistent with observations of saturation) at high IMF magnitudes (> 20 nT), 
although results may not be as accurate as at lower IMF levels.


Dynamic critical cross-over latitudes
-------------------------------------

The TIEGCM uses an imposed electric field model such as the Weimer (2005) or Heelis (1982) electric 
potential models at high magnetic latitudes above a critical magnetic co-latitude called crit(1).  
Before TIEGCM version 1.94, this was set at 15 degrees co-latitude, or at 75 mlat.  Below a second 
magnetic co-latitude called crit(2), the dynamo model calculates the electric field in the TIEGCM.  
This was set at 30 degrees co-latitude, or at 60 mlat.  Between 60 and 75 mlat, there was a linear 
variation of the combined dynamo solution and the imposed high latitude solution.   

In TIEGCM version 1.94, crit(1,2) were revised to change dynamically as a function of the convection 
radius such that crit(1) is 5 degrees larger than the convection radius and crit(2) is set to 
crit(1)+15 degrees.  This ensures that the high-latitude electric field model is used in the TIEGCM 
for the intense sunward flows equatorwards of the convection reversal boundary.  Tests showed that 
vertical ion drifts at the equator, TEC, and NmF2 values on active days are generally improved in 
both the Weimer (2005) and Heelis (1982) convection cases using dynamic crit values.  



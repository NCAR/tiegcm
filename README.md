Changes of the current version compared to public TIEGCM v2\.0 \(2016\):

- The job script is modified to support arbitrary combination of horizontal and vertical resolutions, change of magnetic grids is also supported\.
- The job script, defs\.h and initialization of some altitude\-dependent variables \(xfac/ar\_glbm/aureff/bdriz\) are rewritten to support the extension of the upper boundary\.
- The modeltime of input/output is set to 4 digits \(day/hour/minute/second\) instead of the old 3\-digit format \(day/hour/minute\) to allow a higher cadence\.
- The calculation of N2 mixing ratio, mean molecular mass and scale height is unified in the model to avoid recalculation\. Some artificial caps of N2 mixing ratio are removed\.
- Rewrite the Helium module and it is now included with in all resolutions \(default on\)\.
- The old Fourier filter is replaced with the new ring filter\.
- An additional input parameter NSTEP\_SUB is introduced to control the number of O\+ sub\-cyling\.
- NetCDF4 parallel IO is turned on to reduce the memory usage on root task\.
- ESMF calls are modified to ensure bit\-for\-bit reproducibility\.
- IGRF is updated to the newest version\.
- Magnetospheric coupling module is rewritten to support in-memory MPI data transfer\.
- Simplify the code by removing unused parameters/argmuments/variables in some functions\.
- Some MPI subroutines are rewritten to allow a speed boost\.
- dipmin is set to sin\(dlat\*2\*dtr\) instead of manual setup in different resolutions\.
- Some minor bug fixes\.

Change of physics:
- Modified the coefficients of solar heating \(Astrid Maute\)\.
- Add a scaling factor in accounting for the height variation of equatorward electric field \(elam\) \(Astrid Maute\)\.
- Include the field\-aligned ion drag in the momentum equation \(Jiuhou Lei\)\.
- Change the parameterization scheme of electron heat flux \(fed\) near the equator in settei \(Tong Dang\)\.
- Use a sixth\-order polynomial for thermal electron heating efficiency \(Yihui Cai\)\.
- Electrojet turbulent heating, default off \(Jing Liu\)\.
- Empirical SAPS, default off \(Cheng Sheng\)\.
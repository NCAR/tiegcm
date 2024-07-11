# TIEGCM v3.0

Documentation: [TIEGCM ReadtheDocs](https://tiegcm-docs.readthedocs.io/en/latest/) 

## Code Structures

| Subdirectory | Description           | Summary of Contents                                                |
|--------------|-----------------------|--------------------------------------------------------------------|
| scripts/     | Support scripts       | Job scripts, Make files, Utilities                                 |
| src/         | Source code           | Source files *.F, *.F90, *.h                                       |
| tiegcmrun/   | TIEGCM User Interface | Python program for pre-processes, build routines and job execution |

For a brief set of instructions to build the model and make a short default run, see [documentation](https://tiegcm-docs.readthedocs.io/en/latest/).

Please also see the main TGCM website: [TGCM Website](http://www.hao.ucar.edu/modeling/tgcm)

Additional data may be available on the HAO public FTP site: [HAO Public FTP Site](http://download.hao.ucar.edu/pub/tgcm)

User's Guide, Model Description, and Release documentation are available from the on [TIEGCM ReadtheDocs](https://tiegcm-docs.readthedocs.io/en/latest/).

For any questions or further information, please contact the discussion group email list at tgcmgroup@ucar.edu.

This is a summary of modifications made to the TIEGCM since the release of TIEGCM 2.0 (March 2016).

## New Features and Functional Changes

### Flexible Resolutions
The job script is modified to support arbitrary combinations of horizontal and vertical resolutions. Changes to magnetic grids are also supported.

### Extended Upper Boundary
Job script, defs.h, and initialization of some altitude-dependent variables (xfac, ar_glbm, aureff, bdriz) are rewritten to support the extension of the upper boundary.

### High Cadence Model Time
Model time for input/output is set to 4 digits (day/hour/minute/second) instead of the old 3-digit format (day/hour/minute) to allow higher cadence.

### Unified N2/MBAR/SCHT Calculations
The calculation of N2 mixing ratio, mean molecular mass, and scale height is unified to avoid recalculations. Some artificial caps on N2 mixing ratio are removed.

### Rewritten Helium Module
The Helium module is rewritten and now included with all resolutions (default on). In addition, Helium effects on heating rates, etc., are now accounted correctly throughout the code.

### Ring Filter
The old Fourier filter is replaced with a new ring filter.

### O+ Sub-Cycling
An additional input parameter NSTEP_SUB is introduced to control the number of O+ sub-cycling.

### NetCDF4 Parallel IO
NetCDF4 parallel IO is enabled to reduce memory usage on the root task.

### Bit-for-Bit Reproducibility
ESMF calls are modified to ensure bit-for-bit reproducibility.

### Updated IGRF
The geomagnetic field is updated to the latest International Geomagnetic Reference Field version.

### Rewritten Magnetospheric Coupling Module
The magnetospheric coupling module is rewritten to support in-memory MPI data transfer.

### Code Simplifications
Unused parameters, arguments, and variables are removed from some functions to simplify the code.

### MPI Subroutine Optimizations
Some MPI subroutines are rewritten for a speed boost.

### Consistent Dipmin Calculation
Dipmin is set to sin(dlat*2*dtr) instead of being manually set for different resolutions.

### Miscellaneous Bug Fixes
Various minor bug fixes.

## Changes in Physics

### Solar Heating Coefficients
Modified coefficients of solar heating (Astrid Maute).

### Height Variation of Equatorward Electric Field
A scaling factor is added to account for the height variation of the equatorward electric field (elam) (Astrid Maute).

### Field-Aligned Ion Drag
Field-aligned ion drag is included in the momentum equation (Jiuhou Lei).

### Collision Frequency
In lamdas.F, collision frequency calculation now include all ion species (O+, O2+, N+, N2+ and NO+) instead of only accounting for O+, O2+, and NO+ in previous versions (Haonan Wu).

### N(2D) Transport
Minor species solver now includes N(2D) which was assumed in (photo)-chemical equilibrium, this affects N chemistry at very high altitudes (z>7) (Haonan Wu).

### Electron Heat Flux Parameterization
The parameterization scheme of electron heat flux (fed) near the equator is changed in settei (Tong Dang, Wenbin Wang, Kevin Pham).

### O+ Number Flux Parameterization
The parameterization scheme of O+ (opflux) near the equator is changed in oplus (Haonan Wu, Wenbin Wang)

### Thermal Electron Heating Efficiency
A sixth-order polynomial is used for thermal electron heating efficiency (Yihui Cai).

### Electrojet Turbulent Heating
Electrojet turbulent heating is included, default off (Jing Liu).

### Empirical SAPS
Empirical SAPS is included, default off (Wenbin Wang).

### Eclipse Solar EUV Masking
Support for eclipse solar EUV masking is added (Tong Dang, Jiuhou Lei).

### Lower Boundary Forcing by External Data
Support for lower boundary forcing by external data (SD nudging) (Haonan Wu, Xian Lu).

## Utility Tools

### TIEGCMrun
Tiegcmrun is a Python tool (/tiegcmrun directory) that is used to compile and execute tiegcm in an automated fashion. Tiegcmrun can be executed interactively on the command line. See example of usage under [QuickStart](https://tiegcm-docs.readthedocs.io/en/latest/tiegcm/quickstart.html).

### [TIEGCMpy](https://tiegcmpy.readthedocs.io/en/latest/)
Tiegcmpy is a Python tool ([Tiegcmpy github](https://github.com/NCAR/tiegcmpy)) that is used for post processing and data visualization of TIEGCM outputs. Tiegcmpy can be executed interactively on the command line or as an API in a python script. See example of usage [TIEGCMpy Docs](https://tiegcmpy.readthedocs.io/en/latest/).


Directory Structure
===================

This section describes a typical directory structure the user will be
working with when running the TIEGCM model. The working directory
(:term:`workdir`) is the "root" directory of the user's project. 

The model directory (:term:`modeldir`) is typically a subdirectory under
the working directory, and contains the model source code, supporting
scripts, documentation, and scripts for running tests. 

The data directory (:term:`datadir`) may also be a subdirectory under
the working directory, or it may be on a large temporary disk that is
accessible from the working directory. The data directory contains 
start-up and input data files for running the model.

Working Directory (:term:`workdir`)
-----------------------------------

The user's working directory will typically look something like this
(the :term:`datadir` can be on a large separate disk system)::

                      workdir
                         |
 -----------------------------------------------
    |           |             |            |          
 *.inp       execdir/     modeldir/     datadir/
 *.job          |
 *.out        *.o
              *.mod
              *.nc
              Make*
              exec

Here, ``\*.inp`` are :term:`namelist input` files, ``\*.job`` are 
:term:`job script` s, and ``\*.out`` are are stdout files from model 
runs. The :term:`execdir` is the build/execution
directory (created by the first build), with object code (``\*.o``, ``\*.mod``), 
model :ref:`output history files <output>` (``\*.nc``), make files (``Make\*``), 
and an executable file (``exec``).  Various other files may also be in the execdir. 
The modeldir and datadir directories are described below.

The job script in your working directory contains a shell variable specifying
the :term:`modeldir`, so it knows where to find the source code and supporting scripts
for the build process. The namelist input file also refers to the :term:`datadir` 
for start-up and other data input files (e.g., :ref:`SOURCE <SOURCE>`, 
:ref:`GPI_NCFILE <GPI_NCFILE>`, :ref:`IMF_NCFILE <IMF_NCFILE>`, etc). 
These namelist parameters can use the environment variable :term:`TGCMDATA` to 
specify the :term:`datadir` (see section on :ref:`namelist input files <namelist>`).

Model Directory (:term:`modeldir`)
----------------------------------

The model root directory is what you get when you :ref:`download <download>` the 
model source code tar file. The model directory contains subdirectories with the 
model source code, supporting scripts, documentation, and test scripts::

                                modeldir
                                   |
   ----------------------------------------------------------------------
      |               |                |                 |
     src/          scripts/           doc/             tests/
      |               |                |                 | 
                    Make.*         userguide/         control/
     *.F          linux.job       description/      climatology/
     *.h           ibm.job          release/          dec2006
                 default.inp       diags.table        nov2003/
                tgcm_contents    README.download      whi2008/
                 tgcm_ncdump                           
                    etc                             

:term:`src/` directory contents:

* Fortran source code ``*.F``, ``*.h``. The source code is f90 compliant, and most 
  source files are in fixed-format fortran. There is a single header file, ``defs.h``,
  which contains grid definitions and dimensions.

:term:`scripts/` directory contents:

* **Make.\***: Makefiles containing platform-dependent compiler flags, 
  Make variables, and library locations. For example, 
  :download:`Make.intel_hao64 <_static/Make.intel_hao64>`. These file can
  be copied, renamed, and customized for the user's platform/machine environment.
* **tiegcm-linux.job**: Default model build/execute script for Linux systems.
* **tiegcm-ibm.job**: Default model build/execute script for IBM/AIX systems.
* **tiegcm_default.inp**: Default namelist input file.
* **tgcm_contents**: Utility script to print "contents" of netCDF output history files.
* **tgcm_ncdump**: Utility script to print an "ncdump" of history files, including
  data for scalars and 1-d vectors.

:term:`doc/` directory contents:

* **userguide/**: Directory containing source files for the User's Guide (this document)

* **description/**: Directory containing source files for the 
  :base_url:`TIEGCM Model Description <description/tiegcm_modeldes_6Oct09.pdf>`

* **release/**: Directory containing source files for the 
  :base_url:`Release Documentation <release/html>`

* **diags.table**: :download:`Table of diagnostic fields <_static/diags.table>` that can be 
  saved on secondary history files.

* **README.download**: :download:`Instructions <_static/README.download>` for how to make a 
  quick-start default build and execution of the model after downloading the source code and data.

:term:`tests/` directory contents:

* **README.tests**: Summary of benchmark test runs made with the current version of the model
* **Several directories** (climatology, control, dec2006, nov2003, etc) containing namelist
  input files and job scripts that can be used to reproduce the benchmark runs for
  validation and testing, comparing results from code changes, etc. 
* For more information on benchmark runs made for the current release, please see 
  :base_url:`Release Documentation <release/html>`

Data Directory (:term:`datadir`)
--------------------------------

The public TIEGCM data directory is what you get when you :ref:`download <download>` 
the data tar file. This directory is typically referred to with the environment variable
:term:`TGCMDATA`. Subsequently, after the data download, you may obtain additional needed 
data files from the :term:`NCAR Community Data Portal`. Here is a partial schematic of the 
datadir (where "tiegcmx.xx" is the desired model version)::

                       datadir
                          |
   -----------------------------------------------
              |                      |
          gpi*.nc                tiegcmx.xx/
          gswm*.nc                   |
          imf *.nc          TGCM.tiegcmx.xx.p*.nc
            etc                pcntr*smin*.nc
                               pcntr*smax*.nc 
                                   etc

Files listed on the left side refer to data input files that may be needed when 
running the model in different modes. These are netCDF files, specifically prepared 
for import into the TIEGCM model (they are **not** model start-up SOURCE files). 
These files are version-independent (can be used by (almost) any version of the model). 
They are usually provided to the model as namelist input parameters:

* **gpi\*.nc**
  GeoPhysical Indices data files (3-hourly Kp and F10.7 cm solar flux).
  Namelist Input parameter: :ref:`GPI_NCFILE <GPI_NCFILE>`
  
* **gswm\*.nc**
  `Global Scale Wave Model <http://www.hao.ucar.edu/modeling/gswm/gswm.html>`_
  data files, used to specify tidal perturbations for the lower boundary of
  the TIEGCM. There are 4 separate files for diurnal, semi-diurnal, migrating
  and non-migrating tides. For the namelist input parameters, please see :ref:`GSWM <GSWM>`.

* **imf\*.nc**
  Interplanetary Magnetic Field OMNI data files. Namelist read parameter is
  :ref:`IMF_NCFILE <IMF_NCFILE>`. These files contain data for the BX,BY,BZ 
  components of the IMF, solar wind velocity and solar wind density.

The "tiegcmx.xx/" subdirectory refers to the version of the model that was downloaded. 
This subdirectory contains start-up :ref:`SOURCE <SOURCE>` files from :term:`benchmark runs`
executed by that version of the model (see Section on :ref:`Benchmark Test Runs <tests>`).
These files can be used to remake the benchmark runs for testing and validation.
Here is an example of start-up files provided for benchmark runs made by TIEGCM version 
|version|::

  TGCM.tiegcm1.94.p_dec2006_heelis_gpi_001.nc
  TGCM.tiegcm1.94.p_dec2006_weimer_imf_001.nc
  TGCM.tiegcm1.94.p_nov2003_heelis_gpi_001.nc
  TGCM.tiegcm1.94.p_nov2003_weimer_imf_001.nc
  TGCM.tiegcm1.94.p_whi2008_heelis_gpi_001.nc
  TGCM.tiegcm1.94.p_whi2008_weimer_imf_001.nc
  TGCM.tiegcm1.94.pclim_heelis_001.nc
  pcntr_decsol_smax.nc
  pcntr_decsol_smin.nc
  pcntr_junsol_smax.nc
  pcntr_junsol_smin.nc
  pcntr_mareqx_smax.nc
  pcntr_mareqx_smin.nc
  pcntr_sepeqx_smax.nc
  pcntr_sepeqx_smin.nc


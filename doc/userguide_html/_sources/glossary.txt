
Glossary
--------

.. glossary::

   benchmark runs
     A series of :ref:`validation/test runs <tests>` made with each release 
     of the model.  Namelist files and job scripts used to make these benchmark 
     runs are available in subdirectories under the :term:`tests/` directory of the
     current release (e.g., namelist input files for the "control" test run are
     ``modeldir/tests/control/*.inp``).

   continuation run
     A continuation run continues from the last history of the previous run.
     That history (:ref:`START <START>` time) must be on the first :ref:`OUTPUT <OUTPUT>` file  
     provided by the namelist input file. A continuation run must not specify
     a :ref:`SOURCE <SOURCE>` file or :ref:`SOURCE_START <SOURCE_START>` time.

   datadir
     The directory containing start-up and other input data files required for
     running the model. A minimal set of datadir files are available via 
     :ref:`download <download>`.  The datadir is sometimes referred to by
     the :term:`TGCMDATA` environment variable.  Additional data files are available via the 
     :term:`NCAR Community Data Portal`.

   doc/
     Subdirectory under the :term:`modeldir` containing documentation, e.g., the
     User's Guide, Model Description, Release Notes, etc.

   execdir
     The model execution directory. This is the directory where the model is built 
     and executed. It is typically, but not necessarilly, a subdirectory of your
     working directory :term:`workdir`. When a job script is executed from a working 
     directory, the execdir is created if it does not already exist. During a model run, 
     output history files are written to the execdir.

   history
     A model history records the state of the model at a discrete instant in
     :term:`model time`. One or more histories are stored in netCDF history files.

   initial run
     An initial run is started from a history on a :ref:`SOURCE <SOURCE>` file (see also
     :ref:`SOURCE_START <SOURCE_START>`). Subsequent :term:`continuation runs <continuation run>` 
     do not provide SOURCE or SOURCE_START, but rather search for the :ref:`START <START>` time 
     on the first :ref:`OUTPUT <OUTPUT>` history file provided in the namelist input, and continue
     the run from there.

   job script
     A csh script in the scripts/ directory which, when executed, will build and execute
     the model. The user defines a few shell variables in the job script, such as
     the :term:`modeldir`, and the :term:`namelist input`. For more details, please
     see :ref:`job scripts <jobscript>`.

   model time
     TIEGCM model time is represented by an integer triplet: day,hour,minute, where 
     day is the julian day of the year, and hour is the ut. The variable for model time
     on history files is mtime(3,ntimes). For example, a history file may contain
     24 hourly histories for day 80: mtime = 80,1,0, 80,2,0, ... 81,0,0.
  
   namelist input
     The model reads user specified parameters from the :ref:`namelist input file <namelist>`
     via f90 standard namelist read. Keyword/Value pairs are read from unit 5,
     and are validated by the input module (input.F).
     
   NCAR Community Data Portal
     The `NCAR Community Data Portal <http://cdp.ucar.edu/>`_ is a public data 
     repository with datasets from NCAR, UCAR, UOP, and participating organizations. 
     To browse TIEGCM-related files (mostly netCDF history files for model start-up, 
     or results of :term:`benchmark runs`), click on the "Models" link, then to the
     "Thermospheric General Circulation Models" link, and finally to the desired
     model version. NetCDF Metadata is available without actually downloading files.

   netCDF
     TIEGCM output history files are written in 
     `netCDF <http://www.unidata.ucar.edu/software/netcdf/>`_, a self-describing 
     platform-independent data format written and maintained by the UCAR 
     `Unidata <http://www.unidata.ucar.edu>`_ program.

   resolution
     The TIEGCM can be run in one of two resolutions: 

       * 5   x 5   deg lat x lon, 2 grid levels per scale height (dz = 0.50)
       * 2.5 x 2.5 deg lat x lon, 4 grid levels per scale height (dz = 0.25)

     The resolution is set by the "modelres" shell variable in the TIEGCM 
     :ref:`job script <jobscript>`.

     .. note::

       The 2.5-degree resolution model is available in version |version|, but it is 
       not fully validated or supported by the public release.

   scripts/
     Subdirectory under the :term:`modeldir` containing supporting and utility 
     scripts, including job scripts, the default namelist input file, several
     Make files, etc.

   src/
     Subdirectory under the :term:`modeldir` containing the model source code
     (\*.F, \*.h files).

   tests/
     Subdirectory under the :term:`modeldir`. The tests directory 
     contains subdirectories for :term:`benchmark runs` that were made for 
     the current release.  The subdirectories contain job scripts and namelist input 
     files that can be used to reproduce benchmark runs for testing and validation 
     purposes. For more information, see the section on :ref:`Benchmark Test Runs <tests>`.

   TGCMDATA
     A unix environment variable that refers to the :term:`datadir`. This environment
     variable may be used when referring to data files in the namelist read file, e.g.,
     "GPI_NCFILE = `$TGCMDATA/gpi_xxxxx.nc`". See :ref:`namelist read files <namelist>`.

   modeldir
     The model root directory. This directory typically contains subdirectories
     :term:`src/` (model source code), :term:`scripts/` (utility scripts), 
     :term:`doc/` (documentation), and :term:`tests/` (test runs). The modeldir 
     is available via :ref:`download <download>`, and is typically a subdirectory 
     of the model working directory (:term:`workdir`). 
 
   workdir
     Your local working directory. This will typically contain the model root directory
     :term:`modeldir`, the execution directory :term:`execdir`, and related namelist
     input files, job scripts, stdout files, etc. It may also contain a data subdirectory
     :term:`datadir`.

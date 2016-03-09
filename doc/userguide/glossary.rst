
Glossary
--------

.. glossary::

   benchmark runs
     Selected validation runs made with each release of the model. These runs can
     be made using Python code in the :term:`tgcmrun/` directory. 

   benchmarks/
     Directory in the model root directory containing shell scripts that call
     :term:`tgcmrun/` for making benchmark runs, some utility scripts, and
     a subdirectory postproc/ containing scripts that do post-processing on
     benchmark results. See :ref:`Benchmark Runs <benchmarks>` for more information.
     Benchmark results (plots) for version |tgcm_version| are available here: 
     :base_url:`Release Benchmarks Results <release/html/benchmarks.html>`

   continuation run
     A continuation run continues from the last output history of the previous run.
     That history (:ref:`START <START>` time) must be on the first :ref:`OUTPUT <OUTPUT>` file  
     provided by the namelist input file. A continuation run must not specify
     a :ref:`SOURCE <SOURCE>` file or :ref:`SOURCE_START <SOURCE_START>` time.
     See also :ref:`Continuation Run <continuation_run>`

   diagnostic fields
     A list of diagnostic fields are available to be saved on secondary history files.
     See section :ref:`Saving Diagnostic Fields <diagnostics>`.

   datadir
     Directory containing startup history and data files necessary for running the model.
     This is specified with the :term:`tgcmdata` shell variable in the :term:`job script`.

   doc/
     Subdirectory under the :term:`modeldir` containing documentation, e.g., the
     User's Guide, Model Description, Release Notes, etc.

   ESMF
     "Earth System Modeling Framework". The ESMF library is used in the electro-dynamo 
     code (pdynamo.F in version |version| or later) for regridding between geographic 
     and geomagnetic grids in an mpi environment. This is open software that can
     be downloaded at https://www.earthsystemcog.org/projects/esmf/download/
     If you build the ESMF library, it should be built with the same compiler with 
     which the model is built.

   execdir
     The model execution directory. This is the directory where the model is built 
     and executed. It should be on a large temporary disk, capable of storing
     model object and module code, netCDF output history files, and other data. 
     When a job script is executed from a working directory, the execdir is created 
     if it does not already exist. During a model run, output history files are written 
     to the execdir. The execdir is set in the :term:`job script`. See also
     :ref:`Execution Directory <execdir>`

   geomagnetic coordinates
     The electro-dynamo fields (electric potential, electric field, and ion drift
     velocities) are calculated on a geomagnetic grid, see :ref:`magnetic coordinates <magcoords>` 

   Globus
     `Globus Data Sharing Service <https://www.globus.org>`_ for scientific research. 
     The |tgcm_version| benchmark history files and post-processing are available via Globus. 
     See :ref:`benchmark_history_files` for more information.

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
     the :term:`modeldir`, and the :term:`namelist input`. See example 
     :download:`job script for Linux desktops <../../scripts/tiegcm-linux.job>`, and
     :download:`job script for Super computer <../../scripts/tiegcm-ys.job>`
     See :ref:`jobscript` for more detailed information.

   model time
     TIEGCM model time is represented by an integer triplet: day,hour,minute, where 
     day is the julian day of the year, and hour is the ut. The variable for model time
     on history files is mtime(3,ntimes). For example, a history file may contain
     24 hourly histories for day 80: mtime = 80,1,0, 80,2,0, ... 81,0,0.

   modeldir
     The model root directory. This directory typically contains subdirectories
     :term:`src/` (model source code), :term:`scripts/` (utility scripts), 
     :term:`doc/` (documentation), and :term:`benchmarks/`. The modeldir 
     is available via :ref:`download <download>`, and is typically a subdirectory 
     of the model working directory (:term:`workdir`). 
     See also :ref:`Model Directory <modeldir>`
  
   namelist input
     The model reads user specified parameters from the :ref:`namelist input file <namelist>`
     via f90 standard namelist read. Keyword/Value pairs are read from unit 5,
     and are validated by the input module (input.F). See also :ref:`job scripts <jobscript>`.

   netCDF
     TIEGCM output history files are written in 
     `netCDF <http://www.unidata.ucar.edu/software/netcdf/>`_, a self-describing 
     platform-independent data format written and maintained by the UCAR 
     `Unidata <http://www.unidata.ucar.edu>`_ program.
     
   output
     File to receive stdout output from the model. This file will be created if 
     it does not exist, or overwritten if it does exist.

   resolution
     The TIEGCM can be run in one of two resolutions: 

       * 5   x 5   deg lat x lon, 2 grid levels per scale height (dz = 0.50)
       * 2.5 x 2.5 deg lat x lon, 4 grid levels per scale height (dz = 0.25)

     The resolution is set by the "modelres" shell variable in the TIEGCM 
     :ref:`job script <jobscript>`. See also the section on 
     :ref:`Grid Structure and Resolution <resolution>`.

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

   tgcmrun/
     Subdirectory under the :term:`modeldir`. The tgcmrun directory 
     contains Python code to make :term:`benchmark runs` for the current release.  
     The 'tgcmrun' command may be used to interactively submit selected
     benchmark runs, or tgcmrun can be executed from a shell script using
     command-line options. There are several run_xxxxx shell scripts there 
     demonstrating how to make benchmark runs.

   tgcmdata
     A directory path to start-up and other input data files required for 
     running the model. This should be on a large temporary disk.  tgcmdata is a 
     csh variable optionally specified in the :term:`job script`. If not specified,
     the job script will use the :term:`TGCMDATA` environment variable.  
     See also :ref:`job script shell variables <jobscript>`.

   env var TGCMDATA
     A linux environment variable that refers to the :term:`tgcmdata`. This environment
     variable may be used when referring to data files in the namelist read file, e.g.,
     "GPI_NCFILE = `$TGCMDATA/gpi_xxxxx.nc`". See :ref:`namelist read files <namelist>`.

   tgcmproc_f90
     Post-processor and visualizer for TIEGCM netCDF history files. Written in f90,
     and available at the `TIEGCM download site <http://www.hao.ucar.edu/modeling/tgcm/download.php>`_
     See :ref:`tgcmproc_f90 <tgcmproc_f90>`.

   tgcmproc_idl
     Post-processor and visualizer for TIEGCM netCDF history files. This processor is
     Written in IDL with a GUI, and is available at the 
     `TIEGCM download site <http://www.hao.ucar.edu/modeling/tgcm/download.php>`_
     See :ref:`tgcmproc_idl <tgcmproc_idl>`.

   utproc
     Post-processor and visualizer for TIEGCM netCDF history files. This processor 
     reads time-series history files and makes ut vs pressure and ut vs latitude
     contours.  It is written in IDL with a GUI, and is available at the TGCM download 
     site.  See :ref:`utproc <utproc>`.
 
   workdir
     User-created local working directory. This will typically contain the model root directory
     :term:`modeldir` and related namelist input files, job scripts, stdout files, etc. 
     Because the model source files are critical, this should be on backed-up disk, 
     typically under your home directory.

   Zp
     Vertical log pressure coordinate ln(p0/p) of the |modeluc|. This is the "lev" coordinate
     on the history files. See the chapter on :ref:`Altitude Coordinates the NCAR TIEGCM <altcoords>` 
     for a detailed explanation of the relationship between Zp and Altitude.

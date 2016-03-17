The Make/build process: compile and link
========================================

The |modeluc| model is formally supported on two platform systems: 64-bit Linux Desktop, 
and a Linux cluster supercomputer (e.g., NCAR's |ncarsuper| system. However, the model 
has been built and executed on several other platforms.  The source code is f90 standard 
compliant, and is mostly fixed-format fortran. 

Compilers
---------

The |model| can be built with three compilers on 64-bit Linux Desktop systems,
but only the intel compiler is used on the NCAR supercomputer |ncarsuper|.
The default is intel, since it out-performs the other compilers on both systems.
As of January, 2015, these are the versions of each compiler we are using: 

* Intel

  * On 64-bit desktop: intel ifort 12.0.0
  * On NCAR supercomputer |ncarsuper|: intel ifort 12.1.5

* PGI

  * On 64-bit desktop: PGI pgf90 9.0-4

* GNU gfortran

  * On 64-bit desktop: GNU gfortran 4.4.7 

Each compiler has a makefile in the scripts directory that specifies
compiler-specific flags, library paths, and other parameters necessary
for the build priocess.  These files are provided in the model 
:term:`scripts/` directory, and are included in the main Makefile
at build time:

* | Makefile for Intel compiler on NCAR supercomputer |ncarsuper|:
  | :download:`Make.intel_ys <../../scripts/Make.intel_ys>` (NCAR |ncarsuper|)

* | Makefile for Intel compiler on 64-bit Linux desktop system:
  | :download:`Make.intel_hao64 <../../scripts/Make.intel_hao64>`

* | Makefile for PGI compiler on 64-bit Linux desktop system:
  | :download:`Make.pgi_hao64 <../../scripts/Make.pgi_hao64>` 

* | Makefile for GNU gfortran compiler on 64-bit Linux desktop system:
  | :download:`Make.gfort_hao64 <../../scripts/Make.gfort_hao64>`

One of these files, or the user's own, is specified by the job script 
variable "make" in the :ref:`job script <jobscript>`. The specified 
file is included in the main :download:`Makefile <../../scripts/Makefile>`.
User's outside NCAR are encouraged to copy and rename one of these files, 
and customize it for your own operating system and compiler.

.. _external_libs::

Required External Libraries
---------------------------

External library dependencies are netCDF, MPI, and ESMF.
The MPI implementation is often bundled in with the compiler.
Local paths to these libraries are specified in the compiler-specific
``Make`` files described above.

The Earth System Modeling Framework (ESMF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The electro-dynamo code see :download:`source file pdynamo.F <../../src/pdynamo.F>`) 
in the |modeluc| is calculated on the :ref:`geomagnetic grid <magcoords>`. 
Since the dynamo receives inputs from the neutral atmosphere, which is on 
the :ref:`geographic grid <geocoords>`, there is a need for regridding 
capability between the two grid systems.  The same horizontal geomagnetic 
coordinates are used regardless of the 5-deg or 2.5-deg resolution of the 
geographic grid.

The `Earth System Modeling Framework <https://www.earthsystemcog.org/projects/esmf>`_ (see also 
`Modeling Infractructure for the Geoscience Community <http://www.cisl.ucar.edu/research/2005/esmf.jsp>`_ 
is used in the |modeluc| to perform the grid remapping in an parallel MPI environment,
see :download:`src/esmf.F <../../src/esmf.F>`. To build the |modeluc|, the ESMF library
must be included in the link step. If the ESMF library is not already on your system,
you will need to `download <https://www.earthsystemcog.org/projects/esmf/download/>`_  
and build it, using the same compiler you are using to build the |modeluc|.

HAO is using the following ESMF libraries built on a 64-bit Linux
desktop for each compiler/MPI implementation:

ESMF libraries at HAO for use on Linux desktop systems (these paths are provided
in the scripts/Make.xxxx files described above). The esmf makefiles esmf.mk are
included in the model's main makefile :download:`scripts/Makefile <../../scripts/Makefile>`

 * | For use with the Intel compiler:
   | /home/tgcm/esmf/intel/esmf_6_3_0rp1/lib/libO/Linux.intel.64.intelmpi.default
   | See :download:`Makefile esmf.mk for Intel build <_static/esmf_intel.mk>`

 * | For use with the PGI compiler:
   | /home/tgcm/esmf/pgi-9.04/lib/libO/Linux.pgi.64.mpich.default
   | See :download:`Makefile esmf.mk for PGI build <_static/esmf_pgi.mk>`

 * | For use with the GNU gfortran compiler:
   | /home/tgcm/esmf/gfort/esmf_6_3_0rp1/lib/libO/Linux.gfortran.64.openmpi.default
   | See :download:`Makefile esmf.mk for GNU gfortran build <_static/esmf_gfort.mk>`

 * | For the NCAR Linux cluster |ncarsuper|: esmf-6.3.0r-ncdfio-mpi-O
   | The ESMF library is loaded on |ncarsuper| with the "module load" command,
   | executed by the :download:`job script tiegcm-ys.job <../../scripts/tiegcm-ys.job>`.


netCDF
^^^^^^

  The Network Common Data Form (NetCDF) is a cross-platform, self-describing 
  metadata file format, developed by UNIDATA at UCAR. Please see 
  `NetCDF <http://www.unidata.ucar.edu/software/netcdf>`_ for more information.
  It is necessary to link the netCDF library when the model is built, since
  all data files imported to the model, and all model output history files
  are in NetCDF format. Because NetCDF is platform-independent, all history
  and data files can be used on either linux desktops or the NCAR |ncarsuper|
  system. At HAO on linux desktops, we are using NetCDF version 4.1.1. On the
  NCAR |ncarsuper| system, we are using version 4.3.2.

Build for Debugging
-------------------

  The model can be built with debug flags set in the compiler.
  To do this, simply set debug = TRUE in the :term:`job script`,
  and resubmit (see also :ref:`job scripts section <jobscript>`).

  The debug flags are set in the compiler-specific Make files 
  described above. They can be adjusted there, of course, but
  usually they include floating-point exception traps, and core 
  dumps with traceback. If debug was false in a previous run, 
  the entire code will be rebuilt with the debug flags set, however, 
  it doesn't hurt to go to the execdir and type "gmake clean" before 
  resubmitting. Keep in mind that because optimization is turned off
  when debug flags are set, performance will be destroyed, and the 
  model will run agonizingly slow. 

  Although we do not support the model with MPI turned off, 
  it can also be useful for debugging to run the model with
  only a single MPI task. To do this, set nproc=1 in the linux
  job script, or set #BSUB -n 1 in the |ncarsuper| job script.

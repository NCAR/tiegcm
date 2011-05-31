The Make/build process: compile and link
========================================

The TIEGCM model is formally supported on two platform systems: 64-bit Linux, 
and IBM/AIX. However, the model has been built and executed on several other platforms.
The source code is f90 standard compliant, and is mostly fixed-format fortran. 

Compilers used on Linux systems include 
`Intel's ifort <http://www.cita.utoronto.ca/MISC/computing_guide/docs/Intel_Fortran_v8.1/fcompindex.htm>`_ 
11.1 (with OpenMPI), and `PGI's fortran compiler pgf90 <http://www.pgroup.com/products/index.htm>`_. 
The compiler used on the NCAR IBM/AIX `bluefire <http://www2.cisl.ucar.edu/docs/bluefire-user-guide>`_
system is xlf90.

Library dependencies consist mainly of netCDF and MPI.
The MPI library is often bundled in with the compiler.
Locations of these libraries are specified in "Make.machine"
files, which set platform-specific compile flags and other
parameters for the build process. The following Make.machine 
files are provided in the :term:`scripts/` directory: 

  * :download:`Make.bluefire <_static/Make.bluefire>` (NCAR IBM/AIX machine)
  * :download:`Make.intel_hao64 <_static/Make.intel_hao64>` (ifort compiler on HAO 64-bit Linux desktops)
  * :download:`Make.pgi_hao32 <_static/Make.pgi_hao32>` (pgf90 compiler on HAO 32-bit Linux desktops)
  * :download:`Make.pgi_hao64 <_static/Make.pgi_hao32>` (pgf90 compiler on HAO 64-bit Linux desktops)

One of these files, or the user's own, is specified by the csh 
variable "make" in the :ref:`job script <jobscript>`. The specified 
file is included in the main :download:`Makefile <_static/Makefile>`.
User's outside NCAR are encouraged to copy and rename one of these files, 
and customize it for your own operating system and compiler.

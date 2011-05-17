
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

The user's working directory will typically look something like this::

                      workdir
                         |
 -----------------------------------------------
    |           |             |            |          
 *.inp       execdir      modeldir      datadir
 *.job          |
 *.out        *.o
              *.mod
              *.nc
              Make*
              exec

Where \*.inp are namelist input files, \*.job are job scripts, and \*.out
are stdout files from model runs. The :term:`execdir` is the build/execution
directory (created by the first model run), with object code (\*.o, \*.mod), 
model output history files (\*.nc), make files (Make\*), and an executable file. 
Various other files may also be in the execdir. The modeldir and datadir directories 
are described below.

The job script in your working directory contains a shell variable specifying
the modeldir, so it knows where to find the source code and supporting scripts
for the build process. The namelist input file refers to the datadir for start-up
and other data input files (e.g., SOURCE, GPI_NCFILE, IMF_NCFILE, etc). These
namelist parameters can use the environment variable TGCMDATA to specify the
datadir (see section on :ref:`namelist input files <namelist>`).

Model Directory (:term:`modeldir`)
----------------------------------

The model root directory is what you get when you :ref:`download <download>` the 
model source code tar file.  The model directory contains subdirectories with the 
model source code, supporting scripts, documentation, and test scripts::

                      modeldir
                         |
   -----------------------------------------------
      |          |             |              |
     src      scripts         doc           test
      |          |             |              | 
     *.F      linux.job    userguide        seasons
     *.h       ibm.job    description       climatology
            default.inp  release_notes      dec2006
                etc           etc             etc

The :term:`src/` directory contains all source code, i.e., all \*.F fortran files,
and a single header file with grid definitions and dimensions (defs.h). 
The TIEGCM source code is f90 compliant, and most files are written in fixed-format
Fortran.

The :term:`scripts/` directory contains various csh, perl and other scripts providing 
infrastructure support for the model. These include several Make files, sample 
job scripts, a default namelist input file, and several utility scripts.

The documentation directory :term:`doc/` contains the User's Guide (this document), and 
Model Description (pdf and html files). These docs are also available in the documentation
section of the main `tgcm website <http://www.hao.ucar.edu/modeling/tgcm/>`_ tgcm website.
The User's Guide is written in `Restructured Text <http://docutils.sourceforge.net/rst.html>`_, 
with the `Sphinx <http://sphinx.pocoo.org>`_ markup tool.  Source code (.rst and other files) 
for the User's Guide are available in the "userguide" subdirectory. The Model 
Description is written in `TeX <http://en.wikipedia.org/wiki/TeX>`_, and its source 
code is available in the "description" subdirectory.

Also in the doc directory are Release Notes for the current version, a README
file, and other relevant documents (e.g., table of available diagnostic fields 
from the model).

The :term:`test/` directory contains setup scripts for various test/benchmark runs,
(seasonal, full-year climatology, storm cases, etc.)  Each script produces job scripts
and input files to make the runs for that test.

Data Directory (:term:`datadir`)
--------------------------------

The initial data directory is what you get when you :ref:`download <download>` 
the data tar file. Subsequently, you may obtain additional needed data files
from the :term:`NCAR Community Data Portal`. Here is a partial schematic of the 
datadir (where "tiegcmx.xx" is the desired model version)::

                       datadir
                          |
   ---------------------------------------------
              |                      |
        source *.nc               tiegcmx.xx
        gswm *.nc                    |
        imf *.nc                 seasons/*.nc
          etc                   climatology/*.nc 
                                 dec2006/*.nc
                                    etc

Here, source \*.nc refers to :ref:`SOURCE <SOURCE>` start-up files, gswm \*.nc 
refers to :ref:`GSWM data files <GSWM>`, and imf \*.nc refers to 
:ref:`IMF_NCFILE <IMF_NCFILE>` data files.

"tiegcmx.xx" refers to the version of the model that was downloaded. This
subdirectory contains history file output from test/benchmark runs executed
by that version of the model (see :term:`test/` subdirectory of :term:`modeldir`).




.. _source_section:

Model Source Code
=================

The source code is in the :term:`src/` subdirectory of the model root directory 
(:term:`modeldir`), which is provided in the model :ref:`download <download>` file.


The Academic License Agreement
------------------------------

The TIEGCM :download:`Open Source Academic Research License Agreement <_static/tiegcmlicense.txt>`
specifies the terms and restrictions under which the NCAR/UCAR grants permission to use the
model, including the source code, for research, academic, and non-profit purposes. 

Source Code Flow Diagram
------------------------

A detailed flow diagram and calling tree of the source code structure is available
in single and multi-page pdf files:

.. Warning::

  Some details of these flow charts are out of date with respect to TIEGCM version |version|

* :base_url:`TIEGCM Code Structure (multi-page pdf) <code_structure/tiegcm_codestruct.pdf>`

* :base_url:`TIEGCM Code Structure (single-page pdf) <code_structure/tiegcm_code_poster.pdf>`


.. _modifying_source:

Modifying the Source Code
-------------------------

As a community user, student, research scientist or developer, you may need to modify the model
source code. It is best to do this after building and at least making a default execution 
of the model (see the :ref:`QuickStart <quickstart>` Section). To change one or more 
source files, simply go to the :term:`src/` subdirectory in the model root directory
:term:`modeldir`, and edit the files as necessary. Then return to the working directory 
:term:`workdir` and re-execute the job script. It will recompile the modified files, and 
any other source files that depend on the modified files, and re-execute the model. 
Alternatively, you can enter the execution directory :term:`execdir`, and recompile 
the code by typing "gmake" on the command line, then return to the working directory 
and re-execute the job script.


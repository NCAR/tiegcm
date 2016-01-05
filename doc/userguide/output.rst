Structure of Output History Files
=================================

.. _output:

NetCDF Output Files (:term:`netCDF`)
------------------------------------

TIEGCM history files are output in `netCDF <http://www.unidata.ucar.edu/software/netcdf/>`_, 
a self-describing platform-independent data format written and maintained by the UCAR 
`Unidata <http://www.unidata.ucar.edu>`_ program.

Each netCDF file contains one or more :term:`histories <history>`, i.e., the state of the 
model output fields at a discrete instant in :term:`model time`. Here is an example of the 
metadata content of a sample primary history file: :download:`primary.ncd <_static/primary.ncd>`. 
This example file contains six daily histories (days 355 to 360 of 2002).
This metadata is obtained via the "ncdump" command in the netCDF utilities. This example 
ncdump file contains data values for scalar and singly-dimensioned vectors only. You can
also use the tgcm_ncdump script in the :term:`scripts/` directory.

TIEGCM history files are "CF compliant", i.e., they conform to the 
`NetCDF Climate and Forecast (CF) Metadata Convention <http://cfconventions.org>`_.

Primary and Secondary History Files
-----------------------------------

"Primary" history files contain the "prognostic" fields necessary to restart the model. 
They can be specified in namelist input as the :ref:`SOURCE <SOURCE>` file for starting 
the model in an :term:`initial run` (a :term:`continuation run` does not specify a SOURCE 
file, and is continued from the START time found on the first OUTPUT file).  Typically, 
daily histories are stored on primary history files.

Fields on primary histories necessary for start-up of the TIEGCM are as follows:
TN, UN, VN, O2, O1, N4S, NO, HE, AR, OP, N2D, TI, TE, NE, O2P, OMEGA, Z, POTEN 

"Secondary" history files contain diagnostic fields and/or primary history fields.
Fields to be saved on the secondary history files are listed by the namelist input 
parameter :ref:`SECFLDS <SECFLDS>`. Diagnostics can be saved by calling subroutine addfld 
in the code (addfld.F), or by including one or more of the "standard" diagnostic fields 
available via the :ref:`diagnostics <diagnostics>` module (diags.F). Secondary histories 
are often saved at a higher temporal resolution than primary histories, typically hourly. 
Here is an ncdump of an :download:`example secondary history file <_static/secondary.ncd>`


.. _diagnostics:

Saving Diagnostic Fields
========================

  The diagnostics module (diags.F) in the TIEGCM will calculate and save
  diagnostic fields to the secondary histories. The user can add any subset
  of these fields to the :ref:`SECFLDS <SECFLDS>` parameter list in the
  namelist input file. See the :download:`diagnostics namelist example <_static/diags.inp>`.

Table of Available Diagnostics
------------------------------

  A table of available diagnostic fields for TIEGCM version |version|, is 
  :download:`here <_static/diags.table>`. To request that selected fields be
  saved to the secondary histories, add the "Shortname" of the desired fields
  to the :ref:`SECFLDS <SECFLDS>` parameter list in the namelist input file. 

  Columns of the diagnostics table are described as follows:

    * Field: Field number (not significant)
    * Shortname: Short name of the field (use these in the namelist SECFLDS list)
    * Units: Units of the field (most of these are 
      `udunits <http://www.unidata.ucar.edu/software/udunits>`_ compliant)
    * Levels: The cell position in the vertical dimension on which the field was calculated:

      * lev: Field was calculated on midpoints (geographic grid)
      * ilev: Field was calculated on interfaces (geographic grid)
      * mlev: Field was calculated on midpoints (magnetic grid)
      * imlev: Field was calculated on interfaces (magnetic grid)
      * none: Field is 2d (lat x lon), so no vertical component 

    * Caller: The function or subroutine from which the diags module is called.
    * Longname: Long name of the field.

  A complete list of diagnostics for the SECFLDS namelist parameter is included
  at the end of the diags table. 

Saving Fields/Arrays from the Source
------------------------------------

  Arbitrary 2d and 3d arrays can be saved from the model to secondary histories
  by inserting a call to subroutine addfld (addfld.F) in the source code.
  There are many examples of this in the source code, just grep on "call addfld".
  For more information about how to make calls to addfld, please see comments
  in the addfld.F source file.  

  Here are a couple of examples of addfld calls from near the end of subroutine
  qrj (qrj.F). These calls are inside a latitude loop, where the loop variable
  index is "lat". Normally, in parallel code, subdomains of the field are passed,
  e.g., lon0:lon1 and lat0:lat1::

    call addfld('QO2P'  ,' ',' ',  qo2p(lev0:lev1,lon0:lon1,lat),
   |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
    call addfld('QN2P'  ,' ',' ',  qn2p(lev0:lev1,lon0:lon1,lat),
   |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
    call addfld('QNP'   ,' ',' ',   qnp(lev0:lev1,lon0:lon1,lat),
   |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

  The calling sequence for subroutine addfld is explained in comments at the
  :download:`top of addfld.F <_static/addfld.explain>`.

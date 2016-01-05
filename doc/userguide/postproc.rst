
.. _postproc:

Post-Processing and Visualization
=================================

|modeluc| netCDF history files can be read by any application with access to
the netCDF library, including many freely available software packages developed
for mainpulating or displaying netCDF data 
(see http://www.unidata.ucar.edu/software/netcdf/software.html).
At HAO, we often use the netCDF Operators `NCO <http://nco.sourceforge.net>`_
for file manipulation (subset extracting, concatenation, hyperslabbing, metadata
editing, etc). However for visualization, we typically use one of three post-processors
developed at HAO:

.. _tgcmproc_f90:

tgcmproc_f90
------------

 * :term:`tgcmproc_f90` is a batch-style processor written in fortran 90. 
   This program reads a user namelist
   input file via stdin, and outputs multi-frame plot files (cgm and/or ps), and
   output data files (e.g., ascii, netCDF).
 * Uses the freely available `NCAR Graphics libraries <http://ngwww.ucar.edu/>`_ for 
   basic contouring, making maps at various projections, vector plots, etc.
 * Plots 2d horizontal and vertical slices, and time-dependent plots on the model grid.
 * Calculates a large number of diagnostics from fields on the histories.
 * Custom contouring (setting cmin,cmax,cint)
 * Can interpolate to constant height surfaces.
 * Can be downloaded from the TGCM website, but the f90 code must be compiled, and
   NCAR Graphics libraries must be linked. 

.. _tgcmproc_idl:

tgcmproc_idl
------------

 * :term:`tgcmproc_idl` is an 
   `IDL <http://www.ittvis.com/language/en-US/ProductsServices/IDL.aspx>`_ application
   for browsing and plotting TIEGCM output, with an easy to use Graphical User Interface (GUI).
 * 2d contouring of horizontal and vertical slices, including maps at various projections.
 * Can save images and plots to a variety of image formats.
 * Custom contouring (setting cmin,cmax,cint).
 * Can interpolate to constant height surfaces.
 * Can plot fields on the magnetic grid (tgcmproc_f90 does not do this).
 * Can make and save png movie animations.
 * Should run fine for anybody w/ IDL, but IDL is licensed, and can be expensive.

.. _utproc:

utproc
------

 * :term:`utproc` is an IDL/GUI application that makes time-series contours and images
   including ut vs zp pressure at selected grid lat x lon locations, and ut vs 
   latitude at selected zp pressure surfaces.

These applications are available at the 
`TGCM download page <http://www.hao.ucar.edu/modeling/tgcm/download.php>`_.
Tgcmproc_f90 is best for generating large numbers of plots in a "batch-style" environment,
whereas tgcmproc_idl is best for browsing history files in a GUI interface, and saving plots 
or images as desired. The utproc processor is a hybrid in the sense that a series of plots
can be setup using the GUI, and then created when requested.

At HAO, we also use the NCAR Command Language (`NCL <http://www.ncl.ucar.edu/>`_)
for plotting, analysis, and converting to/from various file formats (GRIB, HDF, etc).
NCL scripts can be used to generate customized plots and images, as well as providing
a variety of analysis and file-manipulation tools.

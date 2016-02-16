
Glossary
--------

.. glossary::

  Heelis
    The Heelis Potential Model is optionally used for high-latitude ion drift velocities
    Namelist input setting: POTENTIAL_MODEL = 'HEELIS'

  Weimer
    The Weimer Potential Model is optionally used for high-latitude ion drift velocities
    Namelist input setting: POTENTIAL_MODEL = 'WEIMER'
 
  GPI
    Geophysical indices data can optionally be used to drive the model, including 3-hourly Kp,
    and F10.7 Solar Flux (average and daily). Namelist input setting: 
    GPI_NCFILE = '$TGCMDATA/tiegcm2.0/gpi_1960001-2015090.nc'

  IMF
    Interplanetary Magnetic Field data from OMNI, used to drive the Weimer potential model.
    Example namelist input setting: IMF_NCFILE = '$TGCMDATA/imf_OMNI_2006001-2006365.nc'

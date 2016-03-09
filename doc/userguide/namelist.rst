
.. index:: namelist input, input

.. _namelist:

Run Control Parameters: The Namelist Input File
===============================================

The :term:`namelist input` file specifies user-provided run control 
parameters for a model run. These parameters include the model startup
file, start and stop times, solar inputs, lower boundary files, and 
several other flags and data files as necessary.  This document describes 
each valid namelist parameter, their valid combinations, and provides 
several example input files for running the model in different modes 
and resolutions.

Example Namelist Input Files
----------------------------

Please refer to the following examples of namelist input files:

.. note::

   Any part of a line in the namelist file following a exclamation mark '!' will
   be treated as a comment (see example files).

.. index:: default ; namelist read file

Example Namelist Input Files:

* The default input files: 

  * :download:`5-degree resolution namelist file <../../scripts/tiegcm_res5.0_default.inp>`
  * :download:`2.5-degree resolution namelist file <../../scripts/tiegcm_res2.5_default.inp>`

* A continuation of the default 5-degree run: 

  * :download:`continuation.inp <_static/namelist_files/continuation.inp>` 

* Saving diagnostics on secondary history files: 

  * :download:`diags.inp <_static/namelist_files/diags.inp>`

.. _namelist_params:

Explanation of Valid Namelist Parameters
-----------------------------------------

Following is a table of valid TIEGCM |version| namelist parameters, and their descriptions. 
Follow the parameter name links to explanations below.

=========================================== ===================== =====================================
Parameter Name                              Data Type and Default Description
=========================================== ===================== =====================================
:ref:`AMIENH,AMIESH <AMIE>`                 string: [none]        Optional AMIE data files
:ref:`AMIE_IBKG <AMIE_IBKG>`                integer: 0            Flag for how to read AMIE data
:ref:`AURORA <AURORA>`                      integer: 1            0/1 flag for auroral parameterization
:ref:`BGRDDATA_NCFILE <BGRDDATA_NCFILE>`    string: [none]        Data file for background lower boundary 
:ref:`BXIMF or BXIMF_TIME <BXIMF>`          real or real array    X-component of the IMF
:ref:`BYIMF or BYIMF_TIME <BYIMF>`          real or real array    Y-component of the IMF
:ref:`BZIMF or BZIMF_TIME <BZIMF>`          real or real array    Z-component of the IMF
:ref:`CALENDAR_ADVANCE <CALENDAR_ADVANCE>`  real: 1               0/1 switch to advance calendar time
:ref:`COLFAC <COLFAC>`                      real: 1.5             O-O+ collision factor
:ref:`CTMT_NCFILE <CTMT_NCFILE>`            string: [none]        Lower boundary data file T,U,V,Z
:ref:`CTPOTEN <CTPOTEN>`                    real:                 Cross-Tail Potential
:ref:`CTPOTEN_TIME <CTPOTEN>`               real: [none]          Time-dependent Cross-Tail Potential
:ref:`CURRENT_KQ <CURRENT_KQ>`              integer: 0            Height-integrated Current Density
:ref:`CURRENT_PG <CURRENT_PG>`              integer: 1            Current due to Plasma Pressure Gradient
:ref:`CALC_HELIUM <CALC_HELIUM>`            integer: 1            0/1 switch for calculation of Helium
:ref:`DYNAMO <DYNAMO>`                      integer: 1            0/1 switch for electro-dynamo
:ref:`EDDY_DIF <EDDY_DIF>`                  integer: 0            0/1 switch for DOY-dependent or constant eddy diffusion
:ref:`ENFORCE_OPFLOOR <ENFORCE_OPFLOOR>`    integer: 1            Gaussian shaped floor for O+
:ref:`F107 or F107_TIME <F107>`             real or real array    Daily F10.7 cm solar flux
:ref:`F107A or F107A_TIME <F107A>`          real or real array    81-day average F10.7 cm solar flux
:ref:`GPI_NCFILE <GPI_NCFILE>`              string: [none]        Geophysical Indices (Kp) data file
:ref:`GSWM data files <GSWM>`               string: [none]        GSWM Model tidal lbc data files
:ref:`HIST <HIST>`                          integer(3)            Primary history write frequency
:ref:`IMF_NCFILE <IMF_NCFILE>`              string: [none]        IMF OMNI data files
:ref:`JOULEFAC <JOULEFAC>`                  real: 1.5             Joule Heating Factor
:ref:`KP or KP_TIME <KP>`                   real or real array    Kp for calc of hpower and ctpoten
:ref:`LABEL <LABEL>`                        string:               Arbitrary string identifying the run
:ref:`MXHIST_PRIM <MXHIST_PRIM>`            integer: 10           Max histories on primary file
:ref:`MXHIST_SECH <MXHIST_SECH>`            integer: 24           Max histories on secondary file
:ref:`OPDIFFCAP <OPDIFFCAP>`                real: 0 [no cap]      Impose maximum O+ diffusion
:ref:`OUTPUT <OUTPUT>`                      string array          Primary history output file(s)
:ref:`POTENTIAL_MODEL <POTENTIAL_MODEL>`    string: [HEELIS]      High-latitude Potential Model
:ref:`POWER or POWER_TIME <POWER>`          real or real array    Hemispheric Power (GW)
:ref:`SABER_NCFILE <SABER_NCFILE>`          string: [none]        SABER data file (T,Z lower boundary condition)
:ref:`SECSTART <SECSTART>`                  integer(3)            Secondary history start time (day,hour,minute)
:ref:`SECSTOP <SECSTOP>`                    integer(3)            Secondary history stop time (day,hour,minute)
:ref:`SECHIST <SECHIST>`                    integer(3)            Secondary history write frequency (day,hour,minute)
:ref:`SECFLDS <SECFLDS>`                    string array          Fields to be stored on secondary histories
:ref:`SECOUT <SECOUT>`                      string array          Secondary history output file(s)
:ref:`SOURCE <SOURCE>`                      string: [none]        Primary SOURCE (start-up) file
:ref:`SOURCE_START <SOURCE_START>`          integer(3)            Model time to start on SOURCE file 
:ref:`START <START>`                        integer(3)            Model start time (day,hour,minute)
:ref:`START_YEAR <START_YEAR>`              integer: 2002         Starting year
:ref:`START_DAY <START_DAY>`                integer: 80           Starting day of year
:ref:`STEP <STEP>`                          integer: [none]       Model time step (seconds)
:ref:`STOP <STOP>`                          integer(3)            Model stop time (day,hour,minute)
:ref:`SWDEN or SWDEN_TIME <SWDEN>`          real or real array    Solar Wind Density
:ref:`SWVEL or SWVEL_TIME <SWVEL>`          real or real array    Solar Wind Velocity
:ref:`TIDE <TIDE>`                          real(10)              Amplitudes and phases of semi-diurnal tide (rarely used)
:ref:`TIDE2 <TIDE2>`                        real(2)               Amplitudes and phases of diurnal tide (rarely used)
:ref:`TIDI_NCFILE <TIDI_NCFILE>`            string: [none]        TIDI data file (U,V lower boundary condition)
=========================================== ===================== =====================================

.. -------------------------------------------------------------------------------------
.. index:: amienh, amiesh, namelist input ; amienh, amiesh
.. _AMIE:
.. describe:: AMIENH, AMIESH

   Data files containing output from the AMIE model, to be imported to the tiegcm.
   AMIENH contains northern hemisphere data, AMIESH contains southern hemisphere data.  
   Contact Gang Lu (ganglu@ucar.edu) for more information.

   | Data type: string
   | Default: [none]

   Example:

    * AMIENH = '$TGCMDATA/amie_apr01_10_2010_nh_ssusi.nc'
    * AMIESH = '$TGCMDATA/amie_apr01_10_2010_sh_ssusi.nc'

   See also: 
   
    * :ref:`AMIE_IBKG <AMIE_IBKG>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: amie_ibkg, namelist input ; amie_ibkg
.. _AMIE_IBKG:
.. describe:: AMIE_IBKG

   Integer flag 0, 1, or 2 for reading real, first, or 24-hour averaged AMIE data

   | Data type: scalar integer
   | Default: 0

   See also: 
   
    * :ref:`AMIE <AMIE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: aurora, namelist input ; aurora
.. _AURORA:
.. describe:: AURORA

   If AURORA > 0 then the auroral parameterization (aurora.F) is called by dynamics
   (dynamics.F), otherwise it is not called. 

   | Data type: scalar integer
   | Default: 1

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: bgrddata_ncfile, namelist input ; bgrddata_ncfile
.. _BGRDDATA_NCFILE:
.. describe:: BGRDDATA_NCFILE

   Data file providing zonal mean climatology of T, U, V using MSIS and HWM empirical models,
   or UARS data. If no input file is specified, a flat lower boundary (u=v=0, Tn=181 K, z=96.4 km) 
   is employed by default. Other zonal mean climatologies can be used by generating and 
   specifying a different file.

   Data type: string
   Default: [none]

   Example:

     * BGRDDATA_NCFILE = '$TGCMDATA/bgndlbc_hwm_msis.nc'
     * BGRDDATA_NCFILE = '$TGCMDATA/bgndlbc_saber_hrdi.nc'

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: bximf, namelist input ; bximf
.. _BXIMF:
.. describe:: BXIMF or BXIMF_TIME

   X-component of the IMF. Can be specified as either a constant (BXIMF), or series of 
   time-dependent values (BXIMF_TIME). If IMF_NCFILE is set and BXIMF is not provided, 
   then BXIMF will be taken from the IMF data file.

   Data type: real or real array

   Examples:
     * BXIMF = 0. ; constant for entire run
     * BXIMF_TIME = 80,0,0,40., 80,1,0,30., 80,5,0,20. ; time series

   See also: 
     * :ref:`BYIMF or BYIMF_TIME <BYIMF>`
     * :ref:`BZIMF or BZIMF_TIME <BZIMF>`
     * :ref:`IMF_NCFILE <IMF_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: byimf, namelist input ; byimf
.. _BYIMF:
.. describe:: BYIMF or BYIMF_TIME

   Y-component of the IMF. Can be specified as either a constant (BYIMF), or series of 
   time-dependent values (BYIMF_TIME). If IMF_NCFILE is set and BYIMF is not provided, 
   then BYIMF will be taken from the IMF data file.

   Data type: real or real array

   Examples:
     * BYIMF = 0. ; constant for entire run
     * BYIMF_TIME = 80,0,0,40., 80,1,0,30., 80,5,0,20. ; time series

   See also: 
     * :ref:`BXIMF or BYIMF_TIME <BXIMF>`
     * :ref:`BZIMF or BZIMF_TIME <BZIMF>`
     * :ref:`IMF_NCFILE <IMF_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: bzimf, namelist input ; bzimf
.. _BZIMF:
.. describe:: BZIMF or BZIMF_TIME

   Z-component of the IMF. Can be specified as either a constant (BZIMF), or series of 
   time-dependent values (BZIMF_TIME). If IMF_NCFILE is set and BZIMF is not provided, 
   then BZIMF will be taken from the IMF data file.

   Data type: real or real array

   Examples:
     * BZIMF = 0. ; constant for entire run
     * BZIMF_TIME = 80,0,0,40., 80,1,0,30., 80,5,0,20. ; time series

   See also: 
     * :ref:`BXIMF or BXIMF_TIME <BXIMF>`
     * :ref:`BYIMF or BYIMF_TIME <BYIMF>`
     * :ref:`IMF_NCFILE <IMF_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: calendar_advance, namelist input ; calendar_advance
.. _CALENDAR_ADVANCE:
.. describe:: CALENDAR_ADVANCE

   Set CALENDAR_ADVANCE=1 to advance calendar time from START_DAY, otherwise 
   calendar time is not advanced. If advancing calendar time, iday (init_module) 
   is incremented every 24 hours, and the sun's declination and longitude is recalculated 
   (see sub advance_day in advance.F and sub sunloc in magfield.F), thereby allowing 
   seasonal change to take place. The earth's orbital eccentricity "sfeps" is also 
   updated as a 6% variation in solar output over a year.

   A run with CALENDAR_ADVANCE=0 is referred to as a "steady-state" run. This is often 
   used to advance the model to a "steady-state" for a given date, prior to a seasonal 
   run with CALENDAR_ADVANCE=1. 

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: colfac, namelist input ; colfac
.. _COLFAC:
.. describe:: COLFAC

   O-O+ Collision Frequency, alias the "Burnside Factor". Default is 1.5, but there 
   have been recommendations for 1.3. COLFAC is used in lamdas.F and oplus.F.

   | Data type: real
   | Default: 1.5 

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: ctmt_ncfile, namelist input ; ctmt_ncfile
.. _CTMT_NCFILE:
.. describe:: CTMT_NCFILE

   Tidal perturbations for lower boundary of Z, T, U, V.  See this reference::

     Ref. Oberheide, J., J. M. Forbes, X. Zhang, and S. L. Bruinsma
       Climatology of upward propagating diurnal and semidiurnal tides in the thermosphere
       J. Geophys. Res., 116, A11306, doi:10.1029/2011JA016784, 2011.

   Note: This is mutually incompatible with GSWM_NCFILE

   | Data type: string
   | Default: [none]

   Examples:
     * CTMT_NCFILE = '$TGCMDATA/ctmt_tides.nc'

   See also:
     * :ref:`GSWM <GSWM>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: ctpoten, namelist input ; ctpoten
.. _CTPOTEN:
.. describe:: CTPOTEN or CTPOTEN_TIME

   Cross-tail (or cross-cap) potential. This is used in the auroral precipitation 
   parameterization. It can be provided either as a single constant (CTPOTEN), or 
   several time-dependent values (CTPOTEN_TIME). If GPI_NCFILE is set and CTPOTEN 
   is not provided, it will be calculated from 3-hourly Kp data read from GPI_NCFILE.

   The time-dependent example below specifies increasing CTPOTEN from model times 
   80,0,0 to 80,1,0, and 80,5,0. Interpolated values will be used between these 
   specified model times.

   Note that if POTENTIAL_MODEL='WEIMER' or 'WEIMER05', then the user is not allowed 
   to provide CTPOTEN because it will be calculated from the Weimer electric potential. 

   | Data type: real or real array

   Examples:
     * CTPOTEN = 60.
     * CTPOTEN_TIME = 80,0,0,60., 80,1,0,65., 80,5,0,100.

   See also:
     * :ref:`POWER or POWER_TIME <POWER>`
     * :ref:`KP or KP_TIME <KP>`
     * :ref:`GPI_NCFILE <GPI_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: current_kq, namelist input ; current_kq
.. _CURRENT_KQ:
.. describe:: CURRENT_KQ

   If CURRENT_KQ=1, then height-integrated current density of current sheet, and
   upward current density at the top of the ionosphere is calculated (default=0)
   (ignored if DYNAMO=0) (see current.F90 to save JQR, JE13D, JE23D, KQPHI, KQLAM)

   | Data type: integer
   | Default: 0

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: eddy_dif, namelist input ; eddy_dif
.. _EDDY_DIF:
.. describe:: EDDY_DIF

   If EDDY_DIF=1, then day-of-year dependent eddy diffusion will be calculated, otherwise
   eddy diffusion will be set to pressure-dependent constants. See cons.F.

   | Data type: integer
   | Default: 0

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: enforce_opfloor, namelist input ; enforce_opfloor
.. _ENFORCE_OPFLOOR:
.. describe:: ENFORCE_OPFLOOR

   A double-Gaussian shaped floor (in latitude and altitude) is applied to O+ at 
   low-to-mid latitudes in the F-region in order to keep the model stable when the 
   ionosphere gets very low in density.

   If set to 1 (the default), the floor is implemented in oplus.F as follows::

       do k=lev0,lev1-1
         opfloor = opmin*exp(-(glat(lat)/90.0)**2/0.3)*
     |             exp(-((zpmid(k)-4.25)/zpmid(nlevp1))**2/0.1)
         do i=lon0,lon1
           if (opout(k,i,lat) < opfloor) opout(k,i,lat) = opfloor
         enddo
       enddo

   | Data type: integer
   | Default: 1

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: current_pg, namelist input ; current_pg
.. _CURRENT_PG:
.. describe:: CURRENT_PG

   If CURRENT_PG=1, current due to plasma pressure gradient and gravity is calculated
   and included as a forcing term in the dynamo equation (default=1) (ignored if DYNAMO=0)

   | Data type: integer
   | Default: 1

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: calc_helium, namelist input ; calc_helium
.. _CALC_HELIUM:
.. describe:: CALC_HELIUM

   If calc_helium=1, helium is calculated as a major composition species.
   If calc_helium=0, helium is zeroed out. If calc_helium=1 and the source history
   does not have helium, then helium will be initialized globally to 0.1154E-5.

   | Data type: integer
   | Default: 1

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: dynamo, namelist input ; dynamo
.. _DYNAMO:
.. describe:: DYNAMO

   Integer switch for electro-dynamo. If DYNAMO=0, then dynamo (pdynamo.F) will not be
   called, and ion drift velocities will be zero.  If DYNAMO=1, then dynamo will be
   called, and ion drift velocities will be calculated.

   | Data type: integer
   | Default: 1

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: f107, namelist input ; f107
.. _F107:
.. describe:: F107 or F107_TIME

   Daily F10.7 cm solar flux. This can be provided either as a single constant (F107), or 
   several time-dependent values (F107_TIME). If GPI_NCFILE is set and F107 is not set, 
   then F107 will be set from the data. The below example of F107_TIME increases the f10.7
   flux from 120 to 150 in the first hour of model time, then to 200 by the fifth hour.
   Values are linearly interpolated at each time-step.

   Data type: real or real array

   Examples:
     * F107 = 120.
     * F107_TIME = 80,0,0,120., 80,1,0,150., 80,5,0,200. 

   See also: 
     * :ref:`F107A <F107A>`
     * :ref:`POTENTIAL_MODEL <POTENTIAL_MODEL>`
     * :ref:`GPI_NCFILE <GPI_NCFILE>`
     * :ref:`IMF_NCFILE <IMF_NCFILE>`
   
   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: f107a, namelist input ; f107a
.. _F107A:
.. describe:: F107A or F107A_TIME

   81-day average F10.7 cm solar flux. This can be provided either as a single constant 
   (F107A), or several time-dependent values (F107A_TIME). If GPI_NCFILE is set and F107A 
   is not set, then F107A will be set from the data. The below example of F107A_TIME
   increases the f10.7a flux from 120 to 130 in 12 hours of model time.

   Data type: real or real array

   Examples:
     * F107A = 120.
     * F107A_TIME = 80,0,0,120., 80,6,0,125., 80,12,0,130. 

   See also: 
     * :ref:`F107 <F107>`
     * :ref:`POTENTIAL_MODEL <POTENTIAL_MODEL>`
     * :ref:`GPI_NCFILE <GPI_NCFILE>`
     * :ref:`IMF_NCFILE <IMF_NCFILE>`
   
   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: gpi, namelist input ; gpi_ncfile
.. _GPI_NCFILE:
.. describe:: GPI_NCFILE

   Specifies a netCDF data file containing 3-hourly Kp and daily F10.7 data to drive 
   high-latitude convection and the auroral precipitation oval.  If GPI_NCFILE is specified, 
   and POTENTIAL_MODEL='HEELIS', then at least one of CTPOTEN,POWER,F107,F107A must **not** 
   be specified. If CTPOTEN or POWER are not specified, they are calculated from the Kp data 
   using empirical relationships (see source file gpi.F).  If F107 or F107A are not specified, 
   the data will be used.

   If GPI_NCFILE is specified when POTENTIAL_MODEL='WEIMER' and IMF_NCFILE is specified,
   then the Weimer model and aurora will be driven by the IMF data, and only F107 and F107A 
   will be read from the GPI data file (F107 is not available on IMF data files).

   If the current model time is not available on the GPI data file, the model will print 
   an error message to stdout, and stop.

   Data Source: Ascii data is obtained from NOAA/NGDC, and an equivalent netCDF data file 
   is written for import to the TGCM models (see code in hao:$TGCMROOT/mkgpi). 

   Datatype: string

   Example:
     * GPI_NCFILE = '$TGCMDATA/gpi_2000001-2009031.nc'

   See also: 
     * :ref:`CTPOTEN or CTPOTEN_TIME <CTPOTEN>`
     * :ref:`POWER or POWER_TIME <POWER>`
     * :ref:`F107 or F107_TIME <F107>`
     * :ref:`IMF_NCFILE <IMF_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: gswm, namelist input ; gswm
.. _GSWM:
.. describe:: GSWM model data files for lbc

   Paths to netCDF data files containing tidal perturbations from the Global Scale Wave Model. 
   If provided, the files will be read and the perturbations will be added to the lower 
   boundary conditions of T,U,V, and Z. If provided, then TIDE and TIDE2 must be zeroed out.

   Warning: As of version |version|, the model is not tuned to use the non-migrating 
   GSWM tidal components.  The default namelist input file specifies migrating diurnal and 
   semi-diurnal tides, but not the non-migrating components. In later releases, non-migrating 
   tides may be supported at the 2.5-deg resolution.

   GSWM files must contain data compatable with the lower boundary of the model (99 km),
   and the horizontal resolution of the model being run (either 5 or 2.5 degrees). 
   See examples below. 

   Datatype: string

   Examples:
    * GSWM files for the 5-degree TIEGCM::

        GSWM_MI_DI_NCFILE   = '$TGCMDATA/gswm_diurn_5.0d_99km.nc'
        GSWM_MI_SDI_NCFILE  = '$TGCMDATA/gswm_semi_5.0d_99km.nc'
        GSWM_NMI_DI_NCFILE  = '$TGCMDATA/gswm_nonmig_diurn_5.0d_99km.nc'
        GSWM_NMI_SDI_NCFILE = '$TGCMDATA/gswm_nonmig_semi_5.0d_99km.nc'

    * GSWM files for 2.5-degree TIEGCM::

        GSWM_MI_DI_NCFILE   = '$TGCMDATA/gswm_diurn_2.5d_99km.nc'
        GSWM_MI_SDI_NCFILE  = '$TGCMDATA/gswm_semi_2.5d_99km.nc'
        GSWM_NMI_DI_NCFILE  = '$TGCMDATA/gswm_nonmig_diurn_2.5d_99km.nc'
        GSWM_NMI_SDI_NCFILE = '$TGCMDATA/gswm_nonmig_semi_2.5d_99km.nc'

   See also:
     * :ref:`TIDE <TIDE>`
     * :ref:`TIDE2 <TIDE2>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: hist, namelist input ; hist
.. _HIST:
.. describe:: HIST

   Primary history write frequency, specified as a model time (day,hour,minute). 
   HIST time must divide evenly into STOP minus START times.

   Examples:
     * HIST = 1,0,0    ;request daily histories
     * HIST = 0,1,0    ;request hourly histories
     * HIST = 0,0,12   ;request 12-minute histories

   See also:
     * :ref:`SECHIST <SECHIST>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: hpower, namelist input ; hpower
.. _POWER:
.. describe:: POWER or POWER_TIME

   Hemispheric Power (GW). This is used in the auroral precipitation parameterization. 
   It can be provided either as a single constant (POWER), or several time-dependent 
   values (POWER_TIME). If GPI_NCFILE is set and POWER is not provided, it will be 
   calculated from 3-hourly Kp data read from GPI_NCFILE.

   The time-dependent example below specifies increasing POWER from model times 
   80,0,0 to 80,1,0, and 80,5,0. Interpolated values will be used between these specified 
   model times.

   Data type: real or real array

   Examples: 
     * POWER = 16.
     * POWER_TIME = 80,0,0,16., 80,1,0,20., 80,5,0,70.

   See also:
     * :ref:`CTPOTEN or CTPOTEN_TIME <CTPOTEN>`
     * :ref:`KP or KP_TIME <KP>`
     * :ref:`GPI_NCFILE <GPI_NCFILE>`

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: saber_ncfile, namelist input ; saber_ncfile
.. _SABER_NCFILE:
.. describe:: SABER_NCFILE

   SABER data file for lower boundary conditions of T and Z 
   (neutral temperature and geopotential height).

   | Data type: string
   | Default: [none]

   See also:
     * :ref:`TIDI_NCFILE <TIDI_NCFILE>`

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: tidi_ncfile, namelist input ; tidi_ncfile
.. _TIDI_NCFILE:
.. describe:: TIDI_NCFILE

   TIDI data file for lower boundary conditions of U and V 
   (zonal and meridional neutral wind).

   | Data type: string
   | Default: [none]

   See also:
     * :ref:`SABER_NCFILE <SABER_NCFILE>`

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: kp, namelist input ; kp
.. _KP:
.. describe:: KP or KP_TIME

   Geomagnetic Activity index. If KP is specified and POWER and/or CTPOTEN are commented,
   then the given KP will be used with empirical formulas to calculate POWER and/or CTPOTEN,
   which are used in the Auroral parameterization.

   KP can be provided as a scalar constant (KP), or as a series of time-dependent values
   (KP_TIME), as in the below examples. KP cannot be set if GPI_NCFILE data file is specified.

   Empirical formula used to calculate POWER from KP (see function hp_from_kp in util.F)::

      if (kp <=7.) hp_from_kp = 16.82*exp(0.32*kp)-4.86
      if (kp > 7.) hp_from_kp = 153.13 + (kp-7.)/(9.-7.)*(300.-153.13)

   Empirical formula used to calculate CTPOTEN from KP (see function ctpoten_from_kp in util.F)::

      ctpoten_from_kp = 15.+15.*kp + 0.8*kp**2

   Examples:
     * KP = 4.0
     * KP_TIME = 80,0,0,4., 80,6,0,4.5, 80,12,0,5.0

   See also:
     * :ref:`CTPOTEN <CTPOTEN>`
     * :ref:`POWER <POWER>`
     * :ref:`GPI_NCFILE <GPI_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: imf, namelist input ; imf_ncfile
.. _IMF_NCFILE:
.. describe:: IMF_NCFILE

   Specifies a netCDF data file containing hourly IMF parameters BX,BY,BZ,SWVEL, and SWDEN.
   This can be set only when POTENTIAL_MODEL='WEIMER'. The data will be used to drive
   the Weimer 2005 potential model. When set, the user must **not** provide at least one
   of the above IMF parameters.  Data will be used for IMF parameters not provided by the
   user. Values (scalar or time-dependent) that are provided by the user will take precedence 
   over the data file. 

   If the current model time is not available on the IMF data file, the model will print 
   an error message to stdout and stop.

   Notes on creation of IMF OMNI data files:

   * IMF data is derived from 1-minute OMNI satellite data available on CDAweb 
     `CDAweb <http://cdaweb.gsfc.nasa.gov/istp_public/>`_.  Our derivation is a multi-step process:   
   * Data gaps in the raw 1-minute OMNI data are linearly interpolated.  If a gap happens to occur 
     at the beginning or end of the time interval, it is set to the next known good data point.
   * Gap-filled data is used to compute a 15 minute trailing average lagged by 5 minutes.
   * Time averaged data is sampled at 5 minutes
   * A data quality flag is calculated for every 5-minute sample point.  The data quality flag is 
     a boolean value set to "1" for all sample points derived from valid (not gap-filled) data.  
     The data quality flag is set to "0" for any sample point that is derived from gap-filled data 
     anywhere in the 15 minute trailing average lagged by 5 minutes.
   * The data quality flag is stored in the NetCDF-formatted IMF input file.  For any variable 
     (ie. "swvel" - solar wind velocity), there exists a mask (ie. "velMask").  Find a complete 
     list of IMF variables with command "ncdump -h [imf-file.nc]".
   * Note:  You should verify the IMF data quality before doing storm simulations.  Known periods
     of invalid IMF data include approximately days 301 to 304 of 2003 (during the "Halloween Storm").

   Example:
     * IMF_NCFILE = '$TGCMDATA/imf_OMNI_2002001-2002365.nc'

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: joulefac, namelist input ; joulefac
.. _JOULEFAC:
.. describe:: JOULEFAC

   Joule heating factor. This factor is multiplied by the joule heating calculation
   (see subroutine qjoule_tn in qjoule.F).

   | Data type: real
   | Default: 1.5

   :ref:`Back to top <namelist_params>`
.. -------------------------------------------------------------------------------------
.. index:: label, namelist input ; label
.. _LABEL:
.. describe:: LABEL

   LABEL may be any string up to 80 characters long, used to identify a run. 
   The LABEL is written to output history files as a global file attribute. 
   This parameter is purely a user convenience, and does not effect the model 
   run in any way. 

   | Data type: string
   | Default: 'tiegcm res=5'

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: mxhist_prim, namelist input ; mxhist_prim
.. _MXHIST_PRIM:
.. describe:: MXHIST_PRIM

   Maximum number of histories to be written to primary OUTPUT files. When this many 
   histories have been written to the current OUTPUT file, the next OUTPUT file is created 
   and it receives subsequent histories. This parameter can be adjusted to control the size 
   of primary OUTPUT files.

   | Data type: integer
   | Default: 10

   Examples: 
     * MXHIST_PRIM = 15 ; allow maximum of 15 histories per primary output file 

   See also: 
     * :ref:`OUTPUT <OUTPUT>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: mxhist_sech, namelist input ; mxhist_sech
.. _MXHIST_SECH:
.. describe:: MXHIST_SECH

   Maximum number of histories to be written to secondary output files (SECOUT). 
   When this many histories have been written to the current SECOUT file, the next SECOUT 
   file is created and it receives subsequent histories. This parameter can be adjusted 
   to control the size of secondary OUTPUT files.

   | Data type: integer
   | Default: 24

   Examples:
    * MXHIST_SECH = 24 ; allow 1 day of hourly histories per file
    * MXHIST_SECH = 48 ; allow 2 days of hourly histories per file

   See also: 
     * :ref:`SECOUT <SECOUT>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: opdiffcap, namelist input ; opdiffcap
.. _OPDIFFCAP:
.. describe:: OPDIFFCAP

   Optional cap on ambipolar diffusion coefficient for O+. This can improve model 
   stability in the topside F-region, but it is only recommended as a last resort 
   since it will change model results. This is a new namelist parameter for |tgcm_version|.
   The default is 0., i.e., no cap. If this is non-zero (provided by the user), then
   it is implemented in subroutine rrk of src/oplus.F.

   | Data type: integer
   | Default: 0

   Examples::
   
   Tests have been made with these values:

     * OPDIFFCAP = 1.5e8
     * OPDIFFCAP = 3.0e8
     * OPDIFFCAP = 6.0e8
     * OPDIFFCAP = 8.0e8

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: output, namelist input ; output
.. _OUTPUT:
.. describe:: OUTPUT

   List of primary history output files. Each file may be an absolute path, or relative 
   to the execution directory. If an initial run (SOURCE is specified), then pre-existing 
   OUTPUT files will be overwritten. If a continuation run (SOURCE is *not* specified), 
   then the first OUTPUT file should contain the source history at START time. In this case, 
   subsequent output histories will be appended to the first OUTPUT file until it is full. 
   As each OUTPUT file is filled (see MXHIST_PRIM), the next OUTPUT file is created and 
   histories are written until it is full, and so on.

   OUTPUT files are usually specified with increasing integers imbedded in the names. See 
   examples below. As a convenience, large sequences of files may be specified in a "short-form", 
   see example 3 below specifying 20 files. By convention, primary history output files 
   may use the letter "p" to indicate primary file series (see all 3 examples below, and 
   contrast with SECOUT). 

   Examples::

     OUTPUT = 'p_myoutput_001.nc'
     OUTPUT = 'myoutput.p001.nc','myoutput.p002.nc','myoutput.p003.nc'
     OUTPUT = 'myoutput_p001.nc','to','myoutput_p020.nc','by','1' 

   See also:
     * :ref:`SECOUT <SECOUT>`
     * :ref:`SOURCE <SOURCE>`
     * :ref:`MXHIST_PRIM <MXHIST_PRIM>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: potential_model, namelist input ; potential_model
.. _POTENTIAL_MODEL:
.. describe:: POTENTIAL_MODEL

   The high-latitude potential model used to calculate electric potential above a specified 
   latitude. This string can have one of two values:

   | POTENTIAL_MODEL = 'HEELIS'
   | POTENTIAL_MODEL = 'WEIMER'

   'HEELIS' is the Rod Heelis model (heelis.F). 'WEIMER' is the Dan Weimer 2005 model 
   (wei05sc.F).

   .. note::

      The Weimer model of high-latitude potential is the intellectual property of Daniel
      Weimer and may not be extracted, distributed, or used for any purpose other
      than as implemented in the TIE-GCM.  For further information concerning this
      model, please contact Dan Weimer (dweimer@vt.edu).

   For a brief discussion of the use of the Weimer 2005 model in TIEGCM, please
   see :ref:`Notes on Weimer05 in TIEGCM <tiegcm_weimer05>`.

   | Data type: string
   | Default: 'HEELIS'

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: secflds, namelist input ; secflds
.. _SECFLDS:
.. describe:: SECFLDS

   List of fields to be saved to secondary histories. These may be either fields that are 
   also saved on primary histories (so-called "prognostic" fields), fields that have been 
   requested via addfld calls in the source code, or fields available via the 
   :ref:`diagnostics module <diagnostics>` (see example below).

   Note the final size of secondary output files is affected by the number of fields specified
   as well as the number of histories on the file.  The file size can be controlled by setting 
   the number of histories allowed on a secondary file :ref:`MXHIST_SECH <MXHIST_SECH>`.

   Data type: one or more character strings

   Examples::

     !
     ! Example list of fields to be written to secondary histories:
     !
       SECFLDS = 'TN' 'UN' 'VN' 'O2' 'O1' 'N2' 'NO' 'N4S' 'HE' 'NE' 'TE' 'TI'
                 'TEC' 'O2P' 'OP' 'OMEGA' 'POTEN' 'UI_ExB' 'VI_ExB' 'WI_ExB' 
                 'DEN' 'QJOULE' 'Z' 'ZG'
     !
     ! This example lists all diagnostic fields available via the diags module
     ! (it is not necessary to call addfld in the code to obtain these fields)
     !
       SECFLDS = 'CO2_COOL','NO_COOL','DEN','HEATING','QJOULE','QJOULE_INTEG',
           'SIGMA_PED','SIGMA_HAL','TEC','UI_ExB','VI_ExB','WI_ExB',
           'LAMDA_PED','LAMDA_HAL','HMF2','NMF2','SCHT','MU_M','O_N2','WN',
           'BX','BY','BZ','BMAG','EX','EY','EZ','ED1','ED2','PHIM2D','N2',
           'CUSP','DRIZZLE','ALFA','NFLUX','EFLUX'

   See also:
     :ref:`MXHIST_SECH <MXHIST_SECH>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: secstart, namelist input ; secstart
.. _SECSTART:
.. describe:: SECSTART

   Secondary history start time, specified as a model time (day,hour,minute). 

   | Data type: 3 integers (day,hour,minute)
   | Valid range: 0-365 for day, 0-23 for hour, 0-59 for minute.

   SECSTART time must follow these rules:
     * It must be a multiple of timestep STEP and less than SECSTOP time.
     * It must be greater than or equal to START time, and less than or equal to STOP time.
     * In the case of an initial run (SOURCE history provided), SECSTART must not be equal
       to START time. This is to avoid zero valued secondary history fields.

   Examples:
     * SECSTART = 80,1,0  ; Start saving secondary histories at model time 80,1,0

   See also:
     * :ref:`SECSTOP <SECSTOP>`
     * :ref:`SECHIST <SECHIST>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: secstop, namelist input ; secstop
.. _SECSTOP:
.. describe:: SECSTOP

   Secondary history stop time, specified as a model time (day,hour,minute). 

   | Data type: 3 integers (day,hour,minute)
   | Valid range: 0-365 for day, 0-23 for hour, 0-59 for minute.

   SECSTOP time must follow these rules:
     * It must be a multiple of timestep STEP and greater than SECSTART time.
     * It must be greater than or equal to START time, and less than or equal to STOP time.

   Examples:
     * SECSTOP = 81,0,0  ; Start saving secondary histories at model time 80,1,0

   See also:
     * :ref:`SECSTART <SECSTART>`
     * :ref:`SECHIST <SECHIST>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: sechist, namelist input ; sechist
.. _SECHIST:
.. describe:: SECHIST

   Secondary history write frequency, specified as a model time (day,hour,minute). 
   SECHIST time must divide evenly into SECSTOP minus SECSTART times.

   | Data type: 3 integers (day,hour,minute)
   | Valid range: 0-365 for day, 0-23 for hour, 0-59 for minute.

   Examples:
     * SECHIST = 0,1,0    ;request hourly histories
     * SECHIST = 0,0,12   ;request 12-minute histories

   See also:
     * :ref:`HIST <HIST>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: secout, namelist input ; secout
.. _SECOUT:
.. describe:: SECOUT

   List of secondary history output files. Secondary histories store diagnostic fields,
   usually at a higher temporal resolution than primary files. Each file may be an 
   absolute path, or relative to the execution directory. Beware that SECOUT will 
   overwrite any pre-existing files with the same names.  As each SECOUT file is filled 
   (see MXHIST_SECH), the next SECOUT file is created and histories are written until 
   it is full, and so on.

   SECOUT files are usually specified with increasing integers imbedded in the names. See 
   examples below. As a convenience, large sequences of files may be specified in a "short-form", 
   see example 3 below specifying 20 files. By convention, secondary history output files 
   may use the letter "s" to indicate secondary file series (see all 3 examples below).

   Examples::

     SECOUT = 's_myoutput_001.nc'
     SECOUT = 'myoutput.s001.nc','myoutput.s002.nc','myoutput.s003.nc'
     SECOUT = 'myoutput_s001.nc','to','myoutput_s020.nc','by','1' 

   See also:
     * :ref:`OUTPUT <OUTPUT>`
     * :ref:`SOURCE <SOURCE>`
     * :ref:`MXHIST_SECH <MXHIST_SECH>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: source, namelist input ; source
.. _SOURCE:
.. describe:: SOURCE

   SOURCE is the start-up history file for an initial run. SOURCE may be a full path or 
   relative to the execution directory. It must be a TIEGCM history with the same grid 
   resolution as the model being run. It does not need to be from the same model version 
   as that being run.

   If SOURCE is specified, then SOURCE_START, the model time of the history to read on the 
   SOURCE file, must also be specified. The code will search the SOURCE file for the 
   SOURCE_START history. If SOURCE is *not* specified, then the run is a continuation run, 
   and the source history is provided in the first OUTPUT file at START time.

   The SOURCE file must be on the local disk. The model will not look for the SOURCE
   history on any archive file system.

   Examples:
     * SOURCE = '$TGCMDATA/tiegcm_res2.5_decsol_smax_prim.nc'

   See also:
     * :ref:`SOURCE_START <SOURCE_START>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: source_start, namelist input ; source_start
.. _SOURCE_START:
.. describe:: SOURCE_START

   This is the model time of the history to read from the SOURCE file. Model time is 
   specified as a 3-integer triplet: day,hour,minute. If SOURCE is specified, then 
   SOURCE_START must also be specified. If the SOURCE_START history is not found on the 
   SOURCE file, the model will stop and print an appropriate error message to stdout.

   | Data type: 3 integers
   | Valid range: 0-365 for day, 0-23 for hour, 0-59 for minute.

   Example: 
     * SOURCE_START = 80,0,0 

   See also:
     * :ref:`SOURCE <SOURCE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: start, namelist input ; start
.. _START:
.. describe:: START

   Model time for start of the run. Model time is a 3-integer triplet: day,hour, minute. 
   If CALENDAR_ADVANCE=0, then START day can be any number between 0 and 365. 
   If CALENDAR_ADVANCE=1, then START day must be the same as START_DAY. If an initial run, 
   START time does not have to be the same as SOURCE_START.

   Data type: 3 integers
   Valid range: 0-365 for day, 0-23 for hour, 0-59 for minute.

   Examples: 
     * START = 80,0,0

   See also:
     * :ref:`SOURCE_START <SOURCE_START>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: start_day, namelist input ; start_day
.. _START_DAY:
.. describe:: START_DAY

   Calendar starting day.

   | Data type: integer
   | Default: 80
   | Valid range: 1 to 365

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: start_year, namelist input ; start_year
.. _START_YEAR:
.. describe:: START_YEAR

   Starting year for the run.

   | Data type: integer
   | Default: 2002

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: step, namelist input ; step
.. _STEP:
.. describe:: STEP

   Model time-step in seconds. Default value is 60 seconds for 5-degree resolution, 
   30 seconds for 2.5-degree resolution. During periods of quiet solar activity, 
   the model can often be run at twice these times. During periods of intense solar 
   activity (e.g., F10.7 > 200, or high magnitude BZ southward), the model may become 
   numerically unstable. In this case, reducing the timestep to as low as 10 seconds 
   may be necessary for the model get through the rough period.

   | Data type: integer
   | Default for 5.0-degree resolution: STEP=60
   | Default for 2.5-degree resolution: STEP=30

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: stop, namelist input ; stop
.. _STOP:
.. describe:: STOP

   Model stop time for the run. Model time is specified as a 3-integer triplet: day,hour,minute.

   | Data type: 3 integers
   | Valid range: 0-365 for day, 0-23 for hour, 0-59 for minute.

   Example: 
     * STOP = 81,0,0

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: swden, namelist input ; swden
.. _SWDEN:
.. describe:: SWDEN or SWDEN_TIME

   Solar Wind Density (#/cm-3). Can be specified as either a constant (SWDEN), or series of 
   time-dependent values (SWDEN_TIME). If IMF_NCFILE is set and SWDEN is not provided, 
   then it SWDEN will be taken from the IMF data file.

   Data type: real or real array

   Examples:
    * SWDEN = 4.0 ; constant for entire run
    * SWDEN_TIME = 80,0,0,2., 80,1,0,3., 80,2,0,4. ; time series

   See also:
     * :ref:`IMF_NCFILE <IMF_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: swvel, namelist input ; swvel
.. _SWVEL:
.. describe:: SWVEL or SWVEL_TIME

   Solar Wind Velocity (Km/s). Can be specified as either a constant (SWVEL), or series of 
   time-dependent values (SWVEL_TIME). If IMF_NCFILE is set and SWVEL is not provided, 
   then it SWVEL will be taken from the IMF data file.

   Data type: real or real array

   Examples:
    * SWVEL = 400. ; constant for entire run
    * SWVEL_TIME = 80,0,0,100., 80,1,0,200., 80,2,0,300. ; time series

   See also:
     * :ref:`IMF_NCFILE <IMF_NCFILE>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: tide, namelist input ; tide
.. _TIDE:
.. describe:: TIDE

   Hough mode amplitudes and phases of the semi-diurnal tide. If GSWM tidal perturbations    
   are specified, TIDE should be set to 0.

     .. note::
        TIDE and TIDE2 should be specified only for experiments where amplitude and phases
        of the tides must be used.  Normally, GSWM tides are specified instead of TIDE,TIDE2.

   Data type: 10 reals

   Example::

       TIDE= 1.9300E+04, 1.5000E+04, 2.3100E+04, 0.7700E+04, 0.1660E+04,
             -2.600E+00,  0.000E+00, -3.300E+00, 4.2000E+00, 5.0000E+00

   See also: 
     * :ref:`GSWM_MI_SDI_NCFILE <GSWM>`
     * :ref:`GSWM_NM_SDI_NCFILE <GSWM>`

   :ref:`Back to top <namelist_params>`

.. -------------------------------------------------------------------------------------
.. index:: tide2, namelist input ; tide2
.. _TIDE2:
.. describe:: TIDE2

   Hough mode amplitudes and phases of the diurnal tide. If GSWM tidal perturbations are 
   specified, TIDE2 should be set to 0.

     .. note::
        TIDE and TIDE2 should be specified only for experiments where amplitude and phases
        of the tides must be used.  Normally, GSWM tides are specified instead of TIDE,TIDE2.

   Data type: 2 floats

   Example::

     TIDE2 = 4.1E+4, -3.7  
 
   See also: 
     * :ref:`GSWM_MI_DI_NCFILE <GSWM>`
     * :ref:`GSWM_NM_DI_NCFILE <GSWM>`

   :ref:`Back to top <namelist_params>`

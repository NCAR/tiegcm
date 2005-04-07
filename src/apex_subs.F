!
! Original file ~bozo/pgms/apex/apxntrp.f copied 2/25/00.
!
      SUBROUTINE APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
     +                  WK,LWK, IST)
C          This contains ENTRYs that determine quantities related to Apex
C          coordinates.  They are designed to speed up coordinate determination
C          by substituting interpolation from previously calculated tables
C          for direct calculation.  They also permit calculating quantities
C	   involving gradients of Apex coordinates, such as base vectors in
C	   the Modified-Apex and Quasi-Dipole systems.  Options allow table 
C	   creation, storage, and read-back in addition to interpolation.
C
C          Tables (arrays) must be prepared prior to coordinate determination.
C          Tables may be created for one date and held only in memory by calling
C
C            APXMKA - make magnetic arrays.
C
C          Or, they may be created and written, then read back later by
C
C            APXWRA - make and write arrays
C            APXRDA - read stored arrays.
C
C          Tables may be created for multiple times (e.g., the DGRF dates)
C          and written to a file using APXWRA.  However, tables are held in
C          memory for only one time due to the potential to exceed memory
C          capacity.  Thus, when reading back multiples times, APXRDA
C          interpolates/extrapolates to a single time.
C
C          Alternatively, If APXWRA is called to make and write tables for one
C          time only, then it is not necessary to call APXRDA because tables
C          have already been initialized for that time.
C
C          Creation of global lookup tables can be time consuming, thus, it
C          is expected that this is done rarely; it might be considered an
C          installation step.  Once the tables are established, programs may
C          reference them by calling APXRDA.  In this case, the grid coordinates
C          (originally specified during the table creation) may be ascertained
C          by calling
C
C            APXGGC - get grid coordinates.
C
C          After the tables have been loaded in memory for a given date, it
C          is possible to interpolate in space to determinine various magnetic
C          parameters using the following calls:
C
C            APXALL  - get Apex coordinates (radius, lat, lon).
C            APXMALL - get Modified Apex coordinates, Quasi-Dipole
C                      coordinates and associated base vectors.
C            APXQ2G  - convert quasi-dipole to geodetic coordinates.
C
C          Details for each ENTRY introduced above are provided below.
C          Following the ENTRY summaries are comments describing EXTERNALS,
C          INSTALLATION SPECIFICS for different computers, and a HISTORY.
C
C          REFERENCE:
C          Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
C          Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
C
C          ALGORITHM:
C          When arrays are created, APXMKA calls subroutine MAKEXYZV, which
C          in turn calls subroutine APEX at each grid point to get the apex
C          radius A, the apex longitude PHIA, and the magnetic potential
C          VMP.  The cosine (CLP) and sine (SLP) of the quasi-dipole
C          latitude are computed from A.  From these can be computed
C          preliminary quantites defined as
C
C            x = cos(quasi-dipole latitude)*cos(apex longitude)
C            y = cos(quasi-dipole latitude)*sin(apex longitude)
C            z = sin(quasi-dipole latitude)
C            v = (VMP/VP)*((RE+ALT)/RE)**2
C
C          where VP is the magnitude of the magnetic potential of the
C          geomagnetic dipole at a radius of RE; ALT is altitude; and RE is
C          the mean Earth radius.  Note that all of these quantities vary
C          smoothly across the apex poles, unlike apex latitude or
C          quasi-dipole latitude, so that they can be linearly interpolated
C          near these poles.  Corresponding values of x,y,z,v for a dipole
C          field on a spherical Earth are computed analytically and
C          subtracted from the above quantities, and the results are put
C          into the 3D arrays X,Y,Z,V.  When APXALL or APXMALL is called,
C          trilinear interpolations (in latitude, longitude, and inverse
C          altitude) are carried out between the grid points.  Gradients are
C          calculated for APXMALL in order to determine the base vectors.
C          Analytic formulas appropriate for a dipole field on a spherical
C          Earth are used to determine the previously removed dipole
C          components of x,y,z,v, and their gradients and these are added
C          back to the interpolated values obtained for the non-dipole
C          components.  Finally, the apex-based coordinates and their base
C          vectors are calculated from x,y,z,v and their gradients.
C
C------------------------------------------------------------------------------
C
C       ENTRY APXMKA and APXWRA:
C          These create gridded arrays of quantities that can later be
C          linearly interpolated to any desired geographic point within
C          their range by APXALL, APXMALL, and APXQ2G.  The spatial extent
C          and resolution of the arrays in latitude, longitude, and altitude
C          can be tailored for the application at hand.  In one extreme,
C          very high interpolation accuracy can be achieved by creating a
C          2x2x2 array with very small spatial increments surrounding the
C          point of interest. (However, since vector quantities are obtained
C          by taking differences of the quantities at adjacent points,
C          computer roundoff errors will limit the benefits of making the
C          increments too small.) In another extreme, a global array from
C          the ground to an altitude of several Earth radii can be created.
C          (In this case, the spatial resolution is limited by the need to
C          maintain reasonable array sizes.) For most purposes, we have
C          found it adequate to use a global grid with dimensions 91,151,6
C          (latitude,longitude,altitude) and altitude from ground to 1274 km
C          as created by GGRID with NVERT = 30 (see file test.f).
C
C          WARNING: execution time to create these look-up tables can be
C          substantial.  It is dependent on the grid dimensions and the
C          number of Epochs to be computed.  For eight epochs and dimensions
C          (91,151,6), it took 76 minutes on a Sun SPARCstation IPX and the
C          resulting file is 15.8 Mbytes.
C
C          Use APXMKA to create the interpolation tables for a single time
C          without writing them.  Otherwise, use APXWRA create the tables
C          and save them in a file.
C
C              CALL APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
C             +            WK,LWK, IST)
C
C              CALL APXWRA (MSGUN, FILNAM,IUN, EPOCH,NEPOCH,
C             +            GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            EPOCH  = Time formatted as UT year and fraction; e.g., 1990.0
C                     is 0 UT 1 Jan 1990.  For APXMKA EPOCH is single valued;
C                     for APXWRA EPOCH contains NEPOCH values.
C            GPLAT,GPLON,GPALT = Grid point latitudes, longitudes, and
C                     altitudes.  NLAT values in GPLAT must be in the range
C                     -90. to +90. degrees North.  NLON values in GPLON
C                     must be in the range -180. to 180. degrees East.  NALT
C                     values in GPALT are km MSL.  Each array must be in
C                     ascending numerical order but they may have arbitrary
C                     spacing (as desired for model grid resolution).
C                     Subroutine GGRID can be used to set up these arrays.
C            NLAT,NLON,NALT = Number of latitudes, longitudes and altitudes
C                     are the single dimensions of the arrays GPLAT,GPLON,
C                     and GPALT.  Subroutine GGRID can be used to select
C                     appropriate values.  In general, increasing the number
C                     of grid points, increases the model resolution.  While
C                     array values at the grid points are exact, values at
C                     locations between grid points are estimated by linear
C                     interpolation.  One can increase the grid resolution to
C                     the point of degraded accuracy very close to the poles
C                     of quantities involving east-west gradients, viz.,
C                     B(1), G, H, W, Bhat(1), D1(1), D2(1), D3, E1, E2, E3,
C                     F1(2), and F2(2).  This is evident with dimensions
C                     301x501x15 for a global grid, 0-1000 km.
C            WK    =  Work array is used internally, i.e., there is no need
C                     to access the contents.  WK should not be altered
C                     between initialization (by APXMKA, APXWRA or APXRDA)
C                     and use (by APXGGC, APXALL, APXMALL, or APXQ2G).
C            LWK   =  Dimension of WK >=  NLAT*NLON*NALT*5 + NLAT+NLON+NALT
C
C          Additional INPUTS for APXWRA only:
C            FILNAM = file name where arrays are stored.
C            IUN    = Fortran unit number to be associated with FILNAM.
C            NEPOCH = Number of times in EPOCH.
C
C          RETURNS:
C            IST = Return status: = 0  okay
C                                 > 0  failed
C
C          Declarations for APXMKA formal arguments:
C            DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT), WK(LWK)
C
C          Additional declarations for APXWRA formal arguments:
C            DIMENSION EPOCH(NEPOCH)
C            CHARACTER*(*) FILNAM
C
C------------------------------------------------------------------------------
C
C       ENTRY APXRDA:
C          Read back tables (previously created by calling APXWRA) in
C          preparation for magnetic coordinate determination (using APXALL,
C          APXMALL, APXQ2G).  If multiple dates were stored, interpolate
C          in time.
C
C              CALL APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            FILNAM = file name where arrays are stored.
C            IUN    = Fortran unit number to be associated with FILNAM.
C            DATE   = Time formatted as UT year and fraction; e.g., 1990.0
C                     is 0 UT 1 Jan 1990.
C            WK,LWK = Same as APXMKA
C
C          RETURNS:
C            IST = Return status: = 0  okay
C                                 = -1 non-fatal date problem
C                                 > 0  failed
C
C          Declarations for APXRDA formal arguments:
C            CHARACTER*(*) FILNAM
C            DIMENSION WK(LWK)
C
C------------------------------------------------------------------------------
C
C       ENTRY APXGGC:
C          Obtain grid coordinates from tables read back using APXRDA.
C
C              CALL APXGGC (MSGUN, WK,LWK, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
C             +            IST)
C
C          INPUTS:
C            MSGUN  = Fortran unit number to write diagnostic messages.
C            WK,LWK = Same as APXMKA
C
C          RETURNS:
C            GPLAT,GPLON,GPALT = Grid point latitudes, longitudes, and
C                                altitudes.  Units are degrees and kilometers.
C            NLAT,NLON,NALT    = Number of latitudes, longitudes and altitudes.
C            IST               = Return status, where IST=0 implies okay.
C
C          Declarations for APXGGC formal arguments:
C            DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT), WK(LWK)
C
C
C------------------------------------------------------------------------------
C
C
C       ENTRY APXALL:
C          Determine Apex coordinates by interpolation from precalculated
C          arrays.  First call APXMKA, APXWRA, or APXRDA to load look-up
C          tables in memory.
C
C              CALL APXALL (GLAT,GLON,ALT, WK, A,ALAT,ALON, IST)
C
C          INPUTS:
C            GLAT = Geographic (geodetic) latitude, degrees, must be within
C                   the grid domain (GPLAT(1) <= GLAT <= GPLAT(NLAT)).
C            GLON = Geographic (geodetic) longitude, degrees, must be within
C                   one revolution of the grid domain:
C                     GPLON(1) <= GLON-360.,GLON, or GLON+360. <= GPLON(NLON))
C            ALT  = Altitude, km
C            WK   = same as entry APXMKA
C          RETURNS:
C            A    = Apex radius, normalized by Req
C            ALAT = Apex latitude, degrees
C            ALON = Apex longitude, degrees
C            IST  = Return status:  okay (0); or failure (1).
C
C          Dimensions of non-scalar arguments to APXALL:
C            WK(LWK)
C
C
C------------------------------------------------------------------------------
C
C       ENTRY APXMALL:
C          Compute Modified Apex coordinates, quasi-dipole coordinates,
C          base vectors and other parameters by interpolation from
C          precalculated arrays.  First call APXMKA, APXWRA, or APXRDA
C          to load look-up tables in memory.
C
C              CALL APXMALL (GLAT,GLON,ALT,HR, WK,
C             +             B,BHAT,BMAG,SI,
C             +             ALON,
C             +             XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,
C             +             XLATQD,F,F1,F2 , IST)
C
C          Reference:  Richmond, A. D., Ionospheric Electrodynamics Using
C          Magnetic Apex Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
C
C          INPUTS:
C            GLAT = Geographic (geodetic) latitude, degrees, must be within
C                   the grid domain (GPLAT(1) <= GLAT <= GPLAT(NLAT)).
C            GLON = Geographic (geodetic) longitude, degrees, must be within
C                   one revolution of the grid domain:
C                     GPLON(1) <= GLON-360.,GLON, or GLON+360. <= GPLON(NLON))
C            ALT  = Altitude, km
C            HR   = Reference altitude, km (used only for Modified Apex coords)
C            WK   = same as entry APXMKA
C          RETURNS:
C            B      = magnetic field components (east, north, up), in nT
C            BHAT   = components (east, north, up) of unit vector along
C                     geomagnetic field direction
C            BMAG   = magnitude of magnetic field, in nT
C            SI     = sin(I)
C            ALON   = apex longitude = modified apex longitude = quasi-dipole
C                     longitude, degrees
C            XLATM  = Modified Apex latitude, degrees
C            VMP    = magnetic potential, in T.m
C            W      = W of reference above, in km**2 /nT  (i.e., 10**15 m**2 /T)
C            D      = D of reference above
C            BE3    = B_e3 of reference above (= Bmag/D), in nT
C            SIM    = sin(I_m) described in reference above
C            D1,D2,D3,E1,E2,E3 = components (east, north, up) of base vectors
C                     described in reference above
C            XLATQD = quasi-dipole latitude, degrees
C            F      = F described in reference above for quasi-dipole
C                     coordinates
C            F1,F2  = components (east, north) of base vectors described in
C                     reference above for quasi-dipole coordinates
C            IST    = Return status:  okay (0); or failure (1).
C
C          Dimensions of non-scalar arguments to APXMALL:
C            GPLAT(NLAT),GPLON(NLON),GPALT(NALT),WK(LWK),
C            B(3),BHAT(3),D1(3),D2(3),D3(3), E1(3),E2(3),E3(3), F1(2),F2(2)
C
C
C------------------------------------------------------------------------------
C
C       ENTRY APXQ2G:
C          Convert from quasi-dipole to geodetic coordinates, APXQ2G
C          (input magnetic, output geodetic) is the functional inverse
C          of APXALL or APXMALL (input geodetic, output magnetic).  First
C          call APXMKA, APXWRA, or APXRDA to load look-up tables in memory.
C
C              CALL APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
C
C          INPUTS:
C            QDLAT = quasi-dipole latitude in degrees
C            QDLON = quasi-dipole longitude in degrees
C            ALT   = altitude in km
C            WK    = same as entry APXMKA
C          RETURNS:
C            GDLAT = geodetic latitude in degrees
C            GDLON = geodetic longitude in degrees
C            IST   = Return status: =  0 okay
C                                   = -1 non-fatal, results are not as
C                                        close as PRECISE
C                                   >  0 failure
C
C          Dimensions of non-scalar arguments to APXQ2G:
C            WK(LWK)
C
C
C------------------------------------------------------------------------------
C
C          EXTERNALS:
C            apex.f   - APEX, etc. source code
C            magfld.f - COFRM, DYPOL, FELDG are the DGRF/IGRF.
C            magloctm.f, subsol.f and cossza.f - are not externals, but they
C                       are related, computing magnetic local time, the sub-
C                       solar point, and the cosine of the solar zenith angle.
C            test.f  -  also not an external, rather it is an example driver
C                       program demonstrating various call sequences.  It
C                       includes subroutine GGRID.
C
C          INSTALLATION SPECIFICS:
C            The only (known) machine dependency is parameter IRLF described
C            below.
C            Compiler default is usually to initialize memory to 0.  If this
C            is not true, an error trap involving KGMA may not work.
C
C          HISTORY:
C            Originally coded 940824 by A. D. Richmond, NCAR.  Components of
C            routines in the package came from previous work by Vincent
C            Wickwar (COFRM, FELDG in magfld.f).
C
C            Modified Sep 95 (Roy Barnes):  Changes were made with the
C            objective to allow the user to completely control definition of
C            the interpolation grid.  While doing this the order of the
C            ENTRYs in the first subroutine and arguments lists were
C            changed:  APXMKA, APXWRA, APXRDA (formerly GETARR) and the other
C            ENTRYs now include a work array (WK) which holds arrays X,Y,Z,V,
C            GPLAT,GPLON and GPALT.  Subroutine SETLIM was removed and the
C            grid creation algorithm based on NVERT originally integral to
C            GETARR and SETLIM has been extracted, but still available as
C            an example (see GGRID in file test.f).  Subroutine TSTDIM has a
C            different role, so it is now CKGP (check grid points).
C            MAKEXYZV was also changed to accomodate explicit grid point
C            arrays instead of computed values from an index, minimum and
C            delta.  Only one format is written now, so that is is possible
C            to concatinate files for different epochs.  This required
C            changing delta time logic from fixed 5 yr epochs to computed
C            differences which may vary.  Also, the ties to DGRF/IGRF dates
C            have been removed.
C
C            Modified Sep 96 (Roy Barnes):  Corrected bug in APXQ2G longitude
C            iteration.  The latitude iteration is now constrained to not
C            step beyond the (partial global) interpolation grid.  NITER was
C            increased from 9 to 14 and code to determine cos(lambda') (CLP)
C            was revised by Art Richmond to reduce truncation errors near
C            the geographic poles.
C
C            Modified Sep 97:  Corrected comments, definition of COLAT.
C
C            Modified Dec 98:  Change GLON input to try +/-360 values
C            before rejecting when out of defined grid (GDLON) range;
C            affects INTRP
C
C            Modified Feb-Mar 99:  (1) Corrected a typo (bad variable name)
C            in diagnostic print in INTRP: GLO -> GLON.  This error was
C            probably introduced Dec 98, so no-one had access to the
C            bad copy and, therefore, no announcement is needed.  (2) Also
C            modified APXMALL and APXQ2G:  When very close to a geographic
C            pole, gradients are recalculated after stepping away; in this
C            situation, the latitude input to INTRP was changed from GLAT
C            to GLATX.  This is affects gradients when within 0.1 degrees
C            of a pole (defined as the larger of GLATLIM, 89.9 deg, or the
C            second largest grid point).  (3) Changed MAKEXYZV to make X,Y,Z,
C            V constant for all longitudes (at each alt) when poleward of
C            POLA; POLA was added to /APXCON/.  (4)  Replaced definition of
C            DYLON in APXQ2G.  (5) Reduced NITER to 10 from 14.  Changes 3-5
C            fix a problem where APXQ2G calculations were failing to satisify
C            PRECISE within NITER iterations at the pole.  (6) Replace
C            XYZ2APX with revised trigonometry to avoid a truncation problem
C            at the magnetic equator.  Most changes were devised by Art
C            Richmond and installed by Roy B.
C
C            Questions about this version should be directed to Roy Barnes
C            NCAR (email: bozo@ucar.edu phone: 303/497-1230).
C
C------------------------------------------------------------------------------

      PARAMETER (XMISS=-32767. , GLATLIM=89.9 , PRECISE=7.6E-11,
     +          DATDMX=1. , DATIMX=2.5 , IRLF=8 )
C          XMISS   = value used to fill returned variables when requested
C                    point is outside array boundaries
C          GLATLIM = Extreme polar latitude allowed before changing east-west
C                    gradient calculation to avoid possible underflow at
C                    poles.  GLATLIM is superseded when the second to last
C                    grid point value is closer to a pole.
C          PRECISE = (angular distance)**2 (radians**2) of precision of
C                    transform from quasi-dipole to geodetic coordinates.
C                    7.6E-11 corresponds to an accuracy of 0.0005 degree.
C          IRLF    = Record length factor required to convert the computed
C                    length of direct access records from words to those
C                    units appropriate to this computer.  IRLF is 1 when
C                    RECL units are words, otherwise it is a function of
C                    the computer word length; e.g.,
C                      IRLF = 1 on a DEC  (32-bit words, RECL units are words)
C                      IRLF = 4 on a PC   (32-bit words, RECL units are bytes)
C                      IRLF = 4 on a Sun  (32-bit words, RECL units are bytes)
C                      IRLF = 8 on a Cray (64-bit words, RECL units are bytes)
C          DATDMX  = maximum time difference (years) allowed between the
C                    requested date and the single epoch in arrays.
C          DATIMX  = maximum time difference (years) allowed between the
C                    requested date and the closest epoch in the stored
C                    arrays (apropos multiple stored dates).

C          Formal argument declarations
      DIMENSION GPLAT(*),GPLON(*),GPALT(*), EPOCH(*), WK(*),
     +          GRADX(3), GRADY(3), GRADZ(3), GRADV(3),
     +          GRCLM(3), CLMGRP(3), RGRLP(3),
     +          B(3),BHAT(3), D1(3),D2(3),D3(3), E1(3),E2(3),E3(3),
     +          F1(2),F2(2)
      CHARACTER*(*) FILNAM

C          Local declarations
      CHARACTER CALNM*7, CMP*10, EDG*5
      SAVE KGMA, GLALMN,GLALMX, NLA,NLO,NAL, LBX,LBY,LBZ,LBV,LLA,LLO,LAL
      DATA KGMA /0/

C          Common APXDIPL is assigned in MAKEXYZV but computed in DYPOL
C            COLAT = geocentric colatitude (degrees) of north geomagnetic pole
C            ELON  = geocentric east longitude (degrees) of north geomagnetic
C                    pole
C            VP    = magnetic potential at 1 RE, south geomagnetic pole
C            CTP   = cos(colat*dtor)
C            STP   = sin(colat*dtor)
      COMMON /APXDIPL/  COLAT,ELON,VP,CTP,STP

C          Common APXCON is assigned here
C            RTOD, DTOR = radians to degrees (180/pi) and inverse
C            RE, REQ    = 6371.2, 6378.160 (Mean and equatorial Earth radius)
C            MSGU       = MSGUN to be passed to subroutines
C            POLA       = Pole angle (deg); when the geographic latitude is
C                         poleward of POLA, X,Y,Z,V are forced to be constant.
C                         for all longitudes at each altitude
!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

      CALNM = 'APXMKA'
      LCN   = 6
      KGMA  = 1
      MSGU  = MSGUN
      NEPOK = 1
      IF (NLAT .LT. 2 .OR. NLON .LT. 2 .OR. NALT .LT. 2) GO TO 9100
      NLA = NLAT
      NLO = NLON
      NAL = NALT
      GO TO 40

      ENTRY APXGGC (MSGUN,WK,LWK, GPLAT,GPLON,GPALT,NLAR,NLOR,NALR,IST)
C          Sep 95 R. Barnes
      CALNM = 'APXGGC'
      LCN   = 6
      MSGU  = MSGUN
      IF (KGMA .LT. 1) GO TO 9300
      J = LLA
      DO 10 I=1,NLA
      GPLAT(I) = WK(J)
   10 J = J + 1
      DO 20 I=1,NLO
      GPLON(I) = WK(J)
   20 J = J + 1
      DO 30 I=1,NAL
      GPALT(I) = WK(J)
   30 J = J + 1
      NLAR = NLA
      NLOR = NLO
      NALR = NAL
      IST = 0
      RETURN

      ENTRY APXWRA (MSGUN, FILNAM,IUN, EPOCH,NEPOCH,
     +             GPLAT,GPLON,GPALT,NLAT,NLON,NALT, WK,LWK, IST)
C          Sep 95 R. Barnes
      CALNM = 'APXWRA'
      LCN   = 6
      KGMA  = 2
      MSGU  = MSGUN
      NEPOK = NEPOCH
      IF (NLAT .LT. 2 .OR. NLON .LT. 2 .OR. NALT .LT. 2) GO TO 9100
      NLA = NLAT
      NLO = NLON
      NAL = NALT
      GO TO 40

      ENTRY APXRDA (MSGUN, FILNAM,IUN, DATE, WK,LWK, IST)
C          Sep 95 R. Barnes
      CALNM = 'APXRDA'
      LCN   = 6
      KGMA  = 3

C          Open the read-back file with a temporary record length, get
C          the grid dimensions from the first values, then close it, (so
C          it can be reopened later with the proper LDR):
      LDR = 7*IRLF
      OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='old',
     +                                                       IOSTAT=IST)
      MSGU = MSGUN
      IF (IST .NE. 0) GO TO 9110
      READ  (IUN,REC=1,IOSTAT=IST) YEAR,COLAT,ELON,VP,NLA,NLO,NAL
      IF (IST .NE. 0) GO TO 9120
      CLOSE (IUN)

   40 RE   = 6371.2
      REQ  = 6378.160
      RTOD = 45./ATAN(1.)
      DTOR = 1./RTOD
      POLA = 90. - SQRT (PRECISE) * RTOD
      LFN = 0
      IF (KGMA .EQ. 1) GO TO 51
      DO 50 I=1,LEN(FILNAM)
      IF (FILNAM(I:I) .EQ. ' ') GO TO 51
   50 LFN = LFN + 1
   51 CONTINUE

C          Save grid dimensions, establish direct access rec length, and
C          determine indices into the work array.  WK contains arrays
C          X,Y,Z,V,temp,GPLAT,GPLON,GPALT where X thru tmp are dimensioned
C          (NLAT,NLON,NALT); tmp is scratch space used during read back.
      NGP = NLA*NLO*NAL
      NGM1= NGP - 1
      LDR = NGP * IRLF
      LBX = 1
      LBY = LBX + NGP
      LBZ = LBY + NGP
      LBV = LBZ + NGP
      LBT = LBV + NGP
      LLA = LBT + NGP
      LLO = LLA + NLA
      LAL = LLO + NLO
      LEG = LAL + NAL-1
      IF (LWK .LT. LEG) GO TO 9130

      IF (KGMA .EQ. 3) GO TO 200

C          Make and optionally write interpolation arrays for NEPOK times
      IF (KGMA .EQ. 2) THEN
	OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='new',
     +       IOSTAT=IST)
	IF (IST .NE. 0) GO TO 9115
      ENDIF

      CALL CKGP (CALNM(:LCN),MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
      IF (IST .NE. 0) RETURN
      I = LLA - 1
      DO 60 J=1,NLAT
      I = I + 1
   60 WK(I) = GPLAT(J)
      DO 70 J=1,NLON
      I = I + 1
   70 WK(I) = GPLON(J)
      DO 80 J=1,NALT
      I = I + 1
   80 WK(I) = GPALT(J)

      IF (NEPOK .LT. 1) GO TO 9140
      J = 1
      DO 100 I=1,NEPOK
      CALL MAKEXYZV (EPOCH(I),NLAT,NLON,NALT,GPLAT,GPLON,GPALT,
     +               WK(LBX),WK(LBY),WK(LBZ),WK(LBV))
      IF (KGMA .EQ. 1) GO TO 100
      WRITE (IUN,REC=J) EPOCH(I),COLAT,ELON,VP,NLAT,NLON,NALT
      WRITE (IUN,REC=J+1) (WK(K),K=LLA,LEG)
      WRITE (IUN,REC=J+2) (WK(K),K=LBX,LBX+NGM1)
      WRITE (IUN,REC=J+3) (WK(K),K=LBY,LBY+NGM1)
      WRITE (IUN,REC=J+4) (WK(K),K=LBZ,LBZ+NGM1)
      WRITE (IUN,REC=J+5) (WK(K),K=LBV,LBV+NGM1)
  100 J = J + 6
      IF (KGMA .EQ. 2) CLOSE (IUN)
      IST = 0
      GO TO 300

C          Read back interpolation arrays.  When arrays for multiple times
C          are available, interpolate using the pair bounding the desired
C          time (DATE).  Make an initial pass only to identify closest
C          available times and issue any diagnostics, then the second pass
C          to read the stored arrays (GPLAT,GPLON,GPALT,X,Y,Z,V) and do
C          the linear interpolation/extrapolation.
  200 OPEN (IUN,FILE=FILNAM,ACCESS='direct',RECL=LDR,STATUS='old',
     +      IOSTAT=IST)
      IF (IST .NE. 0) GO TO 9110

      READ (IUN,REC=1,IOSTAT=IST) TB
      IF (IST .NE. 0) GO TO 9120
      I2 = 1
      TL = TB
      IL = 1
      I = 1
  210 I = I + 6
      READ (IUN,REC=I,IOSTAT=JST) T
C         JST .NE. 0 is assumed to mean read beyond last record
      IF (JST .NE. 0) GO TO 220
      TO = TL
      TL = T
      IL = I
      IF (DATE .GT. TL) GO TO 210

  220 I1 = IL - 6

      IST = 0
      IF (TL .EQ. TB) THEN
	DATD = ABS (DATE-TB)
	IF (DATD .GT. DATDMX) THEN
	  WRITE (MSGU,9150) CALNM(1:LCN),DATE,DATD,TB
	  IF (TB .EQ. 0.) WRITE (MSGU,9155) FILNAM(1:LFN)
	  IST = -1
	ENDIF
	I1 = 1
	I2 = 0
      ELSE IF (DATE .LT. TB) THEN
	WRITE (MSGU,9160) CALNM(1:LCN),DATE,TB,FILNAM(1:LFN)
	IST = -1
      ELSE IF (DATE .GT. TL) THEN
	WRITE (MSGU,9170) CALNM(1:LCN),DATE,TL,FILNAM(1:LFN)
	IST = -1
      ELSE
	DATD = min (DATE-TO,TL-DATE)
	IF (DATD .GT. DATIMX) THEN
	  WRITE (MSGU,9180) CALNM(1:LCN),DATE,TB,TL,FILNAM(1:LFN),DATD
	  IST = -1
	ENDIF
      ENDIF

      READ (IUN,REC=I1) YEAR1,COLAT,ELON,VP,NLAI,NLOI,NALI
      TI = YEAR1
      IF (NLAI.NE.NLA .OR. NLOI.NE.NLO .OR. NALI.NE.NAL) GO TO 9190
      READ (IUN,REC=I1+1) (WK(I),I=LLA,LEG)
      READ (IUN,REC=I1+2) (WK(I),I=LBX,LBX+NGM1)
      READ (IUN,REC=I1+3) (WK(I),I=LBY,LBY+NGM1)
      READ (IUN,REC=I1+4) (WK(I),I=LBZ,LBZ+NGM1)
      READ (IUN,REC=I1+5) (WK(I),I=LBV,LBV+NGM1)
      IF (I2 .EQ. 1) THEN
	READ (IUN,REC=I1+6) YEAR2,COLA2,ELON2,VP2,NLAI,NLOI,NALI
	TI = YEAR2
	IF (NLAI.NE.NLA .OR. NLOI.NE.NLO .OR. NALI.NE.NAL) GO TO 9190
	LE = LBT + NLA+NLO+NAL - 1
	READ (IUN,REC=I1+7) (WK(I),I=LBT,LE)
	J = LLA
	DO 230 I=LBT,LE
	IF (WK(J) .NE. WK(I)) GO TO 9200
  230   J = J + 1
	FRAC = (DATE-YEAR1) / (YEAR2-YEAR1)
	OMF  = 1. - FRAC
	LE = LBT + NGM1
	READ (IUN,REC=I1+8) (WK(I),I=LBT,LE)
	J = LBX
	DO 240 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  240   J = J + 1
	READ (IUN,REC=I1+9) (WK(I),I=LBT,LE)
	DO 250 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  250   J = J + 1
	READ (IUN,REC=I1+10) (WK(I),I=LBT,LE)
	DO 260 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  260   J = J + 1
	READ (IUN,REC=I1+11) (WK(I),I=LBT,LE)
	DO 270 I=LBT,LE
	WK(J) =  OMF*WK(J) + FRAC*WK(I)
  270   J = J + 1
	COLAT = OMF*COLAT + FRAC*COLA2
	ELON  = OMF*ELON  + FRAC*ELON2
	VP    = OMF*VP    + FRAC*VP2
      ENDIF

      CTP  = COS (COLAT*DTOR)
      STP  = SIN (COLAT*DTOR)
      YEAR = DATE
      CLOSE (IUN)

C          Establish for this grid polar latitude limits beyond which east-west
C          gradients are computed differently to avoid potential underflow
  300 GLALMX = max ( GLATLIM,WK(LLA+NLA-2))
      GLALMN = min (-GLATLIM,WK(LLA+1))

      RETURN

C*******************************************************************************

      ENTRY APXMALL (GLAT,GLO ,ALT,HR, WK,                              !Inputs
     +              B,BHAT,BMAG,SI,                                     !Mag Fld
     +              ALON,                                               !Apx Lon
     +              XLATM,VMP,W,D,BE3,SIM,D1,D2,D3,E1,E2,E3,            !Mod Apx
     +              XLATQD,F,F1,F2 , IST)                               !Qsi-Dpl
C          940822 A. D. Richmond, Sep 95 R. Barnes

C          Test to see if WK has been initialized
      CALNM = 'APXMALL'
      LCN   = 7
      IF (KGMA .LT. 1) GO TO 9300

C          Alias geodetic longitude input in case INTRP adjusts it by +/-360
      GLON = GLO

      CALL INTRP (GLAT,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +           NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL),
     +           FX,FY,FZ,FV,
     +           DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,
     +           DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM(1:LCN),IST)

      IF (IST .NE. 0) THEN
	CALL SETMISS (XMISS, XLATM,ALON,VMP,B,BMAG,BE3,SIM,SI,F,D,W,
     +               BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
	RETURN
      ENDIF

      CALL ADPL (GLAT,GLON,CTH,STH,FX,FY,FZ,FV,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
      CALL GRADXYZV (ALT,CTH,STH,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN,
     +          DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)

      IF (GLAT .GT. GLALMX .OR. GLAT .LT. GLALMN) THEN
C          If the point is very close to either the North or South
C          geographic pole, recompute the east-west gradients after
C          stepping a small distance from the pole.
	GLATX = GLALMX
	IF (GLAT .LT. 0.) GLATX = GLALMN
C990225 CALL INTRP (GLAT ,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
	CALL INTRP (GLATX,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +             NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL),
     +             FXDUM,FYDUM,FZDUM,FVDUM,
     +             DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,
     +             DFVDLN,DMXDH,DMYDH,DMZDH,DMVDH, CALNM(1:LCN),IST)
	CALL ADPL (GLATX,GLON,CTH,STH,FXDUM,FYDUM,FZDUM,FVDUM, DMXDTH,
     +            DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
	CALL GRAPXYZV (ALT,CTH,STH, DFXDLN,
     +                DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
      ENDIF

      CALL GRADLPV (HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV,
     +             XLATM,ALON,VMP,GRCLM,CLMGRP,XLATQD,RGRLP,B,CLM,R3_2)
      CALL BASVEC (HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2,
     +            BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
      BE3 = BMAG/D

      IST = 0
      RETURN


C*******************************************************************************

      ENTRY APXALL (GLAT,GLO ,ALT, WK, A,ALAT,ALON, IST)
C          940802 A. D. Richmond, Sep 95 R. Barnes

C          Test to see if WK has been initialized
      CALNM = 'APXALL'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

C          Alias geodetic longitude input in case INTRPSC adjusts it by +/-360
      GLON = GLO

      CALL INTRPSC (GLAT,GLON,ALT, WK(LBX),WK(LBY),WK(LBZ),
     +              NLA,NLO,NAL,   WK(LLA),WK(LLO),WK(LAL),
     +              FX,FY,FZ, CALNM(1:LCN), IST)
      IF (IST .NE. 0) GO TO 600

      CALL ADPLSC (GLAT,GLON,FX,FY,FZ)

      CALL XYZ2APX (ALT,FX,FY,FZ,A,ALAT,ALON,IST)
      IF (IST .EQ. 0) GO TO 601

  600 A    = XMISS
      ALAT = XMISS
      ALON = XMISS
  601 CONTINUE

      RETURN

C*******************************************************************************

      ENTRY APXQ2G (QDLAT,QDLON,ALT, WK, GDLAT,GDLON, IST)
C          940819 A. D. Richmond, Sep 95 R. Barnes, Sep 96 mod A. D. Richmond
C          Input guessed geodetic coordinates (YLAT,YLON) to INTRP and
C          compare the returned magnetic coordinates to those desired.
C          If the guess is not sufficiently close (PRECISE), make another
C          guess by moving in the direction of the gradient of a quantity
C          (DIST2) that approximates the squared angular distance between
C          the returned and desired magnetic coordinates.

C          Test to see if WK has been initialized
      CALNM = 'APXQ2G'
      LCN   = 6
      IF (KGMA .LT. 1) GO TO 9300

C          Determine quasi-cartesian coordinates on a unit sphere of the
C          desired magnetic lat,lon in quasi-dipole coordinates.
      X0 = COS (QDLAT*DTOR) * COS (QDLON*DTOR)
      Y0 = COS (QDLAT*DTOR) * SIN (QDLON*DTOR)
      Z0 = SIN (QDLAT*DTOR)

C          Initial guess:  use centered dipole, convert to geocentric coords
      CALL GM2GC (QDLAT,QDLON,YLAT,YLON)

C          Iterate until (angular distance)**2 (units: radians) is within
C          PRECISE of location (QDLAT,QDLON) on a unit sphere.
!
! 4/00: 10 iters not enough for tgcm14:
!     NITER = 10
      NITER = 20
      DO 1000 ITER=1,NITER
      CALL INTRP (YLAT,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +           NLA,NLO,NAL,    WK(LLA),WK(LLO),WK(LAL),
     +           FX,FY,FZ,FV,
     +           DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,
     +           DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM(1:LCN),IST)
      IF (IST .NE. 0) GO TO 9400
      CALL ADPL (YLAT,YLON,CTH,STH,FX,FY,FZ,FV,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)

      DISTLON = COS(YLAT*DTOR)
      IF (YLAT .GT. GLALMX .OR. YLAT .LT. GLALMN) THEN
	GLATX = GLALMX
	IF (YLAT.LT.0.) GLATX = GLALMN
	DISTLON = COS (GLATX*DTOR)
C990225 CALL INTRP (YLAT ,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
	CALL INTRP (GLATX,YLON,ALT, WK(LBX),WK(LBY),WK(LBZ),WK(LBV),
     +             NLA,NLO,NAL,     WK(LLA),WK(LLO),WK(LAL),
     +             FXDUM,FYDUM,FZDUM,FVDUM,
     +             DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,
     +             DFVDLN,DMXDH,DMYDH,DMZDH,DMVDH, CALNM(1:LCN),IST)
	CALL ADPL (GLATX,YLON,CTH,STH,FXDUM,FYDUM,FZDUM,FVDUM,
     +          DMXDTH,DMYDTH,DMZDTH,DMVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
      ENDIF

C          At this point, FX,FY,FZ are approximate quasi-cartesian
C          coordinates on a unit sphere for the quasi-dipole coordinates
C          corresponding to the geodetic coordinates YLAT, YLON.
C          Normalize the vector length of (FX,FY,FZ) to unity using XNORM
C          so that the resultant vector can be directly compared with the
C          target vector (X0,Y0,Z0).
      XNORM = SQRT(FX*FX + FY*FY + FZ*FZ)
      XDIF = FX/XNORM - X0
      YDIF = FY/XNORM - Y0
      ZDIF = FZ/XNORM - Z0
C          DIST2 = square of distance between normalized (FX,FY,FZ) and
C          X0,Y0,Z0.
      DIST2 = XDIF*XDIF + YDIF*YDIF + ZDIF*ZDIF

      IF (DIST2 .LE. PRECISE) then
!       if (iter > 10)
!    |    write(6,"('Note apxparm: took > 10 iterations: iter=',i3)")
!    |      iter
        GO TO 1001
      endif
C          HGRD2* = one-half of east or north gradient of DIST2 on unit sphere.
      HGRD2E =  (XDIF*DFXDLN + YDIF*DFYDLN + ZDIF*DFZDLN)/DISTLON
      HGRD2N = -(XDIF*DFXDTH + YDIF*DFYDTH + ZDIF*DFZDTH)
      HGRD2  = SQRT(HGRD2E*HGRD2E + HGRD2N*HGRD2N)
C          ANGDIST = magnitude of angular distance to be moved for new guess
C          of YLAT, YLON.
      ANGDIST = DIST2/HGRD2

C          Following spherical trigonometry moves YLAT, YLON to new location,
C          in direction of grad(DIST2), by amount ANGDIST.
      CAL = -HGRD2N/HGRD2
      SAL = -HGRD2E/HGRD2
      COSLM = COS(YLAT*DTOR)
      SLM = SIN(YLAT*DTOR)
      CAD = COS(ANGDIST)
      SAD = SIN(ANGDIST)
      SLP = SLM*CAD + COSLM*SAD*CAL

C          Old code (below) introduced truncation errors near geographic poles
C     SLP = min (SLP,SGLAN) ; sglan = sin(wk(lla+nla-1))*dtor)
C     SLP = max (SLP,SGLAS) ; sglas = sin(wk(lla))      *dtor)
C     CLP = SQRT(1. - SLP*SLP)
C     YLAT = ASIN(SLP)*RTOD

      CLM2 = COSLM*COSLM
      SLM2 = SLM*SLM
      SAD2 = SAD*SAD
      CAL2 = CAL*CAL
      CLP2 = CLM2 + SLM2*SAD2 - 2.*SLM*CAD*COSLM*SAD*CAL -CLM2*SAD2*CAL2
      CLP = SQRT (max(0.,CLP2))
      YLAT = ATAN2(SLP,CLP)*RTOD

C          Restrict latitude iterations to stay within the interpolation grid
C          limits, but let INTRP find any longitude exceedence.  This is only
C          an issue when the interpolation grid does not cover the entire
C          magnetic pole region.
      YLAT = min(YLAT,WK(LLA+NLA-1))
      YLAT = max(YLAT,WK(LLA))

      DYLON = ATAN2 (SAD*SAL,CAD*COSLM-SAD*SLM*CAL)*RTOD
C          Old formula (replaced Mar 99) had truncation problem near poles
C     DYLON = ATAN2 (CLP*SAD*SAL,CAD-SLM*SLP)*RTOD

      YLON  = YLON + DYLON
      IF (YLON .GT.  WK(LLO+NLO-1)) YLON = YLON - 360.
      IF (YLON .LT.  WK(LLO)      ) YLON = YLON + 360.
 1000 CONTINUE

      WRITE (MSGU,'(''APXQ2G: Warning'',I3,'' iterations only reduced th
     +e angular difference to'',/,8X,F8.5,'' degrees ('',F8.5,'' degrees
     + is the test criterion)'')') NITER, SQRT(DIST2  )*RTOD ,
     +                                    SQRT(PRECISE)*RTOD
      EDG = ' '
      IF (YLAT .EQ. WK(LLA+NLA-1)) EDG = 'north'
      IF (YLAT .EQ. WK(LLA))       EDG = 'south'
      IF (EDG .NE. ' ') WRITE (MSGU,'(''        Coordinates are on the '
     +',A,'' edge of the interpolation grid and'',/,''        latitude i
     +s constrained to stay within grid limits when iterating.'')') EDG

      IST = -1
      GO TO 1010

 1001 IST = 0

 1010 GDLAT = YLAT
      GDLON = YLON

      RETURN

C*******************************************************************************

C          Error Trap diagnostics
 9100 WRITE (MSGU,'(A,'':  NLAT,NLON or NALT < 2 '',3I8)')
     +             CALNM(1:LCN),  NLAT,NLON,NALT
      IST = 1
      RETURN
 9110 WRITE (MSGU,'(A,'': Trouble opening old file "'',A,''"'')')
     + CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9115 WRITE (MSGU,'(A,'': Trouble opening new file "'',A,''"'')')
     + CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9120 WRITE (MSGU,'(A,'': Trouble reading first record of '',A)')
     + CALNM(1:LCN), FILNAM(1:LFN)
      RETURN
 9130 WRITE (MSGU,'(A,'': LWK is too small; LWK must be at least'',I5,''
     + but LWK ='',I5)')  CALNM(1:LCN), LEG, LWK
      IST = 1
      RETURN
 9140 WRITE (MSGU,'(A,'':  NEPOCH must be positive; NEPOCH ='',I8)')
     +       CALNM(1:LCN), NEPOK
      IST = 1
      RETURN

 9150 FORMAT (A,': DATE (',F7.2,') differs by',F5.2,' years from the sto
     +red EPOCH (',F7.2,')')
 9155 FORMAT ('        A stored date = 0. implies "',A,'" is incorrectly
     + formatted')
 9160 FORMAT (A,': DATE (',F7.2,') is extrapolated before first EPOCH ('
     +,F7.2,') in "',A,'"')
 9170 FORMAT (A,': DATE (',F7.2,') is extrapolated after last EPOCH (',F
     +7.2,') in "',A,'"')
 9180 FORMAT (A,': DATE (',F7.2,') minimum difference from the nearest s
     +tored ',/,'        EPOCHs (',F7.2,', ',F7.2,') in "',A,'"',/,'
     +   is',F6.2,' years')

 9190 WRITE (MSGU,'(A,'': Dimensions (lat,lon,alt) read from "'',A,''" f
     +or EPOCH '',F7.2,/,''        ('',I5,'','',I5,'','',I5,'') do not m
     +atch ('',I5,'','',I5,'','',I5,'') from the'',/,''        first EPO
     +CH read ('',F7.2,'')'')') CALNM(1:LCN), FILNAM(1:LFN), TI ,
     + NLAI,NLOI,NALI,  NLA,NLO,NAL, TB
      IST = 1
      RETURN
 9200 CMP = 'latitudes'
      LC  = 9
      I1 = LLA - 1
      IT = LBT - 1
      N  = NLA
      IF (I .LT. LLO) GO TO 9201
      CMP = 'longitudes'
      LC   = 10
      I1 = I1 + NLA
      IT = IT + NLA
      N  = NLO
      IF (I .LT. LAL) GO TO 9201
      CMP = 'altitudes'
      LC  = 9
      I1 = I1 + NLO
      IT = IT + NLO
      N  = NAL
 9201 WRITE (MSGU,'(A,'': Grid '',A,'' read from "'',A,''" for EPOCH'',F
     +8.2,'' do not match the'',/,''        first EPOCH ('',F7.2,'')'',/
     +,''        First    Current'',/,(4X,2F10.3))')  CALNM(1:LCN),
     +  CMP(1:LC), FILNAM(1:LFN), TI,TB, (WK(I1+I),WK(IT+I),I=1,N)
      IST = 1
      RETURN

 9300 WRITE(MSGU,'(A,'': Must first load lookup tables by calling APXMKA
     +, APXWRA or APXRDA'')') CALNM(1:LCN)
      IST = 1
      RETURN
 9400 WRITE (MSGU,'(''APXQ2G: INTRP failed (maybe coordinates are not wi
     +thin interpolation grid)'')')
      IST = 1
      RETURN

      END subroutine apxmka

C*******************************************************************************
      SUBROUTINE INTRP (GLAT,GLON,ALT, X,Y,Z,V, NLAT,NLON,NALT,
     +                 GPLAT,GPLON,GPALT, FX,FY,FZ,FV,
     +                 DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,
     +                 DFVDLN,DFXDH,DFYDH,DFZDH,DFVDH, CALNM,IST)
C          Interpolation of x,y,z,v and their derivatives
C          940806 A. D. Richmond
C          INPUTS:
C            GLAT    = latitude, degrees
C            GLON    = longitude, degrees
C            ALT     = altitude, km
C            X,Y,Z,V = gridded arrays
C            NLAT,NLON,NALT = three dimensions of x,y,z,v and respective
C                      dimensions of GP___ arrays
C            GPLAT,GPLON,GPALT = grid point geographic locations
C            CALNM   = Name of calling routine (for error diagnostics)
C          OUTPUT:
C            FX = interpolated value of x
C            FY = interpolated value of y
C            FZ = interpolated value of z
C            FV = interpolated value of v
C            DFXDTH,DFYDTH,DFZDTH,DFVDTH = derivatives of x,y,z,v with
C                  respect to colatitude, in radians-1
C            DFXDLN,DFYDLN,DFZDLN,DFVDLN = derivatives of x,y,z,v with
C                  respect to longitude, in radians-1
C            DFXDH,DFYDH,DFZDH,DFVDH = derivatives of x,y,z,v with
C                  respect to altitude, in km-1
C            IST = Status =  0 = okay.
      DIMENSION X(NLAT,NLON,NALT), Y(NLAT,NLON,NALT), Z(NLAT,NLON,NALT),
     +          V(NLAT,NLON,NALT) , GPLAT(NLAT),GPLON(NLON),GPALT(NALT)
      CHARACTER*(*) CALNM

!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

C          Local declarations
      DATA IO, JO, KO / 1, 1, 1 /
      SAVE IO, JO, KO

      IENTRY = 0
      GO TO 5

C*******************************************************************************
      ENTRY INTRPSC (GLAT,GLON,ALT, X,Y,Z, NLAT,NLON,NALT,
     +               GPLAT,GPLON,GPALT, FX,FY,FZ , CALNM,IST)
C          Interpolation of x,y,z
C          940803 A. D. Richmond.
C          Inputs and outputs:  for definitions, see above.

      IENTRY = 1
    5 IST = 0

      IF (GLAT .LT. GPLAT(1) .OR. GLAT .GT. GPLAT(NLAT)) GO TO 9100
      IF (ALT  .LT. GPALT(1) .OR. ALT  .GT. GPALT(NALT)) GO TO 9300
C          Accept input longitude range +/- one revolution (Dec 98)
      IF (GLON .LT. GPLON(1)   ) GLON = GLON + 360.
      IF (GLON .GT. GPLON(NLON)) GLON = GLON - 360.
      IF (GLON .LT. GPLON(1) .OR. GLON .GT. GPLON(NLON)) GO TO 9200

      I = IO
      IF (GLAT .LE. GPLAT(I)) GO TO 15
   12 I = I + 1
      IF (GPLAT(I) .LT. GLAT) GO TO 12
      I = I - 1
      GO TO 16
   14 I = I - 1
   15 IF (GPLAT(I) .GT. GLAT) GO TO 14
   16 IO = I
      DLAT = GPLAT(I+1) - GPLAT(I)
      XI   = (GLAT - GPLAT(I)) / DLAT

      J = JO
      IF (GLON .LE. GPLON(J)) GO TO 25
   22 J = J + 1
      IF (GPLON(J) .LT. GLON) GO TO 22
      J = J - 1
      GO TO 26
   24 J = J - 1
   25 IF (GPLON(J) .GT. GLON) GO TO 24
   26 JO = J
      DLON = GPLON(J+1) - GPLON(J)
      YJ   = (GLON - GPLON(J)) / DLON

      K = KO
      IF (ALT .LE. GPALT(K)) GO TO 35
   32 K = K + 1
      IF (GPALT(K) .LT. ALT) GO TO 32
      K = K - 1
      GO TO 36
   34 K = K - 1
   35 IF (GPALT(K) .GT. ALT) GO TO 34
   36 KO = K
      HTI  = RE/(RE+ALT)
      DIHT = RE/(RE+GPALT(K+1)) - RE/(RE+GPALT(K))
      ZK  = (HTI - RE/(RE+GPALT(K))) / DIHT

C          For intrpsc:
      IF (IENTRY.EQ.1) THEN
	CALL TRILINS (X(I,J,K),NLAT,NLON,XI,YJ,ZK,FX)
	CALL TRILINS (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,FY)
	CALL TRILINS (Z(I,J,K),NLAT,NLON,XI,YJ,ZK,FZ)
        RETURN
      ENDIF
	
C          For intrp:
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,FX,DFXDN,DFXDE,DFXDD)
      DFXDTH = -DFXDN*RTOD/DLAT
      DFXDLN =  DFXDE*RTOD/DLON
      DFXDH  = -HTI*HTI*DFXDD/(RE*DIHT)
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,FY,DFYDN,DFYDE,DFYDD)
      DFYDTH = -DFYDN*RTOD/DLAT
      DFYDLN =  DFYDE*RTOD/DLON
      DFYDH  = -HTI*HTI*DFYDD/(RE*DIHT)
      CALL TRILIN (Z(I,J,K),NLAT,NLON,XI,YJ,ZK,FZ,DFZDN,DFZDE,DFZDD)
      DFZDTH = -DFZDN*RTOD/DLAT
      DFZDLN =  DFZDE*RTOD/DLON
      DFZDH  = -HTI*HTI*DFZDD/(RE*DIHT)
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,FV,DFVDN,DFVDE,DFVDD)
      DFVDTH = -DFVDN*RTOD/DLAT
      DFVDLN =  DFVDE*RTOD/DLON
      DFVDH  = -HTI*HTI*DFVDD/(RE*DIHT)

      IF (NLAT .LT. 3) RETURN

C          Improve calculation of longitudinal derivatives near poles
      IF (GLAT .GE. DLAT-90.) GO TO 40
      FAC = .5*XI
      OMFAC = 1. - FAC
      XI = XI - 1.
      I = I + 1
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFXDLN = DFXDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFYDLN = DFYDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFVDLN = DFVDLN*OMFAC + FAC*DMDFDE*RTOD/DLON

   40 IF (GLAT .LE. 90.-DLAT) GO TO 50
      FAC = .5*(1.- XI)
      OMFAC = 1. - FAC
      XI = XI + 1.
      I = I - 1
      CALL TRILIN (X(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFXDLN = DFXDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (Y(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFYDLN = DFYDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
      CALL TRILIN (V(I,J,K),NLAT,NLON,XI,YJ,ZK,DMF,DMDFDN,DMDFDE,DMDFDD)
      DFVDLN = DFVDLN*OMFAC + FAC*DMDFDE*RTOD/DLON
   50 RETURN

C          Error trap diagnostics
 9100 WRITE (MSGU,'(A,'':  Latitude out of range; GPLAT(1),GLAT,GPLAT(NL
     +AT)='',3F10.3)')  CALNM,GPLAT(1),GLAT,GPLAT(NLAT)
      IST = 1
      RETURN
 9200 WRITE (MSGU,'(A,'':  Longitude out of range; GPLON(1),GLON,GPLON(N
     +LON)='',3F10.3)') CALNM,GPLON(1),GLON,GPLON(NLON)
      IST = 1
      RETURN
 9300 WRITE (MSGU,'(A,'':  Altitude out of range; GPALT(1),ALT,GPALT(NAL
     +T)='',3F10.3)')   CALNM,GPALT(1),ALT,GPALT(NALT)
      IST = 1
      RETURN

      END subroutine intrp

C*******************************************************************************
      SUBROUTINE TRILIN (U,NLAT,NLON,XI,YJ,ZK,FU,DFUDX,DFUDY,DFUDZ)
C  Trilinear interpolation of u and its derivatives
C 940803 A. D. Richmond
C Inputs:
C   u(1,1,1) = address of lower corner of interpolation box 
C   nlat = first dimension of u from calling routine
C   nlon = second dimension of u from calling routine
C   xi = fractional distance across box in x direction 
C   yj = fractional distance across box in y direction 
C   zk = fractional distance across box in z direction 
C Outputs:
C   fu = interpolated value of u
C   dfudx = interpolated derivative of u with respect to i (x direction)
C   dfudy = interpolated derivative of u with respect to j (y direction)
C   dfudz = interpolated derivative of u with respect to k (z direction)
      DIMENSION U(NLAT,NLON,2)

      IENTRY = 0
      GOTO 5
C*******************************************************************************
      ENTRY TRILINS (U,NLAT,NLON,XI,YJ,ZK,FU)
C  Trilinear interpolation of u only
C 940803 A. D. Richmond
C Inputs and outputs:  see above for definitions
      IENTRY = 1

    5 CONTINUE
      OMXI = 1. - XI
      OMYJ = 1. - YJ
      OMZK = 1. - ZK

      FU = U(1,1,1)*OMXI*OMYJ*OMZK
     2   + U(2,1,1)*XI*OMYJ*OMZK
     3   + U(1,2,1)*OMXI*YJ*OMZK
     4   + U(1,1,2)*OMXI*OMYJ*ZK
     5   + U(2,2,1)*XI*YJ*OMZK
     6   + U(2,1,2)*XI*OMYJ*ZK
     7   + U(1,2,2)*OMXI*YJ*ZK
     8   + U(2,2,2)*XI*YJ*ZK

      IF (IENTRY.NE.0) RETURN

      DFUDX = (U(2,1,1)-U(1,1,1))*OMYJ*OMZK
     2      + (U(2,2,1)-U(1,2,1))*YJ*OMZK
     3      + (U(2,1,2)-U(1,1,2))*OMYJ*ZK
     4      + (U(2,2,2)-U(1,2,2))*YJ*ZK
      DFUDY = (U(1,2,1)-U(1,1,1))*OMXI*OMZK
     2      + (U(2,2,1)-U(2,1,1))*XI*OMZK
     3      + (U(1,2,2)-U(1,1,2))*OMXI*ZK
     4      + (U(2,2,2)-U(2,1,2))*XI*ZK
      DFUDZ = (U(1,1,2)-U(1,1,1))*OMXI*OMYJ
     2      + (U(2,1,2)-U(2,1,1))*XI*OMYJ
     3      + (U(1,2,2)-U(1,2,1))*OMXI*YJ
     4      + (U(2,2,2)-U(2,2,1))*XI*YJ
      RETURN
      END subroutine trilin

C*******************************************************************************
      SUBROUTINE ADPL(GLAT,GLON,CTH,STH,FX,FY,FZ,FV
     1  ,DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN)
C  v is used for vr2n
C  Add-back of pseudodipole component to x,y,z,v and their derivatives.
C  940715 A. D. Richmond
C Inputs:
C   glat = latitude, degrees
C   glon = longitude, degrees
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
C   fv = interpolated value of v
C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
C	colatitude, in radians-1
C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
C	longitude, in radians-1
C Output:
C   cth,sth = cos(colatitude), sin(colatitude)
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
C   fv = interpolated value of v
C   dfxdth,dfydth,dfzdth,dfvdth = derivatives of x,y,z,v with respect to 
C	colatitude, in radians-1
C   dfxdln,dfydln,dfzdln,dfvdln = derivatives of x,y,z,v with respect to 
C	longitude, in radians-1
!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
      CPH = COS((GLON-ELON)*DTOR)
      SPH = SIN((GLON-ELON)*DTOR)
      CTH = SIN(GLAT*DTOR)
      STH = COS(GLAT*DTOR)
      CTM = CTP*CTH + STP*STH*CPH
      FX = FX + STH*CTP*CPH - CTH*STP
      FY = FY + STH*SPH
      FZ = FZ + CTM
      FV = FV - CTM
      DFXDTH = DFXDTH + CTP*CTH*CPH + STP*STH
      DFYDTH = DFYDTH + CTH*SPH
      DFZDTH = DFZDTH - CTP*STH + STP*CTH*CPH
      DFVDTH = DFVDTH + CTP*STH - STP*CTH*CPH
      DFXDLN = DFXDLN - CTP*STH*SPH
      DFYDLN = DFYDLN + STH*CPH
      DFZDLN = DFZDLN - STP*STH*SPH
      DFVDLN = DFVDLN + STP*STH*SPH
      RETURN
      END subroutine adpl

C*******************************************************************************
      SUBROUTINE ADPLSC (GLAT,GLON,FX,FY,FZ)
C  Add-back of pseudodipole component to x,y,z 
C  940801 A. D. Richmond
C Inputs:
C   glat = latitude, degrees
C   glon = longitude, degrees
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
C Output:
C   fx = interpolated value of x
C   fy = interpolated value of y
C   fz = interpolated value of z
!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
      CPH = COS((GLON-ELON)*DTOR)
      SPH = SIN((GLON-ELON)*DTOR)
      CTH = SIN(GLAT*DTOR)
      STH = COS(GLAT*DTOR)
      CTM = CTP*CTH + STP*STH*CPH
      FX = FX + STH*CTP*CPH - CTH*STP
      FY = FY + STH*SPH
      FZ = FZ + CTM
      RETURN
      END subroutine adplsc

C*******************************************************************************
      SUBROUTINE GRADXYZV (ALT,CTH,STH,
     +          DFXDTH,DFYDTH,DFZDTH,DFVDTH,DFXDLN,DFYDLN,DFZDLN,DFVDLN,
     +          DFXDH,DFYDH,DFZDH,DFVDH,GRADX,GRADY,GRADZ,GRADV)
C          Calculates east,north,up components of gradients of x,y,z,v in
C          geodetic coordinates.  All gradients are in inverse km.  Assumes
C          flatness of 1/298.25 and equatorial radius (REQ) of 6378.16 km.
C          940803 A. D. Richmond
      COMMON /APXGEOD/ RHO,DDISDTH
      DIMENSION GRADX(3),GRADY(3),GRADZ(3),GRADV(3)
      IENTRY = 0
      GOTO 5
C*******************************************************************************
      ENTRY GRAPXYZV (ALT,CTH,STH,
     +              DFXDLN,DFYDLN,DFZDLN,DFVDLN,GRADX,GRADY,GRADZ,GRADV)
C          Calculates east component of gradient near pole.
C          940803 A. D. Richmond
C          Inputs and outputs:  see above for definitions
      IENTRY = 1

    5 CONTINUE
      D2 = 40680925.E0 - 272340.E0*CTH*CTH
C          40680925. = req**2 (rounded off)
C          272340.   = req**2 * E2, where E2 = (2. - 1./298.25)/298.25
C                      is square of eccentricity of ellipsoid.
      D = SQRT(D2)
      RHO = STH*(ALT + 40680925.E0/D)
      DDDTHOD = 272340.E0*CTH*STH/D2
      DRHODTH = ALT*CTH + (40680925.E0/D)*(CTH-STH*DDDTHOD)
      DZETDTH =-ALT*STH - (40408585.E0/D)*(STH+CTH*DDDTHOD)
      DDISDTH = SQRT(DRHODTH*DRHODTH + DZETDTH*DZETDTH)
      GRADX(1) = DFXDLN/RHO
      GRADY(1) = DFYDLN/RHO
      GRADZ(1) = DFZDLN/RHO
      GRADV(1) = DFVDLN/RHO

      IF (IENTRY .NE. 0) RETURN

      GRADX(2) = -DFXDTH/DDISDTH
      GRADY(2) = -DFYDTH/DDISDTH
      GRADZ(2) = -DFZDTH/DDISDTH
      GRADV(2) = -DFVDTH/DDISDTH
      GRADX(3) = DFXDH
      GRADY(3) = DFYDH
      GRADZ(3) = DFZDH
      GRADV(3) = DFVDH

      RETURN
      end subroutine gradxyzv

C*******************************************************************************
      SUBROUTINE GRADLPV (HR,ALT,FX,FY,FZ,FV,GRADX,GRADY,GRADZ,GRADV,
     +              XLATM,XLONM,VMP,GRCLM,CLMGRP,QDLAT,RGRLP,B,CLM,R3_2)
C          Uses gradients of x,y,z,v to compute geomagnetic field and
C          gradients of apex latitude, longitude.
C          940819 A. D. Richmond
C          INPUT:
C            HR     = reference altitude
C            ALT    = altitude
C            FX,FY,FZ,FV = interpolated values of x,y,z,v, plus pseudodipole
C                     component
C            GRADX,GRADY,GRADZ,GRADV = interpolated gradients of x,y,z,v,
C                     including pseudodipole components (east,north,up)
C          OUTPUT:
C            XLATM  = modified apex latitude (lambda_m), degrees
C            XLONM  = apex longitude (phi_a), degrees
C            VMP    = magnetic potential, in T.m.
C            GRCLM  = grad(cos(lambda_m)), in km-1
C            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
C            QDLAT  = quasi-dipole latitude, degrees
C            RGRLP  = (re + alt)*grad(lambda')
C            B      = magnetic field, in nT
C            CLM    = cos(lambda_m)
C            R3_2   = ((re + alt)/(re + hr))**(3/2)

      DIMENSION GRADX(3),GRADY(3),GRADZ(3),GRADV(3),
     +          GRCLM(3),CLMGRP(3),RGRLP(3),B(3)

      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA

      RR = RE + HR
      R  = RE + ALT
      RN = R/RE
      SQRROR = SQRT(RR/R)
      R3_2 = 1./SQRROR/SQRROR/SQRROR
      XLONM = ATAN2(FY,FX)
      CPM = COS(XLONM)
      SPM = SIN(XLONM)
      XLONM = RTOD*XLONM
      BO = VP*1.E6
C             1.E6 converts T to nT and km-1 to m-1.
      RN2 = RN*RN
      VMP = VP*FV/RN2
      B(1) = -BO*GRADV(1)/RN2
      B(2) = -BO*GRADV(2)/RN2
      B(3) = -BO*(GRADV(3)-2.*FV/R)/RN2

      X2PY2 = FX*FX + FY*FY
      XNORM = SQRT(X2PY2 + FZ*FZ)
      XLP = ATAN2(FZ,SQRT(X2PY2))
      SLP = SIN(XLP)
      CLP = COS(XLP)
      QDLAT = XLP*RTOD
      CLM = SQRROR*CLP
      IF (CLM.LE.1.) GOTO 5
      WRITE (6,*) 'Stopped in gradlpv because point lies below field lin      
     1e that peaks at reference height.'
      call shutdown('GRADLPV')
    5 XLATM = RTOD*ACOS(CLM)
C  If southern magnetic hemisphere, reverse sign of xlatm
      IF (SLP.LT.0.) XLATM = - XLATM
      DO 10 I=1,3
	GRCLP = CPM*GRADX(I) + SPM*GRADY(I)
	RGRLP(I) = R*(CLP*GRADZ(I) - SLP*GRCLP)
	GRCLM(I) = SQRROR*GRCLP
   10   CLMGRP(I) = SQRROR*(CPM*GRADY(I)-SPM*GRADX(I))
      GRCLM(3) = GRCLM(3) - SQRROR*CLP/(2.*R)
      RETURN
C*******************************************************************************
      ENTRY XYZ2APX (ALT,FX,FY,FZ,A,ALAT,ALON,IERR)
C          Computes apex latitude, longitude.
C          990309 A. D. Richmond
C          INPUT:
C            ALT      = altitude
C            FX,FY,FZ = interpolated values of x,y,z, plus pseudodipole
C                       component
C          OUTPUT:
C            A    = apex radius, normalized by req
C            ALAT = apex latitude, degrees
C            ALON = apex longitude, degrees
C
C          Mod (Mar 99):  Lines 19-30 are changed from the original in order
C          to avoid a truncation error near the magnetic equator.  What was
C          done is to make use of the identity
C
C                  SIN(ALAT/RTOD)**2 + COS(ALAT/RTOD)**2 = 1,
C
C          so that since
C
C                  COS(ALAT/RTOD)**2 = 1/A (Eq. 1),
C
C          then
C
C                  SIN(ALAT/RTOD)**2 = (A-1)/A (Eq. 2)
C
C          Below AM1 = A-1.  If AM1 is less than 1, use Eq. 1;
C          otherwise use Eq. 2.  Mathematically, both equations should
C          give identical results, but numerically it is better to use
C          that function ASIN or ACOS that has the smaller argument.
C          The jump from one equation to the other occurs at ALAT = 45.

      IERR  = 0
      ALON  = ATAN2(FY,FX)*RTOD
      SLP2  = FZ*FZ
      X2PY2 = FX*FX + FY*FY
      XNORM = SLP2 + X2PY2
      SLP2  = SLP2/XNORM
      CLP2  = X2PY2/XNORM
      AM1   = (RE*SLP2 + ALT)/(REQ*CLP2)
      A = 1. + AM1

      IF (AM1.LT.0.) THEN
        IERR = 1
        WRITE (6,*) 'Missing alat returned because point lies below fiel
     1d line that peaks at Earth surface.'
        RETURN
      ELSEIF (AM1.LT.1.) THEN
        ALAT = RTOD*ASIN(SQRT(AM1/A))
      ELSE
        ALAT = RTOD*ACOS(1./SQRT(A))
      ENDIF
C  If southern magnetic hemisphere, reverse sign of alat
      IF (FZ.LT.0.) ALAT = - ALAT
      RETURN
      end subroutine gradlpv

C*******************************************************************************
      SUBROUTINE BASVEC (HR,XLATM,GRCLM,CLMGRP,RGRLP,B,CLM,R3_2,
     +                  BMAG,SIM,SI,F,D,W,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
C          Computes base vectors and other parameters for apex coordinates.
C          Vector components:  east, north, up
C          940801 A. D. Richmond
C          Reference:
C            Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
C            Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
C          INPUTS:
C            HR     = reference altitude
C            XLATM  = modified apex latitude, degrees
C            GRCLM  = grad(cos(lambda_m)), in km-1
C            CLMGRP = cos(lambda_m)*grad(phi_a), in km-1
C            RGRLP  = (re + altitude)*grad(lambda')
C            B      = magnetic field, in nT
C            CLM    = cos(lambda_m)
C            R3_2   = ((re + altitude)/(re + hr))**(3/2)
C          RETURNS:
C            BMAG    = magnitude of magnetic field, in nT
C            SIM     = sin(I_m) of article
C            SI      = sin(I)
C            F       = F of article
C            D       = D of article
C            W       = W of article
C            BHAT    = unit vector along geomagnetic field direction
C            D1...F2 = base vectors of article
!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
      DIMENSION GRCLM(3),CLMGRP(3),RGRLP(3),B(3)
      DIMENSION BHAT(3),D1(3),D2(3),D3(3),E1(3),E2(3),E3(3),F1(2),F2(2)        
      RR = RE + HR
      SIMOSLM = 2./SQRT(4. - 3.*CLM*CLM)
      SIM = SIMOSLM*SIN(XLATM*DTOR)
      BMAG = SQRT(B(1)*B(1) + B(2)*B(2) + B(3)*B(3))
      D1DB = 0.
      D2DB = 0.
      DO 10 I=1,3
        BHAT(I) = B(I)/BMAG
        D1(I) = RR*CLMGRP(I)
        D1DB = D1DB + D1(I)*BHAT(I)
        D2(I) = RR*SIMOSLM*GRCLM(I)
   10   D2DB = D2DB + D2(I)*BHAT(I)
C Ensure that d1,d2 are exactly perpendicular to B:
      DO 15 I=1,3
        D1(I) = D1(I) - D1DB*BHAT(I)
   15   D2(I) = D2(I) - D2DB*BHAT(I)
      E3(1) = D1(2)*D2(3) - D1(3)*D2(2)
      E3(2) = D1(3)*D2(1) - D1(1)*D2(3)
      E3(3) = D1(1)*D2(2) - D1(2)*D2(1)
      D = BHAT(1)*E3(1) + BHAT(2)*E3(2) + BHAT(3)*E3(3)
      DO 20 I=1,3
        D3(I) = BHAT(I)/D
C Following step may be unnecessary, but it ensures that e3 lies along bhat.
        E3(I) = BHAT(I)*D
   20   CONTINUE
      E1(1) = D2(2)*D3(3) - D2(3)*D3(2)
      E1(2) = D2(3)*D3(1) - D2(1)*D3(3)
      E1(3) = D2(1)*D3(2) - D2(2)*D3(1)
      E2(1) = D3(2)*D1(3) - D3(3)*D1(2)
      E2(2) = D3(3)*D1(1) - D3(1)*D1(3)
      E2(3) = D3(1)*D1(2) - D3(2)*D1(1)
      W = RR*RR*CLM*ABS(SIM)/(BMAG*D)
      SI = -BHAT(3)
      F1(1) =  RGRLP(2) 
      F1(2) = -RGRLP(1)
      F2(1) = -D1(2)*R3_2
      F2(2) =  D1(1)*R3_2
      F = F1(1)*F2(2) - F1(2)*F2(1)
      RETURN
      end subroutine basvec

C*******************************************************************************
      SUBROUTINE CKGP (CALNM,MSGUN,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,IST)
C          Check grid point values tests extremes and order of the grid
C          point arrays, producing diagnostics to MSGUN and IST=1 when
C          rules have been broken.
      DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT)
      CHARACTER*(*) CALNM

      IST = 1
      OLAT = -90.
      DO 10 I=1,NLAT
      IF (ABS (GPLAT(I)) .GT.  90.) GO TO 9100
      IF (     GPLAT(I)  .LT. OLAT) GO TO 9200
   10 OLAT = GPLAT(I)
!
! 2/00: reset longitude limits to -270 -> +270, as per communication
!       with Richmond, 2/25/00.
!
!     OLON = -180.
      OLON = -270.
      DO 20 I=1,NLON
!     IF (ABS (GPLON(I)) .GT. 180.) GO TO 9300
      IF (ABS (GPLON(I)) .GT. 270.) GO TO 9300
      IF (     GPLON(I)  .LT. OLON) GO TO 9400
   20 OLON = GPLON(I)

      OALT = 0.
      DO 30 I=1,NALT
      IF (GPALT(I) .LT.   0.) GO TO 9500
      IF (GPALT(I) .LT. OALT) GO TO 9600
   30 OALT = GPALT(I)

      IST = 0

  100 RETURN

 9100 WRITE (MSGUN,'(A,'':  |GPLAT(I)| > 90; I,GPLAT(I)'',I5,F10.3)')
     +            CALNM,                     I,GPLAT(I)
      GO TO 100
 9200 WRITE (MSGUN,'(A,'':  GPLAT(I) < GPLAT(I-1); I,GPLAT(I),GPLAT(I-1)
     +='',I5,2F10.3)')                     CALNM,  I,GPLAT(I),OLAT
      GO TO 100
 9300 WRITE (MSGUN,'(A,'':  |GPLON(I)| > 180; I,GPLON(I)'',I5,F10.3)')
     +                                 CALNM, I,GPLON(I)
      GO TO 100
 9400 WRITE (MSGUN,'(A,'':  GPLON(I) < GPLON(I-1); I,GPLON(I),GPLON(I-1)
     +='',I5,2F10.3)')                    CALNM, I,GPLON(I),OLON
      GO TO 100
 9500 WRITE (MSGUN,'(A,'':  GPALT(I) <  0; I,GPALT(I)'',I5,F10.3)')
     +                                CALNM, I,GPALT(I)
      GO TO 100
 9600 WRITE (MSGUN,'(A,'':  GPALT(I) < GPALT(I-1); I,GPALT(I),GPALT(I-1)
     +='',I5,2F10.3)')                      CALNM, I,GPALT(I),OALT
      GO TO 100
      END subroutine ckgp

C*******************************************************************************
      SUBROUTINE MAKEXYZV (EPOCH,NLAT,NLON,NALT,GPLAT,GPLON,GPALT,
     +                                                          X,Y,Z,V)
C          Sets up grid arrays for later interpolation
C          940822 A. D. Richmond, NCAR
C          INPUT:
C            EPOCH = year and fraction (e.g., 1994.50 for 1994 July 2)
C            NLAT,NLON,NALT = triple dimensions of X,Y,Z,V and respective
C                    single dimensions of GP___ arrays
C            GPLAT,GPLON,GPALT = grid point latitudes, longitudes and altitudes
C          OUTPUT:
C            X = array containing cos(lambda')cos(phi_a) less pseudodipole
C                component
C            Y = array containing cos(lambda')sin(phi_a) less pseudodipole
C                component
C            Z = array containing sin(lambda') less pseudodipole component
C            V = array containing ((magnetic potential)/vp)*((re+height)/re)**2,
C                less pseudodipole component
C
C          Modification (99 Mar):  Make X,Y,Z,V constant near the poles
C          for all GPLON(j) at each height.  Add POLA to APXCON

      DIMENSION X(NLAT,NLON,NALT), Y(NLAT,NLON,NALT), Z(NLAT,NLON,NALT),
     +          V(NLAT,NLON,NALT), GPLAT(NLAT), GPLON(NLON), GPALT(NALT)


C            POLA       = Pole angle (deg); when the geographic latitude is
C                         poleward of POLA, X,Y,Z,V are forced to be constant.
C                         for all longitudes at each altitude.  POLA is defined
C                         in APXMKA (POLA = 90. - SQRT (PRECISE) * RTOD),
C                         which currently makes POLA = 89.995
!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP

      CALL COFRM (EPOCH)
      CALL DYPOL (COLAT,ELON,VP)
      CTP = COS (COLAT*DTOR)
      STP = SIN (COLAT*DTOR)
      REQORE = REQ/RE
      RQORM1 = REQORE-1.

      DO 100 I=1,NLAT
      CT = SIN (GPLAT(I)*DTOR)
      ST = COS (GPLAT(I)*DTOR)
      KPOL = 0
      IF (ABS (GPLAT(I)) .GT. POLA) KPOL = 1

      DO 100 J=1,NLON
      IF (KPOL .EQ. 0) GO TO 20
      IF (J    .EQ. 1) GO TO 20
C          KPOL = 1 (poleward of POLA) and first lon's XYZV are defined
      DO 10 K=1,NALT
      V(I,J,K) = V(I,1,K)
      X(I,J,K) = X(I,1,K)
      Y(I,J,K) = Y(I,1,K)
   10 Z(I,J,K) = Z(I,1,K)
      GO TO 100

   20 CP  = COS ((GPLON(J)-ELON)*DTOR)
      SP  = SIN ((GPLON(J)-ELON)*DTOR)
C           ctm   is pseudodipole component of z
C          -ctm   is pseudodipole component of v
C          stmcpm is pseudodipole component of x
C          stmspm is pseudodipole component of y
      CTM    = CTP*CT + STP*ST*CP
      STMCPM = ST*CTP*CP - CT*STP
      STMSPM = ST*SP

      DO 30 K=1,NALT
      CALL APEX (EPOCH,GPLAT(I),GPLON(J),GPALT(K),
     +           A,ALAT,PHIA,BMAG,XMAG,YMAG,ZDOWN,VMP)
      VNOR = VMP/VP
      RP = 1. + GPALT(K)/RE
      V(I,J,K) = VNOR*RP*RP + CTM
      REQAM1 = REQ*(A-1.)
      SLP = SQRT(max(REQAM1-GPALT(K),0.)/(REQAM1+RE))
C          Reverse sign of slp in southern magnetic hemisphere
      IF (ZDOWN.LT.0.) SLP = -SLP
      CLP = SQRT (RP/(REQORE*A-RQORM1))
      PHIAR = PHIA*DTOR
      X(I,J,K) = CLP*COS (PHIAR) - STMCPM
      Y(I,J,K) = CLP*SIN (PHIAR) - STMSPM
      Z(I,J,K) = SLP - CTM
   30 CONTINUE

  100 CONTINUE

      RETURN
      end subroutine makexyzv

C*******************************************************************************

      SUBROUTINE SETMISS (XMISS
     2 ,XLATM,XLONM,VMP,B,BMAG,BE3,SIM,SI,F,D,W
     3 ,BHAT,D1,D2,D3,E1,E2,E3,F1,F2)
      DIMENSION BHAT(3),D1(3),D2(3),D3(3),E1(3),E2(3),E3(3),F1(2),F2(2)        
     2 ,B(3)
      XLATM = XMISS
      XLONM = XMISS
      VMP = XMISS
      BMAG = XMISS
      BE3 = XMISS
      SIM = XMISS
      SI = XMISS
      F = XMISS
      D = XMISS
      W = XMISS
      DO 5 I=1,3
	B(I) = XMISS
	BHAT(I) = XMISS
	D1(I) = XMISS
	D2(I) = XMISS
	D3(I) = XMISS
	E1(I) = XMISS
	E2(I) = XMISS
    5   E3(I) = XMISS
      DO 6 I=1,2
	F1(I) = XMISS
    6   F2(I) = XMISS
      RETURN
      END subroutine setmiss

C*******************************************************************************
      SUBROUTINE GM2GC (GMLAT,GMLON,GCLAT,GCLON)
C  Converts geomagnetic to geocentric coordinates.
C  940819 A. D. Richmond
C
C  Inputs:
C	gmlat = geomagnetic latitude in degrees
C	gmlon = geomagnetic longitude in degrees
C  Outputs:
C	gclat = geocentric latitude in degrees
C	gclon = geocentric longitude in degrees
C
C  Common/consts/
C	rtod, dtor = 180/pi, pi/180
C	re, req    = 6371.2, 6378.160
!
! bf: separate ints from reals in commons:
      COMMON /APXCON_int/ MSGU
      COMMON /APXCON/ RTOD,DTOR,RE,REQ,POLA
C  Common/dipol/
C       colat = geocentric colatitude of north geomagnetic pole, in degrees
C	elon  = geocentric east longitude of north geomagnetic pole, in degrees
C	vp    = magnetic potential at 1 RE, south geomagnetic pole
C	ctp   = cos(colat*dtor)
C	stp   = sin(colat*dtor)
      COMMON /APXDIPL/ COLAT,ELON,VP,CTP,STP
C
      STM = COS(GMLAT*DTOR)
      CTM = SIN(GMLAT*DTOR)
      CTC = CTP*CTM - STP*STM*COS(GMLON*DTOR)
      CTC = min(CTC,1.)
      CTC = max(CTC,-1.)
      GCLAT = ASIN(CTC)*RTOD
      GCLON = ATAN2(STP*STM*SIN(GMLON*DTOR),CTM-CTP*CTC)
      GCLON = GCLON*RTOD + ELON
      IF (GCLON.LT.-180.) GCLON = GCLON + 360.
      RETURN
      END subroutine gm2gc
!-----------------------------------------------------------------------
!
! Original file ~bozo/pgms/apex/apex.f copied 2/25/00.
!
      SUBROUTINE APEX (DATE,DLAT,DLON,ALT,
     +                A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,V)
C Calculate apex radius, latitude, longitude; and magnetic field and potential
C 940822 A. D. Richmond
C
C          INPUTS:
C            DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
C            DLAT = Geodetic latitude in degrees
C            DLON = Geodetic longitude in degrees
C            ALT = Altitude in km
C
C          RETURNS:
C            A    = (Apex height + REQ)/REQ, where REQ = equatorial Earth
C                    radius.  A is analogous to the L value in invariant
C                    coordinates.
C            ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
C            ALON = Apex longitude (Geomagnetic longitude of apex)
C            BMAG = Magnitude of geomagnetic field in nT
C            XMAG,YMAG,ZMAG = North, east, and downward geomagnetic field components in nT
C            V    = Magnetic potential, in T.m
C
C     COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
C          COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
C          ELON  = East longitude of geomagnetic dipole north pole (deg)
C          VP    = Magnitude, in T.m, of dipole component of magnetic potential at
C                  geomagnetic pole and geocentric radius of 6371.2 km
C          CTP,STP = cosine, sine of COLAT
C
C          MODIFICATIONS:
C          May 1999:  Revise DS calculation in LINAPX to avoid divide by zero.

      PARAMETER (RTOD=5.72957795130823E1, DTOR=1.745329251994330E-2,
     +           RE=6371.2, REQ=6378.160)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

      CALL COFRM (DATE)
      CALL DYPOL (CLATP,POLON,VPOL)
      COLAT = CLATP
      CTP   = COS(CLATP*DTOR)
      STP   = SQRT(1. - CTP*CTP)
      ELON  = POLON
      VP    = VPOL
      CALL LINAPX (DLAT,DLON,ALT,A,ALAT,ALON,XMAG,YMAG,ZMAG,BMAG)
      XMAG = XMAG*1.E5
      YMAG = YMAG*1.E5
      ZMAG = ZMAG*1.E5
      BMAG = BMAG*1.E5
      CALL GD2CART (DLAT,DLON,ALT,X,Y,Z)
      CALL FELDG (3,X/RE,Y/RE,Z/RE,BX,BY,BZ,V)
      RETURN
      END subroutine apex

       SUBROUTINE LINAPX(GDLAT,GLON,ALT
     2 ,A,ALAT,ALON,XMAG,YMAG,ZMAG,F)
C***BEGIN PROLOGUE  LINAPX                                                      
C***DATE WRITTEN   731029   (YYMMDD)                                            
C***AUTHOR  CLARK, W., N.O.A.A. ERL LAB.                                        
C***REVISION DATE  880201   (YYMMDD) Harsh Anand Passi, NCAR
C***LATEST REVISION 940803 A. D. Richmond
C***PURPOSE  Transforms the geographic coordinates to apex coordinates.         
C***DESCRIPTION                                                                 
C     The method used is as follow:                                             
C       1. Calculates step size as a function of the geomagnetic                
C          dipole coordinates of the starting point.                            
C       2. Determine direction of trace                                         
C       3. Convert the geodetic coordinates of the starting point               
C          to the cartesian coordinates for tracing.                            
C       Loop:                                                                   
C       i)   Increment step count, if count > 200,                              
C            assume it is dipole field, call DIPAPX to                          
C            determine Apex coordinates else continue.                          
C       ii)  Get field components in X, Y and Z direction                       
C       iii) Call ITRACE to trace field line.                                   
C       iv)  Test if Apex passed call FNDAPX to determine Apex coordinates      
C            else loop:                                                         
C                                                                               
C   INPUT                                                                       
C     GDLAT  Latitude of starting point (Deg)                                   
C     GLON   Longitude (East=+) of starting point (Deg)                         
C     ALT    Ht. of starting point (Km)                                         
C                                                                               
C  OUTPUT
C          A  (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C             A is analogous to the L value in invariant coordinates.
C       ALAT  Apex Lat. (deg)                                                   
C       ALON  Apex Lon. (deg)                                                   
C       XMAG  North component of magnetic field at starting point
C       YMAG  East component of magnetic field at starting point
C       ZMAG  Down component of magnetic field at starting point
C          F  Magnetic field magnitude at starting point
C
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C       YAPX    Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
C     COLAT   Geographic colatitude of the north pole of the                    
C             earth-centered dipole (Deg).                                      
C     ELON    Geographic longitude of the north pole of the                     
C             earth-centered dipole (Deg).                                      
C     VP      Magnetic potential magnitude at geomagnetic pole (T.m)
C     CTP     cos(COLAT*DTOR)
C     STP     sin(COLAT*DTOR)
C                                                                               
C     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
C     NSTP      Step count. Incremented in subprogram LINAPX.                   
C     Y         Array containing current tracing point cartesian coordinates.   
C     YOLD      Array containing previous tracing point cartesian coordinates.
C     SGN       Determines direction of trace. Set in subprogram LINAPX         
C     DS        Step size (Km) Computed in subprogram LINAPX.                   
C                                                                               
C     /FLDCOMD/ BX, BY, BZ, BB                                                  
C     BX        X comp. of field vector at the current tracing point (Gauss)    
C     BY        Y comp. of field vector at the current tracing point (Gauss)    
C     BZ        Z comp. of field vector at the current tracing point (Gauss)    
C     BB        Magnitude of field vector at the current tracing point (Gauss)  
C                                                                               
C***REFERENCES  Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17             
C                 (1971) GSFC, Greenbelt, Maryland                              
C                                                                               
C***ROUTINES CALLED  GD2CART,CONVRT,FELDG,
C                    ITRACE,DIPAPX,FNDAPX                                       
C                                                                               
C***END PROLOGUE  LINAPX                                                        
      PARAMETER (MAXS = 200,
     +           RTOD=5.72957795130823E1, DTOR=1.745329251994330E-2,
     +           RE=6371.2, REQ=6378.160)
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /APXIN/   YAPX(3,3)
      COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
!
! bf: separate ints and reals into 2 different commons:
!     COMMON /ITRA/    NSTP, Y(3), YP(3),  SGN, DS
      COMMON /ITRA_int/    NSTP
      COMMON /ITRA/  Y(3), YP(3),  SGN, DS

C          Determine step size as a function of geomagnetic dipole
C          coordinates of the starting point
      CALL CONVRT (2,GDLAT,ALT,GCLAT,R)
C     SINGML = .98*SIN(GCLAT*DTOR) + .199*COS(GCLAT*DTOR)*
      SINGML = CTP*SIN(GCLAT*DTOR) +  STP*COS(GCLAT*DTOR)*
     +                                             COS((GLON-ELON)*DTOR)

C          May 99: avoid possible divide by zero (when SINGML = 1.)
      CGML2 = max (0.25,1.-SINGML*SINGML)
      DS = .06*R/CGML2 - 370.
Cold:      Limit DS to its value at 60 deg gm latitude.
Cold: DS = .06*R/(1.-SINGML*SINGML) - 370.
Cold: IF (DS .GT. 1186.) DS = 1186.

C          Initialize YAPX array
      DO 4 J=1,3
      DO 4 I=1,3
    4 YAPX(I,J) = 0.

C          Convert from geodetic to earth centered cartesian coordinates
      CALL GD2CART (GDLAT,GLON,ALT,Y(1),Y(2),Y(3))
      NSTP = 0

C          Get magnetic field components to determine the direction for
C          tracing the field line
      CALL FELDG (1,GDLAT,GLON,ALT,XMAG,YMAG,ZMAG,F)
      SGN = SIGN (1.,-ZMAG)

C          Use cartesian coordinates to get magnetic field components
C          (from which gradients steer the tracing)
   10 IENTRY=2

      CALL FELDG (IENTRY,Y(1)/RE,Y(2)/RE,Y(3)/RE,BX,BY,BZ,BB)

      NSTP = NSTP + 1

C          Quit if too many steps
      IF (NSTP .GE. MAXS) THEN                                                 
        RHO = SQRT(Y(1)*Y(1) + Y(2)*Y(2))
        CALL CONVRT(3,XLAT,HT,RHO,Y(3))
        XLON = RTOD*ATAN2(Y(2),Y(1))
	CALL FELDG (1,XLAT,XLON,HT,BNRTH,BEAST,BDOWN,BABS)
	CALL DIPAPX (XLAT,XLON,HT,BNRTH,BEAST,BDOWN,A,ALON)
        ALAT = -SGN*RTOD*ACOS(SQRT(1./A))
        RETURN                                                                
      END IF                                                                   

C          Find next point using adams algorithm after 7 points
       CALL ITRACE (IAPX)
       IF (IAPX .EQ. 1) GO TO 10

C          Maximum radius just passed.  Find apex coords
      CALL FNDAPX (ALT,ZMAG,A,ALAT,ALON)
      RETURN
      end subroutine linapx

      SUBROUTINE ITRACE (IAPX)
      SAVE
C***BEGIN PROLOGUE  ITRACE                                                      
C***DATE WRITTEN   731029   (YYMMDD)                                            
C***AUTHOR  CLARK, W. N.O.A.A. ERL LAB.                                         
C***REVISION DATE  880201, H. Passi 
C***PURPOSE  Field line integration routine.                                    
C***DESCRIPTION                                                                 
C     It uses 4-point ADAMS formula after initialization.                       
C     First 7 iterations advance point by 3 steps.                              
C                                                                               
C   INPUT                                                                       
C            Passed in through the common blocks ITRA, FLDCOMD.                 
C            See the description below.                                         
C   OUTPUT                                                                      
C    IAPX    Flag set to 2 when APEX passed, otherwise set to 1.                
C                                                                               
C            Passed out through the common block APXIN.                         
C            See the description below.                                         
C                                                                               
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C     /FLDCOMD/ BX, BY, BZ, BB                                                  
C     BX        X comp. of field vector at the current tracing point (Gauss)    
C     BY        Y comp. of field vector at the current tracing point (Gauss)    
C     BZ        Z comp. of field vector at the current tracing point (Gauss)    
C     BB        Magnitude of field vector at the current tracing point (Gauss)  
C                                                                               
C     /ITRA/ NSTP, Y(3), YOLD(3), SGN, DS                                       
C     NSTP      Step count for line integration.                                
C               Incremented in subprogram LINAPX.                               
C     Y         Array containing current tracing point cartesian coordinates.   
C     YOLD      Array containing previous tracing point cartesian coordinates.  
C     SGN       Determines direction of trace. Set in subprogram LINAPX         
C     DS        Integration step size (arc length Km)                           
C               Computed in subprogram LINAPX.                                  
C                                                                               
C***REFERENCES  reference 1                                                     
C                                                                               
C***ROUTINES CALLED  None                                                       
C***COMMON BLOCKS    APXIN,FLDCOMD,ITRA                                         
C***END PROLOGUE  ITRACE                                                        
C **                                                                            
!
! bf: separate ints and reals into 2 different commons:
!     COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
      COMMON /ITRA_int/    NSTP
      COMMON /ITRA/ Y(3), YOLD(3), SGN, DS
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /APXIN/   YAPX(3,3)
      DIMENSION  YP(3, 4)

C          Statement function
      RDUS(D,E,F) = SQRT( D**2 + E**2 + F**2 )
C                                                                               
      IAPX = 1
C          Field line is defined by the following differential equations
C          in cartesian coordinates:
      YP(1, 4) = SGN*BX/BB
      YP(2, 4) = SGN*BY/BB
      YP(3, 4) = SGN*BZ/BB
      IF (NSTP .GT. 7) GO TO 90

C          First seven steps use this block
      DO 80 I = 1, 3
      GO TO (10, 20, 30, 40, 50, 60, 70) NSTP
   10 D2  = DS/2.
      D6  = DS/6.
      D12 = DS/12.
      D24 = DS/24.
      YP(I,1)   = YP(I,4)
      YOLD(I)   = Y(I)
      YAPX(I,1) = Y(I)
      Y(I) = YOLD(I) + DS*YP(I, 1)
      GO TO 80

   20 YP(I, 2) = YP(I,4)
      Y(I) = YOLD(I) + D2*(YP(I,2)+YP(I,1))
      GO TO 80

   30 Y(I) = YOLD(I) + D6*(2.*YP(I,4)+YP(I,2)+3.*YP(I,1))
      GO TO 80

   40 YP(I,2)  = YP(I,4)
      YAPX(I,2)= Y(I)
      YOLD(I)  = Y(I)
      Y(I)     = YOLD(I) + D2*(3.*YP(I,2)-YP(I,1))
      GO TO 80

   50 Y(I) = YOLD(I) + D12*(5.*YP(I,4)+8.*YP(I,2)-YP(I,1))
      GO TO 80

   60 YP(I,3)  = YP(I,4)
      YOLD(I)  = Y(I)
      YAPX(I,3)= Y(I)
      Y(I)     = YOLD(I) + D12*(23.*YP(I,3)-16.*YP(I,2)+5.*YP(I,1))
      GO TO 80

   70 YAPX(I, 1) = YAPX(I, 2)
      YAPX(I, 2) = YAPX(I, 3)
      Y(I) = YOLD(I) + D24*(9.*YP(I,4)+19.*YP(I,3)-5.*YP(I,2)+YP(I,1))
      YAPX(I, 3) = Y(I)
   80 CONTINUE

C          Signal if apex passed
      IF ( NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN
	RC = RDUS( YAPX(1,3), YAPX(2,3), YAPX(3,3))
        RP = RDUS( YAPX(1,2), YAPX(2,2), YAPX(3,2))                             
	IF (RC .LT. RP) IAPX=2
      ENDIF
      RETURN

C          Stepping block for NSTP .gt. 7
   90 DO 110 I = 1, 3                                                           
      YAPX(I, 1) = YAPX(I, 2)
      YAPX(I, 2) = Y(I)
      YOLD(I) = Y(I)
      TERM = 55.*YP(I, 4) - 59.*YP(I, 3) + 37.*YP(I, 2) - 9.*YP(I, 1)
      Y(I) = YOLD(I) + D24*TERM
      YAPX(I, 3) = Y(I)

      DO 100 J = 1, 3
      YP(I, J) = YP(I, J+1)
  100 CONTINUE
  110 CONTINUE                                                                  

      RC = RDUS (   Y(1),    Y(2),    Y(3))
      RP = RDUS (YOLD(1), YOLD(2), YOLD(3))
      IF (RC .LT. RP) IAPX=2

      RETURN
      END subroutine itrace

       SUBROUTINE FNDAPX (ALT,ZMAG,A,ALAT,ALON)
C***BEGIN PROLOGUE  FNDAPX                                                      
C***DATE WRITTEN   731023   (YYMMDD)                                            
C***AUTHOR  CLARK, W., NOAA BOULDER                                             
C***REVISION DATE  940803, A. D. Richmond, NCAR 
C***PURPOSE  Finds apex coords once tracing has signalled that the apex         
C            has been passed.  
C***DESCRIPTION                                                                 
C                                                                               
C     It uses second degree interpolation routine, FINT, to find                
C     apex latitude and apex longtitude.                                        
C   INPUT                                                                       
C     ALT    Altitude of starting point
C     ZMAG   Downward component of geomagnetic field at starting point
C     NMAX   Order of IGRF coefficients being used
C     G      Array of coefficients from COFRM
C   OUTPUT                                                                      
C          A  (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
C             A is analogous to the L value in invariant coordinates.
C       ALAT  Apex Lat. (deg)                                                   
C       ALON  Apex Lon. (deg)                                                   
C                                                                               
C***LONG DESCRIPTION                                                            
C                                                                               
C     COMMON Blocks Used                                                        
C     /APXIN/ YAPX(3,3)                                                         
C     YAPX      Matrix of cartesian coordinates (loaded columnwise)             
C               of the 3 points about APEX. Set in subprogram ITRACE.           
C                                                                               
C   /DIPOLE/COLAT,ELON,VP,CTP,STP                           
C     COLAT   Geocentric colatitude of the north pole of the                    
C             earth-centered dipole (Deg).                                      
C     ELON    Geographic longitude of the north pole of the                     
C             earth-centered dipole (Deg).                                      
C     CTP     cos(COLAT*DTOR)
C     STP     sin(COLAT*DTOR)
C                                                                               
C***ROUTINES CALLED  FINT                                                       
C***COMMON BLOCKS    APXIN,DIPOLE                                         
C***END PROLOGUE  FNDAPX                                                        
C                                                                               
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
      COMMON /APXIN/  YAPX(3,3)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
      DIMENSION Z(3), HT(3), Y(3)

C       ****
C       ****     GET GEODETIC FIELD COMPONENTS
C       ****
      IENTRY = 1
      DO 2 I = 1,3
	RHO = SQRT(YAPX(1,I)**2+YAPX(2,I)**2)
	CALL CONVRT (3,GDLT,HT(I),RHO,YAPX(3,I))
	GDLN = RTOD*ATAN2 (YAPX(2,I),YAPX(1,I))
	CALL FELDG (IENTRY,GDLT,GDLN,HT(I),X,YDUM,Z(I),F)
    2 CONTINUE
C     ****
C     ****     FIND CARTESIAN COORDINATES AT DIP EQUATOR BY INTERPOLATION
C     ****
      DO 3 I = 1,3
	CALL FINT(Z(1),Z(2),Z(3),YAPX(I,1),YAPX(I,2),YAPX(I,3),0.,Y(I))
    3 CONTINUE
C     ****
C     ****     FIND APEX HEIGHT BY INTERPOLATION
C     ****
      CALL FINT(Z(1),Z(2),Z(3),HT(1),HT(2),HT(3),0.,XINTER)
C          Ensure that XINTER is not less than original starting altitude (ALT)
      XINTER = max(ALT,XINTER)
      A = (REQ+XINTER)/(REQ)
C     ****
C     ****     FIND APEX COORDINATES , GIVING ALAT SIGN OF DIP AT
C     ****       STARTING POINT.  ALON IS THE VALUE OF THE GEOMAGNETIC
C     ****       LONGITUDE AT THE APEX.
C     ****
      IF (A.LT.1.) THEN
	WRITE(6,20)
  20    FORMAT (' BOMBED! THIS MAKES A LESS THAN ONE')
	call shutdown('fndapx')
      ENDIF
      RASQ = RTOD*ACOS(SQRT(1./A))
      ALAT = SIGN(RASQ,ZMAG)

C Algorithm for ALON:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let XLON be longitude of apex.
C   Let TE be colatitude of apex.
C   Let ANG be longitude angle from GM to apex.
C   Let TP be colatitude of GM.
C   Let TF be arc length between GM and apex.
C   Let PA = ALON be geomagnetic longitude, i.e., Pi minus angle measured 
C     counterclockwise from arc GM-apex to arc GM-GP.
C   Then, using notation C=cos, S=sin, spherical-trigonometry formulas 
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      XLON = ATAN2(Y(2),Y(1))
      ANG  = XLON-ELON*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      R    = SQRT(Y(1)**2+Y(2)**2+Y(3)**2)
      CTE  = Y(3)/R
      STE  = SQRT(1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END subroutine fndapx 

      SUBROUTINE DIPAPX(GDLAT,GDLON,ALT,BNORTH,BEAST,BDOWN,A,ALON)
C Compute a, alon from local magnetic field using dipole and spherical approx.
C 940501 A. D. Richmond
C Input:
C   GDLAT  = geodetic latitude, degrees
C   GDLON  = geodetic longitude, degrees
C   ALT    = altitude, km
C   BNORTH = geodetic northward magnetic field component (any units)
C   BEAST  = eastward magnetic field component
C   BDOWN  = geodetic downward magnetic field component
C Output:
C   A      = apex radius, 1 + h_A/R_eq
C   ALON   = apex longitude, degrees
C
C Algorithm:
C   Use spherical coordinates.
C   Let GP be geographic pole.
C   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
C   Let G be point at GDLAT,GDLON.
C   Let E be point on sphere below apex of dipolar field line passing through G.
C   Let TD be dipole colatitude of point G, found by applying dipole formula
C     for dip angle to actual dip angle.
C   Let B be Pi plus local declination angle.  B is in the direction 
C     from G to E.
C   Let TG be colatitude of G.
C   Let ANG be longitude angle from GM to G.
C   Let TE be colatitude of E.
C   Let TP be colatitude of GM.
C   Let A be longitude angle from G to E.
C   Let APANG = A + ANG
C   Let PA be geomagnetic longitude, i.e., Pi minus angle measured 
C     counterclockwise from arc GM-E to arc GM-GP.
C   Let TF be arc length between GM and E.
C   Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry formulas 
C     for the functions of the angles are as shown below.  Note: STFCPA,
C     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      PARAMETER (RE=6371.2,REQ=6378.160)
      COMMON/DIPOLE/COLAT,ELON,VP,CTP,STP
      BHOR = SQRT(BNORTH*BNORTH + BEAST*BEAST)
      IF (BHOR.EQ.0.) THEN
	ALON = 0.
	A = 1.E34
	RETURN
      ENDIF
      COTTD = BDOWN*.5/BHOR
      STD = 1./SQRT(1.+COTTD*COTTD)
      CTD = COTTD*STD
      SB = -BEAST/BHOR
      CB = -BNORTH/BHOR
      CTG = SIN(GDLAT*DTOR) 
      STG = COS(GDLAT*DTOR)
      ANG = (GDLON-ELON)*DTOR 
      SANG = SIN(ANG)
      CANG = COS(ANG)
      CTE = CTG*STD + STG*CTD*CB
      STE = SQRT(1. - CTE*CTE)
      SA = SB*CTD/STE
      CA = (STD*STG - CTD*CTG*CB)/STE
      CAPANG = CA*CANG - SA*SANG
      SAPANG = CA*SANG + SA*CANG
      STFCPA = STE*CTP*CAPANG - CTE*STP
      STFSPA = SAPANG*STE
      ALON = ATAN2(STFSPA,STFCPA)*RTOD
      R = ALT + RE
      HA = ALT + R*COTTD*COTTD
      A = 1. + HA/REQ
      RETURN
      end subroutine dipapx

       SUBROUTINE FINT (A1, A2, A3, A4, A5, A6, A7, RESULT)
C***PURPOSE  Second degree interpolation routine                                
C***REFER TO  FNDAPX                                                            
       RESULT = ((A2-A3)*(A7-A2)*(A7-A3)*A4-(A1-A3)*(A7-A1)*(A7-A3)*A5+         
     +  (A1-A2)*(A7-A1)*(A7-A2)*A6)/((A1-A2)*(A1-A3)*(A2-A3))
       RETURN                                                                   
       END subroutine fint 
!
! Original file ~bozo/pgms/apex/magfld.f copied 2/25/00.
!
      SUBROUTINE COFRM (DATE)
C          Assign DGRF/IGRF spherical harmonic coefficients, to degree and
C          order NMAX, for DATE, yyyy.fraction, into array G.  Coefficients
C          are interpolated from the DGRF dates through the current IGRF year.
C          Coefficients for a later DATE are extrapolated using the IGRF
C          initial value and the secular change coefficients.  A warning
C          message is issued if DATE is later than the last recommended
C          (5 yrs later than the IGRF).  An DATE input earlier than the
C          first DGRF (EPOCH(1)), results in a diagnostic and a STOP.
C
C          Output in COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
C             NMAX = Maximum order of spherical harmonic coefficients used
C             GB   = Coefficients for magnetic field calculation
C             GV   = Coefficients for magnetic potential calculation
C             ICHG = Flag indicating when GB,GV have been changed in COFRM
C
C          HISTORY (blame):
C          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
C          (formerly at SRI. Int., currently at Utah State).  Although set
C          up to accomodate second order time derivitives, the IGRF
C          (GTT, HTT) have been zero.  The spherical harmonic coefficients
C          degree and order is defined by NMAX (currently 10).
C
C          Jun 86:  Updated coefficients adding DGRF 1980 & IGRF 1985, which
C          were obtained from Eos Vol. 7, No. 24.  Common block MAG was
C          replaced by MAGCOF, thus removing variables not used in subroutine
C          FELDG.  (Roy Barnes)
C
C          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as described
C          in EOS Vol 73 Number 16 Apr 21 1992.  Other changes were made so
C          future updates should:
C            (1) Increment NDGY;
C            (2) Append to EPOCH the next IGRF year;
C            (3) Append the next DGRF coefficients to G1DIM and H1DIM; and
C            (4) Replace the IGRF initial values (G0, GT) and rates of
C                change indices (H0, HT).
C
C          Apr 94 (Art Richmond): Computation of GV added, for finding
C          magnetic potential.
C
C          Aug 95 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
C          obtained by anonymous ftp geomag.gsfc.nasa.gov (cd pub, mget table*)
C          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov),
C          but, problems are to be reported to baldwin@geomag.gsfc.nasa.gov
 
C          Oct 95 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see
C          email in folder).  Also found bug whereby coefficients were not being
C          updated in FELDG when IENTY did not change. ICHG was added to flag
C          date changes.  Also, a vestigial switch (IS) was removed from COFRM:
C          It was always 0 and involved 3 branch if statements in the main
C          polynomial construction loop (now numbered 200).
 
C          Feb 99 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid
C          possibility of compiler or loader options initializing memory
C          to something else (e.g., indefinite).  Also simplify the algebra
C          in COFRM; this does not effect results.

C          Mar 99 (Barnes):  Removed three branch if's from FELDG and changed
C          statement labels to ascending order

C          Jun 99 (Barnes):  Corrected RTOD definition in GD2CART.

C          May 00 (Barnes):  Replace IGRF 1995 with DGRF 1995, add IGRF
C          2000, and extend the earlier DGRF's backward to 1900.  A complete
C          set of coefficients came from a NGDC web page

C
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)
      DOUBLE PRECISION F,F0
      COMMON /MAGCOF/ G(144),GV(144)
      COMMON /MAGCOF_int/ NMAX,ICHG
!     DATA NMAX,ICHG /10,-99999/
 
!      PARAMETER (NDGY=6 , NYT = NDGY+1 , NGH = 144*NDGY)
      PARAMETER (NDGY=20 , NYT = NDGY+1 , NGH = 144*NDGY)
C          NDGY = Number of DGRF years of sets of coefficients
C          NYT  = Add one for the IGRF set (and point to it).
C          NGH  = Dimension of the equivalenced arrays
      DIMENSION GYR(12,12,NYT) , HYR(12,12,NYT), EPOCH(NYT) ,
     +          G1DIM(NGH)     , H1DIM(NGH) ,
     +          G0(12,12) , GT(12,12) , GTT(12,12) ,
     +          H0(12,12) , HT(12,12) ,  HTT(12,12)
      EQUIVALENCE (GYR(1,1,1),G1DIM(1))  , (HYR(1,1,1),H1DIM(1)) ,
     +            (GYR(1,1,NYT),G0(1,1)) , (HYR(1,1,NYT),H0(1,1))
 
      SAVE DATEL,GYR,HYR,EPOCH,G1DIM,H1DIM,G0,H0,GT,HT,GTT,HTT
C     SAVE DATEL,GYR,HYR,EPOCH,            G0,H0,GT,HT,GTT,HTT
      DATA DATEL /-999./
!      DATA EPOCH /1965. , 1970. , 1975. , 1980. , 1985. , 1990. , 1995./
      DATA EPOCH /1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940,
     +1945, 1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995,
     +2000/
C          D_/Dtime2 coefficients are 0
      DATA GTT/144*0./,HTT/144*0./
C          DGRF g(n,m) for 1965:
C          The "column" corresponds to "n" and
C          the "line" corresponds to "m" as indicated in column 6;
C          e.g., for 1965 g(0,3) = 1297. or g(6,6) = -111.
      DATA (G1DIM(I),I=1,144) /0,
     O  -31543,  -677,  1022,  876, -184,   63,  70,  11,   8, -3,  2*0,
     1   -2298,  2905, -1469,  628,  328,   61, -55,   8,  10, -4,  3*0,
     2            924,  1256,  660,  264,  -11,   0,  -4,   1,  2,  4*0,
     3                   572, -361,    5, -217,  34,  -9, -11, -5,  5*0,
     4                         134,  -86,  -58, -41,   1,  12, -2,  6*0,
     5                               -16,   59, -21,   2,   1,  6,  7*0,
     6                                     -90,  18,  -9,  -2,  4,  8*0,
     7                                            6,   5,   2,  0,  9*0,
     8                                                 8,  -1,  2, 10*0,
     9                                                     -1,  2, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1905:
      DATA (G1DIM(I),I=145,288) /0,
     O  -31464,  -728,  1037,  880, -192,   62,  70,  11,   8, -3,  2*0,
     1   -2298,  2928, -1494,  643,  328,   60, -54,   8,  10, -4,  3*0,
     2           1041,  1239,  653,  259,  -11,   0,  -4,   1,  2,  4*0,
     3                   635, -380,   -1, -221,  33,  -9, -11, -5,  5*0,
     4                         146,  -93,  -57, -41,   1,  12, -2,  6*0,
     5                               -26,   57, -20,   2,   1,  6,  7*0,
     6                                     -92,  18,  -8,  -2,  4,  8*0,
     7                                            6,   5,   2,  0,  9*0,
     8                                                 8,   0,  2, 10*0,
     9                                                     -1,  2, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1910:
      DATA (G1DIM(I),I=289,432) /0,
     O  -31354,  -769,  1058,  884, -201,   62,  71,  11,   8, -3,  2*0,
     1   -2297,  2948, -1524,  660,  327,   58, -54,   8,  10, -4,  3*0,
     2           1176,  1223,  644,  253,  -11,   1,  -4,   1,  2,  4*0,
     3                   705, -400,   -9, -224,  32,  -9, -11, -5,  5*0,
     4                         160, -102,  -54, -40,   1,  12, -2,  6*0,
     5                               -38,   54, -19,   2,   1,  6,  7*0,
     6                                     -95,  18,  -8,  -2,  4,  8*0,
     7                                            6,   5,   2,  0,  9*0,
     8                                                 8,   0,  2, 10*0,
     9                                                     -1,  2, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1915:
      DATA (G1DIM(I),I=433,576) /0,
     O  -31212,  -802,  1084,  887, -211,   61,  72,  11,   8, -3,  2*0,
     1   -2306,  2956, -1559,  678,  327,   57, -54,   8,  10, -4,  3*0,
     2           1309,  1212,  631,  245,  -10,   2,  -4,   1,  2,  4*0,
     3                   778, -416,  -16, -228,  31,  -9, -11, -5,  5*0,
     4                         178, -111,  -51, -38,   2,  12, -2,  6*0,
     5                               -51,   49, -18,   3,   1,  6,  7*0,
     6                                     -98,  19,  -8,  -2,  4,  8*0,
     7                                            6,   6,   2,  0,  9*0,
     8                                                 8,   0,  1, 10*0,
     9                                                     -1,  2, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1920:
      DATA (G1DIM(I),I=577,720) /0,
     O  -31060,  -839,  1111,  889, -221,   61,  73,  11,   8, -3,  2*0,
     1   -2317,  2959, -1600,  695,  326,   55, -54,   7,  10, -4,  3*0,
     2           1407,  1205,  616,  236,  -10,   2,  -3,   1,  2,  4*0,
     3                   839, -424,  -23, -233,  29,  -9, -11, -5,  5*0,
     4                         199, -119,  -46, -37,   2,  12, -2,  6*0,
     5                               -62,   44, -16,   4,   1,  6,  7*0,
     6                                    -101,  19,  -7,  -2,  4,  8*0,
     7                                            6,   6,   2,  0,  9*0,
     8                                                 8,   0,  1, 10*0,
     9                                                     -1,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1925:
      DATA (G1DIM(I),I=721,864) /0,
     O  -30926,  -893,  1140,  891, -230,   61,  73,  11,   8, -3,  2*0,
     1   -2318,  2969, -1645,  711,  326,   54, -54,   7,  10, -4,  3*0,
     2           1471,  1202,  601,  226,   -9,   3,  -3,   1,  2,  4*0,
     3                   881, -426,  -28, -238,  27,  -9, -11, -5,  5*0,
     4                         217, -125,  -40, -35,   2,  12, -2,  6*0,
     5                               -69,   39, -14,   4,   1,  6,  7*0,
     6                                    -103,  19,  -7,  -2,  4,  8*0,
     7                                            6,   7,   2,  0,  9*0,
     8                                                 8,   0,  1, 10*0,
     9                                                     -1,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1930:
      DATA (G1DIM(I),I=865,1008) /0,
     O  -30805,  -951,  1172,  896, -237,   60,  74,  11,   8, -3,  2*0,
     1   -2316,  2980, -1692,  727,  327,   53, -54,   7,  10, -4,  3*0,
     2           1517,  1205,  584,  218,   -9,   4,  -3,   1,  2,  4*0,
     3                   907, -422,  -32, -242,  25,  -9, -12, -5,  5*0,
     4                         234, -131,  -32, -34,   2,  12, -2,  6*0,
     5                               -74,   32, -12,   5,   1,  6,  7*0,
     6                                    -104,  18,  -6,  -2,  4,  8*0,
     7                                            6,   8,   3,  0,  9*0,
     8                                                 8,   0,  1, 10*0,
     9                                                     -2,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1935:
      DATA (G1DIM(I),I=1009,1152) /0,
     O  -30715, -1018,  1206,  903, -241,   59,  74,  11,   8, -3,  2*0,
     1   -2306,  2984, -1740,  744,  329,   53, -53,   7,  10, -4,  3*0,
     2           1550,  1215,  565,  211,   -8,   4,  -3,   1,  2,  4*0,
     3                   918, -415,  -33, -246,  23,  -9, -12, -5,  5*0,
     4                         249, -136,  -25, -33,   1,  11, -2,  6*0,
     5                               -76,   25, -11,   6,   1,  6,  7*0,
     6                                    -106,  18,  -6,  -2,  4,  8*0,
     7                                            6,   8,   3,  0,  9*0,
     8                                                 7,   0,  2, 10*0,
     9                                                     -2,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1940:
      DATA (G1DIM(I),I=1153,1296) /0,
     O  -30654, -1106,  1240,  914, -241,   57,  74,  11,   8, -3,  2*0,
     1   -2292,  2981, -1790,  762,  334,   54, -53,   7,  10, -4,  3*0,
     2           1566,  1232,  550,  208,   -7,   4,  -3,   1,  2,  4*0,
     3                   916, -405,  -33, -249,  20, -10, -12, -5,  5*0,
     4                         265, -141,  -18, -31,   1,  11, -2,  6*0,
     5                               -76,   18,  -9,   6,   1,  6,  7*0,
     6                                    -107,  17,  -5,  -2,  4,  8*0,
     7                                            5,   9,   3,  0,  9*0,
     8                                                 7,   1,  2, 10*0,
     9                                                     -2,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1945:
      DATA (G1DIM(I),I=1297,1440) /0,
     O  -30594, -1244,  1282,  944, -253,   59,  70,  13,   5, -3,  2*0,
     1   -2285,  2990, -1834,  776,  346,   57, -40,   7, -21, 11,  3*0,
     2           1578,  1255,  544,  194,    6,   0,  -8,   1,  1,  4*0,
     3                   913, -421,  -20, -246,   0,  -5, -11,  2,  5*0,
     4                         304, -142,  -25, -29,   9,   3, -5,  6*0,
     5                               -82,   21, -10,   7,  16, -1,  7*0,
     6                                    -104,  15, -10,  -3,  8,  8*0,
     7                                           29,   7,  -4, -1,  9*0,
     8                                                 2,  -3, -3, 10*0,
     9                                                     -4,  5, 11*0,
     A                                                         -2, 13*0/
C          DGRF g(n,m) for 1950:
      DATA (G1DIM(I),I=1441,1584) /0,
     O  -30554, -1341,  1297,  954, -240,   54,  65,  22,   3, -8,  2*0,
     1   -2250,  2998, -1889,  792,  349,   57, -55,  15,  -7,  4,  3*0,
     2           1576,  1274,  528,  211,    4,   2,  -4,  -1, -1,  4*0,
     3                   896, -408,  -20, -247,   1,  -1, -25, 13,  5*0,
     4                         303, -147,  -16, -40,  11,  10, -4,  6*0,
     5                               -76,   12,  -7,  15,   5,  4,  7*0,
     6                                    -105,   5, -13,  -5, 12,  8*0,
     7                                           19,   5,  -2,  3,  9*0,
     8                                                -1,   3,  2, 10*0,
     9                                                      8, 10, 11*0,
     A                                                          3, 13*0/
C          DGRF g(n,m) for 1955:
      DATA (G1DIM(I),I=1585,1728) /0,
     O  -30500, -1440,  1302,  958, -229,   47,  65,  11,   4, -3,  2*0,
     1   -2215,  3003, -1944,  796,  360,   57, -56,   9,   9, -5,  3*0,
     2           1581,  1288,  510,  230,    3,   2,  -6,  -4, -1,  4*0,
     3                   882, -397,  -23, -247,  10, -14,  -5,  2,  5*0,
     4                         290, -152,   -8, -32,   6,   2, -3,  6*0,
     5                               -69,    7, -11,  10,   4,  7,  7*0,
     6                                    -107,   9,  -7,   1,  4,  8*0,
     7                                           18,   6,   2, -2,  9*0,
     8                                                 9,   2,  6, 10*0,
     9                                                      5, -2, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1960:
      DATA (G1DIM(I),I=1729,1872) /0,
     O  -30421, -1555,  1302,  957, -222,   46,  67,  15,   4,  1,  2*0,
     1   -2169,  3002, -1992,  800,  362,   58, -56,   6,   6, -3,  3*0,
     2           1590,  1289,  504,  242,    1,   5,  -4,   0,  4,  4*0,
     3                   878, -394,  -26, -237,  15, -11,  -9,  0,  5*0,
     4                         269, -156,   -1, -32,   2,   1, -1,  6*0,
     5                               -63,   -2,  -7,  10,   4,  4,  7*0,
     6                                    -113,  17,  -5,  -1,  6,  8*0,
     7                                            8,  10,  -2,  1,  9*0,
     8                                                 8,   3, -1, 10*0,
     9                                                     -1,  2, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1965:
      DATA (G1DIM(I),I=1873,2016) /0,
     O  -30334, -1662,  1297,  957, -219,   45,  75,  13,   8, -2,  2*0,
     1   -2119,  2997, -2038,  804,  358,   61, -57,   5,  10, -3,  3*0,
     2           1594,  1292,  479,  254,    8,   4,  -4,   2,  2,  4*0,
     3                   856, -390,  -31, -228,  13, -14, -13, -5,  5*0,
     4                         252, -157,    4, -26,   0,  10, -2,  6*0,
     5                               -62,    1,  -6,   8,  -1,  4,  7*0,
     6                                    -111,  13,  -1,  -1,  4,  8*0,
     7                                            1,  11,   5,  0,  9*0,
     8                                                 4,   1,  2, 10*0,
     9                                                     -2,  2, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1970:
      DATA (G1DIM(I),I=2017,2160) /0,
     O  -30220, -1781,  1287,  952, -216,   43,  72,  14,   8, -3,  2*0,
     1   -2068,  3000, -2091,  800,  359,   64, -57,   6,  10, -3,  3*0,
     2           1611,  1278,  461,  262,   15,   1,  -2,   2,  2,  4*0,
     3                   838, -395,  -42, -212,  14, -13, -12, -5,  5*0,
     4                         234, -160,    2, -22,  -3,  10, -1,  6*0,
     5                               -56,    3,  -2,   5,  -1,  6,  7*0,
     6                                    -112,  13,   0,   0,  4,  8*0,
     7                                           -2,  11,   3,  1,  9*0,
     8                                                 3,   1,  0, 10*0,
     9                                                     -1,  3, 11*0,
     A                                                         -1, 13*0/
C          DGRF g(n,m) for 1975:
      DATA (G1DIM(I),I=2161,2304) /0,
     O  -30100, -1902,  1276,  946, -218,   45,  71,  14,   7, -3,  2*0,
     1   -2013,  3010, -2144,  791,  356,   66, -56,   6,  10, -3,  3*0,
     2           1632,  1260,  438,  264,   28,   1,  -1,   2,  2,  4*0,
     3                   830, -405,  -59, -198,  16, -12, -12, -5,  5*0,
     4                         216, -159,    1, -14,  -8,  10, -2,  6*0,
     5                               -49,    6,   0,   4,  -1,  5,  7*0,
     6                                    -111,  12,   0,  -1,  4,  8*0,
     7                                           -5,  10,   4,  1,  9*0,
     8                                                 1,   1,  0, 10*0,
     9                                                     -2,  3, 11*0,
     A                                                         -1, 13*0/
C          DGRF g(n,m) for 1980:
      DATA (G1DIM(I),I=2305,2448) /0,
     O  -29992, -1997,  1281,  938, -218,   48,  72,  18,   5, -4,  2*0,
     1   -1956,  3027, -2180,  782,  357,   66, -59,   6,  10, -4,  3*0,
     2           1663,  1251,  398,  261,   42,   2,   0,   1,  2,  4*0,
     3                   833, -419,  -74, -192,  21, -11, -12, -5,  5*0,
     4                         199, -162,    4, -12,  -7,   9, -2,  6*0,
     5                               -48,   14,   1,   4,  -3,  5,  7*0,
     6                                    -108,  11,   3,  -1,  3,  8*0,
     7                                           -2,   6,   7,  1,  9*0,
     8                                                -1,   2,  2, 10*0,
     9                                                     -5,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1985:
      DATA (G1DIM(I),I=2449,2592) /0,
     O  -29873, -2072,  1296,  936, -214,   53,  74,  21,   5, -4,  2*0,
     1   -1905,  3044, -2208,  780,  355,   65, -62,   6,  10, -4,  3*0,
     2           1687,  1247,  361,  253,   51,   3,   0,   1,  3,  4*0,
     3                   829, -424,  -93, -185,  24, -11, -12, -5,  5*0,
     4                         170, -164,    4,  -6,  -9,   9, -2,  6*0,
     5                               -46,   16,   4,   4,  -3,  5,  7*0,
     6                                    -102,  10,   4,  -1,  3,  8*0,
     7                                            0,   4,   7,  1,  9*0,
     8                                                -4,   1,  2, 10*0,
     9                                                     -5,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1990:
      DATA (G1DIM(I),I=2593,2736) /0,
     O  -29775, -2131,  1314,  939, -214,   61,  77,  23,   4, -3,  2*0,
     1   -1848,  3059, -2239,  780,  353,   65, -64,   5,   9, -4,  3*0,
     2           1686,  1248,  325,  245,   59,   2,  -1,   1,  2,  4*0,
     3                   802, -423, -109, -178,  26, -10, -12, -5,  5*0,
     4                         141, -165,    3,  -1, -12,   9, -2,  6*0,
     5                               -36,   18,   5,   3,  -4,  4,  7*0,
     6                                     -96,   9,   4,  -2,  3,  8*0,
     7                                            0,   2,   7,  1,  9*0,
     8                                                -6,   1,  3, 10*0,
     9                                                     -6,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF g(n,m) for 1995:
      DATA (G1DIM(I),I=2737,2880) /0,
     O  -29682, -2197,  1329,  941, -210,   66,  78,  24,   4, -3,  2*0,
     1   -1789,  3074, -2268,  782,  352,   64, -67,   4,   9, -4,  3*0,
     2           1685,  1249,  291,  237,   65,   1,  -1,   1,  2,  4*0,
     3                   769, -421, -122, -172,  29,  -9, -12, -5,  5*0,
     4                         116, -167,    2,   4, -14,   9, -2,  6*0,
     5                               -26,   17,   8,   4,  -4,  4,  7*0,
     6                                     -94,  10,   5,  -2,  3,  8*0,
     7                                           -2,   0,   7,  1,  9*0,
     8                                                -7,   0,  3, 10*0,
     9                                                     -6,  3, 11*0,
     A                                                          0, 13*0/
C          DGRF h(n,m) for 1900:
      DATA (H1DIM(I),I=1,144) /13*0,
     1    5922, -1061,  -330,  195, -210,   -9, -45,   8, -20,  2,  3*0,
     2           1121,     3,  -69,   53,   83, -13, -14,  14,  1,  4*0,
     3                   523, -210,  -33,    2, -10,   7,   5,  2,  5*0,
     4                         -75, -124,  -35,  -1, -13,  -3,  6,  6*0,
     5                                 3,   36,  28,   5,  -2, -4,  7*0,
     6                                     -69, -12,  16,   8,  0,  8*0,
     7                                          -22,  -5,  10, -2,  9*0,
     8                                               -18,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1905:
      DATA (H1DIM(I),I=145,288) /13*0,
     1    5909, -1086,  -357,  203, -193,   -7, -46,   8, -20,  2,  3*0,
     2           1065,    34,  -77,   56,   86, -14, -15,  14,  1,  4*0,
     3                   480, -201,  -32,    4, -11,   7,   5,  2,  5*0,
     4                         -65, -125,  -32,   0, -13,  -3,  6,  6*0,
     5                                11,   32,  28,   5,  -2, -4,  7*0,
     6                                     -67, -12,  16,   8,  0,  8*0,
     7                                          -22,  -5,  10, -2,  9*0,
     8                                               -18,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1910:
      DATA (H1DIM(I),I=289,432) /13*0,
     1    5898, -1128,  -389,  211, -172,   -5, -47,   8, -20,  2,  3*0,
     2           1000,    62,  -90,   57,   89, -14, -15,  14,  1,  4*0,
     3                   425, -189,  -33,    5, -12,   6,   5,  2,  5*0,
     4                         -55, -126,  -29,   1, -13,  -3,  6,  6*0,
     5                                21,   28,  28,   5,  -2, -4,  7*0,
     6                                     -65, -13,  16,   8,  0,  8*0,
     7                                          -22,  -5,  10, -2,  9*0,
     8                                               -18,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1915:
      DATA (H1DIM(I),I=433,576) /13*0,
     1    5875, -1191,  -421,  218, -148,   -2, -48,   8, -20,  2,  3*0,
     2            917,    84, -109,   58,   93, -14, -15,  14,  1,  4*0,
     3                   360, -173,  -34,    8, -12,   6,   5,  2,  5*0,
     4                         -51, -126,  -26,   2, -13,  -3,  6,  6*0,
     5                                32,   23,  28,   5,  -2, -4,  7*0,
     6                                     -62, -15,  16,   8,  0,  8*0,
     7                                          -22,  -5,  10, -2,  9*0,
     8                                               -18,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1920:
      DATA (H1DIM(I),I=577,720) /13*0,
     1    5845, -1259,  -445,  220, -122,    0, -49,   8, -20,  2,  3*0,
     2            823,   103, -134,   58,   96, -14, -15,  14,  1,  4*0,
     3                   293, -153,  -38,   11, -13,   6,   5,  2,  5*0,
     4                         -57, -125,  -22,   4, -14,  -3,  6,  6*0,
     5                                43,   18,  28,   5,  -2, -4,  7*0,
     6                                     -57, -16,  17,   9,  0,  8*0,
     7                                          -22,  -5,  10, -2,  9*0,
     8                                               -19,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1925:
      DATA (H1DIM(I),I=721,864) /13*0,
     1    5817, -1334,  -462,  216,  -96,    3, -50,   8, -20,  2,  3*0,
     2            728,   119, -163,   58,   99, -14, -15,  14,  1,  4*0,
     3                   229, -130,  -44,   14, -14,   6,   5,  2,  5*0,
     4                         -70, -122,  -18,   5, -14,  -3,  6,  6*0,
     5                                51,   13,  29,   5,  -2, -4,  7*0,
     6                                     -52, -17,  17,   9,  0,  8*0,
     7                                          -21,  -5,  10, -2,  9*0,
     8                                               -19,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1930:
      DATA (H1DIM(I),I=865,1008) /13*0,
     1    5808, -1424,  -480,  205,  -72,    4, -51,   8, -20,  2,  3*0,
     2            644,   133, -195,   60,  102, -15, -15,  14,  1,  4*0,
     3                   166, -109,  -53,   19, -14,   5,   5,  2,  5*0,
     4                         -90, -118,  -16,   6, -14,  -3,  6,  6*0,
     5                                58,    8,  29,   5,  -2, -4,  7*0,
     6                                     -46, -18,  18,   9,  0,  8*0,
     7                                          -20,  -5,  10, -2,  9*0,
     8                                               -19,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1935:
      DATA (H1DIM(I),I=1009,1152) /13*0,
     1    5812, -1520,  -494,  188,  -51,    4, -52,   8, -20,  2,  3*0,
     2            586,   146, -226,   64,  104, -17, -15,  15,  1,  4*0,
     3                   101,  -90,  -64,   25, -14,   5,   5,  2,  5*0,
     4                        -114, -115,  -15,   7, -15,  -3,  6,  6*0,
     5                                64,    4,  29,   5,  -3, -4,  7*0,
     6                                     -40, -19,  18,   9,  0,  8*0,
     7                                          -19,  -5,  11, -1,  9*0,
     8                                               -19,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1940:
      DATA (H1DIM(I),I=1153,1296) /13*0,
     1    5821, -1614,  -499,  169,  -33,    4, -52,   8, -21,  2,  3*0,
     2            528,   163, -252,   71,  105, -18, -14,  15,  1,  4*0,
     3                    43,  -72,  -75,   33, -14,   5,   5,  2,  5*0,
     4                        -141, -113,  -15,   7, -15,  -3,  6,  6*0,
     5                                69,    0,  29,   5,  -3, -4,  7*0,
     6                                     -33, -20,  19,   9,  0,  8*0,
     7                                          -19,  -5,  11, -1,  9*0,
     8                                               -19,  -2,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1945:
      DATA (H1DIM(I),I=1297,1440) /13*0,
     1    5810, -1702,  -499,  144,  -12,    6, -45,  12, -27,  5,  3*0,
     2            477,   186, -276,   95,  100, -18, -21,  17,  1,  4*0,
     3                   -11,  -55,  -67,   16,   2, -12,  29,-20,  5*0,
     4                        -178, -119,   -9,   6,  -7,  -9, -1,  6*0,
     5                                82,  -16,  28,   2,   4, -6,  7*0,
     6                                     -39, -17,  18,   9,  6,  8*0,
     7                                          -22,   3,   6, -4,  9*0,
     8                                               -11,   1, -2, 10*0,
     9                                                      8,  0, 11*0,
     A                                                         -2, 13*0/
C          DGRF h(n,m) for 1950:
      DATA (H1DIM(I),I=1441,1584) /13*0,
     1    5815, -1810,  -476,  136,    3,   -1, -35,   5, -24, 13,  3*0,
     2            381,   206, -278,  103,   99, -17, -22,  19, -2,  4*0,
     3                   -46,  -37,  -87,   33,   0,   0,  12,-10,  5*0,
     4                        -210, -122,  -12,  10, -21,   2,  2,  6*0,
     5                                80,  -12,  36,  -8,   2, -3,  7*0,
     6                                     -30, -18,  17,   8,  6,  8*0,
     7                                          -16,  -4,   8, -3,  9*0,
     8                                               -17, -11,  6, 10*0,
     9                                                     -7, 11, 11*0,
     A                                                          8, 13*0/
C          DGRF h(n,m) for 1955:
      DATA (H1DIM(I),I=1585,1728) /13*0,
     1    5820, -1898,  -462,  133,   15,   -9, -50,  10, -11, -4,  3*0,
     2            291,   216, -274,  110,   96, -24, -15,  12,  0,  4*0,
     3                   -83,  -23,  -98,   48,  -4,   5,   7, -8,  5*0,
     4                        -230, -121,  -16,   8, -23,   6, -2,  6*0,
     5                                78,  -12,  28,   3,  -2, -4,  7*0,
     6                                     -24, -20,  23,  10,  1,  8*0,
     7                                          -18,  -4,   7, -3,  9*0,
     8                                               -13,  -6,  7, 10*0,
     9                                                      5, -1, 11*0,
     A                                                         -3, 13*0/
C          DGRF h(n,m) for 1960:
      DATA (H1DIM(I),I=1729,1872) /13*0,
     1    5791, -1967,  -414,  135,   16,  -10, -55,  11, -18,  4,  3*0,
     2            206,   224, -278,  125,   99, -28, -14,  12,  1,  4*0,
     3                  -130,    3, -117,   60,  -6,   7,   2,  0,  5*0,
     4                        -255, -114,  -20,   7, -18,   0,  2,  6*0,
     5                                81,  -11,  23,   4,  -3, -5,  7*0,
     6                                     -17, -18,  23,   9,  1,  8*0,
     7                                          -17,   1,   8, -1,  9*0,
     8                                               -20,   0,  6, 10*0,
     9                                                      5,  0, 11*0,
     A                                                         -7, 13*0/
C          DGRF h(n,m) for 1965:
      DATA (H1DIM(I),I=1873,2016) /13*0,
     1    5776, -2016,  -404,  148,   19,  -11, -61,   7, -22,  2,  3*0,
     2            114,   240, -269,  128,  100, -27, -12,  15,  1,  4*0,
     3                  -165,   13, -126,   68,  -2,   9,   7,  2,  5*0,
     4                        -269,  -97,  -32,   6, -16,  -4,  6,  6*0,
     5                                81,   -8,  26,   4,  -5, -4,  7*0,
     6                                      -7, -23,  24,  10,  0,  8*0,
     7                                          -12,  -3,  10, -2,  9*0,
     8                                               -17,  -4,  3, 10*0,
     9                                                      1,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1970:
      DATA (H1DIM(I),I=2017,2160) /13*0,
     1    5737, -2047,  -366,  167,   26,  -12, -70,   7, -21,  1,  3*0,
     2             25,   251, -266,  139,  100, -27, -15,  16,  1,  4*0,
     3                  -196,   26, -139,   72,  -4,   6,   6,  3,  5*0,
     4                        -279,  -91,  -37,   8, -17,  -4,  4,  6*0,
     5                                83,   -6,  23,   6,  -5, -4,  7*0,
     6                                       1, -23,  21,  10,  0,  8*0,
     7                                          -11,  -6,  11, -1,  9*0,
     8                                               -16,  -2,  3, 10*0,
     9                                                      1,  1, 11*0,
     A                                                         -4, 13*0/
C          DGRF h(n,m) for 1975:
      DATA (H1DIM(I),I=2161,2304) /13*0,
     1    5675, -2067,  -333,  191,   31,  -13, -77,   6, -21,  1,  3*0,
     2            -68,   262, -265,  148,   99, -26, -16,  16,  1,  4*0,
     3                  -223,   39, -152,   75,  -5,   4,   7,  3,  5*0,
     4                        -288,  -83,  -41,  10, -19,  -4,  4,  6*0,
     5                                88,   -4,  22,   6,  -5, -4,  7*0,
     6                                      11, -23,  18,  10, -1,  8*0,
     7                                          -12, -10,  11, -1,  9*0,
     8                                               -17,  -3,  3, 10*0,
     9                                                      1,  1, 11*0,
     A                                                         -5, 13*0/
C          DGRF h(n,m) for 1980:
      DATA (H1DIM(I),I=2305,2448) /13*0,
     1    5604, -2129,  -336,  212,   46,  -15, -82,   7, -21,  1,  3*0,
     2           -200,   271, -257,  150,   93, -27, -18,  16,  0,  4*0,
     3                  -252,   53, -151,   71,  -5,   4,   9,  3,  5*0,
     4                        -297,  -78,  -43,  16, -22,  -5,  6,  6*0,
     5                                92,   -2,  18,   9,  -6, -4,  7*0,
     6                                      17, -23,  16,   9,  0,  8*0,
     7                                          -10, -13,  10, -1,  9*0,
     8                                               -15,  -6,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1985:
      DATA (H1DIM(I),I=2449,2592) /13*0,
     1    5500, -2197,  -310,  232,   47,  -16, -83,   8, -21,  1,  3*0,
     2           -306,   284, -249,  150,   88, -27, -19,  15,  0,  4*0,
     3                  -297,   69, -154,   69,  -2,   5,   9,  3,  5*0,
     4                        -297,  -75,  -48,  20, -23,  -6,  6,  6*0,
     5                                95,   -1,  17,  11,  -6, -4,  7*0,
     6                                      21, -23,  14,   9,  0,  8*0,
     7                                           -7, -15,   9, -1,  9*0,
     8                                               -11,  -7,  4, 10*0,
     9                                                      2,  0, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1990:
      DATA (H1DIM(I),I=2593,2736) /13*0,
     1    5406, -2279,  -284,  247,   46,  -16, -80,  10, -20,  2,  3*0,
     2           -373,   293, -240,  154,   82, -26, -19,  15,  1,  4*0,
     3                  -352,   84, -153,   69,   0,   6,  11,  3,  5*0,
     4                        -299,  -69,  -52,  21, -22,  -7,  6,  6*0,
     5                                97,    1,  17,  12,  -7, -4,  7*0,
     6                                      24, -23,  12,   9,  0,  8*0,
     7                                           -4, -16,   8, -2,  9*0,
     8                                               -10,  -7,  3, 10*0,
     9                                                      2, -1, 11*0,
     A                                                         -6, 13*0/
C          DGRF h(n,m) for 1995:
      DATA (H1DIM(I),I=2737,2880) /13*0,
     1    5318, -2356,  -263,  262,   44,  -16, -77,  12, -19,  2,  3*0,
     2           -425,   302, -232,  157,   77, -25, -20,  15,  1,  4*0,
     3                  -406,   98, -152,   67,   3,   7,  11,  3,  5*0,
     4                        -301,  -64,  -57,  22, -21,  -7,  6,  6*0,
     5                                99,    4,  16,  12,  -7, -4,  7*0,
     6                                      28, -23,  10,   9,  0,  8*0,
     7                                           -3, -17,   7, -2,  9*0,
     8                                               -10,  -8,  3, 10*0,
     9                                                      1, -1, 11*0,
     A                                                         -6, 13*0/
C          Initial coefficients g0 (IGRF for 2000):
      DATA G0 /0,
     O  -29615, -2267,  1341,  935, -217,   72,  79,  25,   5, -2,  2*0,
     1   -1728,  3072, -2290,  787,  351,   68, -74,   6,   9, -6,  3*0,
     2           1672,  1253,  251,  222,   74,   0,  -9,   3,  2,  4*0,
     3                   715, -405, -131, -161,  33,  -8,  -8, -3,  5*0,
     4                         110, -169,   -5,   9, -17,   6,  0,  6*0,
     5                               -12,   17,   7,   9,  -9,  4,  7*0,
     6                                     -91,   8,   7,  -2,  1,  8*0,
     7                                           -2,  -8,   9,  2,  9*0,
     8                                                -7,  -4,  4, 10*0,
     9                                                     -8,  0, 11*0,
     A                                                         -1, 13*0/
C          D_/Dtime coefficients gt (IGRF for 2000-2005):
      DATA GT /0,
     O    14.6, -12.4,   0.7, -1.3,  0.0,  1.0,-0.4,-0.3, 0.0,0.0,  2*0,
     1    10.7,   1.1,  -5.4,  1.6, -0.7, -0.4,-0.4, 0.2, 0.0,0.0,  3*0,
     2           -1.1,   0.9, -7.3, -2.1,  0.9,-0.3,-0.3, 0.0,0.0,  4*0,
     3                  -7.7,  2.9, -2.8,  2.0, 1.1, 0.4, 0.0,0.0,  5*0,
     4                        -3.2, -0.8, -0.6, 1.1,-1.0, 0.0,0.0,  6*0,
     5                               2.5, -0.3,-0.2, 0.3, 0.0,0.0,  7*0,
     6                                     1.2, 0.6,-0.5, 0.0,0.0,  8*0,
     7                                         -0.9,-0.7, 0.0,0.0,  9*0,
     8                                              -0.4, 0.0,0.0, 10*0,
     9                                                    0.0,0.0, 11*0,
     A                                                        0.0, 13*0/
C          Initial coefficients h0 (IGRF for 2000):
      DATA H0 /13*0,
     1    5186, -2478,  -227,  272,   44,  -17, -65,  12, -20,  1,  3*0,
     2           -458,   296, -232,  172,   64, -24, -22,  13,  0,  4*0,
     3                  -492,  119, -134,   65,   6,   8,  12,  4,  5*0,
     4                        -304,  -40,  -61,  24, -21,  -6,  5,  6*0,
     5                               107,    1,  15,  15,  -8, -6,  7*0,
     6                                      44, -25,   9,   9, -1,  8*0,
     7                                           -6, -16,   4, -3,  9*0,
     8                                                -3,  -8,  0, 10*0,
     9                                                      5, -2, 11*0,
     A                                                         -8, 13*0/
C          D_/Dtime coefficients ht (IGRF for 2000-2005):
      DATA HT /13*0,
     1   -22.5, -20.6,   6.0,  2.1, -0.1, -0.2, 1.1, 0.1, 0.0,0.0,  3*0,
     2           -9.6,  -0.1,  1.3,  0.6, -1.4, 0.0, 0.0, 0.0,0.0,  4*0,
     3                 -14.2,  5.0,  1.7,  0.0, 0.3, 0.0, 0.0,0.0,  5*0,
     4                         0.3,  1.9, -0.8,-0.1, 0.3, 0.0,0.0,  6*0,
     5                               0.1,  0.0,-0.6, 0.6, 0.0,0.0,  7*0,
     6                                     0.9,-0.7,-0.4, 0.0,0.0,  8*0,
     7                                          0.2, 0.3, 0.0,0.0,  9*0,
     8                                               0.7, 0.0,0.0, 10*0,
     9                                                    0.0,0.0, 11*0,
     A                                                        0.0, 13*0/
 
C          Do not need to load new coefficients if date has not changed
      ICHG = 0
      IF (DATE .EQ. DATEL) GO TO 300
      DATEL = DATE
      ICHG = 1
 
C          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NYT)+5.) WRITE(6,9200) DATE
 
      DO 100 I=1,NYT
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE
 
      TIME = DATE
      T = TIME-EPOCH(IY)
      G(1)  = 0.0
      GV(1) = 0.0
      I = 2
      F0 = -1.0D-5
      DO 200 N=1,NMAX
      F0 = F0 * REAL(N)/2.
      F  = F0 / SQRT(2.0)
      NN = N+1
      MM = 1
      IF (IY .EQ. NYT) THEN
C          Extrapolate coefficients
        G(I) = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F0
      ELSE
C          Interpolate coefficients
        G(I) = (GYR(NN,MM,IY) +
     +         T/5.0 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F0
      ENDIF
      GV(I) = G(I) / REAL(NN)
      I = I+1
      DO 200 M=1,N
      F = F / SQRT( REAL(N-M+1) / REAL(N+M) )
      NN = N+1
      MM = M+1
      I1 = I+1
      IF (IY .EQ. NYT) THEN
C          Extrapolate coefficients
        G(I)  = ((GTT(NN,MM)*T + GT(NN,MM))*T + G0(NN,MM)) * F
        G(I1) = ((HTT(NN,MM)*T + HT(NN,MM))*T + H0(NN,MM)) * F
      ELSE
C          Interpolate coefficients
        G(I)  = (GYR(NN,MM,IY) +
     +          T/5.0 * (GYR(NN,MM,IY+1)-GYR(NN,MM,IY))) * F
        G(I1) = (HYR(NN,MM,IY) +
     +          T/5.0 * (HYR(NN,MM,IY+1)-HYR(NN,MM,IY))) * F
      ENDIF
      RNN = REAL(NN)
      GV(I)  = G(I)  / RNN
      GV(I1) = G(I1) / RNN
  200 I = I+2
 
  300 CONTINUE

      RETURN
 
C          Error trap diagnostics:
 9100 WRITE(6,'('' '',/,
     +'' COFRM:  DATE'',F9.3,'' preceeds DGRF coefficients'',
     +                          '' presently coded.'')') DATE
      call shutdown('mor cod')
 9200 FORMAT(' ',/,
     +' COFRM:  DATE',F9.3,' is after the maximum',
     +                          ' recommended for extrapolation.')
      END
 
      SUBROUTINE DYPOL (COLAT,ELON,VP)
C          Computes parameters for dipole component of geomagnetic field.
C          COFRM must be called before calling DYPOL!
C          940504 A. D. Richmond
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          RETURNS:
C            COLAT = Geocentric colatitude of geomagnetic dipole north pole
C                    (deg)
C            ELON  = East longitude of geomagnetic dipole north pole (deg)
C            VP    = Magnitude, in T.m, of dipole component of magnetic
C                    potential at geomagnetic pole and geocentric radius
C                    of 6371.2 km
 
      PARAMETER (RTOD=5.72957795130823E1 , DTOR=1.745329251994330E-2,
     +           RE=6371.2, REQ=6378.160)
      COMMON /MAGCOF/ G(144),GV(144)
      COMMON /MAGCOF_int/ NMAX,ICHG
 
C          Compute geographic colatitude and longitude of the north pole of
C          earth centered dipole
      GPL = SQRT( G(2  )**2+ G(3  )**2+ G(4  )**2)
      CTP = G(2  )/GPL
      STP = SQRT(1. - CTP*CTP)
      COLAT = (ACOS(CTP))*RTOD
      ELON = ATAN2( G(4  ), G(3  ))*RTOD
 
C          Compute magnitude of magnetic potential at pole, radius Re.
      VP = .2*GPL*RE
C          .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through F0 in COFRM).
 
      RETURN
      END
 
      SUBROUTINE FELDG (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
C          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
C          COFRM must be called to establish coefficients for correct date
C          prior to calling FELDG.
C
C          IENTY is an input flag controlling the meaning and direction of the
C                remaining formal arguments:
C          IENTY = 1
C            INPUTS:
C              GLAT = Latitude of point (deg)
C              GLON = Longitude (east=+) of point (deg)
C              ALT  = Ht of point (km)
C            RETURNS:
C              BNRTH  north component of field vector (Gauss)
C              BEAST  east component of field vector  (Gauss)
C              BDOWN  downward component of field vector (Gauss)
C              BABS   magnitude of field vector (Gauss)
C
C          IENTY = 2
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = X component of field vector (Gauss)
C              BEAST = Y component of field vector (Gauss)
C              BDOWN = Z component of field vector (Gauss)
C              BABS  = Magnitude of field vector (Gauss)
C          IENTY = 3
C            INPUTS:
C              GLAT = X coordinate (in units of earth radii 6371.2 km )
C              GLON = Y coordinate (in units of earth radii 6371.2 km )
C              ALT  = Z coordinate (in units of earth radii 6371.2 km )
C            RETURNS:
C              BNRTH = Dummy variable
C              BEAST = Dummy variable
C              BDOWN = Dummy variable
C              BABS  = Magnetic potential (T.m)
C
C          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(144),GV(144),ICHG
C            NMAX = Maximum order of spherical harmonic coefficients used
C            GB   = Coefficients for magnetic field calculation
C            GV   = Coefficients for magnetic potential calculation
C            ICHG = Flag indicating when GB,GV have been changed
C
C          HISTORY:
C          COFRM and FELDG originated 15 Apr 83 by Vincent B. Wickwar
C          (formerly at SRI. Int., currently at Utah State).
C
C          May 94 (A.D. Richmond): Added magnetic potential calculation
C
C          Oct 95 (Barnes): Added ICHG
 
      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330,
     +           RE=6371.2, REQ=6378.160)
      COMMON /MAGCOF/ GB(144),GV(144)
      COMMON /MAGCOF_int/ NMAX,ICHG
      DIMENSION XI(3),H(144),G(144)
      SAVE IENTYP, G
      DATA IENTYP/-10000/
 
      IF (IENTY .EQ. 1) THEN
        IS=1
        RLAT = GLAT*DTOR
        CT   = SIN(RLAT)
        ST   = COS(RLAT)
        RLON = GLON*DTOR
        CP   = COS(RLON)
        SP   = SIN(RLON)
        CALL GD2CART (GLAT,GLON,ALT,XXX,YYY,ZZZ)
        XXX = XXX/RE
        YYY = YYY/RE
        ZZZ = ZZZ/RE
      ELSE
        IS   = 2
        XXX  = GLAT
        YYY  = GLON
        ZZZ  = ALT
      ENDIF
      RQ    = 1./(XXX**2+YYY**2+ZZZ**2)
      XI(1) = XXX*RQ
      XI(2) = YYY*RQ
      XI(3) = ZZZ*RQ
      IHMAX = NMAX*NMAX+1
      LAST  = IHMAX+NMAX+NMAX
      IMAX  = NMAX+NMAX-1
 
      IF (IENTY .NE. IENTYP .OR. ICHG .EQ. 1) THEN
        IENTYP = IENTY
	ICHG = 0
        IF (IENTY .NE. 3) THEN
	  DO 10 I=1,LAST
   10     G(I) = GB(I)
        ELSE
	  DO 20 I=1,LAST
   20     G(I) = GV(I)
        ENDIF
      ENDIF
 
      DO 30 I=IHMAX,LAST
   30 H(I) = G(I)

      MK = 3
      IF (IMAX .EQ. 1) MK=1

      DO 100 K=1,MK,2
      I  = IMAX
      IH = IHMAX

   60 IL = IH-I
      F = 2./FLOAT(I-K+2)
      X = XI(1)*F
      Y = XI(2)*F
      Z = XI(3)*(F+F)

      I = I-2
      IF (I .LT. 1) GO TO 90
      IF (I .EQ. 1) GO TO 80

      DO 70 M=3,I,2
      IHM = IH+M
      ILM = IL+M
      H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1))
     +                                        -Y*(H(IHM+2)+H(IHM-2))
   70 H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2))
     +                                        +Y*(H(IHM+3)+H(IHM-1))

   80 H(IL+2) = G(IL+2) + Z*H(IH+2) + X*H(IH+4) - Y*(H(IH+3)+H(IH))
      H(IL+1) = G(IL+1) + Z*H(IH+1) + Y*H(IH+4) + X*(H(IH+3)-H(IH))

   90 H(IL)   = G(IL)   + Z*H(IH)   + 2.*(X*H(IH+1)+Y*H(IH+2))
      IH = IL
      IF (I .GE. K) GO TO 60
  100 CONTINUE
 
      S = .5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      T = (RQ+RQ)*SQRT(RQ)
      BXXX = T*(H(3)-S*XXX)
      BYYY = T*(H(4)-S*YYY)
      BZZZ = T*(H(2)-S*ZZZ)
      BABS = SQRT(BXXX**2+BYYY**2+BZZZ**2)
      IF (IS .EQ. 1) THEN            ! (convert back to geodetic)
        BEAST = BYYY*CP-BXXX*SP
        BRHO  = BYYY*SP+BXXX*CP
        BNRTH =  BZZZ*ST-BRHO*CT
        BDOWN = -BZZZ*CT-BRHO*ST
      ELSEIF (IS .EQ. 2) THEN        ! (leave in earth centered cartesian)
        BNRTH = BXXX
        BEAST = BYYY
        BDOWN = BZZZ
      ENDIF
 
C          Magnetic potential computation makes use of the fact that the
C          calculation of V is identical to that for r*Br, if coefficients
C          in the latter calculation have been divided by (n+1) (coefficients
C          GV).  Factor .1 converts km to m and gauss to tesla.
      IF (IENTY.EQ.3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1
 
      RETURN
      END
 
      SUBROUTINE GD2CART (GDLAT,GLON,ALT,X,Y,Z)
C          Convert geodetic to cartesian coordinates by calling CONVRT
C          940503 A. D. Richmond
      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330)
      CALL CONVRT (1,GDLAT,ALT,RHO,Z)
      ANG = GLON*DTOR
      X = RHO*COS(ANG)
      Y = RHO*SIN(ANG)
      RETURN
      END
 
      SUBROUTINE CONVRT (I,GDLAT,ALT,X1,X2)
C          Convert space point from geodetic to geocentric or vice versa.
C
C          I is an input flag controlling the meaning and direction of the
C            remaining formal arguments:
C
C          I = 1  (convert from geodetic to cylindrical)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance above (north of) Earth's equatorial plane (km)
C
C          I = 2  (convert from geodetic to geocentric spherical)
C            INPUTS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C            RETURNS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C
C          I = 3  (convert from cylindrical to geodetic)
C            INPUTS:
C              X1    = Distance from Earth's rotation axis (km)
C              X2    = Distance from Earth's equatorial plane (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C          I = 4  (convert from geocentric spherical to geodetic)
C            INPUTS:
C              X1    = Geocentric latitude (deg)
C              X2    = Geocentric distance (km)
C            RETURNS:
C              GDLAT = Geodetic latitude (deg)
C              ALT   = Altitude above reference ellipsoid (km)
C
C
C          HISTORY:
C          940503 (A. D. Richmond):  Based on a routine originally written
C          by V. B. Wickwar.
C
C          REFERENCE:  ASTRON. J. VOL. 66, p. 15-16, 1961.
 
      PARAMETER (RTOD=57.2957795130823, DTOR=0.01745329251994330 ,
     +  RE=6371.2 , REQ=6378.160 , FLTNVRS=298.25 ,
     +  E2=(2.-1./FLTNVRS)/FLTNVRS , E4=E2*E2 , E6=E4*E2 , E8=E4*E4 ,
     +  OME2REQ = (1.-E2)*REQ ,
     +     A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. ,
     +     A22 =     (                        E6 +     E8)/  32. ,
     +     A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. ,
     +     A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. ,
     +     A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. ,
     +     A43 =                                   15.*E8 / 256. ,
     +     A44 =                                      -E8 /  16. ,
     +     A61 =  3.*(                     4.*E6 +  5.*E8)/1024. ,
     +     A62 = -3.*(                        E6 +     E8)/  32. ,
     +     A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. ,
     +     A81 =                                   -5.*E8 /2048. ,
     +     A82 =                                   64.*E8 /2048. ,
     +     A83 =                                 -252.*E8 /2048. ,
     +     A84 =                                  320.*E8 /2048. )
C          E2 = Square of eccentricity
 
      IF (I .GE. 3) GO TO 300
 
C          Geodetic to geocentric
 
C          Compute RHO,Z
      SINLAT = SIN(GDLAT*DTOR)
      COSLAT = SQRT(1.-SINLAT*SINLAT)
      D      = SQRT(1.-E2*SINLAT*SINLAT)
      Z      = (ALT+OME2REQ/D)*SINLAT
      RHO    = (ALT+REQ/D)*COSLAT
      X1 = RHO
      X2 = Z
      IF (I .EQ. 1) RETURN
 
C          Compute GCLAT,RKM
      RKM   = SQRT(Z*Z + RHO*RHO)
      GCLAT = RTOD*ATAN2(Z,RHO)
      X1 = GCLAT
      X2 = RKM
      RETURN
 
C          Geocentric to geodetic
  300 IF (I .EQ. 3) THEN
         RHO = X1
         Z = X2
         RKM = SQRT(Z*Z+RHO*RHO)
         SCL = Z/RKM
         GCLAT = ASIN(SCL)*RTOD
      ELSEIF (I .EQ. 4) THEN
         GCLAT = X1
         RKM = X2
         SCL = SIN(GCLAT*DTOR)
      ELSE
         RETURN
      ENDIF
 
      RI = REQ/RKM
      A2 = RI*(A21+RI*(A22+RI* A23))
      A4 = RI*(A41+RI*(A42+RI*(A43+RI*A44)))
      A6 = RI*(A61+RI*(A62+RI* A63))
      A8 = RI*(A81+RI*(A82+RI*(A83+RI*A84)))
      CCL = SQRT(1.-SCL*SCL)
      S2CL = 2.*SCL*CCL
      C2CL = 2.*CCL*CCL-1.
      S4CL = 2.*S2CL*C2CL
      C4CL = 2.*C2CL*C2CL-1.
      S8CL = 2.*S4CL*C4CL
      S6CL = S2CL*C4CL+C2CL*S4CL
      DLTCL = S2CL*A2+S4CL*A4+S6CL*A6+S8CL*A8
      GDLAT = DLTCL*RTOD+GCLAT
      SGL = SIN(GDLAT*DTOR)
      ALT = RKM*COS(DLTCL)-REQ*SQRT(1.-E2*SGL*SGL)
      RETURN
      END
!-----------------------------------------------------------------------
      block data blkmagcof
      implicit none
      integer :: nmax,ichg
      real :: g,gv
      COMMON /MAGCOF/ G(144),GV(144)
      COMMON /MAGCOF_int/ NMAX,ICHG
      DATA NMAX,ICHG /10,-99999/
      end block data blkmagcof

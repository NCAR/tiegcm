#if defined(INTERCOMM) || defined(CISMAH)

!-----------------------------------------------------------------------
      subroutine cism_input_analysis (model_hpi)
! This subroutine will:
!  1) calculate MIX energy flux geflux from gflx and geng 
!  2) calculate the Hemispheric power for each hemisphere and define power as the average
! ... Use Association ..................................................
!      use cism_coupling_module,only: gpot,geng,gflx
      use aurora_module,only:	! (nlonp1,0:nlatp1)
     |     gekev, ! electron energy in keV in geographic coordinates (also at poles for mag2geo)
     |     gmwpm2 ! electron energy flux mW/m2 in geographic coordinates (also at poles for mag2geo)

      use cism_coupling_module,only:
     |     gpot ! potential in geographic coordinates, periodic bound.
!     |     geng, ! characteristic energy in geographic coordinates, periodic boundary
!     |     gflx  ! number flux in geographic coordinates, periodic boundary
      use cons_module,only: ! see cons.F
     |  dtr          ! degrees to radians (pi/180)
      use params_module,only: 
     |  nlon,   ! number of geographic longitudes (at 5 deg, nlon==72)
     |  nlonp1, ! nlon+1
     |  nlonp4, ! nlon+4
     |  nlonp2, ! nlon+2
     |  nlev,   ! number of pressure levels (e.g., 28)
     |  nlevp1, ! nlev+1
     |  nlat,   ! number of geographic latitudes (at 5 deg, nlat==36)
     |  nlatp1, ! nlat+1
     |  nlatp2, ! nlat+2
     |  nmlon,  ! number of geomagnetic grid longitudes
     |  nmlonp1,! nmlon+1
     |  nmlat,  ! number of geomagnetic grid latitudes
     |  nmlatp1,! nmlat+1
     |  nmlath, ! (nmlat+1)/2 (index to magnetic equator)
     |  nmlev,  ! number of geomagnetic pressure levels (nmlev==nlevp1+3)
     |  nmlevp1,! number of geomagnetic midpoint levels
     |  nimlevp1! number of geomagnetic interface levels
! ... Constants ........................................................      
      real, parameter:: fac_p2e = 3.204e-9 ! convert from particle to energy flux
                                           ! eflux=flux*fac_p2e*alfa
                                           !(fac_p2e=2*1.602e-9, so 2*alfa in this factor)      

      real, parameter :: Re=6371.e+3  ! radius of the Earth in m at the level of the ionosphere (use 6370e+3 in MIX)

! Args:
      real,intent(out) :: model_hpi(2)

! ... Local variables ..................................................      
      real cism_power   ! hemispheric power   (GW)  (eg., 16.)

      real, dimension(nlonp1,0:nlatp1) :: geflux ! energy flux from CMIT in geog coords, 
                                                 ! but with poles at 0 and nlatp1 (REQUIRED for geo2mag!!)
      real :: areatot,dglat,dglon,ghpi_sh_mix,ghpi_nh_mix
      real, dimension(nlat) :: da,glatdeg
      integer :: nl4,i,jj,j
      real :: pcp_sh_mix,pcp_nh_mix   ! cross-tail potential drop in kV for SH and NH from MIX
      real :: pcp_sh_min,pcp_nh_min,pcp_sh_max,pcp_nh_max  ! min/max kV potentials
      real :: ghpi_sh_max,ghpi_nh_max
      real, dimension(nmlat/4-1) :: magSdeg,magNdeg
      real :: maxefx(2)  ! Peak flux in Southern (1) and Nothern (2) hemispheres
      integer :: i180
! ... Begin ............................................................

!11111111111111111111111111111111111111111111111111111111
!  1) calculate MIX energy flux geflux from gflx and geng 
!11111111111111111111111111111111111111111111111111111111

! Should calculate hemispheric power from MIX energy flux in GW for both hem and store in power
!  Calculate min/max CP and hpi_s,nh_amie in GW
! Can calculate HP in geog coord instead of in magnetic coords
      dglat = 180./nlat  !FIXME: Off by 1 error?
      dglon = 360./nlon  !FIXME: Off by 1 error?
      areatot = 0.
      do j=1,nlat
       glatdeg(j) = -90.0 + dglat/2. + (j-1)*dglat
       !  NOTE:  The HP from this calc at the cos(lat) of the data is ~12% less than the
       !         MIX calcFaceArea python calculations.  Can get good agreement in the MIX grid
       !         if change da(j) to be for the j-1 latitude equatorwards of the j data.  (Won't do)
       da(j) = Re*cos(glatdeg(j)*dtr)*dglon*dtr*(Re*dglat*dtr)
       areatot = areatot + da(j)*nlon
      enddo

! Find peak fluxes in Southern Hemisphere (1) and Northern Hemisphere (2)
      maxefx(1) = 0.
      maxefx(2) = 0.
      do j=1,nlat
       do i=1,nlon
! NOTE: ting_eng_interp is divided by 2 from MIX before entering TIEGCM, but geng=2*ting_eng_interp
! 11/08 EMERY added eflux=flux*fac_p2e*alfa where alfa=mean energy/2  (fac_p2e=2*1.602e-9)
!       geflux(i,j) = gflx(i+2,j)*fac_p2e*geng(i+2,j)
        geflux(i,j) = gmwpm2(i,j)
	if (j .lt. nlat/2) maxefx(1) = max(geflux(i,j),maxefx(1))
	if (j .ge. nlat/2) maxefx(2) = max(geflux(i,j),maxefx(2))
       enddo
       geflux(nlonp1,j) = geflux(1,j)	! periodic point
      enddo

!  Now calculate geflux(i,0) and geflux(i,nlatp1) at across the SH and NH poles for geo2mag
!   (NOTE:  0 never used since ig goes from 1-36 and geo2mag uses 1-37)
      do i=1,nlon
       i180 = i + nlon/2
       if (i180 .gt. nlonp1) i180 = i180 - nlon
       geflux(i,0) = geflux(i180,1)
       geflux(i,nlatp1) = geflux(i180,nlat)
      enddo
     
! 22222222222222222222222222222222222222222222222222222222222222222222222222222222222222
! 2) calculate the Hemispheric power for each hemisphere and define power as the average
! 22222222222222222222222222222222222222222222222222222222222222222222222222222222222222
! Calculate Hemispheric power (HP) for each hemisphere from the pole to ~45 deg (nlat/4)
! SH
      nl4 = nlat/4+1

      ghpi_sh_mix = 0.
      ghpi_sh_max = 0.
      do i=1,nlon
       do j=1,nl4
        ghpi_sh_mix = ghpi_sh_mix + da(j)*geflux(i,j)
        ghpi_sh_max = max(geflux(i,j),ghpi_sh_max)
       enddo
      enddo
	
      ghpi_nh_mix = 0.
      ghpi_nh_max = 0.
      do i=1,nlon
       do jj=1,nl4
	j = nlat - jj + 1
        ghpi_nh_mix = ghpi_nh_mix + da(j)*geflux(i,j)
        ghpi_nh_max = max(geflux(i,j),ghpi_nh_max)
       enddo
      enddo

! Convert from mW/m2 * m2(da) * 1.e-12GW/mW to GW
      ghpi_sh_mix = ghpi_sh_mix*1.e-12  
      ghpi_nh_mix = ghpi_nh_mix*1.e-12  
! 8/13 bae:  Save HP in model_hpi(SH=1,NH=2) to alter power in advance.F for cusp and drizzle 
      model_hpi(1) = ghpi_sh_mix
      model_hpi(2) = ghpi_nh_mix
      cism_power = 0.5*(ghpi_sh_mix + ghpi_nh_mix) !FIXME: shouldn't modify input parameters!


!666666666666666666666666666666666666666666666666
! 6) calculate ctpoten from the average of the NH
!    and SH MIX potential drops
!666666666666666666666666666666666666666666666666
!  This section on finding min/max from gpot can be deleted unless
!    it is desired to use ctpoten from gpot instead of from phihm

!  Find gpot min/max potentials
       pcp_sh_mix = 0.
       pcp_nh_mix = 0.
       pcp_sh_min = 0.
       pcp_nh_min = 0.
       pcp_sh_max = 0.
       pcp_nh_max = 0.
      do i=1,nlon
       do j=1,nl4
        pcp_sh_min = min(gpot(i+2,j),pcp_sh_min)
        pcp_sh_max = max(gpot(i+2,j),pcp_sh_max)
       enddo
       do jj=1,nl4
	j = nlat - jj + 1
        pcp_nh_min = min(gpot(i+2,j),pcp_nh_min)
        pcp_nh_max = max(gpot(i+2,j),pcp_nh_max)
       enddo
      enddo
!  Convert from V to kV
       pcp_sh_min = pcp_sh_min*1.e-3
       pcp_nh_min = pcp_nh_min*1.e-3
       pcp_sh_max = pcp_sh_max*1.e-3
       pcp_nh_max = pcp_nh_max*1.e-3
       pcp_sh_mix = pcp_sh_max - pcp_sh_min
       pcp_nh_mix = pcp_nh_max - pcp_nh_min

      end subroutine cism_input_analysis
!-----------------------------------------------------------------------
#else
!-----------------------------------------------------------------------
      ! Intel Fortran compiler chokes on empty source files.  
      ! This subroutine is empty so this file will have SOMETHING in it
      subroutine cism_input_analysis_null
      end subroutine cism_input_analysis_null
!-----------------------------------------------------------------------
#endif

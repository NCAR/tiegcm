!
      subroutine xray(tn,o2,o1,barm,sco2,sco1,scn2,
     |  lev0,lev1,lon0,lon1,lat)
!
! Calculate production of ions by solar x-rays.
!
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |  p0,expzmid,expzmid_inv,expz,boltz,check_exp
      use qrj_module,only: ! (nlevp1,nlonp4)
     |  qo2p,  ! o2+ ionization
     |  qop,   ! o+  ionization
     |  qn2p,  ! n2+ ionization
     |  qnp    ! n+  ionization
      use init_module,only: sfeps  ! flux variation from orbital excentricity
      use input_module,only: f107  ! 10.7 cm flux (from input and/or gpi)
      use chemrates_module,only: disn2p
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,             ! neutral temperature (deg K)
     |  o2,o1,          ! mass mixing ratios o2,o
     |  barm,           ! mean molecular weight
     |  sco2,sco1,scn2  ! slant-line column integrals from chapman
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: k,i,nlevs,ier
      real,parameter :: 
     |  sf   = 1.0,
     |  epsx = 5.0
      real :: ex
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  tnn,    ! tn at interfaces
     |  o2n,    ! o2 number density at interfaces
     |  o1n,    ! o  number density at interfaces
     |  n2n,    ! n2 number density at interfaces
     |  xo2,    ! o2 ionization from xrays
     |  xo1,    ! o  ionization from xrays
     |  xn2,    ! n2 ionization from xrays
     |  tau,    ! 
!
! This local xnmbar is not f4d(xnmbari), probably because tni 
! (tn at interfaces) boundaries are extrapolated here, whereas 
! they are not in xnmbari (see sub addiag)
     |  xnmbar  ! p0*e(-z)*barm at interfaces
!
! expo() (util.F) is used only if check_exp is true. This will avoid
! NaNS fpe, but will degrade performance. Check_exp is in cons.F.
      real,external :: expo
!
#ifdef VT
!     code = 129 ; state = 'xray' ; activity='ModelCode'
      call vtbegin(129,ier)
#endif

      ex=(0.3+0.5*(f107-70.)/240.)*sfeps ! snoe
!
! Number of levels for addfsech calls:
      nlevs = lev1-lev0+1
!
! o2,o at interface levels (extrapolate boundaries):
      do i=lon0,lon1
        do k=lev0,lev1-2
          tnn(k+1,i) = 0.5*(tn(k,i)+tn(k+1,i))
          o2n(k+1,i) = 0.5*(o2(k,i)+o2(k+1,i))
          o1n(k+1,i) = 0.5*(o1(k,i)+o1(k+1,i))
        enddo ! k=lev0,lev1-2
        tnn(lev0,i) = 1.5*tn(lev0,i)-0.5*tn(lev0+1,i)
        o2n(lev0,i) = 1.5*o2(lev0,i)-0.5*o2(lev0+1,i)
        o1n(lev0,i) = 1.5*o1(lev0,i)-0.5*o1(lev0+1,i)
        tnn(lev1,i) = 1.5*tn(lev1-1,i)-0.5*tn(lev1-2,i)
        o2n(lev1,i) = 1.5*o2(lev1-1,i)-0.5*o2(lev1-2,i)
        o1n(lev1,i) = 1.5*o1(lev1-1,i)-0.5*o1(lev1-2,i)
      enddo
!
! n*mbar:
      do i=lon0,lon1
        do k=lev0,lev1-1
          xnmbar(k,i) = p0*expzmid_inv*expz(k)*barm(k,i)/
     |      (boltz*tnn(k,i))
        enddo ! k=lev0,lev1-1
        xnmbar(lev1,i) = p0*expzmid*expz(lev1-1)*barm(lev1,i)/
     |      (boltz*tnn(lev1,i))
      enddo ! i=lon0,lon1

!     call addfsech('XNMBAR',' ',' ',xnmbar,lon0,lon1,nlevs,nlevs,
!    |  lat)
!
! Convert from mass mixing ratio to number density:
      do i=lon0,lon1
        do k=lev0,lev1
          n2n(k,i) = (1.-o2n(k,i)-o1n(k,i))*rmassinv_n2*xnmbar(k,i)
          o2n(k,i) = o2n(k,i)*rmassinv_o2*xnmbar(k,i)
          o1n(k,i) = o1n(k,i)*rmassinv_o1*xnmbar(k,i)
        enddo
      enddo
!
! Calculate ionization from xrays:
      do i=lon0,lon1
        do k=lev0,lev1
          if (.not.check_exp) then
            tau(k,i) = 2.52e9*ex*epsx*exp(-(4.40e-19*sco2(k,i)+
     |                                      2.00e-19*sco1(k,i)+
     |                                      2.35e-19*scn2(k,i)))
          else
            tau(k,i) = 2.52e9*ex*epsx*expo(-(4.40e-19*sco2(k,i)+
     |                                       2.00e-19*sco1(k,i)+
     |                                       2.35e-19*scn2(k,i)),0)
         
          endif
          xo2(k,i) = 4.40e-19*o2n(k,i)*tau(k,i)
          xo1(k,i) = 2.00e-19*o1n(k,i)*tau(k,i)
          xn2(k,i) = 2.35e-19*n2n(k,i)*tau(k,i)
        enddo
      enddo
!
! Add to total ionization rates:
      do i=lon0,lon1
        do k=lev0,lev1
          qo2p(k,i,lat) = qo2p(k,i,lat)+0.67*xo2(k,i)
          qop (k,i,lat) = qop (k,i,lat)+xo1(k,i)+0.33*xo2(k,i)
          qn2p(k,i,lat) = qn2p(k,i,lat)+0.64*xn2(k,i)
          qnp (k,i,lat) = qnp (k,i,lat)+0.36*xn2(k,i)
          disn2p(k,i,lat) = disn2p(k,i,lat)+0.64*xn2(k,i)
        enddo
      enddo

!     call addfsech('QO2P',' ',' ',qo2p(:,:,lat),lon0,lon1,nlevs,nlevs,
!    |  lat)
!     call addfsech('QOP' ,' ',' ',qop(:,:,lat) ,lon0,lon1,nlevs,nlevs,
!    |  lat)
!     call addfsech('QN2P',' ',' ',qn2p(:,:,lat),lon0,lon1,nlevs,nlevs,
!    |  lat)
!     call addfsech('QNP' ,' ',' ',qnp(:,:,lat) ,lon0,lon1,nlevs,nlevs,
!    |  lat)
!     call addfsech('DISN2P',' ',' ',disn2p(:,:,lat),lon0,lon1,nlevs,
!    |  nlevs,lat)
!
#ifdef VT
!     code = 129 ; state = 'xray' ; activity='ModelCode'
      call vtend(129,ier)
#endif

      end subroutine xray

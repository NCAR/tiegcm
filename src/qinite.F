!
      subroutine qinite(o2,o1,n2,no,xnmbari,vo2,vo1,vn2,
     |  lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate background ionization rates and add to rates from qrj. 
!
      use lbc,only: fb,b
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |  rmassinv_no
      use qrj_module,only:
     |  qo2p,  ! o2+ ionization
     |  qop,   ! o+  ionization
     |  qn2p,  ! n2+ ionization
     |  qnp,   ! n+  ionization
     |  qnop,  ! no+ ionization
     |  qtef   ! n production rate
      use chemrates_module,only: beta9n
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  o2,o1,n2,no,   ! mass mixing ratios o2,o,n2,no
     |  vo2,vo1,vn2,   ! vertical column densities
     |  xnmbari        ! p0*e(-z)*barm at interfaces
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,k,l,nlevs,ier
      real :: sa(3,3),si(3,3)
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  o2i,    ! o2  at interfaces
     |  o1i,    ! o   at interfaces
     |  n2i,    ! n2  at interfaces
     |  tau,    ! 
     |  qbo2,   ! o2 background ionization
     |  qbo1,   ! o  background ionization
     |  qbn2    ! n2 background ionization
!
! al is sometimes adjusted for different conditions ("tuning knob"):
!      real,parameter :: al(3) = (/5.E4, 5.E3, 5.E3/)  !org 1.1.1
! Using these timegcm values allows removal of minimum Ne in elden.F:
      real,parameter :: al(3)=(/1.5E7,1.5E6,1.5E6/) ! timegcm version
!
#ifdef VT
!     code = 128 ; state = 'qinite' ; activity='ModelCode'
      call vtbegin(128,ier)
#endif
!
      sa(:,1) = (/1.6E-18,  0.,        0./)
      sa(:,2) = (/22.0E-18, 10.24E-18, 23.11E-18/)
      sa(:,3) = (/16.0E-18,  8.40E-18, 11.61E-18/)
!
      si(:,1) = (/ 1.0E-18,  0.,        0./)
      si(:,2) = (/22.0E-18, 10.24E-18, 23.11E-18/)
      si(:,3) = (/16.0E-18,  8.40E-18, 11.61E-18/)
!
! Number of levels for addfsech calls:
      nlevs = lev1-lev0+1
!
! o2,o at interface levels:
      do k=lev0,lev1-1
        o2i(k+1,:) = 0.5*(o2(k,lon0:lon1)+o2(k+1,lon0:lon1))
        o1i(k+1,:) = 0.5*(o1(k,lon0:lon1)+o1(k+1,lon0:lon1))
        n2i(k+1,:) = 0.5*(n2(k,lon0:lon1)+n2(k+1,lon0:lon1))
      enddo
!
! Bottom boundary for o2i, o1i:
      do i=lon0,lon1
        o2i(1,i) = .5*((b(i,1,1)+1.)*o2(1,i)+b(i,1,2)*o1(1,i)+fb(i,1))
        o1i(1,i) = .5*(b(i,2,1)*o2(1,i)+(b(i,2,2)+1.)*o1(1,i)+fb(i,2))
!
! btf tiegcm_he 8/21/13: Not sure about n2i bottom boundary. 
! For now, use 1-o2i-o1i since helium is negligible at this height.
!
        n2i(1,i) = 1.-o2i(1,i)-o1i(1,i)
      enddo
!
! n2 at interfaces:
      do k=lev0,lev1
        do i=lon0,lon1                       !!! temporal fix CISM
         if(n2i(k,i)<=0.)n2i(k,i)=1.e-5      !!! temporal fix CISM
        enddo                                !!! temporal fix CISM
      enddo 
!
! Summation over wavelength:
      qbo2 = 0. ! array init
      qbo1 = 0. ! array init
      qbn2 = 0. ! array init
      do l=1,3
        tau(:,:) = 0.
        do i=lon0,lon1
          do k=lev0,lev1
            tau(k,i) = tau(k,i)+
     |        sa(1,l)*vo2(k,i) + sa(2,l)*vo1(k,i) + sa(3,l)*vn2(k,i)
            tau(k,i) = exp(-tau(k,i))
            qbo2(k,i) = qbo2(k,i)+al(l)*si(1,l)*o2i(k,i)*tau(k,i)*
     |        rmassinv_o2
            qbo1(k,i) = qbo1(k,i)+al(l)*si(2,l)*o1i(k,i)*tau(k,i)*
     |        rmassinv_o1
            qbn2(k,i) = qbn2(k,i)+al(l)*si(3,l)*n2i(k,i)*tau(k,i)*
     |        rmassinv_n2
          enddo
        enddo
      enddo ! l=1,3
!
!     call addfsech('TAU' ,' ',' ',tau ,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QBO2',' ',' ',qbo2,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QBO1',' ',' ',qbo1,lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QBN2',' ',' ',qbn2,lon0,lon1,nlevs,nlevs,lat)
!
! Add contributions to o2+, o+, n2+, n+, and no+ ionization rates:
      do i=lon0,lon1
        do k=lev0,lev1
          qo2p(k,i,lat) = qo2p(k,i,lat)+0.67*qbo2(k,i)*xnmbari(k,i)
          qop (k,i,lat) = qop (k,i,lat)+(0.33*qbo2(k,i)+qbo1(k,i))*
     |      xnmbari(k,i)
          qn2p(k,i,lat) = qn2p(k,i,lat)+0.86*qbn2(k,i)*xnmbari(k,i)
          qnp (k,i,lat) = qnp (k,i,lat)+0.14*qbn2(k,i)*xnmbari(k,i)
          qtef(k,i,lat) = qtef(k,i,lat)+
     |                    1.57*(0.86*qbn2(k,i)*xnmbari(k,i))
        enddo
      enddo
!
! NO ionization and add to qnop:
! Bottom boundary:
      do i=lon0,lon1
        qnop(1,i,lat) = qnop(1,i,lat)+beta9n(1,i,lat)*no(1,i)*
     |    xnmbari(1,i)*rmassinv_no
      enddo
!
! qnop at levels 2->top
      do k=2,lev1
        do i=lon0,lon1
          qnop(k,i,lat) = qnop(k,i,lat)+beta9n(k,i,lat)*
     |      0.5*(no(k,i)+no(k-1,i))*xnmbari(k,i)*rmassinv_no
        enddo ! i=lon0,lon1
      enddo ! k=2,lev1
 
!     call addfsech('QO2P',' ',' ',qo2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QOP' ,' ',' ',qop (lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QN2P',' ',' ',qn2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QNP' ,' ',' ',qnp (lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QNOP',' ',' ',qnop(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!
#ifdef VT
!     code = 128 ; state = 'qinite' ; activity='ModelCode'
      call vtend(128,ier)
#endif
      end subroutine qinite

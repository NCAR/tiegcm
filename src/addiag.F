      subroutine addiag(tn,o2,o1,he,vn,vc,barm,xnmbar,xnmbari,xnmbarm,z,
     |  zg,n2,lon0,lon1,lev0,lev1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate needed terms vc, barm, xnmbar[i,m], and Z:
!
      use cons_module,only: cs,rmassinv,rmassinv_he,dz,dzgrav,
     |  freq_semidi,dt,p0,boltz,expz,expzmid,expzmid_inv,zbound
      use init_module,only: iter,igetgswm
      use fields_module,only: tlbc
      use addfld_module,only: addfld
      use lbc,only: z_lbc
      use diags_module,only: mkdiag_SCHT,mkdiag_ZGMID
      implicit none
!
! Input args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn, ! neutral temperature (deg K)
     |  o2, ! molecular oxygen (mmr)
     |  o1, ! atomic oxygen (mmr)
     |  he, ! helium (mmr)
     |  vn  ! meridional wind velocity (cm/s)
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  vc  ,
     |  barm,
     |  z   ,
     |  zg  ,
     |  n2  ,
     |  xnmbar , 
     |  xnmbari,
     |  xnmbarm
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: k,i,j,nlevs,ier
      real :: 
     |  barm1(lon0:lon1),
     |  tni  (lev0:lev1,lon0:lon1), ! tn at interfaces
     |  expzi(lev0:lev1,lon0:lon1), ! e(-z) at interfaces
     |  w1   (lev0:lev1,lon0:lon1)
      complex :: expt
      real :: fmin,fmax
!
#ifdef VT
!     code = 116 ; state = 'addiag' ; activity='ModelCode'
      call vtbegin(116,ier)
#endif
      nlevs = lev1-lev0+1
!
! Save inputs:
!     do j=lat0,lat1
!       call addfld('diag_tn',' ',' ',tn(lev0:lev1-1,lon0:lon1,j),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!       call addfld('diag_o2',' ',' ',o2(lev0:lev1-1,lon0:lon1,j),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!       call addfld('diag_o1',' ',' ',o1(lev0:lev1-1,lon0:lon1,j),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!     enddo
!
! Latitude scan:
      do j=lat0,lat1
!
! Define 3d N2 (fields.F):
!
        do i=lon0,lon1
          n2(:,i,j) = (1.-o2(:,i,j)-o1(:,i,j)-he(:,i,j)) 
        enddo
!       write(6,"('addiag: j=',i4,' n2(:,lon0:lon1,j)=',2e12.4)")
!    |    j,minval(n2(:,lon0:lon1,j)),maxval(n2(:,lon0:lon1,j))
!
! vc = cos(phi)*v
!
        do i=lon0,lon1
          do k=lev0,lev1
            vc(k,i,j) = cs(j)*vn(k,i,j)
          enddo
        enddo
!       call addfld('diag_vc',' ',' ',vc(lev0:lev1-1,lon0:lon1,j),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!
! barm = mean molecular weight (k+1/2):
!
        do i=lon0,lon1
          do k=lev0,lev1
            barm(k,i,j) = 1./ 
     |        (o2(k,i,j)*rmassinv(1)+o1(k,i,j)*rmassinv(2)+
     |         he(k,i,j)*rmassinv_he+n2(k,i,j)*rmassinv(3))
            xnmbarm(k,i,j)=p0*expz(k)*barm(k,i,j)/
     |        (boltz*tn(k,i,j))
          enddo
        enddo
!
! barm1 = barm(k=0) (linear extrapolation)
!
        do i=lon0,lon1
          barm1(i) = 1.5*barm(1,i,j)-0.5*barm(2,i,j)
        enddo
!       call addfld('barm',' ',' ',barm(:,lon0:lon1,j),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! barm(k) = 0.5*(barm(k+1/2)+barm(k-1/2)), k = kmaxp1,2,1
!
        do i=lon0,lon1
          do k=lev1,lev0+1,-1
            barm(k,i,j) = 0.5*(barm(k,i,j)+barm(k-1,i,j))
          enddo
        enddo
!
! barm(1) = barm1
!
        do i=lon0,lon1
          barm(lev0,i,j) = barm1(i)
        enddo
!
! xnmbar = p0*e(-z)*barm/kT at midpoints (used in conversion from mmr to cm3).
! (used by oplus)
        do i=lon0,lon1
          do k=lev0,lev1-1
            xnmbar(k,i,j)=p0*expz(k)*.5*(barm(k,i,j)+barm(k+1,i,j))
     |        /(boltz*tn(k,i,j))
          enddo
        enddo
!
! xnmbari = p0*e(-z)*barm/kT at interfaces (used by qrj and qinite):
        do i=lon0,lon1
!         tni(1,i) = tn(lev1,i,j) ! tn bottom boundary is stored in top slot
          tni(lev0,i) = tlbc(i,j)    ! Lower boundary is in tlbc
          expzi(1,i) = expzmid_inv*expz(1)
          do k=lev0+1,lev1-1
            tni(k,i) = .5*(tn(k-1,i,j)+tn(k,i,j))
            expzi(k,i) = expzmid_inv*expz(k)
          enddo
          tni(lev1,i) = tn(lev1-1,i,j) ! nlevp1 <- nlev
          expzi(lev1,i) = expzmid*expz(lev1-1)
          do k=lev0,lev1
            xnmbari(k,i,j) = p0*expzi(k,i)*barm(k,i,j)/
     |        (boltz*tni(k,i))
          enddo
        enddo
!       call addfld('XNMBARI',' ',' ',xnmbari(:,lon0:lon1,j),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!       call addfld('TNI','TNI from addiag','K',
!    |    tni(:,:),'ilev',lev0,lev1,'lon',lon0,lon1,j)
      enddo ! j=lat0,lat1
!
! Calculate geopotential Z:
!
      do j=lat0,lat1
        z(1,lon0:lon1,j) = z_lbc(lon0:lon1,j)
!
! Complete calculation of geopotential Z:
!
! w1 = barm
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = (barm(k,i,j)+barm(k+1,i,j))*0.5
          enddo
        enddo
!       call addfld('W1a',' ',' ',w1,'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! w1 = tn/w1 (old model comment: s1=s2/s1=(t+t0)/m)
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = tn(k,i,j)/w1(k,i)
          enddo
        enddo
!       call addfld('TNa',' ',' ',tn(:,lon0:lon1,j),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!       call addfld('W1b',' ',' ',w1,'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! w1=(ds*r/g)*w1
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = (dz/dzgrav) * w1(k,i)
          enddo
        enddo
!       call addfld('W1c',' ',' ',w1,'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! Extend Z upward:
        do i=lon0,lon1
          do k=lev0,lev1-1
            z(k+1,i,j) = w1(k,i)+z(k,i,j)
          enddo
        enddo
!       call addfld('ADIAG_Z','Geopotential from addiag',
!    |    'cm',z(:,lon0:lon1,j),'ilev',lev0,lev1,'lon',lon0,lon1,j)
      enddo ! j=lat0,lat1
!
! Calculate diagnostic geopotential zg, with varying gravity:
      call calczg(tn,o2,o1,n2,he,z,zg,lon0,lon1,lev0,lev1,lat0,lat1)
!
! Interpolate Z from interfaces to midpoints for diagnostic:
      do j=lat0,lat1
        call mkdiag_ZGMID('ZGMID',zg(:,lon0:lon1,j),lev0,lev1,
     |    lon0,lon1,j)
      enddo
!
! Calculate scale height diagnostic (using Z here, not ZG):
      call mkdiag_SCHT('SCHT',z(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
#ifdef VT
!     code = 116 ; state = 'addiag' ; activity='ModelCode'
      call vtend(116,ier)
#endif
      end subroutine addiag
!-----------------------------------------------------------------------
      subroutine calczg(tn,o2,o1,n2,he,z,zg,lon0,lon1,lev0,lev1,
     |  lat0,lat1)
!
! Given geopotential z (calculated with the model constant gravity),
!   calculate geopotential zg with varying gravity. This is taken from
!   tgcmproc_f90, routines calchts and glatf in proclat.F.
! ZG will be put on secondary histories, along with the regular Z.
!
      use params_module,only: dz,glat
      use init_module,only: istep
      use cons_module,only: boltz,avo
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn, ! neutral temperature (deg K)
     |  o2, ! molecular oxygen (mmr)
     |  o1, ! atomic oxygen (mmr)
     |  he, ! helium (mmr)
     |  n2, ! molecular nitrogen (mmr)
     |  z   ! geopotential calculated with constant gravity (from addiag)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) :: 
     |  zg  ! output geopotential calculated with varying gravity
!
! Local:
      integer :: i,j,k
      real :: g0,r0,c2
      real,dimension(lev0:lev1) :: xmas,g
      real,parameter :: dgtr=1.74533E-2
!
! Latitude scan:
! 1/20/10 btf: changed float(j) to glat(j) in cos of c2 calculation
      do j=lat0,lat1
        c2 = cos(2.*dgtr*glat(j))
        g0 = 980.616*(1.-.0026373*c2)
        r0 = 2.*g0/(3.085462e-6 + 2.27e-9*c2) ! effective earth radius
!
! Longitude scan:
        do i=lon0,lon1
          g(1)=g0*(r0/(r0+0.5*(z(1,i,j)+z(2,i,j))))**2
          xmas(:) = 1./(o1(:,i,j)/16.+o2(:,i,j)/32.+he(:,i,j)/4.+
     |                  n2(:,i,j)/28.)/avo
!
! Interface levels:
          zg(lev0,i,j) = z(lev0,i,j)
          do k=lev0+1,lev1-1
            zg(k,i,j) = zg(k-1,i,j) + boltz*dz*tn(k-1,i,j) / 
     |        (xmas(k-1)*g(k-1))
            g(k)=g0*(r0/(r0+0.5*(zg(k,i,j)+z(k+1,i,j))))**2
          enddo ! k=lev0+1,lev1-1
          zg(lev1,i,j) = 2.0*zg(lev1-1,i,j)-zg(lev1-2,i,j)
        enddo ! i=lon0,lon1
!
! Save ZG to secondary histories:
! PLEASE DO NOT COMMENT THIS OUT -- ZG IS A MANDATORY FIELD ON SECH HISTORIES
!
        call addfld('ZG','Geometric Height ZG',
     |    'cm',zg(:,lon0:lon1,j),'ilev',lev0,lev1,'lon',lon0,lon1,j)
      enddo ! j=lat0,lat1
      end subroutine calczg

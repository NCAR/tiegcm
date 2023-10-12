      subroutine addiag(tn,o2,o1,he,vn,vc,mbar,barm,xnmbar,xnmbari,
     |  scht,schti,z,zg,n2,lon0,lon1,lev0,lev1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate needed terms vc, barm, xnmbar[i,m], and Z:
!
      use cons_module,only: cs,rmassinv,rmassinv_he,dz,dzgrav,
     |  gask,grav,p0,boltz,expz,expzmid,expzmid_inv
      use fields_module,only: tlbc
      use addfld_module,only: addfld
      use lbc,only: z_lbc
      use diags_module,only: mkdiag_N2,mkdiag_O_N2,mkdiag_ZGMID
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
     |  mbar,
     |  barm,
     |  scht,
     |  schti,
     |  z   ,
     |  zg  ,
     |  n2  ,
     |  xnmbar , 
     |  xnmbari
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
! Provide diagnostics N2 and O/N2 (moved from comp.F):
!
        call mkdiag_N2('N2',n2(:,lon0:lon1,j),lev0,lev1,lon0,lon1,j)
        call mkdiag_O_N2('O_N2',o1(:,lon0:lon1,j),n2(:,lon0:lon1,j),
     |    lev0,lev1,lon0,lon1,j)
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
! mbar = mean molecular weight
!
        do i=lon0,lon1
          do k=lev0,lev1
            mbar(k,i,j) = 1./ 
     |        (o2(k,i,j)*rmassinv(1)+o1(k,i,j)*rmassinv(2)+
     |         he(k,i,j)*rmassinv_he+n2(k,i,j)*rmassinv(3))
            xnmbar(k,i,j)=p0*expz(k)*mbar(k,i,j)/
     |        (boltz*tn(k,i,j))
          enddo
          xnmbar(lev1,i,j)=p0*expz(lev1-1)*mbar(lev1,i,j)/
     |      (boltz*tn(lev1-1,i,j))
        enddo
!
! barm = mean molecular weight (k+1/2):
! barm1 = mbar(k=0) (linear extrapolation)
!
        do i=lon0,lon1
          barm1(i) = 1.5*mbar(1,i,j)-0.5*mbar(2,i,j)
        enddo
!       call addfld('mbar',' ',' ',mbar(:,lon0:lon1,j),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,j)
!
! barm(k) = 0.5*(mbar(k+1/2)+mbar(k-1/2)), k = kmaxp1,2,1
!
        do i=lon0,lon1
          do k=lev1,lev0+1,-1
            barm(k,i,j) = 0.5*(mbar(k,i,j)+mbar(k-1,i,j))
          enddo
        enddo
!
! barm(1) = barm1
!
        do i=lon0,lon1
          barm(lev0,i,j) = barm1(i)
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
!
! Calculate scale height
        do i=lon0,lon1
          do k=lev0,lev1
            scht(k,i,j) = gask*tn(k,i,j)/(mbar(k,i,j)*grav)
            schti(k,i,j) = gask*tni(k,i)/(barm(k,i,j)*grav)
          enddo
          scht(lev1,i,j) =
     |      gask*(2*tn(lev1-1,i,j)-tn(lev1-2,i,j))/(mbar(lev1,i,j)*grav)
        enddo
      enddo ! j=lat0,lat1
!
! Calculate geopotential Z:
!
      do j=lat0,lat1
        z(1,lon0:lon1,j) = z_lbc(lon0:lon1,j)
!
! Complete calculation of geopotential Z:
!
! w1 = mbar
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = mbar(k,i,j)
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
      call calczg(scht,z,zg,lon0,lon1,lev0,lev1,lat0,lat1)
!
! Interpolate Z from interfaces to midpoints for diagnostic:
      do j=lat0,lat1
        call mkdiag_ZGMID('ZGMID',zg(:,lon0:lon1,j),lev0,lev1,
     |    lon0,lon1,j)
      enddo
!
! Calculate scale height diagnostic (using Z here, not ZG):
! As of TIEGCM v3.0, SCHT is f4d, so this call is no longer needed
!     call mkdiag_SCHT('SCHT',z(:,lon0:lon1,lat0:lat1),
!    |  lev0,lev1,lon0,lon1,lat0,lat1)
!
#ifdef VT
!     code = 116 ; state = 'addiag' ; activity='ModelCode'
      call vtend(116,ier)
#endif
      end subroutine addiag
!-----------------------------------------------------------------------
      subroutine calczg(scht,z,zg,lon0,lon1,lev0,lev1,lat0,lat1)
!
! Given geopotential z (calculated with the model constant gravity),
!   calculate geopotential zg with varying gravity. This is taken from
!   tgcmproc_f90, routines calchts and glatf in proclat.F.
! ZG will be put on secondary histories, along with the regular Z.
!
      use params_module,only: dz,glat
      use cons_module,only: grav
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  scht, ! scale height
     |  z   ! geopotential calculated with constant gravity (from addiag)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) :: 
     |  zg  ! output geopotential calculated with varying gravity
!
! Local:
      integer :: i,j,k
      real :: g0,r0,c2
      real,dimension(lev0:lev1) :: g
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
!
! Interface levels:
          zg(lev0,i,j) = z(lev0,i,j)
          do k=lev0+1,lev1-1
            zg(k,i,j) = zg(k-1,i,j) + dz*scht(k-1,i,j)*grav/g(k-1)
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

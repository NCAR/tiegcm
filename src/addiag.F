!
      subroutine addiag(tn,o2,o1,vn,vc,barm,xnmbar,xnmbari,xnmbarm,z,
     |  lon0,lon1,lev0,lev1,lat0,lat1)
!
! Calculate needed terms vc, barm, xnmbar[i,m], and Z:
!
      use cons_module,only: cs,rmassinv,dz,dzgrav,freq_semidi,dt,p0,
     |  boltz,expz,expzmid,expzmid_inv
      use bndry_module,only: zb,zb2,zba,bnd,bnd2,bnda,ci
      use init_module,only: iter
      implicit none
!
! Input args:
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat0,lat1
!     real,dimension(lev0:lev1,lon0  :lon1  ,lat0  :lat1  ),intent(in):: 
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in):: 
     |  tn, ! neutral temperature (deg K)
     |  o2, ! molecular oxygen (mmr)
     |  o1, ! atomic oxygen (mmr)
     |  vn  ! meridional wind velocity (cm/s)
!
! Output args:
!     real,dimension(lev0:lev1,lon0  :lon1  ,lat0  :lat1  ),
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  vc  ,
     |  barm,
     |  z   ,
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
! Latitude scan:
      do j=lat0,lat1
!
! vc = cos(phi)*v
!
        do i=lon0,lon1
          do k=lev0,lev1
            vc(k,i,j) = cs(j)*vn(k,i,j)
          enddo
        enddo
!       call addfsech('diag_vc',' ',' ',vc,lon0,lon1,nlevs,nlevs-1,j)
!
! barm = mean molecular weight (k+1/2):
!
        do i=lon0,lon1
          do k=lev0,lev1
            barm(k,i,j) = 1./
     |        (o2(k,i,j)*rmassinv(1)+o1(k,i,j)*rmassinv(2)+
     |        (1.-o2(k,i,j)-o1(k,i,j))*rmassinv(3))
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
!       call addfsech('barm',' ',' ',barm,lon0,lon1,nlevs,nlevs-1,j)
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
          tni(1,i) = tn(lev1,i,j) ! tn bottom boundary is stored in top slot
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
!       call addfsech('XNMBARI',' ',' ',xnmbari,lon0,lon1,nlevs,nlevs,j)
!
! w1 = barm
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = (barm(k,i,j)+barm(k+1,i,j))*0.5
          enddo
        enddo
!       call addfsech('W1a',' ',' ',w1,lon0,lon1,nlevs,nlevs-1,j)
!
! w1 = tn/w1 (old model comment: s1=s2/s1=(t+t0)/m)
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = tn(k,i,j)/w1(k,i)
          enddo
        enddo
!       call addfsech('TNa',' ',' ',tn,lon0,lon1,nlevs,nlevs-1,j)
!       call addfsech('W1b',' ',' ',w1,lon0,lon1,nlevs,nlevs-1,j)
!
! w1=(ds*r/g)*w1
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            w1(k,i) = (dz/dzgrav) * w1(k,i)
          enddo
        enddo
!       call addfsech('W1c',' ',' ',w1,lon0,lon1,nlevs,nlevs-1,j)
!
! z(1) = zb
        expt=cexp(ci*freq_semidi*dt*iter)
        do i=lon0,lon1
          z(1,i,j) = real(zb(j)*bnd(i)*expt) 
        enddo
!
! Add in effect of (1,1) tidal component to lbc
!
        expt = cexp(ci*.5*freq_semidi*dt*iter)
        do i=lon0,lon1
          z(1,i,j) = z(1,i,j)+real(zb2(j)*bnd2(i)*expt)
        enddo
!
! Add in effect of annual tidal component to lbc
!
        expt = 1.
        do i=lon0,lon1
          z(1,i,j) = z(1,i,j)+real(zba(j)*bnda(i)*expt)
        enddo
!
! z(k+1)=s1(k)+z(k)
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            z(k+1,i,j) = w1(k,i)+z(k,i,j)
          enddo
        enddo
!       call addfsech('ADIAG_Z',' ',' ',z(:,lon0:lon1,j),
!    |    lon0,lon1,nlevs,nlevs-1,j)

!       do k=lev0,lev1-1
!         write(6,"('addiag: k=',i3,' j=',i3,' z(k,lon0:lon1,j)=',
!    |      /,(6e12.4))") k,j,z(k,lon0:lon1,j)
!       enddo ! k=lev0,lev1-1

      enddo ! j=lat0,lat1
!
#ifdef VT
!     code = 116 ; state = 'addiag' ; activity='ModelCode'
      call vtend(116,ier)
#endif
      end subroutine addiag

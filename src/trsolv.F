!
      subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      implicit none
!
! Tri-diagonal solver.
!   a(k,i)*x(k-1,i) + b(k,i)*x(k,i) + c(k,i)*x(k+1,i) = f(k,i)
! This is called from dt, minor, oplus, and settei.
!
! Args:
      integer,intent(in) :: lev0,lev1,k1,k2,lon0,lon1,lonmax
      real,dimension(lev0:lev1,lon0:lon1),intent(in) ::
     |  a, ! input coefficients
     |  b, ! input coefficients
     |  c, ! input coefficients
     |  f  ! input RHS
      real,dimension(lev0:lev1,lon0:lon1),intent(out) ::
     |  x  ! output
!
! Local:
      integer :: k,i,lonbeg,lonend
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  w1,w2,w3  ! work arrays
!
! Set lonbeg,lonend (e.g., 3->74 at 5x resolution):
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==lonmax) lonend = lon1-2
!
! Lower boundary (W(K1)=B(K1), W(2*KF+K1)=F(K1)/W(K1)):
      do i=lonbeg,lonend
        w1(k1,i) = b(k1,i)
        w2(k1,i) = c(k1,i) / w1(k1,i)
        w3(k1,i) = f(k1,i) / w1(k1,i)
      enddo
!
! Set up work arrays:
      do i=lonbeg,lonend
        do k=k1+1,k2
!
! W(K)=B(K)-A(K)*W(KF+K-1)
          w1(k,i) = b(k,i) - a(k,i)*w2(k-1,i)
!
! W(KF+K)=C(K)/W(K):
          w2(k,i) = c(k,i) / w1(k,i)
!
! W(2*KF+K)=(F(K)-A(K)*W(2*KF+K-1))/W(K)
          w3(k,i) = (f(k,i) - a(k,i)*w3(k-1,i)) / w1(k,i)
        enddo ! k=k1+1,k2
      enddo ! i=lonbeg,lonend
!
! Upper boundary (X(K2)=W(2*KF+K2)):
      do i=lonbeg,lonend
        x(k2,i) = w3(k2,i)       
      enddo ! i=lonbeg,lonend
!
! Back substitution:
      do i=lonbeg,lonend
        do k=k2-1,k1,-1
!
! X(K)=W(2*KF+K)-W(KF+K)*X(K+1)
          x(k,i) = w3(k,i) - w2(k,i)*x(k+1,i)
        enddo ! k=k2-1,k1,-1
      enddo ! i=lonbeg,lonend

      end subroutine trsolv

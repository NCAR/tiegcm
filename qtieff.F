!
      subroutine qtieff(lev0,lev1,lon0,lon1,lat)
!
! Calculates electron production of n4s and n2d in qtef:
!
      use qrj_module,only: qtef,qn2p
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Local:
      integer :: k,i,nlevs
!
      do i=lon0,lon1
        do k=lev0,lev1
          qtef(k,i,lat) = qtef(k,i,lat)+1.57*qn2p(k,i,lat)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

      nlevs = lev1-lev0+1
!     call addfsech('QN2P',' ',' ',qn2p(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)
!     call addfsech('QTEF',' ',' ',qtef(lev0:lev1,lon0:lon1,lat),
!    |  lon0,lon1,nlevs,nlevs,lat)

      end subroutine qtieff

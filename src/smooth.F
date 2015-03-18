!
      subroutine smooth(fin,fout,k0,k1,i0,i1,j0,j1,isave)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! 2-part shapiro smoother, called from dt for tn_nm, and from duv for
!   un_nm, vn_nm. 
! Input field fin is defined at full subdomain (including halo cells).
! Output field fout is defined at subdomain (excluding halo cells and
!   periodic points).
!
      use params_module,only: nlonp4
      use cons_module,only: shapiro
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: k0,k1,i0,i1,j0,j1,isave
!
! Input field at full subdomain (with halo):
      real,dimension(k0:k1,i0-2:i1+2,j0-2:j1+2),intent(in) :: fin
!
! Output field at regular subdomain (without halo):
      real,dimension(k0:k1,i0:i1,j0:j1),intent(out) :: fout
!
! Local:
      integer :: k,i,j,ibeg,iend 
      real,dimension(k0:k1,i0-2:i1+2) :: ftmp
!
      ibeg = i0
      if (i0==1) ibeg = 3
      iend = i1
      if (i1==nlonp4) iend = i1-2
!
! Meridional smoothing (differences in latitude): 
      do j=j0,j1
        do i=i0-2,i1+2
          do k=k0,k1
            ftmp(k,i)=fin(k,i,j)-shapiro*(fin(k,i,j+2)+fin(k,i,j-2)-
     |        4.*(fin(k,i,j+1)+fin(k,i,j-1))+6.*fin(k,i,j))
          enddo ! k=k0,k1
        enddo ! i=i0,i1

        if (isave > 0)
     |    call addfld('SMOO1',' ',' ',ftmp(:,i0:i1),
     |      'lev',k0,k1,'lon',i0,i1,j)
!
! Zonal smoothing (differences in longitude):
        do i=ibeg,iend
          do k=k0,k1
            fout(k,i,j) = ftmp(k,i)-shapiro*(ftmp(k,i+2)+ftmp(k,i-2)-
     |        4.*(ftmp(k,i+1)+ftmp(k,i-1))+6.*ftmp(k,i))
          enddo ! k=k0,k1
        enddo ! i=i0,i1
!
! Set periodic points to zero to avoid NaNS is calling routines
! (e.g., dt.F). This was apparently assumed in earlier versions.
        if (i0==1) fout(:,i0:i0+1,j) = 0.
        if (i1==nlonp4) fout(:,i1-1:i1,j) = 0.

        if (isave > 0)
     |    call addfld('SMOO2',' ',' ',fout(:,i0:i1,j),
     |      'lev',k0,k1,'lon',i0,i1,j)

      enddo ! j=j0,j1
      end subroutine smooth

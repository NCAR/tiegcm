!
      subroutine comp_n2d(o2,o1,no,n4s,ne,op,n2p,nop,sco2,xnmbarm,
     |  n2d,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculates psi(n2d) (mmr)  assuming photochemical equilibrium,
! at current latitude (lat).
!
      use cons_module,only: brn2d,rmassinv_o2,rmassinv_o1,rmassinv_no,
     |  rmass_n2d
      use chemrates_module,only: rk3,ra1,ra3,beta2,beta4,beta5,beta6,
     |  beta7,rk10
      use qrj_module,only: qtef
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Input arrays:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  o2,     ! molecular oxygen (mmr)
     |  o1,     ! atomic oxygen (mmr)
     |  no,     ! nitric oxide (mmr)
     |  n4s,    ! N4S (mmr)
     |  ne,     ! electron density (cm3)
     |  op,     ! O+ ion
     |  n2p,    ! N2+ ion
     |  nop,    ! NO+ ion
     |  sco2,   ! O2 line integral (chapman)
     |  xnmbarm ! p0*e(-z)*barm/kT
!
! Output:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  n2d    ! output
!
! Local:
      integer :: k,i
      integer ::i0,i1,nk,nkm1
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  n2d_prod, ! n2d production
     |  n2d_loss  ! n2d production
!
! For addfld:
      i0 = lon0
      i1 = lon1
      nk = lev1-lev0+1
      nkm1 = nk-1
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! Total production of N2D:
          n2d_prod(k,i) = (.5*(qtef(k,i,lat)+qtef(k+1,i,lat))*brn2d)+
     |      rk3(k,i,lat)*n2p(k,i)*xnmbarm(k,i)*o1(k,i)*rmassinv_o1+
     |     (ra1(k,i,lat)*nop(k,i)*0.85+ra3(k,i,lat)*n2p(k,i)*0.9)*
     |      sqrt(ne(k,i)*ne(k+1,i))
!
! Total loss of N2D:
          n2d_loss(k,i) = xnmbarm(k,i)*
     |      (beta2*o2(k,i)*rmassinv_o2+
     |       beta4*o1(k,i)*rmassinv_o1+
     |       beta6*no(k,i)*rmassinv_no)+
     |       beta7+beta5(k,i,lat)*sqrt(ne(k,i)*ne(k+1,i))+
     |       rk10*op(k,i)
!
! Update n2d:
          n2d(k,i)=rmass_n2d*n2d_prod(k,i)/(n2d_loss(k,i)*xnmbarm(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('N2D_PROD' ,' ',' ',n2d_prod(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('N2D_LOSS' ,' ',' ',n2d_loss(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('N2D_UPD'  ,' ',' ',n2d(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)

      end subroutine comp_n2d

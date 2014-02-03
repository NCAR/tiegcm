!
      subroutine cpktkm(tn,o2,o1,n2,he,fcp,fkt,fkm,lev0,lev1,lon0,lon1,
     |  lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Define diagnostics CP, KT, and KM.
!
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |                      rmassinv_he,gask,t0
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_MU_M
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  tn,  ! neutral temperature (deg K)
     |  o2,  ! molecular oxygen (mmr)
     |  o1,  ! atomic oxygen (mmr)
     |  n2,  ! molecular nitrogen (mmr)
     |  he   ! helium (mmr)
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  fcp, ! specific heat at constant pressure (ergs/deg/gm)
     |  fkt, ! molecular diffusion (ergs/cm/deg/sec)
     |  fkm  ! molecular viscosity (gm/cm/sec)
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  fmbar,        ! mean mass
     |  po2,po1,pn2,phe 
      integer :: nlons,nlevs
!
      nlons = lon1-lon0+1
      nlevs = lev1-lev0+1

!     call addfld('tn_cp',' ',' ',tn(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('o2_cp',' ',' ',o2(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('o1_cp',' ',' ',o1(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)

! Modified by EKS (Eric Sutton):
!  Compute po2, po1, phe, and pn2 as the number density mixing ratios.
!  Add the effect of Helium to the cp,km,kt calculations consistent with
!  Banks and Kockarts, 1973.
      do i=lon0,lon1
        do k=lev0,lev1
          fmbar(k,i) = 1./(o2(k,i)*rmassinv_o2 + o1(k,i)*rmassinv_o1 +
     |                 he(k,i)*rmassinv_he + n2(k,i)*rmassinv_n2)
          po2(k,i) = fmbar(k,i)*o2(k,i)*rmassinv_o2
          po1(k,i) = fmbar(k,i)*o1(k,i)*rmassinv_o1
          phe(k,i) = fmbar(k,i)*he(k,i)*rmassinv_he
          pn2(k,i) = fmbar(k,i)*n2(k,i)*rmassinv_n2
          fcp(k,i) = gask*.5*(po2(k,i)*7./32.+pn2(k,i)*7./28.+
     |                        po1(k,i)*5./16.+phe(k,i)*5./4.)
          fkm(k,i) = po2(k,i)*4.03 + pn2(k,i)*3.42 + po1(k,i)*3.9 +
     |               phe(k,i)*3.84
          fkt(k,i) = (po2(k,i)+pn2(k,i))*56. + po1(k,i)*75.9 + 
     |                phe(k,i)*299.
        enddo
      enddo

!     call addfld('fcp1',' ',' ',fcp(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('fkm1',' ',' ',fkm(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('fkt1',' ',' ',fkt(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev0,lev1-1
          fkm(k,i) = fkm(k,i)*(tn(k,i)+.5*(t0(k)+t0(k+1)))**0.69*1.e-6
          fkt(k,i) = fkt(k,i)*(tn(k,i)+.5*(t0(k)+t0(k+1)))**0.69
        enddo
      enddo
      do i=lon0,lon1
        fkm(lev1,i) = 1.e-6*fkm(lev1,i)*(tn(lev1-1,i)+
     |                1.5*t0(lev1)-.5*t0(lev1-1))**0.69
        fkt(lev1,i) =       fkt(lev1,i)*(tn(lev1-1,i)+
     |                1.5*t0(lev1)-.5*t0(lev1-1))**0.69
      enddo

!     call addfld('fcp2',' ',' ',fcp(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('fkm2',' ',' ',fkm(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('fkt2',' ',' ',fkt(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev1,lev0+1,-1
          fcp(k,i) = .5*(fcp(k,i)+fcp(k-1,i))
          fkm(k,i) = .5*(fkm(k,i)+fkm(k-1,i))
          fkt(k,i) = .5*(fkt(k,i)+fkt(k-1,i))
        enddo 
      enddo
      do i=lon0,lon1
        fcp(1,i) = 2.*fcp(1,i)-fcp(2,i)
        fkm(1,i) = 2.*fkm(1,i)-fkm(2,i)
        fkt(1,i) = 2.*fkt(1,i)-fkt(2,i)
      enddo

!     call addfld('CP',' ',' ',fcp(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('KT',' ',' ',fkt(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('KM',' ',' ',fkm(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!
! Save molecular viscosity diagnostic:
      call mkdiag_MU_M('MU_M',fkm(:,lon0:lon1),lev0,lev1,lon0,lon1,lat)

      end subroutine cpktkm

!
      subroutine cpktkm(tn,o2,o1,he,n2,mbar,fcp,fkt,fkm,
     |  lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Define diagnostics CP, KT, and KM.
!
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |                      rmassinv_he,gask
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
     |  he,  ! helium (mmr)
     |  n2,  ! molecular nitrogen (mmr)
     |  mbar ! mean molecular weight
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  fcp, ! specific heat at constant pressure (ergs/deg/gm)
     |  fkt, ! molecular diffusion (ergs/cm/deg/sec)
     |  fkm  ! molecular viscosity (gm/cm/sec)
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  t,po2,po1,phe,pn2,fcp1,fkm1,fkt1
!
!     call addfld('tn_cp',' ',' ',tn(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('o2_cp',' ',' ',o2(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('o1_cp',' ',' ',o1(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)

      t(lev0:lev1-1,:) = tn(lev0:lev1-1,lon0:lon1)
      t(lev1,:) = tn(lev1-1,lon0:lon1)

! Modified by EKS (Eric Sutton):
!  Compute po2, po1, phe, and pn2 as the number density mixing ratios.
!  Add the effect of Helium to the cp,km,kt calculations consistent with
!  Banks and Kockarts, 1973.
      do i=lon0,lon1
        do k=lev0,lev1
          po2(k,i) = o2(k,i)*rmassinv_o2
          po1(k,i) = o1(k,i)*rmassinv_o1
          phe(k,i) = he(k,i)*rmassinv_he
          pn2(k,i) = n2(k,i)*rmassinv_n2
          fcp1(k,i) = gask*.5*mbar(k,i)*
     |      (po2(k,i)*7.*rmassinv_o2+po1(k,i)*5.*rmassinv_o1+
     |       phe(k,i)*5.*rmassinv_he+pn2(k,i)*7.*rmassinv_n2)
          fkm1(k,i) = t(k,i)**0.69*1.e-6*mbar(k,i)*
     |      (po2(k,i)*4.03 + po1(k,i)*3.9 +
     |       phe(k,i)*3.84 + pn2(k,i)*3.42)
          fkt1(k,i) = t(k,i)**0.69*mbar(k,i)*
     |      (po2(k,i)*56.  + po1(k,i)*75.9 +
     |       phe(k,i)*299. + pn2(k,i)*56.)
        enddo
      enddo

!     call addfld('fcp1',' ',' ',fcp1(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('fkm1',' ',' ',fkm1(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('fkt1',' ',' ',fkt1(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev1,lev0+1,-1
          fcp(k,i) = .5*(fcp1(k,i)+fcp1(k-1,i))
          fkm(k,i) = .5*(fkm1(k,i)+fkm1(k-1,i))
          fkt(k,i) = .5*(fkt1(k,i)+fkt1(k-1,i))
        enddo 
      enddo
      do i=lon0,lon1
        fcp(1,i) = 1.5*fcp1(1,i)-.5*fcp1(2,i)
        fkm(1,i) = 1.5*fkm1(1,i)-.5*fkm1(2,i)
        fkt(1,i) = 1.5*fkt1(1,i)-.5*fkt1(2,i)
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

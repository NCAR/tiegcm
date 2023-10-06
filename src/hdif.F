! 
      module hdif_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlevp1,nlonp4,nlat
      use addfld_module,only: addfld
      use mpi_module,only: mpi_timing
      implicit none
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
      contains
!-----------------------------------------------------------------------
!
      subroutine hdif1(tn_nm,un_nm,vn_nm,mbar,fnrh,fkmh,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save global kmh (eddy viscosity) and nrh (M/T) for use in hdif2
!   and hdif3. This routine is called from advance.
!
      use cons_module,only: t0,cs,dlamda,dphi,re_inv
      use mpi_module,only: mp_bndlats_f3d,mp_bndlons_f3d
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn_nm, ! neutral temperature at time n-1
     |  un_nm, ! zonal velocity at time n-1
     |  vn_nm, ! meridional velocity at time n-1
     |  mbar   ! mean molecular weight
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(out)
     |  :: fnrh,fkmh  ! eddy viscosity, M/T
!
! Local:
      integer :: k,i,lat,nlevs,ier
      real :: abcsj,abcsjp,con1,con2,con3,delt,dels
      real :: cp2 = 0.2
!
#ifdef VT
!     code = 119 ; state = 'hdif1' ; activity='ModelCode'
      call vtbegin(119,ier)
#endif
!
! Latitude loop starts at lat0-2.
! kmh and nrh are defined from lat0-2 to lat1+2
! in which lat1+2 is from mp_bndlats
!
      do lat=lat0-2,lat1+1
!
! KMH = eddy viscosity = 2*K0*K0*SQRT(DS*DS+DT*DT)
! cs(lat) = cos(lat)
!
        abcsj = abs(cs(lat))
        abcsjp = abs(cs(lat+1))
        con1 = re_inv*.5/dlamda
        con2 = re_inv/(dphi*(abcsj+abcsjp))
        con3 = 2.*cp2*cp2
!
        do i=lon0-2,lon1+1
          do k=lev0,lev1-1

!           write(6,"('lat=',i2,' i=',i2,' k=',i2,
!    |        ' lat0,1=',2i3,' lon0,1=',2i3,' lev0,1=',2i3)")
!    |        lat,i,k,lat0,lat1,lon0,lon1,lev0,lev1

            delt = 
     |        con1*((un_nm(k,i+1,lat+1)-un_nm(k,i,lat+1))/cs(lat+1)+
     |              (un_nm(k,i+1,lat)-un_nm(k,i,lat))/cs(lat))-
     |        con2*((vn_nm(k,i+1,lat+1)+vn_nm(k,i,lat+1))*abcsjp -
     |              (vn_nm(k,i+1,lat)+vn_nm(k,i,lat))*abcsj)
            dels = 
     |        con1*((vn_nm(k,i+1,lat+1)-vn_nm(k,i,lat+1))/cs(lat+1)+
     |              (vn_nm(k,i+1,lat)-vn_nm(k,i,lat))/cs(lat))+
     |        con2*((un_nm(k,i+1,lat+1)+un_nm(k,i,lat+1))*abcsjp -
     |              (un_nm(k,i+1,lat)+un_nm(k,i,lat))*abcsj)
            fkmh(k,i,lat) = con3*sqrt(dels*dels+delt*delt)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+1
      enddo ! lat=lat0-2,lat1+1
!
! NRH = mbar(k) / (tnm(k)+(t0(k)+t0(k+1))*0.5)
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            fnrh(k,i,lat) = mbar(k,i,lat) /
     |                      (tn_nm(k,i,lat)+((t0(k)+t0(k+1))*.5))
          enddo
        enddo
      enddo
!
! Exchange boundary lats and lons in fkmh:
      nlevs = lev1-lev0+1
      call mp_bndlats_f3d(fkmh,nlevs,lon0,lon1,lat0,lat1,1)
      call mp_bndlons_f3d(fkmh,nlevs,lon0,lon1,lat0,lat1,1,0)
!
!     do lat=lat0,lat1
!         call addfld('FNRH',' ',' ',fnrh(lev0:lev1-1,lon0:lon1,lat),
!    |      'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!         call addfld('FKMH',' ',' ',fkmh(lev0:lev1-1,lon0:lon1,lat),
!    |      'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
#ifdef VT
!     code = 119 ; state = 'hdif1' ; activity='ModelCode'
      call vtend(119,ier)
#endif
      end subroutine hdif1
!-----------------------------------------------------------------------
!
      subroutine hdif2(tn_nm,un_nm,vn_nm,o2_nm,o1_nm,he_nm,fnrh,fkmh,
     |  fkldt,fkldu,fkldv,fkldo2,fkldo1,fkldhe,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save fkldx horizontal diffusion terms for hdif3. fkmh and fnrh are 
! module data, and were calculated by hdif1.
!
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn_nm,  ! neutral temperature at time n-1 (deg K)
     |  un_nm,  ! zonal velocity at time n-1 (cm/s)
     |  vn_nm,  ! meridional velocity at time n-1 (cm/s)
     |  o2_nm,  ! molecular oxygen at time n-1 (mmr)
     |  o1_nm,  ! atomic oxygen at time n-1 (mmr)
     |  he_nm,  ! helium at time n-1 (mmr)
     |  fnrh,   ! eddy viscosity
     |  fkmh    ! M/T
!
! Ouput args:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  fkldt,  ! tn horizontal diffusion term
     |  fkldu,  ! un horizontal diffusion term
     |  fkldv,  ! vn horizontal diffusion term
     |  fkldo1, ! o1 horizontal diffusion term
     |  fkldo2, ! o2 horizontal diffusion term
     |  fkldhe  ! He horizontal diffusion term
!
! Local:
      integer :: k,i,lat,ier
      real :: 
     |  avkmh (lev0:lev1,lon0:lon1), ! average kmh
     |  rhokmh(lev0:lev1,lon0:lon1)  ! nrh*avkmh
!
#ifdef VT
!     code = 120 ; state = 'hdif2' ; activity='ModelCode'
      call vtbegin(120,ier)
#endif
!
! Latitude scan:
      do lat=lat0-1,lat1+1
!       write(6,"('hdif2: lat=',i3)") lat
!
! avkmh = average fkmh:
! Allocation of fkmh in sub allocdata (fields_module.F):
!     allocate(fkmh(nlevp1,lon0-2:lon1+2,lat0-2:lat1+2),stat=ier)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          avkmh(k,i) = (fkmh(k,i-1,lat  )+fkmh(k,i,lat  )+
     |                  fkmh(k,i-1,lat+1)+fkmh(k,i,lat+1))*0.25
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('AVKMH',' ',' ',avkmh(lev0:lev1-1,lon0:lon1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! rhokmh = rho*avg(kmh)
      do i=lon0,lon1
        do k=lev0,lev1-1
          rhokmh(k,i) = avkmh(k,i)*fnrh(k,i,lat)
        enddo
        rhokmh(lev1,i) = 0. ! added to prevent NaNS init fpe
      enddo
!
!     if (lat > 0)
!    |  call addfld('RHOKMH',' ',' ',rhokmh(lev0:lev1-1,lon0:lon1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! CALC RHO*KMH*(L*L(D*D)(PSI)) AT J AND N-1
!
! Define fkldx at i=1,nlonp4 and lat.
! (note prognostic inputs at time n-1 have lon0-2->lon1+2 from bndlons 
!  call in advance)
!
! FOR PSI = U
      call lsqdsq(un_nm(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            un_nm(lev0:lev1,lon0-1:lon1+1,lat),
     |            un_nm(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            avkmh(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1
          fkldu(k,i,lat) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDU',' ',' ',fkldu(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = V
      call lsqdsq(vn_nm(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            vn_nm(lev0:lev1,lon0-1:lon1+1,lat),
     |            vn_nm(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            avkmh(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1
          fkldv(k,i,lat) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDV',' ',' ',fkldv(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = T
      call lsqdsq(tn_nm(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            tn_nm(lev0:lev1,lon0-1:lon1+1,lat),
     |            tn_nm(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            avkmh(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1
          fkldt(k,i,lat) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDT',' ',' ',fkldt(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = O2
      call lsqdsq(o2_nm(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            o2_nm(lev0:lev1,lon0-1:lon1+1,lat),
     |            o2_nm(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            avkmh(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1
          fkldo2(k,i,lat) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDO2',' ',' ',fkldo2(lev0:lev1-1,lon0:lon1,lat)
!    |    ,'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = O
      call lsqdsq(o1_nm(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            o1_nm(lev0:lev1,lon0-1:lon1+1,lat),
     |            o1_nm(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            avkmh(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1
          fkldo1(k,i,lat) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDO1',' ',' ',fkldo1(lev0:lev1-1,lon0:lon1,lat)
!    |    ,'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = He
      call lsqdsq(he_nm(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            he_nm(lev0:lev1,lon0-1:lon1+1,lat),
     |            he_nm(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            avkmh(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1
          fkldhe(k,i,lat) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDHE',' ',' ',fkldhe(lev0:lev1-1,lon0:lon1,lat)
!    |    ,'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      enddo ! lat=lat0-1,lat1+1
!
#ifdef VT
!     code = 120 ; state = 'hdif2' ; activity='ModelCode'
      call vtend(120,ier)
#endif
      end subroutine hdif2
!-----------------------------------------------------------------------
      subroutine hdif_bndlatlons(kldt,kldu,kldv,kldo2,kldo1,kldhe,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Exchange boundary latitudes/longitudes of horizontal diffusion coeffs
! calculated by hdif2. This is called from dynamics, in preparation 
! for hdif3.
!
      use mpi_module,only: mp_bndlats_f3d,mp_bndlons_f3d
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) :: kldt,kldu,kldv,kldo2,kldo1,kldhe
!
! Local:
      integer :: nlevs
      real :: f(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,6)

      nlevs = lev1-lev0+1

      f(:,:,:,1) = kldt
      f(:,:,:,2) = kldu
      f(:,:,:,3) = kldv
      f(:,:,:,4) = kldo2
      f(:,:,:,5) = kldo1
      f(:,:,:,6) = kldhe

      call mp_bndlats_f3d(f,nlevs,lon0,lon1,lat0,lat1,6)
      call mp_bndlons_f3d(f,nlevs,lon0,lon1,lat0,lat1,6,0)

      kldt = f(:,:,:,1)
      kldu = f(:,:,:,2)
      kldv = f(:,:,:,3)
      kldo2 = f(:,:,:,4)
      kldo1 = f(:,:,:,5)
      kldhe = f(:,:,:,6)

      end subroutine hdif_bndlatlons
!-----------------------------------------------------------------------
      subroutine hdif3(cp,fnrh,
     |  kldt,kldu,kldv,kldo2,kldo1,kldhe, ! input
     |  hdt ,hdu ,hdv ,hdo2 ,hdo1,hdhe,   ! output
     |  lev0,lev1,lon0,lon1,lat)
      use mpi_module,only: lat0,lat1
!
! Calculate horizontal diffusion terms for t,u,v,o2,o1,he at current latitude, 
! using coefficients that were output by hdif2. (hdif2 is called from advance, 
! hdif3 is called from dynamics). Sub kld_bndlons has been called prior to
! this routine so tasks have boundary longitudes for kldt,u,v,o2,o1,he.
!
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! 2d input:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: cp
!
! 3d input at full task subdomain:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  fnrh, ! eddy viscosity
     |  kldt,kldu,kldv,kldo2,kldo1,kldhe ! input from hdif2
!
! 2d output for current latitude:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out):: 
     |  hdt ,hdu ,hdv ,hdo2 ,hdo1, hdhe  ! output
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  fnrh_inv, hdout, cpi 
!
! fnrh is eddy viscosity (hdif module data), output by hdif1:
      do i=lon0,lon1
        do k=lev0,lev1-1
          fnrh_inv(k,i) = -1. / fnrh(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! Make hdu from kldu (hdu was nflh):
      call lsqdsq(kldu(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            kldu(lev0:lev1,lon0-1:lon1+1,lat  ),
     |            kldu(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            hdout(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdu(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdu(lev1,i) = 0.
      enddo ! i=lon0,lon1
!     call addfld('HDU',' ',' ',hdu(lev0:lev1-1,lon0:lon1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdv from kldv (hdv was nfph):
      call lsqdsq(kldv(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            kldv(lev0:lev1,lon0-1:lon1+1,lat  ),
     |            kldv(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            hdout(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdv(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdv(lev1,i) = 0.
      enddo ! i=lon0,lon1
!     call addfld('HDV',' ',' ',hdv(lev0:lev1-1,lon0:lon1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdt from kldt (was nqdh):
      call lsqdsq(kldt(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            kldt(lev0:lev1,lon0-1:lon1+1,lat  ),
     |            kldt(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            hdout(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1-1
          cpi(k,i) = .5*(cp(k,i)+cp(k+1,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdt(k,i) = hdout(k,i)*fnrh_inv(k,i)*cpi(k,i)
        enddo ! k=lev0,lev1-1
        hdt(lev1,i) = 0.
      enddo ! i=lon0,lon1
!     call addfld('HDT',' ',' ',hdt(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdo2 from kldo2 (was npsdh):
      call lsqdsq(kldo2(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            kldo2(lev0:lev1,lon0-1:lon1+1,lat  ),
     |            kldo2(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            hdout(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdo2(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdo2(lev1,i) = 0.
      enddo ! i=lon0,lon1
!     call addfld('HDO2',' ',' ',hdo2(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdo1 from kldo1 (was npsdh2):
      call lsqdsq(kldo1(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            kldo1(lev0:lev1,lon0-1:lon1+1,lat  ),
     |            kldo1(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            hdout(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdo1(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdo1(lev1,i) = 0.
      enddo ! i=lon0,lon1
!     call addfld('HDO1',' ',' ',hdo1(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdhe from kldhe:
      call lsqdsq(kldhe(lev0:lev1,lon0-1:lon1+1,lat+1),
     |            kldhe(lev0:lev1,lon0-1:lon1+1,lat  ),
     |            kldhe(lev0:lev1,lon0-1:lon1+1,lat-1),
     |            hdout(lev0:lev1,lon0:lon1),lon0,lon1,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1-1
          hdhe(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdhe(lev1,i) = 0.
      enddo ! i=lon0,lon1
!     call addfld('HDHE',' ',' ',hdhe(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      end subroutine hdif3
!-----------------------------------------------------------------------
      end module hdif_module

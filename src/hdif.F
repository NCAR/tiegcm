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
! An extra latitude is needed in fnrh and fkmh: 
! 10/17/02 bf: changed lon dim from nlonp4+1 to 0:nlonp4+1 
      real,dimension(nlevp1,0:nlonp4+1,-2:nlat) ::
     |  fnrh,  ! eddy viscosity
     |  fkmh   ! M/T
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
      subroutine hdif1(tn_nm,un_nm,vn_nm,barm,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Save global kmh (eddy viscosity) and nrh (M/T) for use in hdif2
!   and hdif3. This routine is called from advance.
!
      use cons_module,only: t0,cs,dlamda,dphi,re_inv
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn_nm, ! neutral temperature at time n-1
     |  un_nm, ! zonal velocity at time n-1
     |  vn_nm, ! meridional velocity at time n-1
     |  barm   ! mean molecular weight
!
! Local:
      integer :: k,ii,i,j,ip1,ip1k,lat,ier
      real :: fmin,fmax
      real :: abcsj,abcsjp,con1,con2,con3,delt,dels
      real :: cp2 = 0.2
!
#ifdef VT
!     code = 119 ; state = 'hdif1' ; activity='ModelCode'
      call vtbegin(119,ier)
#endif
!
! Latitude loop starts at lat0-3. Inputs are referenced only through
! lat+2, but kmh and nrh are defined from lat0-3 to lat1.
!
      do lat=lat0-3,lat1
!
! KMH = eddy viscosity = 2*K0*K0*SQRT(DS*DS+DT*DT)
! cs(lat) = cos(lat)
!
        abcsj = abs(cs(lat+1))
        abcsjp = abs(cs(lat+2))
        con1 = re_inv*.5/dlamda
        con2 = re_inv/(dphi*(abcsj+abcsjp))
        con3 = 2.*cp2*cp2
!
!  ii = 1 -> nlev*nlonp4
!   i = 1 -> 75,76
! ip1 = 2 -> 76, 1
!
        do i=lon0,lon1
          do k=lev0,lev1-1
            ip1 = i+1

!           if (i==nlonp4) ip1 = 1
!           if (ip1==1) then

! 10/17/02 bf: if i==nlonp4, then current task does not have i==1, so try
!   using ip1==nlonp4-3 instead (it should be same as i==1 because of
!   periodic points). 
!
            if (i==nlonp4) ip1 = 9999
            if (ip1==9999) then
              ip1k = k+1
              ip1 = nlonp4-3
            else
              ip1k = k
            endif

!           write(6,"('lat=',i2,' i=',i2,' k=',i2,' ip1k=',i2,' ip1=',
!    |        i2,' lat0,1=',2i3,' lon0,1=',2i3,' lev0,1=',2i3)")
!    |        lat,i,k,ip1k,ip1,lat0,lat1,lon0,lon1,lev0,lev1

            delt = 
     |        con1*((un_nm(ip1k,ip1,lat+2)-un_nm(k,i,lat+2))/cs(lat+2)+
     |              (un_nm(ip1k,ip1,lat+1)-un_nm(k,i,lat+1))/cs(lat+1))-
     |        con2*((vn_nm(ip1k,ip1,lat+2)+vn_nm(k,i,lat+2))*abcsjp -
     |              (vn_nm(ip1k,ip1,lat+1)+vn_nm(k,i,lat+1))*abcsj)
            dels = 
     |        con1*((vn_nm(ip1k,ip1,lat+2)-vn_nm(k,i,lat+2))/cs(lat+2)+
     |              (vn_nm(ip1k,ip1,lat+1)-vn_nm(k,i,lat+1))/cs(lat+1))+
     |        con2*((un_nm(ip1k,ip1,lat+2)+un_nm(k,i,lat+2))*abcsjp -
     |              (un_nm(ip1k,ip1,lat+1)+un_nm(k,i,lat+1))*abcsj)
            fkmh(k,i,lat) = con3*sqrt(dels*dels+delt*delt)
          enddo ! i=lon0,lon1
        enddo ! k=lev0,lev1-1
!
! NRH = ((barm(k)+barm(k+1))*0.5) / (tnm(k)+(t0(k)+t0(k+1))*0.5)
        do i=lon0,lon1
          do k=lev0,lev1-1
            fnrh(k,i,lat) = ((barm(k,i,lat+1)+barm(k+1,i,lat+1))*0.5) /
     |                      (tn_nm(k,i,lat+1)+((t0(k)+t0(k+1))*.5))
          enddo
        enddo
      enddo ! lat=lat0-3,lat1
#ifdef MPI
!
! Exchange boundary lats and lons in fkmh:
      call mp_bndlats_kmh(lev0,lev1,lon0,lon1,lat0,lat1)
      call mp_bndlons_kmh(lev0,lev1,lon0,lon1,lat0,lat1)
#endif
!     do j=lat0,lat1
!         call addfld('FNRH',' ',' ',fnrh(lev0:lev1-1,lon0:lon1,j),
!    |      'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!         call addfld('FKMH',' ',' ',fkmh(lev0:lev1-1,lon0:lon1,j),
!    |      'lev',lev0,lev1-1,'lon',lon0,lon1,j)
!     enddo ! j=lat0,lat1
!
#ifdef VT
!     code = 119 ; state = 'hdif1' ; activity='ModelCode'
      call vtend(119,ier)
#endif
      end subroutine hdif1
!-----------------------------------------------------------------------
!
      subroutine hdif2(tn_nm,un_nm,vn_nm,o2_nm,o1_nm,he_nm,
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
     |  he_nm   ! helium at time n-1 (mmr)
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
      integer :: k,i,ip1,ip1k,lonbeg,lonend,lat,ier
      real :: 
     |  avkmh (lev0:lev1,lon0:lon1), ! average kmh
     |  rhokmh(lev0:lev1,lon0:lon1)  ! nrh*avkmh
      real :: fmin,fmax
!
#ifdef VT
!     code = 120 ; state = 'hdif2' ; activity='ModelCode'
      call vtbegin(120,ier)
#endif
!
! Latitude scan:
      do lat=lat0-2,lat1
!       write(6,"('hdif2: lat=',i3)") lat
!
! avkmh = average fkmh:
! Allocation of fkmh in sub init_kmh_nrh (fields_module.F):
!     allocate(fkmh%data3d(nlevp1,nlonp4+1,-2:nlat),stat=ier)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          avkmh(k,i) = ((fkmh(k,i-1,lat-1)+fkmh(k,i,lat-1)+
     |                   fkmh(k,i-1,lat)) +fkmh(k,i,lat))*0.25
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
! CALC RHO*KMH*(L*L(D*D)(PSI)) AT J+1 AND N-1
!
! Define fkldx at i=2,nlonp4-1 and lat+1.
! (note pronostic inputs at time n-1 have lon0-2->lon1+2 from bndlons 
!  call in advance)
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 2
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-1
!
! FOR PSI = U
      call lsqdsq(un_nm(lev0:lev1,lonbeg-1:lonend+1,lat+2),
     |            un_nm(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            un_nm(lev0:lev1,lonbeg-1:lonend+1,lat),
     |            avkmh(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat+1)
      do i=lonbeg,lonend
        do k=lev0,lev1
          fkldu(k,i,lat+1) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDU',' ',' ',fkldu(lev0:lev1-1,lon0:lon1,lat+1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = V
      call lsqdsq(vn_nm(lev0:lev1,lonbeg-1:lonend+1,lat+2),
     |            vn_nm(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            vn_nm(lev0:lev1,lonbeg-1:lonend+1,lat),
     |            avkmh(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat+1)
      do i=lonbeg,lonend
        do k=lev0,lev1
          fkldv(k,i,lat+1) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDV',' ',' ',fkldv(lev0:lev1-1,lon0:lon1,lat+1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = T
      call lsqdsq(tn_nm(lev0:lev1,lonbeg-1:lonend+1,lat+2),
     |            tn_nm(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            tn_nm(lev0:lev1,lonbeg-1:lonend+1,lat),
     |            avkmh(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat+1)
      do i=lonbeg,lonend
        do k=lev0,lev1
          fkldt(k,i,lat+1) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDT',' ',' ',fkldt(lev0:lev1-1,lon0:lon1,lat+1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = O2
      call lsqdsq(o2_nm(lev0:lev1,lonbeg-1:lonend+1,lat+2),
     |            o2_nm(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            o2_nm(lev0:lev1,lonbeg-1:lonend+1,lat),
     |            avkmh(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat+1)
      do i=lonbeg,lonend
        do k=lev0,lev1
          fkldo2(k,i,lat+1) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDO2',' ',' ',fkldo2(lev0:lev1-1,lon0:lon1,lat+1)
!    |    ,'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = O
      call lsqdsq(o1_nm(lev0:lev1,lonbeg-1:lonend+1,lat+2),
     |            o1_nm(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            o1_nm(lev0:lev1,lonbeg-1:lonend+1,lat),
     |            avkmh(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat+1)
      do i=lonbeg,lonend
        do k=lev0,lev1
          fkldo1(k,i,lat+1) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDO1',' ',' ',fkldo1(lev0:lev1-1,lon0:lon1,lat+1)
!    |    ,'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! FOR PSI = He
      call lsqdsq(he_nm(lev0:lev1,lonbeg-1:lonend+1,lat+2),
     |            he_nm(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            he_nm(lev0:lev1,lonbeg-1:lonend+1,lat),
     |            avkmh(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat+1)
      do i=lonbeg,lonend
        do k=lev0,lev1
          fkldhe(k,i,lat+1) = avkmh(k,i)*rhokmh(k,i)
        enddo
      enddo
!     if (lat > 0)
!    |  call addfld('FKLDHE',' ',' ',fkldhe(lev0:lev1-1,lon0:lon1,lat+1)
!    |    ,'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      enddo ! lat=lat0-2,lat1
!
#ifdef VT
!     code = 120 ; state = 'hdif2' ; activity='ModelCode'
      call vtend(120,ier)
#endif
      end subroutine hdif2
!-----------------------------------------------------------------------
      subroutine hdif_bndlons(kldt,kldu,kldv,kldo2,kldo1,kldhe,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Exchange boundary longitudes of horizontal diffusion coeffs
! calculated by hdif2. This is called from dynamics, in preparation 
! for hdif3.
!
      use fields_module,only: f3d
#ifdef MPI
      use mpi_module,only: mp_bndlons_f3d
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) :: kldt,kldu,kldv,kldo2,kldo1,kldhe
!
! Local:
      integer :: i,nlevs
      real :: f(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,16)
!
#ifdef MPI
      nlevs = lev1-lev0+1
      do i=1,16
        f(:,:,:,i) = f3d(i)%data(:,:,:)
      enddo
      call mp_bndlons_f3d(f,nlevs,lon0,lon1,lat0,lat1,16,0)
      do i=1,16
        f3d(i)%data(:,:,:) = f(:,:,:,i)
      enddo
#endif
      end subroutine hdif_bndlons
!-----------------------------------------------------------------------
      subroutine hdif3(cp,
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
     |  kldt,kldu,kldv,kldo2,kldo1,kldhe ! input from hdif2
!
! 2d output for current latitude:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out):: 
     |  hdt ,hdu ,hdv ,hdo2 ,hdo1, hdhe  ! output
!
! Local:
      integer :: k,i,lonbeg,lonend
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  fnrh_inv, hdout, cpi 
!
! fnrh is eddy viscosity (hdif module data), output by hdif1:
      do i=lon0,lon1
        do k=lev0,lev1-1
          fnrh_inv(k,i) = -1. / fnrh(k,i,lat-1)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
! At 5 deg res (nlonp4=76), lonbeg,end will be 3,74
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
! Make hdu from kldu (hdu was nflh):
      call lsqdsq(kldu(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            kldu(lev0:lev1,lonbeg-1:lonend+1,lat  ),
     |            kldu(lev0:lev1,lonbeg-1:lonend+1,lat-1),
     |            hdout(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat)
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          hdu(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdu(lev1,i) = 0.
      enddo ! i=lonbeg,lonend
!     call addfld('HDU',' ',' ',hdu(lev0:lev1-1,lon0:lon1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdv from kldv (hdv was nfph):
      call lsqdsq(kldv(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            kldv(lev0:lev1,lonbeg-1:lonend+1,lat  ),
     |            kldv(lev0:lev1,lonbeg-1:lonend+1,lat-1),
     |            hdout(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat)
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          hdv(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdv(lev1,i) = 0.
      enddo ! i=lonbeg,lonend
!     call addfld('HDV',' ',' ',hdv(lev0:lev1-1,lon0:lon1),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdt from kldt (was nqdh):
      call lsqdsq(kldt(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            kldt(lev0:lev1,lonbeg-1:lonend+1,lat  ),
     |            kldt(lev0:lev1,lonbeg-1:lonend+1,lat-1),
     |            hdout(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat)
      do i=lon0,lon1
        do k=lev0,lev1-1
          cpi(k,i) = .5*(cp(k,i)+cp(k+1,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          hdt(k,i) = hdout(k,i)*fnrh_inv(k,i)*cpi(k,i)
        enddo ! k=lev0,lev1-1
        hdt(lev1,i) = 0.
      enddo ! i=lonbeg,lonend
!     call addfld('HDT',' ',' ',hdt(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdo2 from kldo2 (was npsdh):
      call lsqdsq(kldo2(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            kldo2(lev0:lev1,lonbeg-1:lonend+1,lat  ),
     |            kldo2(lev0:lev1,lonbeg-1:lonend+1,lat-1),
     |            hdout(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat)
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          hdo2(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdo2(lev1,i) = 0.
      enddo ! i=lonbeg,lonend
!     call addfld('HDO2',' ',' ',hdo2(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdo1 from kldo1 (was npsdh2):
      call lsqdsq(kldo1(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            kldo1(lev0:lev1,lonbeg-1:lonend+1,lat  ),
     |            kldo1(lev0:lev1,lonbeg-1:lonend+1,lat-1),
     |            hdout(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat)
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          hdo1(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdo1(lev1,i) = 0.
      enddo ! i=lonbeg,lonend
!     call addfld('HDO1',' ',' ',hdo1(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Make hdhe from kldhe:
      call lsqdsq(kldhe(lev0:lev1,lonbeg-1:lonend+1,lat+1),
     |            kldhe(lev0:lev1,lonbeg-1:lonend+1,lat  ),
     |            kldhe(lev0:lev1,lonbeg-1:lonend+1,lat-1),
     |            hdout(lev0:lev1,lonbeg:lonend),lonbeg,lonend,
     |            lev0,lev1,lat)
      do i=lonbeg,lonend
        do k=lev0,lev1-1
          hdhe(k,i) = hdout(k,i)*fnrh_inv(k,i)
        enddo ! k=lev0,lev1-1
        hdhe(lev1,i) = 0.
      enddo ! i=lonbeg,lonend
!     call addfld('HDHE',' ',' ',hdhe(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
      end subroutine hdif3
!-----------------------------------------------------------------------
      subroutine hdif_periodic(hdt,hdu,hdv,hdo2,hdo1,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
#ifdef MPI
      use mpi_module,only: mp_periodic_f3d
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) :: hdt,hdu,hdv,hdo2,hdo1
!
! Local:
      real :: ftmp(lev0:lev1,lon0:lon1,lat0:lat1,5)
#ifdef MPI
      ftmp(:,:,:,1) = hdt (:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,2) = hdu (:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,3) = hdv (:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,4) = hdo2(:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,5) = hdo1(:,lon0:lon1,lat0:lat1)

      call mp_periodic_f3d(ftmp,lev0,lev1,lon0,lon1,lat0,lat1,5)

      hdt (:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,1)
      hdu (:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,2)
      hdv (:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,3)
      hdo2(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,4)
      hdo1(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,5)
#endif
      end subroutine hdif_periodic
!-----------------------------------------------------------------------
#ifdef MPI
      subroutine mp_bndlats_kmh(lev0,lev1,lon0,lon1,lat0,lat1)
      use mpi_module,only: mytidi,mytidj,itask_table,handle_mpi_err,
     |  time_bndlats_kmh
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
#include <mpif.h>
!
! Local:
      integer ::
     |  irstat(MPI_STATUS_SIZE)      ! mpi receive status
      real,allocatable :: sndbuf(:,:,:),rcvbuf(:,:,:)
      integer :: ier,nlons,len,jprev,jnext,ireqrecv,ireqsend
      real :: starttime,endtime
#ifdef VT
!     call vtsymdef(111, 'mp_bndlats_kmh','Communication',ier)
      call vtbegin(111,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
!
      nlons = lon1-lon0+1
      allocate(sndbuf(nlevp1,nlons,3),stat=ier)
      if (ier /= 0) write(6,"('>>> mp_bndlats_kmh: error allocating',
     |  ' sndbuf: nlons=',i3,' nlevp1=',i3)") nlons,nlevp1
      allocate(rcvbuf(nlevp1,nlons,3),stat=ier)
      if (ier /= 0) write(6,"('>>> mp_bndlats_kmh: error allocating',
     |  ' rcvbuf: nlons=',i3,' nlevp1=',i3)") nlons,nlevp1
      len = nlons*nlevp1*3
!
      jprev = itask_table(mytidi,mytidj-1) ! task to south
      jnext = itask_table(mytidi,mytidj+1) ! task to north
!
! Send lat1-1,lat1 to jnext:
      sndbuf(:,:,:) = fkmh(:,lon0:lon1,lat1-2:lat1)
      call mpi_isend(sndbuf,len,MPI_REAL8,jnext,1,MPI_COMM_WORLD,
     |  ireqsend,ier)      
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlats_kmh send to jnext')
!
! Receive lat0-2,lat0-1 from jprev:
      call mpi_irecv(rcvbuf,len,MPI_REAL8,jprev,1,MPI_COMM_WORLD,
     |  ireqrecv,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlats_kmh recv from jprev')
!
! Wait for completions:
      call mpi_wait(ireqsend,irstat,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_bndlats_kmh wait for send')
      call mpi_wait(ireqrecv,irstat,ier)
      if (ier /= 0)
     |  call handle_mpi_err(ier,'mp_bndlats_kmh wait for recv')
!
! Define array from receive buffer:
      if (lat0 > 1) ! i.e., mytaskj > 0
     |  fkmh(:,lon0:lon1,lat0-3:lat0-1) = rcvbuf(:,:,:)
!
! Release local buffer space:
      deallocate(sndbuf)
      deallocate(rcvbuf)
!
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_bndlats_kmh=time_bndlats_kmh+(endtime-starttime)
      endif
#ifdef VT
!     call vtsymdef(111, 'mp_bndlats_kmh','Communication',ier)
      call vtend(111,ier)
#endif
      end subroutine mp_bndlats_kmh
!-----------------------------------------------------------------------
      subroutine mp_bndlons_kmh(lev0,lev1,lon0,lon1,lat0,lat1)
      use mpi_module,only: mytidi,mytidj,itask_table,ntaski,
     |  handle_mpi_err,time_bndlons_kmh
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
#include <mpif.h>
!
! Each task receives fkmh(:,lon0-1,lat0-3:lat1) from task iprev, 
!   and sends fkmh%data(:,lon1,lat0-3:lat1) to task inext.
! lon1==nlonp4 is sent to task mytidi==0, and lon0-1==0 is received 
!   from task ntaski-1
!
! Local:
      integer ::
     |  irstat(MPI_STATUS_SIZE)      ! mpi receive status
      integer :: j,jj,ier,len,iprev,inext,isend,irecv,nlats
      real,allocatable :: rcvbuf(:,:),sndbuf(:,:)
      real :: fmin,fmax
      real :: starttime,endtime
#ifdef VT
!     call vtsymdef(112, 'mp_bndlons_kmh','Communication',ier)
      call vtbegin(112,ier)
#endif
      if (mpi_timing) starttime = mpi_wtime()
!
! Allocate send and receive buffers:
      nlats = lat1-lat0+1+3
      allocate(rcvbuf(nlevp1,nlats),stat=ier)
      if (ier /= 0) write(6,"('>>> mp_bndlons_kmh: error allocating',
     |  ' rcvbuf: nlevp1=',i3,' nlats=',i3)") nlevp1,nlats
      allocate(sndbuf(nlevp1,nlats),stat=ier)
      if (ier /= 0) write(6,"('>>> mp_bndlons_kmh: error allocating',
     |  ' sndbuf: nlevp1=',i3,' nlats=',i3)") nlevp1,nlats
      len = nlevp1*nlats
!
! Receive lon0-1 from iprev:
! (if lon0==1, receive lon0-1 from task with mytidi==ntask-1)
      iprev = itask_table(mytidi-1,mytidj)
      if (lon0==1) iprev = itask_table(ntaski-1,mytidj) 
      call mpi_irecv(rcvbuf,len,MPI_REAL8,iprev,1,MPI_COMM_WORLD,
     |  irecv,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_kmh recv from iprev')
!
! Send lon1 to inext:
! (if lon1==nlonp4, send lon1 to task with mytidi==0)
      jj = 0
      do j=lat0-3,lat1
        jj = jj+1
        sndbuf(:,jj) = fkmh(:,lon1,j)
      enddo
      inext = itask_table(mytidi+1,mytidj)
      if (lon1==nlonp4) inext = itask_table(0,mytidj)

!     call fminmax(sndbuf(:,:),nlevp1*nlats,fmin,fmax)
!     write(6,"('Send lon1=',i2,' to inext=',i2,' sndbuf min,max=',
!    |  2e12.4)") lon1,inext,fmin,fmax

      call mpi_isend(sndbuf,len,MPI_REAL8,inext,1,MPI_COMM_WORLD,
     |  isend,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_kmh send to inext')
!
! Wait for completions:
      call mpi_wait(isend,irstat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_kmh wait for send')
      call mpi_wait(irecv,irstat,ier)
      if (ier /= 0) 
     |  call handle_mpi_err(ier,'mp_bndlons_kmh wait for recv')
!
! Copy lon0-1 from rcvbuf:
      call fminmax(rcvbuf(:,:),nlevp1*nlats,fmin,fmax)
!     write(6,"('Recv lon0-1=',i2,' from iprev=',i2,' rcvbuf min,max=',
!    |  2e12.4)") lon0-1,iprev,fmin,fmax

      jj = 0
      do j=lat0-3,lat1
        jj = jj+1
        fkmh(:,lon0-1,j) = rcvbuf(:,jj)
      enddo
!
! Release local buffer space:
      deallocate(sndbuf)
      deallocate(rcvbuf)
!
      if (mpi_timing) then
        endtime = mpi_wtime()
        time_bndlons_kmh=time_bndlons_kmh+(endtime-starttime)
      endif
#ifdef VT
!     call vtsymdef(112, 'mp_bndlons_kmh','Communication',ier)
      call vtend(112,ier)
#endif
      end subroutine mp_bndlons_kmh
!
! #endif MPI:
#endif
      end module hdif_module

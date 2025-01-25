      module hdif_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! tiegcmlicense.txt.
!
      implicit none
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! 2025/01/24 Haonan Wu: The previous version doesn't fill halo points.
! Reorganized so that halo points are filled before horizontal derivatives.
!
! Reference:
!   Washington and Williamson (1977),
!   A Description of the NCAR Global Circulation Models,
!     III. Numerical Approximations
!       C. Discrete Equations
!         3. Physical Processes
!           a. Subgrid Scale Diffusion
!
      contains
!-----------------------------------------------------------------------
      subroutine hdif(cp,mbar,tn_nm,un_nm,vn_nm,o2_nm,o1_nm,he_nm,
     |  hdt,hdu,hdv,hdo2,hdo1,hdhe,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate horizontal diffusion terms for t,u,v,o2,o1,he.
!
      use params_module,only: nlat,nlatp1,nlatp2
      use cons_module,only: cs,dlamda,dphi,re
      use mpi_module,only: mp_polelats_f3d,mp_bndlats_f3d,mp_bndlons_f3d
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)
     |  :: cp, ! specific heat
     |  mbar,  ! mean molecular weight
     |  tn_nm, ! neutral temperature at time n-1 (deg K)
     |  un_nm, ! zonal velocity at time n-1 (cm/s)
     |  vn_nm, ! meridional velocity at time n-1 (cm/s)
     |  o2_nm, ! molecular oxygen at time n-1 (mmr)
     |  o1_nm, ! atomic oxygen at time n-1 (mmr)
     |  he_nm  ! helium at time n-1 (mmr)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(out)
     |  :: hdt, hdu, hdv, hdo2, hdo1, hdhe  ! output
!
! Local:
      real,parameter :: cp2 = 0.2
      integer :: k,i,lat,n,nlevs,ier
      real :: con1,con2,con3
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  avkmh, ! average kmh
     |  nrh,   ! eddy viscosity, M/T
     |  cpi
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1,6) ::
     |  kldout,hdout
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  du,dv,um,vm,dup,dvp,ump,vmp,delt,dels,kmh,
     |  kldt,  ! tn horizontal diffusion term
     |  kldu,  ! un horizontal diffusion term
     |  kldv,  ! vn horizontal diffusion term
     |  kldo1, ! o1 horizontal diffusion term
     |  kldo2, ! o2 horizontal diffusion term
     |  kldhe  ! He horizontal diffusion term
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,2) :: f2
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,4) :: f4
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,6) :: f6
!
#ifdef VT
!     code = 119 ; state = 'hdif' ; activity='ModelCode'
      call vtbegin(119,ier)
#endif
!
      nlevs = lev1-lev0+1
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+1
          do k=lev0,lev1-1
            f4(k,i,lat,1) = (un_nm(k,i+1,lat)-un_nm(k,i,lat))/cs(lat)
            f4(k,i,lat,2) = (vn_nm(k,i+1,lat)-vn_nm(k,i,lat))/cs(lat)
            f4(k,i,lat,3) =
     |        (un_nm(k,i+1,lat)+un_nm(k,i,lat))*abs(cs(lat))
            f4(k,i,lat,4) =
     |        (vn_nm(k,i+1,lat)+vn_nm(k,i,lat))*abs(cs(lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+1
      enddo ! lat=lat0-2,lat1+2
!
! Exchange boundary lons in f4:
      call mp_bndlons_f3d(f4,nlevs,lon0,lon1,lat0,lat1,4,0)
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            du(k,i,lat) = f4(k,i,lat,1)
            dv(k,i,lat) = f4(k,i,lat,2)
            um(k,i,lat) = f4(k,i,lat,3)
            vm(k,i,lat) = f4(k,i,lat,4)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      do lat=lat0-2,lat1+1
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            f4(k,i,lat,1) = du(k,i,lat+1)
            f4(k,i,lat,2) = dv(k,i,lat+1)
            f4(k,i,lat,3) = um(k,i,lat+1)
            f4(k,i,lat,4) = vm(k,i,lat+1)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+1
!
! Exchange boundary lats in f4:
      call mp_bndlats_f3d(f4,nlevs,lon0,lon1,lat0,lat1,4)
!
! Extrapolate north pole:
      if (lat1 == nlat) then
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            do n=1,4
              f4(k,i,nlatp2,n) = 2*f4(k,i,nlatp1,n)-f4(k,i,nlat,n)
            enddo ! n=1,4
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      endif
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            dup(k,i,lat) = f4(k,i,lat,1)
            dvp(k,i,lat) = f4(k,i,lat,2)
            ump(k,i,lat) = f4(k,i,lat,3)
            vmp(k,i,lat) = f4(k,i,lat,4)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      con1 = 1/(2*re*dlamda)
      do lat=lat0-2,lat1+1
        con2 = 1/(re*dphi*(abs(cs(lat))+abs(cs(lat+1))))
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            f2(k,i,lat,1) =
     |        con1*(dup(k,i,lat)+du(k,i,lat))-
     |        con2*(vmp(k,i,lat)-vm(k,i,lat))
            f2(k,i,lat,2) =
     |        con1*(dvp(k,i,lat)+dv(k,i,lat))+
     |        con2*(ump(k,i,lat)-um(k,i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+1
!
! Exchange boundary lats in f2:
      call mp_bndlats_f3d(f2,nlevs,lon0,lon1,lat0,lat1,2)
!
! Extrapolate north pole
      if (lat1 == nlat) then
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            do n=1,2
              f2(k,i,nlatp2,n) = 2*f2(k,i,nlatp1,n)-f2(k,i,nlat,n)
            enddo ! n=1,2
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      endif
!
! KMH = eddy viscosity = 2*K0*K0*SQRT(DS*DS+DT*DT)
      con3 = 2*cp2**2
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            delt(k,i,lat) = f2(k,i,lat,1)
            dels(k,i,lat) = f2(k,i,lat,2)
            kmh(k,i,lat) = con3*sqrt(dels(k,i,lat)**2+delt(k,i,lat)**2)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
!
! avkmh = average kmh:
            avkmh(k,i,lat) = (kmh(k,i-1,lat  )+kmh(k,i,lat  )+
     |                        kmh(k,i-1,lat+1)+kmh(k,i,lat+1))*0.25
!
! NRH = mbar(k) / tnm(k)
            nrh(k,i,lat) = mbar(k,i,lat)/tn_nm(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! CALC RHO*KMH*(L*L(D*D)(PSI)) AT J AND N-1
!
! Define kldx at i=1,nlonp4 and lat.
! (note prognostic inputs at time n-1 have lon0-2->lon1+2 from bndlons
!  call in advance)
!
! FOR PSI = T,U,V,O2,O,He
      call lsqdsq(tn_nm,kldout(:,:,:,1),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(un_nm,kldout(:,:,:,2),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(vn_nm,kldout(:,:,:,3),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(o2_nm,kldout(:,:,:,4),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(o1_nm,kldout(:,:,:,5),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(he_nm,kldout(:,:,:,6),lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            do n=1,6
              f6(k,i,lat,n) = kldout(k,i,lat,n)*
     |          avkmh(k,i,lat)*nrh(k,i,lat)
            enddo ! n=1,6
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Exchange boundary latitudes/longitudes of horizontal diffusion coeffs.
      call mp_polelats_f3d(f6(:,lon0:lon1,:,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,6,(/1.,1.,1.,1.,1.,1./))
      call mp_bndlats_f3d(f6,nlevs,lon0,lon1,lat0,lat1,6)
      call mp_bndlons_f3d(f6,nlevs,lon0,lon1,lat0,lat1,6,0)
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            kldt(k,i,lat) = f6(k,i,lat,1)
            kldu(k,i,lat) = f6(k,i,lat,2)
            kldv(k,i,lat) = f6(k,i,lat,3)
            kldo2(k,i,lat) = f6(k,i,lat,4)
            kldo1(k,i,lat) = f6(k,i,lat,5)
            kldhe(k,i,lat) = f6(k,i,lat,6)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! Calculate horizontal diffusion terms for t,u,v,o2,o1,he.
!
! Make hdx from kldx:
      call lsqdsq(kldt,hdout(:,:,:,1),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(kldu,hdout(:,:,:,2),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(kldv,hdout(:,:,:,3),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(kldo2,hdout(:,:,:,4),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(kldo1,hdout(:,:,:,5),lev0,lev1,lon0,lon1,lat0,lat1)
      call lsqdsq(kldhe,hdout(:,:,:,6),lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            cpi(k,i,lat) = .5*(cp(k,i,lat)+cp(k+1,i,lat))
            hdt(k,i,lat) = -hdout(k,i,lat,1)/nrh(k,i,lat)*cpi(k,i,lat)
            hdu(k,i,lat) = -hdout(k,i,lat,2)/nrh(k,i,lat)
            hdv(k,i,lat) = -hdout(k,i,lat,3)/nrh(k,i,lat)
            hdo2(k,i,lat) = -hdout(k,i,lat,4)/nrh(k,i,lat)
            hdo1(k,i,lat) = -hdout(k,i,lat,5)/nrh(k,i,lat)
            hdhe(k,i,lat) = -hdout(k,i,lat,6)/nrh(k,i,lat)
          enddo ! k=lev0,lev1-1
          hdt(lev1,i,lat) = 0.
          hdu(lev1,i,lat) = 0.
          hdv(lev1,i,lat) = 0.
          hdo2(lev1,i,lat) = 0.
          hdo1(lev1,i,lat) = 0.
          hdhe(lev1,i,lat) = 0.
        enddo ! i=lon0,lon1
!
!       call addfld('DELT',' ',' ',delt(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('DELS',' ',' ',dels(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('KMH',' ',' ',kmh(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('AVKMH',' ',' ',avkmh(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('NRH',' ',' ',nrh(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('KLDT',' ',' ',kldt(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('KLDU',' ',' ',kldu(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('KLDV',' ',' ',kldv(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('KLDO2',' ',' ',kldo2(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('KLDO1',' ',' ',kldo1(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('KLDHE',' ',' ',kldhe(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('HDT',' ',' ',hdt(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDU',' ',' ',hdu(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDV',' ',' ',hdv(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDO2',' ',' ',hdo2(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDO1',' ',' ',hdo1(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDHE',' ',' ',hdhe(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
#ifdef VT
!     code = 119 ; state = 'hdif' ; activity='ModelCode'
      call vtend(119,ier)
#endif
      end subroutine hdif
!-----------------------------------------------------------------------
      subroutine lsqdsq(fj,fout,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Subroutine to calculate L**2DEL**2(PSI)
!
      use cons_module,only: cs
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)
     |  :: fj
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: fout
!
! Local:
      real :: dlam,dpsi,dpsip,dpsim
      integer :: k,i,lat
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
!
! fout = PSI(I+1)+PSI(I-1)-2*PSI(I)
            dlam = fj(k,i+1,lat)+fj(k,i-1,lat)-2*fj(k,i,lat)
!
! DPSIP=(PSI(J+1)-PSI(J))*.5*(CS(J+1)+CS(J))
            dpsip = (fj(k,i,lat+1)-fj(k,i,lat))*
     |        (0.5*(cs(lat+1)+cs(lat)))
!
! DPSIM=(PSI(J)-PSI(J-1))*.5*(CS(J)+CS(J-1))
            dpsim = (fj(k,i,lat)-fj(k,i,lat-1))*
     |        (0.5*(cs(lat)+cs(lat-1)))
!
! DPSI = (DPSIP-DPSIM)/CS(J)
            dpsi = (dpsip-dpsim)/cs(lat)
!
! S = DLAM + DPSI = (L**2)(D**2)PSI
            fout(k,i,lat) = dlam + dpsi
          enddo
        enddo
      enddo
!
      end subroutine lsqdsq
!-----------------------------------------------------------------------
      end module hdif_module

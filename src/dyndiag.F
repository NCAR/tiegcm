      module dyndiag_module
      use params_module,only: nlonp1,nlonp4,nlat,nlevp1,nlatp1
      use addfld_module,only: addfld
      implicit none
!
      real, dimension(:,:),allocatable ::
     |  sigp,    ! height integrated pedersen conductivity
     |  sigh,    ! height integrated hall conductivity
     |  qwind,   ! joule heating associated with neutral winds
     |  qamie,   ! joule heating from amie
     |  wtot,    ! total electromagnetic power
     |  work,    ! mechanical work
     |  fwindu,fwindv, ! horizontal current
     |  famieu,famiev, ! horizontal current from amie
     |  tec      ! total electron content
      real, dimension(:,:,:),allocatable ::
     |  qwind_sec,qamie_sec,work_sec,wtot_sec,fwindu_sec,
     |  fwindv_sec,famieu_sec,famiev_sec,
     |  tec_sec,ped_sec,hall_sec,qjoule_sec
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_dyndiag(lond0,lond1,latd0,latd1)

      integer,intent(in) :: lond0,lond1,latd0,latd1

      allocate(sigp(lond0:lond1,latd0:latd1))
      allocate(sigh(lond0:lond1,latd0:latd1))
      allocate(qwind(lond0:lond1,latd0:latd1))
      allocate(qamie(lond0:lond1,latd0:latd1))
      allocate(wtot(lond0:lond1,latd0:latd1))
      allocate(work(lond0:lond1,latd0:latd1))
      allocate(fwindu(lond0:lond1,latd0:latd1))
      allocate(fwindv(lond0:lond1,latd0:latd1))
      allocate(famieu(lond0:lond1,latd0:latd1))
      allocate(famiev(lond0:lond1,latd0:latd1))
      allocate(tec(lond0:lond1,latd0:latd1))
      allocate(qwind_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(qamie_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(work_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(wtot_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(fwindu_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(fwindv_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(famieu_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(famiev_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(tec_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(ped_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(hall_sec(lond0:lond1,latd0:latd1,nlevp1))
      allocate(qjoule_sec(lond0:lond1,latd0:latd1,nlevp1))

      end subroutine alloc_dyndiag
!-----------------------------------------------------------------------
      subroutine dyndiag(sigma1,sigma2,qji_tn,z,un,vn,ui,vi,
     |  lev0,lev1,lon0,lon1,lat)
      use magfield_module,only: dipmag,sn2dec,csdec,sndec,bmod
!
! Calculate dynamo related diagnostics for secondary histories.
! (this sub called from inside lat loop in dynamics, after lamdas)
!
! Args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  sigma1,  ! pederson conductivity (from lamdas)
     |  sigma2,  ! hall conductivity (from lamdas)
     |  qji_tn,  ! ion Joule heating for Tn (from qjoule)
     |  z,       ! geopotential from addiag
     |  un,vn,   ! neutral velocities
     |  ui,vi    ! ion velocities
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Local:
      integer :: i,k
      real,dimension(nlonp4) :: tm1,tm2,tm3,sndip 
      real :: qamie0,qwind0,wtot0,work0,famieu0,famiev0,
     |        fwindu0,fwindv0
      real :: bmod_lat(nlonp4) ! local bmod at current lat
!
! Get height integral of Ped, Hall conductance
      sigp(:,lat) = 0.
      sigh(:,lat) = 0.

      do i=lon0,lon1
        do k=lev0,lev1-1
!
!  Times ht integ by 1.e-2 because appears off for siemens (1.e-2cm/m?)
          sigp(i,lat) = sigp(i,lat) +
     |      (z(k+1,i)-z(k,i)) * sigma1(k,i) * 1.e-2
          sigh(i,lat) = sigh(i,lat) +
     |      (z(k+1,i)-z(k,i)) * sigma2(k,i) * 1.e-2
          ped_sec(i,lat,k) = sigma1(k,i)
          hall_sec(i,lat,k) = sigma2(k,i)
          qjoule_sec(i,lat,k) = qji_tn(k,i)
        enddo
        ped_sec(i,lat,nlevp1) = sigp(i,lat)
        hall_sec(i,lat,nlevp1) = sigh(i,lat)
      enddo

      call addfld('SIGMA1',' ',' ',sigma1(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      call addfld('SIGMA2',' ',' ',sigma2(lev0:lev1-1,lon0:lon1),
     |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

!
! Calculate the Joule heating due to wind and ion drift
!  qwind = Joule Heting due to wind
!  qamie = Joule Heating due to ion drift
!  work = the mechanical work
!  wtot = total electromagnetic power
!
      do i=lon0,lon1
        sndip(i) = sin(dipmag(i,lat))
!       if (abs(sndip(i)) .lt. 1.e-3)
!    |     print *, 'sndip is zero at i= ',i,sndip(i)
!
! term1=(bx**2 + bz**2)
! term2=(by**2 + bz**2)
! term3=bxby
!
        tm1(i) = sn2dec(i,lat)+(1.-sn2dec(i,lat))*sndip(i)**2
        tm2(i) = (1.-sn2dec(i,lat))+sn2dec(i,lat)*sndip(i)**2
        tm3(i) = sndec(i,lat)*csdec(i,lat) *(1.-sndip(i)**2)
      enddo
!
! Init:
      qwind(:,lat) = 0.
      qamie(:,lat) = 0.
      wtot(:,lat) = 0.
      work(:,lat) = 0.
      fwindu(:,lat) = 0.
      fwindv(:,lat) = 0.
      famieu(:,lat) = 0.
      famiev(:,lat) = 0.
!
! Set local global bmod_lat(nlonp4) (at current latitude) to avoid 
! out-of-bounds in bmod(nlonp1) in longitude dimension:
      bmod_lat(3:nlonp1+1) = bmod(1:nlonp1-1,lat) ! at 5 deg: 3-74 <- 1-72
      bmod_lat(1:2) = bmod_lat(nlonp1+1:nlonp1+2) ! periodic
      bmod_lat(nlonp4-1:nlonp4) = bmod_lat(3:4)   ! periodic
!
! Longitude, level loops:
      do k=lev0,lev1-1
        do i=lon0,lon1
!
         sndip(i) = sin(dipmag(i,lat))
         if (abs(sndip(i)) .gt. 0.1) then
!
! qwind is the Joule heating associated with nuetral wind terms
! qamie=sigp{(Bx^2+Bz^2)*Wx^2 + (By^2+Bz^2)*Wy^2-2.*Bx*By*Wx*Wy
!        -2.*B^2*Vx*Wx - 2.*B^2*Vy*Wy}
!
          qwind0 = sigma1(k,i)*bmod_lat(i)**2*(tm2(i)*
     |      un(k,i)**2+tm1(i)*vn(k,i)**2-2.*tm3(i)*un(k,i)*
     |      vn(k,i)-2.*(un(k,i)*ui(k,i)+vn(k,i)*vi(k,i)))
!
! qamie is Joule heating without neutral wind
! qamie=sigp{(Bx^2+Bz^2)*Vx^2+(By^2+Bz^2)*Vy^2+2.*Bx*By*Vx*Vy}/
!       sin(dip)^2
!
          qamie0 = sigma1(k,i)*(tm1(i)*ui(k,i)**2+tm2(i)*
     |      vi(k,i)**2+2.*tm3(i)*ui(k,i)*vi(k,i))
          qamie0 = qamie0*bmod_lat(i)**2 / sndip(i)**2
!
! wtot = total electric power
! wtot = sigp (E + UxB).E = Qamie + sigp*{(UxB).E}
!
          wtot0 = - sigma1(k,i)*bmod_lat(i)**2*(un(k,i)*ui(k,i)+
     |      vn(k,i)*vi(k,i)) + sigma2(k,i)*bmod_lat(i)**2*
     |      (tm3(i)*(un(k,i)*ui(k,i)-vn(k,i)*vi(k,i))-tm1(i)*
     |      vn(k,i)*ui(k,i)+tm2(i)*un(k,i)*vi(k,i))/sndip(i)
          wtot0 = qamie0 + wtot0
C Calculate the mechanical work done by wind
               work0 = wtot0 - qwind0 - qamie0
 
               qwind_sec(i,lat,k) = qwind0 * 1.e-9
               qamie_sec(i,lat,k) = qamie0 * 1.e-9
               work_sec(i,lat,k) = work0 * 1.e-9
               wtot_sec(i,lat,k) = wtot0 * 1.e-9

!
! To convert to mW/m^2 by mutipling by e-11 (since sigp (mho/m),height(cm)
!  B(in Gauss), and velocity (cm/s)
!
          qwind(i,lat) = qwind(i,lat)+(z(k+1,i)-z(k,i))*qwind0*1.e-11
          qamie(i,lat) = qamie(i,lat)+(z(k+1,i)-z(k,i))*qamie0*1.e-11
          wtot(i,lat) = wtot(i,lat)  +(z(k+1,i)-z(k,i))*wtot0*1.e-11
          work(i,lat) = work(i,lat)  +(z(k+1,i)-z(k,i))*work0*1.e-11
!
! Calculate horizontal currents:
          fwindu0 = -tm1(i)*vn(k,i)*sigma1(k,i)+tm3(i)*un(k,i)*
     |      sigma1(k,i) + sndip(i)*un(k,i)*sigma2(k,i)
          fwindv0 = tm2(i)*un(k,i)*sigma1(k,i)+sndip(i)*vn(k,i)*
     |      sigma2(k,i) - tm3(i)*vn(k,i)*sigma1(k,i)
!
! Convert wind velocity from cm/s to m/s by multipling 1.e-2 and
! Convert B from Gauss to T multipling 1.e-4 and 1.e-2 from cm to m
!
          fwindu0 = fwindu0*bmod_lat(i)*1.e-8/sndip(i)
          fwindv0 = fwindv0*bmod_lat(i)*1.e-8/sndip(i)
          famieu0 = -tm1(i)*sigma2(k,i)*
     |       ui(k,i) + sndip(i)*sigma1(k,i)*vi(k,i)
     |       - tm3(i)*sigma2(k,i)*vi(k,i)
          famiev0 = -tm2(i)*sigma2(k,i)*
     |       vi(k,i) - sndip(i)*sigma1(k,i)*ui(k,i)
     |       - tm3(i)*sigma2(k,i)*ui(k,i)
!
! Convert wind velocity from cm/s to m/s by multipling 1.e-2 and
! Convert B from Gauss to T multipling 1.e-4 and 1.e-2 from cm to m
!
          famieu0 = famieu0*bmod_lat(i)*1.e-8/sndip(i)**2
          famiev0 = famiev0*bmod_lat(i)*1.e-8/sndip(i)**2
          fwindu(i,lat) = fwindu(i,lat)+(z(k+1,i)-z(k,i))*fwindu0
          fwindv(i,lat) = fwindv(i,lat)+(z(k+1,i)-z(k,i))*fwindv0
          famieu(i,lat) = famieu(i,lat)+(z(k+1,i)-z(k,i))*famieu0
          famiev(i,lat) = famiev(i,lat)+(z(k+1,i)-z(k,i))*famiev0

          fwindu_sec(i,lat,k) = fwindu0
          fwindv_sec(i,lat,k) = fwindv0
          famieu_sec(i,lat,k) = famieu0
          famiev_sec(i,lat,k) = famiev0
!
         endif
!
        enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1-1

      qwind_sec(:,lat,nlevp1) = qwind(:,lat)
      qamie_sec(:,lat,nlevp1) = qamie(:,lat)
      work_sec(:,lat,nlevp1) = work(:,lat)
      wtot_sec(:,lat,nlevp1) = wtot(:,lat)
      fwindu_sec(:,lat,nlevp1) = fwindu(:,lat)
      fwindv_sec(:,lat,nlevp1) = fwindv(:,lat)
      famieu_sec(:,lat,nlevp1) = famieu(:,lat)
      famiev_sec(:,lat,nlevp1) = famiev(:,lat)
!     if (lat==25) write(6,"('dyndiag>>> ped_sec,',
!    |  'qamie_sec = ',/,(6g12.2))")
!    |   ped_sec(lon0,lat,:),qamie_sec(lon0,lat,:)
      call addfld('QAMIE1',' ',' ',qamie_sec(lon0:lon1,lat,:),
     |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
      call addfld('QWIND1',' ',' ',qwind_sec(lon0:lon1,lat,:),
     |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
      end subroutine dyndiag
!-----------------------------------------------------------------------
      subroutine dyndiag_sech(lon0,lon1,lev0,lev1,lat)
!
! Save 2d (lon,lat) diagnostics to secondary histories:
! This is called from dynamics after lamdas lat loop.
!
! Note total fwind and famie are not calculated here (as they were in
!   tgcm15) because of the need for lat-1,lat+1,lon-1,lon+1. These can 
!   be calculated by a post-processor after reading the u and v components 
!   from the secondary histories.
!
      integer,intent(in) :: lon0,lon1,lev0,lev1,lat

!     call addfld('PEDERSEN',' ',' ',ped_sec(lon0:lon1,lat,:),
!    |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
!     call addfld('HALL'    ,' ',' ',hall_sec(lon0:lon1,lat,:),
!    |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
      call addfld('QWIND'   ,' ',' ',qwind_sec(lon0:lon1,lat,:),
     |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
      call addfld('QAMIE'   ,' ',' ',qamie_sec(lon0:lon1,lat,:),
     |  'lon',lon0,lon1,'lev',lev0,lev1,lat)

!     call addfld('FWINDU'   ,' ',' ',fwindu(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
!     call addfld('FWINDV'   ,' ',' ',fwindv(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
!     call addfld('FAMIEU'   ,' ',' ',famieu(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
!     call addfld('FAMIEV'   ,' ',' ',famiev(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lev',lev0,lev1,lat)
      end subroutine dyndiag_sech
!-----------------------------------------------------------------------
      subroutine dyndiag_bndlatlon

      use mpi_module,only: lon0,lon1,lat0,lat1,
     |  mp_bndlats_f2d,mp_bndlons_f2d,
     |  mp_bndlats_f3d,mp_bndlons_f3d

      integer,parameter :: nf2d = 11, nf3d = 12
      integer :: k
      real,dimension(lon0-2:lon1+2,lat0-2:lat1+2,nf2d) :: f2d
      real,dimension(nlevp1,lon0-2:lon1+2,lat0-2:lat1+2,nf3d) :: f3d

      f2d(:,:,1) = sigp
      f2d(:,:,2) = sigh
      f2d(:,:,3) = qwind
      f2d(:,:,4) = qamie
      f2d(:,:,5) = wtot
      f2d(:,:,6) = work
      f2d(:,:,7) = fwindu
      f2d(:,:,8) = fwindv
      f2d(:,:,9) = famieu
      f2d(:,:,10) = famiev
      f2d(:,:,11) = tec

      call mp_bndlats_f2d(f2d,lon0,lon1,lat0,lat1,nf2d)
      call mp_bndlons_f2d(f2d,lon0,lon1,lat0,lat1,nf2d)

      sigp = f2d(:,:,1)
      sigh = f2d(:,:,2)
      qwind = f2d(:,:,3)
      qamie = f2d(:,:,4)
      wtot = f2d(:,:,5)
      work = f2d(:,:,6)
      fwindu = f2d(:,:,7)
      fwindv = f2d(:,:,8)
      famieu = f2d(:,:,9)
      famiev = f2d(:,:,10)
      tec = f2d(:,:,11)

      do k = 1,nlevp1
        f3d(k,:,:,1) = qwind_sec(:,:,k)
        f3d(k,:,:,2) = qamie_sec(:,:,k)
        f3d(k,:,:,3) = work_sec(:,:,k)
        f3d(k,:,:,4) = wtot_sec(:,:,k)
        f3d(k,:,:,5) = fwindu_sec(:,:,k)
        f3d(k,:,:,6) = fwindv_sec(:,:,k)
        f3d(k,:,:,7) = famieu_sec(:,:,k)
        f3d(k,:,:,8) = famiev_sec(:,:,k)
        f3d(k,:,:,9) = tec_sec(:,:,k)
        f3d(k,:,:,10) = ped_sec(:,:,k)
        f3d(k,:,:,11) = hall_sec(:,:,k)
        f3d(k,:,:,12) = qjoule_sec(:,:,k)
      enddo

      call mp_bndlats_f3d(f3d,nlevp1,lon0,lon1,lat0,lat1,nf3d)
      call mp_bndlons_f3d(f3d,nlevp1,lon0,lon1,lat0,lat1,nf3d,0)

      do k = 1,nlevp1
        qwind_sec(:,:,k) = f3d(k,:,:,1)
        qamie_sec(:,:,k) = f3d(k,:,:,2)
        work_sec(:,:,k) = f3d(k,:,:,3)
        wtot_sec(:,:,k) = f3d(k,:,:,4)
        fwindu_sec(:,:,k) = f3d(k,:,:,5)
        fwindv_sec(:,:,k) = f3d(k,:,:,6)
        famieu_sec(:,:,k) = f3d(k,:,:,7)
        famiev_sec(:,:,k) = f3d(k,:,:,8)
        tec_sec(:,:,k) = f3d(k,:,:,9)
        ped_sec(:,:,k) = f3d(k,:,:,10)
        hall_sec(:,:,k) = f3d(k,:,:,11)
        qjoule_sec(:,:,k) = f3d(k,:,:,12)
      enddo

      end subroutine dyndiag_bndlatlon
!-----------------------------------------------------------------------
      end module dyndiag_module

!
      subroutine lamdas(tn,barm,o2,o1,ti,te,o2p,op,nop, lxx,lyy,lxy,lyx,
     |  ped_out,hall_out,lev0,lev1,lon0,lon1,lat)
!
! Compute ion drag coefficients lxx,lyy,lxy, and lyx.
!
      use params_module,only: nlonp4
      use magfield_module,only: bmod2,sn2dec,sncsdc
      use cons_module,only: dipmin,rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |  p0,expz,rmass_o2,rmass_o1,boltz,avo
      use magfield_module,only: dipmag
      use input_module,only: colfac
      implicit none
!
! Input args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,   ! neutral temperature (deg K)
     |  barm, ! mean molecular mass
     |  o2,   ! molecular oxygen (mmr)
     |  o1,   ! atomic oxygen (mmr)
     |  ti,   ! ion temperature (deg K)
     |  te,   ! electron temperature (deg K)
     |  o2p,  ! O2+
     |  op,   ! O+
     |  nop   ! NO+
!
! Output args:
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  lxx,     ! lamda XX term
     |  lyy,     ! lamda YY term
     |  lxy,     ! lamda XY term
     |  lyx,     ! lamda YX term
     |  ped_out, ! pedersen conductivity (will be input to dynamo)
     |  hall_out ! hall conductivity (will be input to dynamo)
!
! Local:
      integer :: k,i,l,nlevs,lonbeg,lonend
      integer ::i0,i1,nk,nkm1 ! for addfsech
      real :: sqrt_te,xn2tmp
      real,parameter :: qe=1.602e-19
      real,parameter :: rmass_nop = 30.
      real,dimension(lev0:lev1,lon0:lon1) :: ped_plt,hall_plt
!
! (lon):
      real,dimension(lon0:lon1) :: 
     |  bgauss,omega_i,omega_e,sindip,sindip2
!
! (lev,lon):
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  xnmbarm,    ! for conversion to volume density
     |  vo2,vo,vn2, ! inputs to sub nufac (these were v1,v2,v3)
     |  vt,         ! input to sub nufac ti+tn (was v6)
     |  xn2,        ! N2 (1-o2-o)
     |  rnu_omega,  ! nu/omega_e (was v4)
     |  pedxions,hallxions, ! ped*ions, hall*ions
     |  wk1,wk2,wk3 ! work arrays
!
! (lev,lon,3):
      real,dimension(lev0:lev1,lon0:lon1,3) ::
     |  rmo2,rmo,rmn2, ! numerical factors output by sub nufac
     |  ped,hall,      ! pedersen and hall factors
     |  xions          ! holding array for ions o2p,op,nop
!
! Convenience ints for addfsech calls:
      i0 = lon0
      i1 = lon1
      nk = lev1-lev0+1
      nkm1 = nk-1
      nlevs = nk

!     call addfsech('O2P_LAM',' ',' ',o2p(:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('OP_LAM' ,' ',' ',op (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('NOP_LAM',' ',' ',nop(:,i0:i1),i0,i1,nk,nkm1,lat)
!
! Set local needs:
      do i=lon0,lon1
        bgauss(i) = bmod2(i,lat)
        omega_i(i) = 9.6489E3*bgauss(i)
        omega_e(i) = 1.7588028E7*bgauss(i)
!
! dipmag is in magfield module, dipmin is in cons module:
        if (abs(dipmag(i,lat)) >= dipmin) then
          sindip(i) = sin(dipmag(i,lat))
        else
          sindip(i) = sin(dipmin)
        endif
        sindip2(i) = sindip(i)**2
      enddo ! i=lon0,lon1
!
! vo2,vo,vn2:
      do i=lon0,lon1
        do k=lev0,lev1-1
          vo2(k,i) = o2(k,i)*rmassinv_o2
          vo (k,i) = o1(k,i)*rmassinv_o1
          xn2(k,i) = (1.-o2(k,i)-o1(k,i))
          vn2(k,i) = xn2(k,i)*rmassinv_n2
          if (vn2(k,i) < 1.e-20) vn2(k,i) = 1.e-20
          vt (k,i) = 0.5*(ti(k,i)+tn(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('VO2',' ',' ',vo2,i0,i1,nk,nkm1,lat)
!     call addfsech('VO' ,' ',' ',vo ,i0,i1,nk,nkm1,lat)
!     call addfsech('VN2',' ',' ',vn2,i0,i1,nk,nkm1,lat)
!     call addfsech('VT' ,' ',' ',vt ,i0,i1,nk,nkm1,lat)
!
! Generate numerical factors in rmo2,rmo,rmn2 (lev0:lev1,lon0:lon1,3):
! (this was sub new.F in earlier versions)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! O2:
          rmo2(k,i,1)=2.59E-11*sqrt(vt(k,i))*       ! O2+ - O2
     |                (1.-0.073*alog10(vt(k,i)))**2
          rmo2(k,i,2)=6.64E-10                      ! O+  - O2
          rmo2(k,i,3)=4.27E-10                      ! NO+ - O2
!
! O:
          rmo(k,i,1)=2.31E-10                       ! O2+ - O
          rmo(k,i,2)=3.67e-11*sqrt(vt(k,i))*        ! O+  - O
     |      (1.-0.064*alog10(vt(k,i)))**2*colfac    
          rmo(k,i,3)=2.44E-10                       ! NO+ - O
!
! N2:
          rmn2(k,i,1)=4.13E-10                      ! O2+ - N2
          rmn2(k,i,2)=6.82E-10                      ! O+  - N2
          rmn2(k,i,3)=4.34E-10                      ! NO+ - N2
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('RMO2_1',' ',' ',rmo2(lev0:lev1,lon0:lon1,1),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMO2_2',' ',' ',rmo2(lev0:lev1,lon0:lon1,2),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMO2_3',' ',' ',rmo2(lev0:lev1,lon0:lon1,3),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMO_1' ,' ',' ',rmo(lev0:lev1,lon0:lon1,1) ,
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMO_2' ,' ',' ',rmo(lev0:lev1,lon0:lon1,2) ,
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMO_3' ,' ',' ',rmo(lev0:lev1,lon0:lon1,3) ,
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMN2_1',' ',' ',rmn2(lev0:lev1,lon0:lon1,1),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMN2_2',' ',' ',rmn2(lev0:lev1,lon0:lon1,2),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('RMN2_3',' ',' ',rmn2(lev0:lev1,lon0:lon1,3),
!    |  i0,i1,nk,nkm1,lat)
!
! Multiply by major species number densities:
      do l=1,3
        do i=lon0,lon1
          do k=lev0,lev1-1
            rmo2(k,i,l) = rmo2(k,i,l)*vo2(k,i)
            rmo (k,i,l) = rmo (k,i,l)*vo (k,i)
            rmn2(k,i,l) = rmn2(k,i,l)*vn2(k,i)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! l=1,3
!
! vo2,vo,vn2 = sum(rmo2(l)+rmo(l)+rmn2(l))
      vo2(:,:) = 0.
      vo (:,:) = 0.
      vn2(:,:) = 0.
      do i=lon0,lon1
        do k=lev0,lev1-1
          vo2(k,i) = vo2(k,i)+rmo2(k,i,1)+rmo(k,i,1)+rmn2(k,i,1)
          vo (k,i) = vo (k,i)+rmo2(k,i,2)+rmo(k,i,2)+rmn2(k,i,2)
          vn2(k,i) = vn2(k,i)+rmo2(k,i,3)+rmo(k,i,3)+rmn2(k,i,3)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('VO2',' ',' ',vo2,i0,i1,nk,nkm1,lat)
!     call addfsech('VO' ,' ',' ',vo ,i0,i1,nk,nkm1,lat)
!     call addfsech('VN2',' ',' ',vn2,i0,i1,nk,nkm1,lat)
!
! Convert to number densities:
      do i=lon0,lon1
        do k=lev0,lev1-1
          xnmbarm(k,i) = p0*expz(k)*.5*(barm(k,i)+barm(k+1,i))/
     |      (boltz*tn(k,i))
          vo2(k,i) = vo2(k,i)*xnmbarm(k,i)*rmass_o2/omega_i(i)
          vo (k,i) = vo (k,i)*xnmbarm(k,i)*rmass_o1/omega_i(i)
!
! N2 is a minor contributor to condumctivity in lower thermosphere,
! and NO+ is the major ion. Using rmass_nop instead of rmass_n2.
          vn2(k,i) = vn2(k,i)*xnmbarm(k,i)*rmass_nop/omega_i(i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('XNMBARM',' ',' ',xnmbarm,i0,i1,nk,nkm1,lat)
!     call addfsech('VO2',' ',' ',vo2,i0,i1,nk,nkm1,lat)
!     call addfsech('VO' ,' ',' ',vo ,i0,i1,nk,nkm1,lat)
!     call addfsech('VN2',' ',' ',vn2,i0,i1,nk,nkm1,lat)
!
! rnu_omega = nu/omega_e
      do i=lon0,lon1
        do k=lev0,lev1-1
          sqrt_te = sqrt(te(k,i))
          xn2tmp = xn2(k,i)*rmassinv_n2
          if (xn2tmp < 1.e-20) xn2tmp = 1.e-20
!
! Original Banks & Kockards:
          rnu_omega(k,i) = xnmbarm(k,i)*(2.33e-11*xn2tmp*te(k,i)*
     |      (1.-1.21e-4*te(k,i))+1.82e-10*o2(k,i)*rmassinv_o2*
     |      sqrt_te*(1.+3.6e-2*sqrt_te)+8.9e-11*o1(k,i)*rmassinv_o1*
     |      sqrt_te*(1.+5.7e-4*te(k,i)))
          rnu_omega(k,i) = rnu_omega(k,i)/omega_e(i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('NU_OMEGA',' ',' ',rnu_omega,i0,i1,nk,nkm1,lat)
!
! Pedersen and Hall factors:
      do i=lon0,lon1
        do k=lev0,lev1-1
          ped(k,i,1) = rnu_omega(k,i)/(1.+rnu_omega(k,i)**2)+
     |      vo2(k,i)/(1.+vo2(k,i)**2)
          ped(k,i,2) = rnu_omega(k,i)/(1.+rnu_omega(k,i)**2)+
     |      vo (k,i)/(1.+vo (k,i)**2)
          ped(k,i,3) = rnu_omega(k,i)/(1.+rnu_omega(k,i)**2)+
     |      vn2(k,i)/(1.+vn2(k,i)**2)
!
          hall(k,i,1) = 1./(1.+rnu_omega(k,i)**2)-1./(1.+vo2(k,i)**2)
          hall(k,i,2) = 1./(1.+rnu_omega(k,i)**2)-1./(1.+vo (k,i)**2)
          hall(k,i,3) = 1./(1.+rnu_omega(k,i)**2)-1./(1.+vn2(k,i)**2)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('PED_1',' ',' ',ped(lev0:lev1,lon0:lon1,1),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('PED_2',' ',' ',ped(lev0:lev1,lon0:lon1,2),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('PED_3',' ',' ',ped(lev0:lev1,lon0:lon1,3),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('HALL_1',' ',' ',hall(lev0:lev1,lon0:lon1,1),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('HALL_2',' ',' ',hall(lev0:lev1,lon0:lon1,2),
!    |  i0,i1,nk,nkm1,lat)
!     call addfsech('HALL_3',' ',' ',hall(lev0:lev1,lon0:lon1,3),
!    |  i0,i1,nk,nkm1,lat)
!
! pedxions = ped*ions, hallxions = hall*ions
      pedxions(:,:) = 0.  
      hallxions(:,:) = 0.  
      xions(:,lon0:lon1,1) = o2p(:,lon0:lon1)
      xions(:,lon0:lon1,2) = op (:,lon0:lon1)
      xions(:,lon0:lon1,3) = nop(:,lon0:lon1)
      do l=1,3
        do i=lon0,lon1
          do k=lev0,lev1-1
            pedxions (k,i) = pedxions (k,i)+ped (k,i,l)*xions(k,i,l)
            hallxions(k,i) = hallxions(k,i)+hall(k,i,l)*xions(k,i,l)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! l=1,3
!
! small (3e-8) "diamond diff" ~tgcm15 in pedxions:
!     call addfsech('XNMBARM',' ',' ' ,xnmbarm  ,i0,i1,nk,nkm1,lat)
!     call addfsech('PEDXIONS',' ',' ',pedxions ,i0,i1,nk,nkm1,lat)
!     call addfsech('HALXIONS',' ',' ',hallxions,i0,i1,nk,nkm1,lat)
!
! Set lxx, lyy, lxy:
      do i=lon0,lon1
        do k=lev0,lev1-1
          wk1(k,i) = pedxions(k,i)*bgauss(i)*1.e-1*qe*avo/xnmbarm(k,i)
          wk2(k,i) = wk1(k,i)*sindip2(i)
          wk3(k,i) = hallxions(k,i)*bgauss(i)*1.e-1*qe*avo/xnmbarm(k,i)
        enddo ! k=lev0,lev1-1
        do k=lev0,lev1-2
          lxx(k+1,i) = sqrt(wk1(k,i)*wk1(k+1,i))
          lyy(k+1,i) = sqrt(wk2(k,i)*wk2(k+1,i))
          lxy(k+1,i) = sqrt(wk3(k,i)*wk3(k+1,i))
        enddo ! k=lev0,lev1-2
!
! Lower boundary:
        lxx(lev0,i) = sqrt(wk1(lev0,i)**3/wk1(lev0+1,i))
        lyy(lev0,i) = sqrt(wk2(lev0,i)**3/wk2(lev0+1,i))
        lxy(lev0,i) = sqrt(wk3(lev0,i)**3/wk3(lev0+1,i))
!
! Upper boundary:
        lxx(lev1,i) = sqrt(wk1(lev1-1,i)**3/wk1(lev1-2,i))
        lyy(lev1,i) = sqrt(wk2(lev1-1,i)**3/wk2(lev1-2,i))
        lxy(lev1,i) = sqrt(wk3(lev1-1,i)**3/wk3(lev1-2,i))
      enddo ! i=lon0,lon1
!
! lxy = lxy*sindip:
      do i=lon0,lon1
        do k=lev0,lev1
          lxy(k,i) = lxy(k,i)*sindip(i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('LXX',' ',' ',lxx(:,i0:i1),i0,i1,nk,nk,lat)
!     call addfsech('LYY',' ',' ',lyy(:,i0:i1),i0,i1,nk,nk,lat)
!     call addfsech('LXY',' ',' ',lxy(:,i0:i1),i0,i1,nk,nk,lat)
!     call addfsech('LYX',' ',' ',lxy(:,i0:i1),i0,i1,nk,nk,lat)

!
! Define 3-d fields for dynamo (in old dynamo.h) (do this later)

!
! Rotate lamdas (formerly sub rotate.F):
      do i=lon0,lon1
        do k=lev0,lev1
          wk1(k,i) = lxx(k,i)*(1.-sn2dec(i,lat))+lyy(k,i)*sn2dec(i,lat)
          wk2(k,i) = lyy(k,i)*(1.-sn2dec(i,lat))+lxx(k,i)*sn2dec(i,lat)
          lyx(k,i) = lxy(k,i)-(lyy(k,i)-lxx(k,i))*sncsdc(i,lat)
          lxy(k,i) = lxy(k,i)+(lyy(k,i)-lxx(k,i))*sncsdc(i,lat)
          lxx(k,i) = wk1(k,i)
          lyy(k,i) = wk2(k,i)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1
!
! Periodic points:
! Probably not necessary:
!     call periodic_f2d(lxx,lon0,lon1,nlevs)
!     call periodic_f2d(lyy,lon0,lon1,nlevs)
!     call periodic_f2d(lxy,lon0,lon1,nlevs)
!     call periodic_f2d(lyx,lon0,lon1,nlevs)

! (For tgcm15, these are from rotate.F, not lamdas.F)
!     call addfsech('LXX',' ',' ',lxx(:,i0:i1),i0,i1,nk,nk,lat)
!     call addfsech('LYY',' ',' ',lyy(:,i0:i1),i0,i1,nk,nk,lat)
!     call addfsech('LXY',' ',' ',lxy(:,i0:i1),i0,i1,nk,nk,lat)
!     call addfsech('LYX',' ',' ',lyx(:,i0:i1),i0,i1,nk,nk,lat)
!
! Return pedersen and hall conductivities for use in the dynamo:
! This is not strictly necessary unless user input flag DYNAMO is set,
!   but do it anyway, to save the 3d ped and hall fields.
! Note tgcm15 (lamdas.F) assigned do i=1,imax+1: f(i) <= f(i+2),
!   so to match that ped_out and hall_out would need to be rotated
!   2 lon points to the "left", i.e., do i=1,nlon+1: f(i) = f(i+2)
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2

      do i=lonbeg,lonend
        do k=lev0,lev1
          ped_out (k,i) = pedxions (k,i)*1.e10*qe/bgauss(i)
          hall_out(k,i) = hallxions(k,i)*1.e10*qe/bgauss(i)
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     call addfsech('SPED_LAM',' ',' ',ped_out(:,i0:i1),i0,i1,nk,nk,
!    |  lat)
!     call addfsech('SHAL_LAM',' ',' ',hall_out(:,i0:i1),i0,i1,nk,nk,
!    |  lat)
     
!     real,dimension(lev0:lev1,lon0:lon1) :: ped_plt,hall_plt
!     do i=lon0,lon1
!       ped_plt(:,i)  = ped_out(:,i)
!       hall_plt(:,i) = hall_out(:,i)
!     enddo
!     call addfsech('SIGMAPED',' ',' ',ped_plt,i0,i1,nk,nk,lat)
!     call addfsech('SIGMAHAL',' ',' ',hall_plt,i0,i1,nk,nk,lat)

      end subroutine lamdas

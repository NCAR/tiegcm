!
      module oplus_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! tiegcmlicense.txt.
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
      implicit none
!
! 2024/09 Haonan Wu: rewrite the whole oplus module
! to split calculations independent of O+ sub-cycling
! this will speed up the model a little
! when the number of O+ sub-cycling is large
!
      contains
!-----------------------------------------------------------------------
      subroutine oplus(tn,te,ti,o2,o1,he,n2,n2d,ne,un,vn,w,ui,vi,wi,
     |  xnmbar,scht,op,optm1,opout,optm1out,xiop2p,xiop2d,Fe,Fn,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Update O+ ion at 3d task subdomain.
! Outputs are opout, optm1out, xiop2p, xiop2d, Fe, and Fn,
! all other args are input.
!
      use params_module,only: rp
      use input_module,only: nstep_sub
      use mpi_module,only: mp_polelats_f3d,mp_bndlats_f3d,mp_bndlons_f3d
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn, te, ti, ! neutral, electron, and ion temperatures (deg K)
     |  o2, o1, he, ! O2, O, He mass mixing ratios
     |  n2,         ! N2 mass mixing ratio
     |  n2d,        ! N(2D) mass mixing ratio
     |  ne,         ! electron density
     |  un,vn,w,    ! neutral wind velocities (zonal, meridional, omega)
     |  ui,vi,wi,   ! zonal, meridional, and vertical ion velocities
     |  xnmbar,     ! p0*e(-z)*mbar/kT
     |  scht,       ! scale height
     |  op,         ! O+ ion
     |  optm1       ! O+ at time n-1
!
! Output fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  opout,    ! O+ output for next timestep
     |  optm1out, ! O+ output for time n-1
     |  xiop2p,xiop2d,Fe,Fn
!
! Local:
      integer :: k,i,lat,nlevs,istep
      real(rp),dimension(lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  dvb,ubcrhs,lbcrhs
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  tp,dj,bdotu,op_prod,op_loss,diffp,diffq,diffr,
     |  driftp,driftq,driftr,windp,windq,windr,
     |  p_coeff,q_coeff,r_coeff,optm1_smooth,
     |  op_sub,optm1_sub,opout_sub,optm1out_sub,
     |  diffj,diffexp,vdotn_h,bdotdh_bvel
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,2) ::
     |  tmpf3d
!
! Number of pressure levels (this will equal nlevp1):
      nlevs = lev1-lev0+1 ! for bndlons calls
!
      call prep_oplus(tn,te,ti,
     |  o2,o1,he,n2,n2d,ne,un,vn,w,wi,xnmbar,scht,
     |  dvb,ubcrhs,lbcrhs,tp,dj,bdotu,xiop2p,xiop2d,
     |  op_prod,op_loss,diffp,diffq,diffr,
     |  driftp,driftq,driftr,windp,windq,windr,
     |  p_coeff,q_coeff,r_coeff,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1
            op_sub(k,i,lat) = op(k,i,lat)
            optm1_sub(k,i,lat) = optm1(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      do istep=1,nstep_sub
        call smooth_oplus(optm1_sub,optm1_smooth,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
!
        call iterate_oplus(dvb,ubcrhs,lbcrhs,tp,dj,bdotu,
     |    optm1_smooth,op_prod,ui,vi,scht,p_coeff,q_coeff,r_coeff,
     |    op_sub,opout_sub,diffj,diffexp,vdotn_h,bdotdh_bvel,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
!
        if (istep == nstep_sub)
     |    call calc_terms(xnmbar,diffj,Fe,Fn,
     |    opout_sub,optm1_smooth,op_prod,op_loss,
     |    diffexp,diffp,diffq,diffr,
     |    vdotn_h,driftp,driftq,driftr,
     |    bdotdh_bvel,windp,windq,windr,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
!
! Filter updated O+:
        call filter_op(opout_sub(:,lon0:lon1,lat0:lat1),
     |    lev0,lev1,lon0,lon1,lat0,lat1,'OPLUS')
!
!       do lat=lat0,lat1
!         call addfld('OP_FILT',' ',' ',opout_sub(:,lon0:lon1,lat),
!    |      'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       enddo ! lat=lat0,lat1
!
        call post_oplus(op_sub,optm1_sub,opout_sub,optm1out_sub,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
!
        do lat=lat0,lat1
          do i=lon0,lon1
            do k=lev0,lev1
              tmpf3d(k,i,lat,1) = opout_sub(k,i,lat)
              tmpf3d(k,i,lat,2) = optm1out_sub(k,i,lat)
            enddo ! k=lev0,lev1
          enddo ! i=lon0,lon1
        enddo ! lat=lat0,lat1
        call mp_polelats_f3d(tmpf3d(:,lon0:lon1,:,:),
     |    lev0,lev1,lon0,lon1,lat0,lat1,2,(/1.,1./))
        call mp_bndlats_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2)
        call mp_bndlons_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2,0)
        do lat=lat0-2,lat1+2
          do i=lon0-2,lon1+2
            do k=lev0,lev1
              opout_sub(k,i,lat) = tmpf3d(k,i,lat,1)
              optm1out_sub(k,i,lat) = tmpf3d(k,i,lat,2)
            enddo ! k=lev0,lev1
          enddo ! i=lon0-2,lon1+2
        enddo ! lat=lat0-2,lat1+2
!
        do lat=lat0-2,lat1+2
          do i=lon0-2,lon1+2
            do k=lev0,lev1
              op_sub(k,i,lat) = opout_sub(k,i,lat)
              optm1_sub(k,i,lat) = optm1out_sub(k,i,lat)
            enddo ! k=lev0,lev1
          enddo ! i=lon0-2,lon1+2
        enddo ! lat=lat0-2,lat1+2
      enddo ! istep=1,nstep_sub
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1
            opout(k,i,lat) = opout_sub(k,i,lat)
            optm1out(k,i,lat) = optm1out_sub(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      end subroutine oplus
!-----------------------------------------------------------------------
      subroutine prep_oplus(tn,te,ti,
     |  o2,o1,he,n2,n2d,ne,un,vn,w,wi,xnmbar,scht,
     |  dvb,ubcrhs,lbcrhs,tp,dj,bdotu,xiop2p,xiop2d,
     |  op_prod,op_loss,diffp,diffq,diffr,
     |  driftp,driftq,driftr,windp,windq,windr,
     |  p_coeff,q_coeff,r_coeff,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate terms that do not change with O+ sub-cycling
!
      use params_module,only: nlonp4,dz,zpmid,spval,rp
      use cons_module,only: rmass_op,gask,grav,re,cs,dphi,dlamda,
     |  boltz,dtx2inv,pi,dtr,rmassinv_o2,rmassinv_o1,
     |  rmassinv_he,rmassinv_n2,rmassinv_n2d
      use input_module,only: nstep_sub,colfac,
     |  opdiffcap,opdiffrate,opdifflev
      use chapman_module,only: chi     ! was t2 in old sub opflux
      use qrj_module,only:
     |  qop2p, ! O+(2p) ionization from qrj, used in xiop2p
     |  qop2d, ! O+(2d) ionization from qrj, used in xiop2d
     |  qop    ! O+ ionization from qrj
      use chemrates_module,only: ! needed chemical reaction rates
     |  rk1 ,rk2 ,rk10,rk16,rk17,rk18,rk19,rk20,
     |  rk21,rk22,rk23,rk24,rk25,rk26,rk27
      use magfield_module,only: bx,by,bz,bmod2, ! (nlonp4,-1:nlat+2)
     |  rlatm
      use mpi_module,only: mp_periodic_f2d,
     |  mp_polelats_f3d,mp_bndlats_f3d,mp_bndlons_f3d
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_BXYZ,mkdiag_BMAG
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tn, te, ti, ! neutral, electron, and ion temperatures (deg K)
     |  o2, o1, he, ! O2, O, He mass mixing ratios
     |  n2,         ! N2 mass mixing ratio
     |  n2d,        ! N(2D) mass mixing ratio
     |  ne,         ! electron density
     |  un,vn,w,    ! neutral wind velocities (zonal, meridional, omega)
     |  wi,         ! vertical ion velocities
     |  xnmbar,     ! p0*e(-z)*mbar/kT
     |  scht        ! scale height
!
! Output fields (full 2d task subdomain):
      real(rp),dimension(lon0-2:lon1+2,lat0-2:lat1+2),intent(out) ::
     |  dvb,        ! output of sub divb
     |  ubcrhs,lbcrhs
!
! Output fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  tp,          ! Plasma temperature 0.5*(te+ti)
     |  dj,          ! diffusion coefficients (s13,s14,s15)
     |  bdotu,       ! was s7 (B.U)
     |  op_prod,     ! production rate
     |  op_loss,     ! was s13
     |  xiop2p,xiop2d,
     |  diffp,diffq,diffr,      ! tridiagonal coefficients related to diffusion
     |  driftp,driftq,driftr,   ! tridiagonal coefficients related to ion drift
     |  windp,windq,windr,      ! tridiagonal coefficients related to wind
     |  p_coeff,q_coeff,r_coeff ! coefficients for tridiagonal solver (s1,s2,s3)
!
! Local:
      real(rp),parameter ::
     |  phid =  2.0e8,
     |  phin = -2.0e8,
     |  ppolar = 0.,
     |  mp = 1.6726e-24 ! proton mass (g)
      integer :: k,i,lat,lonbeg,lonend,nlevs
      real(rp) :: a,djmin,gmr
      real(rp),dimension(lon0:lon1,lat0:lat1) ::
     |  phi,   opflux, ! upward number flux of O+ (returned by sub oplus_flux) (t7)
     |  dbxdx, dbycdy, ! lon/lat derivatives of unit vector B
     |  dvb_h,         ! horizontal divergence of unit vector B
     |  ubca,  ubcb    ! O+ upper boundary condition (were t2,t3)
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  bdotdh_djbz,   ! (B(H).DEL(H))(D*BZ)
     |  hdz,hdzi,      ! was s15
     |  tp1,           ! 0.5*(te+ti)
     |  djint,         ! was s11
     |  tphdz1,tphdz0, ! were s13,s12 (using gmr)
     |  bdotdh_djbz2,  ! (B(H).DEL(H))(D*BZ) / (D*BZ)
     |  divbz,         ! was s7 (DIV(B)+(DH*D*BZ)/(D*BZ)
     |  hdzmbz,hdzpbz, ! were s10,s9
     |  djint_tphdz0,djint_tphdz1,
     |  djint_tphdz0_1,djint_tphdz1_1,
     |  bdotu_m1,bdotu_p1,wii
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  tr,            ! Reduced temperature 0.5*(tn+ti)
     |  vni,           ! O+ collision frequency
     |  djbz,          ! D*BZ
     |  o2_cm3,o1_cm3,he_cm3,n2_cm3,n2d_cm3
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,2) ::
     |  tmpf3d
!
! Number of pressure levels (this will equal nlevp1):
      nlevs = lev1-lev0+1 ! for bndlons calls
!
! Calculate O+ number flux in opflux for sub oplus (was sub opflux).
!
! Change upper boundary flux as electron heat flux in settei
      do lat=lat0,lat1
        do i=lon0,lon1
          if (abs(rlatm(i,lat)) >= pi/4.5) then
            a = 1.
          elseif (abs(rlatm(i,lat)) <= pi/18.) then
            a = 0.
          else
            a = .5*(1.-cos(abs(rlatm(i,lat))*6.-pi/3.))
          endif
!
          if (chi(i,lat) >= 100.*dtr) then
            phi(i,lat) = phin
          elseif (chi(i,lat) <= 80.*dtr) then
            phi(i,lat) = phid
          else
            phi(i,lat) = .5*(phid+phin+(phid-phin)*cos(chi(i,lat)*9.))
          endif
!
          opflux(i,lat) = phi(i,lat)*a
!
! Add ppolar if magnetic latitude >= 60 degrees:
          if (abs(rlatm(i,lat)) >= pi/3.)
     |      opflux(i,lat) = opflux(i,lat)+ppolar
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
!     call addfld('OPFLUX',' ',' ',opflux,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Evaluate divergence of B, the unit magnetic field vector.
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
      dvb = 0.
      do lat=lat0,lat1
        do i=lonbeg,lonend
          dbxdx(i,lat) = (bx(i+1,lat)-bx(i-1,lat))/(2.*dlamda)
          dbycdy(i,lat) =
     |      (cs(lat+1)*by(i,lat+1)-
     |       cs(lat-1)*by(i,lat-1))/(2.*dphi)
          dvb_h(i,lat) = (dbxdx(i,lat)+dbycdy(i,lat))/(re*cs(lat))
          dvb(i,lat) = dvb_h(i,lat)+2.*bz(i,lat)/re
        enddo ! i=lonbeg,lonend
      enddo ! lat=lat0,lat1
!
      call mp_periodic_f2d(dvb(lon0:lon1,lat0:lat1),
     |  lon0,lon1,lat0,lat1,1)
!
!     call addfld('OP_DIVB',' ',' ',dvb(lon0:lon1,lat0:lat1),
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
      call mkdiag_BXYZ('BX',bx(lon0:lon1,lat0:lat1),lon0,lon1,lat0,lat1)
      call mkdiag_BXYZ('BY',by(lon0:lon1,lat0:lat1),lon0,lon1,lat0,lat1)
      call mkdiag_BXYZ('BZ',bz(lon0:lon1,lat0:lat1),lon0,lon1,lat0,lat1)
      call mkdiag_BMAG('BMAG',bmod2(lon0:lon1,lat0:lat1),
     |  lon0,lon1,lat0,lat1)
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
!
! Plasma temperature:
            tp(k,i,lat) = 0.5*(te(k,i,lat)+ti(k,i,lat))
!
! Set reduced temperature (average of tn and ti)
!
! 1/2/16 btf:
!   This was previously set incorrectly as tp = 0.5*(te+ti))
!   (see also plasma temp tp). Correcting this resulted in a
!   40% change in ambipolar diffusion DJ, but only 10% or so
!   change in O+ and Ne, and very small differences in NMF2,HMF2,TEC.
!   Wenbin feels there should be larger changes in these fields
!   as a consequence of fixing this bug.
            tr(k,i,lat) = 0.5*(tn(k,i,lat)+ti(k,i,lat))
          enddo ! k=lev0,lev1-1
!
          tp(lev1,i,lat) = 2.*tp(lev1-1,i,lat)-tp(lev1-2,i,lat)
          tr(lev1,i,lat) = 2.*tr(lev1-1,i,lat)-tr(lev1-2,i,lat)
!
          do k=lev0,lev1
            o2_cm3(k,i,lat) = xnmbar(k,i,lat)*o2(k,i,lat)*rmassinv_o2
            o1_cm3(k,i,lat) = xnmbar(k,i,lat)*o1(k,i,lat)*rmassinv_o1
            he_cm3(k,i,lat) = xnmbar(k,i,lat)*he(k,i,lat)*rmassinv_he
            n2_cm3(k,i,lat) = xnmbar(k,i,lat)*n2(k,i,lat)*rmassinv_n2
            n2d_cm3(k,i,lat) = xnmbar(k,i,lat)*n2d(k,i,lat)*rmassinv_n2d
!
! 8/28/13 btf: Use n2=1-o2-o-he, and include O+/He collision rate from Wenbin:
            vni(k,i,lat) = 1e-10*
     |        (6.64 *o2_cm3(k,i,lat)+
     |         1.32 *he_cm3(k,i,lat)+
     |         6.82 *n2_cm3(k,i,lat)+
     |         0.367*o1_cm3(k,i,lat)*sqrt(tr(k,i,lat))*
     |         (1.-0.064*log10(tr(k,i,lat)))**2*colfac)
!
! ambipolar diffusion coefficient
            dj(k,i,lat) = boltz/(mp*rmass_op*vni(k,i,lat))
          enddo ! k=lev0,lev1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! 1/30/16 btf:
! Cap ambipolar diffusion coefficient. Namelist parameter OPDIFFCAP.
! This was tested with various values (1.5e8, 3e8, 6e8, 8e8),
! and was found to improve numerical stability in some storm cases,
! for example the November, 2003 and July, 2000 storms, with Weimer
! potential model and IMF/OMNI data. Both of these cases generally
! will not complete if timestep is longer than 10 sec. The Nov, 2003
! may succeed with opdiffcap turned off, but the July, 2000
! "Bastille day storm" will succeed only with step=10 and opdiffcap=6.e8.
!
! 2024/08 Haonan Wu: Change the original constant O+ diffusion cap
!   to an altitude-dependent cap (use logistic function for now)
      if (opdiffcap/=0. .and.
     |    opdiffrate/=0. .and.
     |    opdifflev/=spval) then
        do lat=lat0-2,lat1+2
          do i=lon0-2,lon1+2
            do k=lev0,lev1
              djmin = opdiffcap/(1+exp(opdiffrate*(zpmid(k)-opdifflev)))
              if (dj(k,i,lat) > djmin) dj(k,i,lat) = djmin
            enddo ! k=lev0,lev1
          enddo ! i=lon0-2,lon1+2
        enddo ! lat=lat0-2,lat1+2
      endif
!
      do lat=lat0,lat1
        do i=lon0,lon1
!
! bdotu = B.U (s7)
          do k=lev0,lev1-1
            bdotu(k,i,lat) =
     |        bx(i,lat)*un(k,i,lat)+by(i,lat)*vn(k,i,lat)+
     |        scht(k,i,lat)*bz(i,lat)*0.5*(w(k,i,lat)+w(k+1,i,lat))
          enddo ! k=lev0,lev1-1
          bdotu(lev1,i,lat) = 2.*bdotu(lev1-1,i,lat)-bdotu(lev1-2,i,lat)
!
          do k=lev0,lev1
            djbz(k,i,lat) = dj(k,i,lat)*bz(i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Fill up halo points of bdotu, djbz (used in bdotdh calls)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            tmpf3d(k,i,lat,1) = bdotu(k,i,lat)
            tmpf3d(k,i,lat,2) = djbz(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
      call mp_polelats_f3d(tmpf3d(:,lon0:lon1,:,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,2,(/1.,1./))
      call mp_bndlats_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2)
      call mp_bndlons_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2,0)
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1
            bdotu(k,i,lat) = tmpf3d(k,i,lat,1)
            djbz(k,i,lat) = tmpf3d(k,i,lat,2)
          enddo ! k=lev0,lev1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      call bdotdh(djbz,bdotdh_djbz,lev0,lev1,lon0,lon1,lat0,lat1)
!
! gmr = G*M(O+)/(2.*R)
      gmr = grav*rmass_op/(2.*gask)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            hdz(k,i,lat) = 1./(scht(k,i,lat)*dz)
          enddo ! k=lev0,lev1
!
          do k=lev0,lev1-1
            tp1(k+1,i,lat) = tp(k,i,lat)
            bdotu_m1(k+1,i,lat) = bdotu(k,i,lat)
            bdotu_p1(k,i,lat) = bdotu(k+1,i,lat)
            hdzi(k+1,i,lat) = 0.5*(hdz(k,i,lat)+hdz(k+1,i,lat))
            djint(k+1,i,lat) = 0.5*(dj(k,i,lat)+dj(k+1,i,lat))
            wii(k,i,lat) = 0.5*(wi(k,i,lat)+wi(k+1,i,lat))
          enddo ! k=lev0,lev1-1
!
          tp1(lev0,i,lat) = 2.*tp(lev0,i,lat)-tp(lev0+1,i,lat)
          bdotu_m1(lev0,i,lat) =
     |      2.*bdotu(lev0,i,lat)-bdotu(lev0+1,i,lat)
          bdotu_p1(lev1,i,lat) =
     |      2.*bdotu(lev1,i,lat)-bdotu(lev1-1,i,lat)
          hdzi(lev0,i,lat) = 1.5*hdz(lev0,i,lat)-0.5*hdz(lev0+1,i,lat)
          djint(lev0,i,lat) = 1.5*dj(lev0,i,lat)-0.5*dj(lev0+1,i,lat)
          wii(lev1,i,lat) = 1.5*wi(lev1,i,lat)-0.5*wi(lev1-1,i,lat)
!
          do k=lev0,lev1
            tphdz1(k,i,lat) = 2.*tp (k,i,lat)*hdzi(k,i,lat)+gmr ! s13
            tphdz0(k,i,lat) = 2.*tp1(k,i,lat)*hdzi(k,i,lat)-gmr ! s12
!
! divbz = (DIV(B)+(DH*D*BZ)/(D*BZ) (was s7)
            bdotdh_djbz2(k,i,lat) = bdotdh_djbz(k,i,lat)/djbz(k,i,lat)
            divbz(k,i,lat) =
     |        (dvb(i,lat)+bdotdh_djbz2(k,i,lat))/bz(i,lat)
!
! hdzmbz = (1./(H*DZ)-(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s10)
! hdzpbz = (1./(H*DZ)+(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s9 )
            hdzmbz(k,i,lat) =
     |        (hdz(k,i,lat)-0.5*divbz(k,i,lat))*bz(i,lat)**2
            hdzpbz(k,i,lat) =
     |        (hdz(k,i,lat)+0.5*divbz(k,i,lat))*bz(i,lat)**2
!
            djint_tphdz0(k,i,lat) = djint(k,i,lat)*tphdz0(k,i,lat)
            djint_tphdz1(k,i,lat) = djint(k,i,lat)*tphdz1(k,i,lat)
          enddo ! k=lev0,lev1
!
          do k=lev0,lev1-1
            djint_tphdz0_1(k,i,lat) = djint_tphdz0(k+1,i,lat)
            djint_tphdz1_1(k,i,lat) = djint_tphdz1(k+1,i,lat)
!
! xiop2p and xiop2d are outputs:
            xiop2p(k,i,lat) =
     |        0.5*(qop2p(k,i,lat)+qop2p(k+1,i,lat))/
     |        ((rk16+rk17)*n2_cm3(k,i,lat)+
     |         rk18*o1_cm3(k,i,lat)+
     |         (rk19(k,i,lat)+rk20(k,i,lat))*ne(k,i,lat)+
     |         rk21+rk22)
            xiop2d(k,i,lat) =
     |        (0.5*(qop2d(k,i,lat)+qop2d(k+1,i,lat))+
     |         rk20(k,i,lat)*xiop2p(k,i,lat)*ne(k,i,lat)+
     |         rk22*xiop2p(k,i,lat))/
     |        (rk23*n2_cm3(k,i,lat)+
     |         rk24*o1_cm3(k,i,lat)+
     |         rk25(k,i,lat)*ne(k,i,lat)+
     |         rk26*o2_cm3(k,i,lat)+
     |         rk27)
!
! Sources and sinks
            op_prod(k,i,lat) =
     |        0.5*(qop(k,i,lat)+qop(k+1,i,lat))+
     |        rk18*xiop2p(k,i,lat)*o1_cm3(k,i,lat)+
     |        rk19(k,i,lat)*xiop2p(k,i,lat)*ne(k,i,lat)+
     |        rk21*xiop2p(k,i,lat)+
     |        rk24*xiop2d(k,i,lat)*o1_cm3(k,i,lat)+
     |        rk25(k,i,lat)*xiop2d(k,i,lat)*ne(k,i,lat)+
     |        rk27*xiop2d(k,i,lat)
            op_loss(k,i,lat) =
     |        rk1(k,i,lat)*o2_cm3(k,i,lat)+
     |        rk2(k,i,lat)*n2_cm3(k,i,lat)+
     |        rk10*n2d_cm3(k,i,lat)
          enddo ! k=lev0,lev1-1
!
          djint_tphdz0_1(lev1,i,lat) =
     |      2.*djint_tphdz0(lev1,i,lat)-djint_tphdz0(lev1-1,i,lat)
          djint_tphdz1_1(lev1,i,lat) =
     |      2.*djint_tphdz1(lev1,i,lat)-djint_tphdz1(lev1-1,i,lat)
          xiop2p(lev1,i,lat) =
     |      2.*xiop2p(lev1-1,i,lat)-xiop2p(lev1-2,i,lat)
          xiop2d(lev1,i,lat) =
     |      2.*xiop2d(lev1-1,i,lat)-xiop2d(lev1-2,i,lat)
          op_prod(lev1,i,lat) =
     |      2.*op_prod(lev1-1,i,lat)-op_prod(lev1-2,i,lat)
          op_loss(lev1,i,lat) =
     |      2.*op_loss(lev1-1,i,lat)-op_loss(lev1-2,i,lat)
!
          do k=lev0,lev1
!
! Begin coefficients p_coeff, q_coeff, r_coeff  (s1,s2,s3)
            diffp(k,i,lat) =
     |          hdzmbz(k,i,lat)*djint_tphdz0  (k,i,lat)
            diffq(k,i,lat) =
     |        -(hdzpbz(k,i,lat)*djint_tphdz0_1(k,i,lat)+
     |          hdzmbz(k,i,lat)*djint_tphdz1  (k,i,lat))
            diffr(k,i,lat) =
     |          hdzpbz(k,i,lat)*djint_tphdz1_1(k,i,lat)
!
! Continue coefficients with vertical ion velocity:
            driftp(k,i,lat) = wii(k,i,lat)*0.5*hdz(k,i,lat)
            driftq(k,i,lat) = -wii(k,i,lat)*6./re
            driftr(k,i,lat) = -wii(k,i,lat)*0.5*hdz(k,i,lat)
!
! Continue coefficients with vertical neutral wind:
            windp(k,i,lat) =
     |        bz(i,lat)*bdotu_m1(k,i,lat)*0.5*hdz(k,i,lat)
            windq(k,i,lat) = -bdotu(k,i,lat)*dvb(i,lat)
            windr(k,i,lat) =
     |        -bz(i,lat)*bdotu_p1(k,i,lat)*0.5*hdz(k,i,lat)
          enddo ! k=lev0,lev1
!
! Additions to Q coefficients, Sinks:
          do k=lev0,lev1
            p_coeff(k,i,lat) =
     |        diffp(k,i,lat)+driftp(k,i,lat)+windp(k,i,lat)
            q_coeff(k,i,lat) =
     |        diffq(k,i,lat)+driftq(k,i,lat)+windq(k,i,lat)-
     |        dtx2inv*nstep_sub-op_loss(k,i,lat)
            r_coeff(k,i,lat) =
     |        diffr(k,i,lat)+driftr(k,i,lat)+windr(k,i,lat)
          enddo ! k=lev0,lev1
!
! Upper boundary condition for O+:
          ubcb(i,lat) = bz(i,lat)**2*djint_tphdz0(lev1,i,lat) ! t3
          ubca(i,lat) = bz(i,lat)**2*djint_tphdz1(lev1,i,lat) ! t2
!
! Q = Q+B/A*R
          q_coeff(lev1-1,i,lat) = q_coeff(lev1-1,i,lat)+
     |      ubcb(i,lat)/ubca(i,lat)*r_coeff(lev1-1,i,lat)
!
! F = F-R/A*PHI
          ubcrhs(i,lat) = -opflux(i,lat)*
     |      r_coeff(lev1-1,i,lat)/ubca(i,lat)
          r_coeff(lev1-1,i,lat) = 0.
!
! Lower boundary condition N(O+) = Q/L:
          q_coeff(lev0,i,lat) = q_coeff(lev0,i,lat)-p_coeff(lev0,i,lat)
          lbcrhs(i,lat) = 2.*p_coeff(lev0,i,lat)*qop(lev0,i,lat)/
     |      (1.5*op_loss(lev0,i,lat)-0.5*op_loss(lev0+1,i,lat))
          p_coeff(lev0,i,lat) = 0.
        enddo ! i=lon0,lon1
!
!       call addfld('TR','TR: reduced temperature (0.5*(tn+ti))',' ',
!    |    tr(:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('DJ','DJ: Ambipolar Diffusion of O+',' ',
!    |    dj(:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('BDOTU' ,' ',' ',bdotu(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDZ',' ',' ',hdz(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('TP1',' ',' ',tp1(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('DJINT',' ',' ',djint(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('TPHDZ1',' ',' ',tphdz1(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('TPHDZ0',' ',' ',tphdz0(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('DIVBZ',' ',' ',divbz(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDZMBZ' ,' ',' ',hdzmbz(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('HDZPBZ' ,' ',' ',hdzpbz(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('XIOP2P',' ',' ',xiop2p(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('XIOP2D',' ',' ',xiop2d(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_LOSS_COEF',' ',' ',op_loss(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_PROD',' ',' ',op_prod(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('P_COEFF',' ',' ',
!    |    p_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('Q_COEFF',' ',' ',
!    |    q_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('R_COEFF',' ',' ',
!    |    r_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      end subroutine prep_oplus
!-----------------------------------------------------------------------
      subroutine smooth_oplus(optm1,optm1_smooth,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Shapiro smoother for O+
!
      use params_module,only: rp
      use cons_module,only: shapiro
      use input_module,only: nstep_sub
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  optm1      ! O+ at time n-1
!
! Output fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  optm1_smooth ! op at time n-1, with shapiro smoother (was s1)
!
! Local:
      integer :: k,i,lat
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0:lat1) ::
     |  optm1_smooth1 ! op at time n-1, with shapiro smoother (was s1)
!
! Shapiro smoother: optm1 is O+ at time n-1 (optm1_smooth was s1)
! optm1_smooth will be used in explicit terms below.
      do lat=lat0,lat1
        do i=lon0-2,lon1+2
          do k=lev0,lev1
            optm1_smooth1(k,i,lat) = optm1(k,i,lat)-
     |        shapiro/nstep_sub*
     |        (optm1(k,i,lat+2)+optm1(k,i,lat-2)-
     |        4.*(optm1(k,i,lat+1)+optm1(k,i,lat-1))+
     |        6.*optm1(k,i,lat))
          enddo ! k=lev0,lev1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0,lat1
!
! Boundary longitudes for optm1_smooth1 were obtained after first
! latitude scan above.
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            optm1_smooth(k,i,lat) = optm1_smooth1(k,i,lat)-
     |        shapiro/nstep_sub*
     |        (optm1_smooth1(k,i+2,lat)+optm1_smooth1(k,i-2,lat)-
     |        4.*(optm1_smooth1(k,i+1,lat)+optm1_smooth1(k,i-1,lat))+
     |        6.*optm1_smooth1(k,i,lat))
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
!     do lat=lat0,lat1
!       call addfld('OPTM1',' ',' ',optm1(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OPTM1_SM0' ,' ',' ',
!    |    optm1_smooth1(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OPTM1_SM1' ,' ',' ',
!    |    optm1_smooth(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
      end subroutine smooth_oplus
!-----------------------------------------------------------------------
      subroutine iterate_oplus(dvb,ubcrhs,lbcrhs,tp,dj,bdotu,
     |  optm1_smooth,op_prod,ui,vi,scht,p_coeff,q_coeff,r_coeff,
     |  op,opout,diffj,diffexp,vdotn_h,bdotdh_bvel,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Update O+ ion at 3d task subdomain.
! Outputs are opout, diffj, diffexp, vdotn_h, bdotdh_bvel,
! all other args are input.
!
      use params_module,only: nlonp4,dz,rp
      use cons_module,only: rmass_op,gask,grav,re,cs,dphi,dlamda,dtx2inv
      use input_module,only: nstep_sub
      use magfield_module,only: bz,bmod2 ! (nlonp4,-1:nlat+2)
      use mpi_module,only: mp_periodic_f3d,mp_polelats_f3d,
     |  mp_bndlats_f3d,mp_bndlons_f3d
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 2d task subdomain):
      real(rp),dimension(lon0-2:lon1+2,lat0-2:lat1+2),intent(in) ::
     |  dvb,        ! output of sub divb
     |  ubcrhs,lbcrhs
!
! Input fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  tp,          ! Plasma temperature (te+ti)
     |  dj,          ! diffusion coefficients (s13,s14,s15)
     |  bdotu,       ! was s7 (B.U)
     |  optm1_smooth,! op at time n-1, with shapiro smoother (was s1)
     |  op_prod,
     |  ui,vi,       ! zonal, and meridional ion velocities
     |  scht,        ! scale height
     |  p_coeff,q_coeff,r_coeff, ! coefficients for tridiagonal solver (s1,s2,s3)
     |  op           ! O+ ion
!
! Output fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  opout,    ! O+ output for next timestep
     |  diffj,    ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+) (s7,s8,s9)
     |  diffexp,vdotn_h,bdotdh_bvel
!
! Local:
      integer :: k,i,lat,ier,nlevs
      real(rp) :: mgr
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  uii,vii,dtpopdz,dbdotdh_opjdz,
     |  bdotdh_op,        ! (b(h)*del(h))*phi
     |  bdotdh_diff,      ! (b(h)*del(h))*phi
     |  bdzdvb_op,        ! was s7
     |  explicit          ! was s4
!
! Local fields at 3d subdomain (must be 3d to bridge latitude scans):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  tpop,bvel,opb,bdotdh_opj ! (b(h)*del(h))*phi
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,2) ::
     |  tmpf3d
!
#ifdef VT
!     code = 113 ; state = 'oplus' ; activity='ModelCode'
      call vtbegin(113,ier)
#endif
!
! Number of pressure levels (this will equal nlevp1):
      nlevs = lev1-lev0+1 ! for bndlons calls
!
      mgr = rmass_op*grav/gask
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1
!
! Plasma temperature times O+
            tpop(k,i,lat) = tp(k,i,lat)*op(k,i,lat)
!
! bvel @ j   = (B.U)*N(O+)      (J) (was s7)
            bvel(k,i,lat) = bdotu(k,i,lat)*op(k,i,lat)
          enddo ! k=lev0,lev1
!
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
          do k=lev0+1,lev1-1
            dtpopdz(k,i,lat) = (tpop(k+1,i,lat)-tpop(k-1,i,lat))/(2.*dz)
          enddo ! k=lev0+1,lev1-1
!
! Upper and lower boundaries:
          dtpopdz(lev1,i,lat) = (tpop(lev1,i,lat)-tpop(lev1-1,i,lat))/dz
          dtpopdz(lev0,i,lat) = (tpop(lev0+1,i,lat)-tpop(lev0,i,lat))/dz
!
          do k=lev0,lev1
            diffj(k,i,lat) = 2.*dtpopdz(k,i,lat)/scht(k,i,lat)+
     |        mgr*op(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      call bdotdh(bvel,bdotdh_bvel(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! bdotdh_op = (B(H).DEL(H))*(D/(H*DZ)*TP+M*G/R)*N(O+)
      call bdotdh(diffj,bdotdh_op,lev0,lev1,lon0,lon1,lat0,lat1)
!
! bdotdh_opj   = (B(H).DEL(H))*2.*TP*N(O+)      (J)
      call bdotdh(2.*tpop,bdotdh_opj(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            bdotdh_opj(k,i,lat) = bdotdh_opj(k,i,lat)*dj(k,i,lat)
            opb(k,i,lat) = op(k,i,lat)/bmod2(i,lat)**2
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            tmpf3d(k,i,lat,1) = bdotdh_opj(k,i,lat)
            tmpf3d(k,i,lat,2) = opb(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
      call mp_polelats_f3d(tmpf3d(:,lon0:lon1,:,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,2,(/1.,1./))
      call mp_bndlats_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2)
      call mp_bndlons_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2,0)
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1
            bdotdh_opj(k,i,lat) = tmpf3d(k,i,lat,1)
            opb(k,i,lat) = tmpf3d(k,i,lat,2)
          enddo ! k=lev0,lev1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! bdotdh_opj = (B(H).DEL(H))*D*(B(H).DEL(H))*2.*TP*N(O+)   (J)
! Note bdotdh_opj longitude dimension is lon-2:lon+2. bdotdh_diff is returned.
! (periodic points apparently not necessary for bdotdh_diff)
      call bdotdh(bdotdh_opj,bdotdh_diff,lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            uii(k,i,lat) = 0.5*(ui(k,i,lat)+ui(k+1,i,lat))
            vii(k,i,lat) = 0.5*(vi(k,i,lat)+vi(k+1,i,lat))
          enddo ! k=lev0,lev1-1
          uii(lev1,i,lat) = 1.5*ui(lev1,i,lat)-0.5*ui(lev1-1,i,lat)
          vii(lev1,i,lat) = 1.5*vi(lev1,i,lat)-0.5*vi(lev1-1,i,lat)
!
! bdzdvb_op = (BZ*D/(H*DZ)+DIV(*B))*S2
!
! Evaluates (bz*d/(h*dz)+divb)*phi
          do k=lev0+1,lev1-1
            dbdotdh_opjdz(k,i,lat) =
     |        (bdotdh_opj(k+1,i,lat)-bdotdh_opj(k-1,i,lat))/(2.*dz)
          enddo ! k=lev0+1,lev1-1
!
! Upper and lower boundaries:
          dbdotdh_opjdz(lev1,i,lat) =
     |      (bdotdh_opj(lev1,i,lat)-bdotdh_opj(lev1-1,i,lat))/dz
          dbdotdh_opjdz(lev0,i,lat) =
     |      (bdotdh_opj(lev0+1,i,lat)-bdotdh_opj(lev0,i,lat))/dz
!
          do k=lev0,lev1
            bdzdvb_op(k,i,lat) =
     |        bz(i,lat)*dbdotdh_opjdz(k,i,lat)/scht(k,i,lat)+
     |        dvb(i,lat)*bdotdh_opj(k,i,lat)
!
! bdotdh_op = (B(H).DEL(H))*(D/(H*DZ)*TP+M*G/R)*N(O+)
! then bdotdh_op = d*bz*bdotdh_op
            bdotdh_op(k,i,lat) =
     |        dj(k,i,lat)*bz(i,lat)*bdotdh_op(k,i,lat)
!
            diffexp(k,i,lat) =
     |        bdzdvb_op(k,i,lat)+bdotdh_diff(k,i,lat)+bdotdh_op(k,i,lat)
!
! Note if input flag DYNAMO<=0, then ui,vi,wi velocities will be zero.
            vdotn_h(k,i,lat) = bmod2(i,lat)**2/(2.*re)*
     |        (uii(k,i,lat)/(cs(lat)*dlamda)*
     |         (opb(k,i+1,lat)-opb(k,i-1,lat))+
     |         vii(k,i,lat)/dphi*
     |         (opb(k,i,lat+1)-opb(k,i,lat-1)))
!
! Collect explicit terms:
! Sum O+ at time n-1 to explicit terms: N(O+)/(2*DT) (N-1) (was s4)
! Add source term to RHS (explicit terms):
            explicit(k,i,lat) =
     |        vdotn_h(k,i,lat)+bdotdh_bvel(k,i,lat)-
     |        diffexp(k,i,lat)-op_prod(k,i,lat)-
     |        optm1_smooth(k,i,lat)*dtx2inv*nstep_sub
          enddo ! k=lev0,lev1
!
! Upper boundary condition for O+ (F = F-R/A*PHI):
          explicit(lev1-1,i,lat) = explicit(lev1-1,i,lat)-ubcrhs(i,lat)
!
! Lower boundary condition N(O+) = Q/L:
          explicit(lev0,i,lat) = explicit(lev0,i,lat)-lbcrhs(i,lat)
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Tridiagonal solver returns updated O+ in opout (all other args are input):
      do lat=lat0,lat1
        call trsolv(
     |    p_coeff(:,lon0:lon1,lat),
     |    q_coeff(:,lon0:lon1,lat),
     |    r_coeff(:,lon0:lon1,lat),
     |    explicit(:,lon0:lon1,lat),
     |    opout(:,lon0:lon1,lat),
     |    lev0,lev1,lev0,lev1-1,lon0,lon1,nlonp4,lat,0)
!
! Upper boundary:
        do i=lon0,lon1
          opout(lev1,i,lat) = 2.*opout(lev1-1,i,lat)-opout(lev1-2,i,lat)
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Periodic points for outputs:
      call mp_periodic_f3d(opout(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1,1)
!
!     do lat=lat0,lat1
!       call addfld('TPJ','TP after times OP',' ',tpop(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('BVEL_J',' ',' ',bvel(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('DIFFJ','DIFFJ after diffus',' ',
!    |    diffj(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('BDOTDH_1',' ',' ',bdotdh_op(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('BDOTDH_2',' ',' ',bdotdh_opj(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('BDOT_DIF',' ',' ',bdotdh_diff(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('BDZDVB',' ',' ',bdzdvb_op(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('EXPLIC',' ',' ',explicit(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('OP_SOLV',' ',' ',opout(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
#ifdef VT
!     code = 113 ; state = 'oplus' ; activity='ModelCode'
      call vtend(113,ier)
#endif
!
      end subroutine iterate_oplus
!-----------------------------------------------------------------------
      subroutine post_oplus(op,optm1,opout,optm1out,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Post-processing time smoothing.
! This is called after filter_op
!
      use params_module,only: zpmid,spval,rp  ! for O+ minimum
      use cons_module,only: rtd,dtsmooth,dtsmooth_div2
      use input_module,only: opfloor,oprate,oplev,oplatwidth
      use magfield_module,only: rlatm
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  op,         ! O+ ion
     |  optm1       ! O+ at time n-1
!
! Output fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) ::
     |  opout     ! O+ output for next timestep
!
! Output fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  optm1out  ! O+ output for time n-1
!
! Local:
      integer :: k,i,lat
!
! For O+ "smooth floor". The floor is a
! spatial gaussian based on the O+ minimum (opmin below).
!
! 1/21/16 btf: Increased O+ minimum.
! (This alleviates numerical instability that can develop
!  near the equator at the top of the model in some cases,
!  e.g., it allows the 2.5-deg June solstice solar max
!  benchmark to run with step=30 rather than step=20, and
!  may improve stability in other cases as well.)
!
! 2024/08 Haonan Wu: opmin is determined from input parameters
! OPFLOOR, OPRATE, OPLEV, and OPLATWIDTH to be more flexible
      real(rp) :: opmin
!
! 12/4/14 btf: Enforce O+ minimum.
! Opfloor is Stan's "smooth floor" (product of two Gaussians,
!   dependent on latitude and pressure level):
!
! 2024/08 Haonan Wu:
!   Change the altitude distribution of O+ minimum from Gaussian to logistic,
!   also change the latitude distribution based on magnetic latitudes
!   instead of geographic latitudes
      if (opfloor/=0. .and.
     |    oprate/=0. .and.
     |    oplev/=spval .and.
     |    oplatwidth/=0.) then
        do lat=lat0,lat1
          do i=lon0,lon1
            do k=lev0,lev1
              opmin = opfloor*exp(-(rlatm(i,lat)*rtd/oplatwidth)**2)/
     |          (1+exp(-oprate*(zpmid(k)-oplev)))
              if (opout(k,i,lat) < opmin) opout(k,i,lat) = opmin
            enddo ! k=lev0,lev1
          enddo ! i=lon0,lon1
        enddo ! lat=lat0,lat1
      endif
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
!
! Ensure global non-negative O+:
! This should already been satisfied if O+ minimum is applied above
            if (opout(k,i,lat) < 0.) opout(k,i,lat) = 0.
!
! Time smoothing:
!
! optm1out(k,i,lat): New O+ at current latitude and time n-1.
! op(k,i,lat)      : O+ at current latitude and time.
! optm1(k,i,lat)   : O+ at current latitude and time n-1.
! opout(k,i,lat)   : New O+ at current latitude and time.
            optm1out(k,i,lat) = dtsmooth*op(k,i,lat)+
     |        dtsmooth_div2*(optm1(k,i,lat)+opout(k,i,lat))
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Save outputs on secondary history for diagnostics:
!     do lat=lat0,lat1
!       call addfld('OPOUT',' ',' ',opout(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OPOUTM1',' ',' ',optm1out(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
      end subroutine post_oplus
!-----------------------------------------------------------------------
      subroutine bdotdh(phij,ans,lev0,lev1,lon0,lon1,lat0,lat1)
!
      use params_module,only: rp
      use cons_module,only: re,dphi,dlamda,cs
      use magfield_module,only: bx,by
!
! Evaluates ans = (b(h)*del(h))*phi
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) :: phij
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) ::
     |  ans
!
! Local:
      integer :: k,i,lat
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0:lat1) :: dphidx,dphidy
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Boundary longitudes i-1 and i+1 must have been set before this routine is
! called (e.g., call mp_bndlons_f3d).
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            dphidx(k,i,lat) =
     |        (phij(k,i+1,lat)-phij(k,i-1,lat))/(2.*dlamda)
            dphidy(k,i,lat) =
     |        (phij(k,i,lat+1)-phij(k,i,lat-1))/(2.*dphi)
            ans(k,i,lat) =
     |        (bx(i,lat)*dphidx(k,i,lat)/cs(lat)+
     |         by(i,lat)*dphidy(k,i,lat))/re
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
      end subroutine bdotdh
!-----------------------------------------------------------------------
      subroutine calc_terms(xnmbar,diffj,Fe,Fn,
     |  opout,optm1_smooth,op_prod,op_loss,
     |  diffexp,diffp,diffq,diffr,
     |  vdotn_h,driftp,driftq,driftr,
     |  bdotdh_bvel,windp,windq,windr,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Post-processing diagnostic term analysis.
! This is called only at the last sub-cycling step
!
      use params_module,only: rp
      use cons_module,only: rmass_op,dtx2inv
      use input_module,only: nstep_sub
      use magfield_module,only: dipmag,sndec,csdec
      use addfld_module,only: addfld
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  xnmbar,diffj,
     |  opout,        ! O+ output for next timestep
     |  optm1_smooth, ! op at time n-1, with shapiro smoother (was s1)
     |  op_prod,op_loss,
     |  diffexp,diffp,diffq,diffr,
     |  vdotn_h,driftp,driftq,driftr,
     |  bdotdh_bvel,windp,windq,windr
!
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) :: Fe,Fn
!
! Local:
      integer :: k,i,lat
      real(rp),dimension(lon0:lon1,lat0:lat1) :: proj_e,proj_n
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  wd,dopdt,op_loss_out,diffsum,driftsum,windsum
!
      do lat=lat0,lat1
        do i=lon0,lon1
          proj_e(i,lat) =
     |      sin(dipmag(i,lat))*cos(dipmag(i,lat))*sndec(i,lat)
          proj_n(i,lat) =
     |      sin(dipmag(i,lat))*cos(dipmag(i,lat))*csdec(i,lat)
!
          do k=lev0,lev1
            wd(k,i,lat) = 3.53*1.42E6*rmass_op/
     |        xnmbar(k,i,lat)*diffj(k,i,lat)
            Fe(k,i,lat) = wd(k,i,lat)*proj_e(i,lat)
            Fn(k,i,lat) = wd(k,i,lat)*proj_n(i,lat)
!
            dopdt(k,i,lat) = dtx2inv*nstep_sub*
     |        (opout(k,i,lat)-optm1_smooth(k,i,lat))
            op_loss_out(k,i,lat) = op_loss(k,i,lat)*opout(k,i,lat)
          enddo ! k=lev0,lev1
!
          do k=lev0+1,lev1-1
            diffsum (k,i,lat) = diffexp(k,i,lat)+
     |        diffp (k,i,lat)*opout(k-1,i,lat)+
     |        diffq (k,i,lat)*opout(k  ,i,lat)+
     |        diffr (k,i,lat)*opout(k+1,i,lat)
            driftsum(k,i,lat) = -vdotn_h(k,i,lat)+
     |        driftp(k,i,lat)*opout(k-1,i,lat)+
     |        driftq(k,i,lat)*opout(k  ,i,lat)+
     |        driftr(k,i,lat)*opout(k+1,i,lat)
            windsum (k,i,lat) = -bdotdh_bvel(k,i,lat)+
     |        windp (k,i,lat)*opout(k-1,i,lat)+
     |        windq (k,i,lat)*opout(k  ,i,lat)+
     |        windr (k,i,lat)*opout(k+1,i,lat)
          enddo ! k=lev0+1,lev1-1
        enddo ! i=lon0,lon1
!
!       call addfld('PARDRAG_U','','',Fe(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('PARDRAG_V','','',Fn(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('DOPDT','O+ changing rate','cm-3s-1',
     |    dopdt(:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('OP_PROD','O+ production','cm-3s-1',
     |    op_prod(:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('OP_LOSS','O+ loss','cm-3s-1',
     |    -op_loss_out(:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('OP_DIFF',
     |    'O+ transport due to ambipolar diffusion','cm-3s-1',
     |    diffsum(lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_DRIFT',
     |    'O+ transport due to ion drift','cm-3s-1',
     |     driftsum(lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_WIND',
     |    'O+ transport due to neutral wind','cm-3s-1',
     |    windsum(lev0+1:lev1-1,:,lat),
     |    'lev',lev0+1,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      end subroutine calc_terms
!-----------------------------------------------------------------------
      subroutine filter_op(opout,lev0,lev1,lon0,lon1,lat0,lat1,name)
!
! Filter updated O+. This is called once per timestep.
!
      use params_module,only: nlonp4,rp
      use mpi_module,only: mytidi,mp_gatherlons_f3d,mp_scatterlons_f3d
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real(rp),dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(inout) ::
     |  opout
      character(len=*),intent(in) :: name
!
! Local:
      integer :: k,i,lat,ier
      real(rp),dimension(nlonp4,lev0:lev1) :: op_ik
      real(rp),dimension(lev0:lev1,nlonp4,lat0:lat1) :: op_kij
!
#ifdef VT
!     code = 124 ; state = 'filter_op' ; activity='Filtering'
      call vtbegin(124,ier)
#endif
!
! Define lons in op_kij from current task subdomain opout:
      op_kij = 0.
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            op_kij(k,i,lat) = opout(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each lat-row) for global fft. (i.e., tasks with mytidi==0
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(op_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)
!
! Only leftmost tasks at each lat-row of tasks does the global filtering:
        if (mytidi==0) then
!
! Loop through subdomain latitudes:
          latscan: do lat=lat0,lat1
!
! Define 2d array at all longitudes for filter:
            do i=1,nlonp4
              do k=lev0,lev1
                op_ik(i,k) = op_kij(k,i,lat)
              enddo ! k=lev0,lev1
            enddo ! i=1,nlonp4
!
! Do the filtering:
            call ringfilter(op_ik,lev0,lev1,3,name,lat)
!
! Return filtered array to op_kij:
            do i=1,nlonp4
              do k=lev0,lev1
                op_kij(k,i,lat) = op_ik(i,k)
              enddo ! k=lev0,lev1
            enddo ! i=1,nlonp4
          enddo latscan ! lat=lat0,lat1
        endif ! mytidi==0
!
! Now leftmost task at each lat-row must redistribute filtered data
! back to other tasks in the lat-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(op_kij,lev0,lev1,lon0,lon1,lat0,lat1,1,
     |  name)
!
! Return filtered array to opout at task subdomain:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            opout(k,i,lat) = op_kij(k,i,lat)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
#ifdef VT
!     code = 124 ; state = 'filter_op' ; activity='Filtering'
      call vtend(124,ier)
#endif
!
      end subroutine filter_op
!-----------------------------------------------------------------------
      end module oplus_module

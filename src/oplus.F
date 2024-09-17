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
      use params_module,only: zpmid,spval  ! for O+ minimum
      use cons_module,only: rtd
      use input_module,only: nstep_sub,opfloor,oprate,oplev,oplatwidth
      use magfield_module,only: rlatm
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
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
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
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  opout,    ! O+ output for next timestep
     |  optm1out, ! O+ output for time n-1
     |  xiop2p,xiop2d,Fe,Fn
!
! Local:
      integer :: k,i,lat,nlevs,istep
      real,dimension(lon0-2:lon1+2,lat0-2:lat1+2) :: dvb,ubcrhs,lbcrhs
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  dj,tp,bdotu,op_prod,op_loss,diffp,diffq,diffr,
     |  driftp,driftq,driftr,windp,windq,windr,
     |  p_coeff,q_coeff,r_coeff,optm1_smooth,
     |  op_sub,optm1_sub,opout_sub,optm1out_sub,
     |  diffj,diffexp,vdotn_h,bdotdh_bvel
!
! For O+ "smooth floor". The floor is a
! spatial gaussian based on the O+ minimum (opmin below).
!
! 1/21/16 btf: O+ minimum increased from 1000 to 3000.
! (This alleviates numerical instability that can develop
!  near the equator at the top of the model in some cases,
!  e.g., it allows the 2.5-deg June solstice solar max
!  benchmark to run with step=30 rather than step=20, and
!  may improve stability in other cases as well.)
!
      real :: opmin
!
! Number of pressure levels (this will equal nlevp1):
      nlevs = lev1-lev0+1 ! for bndlons calls
!
      call prep_oplus(tn,te,ti,
     |  o2,o1,he,n2,n2d,ne,un,vn,w,wi,xnmbar,scht,
     |  dvb,ubcrhs,lbcrhs,dj,tp,bdotu,xiop2p,xiop2d,
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
        call iterate_oplus(dvb,ubcrhs,lbcrhs,dj,tp,bdotu,
     |    optm1_smooth,op_prod,ui,vi,scht,p_coeff,q_coeff,r_coeff,
     |    op_sub,opout_sub,diffj,diffexp,vdotn_h,bdotdh_bvel,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
!
        if (istep == nstep_sub) then
          call calc_terms(xnmbar,diffj,Fe,Fn,
     |      opout_sub,optm1_smooth,op_prod,op_loss,
     |      diffexp,diffp,diffq,diffr,
     |      vdotn_h,driftp,driftq,driftr,
     |      bdotdh_bvel,windp,windq,windr,
     |      lev0,lev1,lon0,lon1,lat0,lat1)
        endif
!
! Filter updated O+:
!
        call filter_op(opout_sub(:,lon0:lon1,lat0:lat1),
     |    lev0,lev1,lon0,lon1,lat0,lat1,'OPLUS')
!
!       do lat=lat0,lat1
!         call addfld('OP_FILT',' ',' ',
!    |      opout_sub(lev0:lev1-1,lon0:lon1,lat),
!    |      'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       enddo ! lat=lat0,lat1
!
! 12/4/14 btf: Enforce O+ minimum.
! Opfloor is Stan's "smooth floor" (product of two Gaussians,
!   dependent on latitude and pressure level) (opmin=1000.0):
!
! 2024/08 Haonan Wu:
!   Change the altitude distribution of O+ minimum from Gaussian to logistic,
!   also change the latitude distribution based on magnetic latitudes
!   instead of geographic latitudes
!
        if (opfloor/=0. .and. oprate/=0. .and.
     |      oplev/=spval .and. oplatwidth/=0.) then
          do lat=lat0,lat1
            do i=lon0,lon1
              do k=lev0,lev1-1
                opmin = opfloor*exp(-(rlatm(i,lat)*rtd/oplatwidth)**2)/
     |            (1+exp(-oprate*(zpmid(k)-oplev)))
                if (opout_sub(k,i,lat) < opmin) then
                  opout_sub(k,i,lat) = opmin
                endif ! opout_sub < opmin
              enddo ! k=lev0,lev1-1
            enddo ! i=lon0,lon1
          enddo ! lat=lat0,lat1
        endif
!
        call mp_polelats_f3d(opout_sub(:,lon0:lon1,:),
     |    lev0,lev1,lon0,lon1,lat0,lat1,1,(/1./))
        call mp_bndlats_f3d(opout_sub,nlevs,lon0,lon1,lat0,lat1,1)
        call mp_bndlons_f3d(opout_sub,nlevs,lon0,lon1,lat0,lat1,1,0)
!
        call post_oplus(op_sub,optm1_sub,opout_sub,optm1out_sub,
     |    lev0,lev1,lon0,lon1,lat0,lat1)
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
     |  dvb,ubcrhs,lbcrhs,dj,tp,bdotu,xiop2p,xiop2d,
     |  op_prod,op_loss,diffp,diffq,diffr,
     |  driftp,driftq,driftr,windp,windq,windr,
     |  p_coeff,q_coeff,r_coeff,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate terms that are not changing across O+ sub-cycling
!
      use params_module,only: nlonp4,dz,zpmid,spval
      use cons_module,only: rmass_op,gask,grav,re,cs,dphi,dlamda,
     |  dtx2inv,pi,dtr,rmassinv_o2,rmassinv_o1,
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
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
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
      real,dimension(lon0-2:lon1+2,lat0-2:lat1+2),intent(out) ::
     |  dvb,        ! output of sub divb
     |  ubcrhs,lbcrhs
!
! Output fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  dj,          ! diffusion coefficients (s13,s14,s15)
     |  tp,          ! Plasma temperature 0.5*(te+ti)
     |  bdotu,       ! was s7 (B.U)
     |  xiop2p,xiop2d,
     |  op_prod,op_loss,diffp,diffq,diffr,
     |  driftp,driftq,driftr,windp,windq,windr,
     |  p_coeff,q_coeff,r_coeff ! coefficients for tridiagonal solver (s1,s2,s3)
!
! Local:
      real,parameter ::
     |  phid =  2.0e8,
     |  phin = -2.0e8,
     |  ppolar = 0.
      integer :: k,i,lat,lonbeg,lonend,nlevs
      real :: gmr,a,fed,fen,dbxdx,dbydy,tr,djmin,
     |  o2_cm3,o1_cm3,n2_cm3,he_cm3,n2d_cm3,
     |  ubca, ubcb ! O+ upper boundary condition (were t2,t3)
      real,dimension(lev0:lev1) :: bdotu1
      real,dimension(lon0:lon1,lat0:lat1) ::
     |  opflux,    ! upward number flux of O+ (returned by sub oplus_flux) (t7)
     |  dvb_bz
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  hdz,              ! was s15
     |  tphdz1,tphdz0,    ! were s13,s12 (using gmr)
     |  djint,            ! was s11
     |  divbz,            ! was s7 (DIV(B)+(DH*D*BZ)/(D*BZ)
     |  hdzmbz,hdzpbz,    ! were s10,s9
     |  tp1,              ! 0.5*(te+ti)
     |  hdzi,bdotdh_djbz,
     |  djint_tphdz0,djint_tphdz1,
     |  djint_tphdz0_1,djint_tphdz1_1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) :: djbz
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,2) :: tmpf3d
!
! Number of pressure levels (this will equal nlevp1):
      nlevs = lev1-lev0+1 ! for bndlons calls
!
! Save inputs to secondary history file:
!
!     do lat=lat0,lat1
!       call addfld('TE_OP',' ',' ',te(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('TI_OP',' ',' ',ti(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('N2D_OP',' ',' ',n2d(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('NE_OP',' ',' ',ne(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('WI_OP',' ',' ',wi(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
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
          fed = phid*a
          fen = phin*a
          if (chi(i,lat) >= 100.*dtr) then
            opflux(i,lat) = fen
          elseif (chi(i,lat) <= 80.*dtr) then
            opflux(i,lat) = fed
          else
            opflux(i,lat) = .5*(fed+fen+(fed-fen)*cos(chi(i,lat)*9.))
          endif
!
! Add ppolar if magnetic latitude >= 60 degrees:
          if (abs(rlatm(i,lat)) >= pi/3.)
     |      opflux(i,lat) = opflux(i,lat)+ppolar
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!     call addfld('OPFLUX',' ',' ',opflux,
!    |  'lon',lon0,lon1,'lat',lat0,lat1,0)
!
! Evaluate divergence of B, the unit magnetic field vector.
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp4-2
!
      dvb = 0.
      do lat=lat0,lat1
        do i=lonbeg,lonend
          dbxdx = (bx(i+1,lat)-bx(i-1,lat))/(2.*dlamda)
          dbydy =
     |      (cs(lat+1)*by(i,lat+1)-cs(lat-1)*by(i,lat-1))/(2.*dphi)
          dvb(i,lat) = ((dbxdx+dbydy)/cs(lat)+2.*bz(i,lat))/re
        enddo ! i=lonbeg,lonend
      enddo ! lat=lat0,lat1
      call mp_periodic_f2d(dvb(lon0:lon1,lat0:lat1),
     |  lon0,lon1,lat0,lat1,1)
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
! Set reduced temperature (average of tn and ti)
!
! 1/2/16 btf:
!   This was previously set incorrectly as tp = 0.5*(te+ti))
!   (see also plasma temp tp). Correcting this resulted in a
!   40% change in ambipolar diffusion DJ, but only 10% or so
!   change in O+ and Ne, and very small differences in NMF2,HMF2,TEC.
!   Wenbin feels there should be larger changes in these fields
!   as a consequence of fixing this bug.
!
            tr = 0.5*(tn(k,i,lat)+ti(k,i,lat))
            o2_cm3 = xnmbar(k,i,lat)*o2(k,i,lat)*rmassinv_o2
            o1_cm3 = xnmbar(k,i,lat)*o1(k,i,lat)*rmassinv_o1
            he_cm3 = xnmbar(k,i,lat)*he(k,i,lat)*rmassinv_he
            n2_cm3 = xnmbar(k,i,lat)*n2(k,i,lat)*rmassinv_n2
!
! 8/28/13 btf: Use n2=1-o2-o-he, and include O+/He collision rate from Wenbin:
!
            dj(k,i,lat) = 1.42E17/(18.1*o2_cm3+3.6*he_cm3+18.6*n2_cm3+
     |        o1_cm3*sqrt(tr)*(1.-0.064*log10(tr))**2*colfac)
          enddo ! k=lev0,lev1-1
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
!
      if (opdiffcap/=0. .and. opdiffrate/=0. .and. opdifflev/=spval)
     |  then
        do lat=lat0-2,lat1+2
          do i=lon0-2,lon1+2
            do k=lev0,lev1-1
              djmin = opdiffcap/(1+exp(opdiffrate*(zpmid(k)-opdifflev)))
              if (dj(k,i,lat) > djmin) dj(k,i,lat) = djmin
            enddo ! k=lev0,lev1-1
          enddo ! i=lon0-2,lon1+2
        enddo ! lat=lat0-2,lat1+2
      endif

!     do lat=lat0,lat1
!       call addfld('DJ','DJ: Ambipolar Diffusion of O+',' ',
!    |    dj(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! Plasma temperature:
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            tp(k,i,lat) = 0.5*(te(k,i,lat)+ti(k,i,lat))
          enddo ! k=lev0,lev1-1
          tp(lev1,i,lat) = 2.*tp(lev1-1,i,lat)-tp(lev1-2,i,lat)
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! bdotu = B.U (s7)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            bdotu(k,i,lat) =
     |        bx(i,lat)*un(k,i,lat)+by(i,lat)*vn(k,i,lat)+
     |        bz(i,lat)*scht(k,i,lat)*0.5*(w(k,i,lat)+w(k+1,i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!       call addfld('BDOTU' ,' ',' ',bdotu(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            djbz(k,i,lat) = dj(k,i,lat)*bz(i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Fill up halo points of bdotu, djbz (used in bdotdh calls)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            tmpf3d(k,i,lat,1) = bdotu(k,i,lat)
            tmpf3d(k,i,lat,2) = djbz(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
      call mp_polelats_f3d(tmpf3d(:,lon0:lon1,:,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,2,(/1.,1./))
      call mp_bndlats_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2)
      call mp_bndlons_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2,0)
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            bdotu(k,i,lat) = tmpf3d(k,i,lat,1)
            djbz(k,i,lat) = tmpf3d(k,i,lat,2)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
      do lat=lat0,lat1
        do i=lon0,lon1
          dvb_bz(i,lat) = dvb(i,lat)/bz(i,lat)
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            hdz(k,i,lat) = 1./(scht(k,i,lat)*dz)
          enddo ! k=lev0,lev1
          do k=lev0,lev1-1
            tp1(k+1,i,lat) = tp(k,i,lat)
          enddo ! k=lev0,lev1-1
          tp1(lev0,i,lat) = 2.*tp(lev0,i,lat)-tp(lev0+1,i,lat)
        enddo ! i=lon0,lon1

!       call addfld('TP1',' ',' ',tp1(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('HDZ',' ',' ',hdz(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-2
            hdzi(k+1,i,lat) = 0.5*(hdz(k,i,lat)+hdz(k+1,i,lat))
          enddo ! k=lev0,lev1-2
!
! Upper and lower boundaries:
          hdzi(lev0,i,lat) = 1.5*hdz(lev0  ,i,lat)-0.5*hdz(lev0+1,i,lat)
          hdzi(lev1,i,lat) = 1.5*hdz(lev1-1,i,lat)-0.5*hdz(lev1-2,i,lat)
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! gmr = G*M(O+)/(2.*R)
      gmr = grav*rmass_op/(2.*gask)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            tphdz1(k,i,lat) = 2.*tp (k,i,lat)*hdzi(k,i,lat)+gmr ! s13
            tphdz0(k,i,lat) = 2.*tp1(k,i,lat)*hdzi(k,i,lat)-gmr ! s12
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!       call addfld('TPHDZ1',' ',' ',tphdz1(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('TPHDZ0',' ',' ',tphdz0(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! djint = dj at interfaces:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-2
            djint(k+1,i,lat) = 0.5*(dj(k,i,lat)+dj(k+1,i,lat))
          enddo ! k=lev0,lev1-2
          djint(lev0,i,lat) = 1.5*dj(lev0  ,i,lat)-0.5*dj(lev0+1,i,lat)
          djint(lev1,i,lat) = 1.5*dj(lev1-1,i,lat)-0.5*dj(lev1-2,i,lat)
        enddo ! i=lon0,lon1
!       call addfld('DJINT' ,' ',' ',djint(:,:,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      call bdotdh(djbz,bdotdh_djbz,lev0,lev1,lon0,lon1,lat0,lat1)
!
! divbz = (DIV(B)+(DH*D*BZ)/(D*BZ) (was s7)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            divbz(k,i,lat) = dvb_bz(i,lat)+
     |        bdotdh_djbz(k,i,lat)/(djbz(k,i,lat)*bz(i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!       call addfld('DIVBZ' ,' ',' ',divbz(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! hdzmbz = (1./(H*DZ)-(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s10)
! hdzpbz = (1./(H*DZ)+(DIV(B)+DH*D*BZ/(D*BZ))/(2*BZ))*BZ**2 (was s9 )
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            hdzmbz(k,i,lat) =
     |        (hdz(k,i,lat)-0.5*divbz(k,i,lat))*bz(i,lat)**2
            hdzpbz(k,i,lat) =
     |        (hdz(k,i,lat)+0.5*divbz(k,i,lat))*bz(i,lat)**2
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('HDZMBZ' ,' ',' ',hdzmbz(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('HDZPBZ' ,' ',' ',hdzpbz(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1
            djint_tphdz0(k,i,lat) = djint(k,i,lat)*tphdz0(k,i,lat)
            djint_tphdz1(k,i,lat) = djint(k,i,lat)*tphdz1(k,i,lat)
          enddo ! k=lev0,lev1
          do k=lev0,lev1-1
            djint_tphdz0_1(k,i,lat) = djint_tphdz0(k+1,i,lat)
            djint_tphdz1_1(k,i,lat) = djint_tphdz1(k+1,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Begin coefficients p_coeff, q_coeff, r_coeff  (s1,s2,s3)
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            diffp(k,i,lat) =   hdzmbz(k,i,lat)*djint_tphdz0  (k,i,lat)
            diffq(k,i,lat) = -(hdzpbz(k,i,lat)*djint_tphdz0_1(k,i,lat)+
     |                         hdzmbz(k,i,lat)*djint_tphdz1  (k,i,lat))
            diffr(k,i,lat) =   hdzpbz(k,i,lat)*djint_tphdz1_1(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('P_COEFF0',' ',' ',
!    |    diffp(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('Q_COEFF0',' ',' ',
!    |    diffq(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('R_COEFF0',' ',' ',
!    |    diffr(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Continue coefficients with vertical ion velocity:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            driftp(k,i,lat) =
     |        0.5*(wi(k,i,lat)+wi(k+1,i,lat))*0.5*hdz(k,i,lat)
            driftq(k,i,lat) = -0.5*(wi(k,i,lat)+wi(k+1,i,lat))*6./re
            driftr(k,i,lat) =
     |        -0.5*(wi(k,i,lat)+wi(k+1,i,lat))*0.5*hdz(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('P_COEFF1',' ',' ',
!    |    driftp(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('Q_COEFF1',' ',' ',
!    |    driftq(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('R_COEFF1',' ',' ',
!    |    driftr(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Continue coefficients with vertical neutral wind:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-2
            bdotu1(k+1) = bdotu(k,i,lat)
          enddo ! k=lev0,lev1-2
          bdotu1(lev0) = 2.*bdotu(lev0,i,lat)-bdotu(lev0+1,i,lat)
          do k=lev0,lev1-1
            windp(k,i,lat) = bz(i,lat)*bdotu1(k)*0.5*hdz(k,i,lat)
          enddo ! k=lev0,lev1-1
!
          do k=lev0,lev1-1
            windq(k,i,lat) = -bdotu(k,i,lat)*dvb(i,lat)
          enddo ! k=lev0,lev1-1
!
          do k=lev0,lev1-2
            bdotu1(k) = bdotu(k+1,i,lat)
          enddo ! k=lev0,lev1-2
          bdotu1(lev1-1) = 2.*bdotu(lev1-1,i,lat)-bdotu(lev1-2,i,lat)
          do k=lev0,lev1-1
            windr(k,i,lat) = -bz(i,lat)*bdotu1(k)*0.5*hdz(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('P_COEFF2',' ',' ',
!    |    windp(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('Q_COEFF2',' ',' ',
!    |    windq(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('R_COEFF2',' ',' ',
!    |    windr(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Additions to Q coefficients:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            p_coeff(k,i,lat) =
     |        diffp(k,i,lat)+driftp(k,i,lat)+windp(k,i,lat)
            q_coeff(k,i,lat) =
     |        diffq(k,i,lat)+driftq(k,i,lat)+windq(k,i,lat)-
     |        dtx2inv*nstep_sub
            r_coeff(k,i,lat) =
     |        diffr(k,i,lat)+driftr(k,i,lat)+windr(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Upper boundary condition for O+:
      do lat=lat0,lat1
        do i=lon0,lon1
          ubcb = -bz(i,lat)**2*djint_tphdz0(lev1,i,lat) ! t3
          ubca = -bz(i,lat)**2*djint_tphdz1(lev1,i,lat) ! t2
!
! Q = Q+B/A*R
          q_coeff(lev1-1,i,lat) = q_coeff(lev1-1,i,lat)+
     |      ubcb/ubca*r_coeff(lev1-1,i,lat)
!
! F = F -R/A*PHI
          ubcrhs(i,lat) = opflux(i,lat)*r_coeff(lev1-1,i,lat)/ubca
          r_coeff(lev1-1,i,lat) = 0.
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1

!     do lat=lat0,lat1
!       call addfld('P_COEFF3',' ',' ',
!    |    p_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('Q_COEFF3',' ',' ',
!    |    q_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('R_COEFF3',' ',' ',
!    |    r_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

!       call addfld('QOP2P_OP',' ',' ',qop2p(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('QOP2D_OP',' ',' ',qop2d(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RK20',' ',' ',rk20(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('RK25',' ',' ',rk25(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! Sources and sinks (xiop2p and xiop2d are outputs):
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            o2_cm3 = xnmbar(k,i,lat)*o2(k,i,lat)*rmassinv_o2
            o1_cm3 = xnmbar(k,i,lat)*o1(k,i,lat)*rmassinv_o1
            n2_cm3 = xnmbar(k,i,lat)*n2(k,i,lat)*rmassinv_n2
            n2d_cm3 = xnmbar(k,i,lat)*n2d(k,i,lat)*rmassinv_n2d
!
            xiop2p(k,i,lat) =
     |        0.5*(qop2p(k,i,lat)+qop2p(k+1,i,lat))/
     |        ((rk16+rk17)*n2_cm3+rk18*o1_cm3+
     |        (rk19(k,i,lat)+rk20(k,i,lat))*ne(k,i,lat)+rk21+rk22)
!
            xiop2d(k,i,lat) =
     |        (0.5*(qop2d(k,i,lat)+qop2d(k+1,i,lat))+
     |        (rk20(k,i,lat)*ne(k,i,lat)+rk22)*xiop2p(k,i,lat))/
     |        (rk23*n2_cm3+rk24*o1_cm3+rk26*o2_cm3+
     |        rk25(k,i,lat)*ne(k,i,lat)+rk27)
!
            op_loss(k,i,lat) =
     |        rk1(k,i,lat)*o2_cm3+rk2(k,i,lat)*n2_cm3+rk10*n2d_cm3
!
            q_coeff(k,i,lat) = q_coeff(k,i,lat)-op_loss(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('XIOP2P',' ',' ',xiop2p(lev0:lev1,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('XIOP2D',' ',' ',xiop2d(lev0:lev1,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_LOSS_COEF',' ',' ',
!    |    op_loss(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('OP_QOP',' ',' ',qop(lev0:lev1,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_NE' ,' ',' ',ne(lev0:lev1,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_O1' ,' ',' ',o1(lev0:lev1,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_TN' ,' ',' ',tn(lev0:lev1,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_RK19' ,' ',' ',rk19(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('OP_RK25' ,' ',' ',rk25(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Add source term to RHS (explicit terms):
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            op_prod(k,i,lat) =
     |        0.5*(qop(k,i,lat)+qop(k+1,i,lat))+
     |        (rk19(k,i,lat)*ne(k,i,lat)+rk21)*xiop2p(k,i,lat)+
     |        (rk25(k,i,lat)*ne(k,i,lat)+rk27)*xiop2d(k,i,lat)+
     |        (rk18*xiop2p(k,i,lat)+rk24*xiop2d(k,i,lat))*o1_cm3
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
! Lower boundary condition N(O+) = Q/L:
      do lat=lat0,lat1
        do i=lon0,lon1
          q_coeff(lev0,i,lat) = q_coeff(lev0,i,lat)-p_coeff(lev0,i,lat)
          lbcrhs(i,lat) = 2.*p_coeff(lev0,i,lat)*qop(lev0,i,lat)/
     |      (1.5*op_loss(lev0,i,lat)-0.5*op_loss(lev0+1,i,lat))
          p_coeff(lev0,i,lat) = 0.
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1

!     do lat=lat0,lat1
!       call addfld('P_COEFF',' ',' ',
!    |    p_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('Q_COEFF',' ',' ',
!    |    q_coeff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('R_COEFF',' ',' ',
!    |    r_coeff(lev0:lev1-2,lon0:lon1,lat),
!    |    'lev',lev0,lev1-2,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1

      end subroutine prep_oplus
!-----------------------------------------------------------------------
      subroutine smooth_oplus(optm1,optm1_smooth,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Shapiro smoother for O+
!
      use cons_module,only: shapiro
      use input_module,only: nstep_sub
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  optm1      ! O+ at time n-1
!
! Output fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  optm1_smooth ! op at time n-1, with shapiro smoother (was s1)
!
! Local:
      integer :: k,i,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0:lat1) ::
     |  optm1_smooth1 ! op at time n-1, with shapiro smoother (was s1)
!
! Save inputs to secondary history file:
!
!     do lat=lat0,lat1
!       call addfld('OPTM1',' ',' ',optm1(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! Shapiro smoother: optm1 is O+ at time n-1 (optm1_smooth was s1)
! optm1_smooth will be used in explicit terms below.
      do lat=lat0,lat1
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            optm1_smooth1(k,i,lat) = optm1(k,i,lat)-
     |        shapiro/nstep_sub*
     |        (optm1(k,i,lat+2)+optm1(k,i,lat-2)-
     |        4.*(optm1(k,i,lat+1)+optm1(k,i,lat-1))+
     |        6.*optm1(k,i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
!       call addfld('OPTM1_SM0' ,' ',' ',
!    |    optm1_smooth1(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Boundary longitudes for optm1_smooth1 were obtained after first
! latitude scan above.
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            optm1_smooth(k,i,lat) = optm1_smooth1(k,i,lat)-
     |        shapiro/nstep_sub*
     |        (optm1_smooth1(k,i+2,lat)+optm1_smooth1(k,i-2,lat)-
     |        4.*(optm1_smooth1(k,i+1,lat)+optm1_smooth1(k,i-1,lat))+
     |        6.*optm1_smooth1(k,i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!       call addfld('OPTM1_SM1' ,' ',' ',
!    |    optm1_smooth(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      end subroutine smooth_oplus
!-----------------------------------------------------------------------
      subroutine iterate_oplus(dvb,ubcrhs,lbcrhs,dj,tp,bdotu,
     |  optm1_smooth,op_prod,ui,vi,scht,p_coeff,q_coeff,r_coeff,
     |  op,opout,diffj,diffexp,vdotn_h,bdotdh_bvel,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Update O+ ion at 3d task subdomain.
! Outputs are opout, diffj, diffexp, vdotn_h, bdotdh_bvel,
! all other args are input.
!
      use params_module,only: nlonp4,dz
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
      real,dimension(lon0-2:lon1+2,lat0-2:lat1+2),intent(in) ::
     |  dvb,        ! output of sub divb
     |  ubcrhs,lbcrhs
!
! Input fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  dj,          ! diffusion coefficients (s13,s14,s15)
     |  tp,          ! Plasma temperature (te+ti)
     |  bdotu,       ! was s7 (B.U)
     |  optm1_smooth,! op at time n-1, with shapiro smoother (was s1)
     |  op_prod,
     |  ui,vi,       ! zonal, and meridional ion velocities
     |  scht,        ! scale height
     |  p_coeff,q_coeff,r_coeff, ! coefficients for tridiagonal solver (s1,s2,s3)
     |  op           ! O+ ion
!
! Output fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  opout,    ! O+ output for next timestep
     |  diffj,    ! (D/(H*DZ)*2.*TP+M*G/R)*N(O+) (s7,s8,s9)
     |  diffexp,vdotn_h,bdotdh_bvel
!
! Local:
      integer :: k,i,lat,ier,nlevs
      real :: mgr
      real,dimension(lev0:lev1) :: dtpop,dbdotdh_opj
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  bdzdvb_op,        ! was s7
     |  explicit,         ! was s4
     |  bdotdh_op,        ! (b(h)*del(h))*phi
     |  bdotdh_diff       ! (b(h)*del(h))*phi
!
! Local fields at 3d subdomain (must be 3d to bridge latitude scans):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  bvel,tpop,op_bmod2,bdotdh_opj ! (b(h)*del(h))*phi
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2,2) :: tmpf3d
!
#ifdef VT
!     code = 113 ; state = 'oplus' ; activity='ModelCode'
      call vtbegin(113,ier)
#endif
!
! Number of pressure levels (this will equal nlevp1):
      nlevs = lev1-lev0+1 ! for bndlons calls
!
! Save inputs to secondary history file:
!
!     do lat=lat0,lat1
!       call addfld('OP_OPLUS',' ',' ',op(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('UI_OP',' ',' ',ui(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VI_OP',' ',' ',vi(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! bvel @ j   = (B.U)*N(O+)      (J) (was s7)
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            bvel(k,i,lat) = bdotu(k,i,lat)*op(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!     do lat=lat0,lat1
!       call addfld('BVEL_J',' ',' ',bvel(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! Plasma temperature times O+
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            tpop(k,i,lat) = tp(k,i,lat)*op(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!     do lat=lat0,lat1
!       call addfld('TPJ','TP after times OP',' ',
!    |    tpop(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! Evaluates ans = (d/(h*dz)*tp+m*g/r)*en
!
      mgr = rmass_op*grav/gask
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0+1,lev1-2
            dtpop(k) = (tpop(k+1,i,lat)-tpop(k-1,i,lat))/2.
          enddo ! k=lev0+1,lev1-2
!
! Upper and lower boundaries:
          dtpop(lev1-1) = tpop(lev1-1,i,lat)-tpop(lev1-2,i,lat)
          dtpop(lev0) = tpop(lev0+1,i,lat)-tpop(lev0,i,lat)
!
          do k=lev0,lev1-1
            diffj(k,i,lat) = 1./(scht(k,i,lat)*dz)*2.*dtpop(k)+
     |        mgr*op(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2

!     do lat=lat0,lat1
!       call addfld('DIFFJ','DIFFJ after diffus',' ',
!    |    diffj(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! bdotdh_op = (B(H).DEL(H))*(D/(H*DZ)*TP+M*G/R)*N(O+)
! then bdotdh_op = d*bz*bdotdh_op
!
      call bdotdh(diffj,bdotdh_op,lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            bdotdh_op(k,i,lat) =
     |        dj(k,i,lat)*bz(i,lat)*bdotdh_op(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('BDOTDH_1',' ',' ',
!    |    bdotdh_op(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! bdotdh_opj   = (B(H).DEL(H))*2.*TP*N(O+)      (J)
!
      call bdotdh(2.*tpop,bdotdh_opj(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            bdotdh_opj(k,i,lat) = bdotdh_opj(k,i,lat)*dj(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('BDOTDH_2',' ',' ',
!    |    bdotdh_opj(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            op_bmod2(k,i,lat) = op(k,i,lat)/bmod2(i,lat)**2
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            tmpf3d(k,i,lat,1) = bdotdh_opj(k,i,lat)
            tmpf3d(k,i,lat,2) = op_bmod2(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
      call mp_polelats_f3d(tmpf3d(:,lon0:lon1,:,:),
     |  lev0,lev1,lon0,lon1,lat0,lat1,2,(/1.,1./))
      call mp_bndlats_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2)
      call mp_bndlons_f3d(tmpf3d,nlevs,lon0,lon1,lat0,lat1,2,0)
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            bdotdh_opj(k,i,lat) = tmpf3d(k,i,lat,1)
            op_bmod2(k,i,lat) = tmpf3d(k,i,lat,2)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! bdotdh_opj = (B(H).DEL(H))*D*(B(H).DEL(H))*2.*TP*N(O+)   (J)
! Note bdotdh_opj longitude dimension is lon-2:lon+2. bdotdh_diff
!   is returned. (periodic points apparently not necessary for
!   bdotdh_diff)
!
      call bdotdh(bdotdh_opj,bdotdh_diff,lev0,lev1,lon0,lon1,lat0,lat1)

!     do lat=lat0,lat1
!       call addfld('BDOT_J'  ,' ',' ',
!    |    bdotdh_opj(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('BDOT_DIF',' ',' ',
!    |    bdotdh_diff(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
! bdzdvb_op = (BZ*D/(H*DZ)+DIV(*B))*S2
!
! Evaluates (bz*d/(h*dz)+divb)*phi
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0+1,lev1-2
            dbdotdh_opj(k) =
     |        (bdotdh_opj(k+1,i,lat)-bdotdh_opj(k-1,i,lat))/2.
          enddo ! k=lev0+1,lev1-2
!
! Upper and lower boundaries:
          dbdotdh_opj(lev1-1) =
     |      bdotdh_opj(lev1-1,i,lat)-bdotdh_opj(lev1-2,i,lat)
          dbdotdh_opj(lev0) =
     |      bdotdh_opj(lev0+1,i,lat)-bdotdh_opj(lev0,i,lat)
!
          do k=lev0,lev1-1
            bdzdvb_op(k,i,lat) =
     |        bz(i,lat)/(scht(k,i,lat)*dz)*dbdotdh_opj(k)+
     |        dvb(i,lat)*bdotdh_opj(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('BDZDVB',' ',' ',bdzdvb_op(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Collect explicit terms:
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            diffexp(k,i,lat) =
     |        bdzdvb_op(k,i,lat)+bdotdh_diff(k,i,lat)+bdotdh_op(k,i,lat)
            explicit(k,i,lat) = -diffexp(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('EXPLIC0',' ',' ',explicit(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1

!     do lat=lat0,lat1
!       call addfld('UI_VEL',' ',' ',ui(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('VI_VEL',' ',' ',vi(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
!
      call bdotdh(bvel,bdotdh_bvel(:,lon0:lon1,lat0:lat1),
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Note if input flag DYNAMO<=0, then ui,vi,wi velocities will be zero.
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            vdotn_h(k,i,lat) = bmod2(i,lat)**2/(2.*re)*
     |        (1./(cs(lat)*dlamda)*
     |        0.5*(ui(k,i,lat)+ui(k+1,i,lat))*
     |        (op_bmod2(k,i+1,lat)-op_bmod2(k,i-1,lat))+
!
     |        1./dphi*
     |        0.5*(vi(k,i,lat)+vi(k+1,i,lat))*
     |        (op_bmod2(k,i,lat+1)-op_bmod2(k,i,lat-1)))
!
            explicit(k,i,lat) = explicit(k,i,lat)+
     |        vdotn_h(k,i,lat)+bdotdh_bvel(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1

!       call addfld('EXPLIC1',' ',' ',explicit(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Sum O+ at time n-1 to explicit terms: N(O+)/(2*DT) (N-1) (was s4)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            explicit(k,i,lat) = explicit(k,i,lat)-
     |        optm1_smooth(k,i,lat)*dtx2inv*nstep_sub
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1

!     call addfld('OPTM1_SM1' ,' ',' ',
!    |  optm1_smooth(lev0:lev1-1,lon0:lon1,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('EXPLIC2' ,' ',' ',explicit(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Upper boundary condition for O+:
      do lat=lat0,lat1
        do i=lon0,lon1
!
! F = F -R/A*PHI
          explicit(lev1-1,i,lat) = explicit(lev1-1,i,lat)-ubcrhs(i,lat)
        enddo ! i=lon0,lon1

!       call addfld('EXPLIC3',' ',' ',explicit(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Add source term to RHS (explicit terms):
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            explicit(k,i,lat) = explicit(k,i,lat)-op_prod(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!       call addfld('EXPLIC4',' ',' ',explicit(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Lower boundary condition N(O+) = Q/L:
      do lat=lat0,lat1
        do i=lon0,lon1
          explicit(lev0,i,lat) = explicit(lev0,i,lat)-lbcrhs(i,lat)
        enddo ! i=lon0,lon1

!       call addfld('EXPLIC5',' ',' ',explicit(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
! Tridiagonal solver returns updated O+ in opout (all other args are input):
!     subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax,lat,
!    |  idebug)
!
      do lat=lat0,lat1
        call trsolv(
     |    p_coeff(:,lon0:lon1,lat),
     |    q_coeff(:,lon0:lon1,lat),
     |    r_coeff(:,lon0:lon1,lat),
     |    explicit(:,lon0:lon1,lat),
     |    opout(:,lon0:lon1,lat),
     |    lev0,lev1,lev0,lev1-1,lon0,lon1,nlonp4,lat,0)
!
!       call addfld('OP_SOLV',' ',' ',opout(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      enddo ! lat=lat0,lat1
!
! Periodic points for outputs:
        call mp_periodic_f3d(opout(:,lon0:lon1,lat0:lat1),
     |    lev0,lev1,lon0,lon1,lat0,lat1,1)
!
#ifdef VT
!     code = 113 ; state = 'oplus' ; activity='ModelCode'
      call vtend(113,ier)
#endif
      end subroutine iterate_oplus
!-----------------------------------------------------------------------
      subroutine post_oplus(op,optm1,opout,optm1out,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! Post-processing time smoothing.
! This is called after filter_op
!
      use cons_module,only: dtsmooth,dtsmooth_div2
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1, ! first,last pressure  indices for current task (bot->top)
     |  lon0,lon1, ! first,last longitude indices for current task (W->E)
     |  lat0,lat1  ! first,last latitude  indices for current task (S->N)
!
! Input fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  op,         ! O+ ion
     |  optm1       ! O+ at time n-1
!
! Output fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(inout) ::
     |  opout     ! O+ output for next timestep
!
! Output fields (full 3d task subdomain):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  optm1out  ! O+ output for time n-1
!
! Local:
      integer :: k,i,lat
!
! Time smoothing:
!
! optm1out(k,i,lat): New O+ at current latitude and time n-1.
! op(k,i,lat)      : O+ at current latitude and time.
! optm1(k,i,lat)   : O+ at current latitude and time n-1.
! opout(k,i,lat)   : New O+ at current latitude and time.
!
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            optm1out(k,i,lat) = dtsmooth*op(k,i,lat)+
     |        dtsmooth_div2*(optm1(k,i,lat)+opout(k,i,lat))
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! Upper boundary:
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          opout(lev1,i,lat) = 0.
          optm1out(lev1,i,lat) = 0.
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! Insure global non-negative O+:
      do lat=lat0-2,lat1+2
        do i=lon0-2,lon1+2
          do k=lev0,lev1-1
            if (opout(k,i,lat)    < 0.) opout(k,i,lat)    = 0.
            if (optm1out(k,i,lat) < 0.) optm1out(k,i,lat) = 0.
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0-2,lon1+2
      enddo ! lat=lat0-2,lat1+2
!
! Save outputs on secondary history for diagnostics:
!     do lat=lat0,lat1
!       call addfld('OPOUT',' ',' ',opout(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('OPOUTM1',' ',' ',
!    |    optm1out(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     enddo ! lat=lat0,lat1
      end subroutine post_oplus
!-----------------------------------------------------------------------
      subroutine bdotdh(phij,ans,lev0,lev1,lon0,lon1,lat0,lat1)
      use cons_module,only: re,dphi,dlamda,cs
      use magfield_module,only: bx,by
!
! Evaluates ans = (b(h)*del(h))*phi
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) :: phij
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(out) :: ans
!
! Local:
      integer :: k,i,lat
      real :: dphidx,dphidy
!
! Note phij longitude dimension is lon0-2:lon1+2 (only i-1 and i+1 are used).
! Boundary longitudes i-1 and i+1 must have been set before this routine is
! called (e.g., call mp_bndlons_f3d).
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            dphidx = (phij(k,i+1,lat)-phij(k,i-1,lat))/(2.*dlamda)
            dphidy = (phij(k,i,lat+1)-phij(k,i,lat-1))/(2.*dphi)
            ans(k,i,lat) = 1./re*
     |        (bx(i,lat)*dphidx/cs(lat)+by(i,lat)*dphidy)
          enddo ! k=lev0,lev1-1
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
      use cons_module,only: dtx2inv
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
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) ::
     |  xnmbar,diffj,
     |  opout,         ! O+ output for next timestep
     |  optm1_smooth, ! op at time n-1, with shapiro smoother (was s1)
     |  op_prod,op_loss,
     |  diffexp,diffp,diffq,diffr,
     |  vdotn_h,driftp,driftq,driftr,
     |  bdotdh_bvel,windp,windq,windr
!
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) :: Fe,Fn
!
! Local:
      integer :: k,i,lat
      real :: wd
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  dopdt,op_loss_out,diffsum,driftsum,windsum
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            wd = 16*3.53*1.42E6/xnmbar(k,i,lat)*diffj(k,i,lat)*
     |        sin(dipmag(i,lat))*cos(dipmag(i,lat))
            Fe(k,i,lat) = wd*sndec(i,lat)
            Fn(k,i,lat) = wd*csdec(i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1
!       call addfld('PARDRAG_U','','',Fe(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!       call addfld('PARDRAG_V','','',Fn(lev0:lev1-1,lon0:lon1,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        do i=lon0,lon1
          do k=lev0,lev1-1
            dopdt(k,i,lat) = dtx2inv*nstep_sub*
     |        (opout(k,i,lat)-optm1_smooth(k,i,lat))
            op_loss_out(k,i,lat) = op_loss(k,i,lat)*opout(k,i,lat)
          enddo ! k=lev0,lev1-1
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
        call addfld('DOPDT','O+ changing rate','cm-3s-1',
     |    dopdt(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_PROD','O+ production','cm-3s-1',
     |    op_prod(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
        call addfld('OP_LOSS','O+ loss','cm-3s-1',
     |    -op_loss_out(lev0:lev1-1,:,lat),
     |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
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
! Filter updated O+. This is called from outside latitude loop,
! i.e., once per timestep.
!
      use params_module,only: nlonp4
      use mpi_module,only: mytidi,mp_gatherlons_f3d,mp_scatterlons_f3d
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1),intent(inout) ::
     |  opout
      character(len=*),intent(in) :: name
!
! Local:
      integer :: k,i,lat,ier
      real,dimension(nlonp4,lev0:lev1) :: op_ik
      real,dimension(lev0:lev1,nlonp4,lat0:lat1) :: op_kij
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
      end subroutine filter_op
!-----------------------------------------------------------------------
      end module oplus_module

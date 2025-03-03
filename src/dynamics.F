!
      subroutine dynamics(nstep,istep)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Call dynamics layer of the model for current time step.
! Routines that do not have to do message-passing are generally called
!   from inside a latitude scan in this routine. Routines that do have
!   to make mpi calls are called once from this routine, and latitude
!   scan(s) are inside the called routine.
!
! Include the option of setting the lower boundary by using GSWM output
! this affects subroutine dt and duv
!
! am 6/06: change call to qjoule_ti to include ion and neutral density
!          factors, include vertical component
!          change call to qjoule_tn to include vertical velocity
!
      use params_module,only: nlev,nlevp1
      use cons_module,only: rmass_o2,rmass_o1,rmass_n2,step=>dt
      use fields_module
      use mpi_module,only: lon0,lon1,lat0,lat1,distribute_1d,mp_close,
     |  mp_polelats_f3d,mp_bndlats_f3d,mp_bndlons_f3d
      use chemrates_module,only: chemrates_tdep,
     |  rk16,rk17,rk18,rk21,rk22,rk23,rk24,rk26,rk27,
     |  rk1,rk2,rk3,rk19,rk20,rk25,ra1,ra2,ra3,beta9
      use qrj_module,only: qrj,qtotal,qop,qo2p,qn2p,qnp,qnop,qop2p,qop2d
      use chapman_module,only: chapman
      use oplus_module,only: oplus
      use input_module,only: iaurora=>aurora,iamie,nstep_sub
      use aurora_module,only: aurora,cusp,drizl,alfa,flux,eflux
      use n2d_module,only: comp_n2d,minor_n2d
      use n4s_module,only: comp_n4s,minor_n4s
      use no_module,only: comp_no,minor_no
      use ar_module,only: comp_ar,minor_ar
      use hdif_module,only: hdif3,hdif_bndlatlons
      use timing_module,only: timer,timing
      use addfld_module,only: addfld
      use qjoule,only: qjoule_tn,qjoule_ti
      use diags_module,only: mkdiag_CUSP,mkdiag_DRIZZLE,mkdiag_ALFA,
     |  mkdiag_NFLUX,mkdiag_EFLUX,mkdiag_TEC,mkdiag_HNMF2,
     |  mkdiag_QJOULE_INTEG,mkdiag_UI,mkdiag_VI,mkdiag_WI
      use dyndiag_module,only: dyndiag
      use ionden_module,only: calculate_prodloss,ionden_local

      implicit none
!
! Args:
      integer,intent(in) :: nstep,istep
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: k,i,lat,ier
      integer :: i0,i1,nk,nkm1,nlats,k0,k1
      logical,parameter :: debug=.false.  ! add prints to stdout
!
! If diags is set true, the following fields will be written to
! secondary histories:
!
! SECFLDS = 'CP','KT','KM','W_OMEGA','UI_ExB','VI_ExB','WI_ExB',
!           'SCO2','SCO1','SCN2','QTOTAL','QOP','QO2P','QN2P','QNP',
!           'QNOP','OPLUS','XIOP2P','XIOP2D','NP_diag','N2P_diag',
!           'NOP_diag','O2P_diag','ELDEN','LXX','LYY','LXY','LYX',
!           'LAM1','SIGMAPED','SIGMAHAL','CMP_N4S','CMP_NO','QJI_TI',
!           'T_ELEC','T_ION','COOL_IMPLICIT','COOL_EXPLICIT',
!           'QJI_TN','TN_UPD','UN_UPD','VN_UPD','O2_UPD','O1_UPD'
!
      logical,parameter :: diags=.true. ! put fields on sec hist
      real :: time0,time1,time1_qrj,time1_cmpminor,time1_cmpmajor,
     |  o2_cm3,o1_cm3,n2_cm3,
     |  beta9m,qopm,qo2pm,qnpm,qn2pm,qnopm,qop2pm,qop2dm,
     |  o2p_prod,np_prod,n2p_prod,nop_prod,op2p_prod,op2d_prod,
     |  o2p_loss,np_loss,n2p_loss,nop_loss,op2p_loss,op2d_loss
      real,dimension(lon0:lon1,lat0:lat1) ::
     |  tec ! diagnostic total electron content
      real,dimension(nlevp1) :: nem
      real,dimension(nlevp1,lond0:lond1,latd0:latd1) :: op_prod,op_loss
      logical,external :: time2print
!
      if (debug) write(6,"(/,'Enter dynamics.')")
      i0=lon0 ; i1=lon1 ; nk=nlevp1 ; nkm1=nk-1 ! for addfld
      k0=1 ; k1=nlevp1
      nlats = lat1-lat0+1
!
! Calculate specific heat and molecular viscosity.
! f4d itp fields are input, f3d fields are output.
!
      do lat=lat0,lat1
        call cpktkm(
     |    tn(levd0,lond0,lat,itp),  ! 4d input
     |    o2(levd0,lond0,lat,itp),  ! 4d input
     |    o1(levd0,lond0,lat,itp),  ! 4d input
     |    he(levd0,lond0,lat,itp),  ! 4d input
     |    n2(levd0,lond0,lat),      ! 3d input
     |    mbar(levd0,lond0,lat,itp),! 4d input
     |    cp(levd0,lond0,lat),      ! 3d output
     |    kt(levd0,lond0,lat),      ! 3d output
     |    km(levd0,lond0,lat),      ! 3d output
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after cpktkm')")
      if (diags) then
        do lat=lat0,lat1
          call addfld('CP','CP',' ',cp(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('KT','KT',' ',kt(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('KM','KM',' ',km(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Calculate omega for vertical velocity:
! (latitude scan is inside swdot)
!
      if (debug) write(6,"('dynamics before swdot')")
      call swdot(
     |  un(levd0,lond0,latd0,itp), ! un input
     |  vc(levd0,lond0,latd0,itp), ! vc input
     |  z (levd0,lond0,latd0,itc), ! input for WN diag
     |  w (levd0,lond0,latd0,itc), ! omega output
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after swdot')")
      if (diags) then
        do lat=lat0,lat1
          call addfld('W_OMEGA','VERTICAL MOTION (DZp/DT)',
     |    's-1',w(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Calculate ExB ion drift velocities from electric field:
!
      call ionvel(
     |  z (levd0,lond0,latd0,itc), ! geopotential input
     |  ui(levd0,lond0,latd0),     ! U ion velocity output
     |  vi(levd0,lond0,latd0),     ! V ion velocity output
     |  wi(levd0,lond0,latd0),     ! W ion velocity output
     |  Etot(levd0,lond0,latd0),   ! total eletric field
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after ionvel')")
!
! mkdiags_UI,VI,WI are called from ionvel.
!     if (diags) then
!       do lat=lat0,lat1
!         call addfld('UI_ExB','Zonal ExB Velocity','cm/s',
!    |      ui(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
!         call addfld('VI_ExB','Meridional ExB Velocity','cm/s',
!    |      vi(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
!         call addfld('WI_ExB','Vertical ExB Velocity','cm/s',
!    |      wi(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
!       enddo
!     endif
!
! Calculate column densities and line integrals of o2, o, and n2.
! Prognostics z,tn,o2,o1,barm are input, diagnostics vo2,vo1,vn2,
!   sco2,sco1,scn2 are output.
!
      do lat=lat0,lat1
        call chapman(
     |    z   (levd0,lond0,lat,itc), ! updated Z from addiag
     |    tn  (levd0,lond0,lat,itp), ! 4d in
     |    o2  (levd0,lond0,lat,itp), ! 4d in
     |    o1  (levd0,lond0,lat,itp), ! 4d in
     |    n2  (levd0,lond0,lat),     ! 3d in
     |    barm(levd0,lond0,lat,itp), ! 4d in
     |    vo2 (levd0,lond0,lat),     ! 3d out
     |    vo1 (levd0,lond0,lat),     ! 3d out
     |    vn2 (levd0,lond0,lat),     ! 3d out
     |    sco2(levd0,lond0,lat),     ! 3d out
     |    sco1(levd0,lond0,lat),     ! 3d out
     |    scn2(levd0,lond0,lat),     ! 3d out
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after chapman')")
      if (diags) then
!
! Post-processors cannot read 1.e80, so change these to
! 1.e36 if saving on secondary histories (then change back
! for rest of model):
        do lat=lat0,lat1
          call addfld('SCO2','Chapman slant-column O2',' ',
     |      sco2(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('SCO1','Chapman slant-column O1',' ',
     |      sco1(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('SCN2','Chapman slant-column N2',' ',
     |      scn2(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Calculate temperature dependent reaction rates (chemrates_module): 
!
      do lat=lat0,lat1
        call chemrates_tdep(
     |    tn  (levd0,lond0,lat,itp),
     |    te  (levd0,lond0,lat,itp),
     |    ti  (levd0,lond0,lat,itp),
     |    sco2(levd0,lond0,lat),
     |    vo2 (levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
      enddo
      if (debug) write(6,"('dynamics after chemrates_tdep')")
!
! Latitude scan for qrj, qinite, and aurora.
!
#ifdef VT
!     code = 122 ; state = 'ions' ; activity='ModelCode'
      call vtbegin(122,ier)
#endif
      time1_qrj = 0.
!
      do lat=lat0,lat1
!
! Calculate ionization and heating rates:
!     subroutine qrj(sco2,sco1,scn2,tn,no,o2,o1,n4s,xnmbari,
!    |  lev0,lev1,lon0,lon1,lat)
!
!       if (mod(istep,nstepqrj)==0) then
          if (debug) write(6,"('dynamics call qrj: lat=',i3)") lat
          call timer(time0,time1,'QRJ',0,0) ! start qrj timing for current lat
          call qrj(
     |      sco2(levd0,lond0,lat),
     |      sco1(levd0,lond0,lat),
     |      scn2(levd0,lond0,lat),
     |      tn  (levd0,lond0,lat,itp),
     |      no  (levd0,lond0,lat,itp),
     |      o2  (levd0,lond0,lat,itp),
     |      o1  (levd0,lond0,lat,itp),
     |      he  (levd0,lond0,lat,itp),
     |      n2  (levd0,lond0,lat),
     |      xnmbari(levd0,lond0,lat,itp),
     |      1,nlevp1,lon0,lon1,lat)
          call timer(time0,time1,'QRJ',1,0) ! end qrj timing for current lat
          time1_qrj = time1_qrj+time1
          if (debug) write(6,"('dynamics after qrj: lat=',i3)") lat
!       endif
        if (diags) then
          call addfld('QTOTAL','Total Heating',' ',
     |      qtotal(:,i0:i1,lat),'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('QOP'   ,'QO+ from QRJ',' ',
     |      qop(:,i0:i1,lat),'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('QO2P'  ,'QO2+ from QRJ',' ',
     |      qo2p(:,i0:i1,lat),'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('QN2P'  ,'QN2+ from QRJ',' ',
     |      qn2p(:,i0:i1,lat),'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('QNP'   ,'QN+ from QRJ',' ',
     |      qnp(:,i0:i1,lat),'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('QNOP'  ,'QNO+ from QRJ',' ',
     |      qnop(:,i0:i1,lat),'ilev',k0,k1,'lon',i0,i1,lat)
        endif
!
! Calculate background (nightime) ionization:
        call qinite(
     |    o2     (levd0,lond0,lat,itp),
     |    o1     (levd0,lond0,lat,itp),
     |    he     (levd0,lond0,lat,itp),
     |    n2     (levd0,lond0,lat),
     |    no     (levd0,lond0,lat,itp),
     |    xnmbari(levd0,lond0,lat,itp),
     |    vo2    (levd0,lond0,lat),
     |    vo1    (levd0,lond0,lat),
     |    vn2    (levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after qinite: lat=',i3)") lat
      enddo ! lat=lat0,lat1 for qrj, qinite
!
! Add aurora particle precip and ionization:
      if (iaurora > 0) then
        call aurora(
     |    o2  (levd0,lond0,latd0,itp),
     |    o1  (levd0,lond0,latd0,itp),
     |    he  (levd0,lond0,latd0,itp),
     |    n2  (levd0,lond0,latd0),
     |    xnmbar(levd0,lond0,latd0,itp),
     |    scht(levd0,lond0,latd0,itp),
     |    1,nlevp1,lon0,lon1,lat0,lat1)
        if (debug) write(6,"('dynamics after aurora')")
!
! Save auroral diagnostics (see diags.F):
!
        call mkdiag_CUSP('CUSP',cusp(lon0:lon1,lat0:lat1),
     |    lon0,lon1,lat0,lat1)
        call mkdiag_DRIZZLE('DRIZZLE',drizl(lon0:lon1,lat0:lat1),
     |    lon0,lon1,lat0,lat1)
        call mkdiag_ALFA('ALFA',alfa(lon0:lon1,lat0:lat1),
     |    lon0,lon1,lat0,lat1)
        call mkdiag_NFLUX('NFLUX',flux(lon0:lon1,lat0:lat1),
     |    lon0,lon1,lat0,lat1)
        call mkdiag_EFLUX('EFLUX',eflux(lon0:lon1,lat0:lat1),
     |    lon0,lon1,lat0,lat1)
!
      endif
      if (debug) write(6,"('dynamics after qrj,qinite,aurora')")
!
! Report timing from above latitude loop:
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in QRJ =          ',
     |    f12.2,' Dynamics: step ',i5)") time1_qrj,istep
!
! Save ion drifts if requested:
      call mkdiag_UI('UI_ExB',ui(:,lon0:lon1,lat0:lat1),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      call mkdiag_VI('VI_ExB',vi(:,lon0:lon1,lat0:lat1),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      call mkdiag_WI('WI_ExB',wi(:,lon0:lon1,lat0:lat1),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
!
#ifdef VT
!     code = 122 ; state = 'ions' ; activity='ModelCode'
      call vtend(122,ier)
#endif
!
! Advance O+:
! Oplus makes several calls to mp_bndlons_f3d, mp_periodic_f3d,
!   mp_gatherlons_f3d, and mp_scatterlons_f3.
! Full task subdomains (including ghost cells) are passed.
!
      call timer(time0,time1,'OPLUS',0,0) ! start oplus timing
!
! Make latitudes j=-1,0 from j=2,1, and j=nlat+1,nlat+2 from nlat-1,nlat for n2, 
! then define 2d halos points.  Note that mp_geo_halos_f3d and mp_polelats_f3d 
! do not define "outside corner" halo points, but this should be ok since 
! cross-derivatives are not performed.
!
      call mp_polelats_f3d(n2(:,lon0:lon1,:),   ! 3rd dim is lat0-2:lat1+2
     |  1,nlevp1,lon0,lon1,lat0,lat1,1,(/1./))  ! last arg means no change in polesign
      call mp_bndlats_f3d(n2,nk,lon0,lon1,lat0,lat1,1)
      call mp_bndlons_f3d(n2,nk,lon0,lon1,lat0,lat1,1,0)
!
      if (istep == 1) then
        do concurrent (lat=lat0:lat1, i=lon0:lon1, k=1:nlev)
          o2_cm3 = xnmbar(k,i,lat,itp)*o2(k,i,lat,itp)/rmass_o2
          o1_cm3 = xnmbar(k,i,lat,itp)*o1(k,i,lat,itp)/rmass_o1
          n2_cm3 = xnmbar(k,i,lat,itp)*n2(k,i,lat)/rmass_n2
          nem(k) = sqrt(ne(k,i,lat,itp)*ne(k+1,i,lat,itp))
          qop2pm = (qop2p(k,i,lat)+qop2p(k+1,i,lat))/2
          qop2dm = (qop2d(k,i,lat)+qop2d(k+1,i,lat))/2
          op2p(k,i,lat,itc) = qop2pm/
     |      ((rk16+rk17)*n2_cm3+rk18*o1_cm3+
     |       (rk19(k,i,lat)+rk20(k,i,lat))*nem(k)+rk21+rk22)
          op2d(k,i,lat,itc) = (qop2dm+
     |       rk20(k,i,lat)*op2p(k,i,lat,itc)*nem(k)+
     |       rk22*op2p(k,i,lat,itc))/
     |      (rk23*n2_cm3+rk24*o1_cm3+rk25(k,i,lat)*nem(k)+
     |       rk26*o2_cm3+rk27)
        enddo
        do lat=lat0,lat1
          call elden(
     |      xnmbar(levd0,lond0,lat,itp),
     |      op    (levd0,lond0,lat,itp), ! O+ at current time step
     |      op    (levd0,lond0,lat,itp), ! O+ updated by sub oplus
     |      o2    (levd0,lond0,lat,itp),
     |      o1    (levd0,lond0,lat,itp),
     |      n2    (levd0,lond0,lat),
     |      n2d   (levd0,lond0,lat,itp),
     |      no    (levd0,lond0,lat,itp),
     |      n4s   (levd0,lond0,lat,itp),
     |      op2p  (levd0,lond0,lat,itc),  ! from oplus
     |      op2d  (levd0,lond0,lat,itc),  ! from oplus
     |      np    (levd0,lond0,lat,itc),  ! N+  output
     |      n2p   (levd0,lond0,lat,itc),  ! N2+ output
     |      nop   (levd0,lond0,lat,itc),  ! NO+ output
     |      o2p   (levd0,lond0,lat,itc),  ! O2+ output
     |      ne    (levd0,lond0,lat,itc),  ! electron density output
     |      1,nlevp1,lon0,lon1,lat)
        enddo
        do concurrent (lat=lat0:lat1, i=lon0:lon1, k=1:nlev)
          np(k,i,lat,itp) = np(k,i,lat,itc)
          n2p(k,i,lat,itp) = n2p(k,i,lat,itc)
          nop(k,i,lat,itp) = nop(k,i,lat,itc)
          op2p(k,i,lat,itp) = op2p(k,i,lat,itc)
          op2d(k,i,lat,itp) = op2d(k,i,lat,itc)
        enddo
      endif
      do concurrent (lat=lat0:lat1, i=lon0:lon1, k=1:nlev)
        nem(k) = sqrt(ne(k,i,lat,itp)*ne(k+1,i,lat,itp))
        beta9m = (beta9(k,i,lat)+beta9(k+1,i,lat))/2
        qopm = (qop(k,i,lat)+qop(k+1,i,lat))/2
        qo2pm = (qo2p(k,i,lat)+qo2p(k+1,i,lat))/2
        qnpm = (qnp(k,i,lat)+qnp(k+1,i,lat))/2
        qn2pm = (qn2p(k,i,lat)+qn2p(k+1,i,lat))/2
        qnopm = (qnop(k,i,lat)+qnop(k+1,i,lat))/2
        qop2pm = (qop2p(k,i,lat)+qop2p(k+1,i,lat))/2
        qop2dm = (qop2d(k,i,lat)+qop2d(k+1,i,lat))/2
        call calculate_prodloss(
     |    xnmbar(k,i,lat,itp),
     |    o2    (k,i,lat,itp),
     |    o1    (k,i,lat,itp),
     |    n2    (k,i,lat),
     |    n2d   (k,i,lat,itp),
     |    n4s   (k,i,lat,itp),
     |    no    (k,i,lat,itp),
     |    nem   (k),
     |    op    (k,i,lat,itp),
     |    o2p   (k,i,lat,itp),
     |    np    (k,i,lat,itp),
     |    n2p   (k,i,lat,itp),
     |    op2p  (k,i,lat,itp),
     |    op2d  (k,i,lat,itp),
     |    rk1   (k,i,lat),
     |    rk2   (k,i,lat),
     |    rk3   (k,i,lat),
     |    rk19  (k,i,lat),
     |    rk20  (k,i,lat),
     |    rk25  (k,i,lat),
     |    ra1   (k,i,lat),
     |    ra2   (k,i,lat),
     |    ra3   (k,i,lat),
     |    beta9m,qopm,qo2pm,qnpm,qn2pm,qnopm,qop2pm,qop2dm,
     |    op_prod(k,i,lat),o2p_prod,np_prod,n2p_prod,
     |    nop_prod,op2p_prod,op2d_prod,
     |    op_loss(k,i,lat),o2p_loss,np_loss,n2p_loss,
     |    nop_loss,op2p_loss,op2d_loss)
        o2p(k,i,lat,itc) = ionden_local(step,
     |    o2p(k,i,lat,itp),o2p_prod,o2p_loss)
        np(k,i,lat,itc) = ionden_local(step,
     |    np(k,i,lat,itp),np_prod,np_loss)
        n2p(k,i,lat,itc) = ionden_local(step,
     |    n2p(k,i,lat,itp),n2p_prod,n2p_loss)
        nop(k,i,lat,itc) = ionden_local(step,
     |    nop(k,i,lat,itp),nop_prod,nop_loss)
        op2p(k,i,lat,itc) = ionden_local(step,
     |    op2p(k,i,lat,itp),op2p_prod,op2p_loss)
        op2d(k,i,lat,itc) = ionden_local(step,
     |    op2d(k,i,lat,itp),op2d_prod,op2d_loss)
      enddo
!
      call oplus(
     |  tn    (levd0,lond0,latd0,itp),
     |  te    (levd0,lond0,latd0,itp),
     |  ti    (levd0,lond0,latd0,itp),
     |  o2    (levd0,lond0,latd0,itp),
     |  o1    (levd0,lond0,latd0,itp),
     |  he    (levd0,lond0,latd0,itp),
     |  n2    (levd0,lond0,latd0),
     |  n2d   (levd0,lond0,latd0,itp),
     |  un    (levd0,lond0,latd0,itp),
     |  vn    (levd0,lond0,latd0,itp),
     |  w     (levd0,lond0,latd0,itp),
     |  ui    (levd0,lond0,latd0),
     |  vi    (levd0,lond0,latd0),
     |  wi    (levd0,lond0,latd0),
     |  xnmbar(levd0,lond0,latd0,itp),
     |  scht  (levd0,lond0,latd0,itp),
     |  op    (levd0,lond0,latd0,itp),
     |  op_nm (levd0,lond0,latd0,itp),
     |  op    (levd0,lond0,latd0,itc), ! out
     |  op_nm (levd0,lond0,latd0,itc), ! out
     |  op_prod(levd0,lond0,latd0),     ! out
     |  op_loss(levd0,lond0,latd0),     ! out
     |  Fe    (levd0,lond0,latd0),
     |  Fn    (levd0,lond0,latd0),
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after oplus')")
      call timer(time0,time1,'OPLUS',1,0) ! end oplus timing
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in OPLUS =        ',
     |    f12.2,' Dynamics: step ',i5)") time1,istep
      if (diags) then
        do lat=lat0,lat1
          call addfld('OPLUS' ,'O+ Ion',' ',op(:,i0:i1,lat,itc),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('XIOP2P','O+(2P)',' ',op2p(:,i0:i1,lat,itc),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('XIOP2D','O+(2D)',' ',op2d(:,i0:i1,lat,itc),
     |      'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Electron density:
!
      do concurrent (lat=lat0:lat1, i=lon0:lon1)
        do k=1,nlev
          nem(k) = max(op(k,i,lat,itc)+o2p(k,i,lat,itc)+
     |      np(k,i,lat,itc)+n2p(k,i,lat,itc)+nop(k,i,lat,itc),1.0)
        enddo
        do k=2,nlev
          ne(k,i,lat,itc) = sqrt(nem(k-1)*nem(k))
        enddo
        ne(1,i,lat,itc) = sqrt(nem(1)**3/nem(2))
        ne(nlevp1,i,lat,itc) = sqrt(nem(nlev)**3/nem(nlev-1))
      enddo
      if (debug) write(6,"('dynamics after oplus')")
!
! Save diagnostic tec (total electron content):
      call mkdiag_TEC('TEC',tec,z(:,i0:i1,lat0:lat1,itc),
     |  ne(:,i0:i1,lat0:lat1,itc),k0,k1,i0,i1,lat0,lat1)
!
! Save diagnostics hmf2, nmf2:
      call mkdiag_HNMF2('HMF2',z(:,i0:i1,lat0:lat1,itc),
     |  ne(:,i0:i1,lat0:lat1,itc),k0,k1,i0,i1,lat0,lat1)
      call mkdiag_HNMF2('NMF2',z(:,i0:i1,lat0:lat1,itc),
     |  ne(:,i0:i1,lat0:lat1,itc),k0,k1,i0,i1,lat0,lat1)
      call mkdiag_HNMF2('FOF2',z(:,i0:i1,lat0:lat1,itc),
     |  ne(:,i0:i1,lat0:lat1,itc),k0,k1,i0,i1,lat0,lat1)
!
! Ions and electron density are at interface levels (ilev):
      do lat=lat0,lat1
        if (diags) then
          call addfld('NP_diag' ,'N+ Ion' ,'cm^3',np(:,i0:i1,lat,itc),
     |      'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('N2P_diag','N2+ Ion','cm^3',n2p(:,i0:i1,lat,itc),
     |      'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('NOP_diag','NO+ Ion','cm^3',nop(:,i0:i1,lat,itc),
     |      'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('O2P_diag','O2+ Ion','cm^3',o2p(:,i0:i1,lat,itc),
     |      'ilev',k0,k1,'lon',i0,i1,lat)
          call addfld('ELDEN' ,'Electron Density','cm-3',
     |      ne (:,i0:i1,lat,itc),'ilev',k0,k1,'lon',i0,i1,lat)
        endif
      enddo ! lat=lat0,lat1
!
! Lamda ion drag coefficients:
!
      time1_cmpminor = 0.
      do lat=lat0,lat1
        call lamdas(
     |    tn  (levd0,lond0,lat,itp),
     |    xnmbar(levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    he  (levd0,lond0,lat,itp),
     |    n2  (levd0,lond0,lat),
     |    ti  (levd0,lond0,lat,itp),
     |    te  (levd0,lond0,lat,itp),
     |    o2p (levd0,lond0,lat,itp),
     |    op  (levd0,lond0,lat,itp),
     |    np  (levd0,lond0,lat,itp),
     |    n2p (levd0,lond0,lat,itp),
     |    nop (levd0,lond0,lat,itp),
     |    Etot(levd0,lond0,lat),
     |    lxx (levd0,lond0,lat), ! LXX output
     |    lyy (levd0,lond0,lat), ! LYY output
     |    lxy (levd0,lond0,lat), ! LXY output
     |    lyx (levd0,lond0,lat), ! LYX output
     |    lam1(levd0,lond0,lat), ! LAM1 output
     |    ped (levd0,lond0,lat), ! pedersen conductivity output
     |    hall(levd0,lond0,lat), ! hall conductivity output
     |    Q1  (levd0,lond0,lat), ! electron ohmic heating output
     |    Q2  (levd0,lond0,lat), ! electrojet turbulent heating output
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after lamdas: lat=',i3)") lat
        if (diags) then
          call addfld('LXX','Lamda XX',' ',lxx(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('LYY','Lamda YY',' ',lyy(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('LXY','Lamda XY',' ',lxy(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('LYX','Lamda YX',' ',lyx(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('LAM1','LAM1','1/s',lam1(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('SIGMAPED','Pedersen Conductivity','S/m',
     |      ped(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('SIGMAHAL','Hall Conductivity','S/m',
     |      hall(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
        endif
      enddo ! lat=lat0,lat1
!
! Advance n2d:
!
      call timer(time0,time1,'CMPMINOR',0,0) ! start comp minor timing
      do lat=lat0,lat1
        call comp_n2d(
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    no  (levd0,lond0,lat,itp),
     |    ne  (levd0,lond0,lat,itp),
     |    op  (levd0,lond0,lat,itp),
     |    n2p (levd0,lond0,lat,itp),
     |    nop (levd0,lond0,lat,itp),
     |    xnmbar(levd0,lond0,lat,itp),
     |    xnmbari(levd0,lond0,lat,itp),
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after comp_n2d')")
!
! Minor_n2d calls sub minor, which has 3d mpi calls, including gather/scatter
! for fft filtering. Full 3d task subdomains are passed.
!
      call minor_n2d(
     |    tn     (levd0,lond0,latd0,itp),  ! TN (deg K)
     |    o2     (levd0,lond0,latd0,itp),  ! O2 (mmr)
     |    o1     (levd0,lond0,latd0,itp),  ! O (mmr)
     |    he     (levd0,lond0,latd0,itp),  ! HE (mmr)
     |    mbar   (levd0,lond0,latd0,itp),  ! mean molecular weight
     |    barm   (levd0,lond0,latd0,itp),  ! mean molecular weight
     |    xnmbar (levd0,lond0,latd0,itp),  ! p0*e(-z)*mbar/kT
     |    w      (levd0,lond0,latd0,itc),  ! vertical wind
     |    n2d    (levd0,lond0,latd0,itp),  ! n2d from previous step
     |    n2d_nm (levd0,lond0,latd0,itp),  ! n2d at time n-1
     |    n2d    (levd0,lond0,latd0,itc),  ! output n2d
     |    n2d_nm (levd0,lond0,latd0,itc),  ! output n2d at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after minor_n2d')")
      call timer(time0,time1,'CMPMINOR',1,0) ! suspend comp minor timing
      time1_cmpminor = time1_cmpminor+time1
      if (diags) then
        do lat=lat0,lat1
          call addfld('CMP_N2D',' ',' ',n2d(:,i0:i1,lat,itc),
     |      'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Advance n4s:
!
      call timer(time0,time1,'CMPMINOR',0,0) ! resume comp minor timing
      do lat=lat0,lat1
        call comp_n4s(
     |    o2     (levd0,lond0,lat,itp),  ! O2 (mmr)
     |    o1     (levd0,lond0,lat,itp),  ! O (mmr)
     |    n2     (levd0,lond0,lat),      ! N2 (mmr)
     |    xnmbar (levd0,lond0,lat,itp),  ! p0*e(-z)*mbar/kT
     |    xnmbari(levd0,lond0,lat,itp),  ! p0*e(-z)*barm/kT
     |    no     (levd0,lond0,lat,itp),  ! NO
     |    n2d    (levd0,lond0,lat,itc),  ! N(2D) from comp_n2d
     |    ne     (levd0,lond0,lat,itp),  ! Ne
     |    o2p    (levd0,lond0,lat,itp),  ! O2+
     |    op     (levd0,lond0,lat,itp),  ! O+
     |    n2p    (levd0,lond0,lat,itp),  ! N2+ from elden
     |    np     (levd0,lond0,lat,itp),  ! N+  from elden
     |    nop    (levd0,lond0,lat,itp),  ! NO+ from elden
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1 
      if (debug) write(6,"('dynamics after comp_n4s')")
!
! Minor_n4s calls sub minor, which has 3d mpi calls, including gather/scatter
! for fft filtering. Full 3d task subdomains are passed.
!
      call minor_n4s(
     |    tn     (levd0,lond0,latd0,itp),  ! TN (deg K)
     |    o2     (levd0,lond0,latd0,itp),  ! O2 (mmr)
     |    o1     (levd0,lond0,latd0,itp),  ! O (mmr)
     |    he     (levd0,lond0,latd0,itp),  ! HE (mmr)
     |    mbar   (levd0,lond0,latd0,itp),  ! mean molecular weight
     |    barm   (levd0,lond0,latd0,itp),  ! mean molecular weight
     |    xnmbar (levd0,lond0,latd0,itp),  ! p0*e(-z)*mbar/kT
     |    w      (levd0,lond0,latd0,itc),  ! vertical wind
     |    n4s    (levd0,lond0,latd0,itp),  ! n4s from previous step
     |    n4s_nm (levd0,lond0,latd0,itp),  ! n4s at time n-1
     |    n4s    (levd0,lond0,latd0,itc),  ! output n4s
     |    n4s_nm (levd0,lond0,latd0,itc),  ! output n4s at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after minor_n4s')")
      call timer(time0,time1,'CMPMINOR',1,0) ! suspend comp minor timing
      time1_cmpminor = time1_cmpminor+time1
      if (diags) then
        do lat=lat0,lat1
          call addfld('CMP_N4S','N4S density','mmr',
     |      n4s(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
! 
! Advance NO: 
!
      call timer(time0,time1,'CMPMINOR',0,0) ! resume comp minor timing
      do lat=lat0,lat1
        call comp_no(
     |    o2     (levd0,lond0,lat,itp),
     |    o1     (levd0,lond0,lat,itp),
     |    n2     (levd0,lond0,lat),
     |    xnmbar (levd0,lond0,lat,itp),
     |    xnmbari(levd0,lond0,lat,itp),
     |    n4s    (levd0,lond0,lat,itp),
     |    n2d    (levd0,lond0,lat,itc),
     |    o2p    (levd0,lond0,lat,itp),
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1 
      if (debug) write(6,"('dynamics after comp_no')")
!
! Minor_no calls minor, which has 3d mpi calls:
      call minor_no(
     |    tn   (levd0,lond0,latd0,itp),  ! TN (deg K)
     |    o2   (levd0,lond0,latd0,itp),  ! O2 (mmr)
     |    o1   (levd0,lond0,latd0,itp),  ! O (mmr)
     |    he   (levd0,lond0,latd0,itp),  ! HE (mmr)
     |    mbar (levd0,lond0,latd0,itp),  ! mean molecular weight
     |    barm (levd0,lond0,latd0,itp),  ! mean molecular weight
     |    xnmbar(levd0,lond0,latd0,itp), ! p0*e(-z)*mbar/kT
     |    w    (levd0,lond0,latd0,itc),  ! vertical wind
     |    no   (levd0,lond0,latd0,itp),  ! current no
     |    no_nm(levd0,lond0,latd0,itp),  ! no at time n-1
     |    no   (levd0,lond0,latd0,itc),  ! output no
     |    no_nm(levd0,lond0,latd0,itc),  ! output no at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after minor_no')")
      call timer(time0,time1,'CMPMINOR',1,0) ! suspend comp minor timing
      time1_cmpminor = time1_cmpminor+time1
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in CMPMINOR =     ',f12.2,
     |    ' Dynamics: step ',i5)") time1_cmpminor,istep
      if (diags) then
        do lat=lat0,lat1
          call addfld('CMP_NO','NO density','mmr',
     |      no(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Advance argon:
!
      call timer(time0,time1,'CMPMINOR',0,0) ! resume comp minor timing
      do lat=lat0,lat1
        call comp_ar(lon0,lon1,lat)
        if (debug) write(6,"('dynamics after comp_ar: lat=',i3,
     |    ' istep=',i5)") lat,istep
      enddo ! lat=lat0,lat1
!
      call minor_ar(
     |  tn     (levd0,lond0,latd0,itp),   ! tn (deg K)
     |  o2     (levd0,lond0,latd0,itp),   ! o2 (mmr)
     |  o1     (levd0,lond0,latd0,itp),   ! o1 (mmr)
     |  he     (levd0,lond0,latd0,itp),   ! he (mmr)
     |  mbar   (levd0,lond0,latd0,itp),   ! mean molecular weight
     |  barm   (levd0,lond0,latd0,itp),   ! mean molecular weight
     |  xnmbar (levd0,lond0,latd0,itp),   ! p0*e(-z)*mbar/kT
     |  w      (levd0,lond0,latd0,itc),   ! vertical wind
     |  ar     (levd0,lond0,latd0,itp),   ! ar from previous step
     |  ar_nm  (levd0,lond0,latd0,itp),   ! ar at time n-1
     |  ar     (levd0,lond0,latd0,itc),   ! ar output
     |  ar_nm  (levd0,lond0,latd0,itc),   ! ar output at time n-1
     |  1,nlevp1,lon0,lon1,lat0,lat1)
!
      if (debug) write(6,"('dynamics after minor_ar: istep=',i5)") istep
      call timer(time0,time1,'CMPMINOR',1,0) ! suspend comp minor timing
      time1_cmpminor = time1_cmpminor+time1
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in CMPMINOR =     ',f12.2,
     |    ' Dynamics: step ',i5)") time1_cmpminor,istep
      if (diags) then
        do lat=lat0,lat1
          call addfld('CMP_AR',' ',' ',ar(:,i0:i1,lat,itc),
     |      'lev',1,nk,'lon',i0,i1,lat)
        enddo
      endif

!
! Calculate additions to neutral gas heating and O2 dissociation
! due to N, NO chemistry:
!
      do lat=lat0,lat1
        call qjnno(
     |    o2   (levd0,lond0,lat,itp),
     |    o1   (levd0,lond0,lat,itp),
     |    ne   (levd0,lond0,lat,itp),
     |    no   (levd0,lond0,lat,itp),
     |    n4s  (levd0,lond0,lat,itp),
     |    n2d  (levd0,lond0,lat,itc), ! from comp_n2d
     |    xnmbar(levd0,lond0,lat,itp),
     |    1,nlevp1,lon0,lon1,lat)
      enddo
      if (debug) write(6,"('dynamics after qjnno')")

!
! Calculate ion chemistry contribution to neutral gas heating
!   and O2 dissociation:
!
      do lat=lat0,lat1
        call qjion(
     |    o2    (levd0,lond0,lat,itp),
     |    o1    (levd0,lond0,lat,itp),
     |    n2    (levd0,lond0,lat),
     |    o2p   (levd0,lond0,lat,itp),
     |    op    (levd0,lond0,lat,itp),
     |    n4s   (levd0,lond0,lat,itp),
     |    n2d   (levd0,lond0,lat,itp),
     |    no    (levd0,lond0,lat,itp),
     |    ne    (levd0,lond0,lat,itp),
     |    xnmbar(levd0,lond0,lat,itp),
     |    xnmbari(levd0,lond0,lat,itp),
     |    n2p   (levd0,lond0,lat,itp), ! from elden
     |    np    (levd0,lond0,lat,itp), ! from elden
     |    nop   (levd0,lond0,lat,itp), ! from elden
     |    op2p  (levd0,lond0,lat,itp), ! from oplus 
     |    op2d  (levd0,lond0,lat,itp), ! from oplus
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after qjion')")
!
! Electron and Ion temperatures:
!
      do lat=lat0,lat1
!
! Joule heating for ion temperature:
        call qjoule_ti(
     |    un (levd0,lond0,lat,itp), ! zonal neutral velocity input
     |    vn (levd0,lond0,lat,itp), ! meridional velocity input
     |    w  (levd0,lond0,lat,itp), ! vertical velocity input
     |    ui (levd0,lond0,lat),     ! zonal ion velocity input (ExB)
     |    vi (levd0,lond0,lat),     ! meridional ion velocity input (ExB)
     |    wi (levd0,lond0,lat),     ! vertical ion velocity input (ExB)
     |    lam1(levd0,lond0,lat),    ! lambda_1 ion drag coefficient [1/s]
     |    op (levd0,lond0,lat,itp), ! O+ 1/cm^3
     |    np (levd0,lond0,lat,itp), ! N+  1/cm^3
     |    n2p(levd0,lond0,lat,itp), ! N2+ 1/cm^3
     |    nop(levd0,lond0,lat,itp), ! NO+  1/cm^3
     |    o2p(levd0,lond0,lat,itp), ! O2+  1/cm^3
     |    mbar(levd0,lond0,lat,itp),! mean molecular weight
     |    scht(levd0,lond0,lat,itp),! scale height
     |    qji_ti(levd0,lond0,lat),  ! ion joule heating output
     |    1,nlevp1,lon0,lon1)
        if (debug) write(6,"('dynamics after qjoule_ti: lat=',i3)")lat
        if (diags) then
          call addfld('QJI_TI','Joule Heating for TI',' ',
     |      qji_ti(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
        endif
!
        call settei(
     |    tn  (levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    he  (levd0,lond0,lat,itp),
     |    n2  (levd0,lond0,lat),
     |    ne  (levd0,lond0,lat,itp),
     |    te  (levd0,lond0,lat,itp),
     |    ti  (levd0,lond0,lat,itp),
     |    op  (levd0,lond0,lat,itp),
     |    o2p (levd0,lond0,lat,itp),
     |    np  (levd0,lond0,lat,itp),
     |    n2p (levd0,lond0,lat,itp),
     |    nop (levd0,lond0,lat,itp),
     |    xnmbar(levd0,lond0,lat,itp),
     |    xnmbari(levd0,lond0,lat,itp),
     |    scht(levd0,lond0,lat,itp),
     |    schti(levd0,lond0,lat,itp),
     |    qji_ti(levd0,lond0,lat),
     |    Q1  (levd0,lond0,lat),
     |    Q2  (levd0,lond0,lat),
     |    te  (levd0,lond0,lat,itc), ! output
     |    ti  (levd0,lond0,lat,itc), ! output
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after settei: lat=',i3)") lat
        if (diags) then
          call addfld('T_ELEC','Electron Temperature','K',
     |      te(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('T_ION' ,'Ion Temperature','K',
     |      ti(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
        endif
      enddo ! lat=lat0,lat1
!
! Implicit and explicit cooling terms:
      do lat=lat0,lat1
        call newton(
     |    tn  (levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    n2  (levd0,lond0,lat),
     |    no  (levd0,lond0,lat,itp), 
     |    barm(levd0,lond0,lat,itp), 
     |    xnmbar(levd0,lond0,lat,itp),
     |    cp  (levd0,lond0,lat),
     |    cool_implicit(levd0,lond0,lat), ! output implicit cooling
     |    cool_explicit(levd0,lond0,lat), ! output explicit cooling
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after newton: lat=',i3)") lat
!
! Add O(3P) contribution to cooling terms:
        call newto3p(
     |    tn  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    cp  (levd0,lond0,lat),
     |    cool_implicit(levd0,lond0,lat), ! output implicit cooling
     |    cool_explicit(levd0,lond0,lat), ! output explicit cooling
     |    1,nlevp1,lon0,lon1)
        if (debug) write(6,"('dynamics after newto3p: lat=',i3)") lat
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after newton,newto3p')")
      if (diags) then
        do lat=lat0,lat1
          call addfld('COOL_IMP',' ',' ',cool_implicit(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('COOL_EXP',' ',' ',cool_explicit(:,i0:i1,lat),
     |      'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Calculate horizontal diffusion terms for t,u,v,o2,o1,he, using
! coefficients from hdif2 (hdif2 was called from advance).
!
! Make boundary latitudes/longitudes available for inputs to hdif3:
      call hdif_bndlatlons(kldt,kldu,kldv,kldo2,kldo1,kldhe,
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after hdif_bndlatlons')")
!
! Make hdt,u,v,o2,o1,he:
      do lat=lat0,lat1
        call hdif3(
     |    cp(levd0,lond0,lat),     ! specific heat input
     |    fnrh (levd0,lond0,latd0),! 3d input
     |    kldt (levd0,lond0,latd0),! 3d input from hdif2 for tn
     |    kldu (levd0,lond0,latd0),! 3d input from hdif2 for un
     |    kldv (levd0,lond0,latd0),! 3d input from hdif2 for vn
     |    kldo2(levd0,lond0,latd0),! 3d input from hdif2 for o2
     |    kldo1(levd0,lond0,latd0),! 3d input from hdif2 for o1
     |    kldhe(levd0,lond0,latd0),! 3d input from hdif2 for he
     |    hdt (levd0,lond0,lat),   ! 2d tn output at current lat
     |    hdu (levd0,lond0,lat),   ! 2d un output at current lat
     |    hdv (levd0,lond0,lat),   ! 2d vn output at current lat
     |    hdo2(levd0,lond0,lat),   ! 2d o2 output at current lat
     |    hdo1(levd0,lond0,lat),   ! 2d o1 output at current lat
     |    hdhe(levd0,lond0,lat),   ! 2d he output at current lat
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after hdif3')")
!
! Calculate joule heating for tn (3d field qji_tn):
      do lat=lat0,lat1
        call qjoule_tn(
     |    un (levd0,lond0,lat,itp), ! zonal neutral velocity input
     |    vn (levd0,lond0,lat,itp), ! meridional velocity input
     |    w  (levd0,lond0,lat,itp), ! vertical velocity input
     |    ui (levd0,lond0,lat),     ! zonal ion velocity input
     |    vi (levd0,lond0,lat),     ! meridional ion velocity input
     |    wi (levd0,lond0,lat),     ! vertical ion velocity input
     |    lam1(levd0,lond0,lat),    ! lambda_1 ion drag coefficient
     |    scht(levd0,lond0,lat,itp),! scale height
     |    qji_tn(levd0,lond0,lat),  ! ion joule heating output
     |    1,nlevp1,lon0,lon1,lat)
        if (diags) then
          call addfld('QJI_TN','Joule Heating for TN',' ',
     |      qji_tn(:,i0:i1,lat),'lev',k0,k1,'lon',i0,i1,lat)
        endif
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after qjoule_tn')")
!
! Save diagnostic Height-integrated Joule Heating:
      call mkdiag_QJOULE_INTEG('QJOULE_INTEG',
     |  qji_tn(:,i0:i1,lat0:lat1),k0,k1,i0,i1,lat0,lat1)
!
! Calculate 2d (lon,lat) diagnostics for amie run:
!
      if (iamie > 0) then
        do lat=lat0,lat1
          call dyndiag(
     |      ped (levd0,lond0,lat),     ! pedersen conductivity
     |      hall(levd0,lond0,lat),     ! hall conductivity
     |      qji_tn(levd0,lond0,lat),   ! ion Joule heating for tn
     |      z   (levd0,lond0,lat,itc), ! updated Z from addiag
     |      un  (levd0,lond0,lat,itp), ! neutral zonal velocity
     |      vn  (levd0,lond0,lat,itp), ! neutral meridional velocity
     |      ui  (levd0,lond0,lat),     ! ion zonal velocity
     |      vi  (levd0,lond0,lat),     ! ion meridional velocity
     |      1,nlevp1,lon0,lon1,lat)
        enddo ! lat=lat0,lat1
      endif
!
! Advance neutral temperature:
!
      call timer(time0,time1,'DT',0,0) ! start dt timing
      call dt(
     |  tn	(levd0,lond0,latd0,itp),
     |  tn_nm   (levd0,lond0,latd0,itp),
     |  un	(levd0,lond0,latd0,itp),
     |  vn	(levd0,lond0,latd0,itp),
     |  o1	(levd0,lond0,latd0,itp),
     |  mbar  (levd0,lond0,latd0,itp),
     |  xnmbar(levd0,lond0,latd0,itp),
     |  scht  (levd0,lond0,latd0,itp),
     |  schti (levd0,lond0,latd0,itp),
     |  cp	(levd0,lond0,latd0),	  ! specific heat	 (cpktkm.F)
     |  kt	(levd0,lond0,latd0),	  ! molecular diffusion  (cpktkm.F)
     |  km	(levd0,lond0,latd0),	  ! molecular viscosity  (cpktkm.F)
     |  hdt	(levd0,lond0,latd0),	  ! horizontal diffusion (hdif.F)
     |  qji_tn  (levd0,lond0,latd0),      ! joule heating (sub qjoule_tn)
     |  cool_implicit(levd0,lond0,latd0), ! implicit cooling (newton.F)
     |  cool_explicit(levd0,lond0,latd0), ! explicit cooling (newton.F)
     |  w	(levd0,lond0,latd0,itc),  ! updated W (swdot.F)
     |  tn	(levd0,lond0,latd0,itc),  ! output updated tn
     |  tn_nm (levd0,lond0,latd0,itc),    ! output updated tn at time n-1
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after dt')")
      call timer(time0,time1,'DT',1,0) ! end dt timing
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in DT =           ',f12.2,
     |    ' Dynamics: step ',i5)") time1,istep
      if (diags) then
        do lat=lat0,lat1
          call addfld('TN_UPD','Updated TN','K',
     |      tn(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Advance neutral velocities:
!
      call timer(time0,time1,'DUV',0,0) ! start duv timing
      call duv(
     |  tn    (levd0,lond0,latd0,itp), ! neutral temperature
     |  tn    (levd0,lond0,latd0,itc), ! updated neutral temperature (dt.F)
     |  tn_nm (levd0,lond0,latd0,itp), ! tn at time n-1
     |  un    (levd0,lond0,latd0,itp), ! zonal velocity
     |  vn    (levd0,lond0,latd0,itp), ! meridional velocity
     |  un_nm (levd0,lond0,latd0,itp), ! zonal velocity at time n-1
     |  vn_nm (levd0,lond0,latd0,itp), ! meridional velocity at time n-1
     |  w     (levd0,lond0,latd0,itc), ! updated vertical velocity (swdot.F)
     |  mbar  (levd0,lond0,latd0,itp), ! mean molecular weight
     |  scht  (levd0,lond0,latd0,itp), ! scale height
     |  schti (levd0,lond0,latd0,itp), ! scale height
     |  z     (levd0,lond0,latd0,itp), ! geopotential height
     |  hdu   (levd0,lond0,latd0),     ! horizontal diffusion of U (hdif3.F)
     |  hdv   (levd0,lond0,latd0),     ! horizontal diffusion of V (hdif3.F)
     |  ui    (levd0,lond0,latd0),     ! zonal ion velocity input
     |  vi    (levd0,lond0,latd0),     ! meridional ion velocity input
     |  lxx   (levd0,lond0,latd0),     ! xx ion drag coefficient
     |  lyy   (levd0,lond0,latd0),     ! yy ion drag coefficient
     |  lxy   (levd0,lond0,latd0),     ! xy ion drag coefficient
     |  lyx   (levd0,lond0,latd0),     ! yx ion drag coefficient
     |  km    (levd0,lond0,latd0),     ! molecular viscosity  (cpktkm.F)
     |  Fe    (levd0,lond0,latd0),     ! eastward field-aligned ion drag
     |  Fn    (levd0,lond0,latd0),     ! northward field-aligned ion drag
     |  un    (levd0,lond0,latd0,itc), ! output updated un
     |  un_nm (levd0,lond0,latd0,itc), ! output updated un at time n-1
     |  vn    (levd0,lond0,latd0,itc), ! output updated vn
     |  vn_nm (levd0,lond0,latd0,itc), ! output updated vn at time n-1
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after duv')")
      call timer(time0,time1,'DUV',1,0) ! end duv timing
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in DUV =          ',f12.2,
     |    ' Dynamics: step ',i5)") time1,istep
      if (diags) then
        do lat=lat0,lat1
          call addfld('UN_UPD','Updated UN','cm/s',
     |      un(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('VN_UPD','Updated VN','cm/s',
     |      vn(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif
!
! Sources and sinks for major species composition O2, O:
!
      do lat=lat0,lat1
        call comp_o2o(
     |    o1   (levd0,lond0,lat,itp), ! atomic oxygen (mmr)
     |    mbar (levd0,lond0,lat,itp), ! mean molecular weight
     |    xnmbar(levd0,lond0,lat,itp),! p0*e(-z)*mbar/kT
     |    op   (levd0,lond0,lat,itp), ! O+
     |    no   (levd0,lond0,lat,itp), ! NO
     |    n4s  (levd0,lond0,lat,itp), ! N(4S)
     |    n2d  (levd0,lond0,lat,itp), ! N(2D)
     |    o2p  (levd0,lond0,lat,itp), ! O2+
     |    ne   (levd0,lond0,lat,itp), ! electron density
     |    n2p  (levd0,lond0,lat,itp), ! N2+
     |    np   (levd0,lond0,lat,itp), ! N+
     |    nop  (levd0,lond0,lat,itp), ! NO+
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after comp_o2o')")
!
! Advance O2, O, HE:
!
      time1_cmpmajor = 0.
      call timer(time0,time1,'CMPMAJOR',0,0) ! start comp major timing
      call comp(
     |  tn    (levd0,lond0,latd0,itp), ! neutral temperature
     |  o2    (levd0,lond0,latd0,itp), ! O2 (mmr)
     |  o2_nm (levd0,lond0,latd0,itp), ! O2 (mmr) at time n-1
     |  o1    (levd0,lond0,latd0,itp), ! O1 (mmr)
     |  o1_nm (levd0,lond0,latd0,itp), ! O1 (mmr) at time n-1
     |  he    (levd0,lond0,latd0,itp), ! He (mmr)
     |  he_nm (levd0,lond0,latd0,itp), ! He (mmr) at time n-1
     |  w     (levd0,lond0,latd0,itp), ! vertical velocity
     |  mbar  (levd0,lond0,latd0,itp), ! mean molecular weight
     |  xnmbar(levd0,lond0,latd0,itp),
     |  hdo2  (levd0,lond0,latd0),     ! O2 horizontal diffusion (hdif3)
     |  hdo1  (levd0,lond0,latd0),     ! O  horizontal diffusion (hdif3)
     |  hdhe  (levd0,lond0,latd0),     ! He horizontal diffusion (hdif3)
     |  o2    (levd0,lond0,latd0,itc), ! output: O2 updated for current step 
     |  o2_nm (levd0,lond0,latd0,itc), ! output: O2 updated for previous step
     |  o1    (levd0,lond0,latd0,itc), ! output: O1 updated for current step
     |  o1_nm (levd0,lond0,latd0,itc), ! output: O1 updated for previous step
     |  he    (levd0,lond0,latd0,itc), ! output: He updated for current step
     |  he_nm (levd0,lond0,latd0,itc), ! output: He updated for previous step
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after comp')")
      call timer(time0,time1,'CMPMAJOR',1,0)
      if (timing%level >= 2.and.time2print(nstep,istep))
     |  write(6,"('Time in CMPMAJOR =     ',f12.2,
     |    ' Dynamics: step ',i5)") time1_cmpmajor,istep
      if (diags) then
        do lat=lat0,lat1
          call addfld('O2_UPD','Updated O2','mmr',
     |      o2(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
          call addfld('O1_UPD','Updated O1','mmr',
     |      o1(:,i0:i1,lat,itc),'lev',k0,k1,'lon',i0,i1,lat)
        enddo
      endif

      if (debug) write(6,"('dynamics returning')")
      end subroutine dynamics

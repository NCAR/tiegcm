!
      subroutine dynamics(istep,do_rtc)
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
      use fields_module
      use mpi_module,only: lon0,lon1,lat0,lat1,mytidi
#ifdef MPI
      use mpi_module,only: distribute_1d,mp_close
#endif
      use bndry_module,only: bndcmp
      use chemrates_module,only: chemrates_tdep,disn2p
      use qrj_module,only: qrj,qtotal,qop,qo2p,qn2p,qnp,qnop
      use chapman_module,only: chapman
      use oplus_module,only: oplus
      use input_module,only: iaurora=>aurora
      use aurora_module,only: aurora
      use n4s_module,only: comp_n4s,minor_n4s
      use no_module,only: comp_no,minor_no
      use hdif_module,only: hdif3,hdif_bndlons,hdif_periodic
      use init_module,only: igetgswmdi,igetgswmsdi, ! GSWM input lower boundary
     |   igetgswmnmidi,igetgswmnmisdi
      implicit none
!
! Args:
      integer,intent(in) :: istep
      logical,intent(in) :: do_rtc
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,n,lat,ier
      integer :: i0,i1,nk,nkm1,nlats
      logical,parameter :: debug=.false. ! add prints to stdout
      logical,parameter :: diags=.false. ! put fields on sec hist
      real :: time0,tsec,rtc0_qrj,rtc_qrj
!
      if (debug) write(6,"(/,'Enter dynamics.')")
      i0=lon0 ; i1=lon1 ; nk=nlevp1 ; nkm1=nk-1 ! for addfsech
      nlats = lat1-lat0+1
!
      call bndcmp
      if (debug) write(6,"('dynamics after bndcmp')")
!
! Calculate specific heat and molecular viscosity.
! f4d itp fields are input, f3d fields are output.
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call cpktkm(
     |    tn(levd0,lond0,lat,itp), ! 4d input
     |    o2(levd0,lond0,lat,itp), ! 4d input
     |    o1(levd0,lond0,lat,itp), ! 4d input
     |    cp(levd0,lond0,lat),     ! 3d output
     |    kt(levd0,lond0,lat),     ! 3d output
     |    km(levd0,lond0,lat),     ! 3d output
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after cpktkm')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for cpktkm=',f8.2)")
     |    istep,tsec 
      endif
      if (diags) then
        do lat=lat0,lat1
          call addfsech('CP' ,' ',' ',cp(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('KT' ,' ',' ',kt(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('KM' ,' ',' ',km(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
        enddo
      endif
!
! Calculate omega for vertical velocity:
! (latitude scan is inside swdot)
!
      if (do_rtc) call timer(time0,tsec,'begin')
      if (debug) write(6,"('dynamics before swdot')")
      call swdot(
     |  un(levd0,lond0,latd0,itp), ! un input
     |  vc(levd0,lond0,latd0,itp), ! vc input
     |  w (levd0,lond0,latd0,itc), ! omega output
     |  1,nlevp1,lon0,lon1,lat0,lat1,lat0,lat1)
      if (debug) write(6,"('dynamics after swdot')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for swdot=',f8.2)")
     |    istep,tsec 
      endif
      if (diags) then
        do lat=lat0,lat1
          call addfsech('OMEGA' ,' ',' ',w(:,i0:i1,lat,itc),
     |      i0,i1,nk,nkm1,lat)
        enddo
      endif
!
! Calculate ExB ion drift velocities from electric field:
!
      if (do_rtc) call timer(time0,tsec,'begin')
      call ionvel(
     |  z (levd0,lond0,latd0,itc), ! geopotential input
     |  ui(levd0,lond0,latd0),     ! U ion velocity output
     |  vi(levd0,lond0,latd0),     ! V ion velocity output
     |  wi(levd0,lond0,latd0),     ! W ion velocity output
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after ionvel')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for ionvel=',f8.2)")
     |    istep,tsec 
      endif
      if (diags) then
        do lat=lat0,lat1
          call addfsech('UIVEL' ,' ',' ',ui(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('VIVEL' ,' ',' ',vi(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('WIVEL' ,' ',' ',wi(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
        enddo
      endif
!
! Calculate column densities and line integrals of o2, o, and n2.
! Prognostics z,tn,o2,o1,barm are input, diagnostics vo2,vo1,vn2,
!   sco2,sco1,scn2 are output.
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call chapman(
     |    z   (levd0,lond0,lat,itc), ! updated Z from addiag
     |    tn  (levd0,lond0,lat,itp), ! 4d in
     |    o2  (levd0,lond0,lat,itp), ! 4d in
     |    o1  (levd0,lond0,lat,itp), ! 4d in
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
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for chapman=',f8.2)")
     |    istep,tsec 
      endif
      if (diags) then
        do lat=lat0,lat1
          call addfsech('SCO2' ,' ',' ',sco2(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('SCO1' ,' ',' ',sco1(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('SCN2' ,' ',' ',scn2(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
        enddo
      endif
!
! Calculate temperature dependent reaction rates (chemrates_module): 
!
      if (do_rtc) call timer(time0,tsec,'begin')
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
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for chemrates=',
     |    f8.2)") istep,tsec 
      endif
!
! Latitude scan for qrj, qinite, xray and aurora.
!
      if (do_rtc) then
        call timer(time0,tsec,'begin')
        rtc_qrj = 0.
      endif
!
#ifdef VT
!     code = 122 ; state = 'ions' ; activity='ModelCode'
      call vtbegin(122,ier)
#endif
      do lat=lat0,lat1
!
! Calculate ionization and heating rates:
! Note timer has barrier calls, and timing for qrj is
!   inside j-loop.
!
        if (do_rtc) call timer(rtc0_qrj,tsec,'begin')
!       if (debug) write(6,"('dynamics call qrj: lat=',i3)") lat
        call qrj(
     |    sco2(levd0,lond0,lat),
     |    sco1(levd0,lond0,lat),
     |    scn2(levd0,lond0,lat),
     |    tn  (levd0,lond0,lat,itp),
     |    no  (levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    n4s (levd0,lond0,lat,itp),
     |    xnmbari(levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after qrj: lat=',i3)") lat
        if (do_rtc) then
          call timer(rtc0_qrj,tsec,'end')
          rtc_qrj = rtc_qrj+tsec
        endif
        if (diags) then
          call addfsech('QTOTAL' ,' ',' ',qtotal(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('QOP'    ,' ',' ',qop(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('QO2P'   ,' ',' ',qo2p(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('QN2P'   ,' ',' ',qn2p(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('QNP'    ,' ',' ',QNP(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('QNOP'   ,' ',' ',qnop(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
        endif
!
! Calculate background (nightime) ionization:
        call qinite(
     |    o2     (levd0,lond0,lat,itp),
     |    o1     (levd0,lond0,lat,itp),
     |    no     (levd0,lond0,lat,itp),
     |    xnmbari(levd0,lond0,lat),
     |    vo2    (levd0,lond0,lat),
     |    vo1    (levd0,lond0,lat),
     |    vn2    (levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after qinite: lat=',i3)") lat
!
! Ionization from xrays:
        call xray(
     |    tn  (levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    barm(levd0,lond0,lat,itp),
     |    sco2(levd0,lond0,lat),
     |    sco1(levd0,lond0,lat),
     |    scn2(levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after xray: lat=',i3)") lat
!
! Add aurora particle precip and ionization:
        if (iaurora > 0) then
          call aurora(
     |      tn  (levd0,lond0,lat,itp),
     |      o2  (levd0,lond0,lat,itp),
     |      o1  (levd0,lond0,lat,itp),
     |      barm(levd0,lond0,lat,itp),
     |      1,nlevp1,lon0,lon1,lat)
          if (debug) write(6,"('dynamics after aurora: lat=',i3)") lat
        endif
      enddo ! lat=lat0,lat1 for qrj, qinite, xray, and aurora
      if (debug) write(6,"('dynamics after qrj,qinite,xray,aurora')")
!
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for qrj+qinite+',
     |    'xray+aurora=',f8.2)") istep,tsec 
        write(6,"('Dynamics step ',i4,': rtc time for qrj only=',
     |    f8.2)") istep,rtc_qrj
      endif
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
      if (do_rtc) call timer(time0,tsec,'begin')
        call oplus(
     |    tn    (levd0,lond0,latd0,itp),
     |    te    (levd0,lond0,latd0,itp),
     |    ti    (levd0,lond0,latd0,itp),
     |    o2    (levd0,lond0,latd0,itp),
     |    o1    (levd0,lond0,latd0,itp),
     |    n2d   (levd0,lond0,latd0,itp),
     |    ne    (levd0,lond0,latd0,itp),
     |    un    (levd0,lond0,latd0,itp),
     |    vn    (levd0,lond0,latd0,itp),
     |    w     (levd0,lond0,latd0,itp),
     |    barm  (levd0,lond0,latd0,itp),
     |    ui    ,
     |    vi    ,
     |    wi    ,
     |    xnmbar,
     |    op    (levd0,lond0,latd0,itp),
     |    op_nm (levd0,lond0,latd0,itp),
     |    op    (levd0,lond0,latd0,itc), ! out
     |    op_nm (levd0,lond0,latd0,itc), ! out
     |    xiop2p,                        ! out
     |    xiop2d,                        ! out
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after oplus')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for oplus=',
     |    f8.2)") istep,tsec 
      endif
      if (diags) then
        do lat=lat0,lat1
          call addfsech('OPLUS' ,' ',' ',op(:,i0:i1,lat,itc),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('XIOP2P',' ',' ',xiop2p(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('XIOP2D',' ',' ',xiop2d(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
        enddo
      endif
!
! Electron density:
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call elden(
     |    tn    (levd0,lond0,lat,itp),
     |    barm  (levd0,lond0,lat,itp),
     |    op    (levd0,lond0,lat,itp), ! O+ at current time step
     |    op    (levd0,lond0,lat,itc), ! O+ updated by sub oplus
     |    o2    (levd0,lond0,lat,itp),
     |    o1    (levd0,lond0,lat,itp),
     |    n2d   (levd0,lond0,lat,itp),
     |    no    (levd0,lond0,lat,itp),
     |    n4s   (levd0,lond0,lat,itp),
     |    xiop2p(levd0,lond0,lat),      ! from oplus
     |    xiop2d(levd0,lond0,lat),      ! from oplus
     |    nplus (levd0,lond0,lat),      ! N+  output
     |    n2p   (levd0,lond0,lat),      ! N2+ output
     |    nop   (levd0,lond0,lat),      ! NO+ output
     |    o2p   (levd0,lond0,lat,itc),  ! O2+ output
     |    ne    (levd0,lond0,lat,itc),  ! electron density output
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after oplus')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for elden=',
     |    f8.2)") istep,tsec 
      endif
      do lat=lat0,lat1
        if (diags) then
          call addfsech('NPLUS'  ,' ',' ',nplus(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('ELD_N2P',' ',' ',n2p(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('ELD_NOP',' ',' ',nop(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('ELD_O2P',' ',' ',o2p(:,i0:i1,lat,itc),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('ELD_NE' ,' ',' ',ne (:,i0:i1,lat,itc),
     |      i0,i1,nk,nkm1,lat)
        endif
      enddo ! lat=lat0,lat1

!
! Lamda ion drag coefficients:
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call lamdas(
     |    tn  (levd0,lond0,lat,itp),
     |    barm(levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    ti  (levd0,lond0,lat,itp),
     |    te  (levd0,lond0,lat,itp),
     |    o2p (levd0,lond0,lat,itp),
     |    op  (levd0,lond0,lat,itp),
     |    nop (levd0,lond0,lat),
     |    lxx (levd0,lond0,lat), ! LXX output
     |    lyy (levd0,lond0,lat), ! LYY output
     |    lxy (levd0,lond0,lat), ! LXY output
     |    lyx (levd0,lond0,lat), ! LYX output
     |    ped (levd0,lond0,lat), ! pedersen conductivity output
     |    hall(levd0,lond0,lat), ! hall conductivity output
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after lamdas: lat=',i3)") lat
        if (diags) then
          call addfsech('LXX',' ',' ',lxx(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('LYY',' ',' ',lyy(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('LXY',' ',' ',lxy(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('LYX',' ',' ',lyx(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('SIGMAPED',' ',' ',ped(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('SIGMAHAL',' ',' ',hall(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
        endif
!
! Add qn2p to qtef:
        call qtieff(1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after qtieff: lat=',i3)") lat
!
! Calculate n2d:
        call comp_n2d(
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    no  (levd0,lond0,lat,itp),
     |    n4s (levd0,lond0,lat,itp),
     |    ne  (levd0,lond0,lat,itp),
     |    op  (levd0,lond0,lat,itp),
     |    n2p (levd0,lond0,lat),
     |    nop (levd0,lond0,lat),
     |    sco2(levd0,lond0,lat),
     |    xnmbarm(levd0,lond0,lat),
     |    n2d (levd0,lond0,lat,itc), ! N2D output
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after comp_n2d: lat=',i3)") lat
      enddo ! lat=lat0,lat1
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for lamdas+qtieff+',
     |    'comp_n2d=',f8.2)") istep,tsec 
      endif
      if (diags) then
        call addfsech('CMP_N2D',' ',' ',n2d(:,i0:i1,lat,itc),
     |    i0,i1,nk,nkm1,lat)
      endif
!
! Advance n4s:
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call comp_n4s(
     |    tn     (levd0,lond0,lat,itp),  ! TN (deg K)
     |    o2     (levd0,lond0,lat,itp),  ! O2 (mmr)
     |    o1     (levd0,lond0,lat,itp),  ! O (mmr)
     |    barm   (levd0,lond0,lat,itp),  ! mean mol weight
     |    xnmbarm(levd0,lond0,lat),      ! p0*e(-z)*barm
     |    no     (levd0,lond0,lat,itp),  ! NO
     |    n2d    (levd0,lond0,lat,itc),  ! N(2D) from comp_n2d
     |    ne     (levd0,lond0,lat,itp),  ! Ne
     |    o2p    (levd0,lond0,lat,itp),  ! O2+
     |    op     (levd0,lond0,lat,itp),  ! O+
     |    n2p    (levd0,lond0,lat),      ! N2+ from elden
     |    nplus  (levd0,lond0,lat),      ! N+  from elden
     |    nop    (levd0,lond0,lat),      ! NO+ from elden
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1 
      if (debug) write(6,"('dynamics after comp_n4s')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for comp_n4s=',f8.2)")
     |    istep,tsec 
      endif
!
! Minor_n4s calls sub minor, which has 3d mpi calls, including gather/scatter
! for fft filtering. Full 3d task subdomains are passed.
!
      if (do_rtc) call timer(time0,tsec,'begin')
      call minor_n4s(
     |    tn     (levd0,lond0,latd0,itp),  ! TN (deg K)
     |    o2     (levd0,lond0,latd0,itp),  ! O2 (mmr)
     |    o1     (levd0,lond0,latd0,itp),  ! O (mmr)
     |    n4s    (levd0,lond0,latd0,itp),  ! n4s from previous step
     |    n4s_nm (levd0,lond0,latd0,itp),  ! n4s at time n-1
     |    n4s    (levd0,lond0,latd0,itc),  ! output n4s
     |    n4s_nm (levd0,lond0,latd0,itc),  ! output n4s at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after minor_n4s')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for minor_n4s=',
     |    f8.2)") istep,tsec 
      endif
      if (diags) then
        call addfsech('CMP_N2D',' ',' ',n4s(:,i0:i1,lat,itc),
     |    i0,i1,nk,nkm1,lat)
      endif
! 
! Advance NO: 
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call comp_no(
     |    tn     (levd0,lond0,lat,itp),
     |    o2     (levd0,lond0,lat,itp),
     |    o1     (levd0,lond0,lat,itp),
     |    barm   (levd0,lond0,lat,itp),
     |    xnmbarm(levd0,lond0,lat),
     |    n4s    (levd0,lond0,lat,itp),
     |    n2d    (levd0,lond0,lat,itc),
     |    ne     (levd0,lond0,lat,itp),
     |    o2p    (levd0,lond0,lat,itp),
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1 
      if (debug) write(6,"('dynamics after comp_no')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for comp_no=',f8.2)")
     |    istep,tsec 
      endif
!
! Minor_no calls minor, which has 3d mpi calls:
      if (do_rtc) call timer(time0,tsec,'begin')
      call minor_no(
     |    tn   (levd0,lond0,latd0,itp),  ! TN (deg K)
     |    o2   (levd0,lond0,latd0,itp),  ! O2 (mmr)
     |    o1   (levd0,lond0,latd0,itp),  ! O (mmr)
     |    no   (levd0,lond0,latd0,itp),  ! current no
     |    no_nm(levd0,lond0,latd0,itp),  ! no at time n-1
     |    no   (levd0,lond0,latd0,itc),  ! output no
     |    no_nm(levd0,lond0,latd0,itc),  ! output no at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after minor_no')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for minor_no=',f8.2)")
     |    istep,tsec 
      endif
      if (diags) then
        call addfsech('CMP_NO',' ',' ',no(:,i0:i1,lat,itc),
     |    i0,i1,nk,nkm1,lat)
      endif
!
! Calculate additions to neutral gas heating and O2 dissociation
! due to N, NO chemistry:
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call qjnno(
     |    tn   (levd0,lond0,lat,itp),
     |    o2   (levd0,lond0,lat,itp),
     |    o1   (levd0,lond0,lat,itp),
     |    ne   (levd0,lond0,lat,itp),
     |    no   (levd0,lond0,lat,itp),
     |    n4s  (levd0,lond0,lat,itp),
     |    n2d  (levd0,lond0,lat,itc), ! from comp_n2d
     |    xnmbarm(levd0,lond0,lat),
     |    1,nlevp1,lon0,lon1,lat)
      enddo
      if (debug) write(6,"('dynamics after qjnno')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for qjnno=',f8.2)")
     |    istep,tsec 
      endif
!
! Calculate ion chemistry contribution to neutral gas heating
!   and O2 dissociation:
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
        call qjion(
     |    tn    (levd0,lond0,lat,itp),
     |    o2    (levd0,lond0,lat,itp),
     |    o1    (levd0,lond0,lat,itp),
     |    o2p   (levd0,lond0,lat,itp),
     |    op    (levd0,lond0,lat,itp),
     |    n4s   (levd0,lond0,lat,itp),
     |    n2d   (levd0,lond0,lat,itp),
     |    no    (levd0,lond0,lat,itp),
     |    ne    (levd0,lond0,lat,itp),
     |    barm  (levd0,lond0,lat,itp),
     |    n2p   (levd0,lond0,lat),     ! from elden
     |    nplus (levd0,lond0,lat),     ! from elden
     |    nop   (levd0,lond0,lat),     ! from elden
     |    xiop2p(levd0,lond0,lat),     ! from oplus 
     |    xiop2d(levd0,lond0,lat),     ! from oplus
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after qjion')")
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for qjion=',f8.2)")
     |    istep,tsec 
      endif
!
! Electron and Ion temperatures:
!
      if (do_rtc) call timer(time0,tsec,'begin')
      do lat=lat0,lat1
!
! Joule heating for ion temperature:
        call qjoule_ti(
     |    un (levd0,lond0,lat,itp), ! zonal neutral velocity input
     |    vn (levd0,lond0,lat,itp), ! meridional velocity input
     |    ui (levd0,lond0,lat),     ! zonal ion velocity input
     |    vi (levd0,lond0,lat),     ! meridional ion velocity input
     |    lxx(levd0,lond0,lat),     ! xx ion drag coefficient
     |    lyy(levd0,lond0,lat),     ! yy ion drag coefficient
     |    lxy(levd0,lond0,lat),     ! xy ion drag coefficient
     |    lyx(levd0,lond0,lat),     ! yx ion drag coefficient
     |    qji_ti(levd0,lond0,lat),  ! ion joule heating output
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after qjoule_ti: lat=',i3)")lat
!
        call settei(
     |    tn  (levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    ne  (levd0,lond0,lat,itp),
     |    te  (levd0,lond0,lat,itp),
     |    ti  (levd0,lond0,lat,itp),
     |    op  (levd0,lond0,lat,itp),
     |    o2p (levd0,lond0,lat,itp),
     |    nop (levd0,lond0,lat),
     |    barm(levd0,lond0,lat,itp),
     |    qji_ti(levd0,lond0,lat),
     |    te  (levd0,lond0,lat,itc), ! output
     |    ti  (levd0,lond0,lat,itc), ! output
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after settei: lat=',i3)") lat
      enddo ! lat=lat0,lat1
      if (do_rtc) then
        call timer(time0,tsec,'end')
        write(6,"('Dynamics step ',i4,': rtc time for ',
     |    'qjoule_ti+settei=',f8.2)") istep,tsec 
      endif
      if (diags) then
        do lat=lat0,lat1
          call addfsech('TE_UPD',' ',' ',te(:,i0:i1,lat,itc),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('TI_UPD',' ',' ',ti(:,i0:i1,lat,itc),
     |      i0,i1,nk,nkm1,lat)
        enddo
      endif
!
! Implicit and explicit cooling terms:
      do lat=lat0,lat1
        call newton(
     |    tn  (levd0,lond0,lat,itp),
     |    o2  (levd0,lond0,lat,itp),
     |    o1  (levd0,lond0,lat,itp),
     |    no  (levd0,lond0,lat,itp),
     |    barm(levd0,lond0,lat,itp),
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
     |    1,nlevp1,lon0,lon1,lat)
        if (debug) write(6,"('dynamics after newto3p: lat=',i3)") lat
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after newton,newto3p')")
      if (diags) then
        do lat=lat0,lat1
          call addfsech('COOL_IMP',' ',' ',cool_implicit(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
          call addfsech('COOL_EXP',' ',' ',cool_explicit(:,i0:i1,lat),
     |      i0,i1,nk,nkm1,lat)
        enddo
      endif
!
! Calculate horizontal diffusion terms for t,u,v,o2,o1, using
! coefficients from hdif2 (hdif2 was called from advance).
!
! Make boundary longitudes available for inputs to hdif3:
      call hdif_bndlons(kldt,kldu,kldv,kldo2,kldo1,
     |    1,nlevp1,lon0,lon1,lat0,lat1)
!
! Make hdt,u,v,o2,o1:
      do lat=lat0,lat1
        call hdif3(
     |    cp(levd0,lond0,lat),     ! specific heat input
     |    kldt ,                   ! 3d input from hdif2 for tn
     |    kldu ,                   ! 3d input from hdif2 for un
     |    kldv ,                   ! 3d input from hdif2 for vn
     |    kldo2,                   ! 3d input from hdif2 for o2
     |    kldo1,                   ! 3d input from hdif2 for o1
     |    hdt (levd0,lond0,lat),   ! 2d tn output at current lat
     |    hdu (levd0,lond0,lat),   ! 2d un output at current lat
     |    hdv (levd0,lond0,lat),   ! 2d vn output at current lat
     |    hdo2(levd0,lond0,lat),   ! 2d o2 output at current lat
     |    hdo1(levd0,lond0,lat),   ! 2d o1 output at current lat
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after hdif3')")
!
! Periodic points for horizontal diffusion terms output from hdif3:
! This may not be necessary.
!     call hdif_periodic(hdt,hdu,hdv,hdo2,hdo1,
!    |    1,nlevp1,lon0,lon1,lat0,lat1)
!
! Calculate joule heating for tn (3d field qji_tn):
      do lat=lat0,lat1
        call qjoule_tn(
     |    un (levd0,lond0,lat,itp), ! zonal neutral velocity input
     |    vn (levd0,lond0,lat,itp), ! meridional velocity input
     |    ui (levd0,lond0,lat),     ! zonal ion velocity input
     |    vi (levd0,lond0,lat),     ! meridional ion velocity input
     |    lxx(levd0,lond0,lat),     ! xx ion drag coefficient
     |    lyy(levd0,lond0,lat),     ! yy ion drag coefficient
     |    lxy(levd0,lond0,lat),     ! xy ion drag coefficient
     |    lyx(levd0,lond0,lat),     ! yx ion drag coefficient
     |    qji_tn(levd0,lond0,lat),  ! ion joule heating output
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after qjoule_tn')")
!    
      if( igetgswmdi==0.and.igetgswmsdi==0.and.      ! no GSWM input used 
     |  igetgswmnmidi==0.and.igetgswmnmisdi==0) then
!
! Advance neutral temperature:
!
        call dt(
     |    tn	(levd0,lond0,latd0,itp),
     |    tn_nm (levd0,lond0,latd0,itp),
     |    un	(levd0,lond0,latd0,itp),
     |    vn	(levd0,lond0,latd0,itp),
     |    o2	(levd0,lond0,latd0,itp),
     |    o1	(levd0,lond0,latd0,itp),
     |    barm  (levd0,lond0,latd0,itp),    ! mean molecular weight
     |    cp	(levd0,lond0,latd0),	    ! specific heat	   (cpktkm.F)
     |    kt	(levd0,lond0,latd0),	    ! molecular diffusion  (cpktkm.F)
     |    km	(levd0,lond0,latd0),	    ! molecular viscosity  (cpktkm.F)
     |    hdt	(levd0,lond0,latd0),	    ! horizontal diffusion (hdif.F)
     |    qji_tn(levd0,lond0,latd0),	    ! joule heating (sub qjoule_tn)
     |    cool_implicit(levd0,lond0,latd0), ! implicit cooling (newton.F)
     |    cool_explicit(levd0,lond0,latd0), ! explicit cooling (newton.F)
     |    w	(levd0,lond0,latd0,itc),    ! updated W (swdot.F)
     |    tn	(levd0,lond0,latd0,itc),    ! output updated tn
     |    tn_nm (levd0,lond0,latd0,itc),    ! output updated tn at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
        if (debug) write(6,"('dynamics after dt')")
!
! Advance neutral velocities:
!
        call duv(
     |    tn	(levd0,lond0,latd0,itp), ! neutral temperature
     |    tn	(levd0,lond0,latd0,itc), ! updated neutral temperature (dt.F)
     |    tn_nm (levd0,lond0,latd0,itp), ! tn at time n-1
     |    un	(levd0,lond0,latd0,itp), ! zonal velocity
     |    vn	(levd0,lond0,latd0,itp), ! meridional velocity
     |    un_nm (levd0,lond0,latd0,itp), ! zonal velocity at time n-1
     |    vn_nm (levd0,lond0,latd0,itp), ! meridional velocity at time n-1
     |    w	(levd0,lond0,latd0,itc), ! updated vertical velocity (swdot.F)
     |    barm  (levd0,lond0,latd0,itp), ! mean molecular weight
     |    z	(levd0,lond0,latd0,itp), ! geopotential height
     |    hdu	(levd0,lond0,latd0),	 ! horizontal diffusion of U (hdif3.F)
     |    hdv	(levd0,lond0,latd0),	 ! horizontal diffusion of V (hdif3.F)
     |    ui	(levd0,lond0,latd0),	 ! zonal ion velocity input
     |    vi	(levd0,lond0,latd0),	 ! meridional ion velocity input
     |    lxx	(levd0,lond0,latd0),	 ! xx ion drag coefficient
     |    lyy	(levd0,lond0,latd0),	 ! yy ion drag coefficient
     |    lxy	(levd0,lond0,latd0),	 ! xy ion drag coefficient
     |    lyx	(levd0,lond0,latd0),	 ! yx ion drag coefficient
     |    km	(levd0,lond0,latd0),	 ! molecular viscosity  (cpktkm.F)
     |    un	(levd0,lond0,latd0,itc), ! output updated un
     |    un_nm (levd0,lond0,latd0,itc), ! output updated un at time n-1
     |    vn	(levd0,lond0,latd0,itc), ! output updated vn
     |    vn_nm (levd0,lond0,latd0,itc), ! output updated vn at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
        if (debug) write(6,"('dynamics after duv')")
!     
      else  				! GSWM input used 
!
! Advance neutral temperature:
!
        call dt_gswm(
     |    tn	(levd0,lond0,latd0,itp),
     |    tn_nm (levd0,lond0,latd0,itp),
     |    un	(levd0,lond0,latd0,itp),
     |    vn	(levd0,lond0,latd0,itp),
     |    o2	(levd0,lond0,latd0,itp),
     |    o1	(levd0,lond0,latd0,itp),
     |    barm  (levd0,lond0,latd0,itp),    ! mean molecular weight
     |    cp	(levd0,lond0,latd0),	    ! specific heat	   (cpktkm.F)
     |    kt	(levd0,lond0,latd0),	    ! molecular diffusion  (cpktkm.F)
     |    km	(levd0,lond0,latd0),	    ! molecular viscosity  (cpktkm.F)
     |    hdt	(levd0,lond0,latd0),	    ! horizontal diffusion (hdif.F)
     |    qji_tn(levd0,lond0,latd0),	    ! joule heating (sub qjoule_tn)
     |    cool_implicit(levd0,lond0,latd0), ! implicit cooling (newton.F)
     |    cool_explicit(levd0,lond0,latd0), ! explicit cooling (newton.F)
     |    w	(levd0,lond0,latd0,itc),    ! updated W (swdot.F)
     |    tn	(levd0,lond0,latd0,itc),    ! output updated tn
     |    tn_nm (levd0,lond0,latd0,itc),    ! output updated tn at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
        if (debug) write(6,"('dynamics after dt_gswm')")
!
! Advance neutral velocities:
!
        call duv_gswm(
     |    tn	(levd0,lond0,latd0,itp), ! neutral temperature
     |    tn	(levd0,lond0,latd0,itc), ! updated neutral temperature (dt.F)
     |    tn_nm (levd0,lond0,latd0,itp), ! tn at time n-1
     |    un	(levd0,lond0,latd0,itp), ! zonal velocity
     |    vn	(levd0,lond0,latd0,itp), ! meridional velocity
     |    un_nm (levd0,lond0,latd0,itp), ! zonal velocity at time n-1
     |    vn_nm (levd0,lond0,latd0,itp), ! meridional velocity at time n-1
     |    w	(levd0,lond0,latd0,itc), ! updated vertical velocity (swdot.F)
     |    barm  (levd0,lond0,latd0,itp), ! mean molecular weight
     |    z	(levd0,lond0,latd0,itp), ! geopotential height
     |    hdu	(levd0,lond0,latd0),	 ! horizontal diffusion of U (hdif3.F)
     |    hdv	(levd0,lond0,latd0),	 ! horizontal diffusion of V (hdif3.F)
     |    ui	(levd0,lond0,latd0),	 ! zonal ion velocity input
     |    vi	(levd0,lond0,latd0),	 ! meridional ion velocity input
     |    lxx	(levd0,lond0,latd0),	 ! xx ion drag coefficient
     |    lyy	(levd0,lond0,latd0),	 ! yy ion drag coefficient
     |    lxy	(levd0,lond0,latd0),	 ! xy ion drag coefficient
     |    lyx	(levd0,lond0,latd0),	 ! yx ion drag coefficient
     |    km	(levd0,lond0,latd0),	 ! molecular viscosity  (cpktkm.F)
     |    un	(levd0,lond0,latd0,itc), ! output updated un
     |    un_nm (levd0,lond0,latd0,itc), ! output updated un at time n-1
     |    vn	(levd0,lond0,latd0,itc), ! output updated vn
     |    vn_nm (levd0,lond0,latd0,itc), ! output updated vn at time n-1
     |    1,nlevp1,lon0,lon1,lat0,lat1)
        if (debug) write(6,"('dynamics after duv_gswm')")
      
      endif	! endif for GSWM input
!
! Sources and sinks for major species composition O2, O:
!     subroutine comp_o2o(tn,o2,o1,barm,op,no,n4s,n2d,o2p,ne,
!    |  n2p,nplus,nop,lev0,lev1,lon0,lon1,lat)
!
      do lat=lat0,lat1
        call comp_o2o(
     |    tn   (levd0,lond0,lat,itp), ! neutral temperature
     |    o2   (levd0,lond0,lat,itp), ! molecular oxygen (mmr)
     |    o1   (levd0,lond0,lat,itp), ! atomic oxygen (mmr)
     |    barm (levd0,lond0,lat,itp), ! mean molecular weight
     |    op   (levd0,lond0,lat,itp), ! O+
     |    no   (levd0,lond0,lat,itp), ! NO
     |    n4s  (levd0,lond0,lat,itp), ! N(4S)
     |    n2d  (levd0,lond0,lat,itp), ! N(2D)
     |    o2p  (levd0,lond0,lat,itp), ! O2+
     |    ne   (levd0,lond0,lat,itp), ! electron density
     |    n2p  (levd0,lond0,lat),     ! N2+
     |    nplus(levd0,lond0,lat),     ! N+
     |    nop  (levd0,lond0,lat),     ! NO+
     |    1,nlevp1,lon0,lon1,lat)
      enddo ! lat=lat0,lat1
      if (debug) write(6,"('dynamics after comp_o2o')")
!
! Advance O2, O:
!     subroutine comp(tn,o2,o2_nm,o1,o1_nm,un,vn,w,hdo2,hdo1,
!    |                o2_upd,o2nm_upd,o1_upd,o1nm_upd,
!    |                lev0,lev1,lon0,lon1,lat0,lat1)
!
      call comp(
     |  tn    (levd0,lond0,latd0,itp), ! neutral temperature
     |  o2    (levd0,lond0,latd0,itp), ! O2 (mmr)
     |  o2_nm (levd0,lond0,latd0,itp), ! O2 (mmr) at time n-1
     |  o1    (levd0,lond0,latd0,itp), ! O1 (mmr)
     |  o1_nm (levd0,lond0,latd0,itp), ! O1 (mmr) at time n-1
     |  un    (levd0,lond0,latd0,itp), ! zonal velocity
     |  vn    (levd0,lond0,latd0,itp), ! meridional velocity
     |  w     (levd0,lond0,latd0,itp), ! vertical velocity
     |  hdo2  (levd0,lond0,latd0),     ! O2 horizontal diffusion (hdif3)
     |  hdo1  (levd0,lond0,latd0),     ! O  horizontal diffusion (hdif3)
     |  o2    (levd0,lond0,latd0,itc), ! output: O2 updated for current step 
     |  o2_nm (levd0,lond0,latd0,itc), ! output: O2 updated for previous step
     |  o1    (levd0,lond0,latd0,itc), ! output: O1 updated for current step
     |  o1_nm (levd0,lond0,latd0,itc), ! output: O1 updated for previous step
     |  1,nlevp1,lon0,lon1,lat0,lat1)
      if (debug) write(6,"('dynamics after comp')")

      if (debug) write(6,"('dynamics returning')")
      end subroutine dynamics

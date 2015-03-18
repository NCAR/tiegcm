!
      module chemrates_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlevp1,nlonp4
!
! Chemical reaction rates.
!
! Temperature independent rates are real parameters (constants),
!   set at compile time.
! Temperature dependent rates are dimensioned (nlev,nlon) and are
!   updated by sub rates_tdep at every latitude of every time step. 
! See refs:
! Roble,R.G., 1995. Energetics of the Mesosphere and Thermosphere,
!   Geophysical Monograph 87, American Geophysical Union
! R.G. Roble, E.C. Ridley
!   An auroral model for the NCAR thermospheric general circulation model,
!   Annales Geophysicae,5A, (6), 369-382, 1987. 
!
      implicit none 
!
! Temperature independent reaction rate constants:
      real,parameter ::
!
! Ion chemistry:
     |  rk4   = 1.0E-10,  ! O2+ + N4S -> NO+ + O   + 4.21 eV
     |  rk5   = 4.4E-10,  ! O2+ + NO  -> NO+ + O2  + 2.813 eV
     |  rk6   = 4.0E-10,  ! N+  + O2  -> O2+ + N4S + 2.486 eV
     |  rk7   = 2.0E-10,  ! N+  + O2  -> NO+ + O   + 6.669 eV
     |  rk8   = 1.0E-12,  ! N+  + O   -> O+  + N   + 0.98 eV
     |  rk9   = 6.0E-11,  ! N2+ + O2  -> O2+ + N2  + 3.52 eV
     |  rk10  = 1.3E-10,  ! O+  + N2D -> N+  + O   + 1.45 eV
     |  rk16  = 4.8E-10,  ! O+(2P) + N2 -> N2+ + O  + 3.02 eV
     |  rk17  = 1.0E-10,  ! O+(2P) + N2 -> N+  + NO + 0.70 eV
     |  rk18  = 4.0E-10,  ! O+(2P) + O  -> O+  + O  + 5.2 eV
     |  rk21  = 0.047,    ! O+(2P)      -> O+  + hv (2470A)
     |  rk22  = 0.171,    ! O+(2P)      -> O+  + hv (7320A)
     |  rk23  = 8.E-10,   ! O+(2D)      -> N2+O     + 1.33 eV
     |  rk24  = 5.0E-12,  ! O+(2D) + O  -> O+(4S) +e+ 3.31 eV
     |  rk26  = 7.E-10,   ! O+(2D) + O2 -> O2+O     + 4.865 eV
     |  rk27  = 7.7E-5,
!
! Neutral chemistry:
     |  beta2 = 5.0E-12,  ! N2D + O2 -> NO  + O1D + 1.84 eV
     |  beta4 = 7.0E-13,  ! N2D + O  -> N4S + O   + 2.38 eV
     |  beta6 = 7.0E-11,  ! N2D + NO -> N2  + O   + 5.63 eV
     |  beta7 = 1.06E-5   ! N2D      -> N4S + hv
!
! Temperature dependent chemical reaction rate coefficients. 
! These are set by rates_tdep for every latitude at every timestep. 
!
!     real,dimension(nlevp1,nlonp4) :: 
!
      real,dimension(:,:,:),allocatable :: ! (nlevp1,lon0:lon1,lat0:lat1)
!
! Ion chemistry:
     |  rk1, ! O+  + O2 -> O + O2+ + 1.555 eV
     |  rk2, ! O+  + N2 -> NO+ + N4s + 1.0888 eV
     |  rk3, ! N2+ + O  -> NO+ + N2D + 0.70 eV
     |  ra1, ! NO+ + e  -> (20% N4S + O + 2.75 eV), (80% N2D + O + 0.38 eV)
     |  ra2, ! O2+ + e  -> (15% O   + O + 6.95 eV), (85% O   + O + 4.98 eV)
     |  ra3, ! N2+ + e  -> (10% N4S + N4S + 5.82 eV), (90% N4S + N2D + 3.44 eV)
!
! Neutral chemistry:
     |  beta1, ! N4S + O2 -> NO  + O + 1.4  eV
     |  beta3, ! N4S + NO -> N2  + O + 3.25 eV
     |  beta5, ! N2D + e  -> N4S + e + 2.38 eV
     |  beta8, ! NO  + hv -> N4S + O
     |  beta9, ! NO  + hv/Ly-a -> NO+ + e
     |  beta9n,
     |  beta17,

     |  rk19,  ! O+(2P) + e -> O+(4S) + e + 5.0  eV
     |  rk20,  ! O+(2P) + e -> O+(2O) + e + 1.69 eV
     |  rk25,  ! O+(2D) + e -> O+(4S) + e + 3.31 eV
     |  rkm12, ! O + O + N2 -> O2 + N2
     |  tvib
!
! Matrices for O2,O (set by comp_o2o.F, used by comp.F) 
      real,allocatable :: fs(:,:,:,:,:) ! (i,k,3,0:3,j)
!
      contains
!-----------------------------------------------------------------------
      subroutine alloc_tdep
!
! Allocate temperature-dependent reaction rates for task subdomain: 
! Called once per run from init_fields.
!
      use mpi_module,only: lon0,lon1,lat0,lat1
      integer :: istat
!
      allocate(rk1(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk1: stat=',i3)") istat
      allocate(rk2(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk2: stat=',i3)") istat
      allocate(rk3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk3: stat=',i3)") istat
      allocate(ra1(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' ra1: stat=',i3)") istat
      allocate(ra2(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' ra2: stat=',i3)") istat
      allocate(ra3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' ra3: stat=',i3)") istat
      allocate(beta1(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' beta1: stat=',i3)") istat
      allocate(beta3(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' beta3: stat=',i3)") istat
      allocate(beta5(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' beta5: stat=',i3)") istat
      allocate(beta8(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' beta8: stat=',i3)") istat
      allocate(beta9(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' beta9: stat=',i3)") istat
      allocate(beta9n(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' beta9n: stat=',i3)") istat
      allocate(beta17(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' beta17: stat=',i3)") istat
      allocate(rk19(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk19: stat=',i3)") istat
      allocate(rk20(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk20: stat=',i3)") istat
      allocate(rk25(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rk25: stat=',i3)") istat
      allocate(rkm12(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' rkm12: stat=',i3)") istat
      allocate(tvib(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' tvib: stat=',i3)") istat
!
      allocate(fs(lon0:lon1,nlevp1,3,0:3,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_tdep: error allocating',
     |  ' fs: stat=',i3)") istat

      end subroutine alloc_tdep
!-----------------------------------------------------------------------
      subroutine chemrates_tdep(tn,te,ti,fno2,fnvo2,lev0,lev1,lon0,lon1,
     |  lat)
!
! Calculate temperature-dependent reaction rates (called at each latitude)
!
      use input_module,only: f107  ! 10.7 cm flux (from input and/or gpi)
      use init_module,only: sfeps  ! flux variation from orbital excentricity
      use cons_module,only: check_exp
!
! Args:
      integer,intent(in) ::
     |  lev0,lev1,           ! first and last level indices, this task
     |  lon0,lon1,           ! first and last longitude indices, this task
     |  lat                  ! latitude index
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,     ! neutral temperature (deg K)
     |  te,     ! electron temperature (deg K)
     |  ti,     ! ion temperature (deg K)
     |  fno2,   ! o2 line integral (see chapman.F)
     |  fnvo2   ! o2 column number density (see chapman.F)
!
! Local:
      integer :: k,i
      real :: ti1(lev0:lev1,lon0:lon1), 
     |        ti2(lev0:lev1,lon0:lon1), 
     |        ti3(lev0:lev1,lon0:lon1),
     |        etvib(lev0:lev1,lon0:lon1)
!
! expo() (util.F) is used only if check_exp is true. This will avoid
! NaNS fpe, but will degrade performance. Check_exp is in cons.F.
      real,external :: expo
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!         write(6,"('chemrates_tdep: k=',i2,' i=',i2,' lat=',i2)")
!    |      k,i,lat
!
! ti1 = T1/300. (S15)    ti2 = T2/300. (S14)   ti3 = TR/300. (S13)
!
          ti1(k,i) = (0.667* ti(k,i)+0.333 *tn(k,i))/300.
          ti2(k,i) = (0.6363*ti(k,i)+0.3637*tn(k,i))/300.
          ti3(k,i) = .5*(ti(k,i)+tn(k,i))/300.
!
! Rate coefficient for O+ + O2 -> O + O2+ + 1.555eV (HIERL)
          rk1(k,i,lat)=1.6E-11*ti1(k,i)**(-0.52)+5.5E-11*
     |      exp(-22.85/ti1(k,i))
!
! k2: O+  + N2 -> NO+ + N4s + 1.0888 eV
! (assumes t2/300 in ti2)
          rk2(k,i,lat)=(8.6E-14*ti2(k,i)-5.92E-13)*ti2(k,i)+1.533E-12
!
! Rate coefficient for O+ + N2 -> N + NO+
! See also sub altv (altv.f). Tvib is output of sub altv, if it
! is called (see dynamics.f).
!
          tvib(k,i,lat) = tn(k,i)            ! tvib in crates_tdep.h
          etvib(k,i) = exp(-3353./tvib(k,i,lat)) ! etvib is local
!
! rk2: O+  + N2 -> NO+ + N4s + 1.0888 eV
          rk2(k,i,lat) = (((((270.*etvib(k,i)+220.)*etvib(k,i)+85.)*
     |      etvib(k,i)+38.)*etvib(k,i)+1.)*rk2(k,i,lat)*etvib(k,i)+
     |      rk2(k,i,lat))*(1.-etvib(k,i))
!
! rk3: N2+ + O  -> NO+ + N2D + 0.70 eV
          if (300.*ti3(k,i) >= 1500.) then
            rk3(k,i,lat)=5.2E-11*ti3(k,i)**0.2
          else
            rk3(k,i,lat)=1.4E-10*ti3(k,i)**(-0.44)
          endif
!
! rk19: O+(2P) + e -> O+(4S) + e + 5.0  eV
          rk19(k,i,lat) = 4.0E-8*sqrt(300./te(k,i))
!
! rk20: O+(2P) + e -> O+(2O) + e + 1.69 eV
          rk20(k,i,lat) = 1.5E-7*sqrt(300./te(k,i))
!
! rk25: O+(2D) + e -> O+(4S) + e + 3.31 eV
          rk25(k,i,lat) = 6.6E-8*sqrt(300./te(k,i))
!
! rkm12: O + O + N2 -> O2 + N2 (loss of O to O2)
! (same for O + O + O2 -> O2 + O2, Dickinson, et.al., 1984)
          rkm12(k,i,lat) = 9.59E-34*exp(480./tn(k,i))
!
! ra1: NO+ + e  -> (20% N4S + O + 2.75 eV), (80% N2D + O + 0.38 eV)
          ra1(k,i,lat)=4.2E-7*(300./te(k,i))**0.85
!
! ra2: O2+ + e  -> (15% O   + O + 6.95 eV), (85% O   + O + 4.98 eV)
          if (te(k,i) >= 1200.) then
            ra2(k,i,lat)=1.6E-7*(300./te(k,i))**0.55
          else
! ra2: notes replace 2.7e-7 with 1.95e-7
            ra2(k,i,lat)=2.7E-7*(300./te(k,i))**0.7 
          endif
!
! ra3: N2+ + e  -> (10% N4S + N4S + 5.82 eV), (90% N4S + N2D + 3.44 eV)
          ra3(k,i,lat)=1.8E-7*(300./te(k,i))**0.39
!
! beta1: N4S + O2 -> NO  + O + 1.4  eV
          beta1(k,i,lat) = 1.5E-11*exp(-3600./tn(k,i))
!
! beta3: N4S + NO -> N2  + O + 3.25 eV
          beta3(k,i,lat) =3.4e-11*sqrt(tn(k,i)/300.)
!
! beta5: N2D + e  -> N4S + e + 2.38 eV
          beta5(k,i,lat) = 3.6E-10*sqrt(te(k,i)/300.)
!
! beta8: NO  + hv -> N4S + O
          if (.not.check_exp) then
            beta8(k,i,lat) = 4.5E-6*(1.+0.11*(f107-65.)/165.)*
     |        exp(-1.E-8*fno2(k,i)**0.38)*sfeps
          else
            beta8(k,i,lat) = 4.5E-6*(1.+0.11*(f107-65.)/165.)*
     |        expo(-1.E-8*fno2(k,i)**0.38,0)*sfeps
          endif
!
! beta9: NO  + hv/Ly-a -> NO+ + e
          if (.not.check_exp) then
            beta9(k,i,lat)=2.91E11*(1.+0.2*(f107-65.)/100.)*2.E-18*
     |        exp(-8.E-21*fno2(k,i))*sfeps
          else
            beta9(k,i,lat)=2.91E11*(1.+0.2*(f107-65.)/100.)*2.E-18*
     |        expo(-8.E-21*fno2(k,i),0)*sfeps
          endif
          beta9n(k,i,lat)=5.E9*(1.+0.2*(f107-65.)/100.)*2.E-18*
     |      exp(-8.E-21*fnvo2(k,i))*sfeps
          beta17(k,i,lat)=1.E-32*sqrt(300./te(k,i))
        enddo ! i=lon0,lon1
      enddo ! k=lev0,lev1
!
! Upper boundary:
      do i=lon0,lon1

!       write(6,"('chemrates: i=',i3,' lon0,1=',2i3,' f107=',e12.4,
!    |    ' fno2(lev1,i)=',e12.4)") i,lon0,lon1,f107,fno2(lev1,i)

        beta8(lev1,i,lat)=4.5E-6*(1.+0.11*(f107-65.)/165.)*
     |    exp(-1.E-8*fno2(lev1,i)**0.38)
        beta9(lev1,i,lat)=2.91E11*(1.+0.2*(f107-65.)/100.)*2.E-18*
     |    exp(-8.E-21*fno2(lev1,i))
        beta9n(lev1,i,lat)=5.E9*(1.+0.2*(f107-65.)/100.)*2.E-18*
     |    exp(-8.E-21*fnvo2(lev1,i))
      enddo
      end subroutine chemrates_tdep
      end module chemrates_module

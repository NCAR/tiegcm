!
      subroutine dt(tn,tn_nm,un,vn,o2,o1,n2,barm,cp,kt,km,hdt,qji_tn,
     |  cool_imp,cool_exp,w_upd,tn_upd,tn_nm_upd,
     |  lev0,lev1,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Advance neutral temperature at current latitude: 
! 4/4/05 btf: added gswm lbc option (eliminated dt_gswm.F)
!
      use params_module,only: nlonp4,dz,nlat,spval
      use input_module,only: step
      use init_module,only: iter,igetgswm,iday
      use cons_module,only: freq_semidi,tbound,shapiro,dtx2inv,expz,
     |  rmassinv_o2,rmassinv_o1,rmassinv_n2,tsurplus,p0,boltz,avo,grav,
     |  gask,expzmid,expzmid_inv,dift,dtsmooth,dtsmooth_div2,kut
      use lbc,only: t_lbc
      use qrj_module,only: qtotal ! qtotal(nlevp1,lon0:lon1,lat0:lat1)
      use chemrates_module,only: rkm12
      use fields_module,only: tlbc,ulbc,vlbc,tlbc_nm
      use addfld_module
      use diags_module,only: mkdiag_DEN,mkdiag_HEAT
#ifdef MPI
      use mpi_module,only: mp_bndlons_f3d, mp_periodic_f3d
#endif
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Full subdomains:
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),intent(in)::
     |  tn,    ! neutral temperature (deg K)
     |  tn_nm, ! neutral temperature, time n-1 
     |  un,    ! neutral zonal velocity (cm/sec)
     |  vn,    ! neutral zonal velocity (cm/sec)
     |  o2,    ! molecular oxygen (mmr)
     |  o1,    ! atomic oxygen (mmr)
     |  n2,    ! molecular nitrogen (mmr)
     |  barm,  ! mean molecular weight
     |  cp,    ! specific heat (ergs/deg/gm)           (sub cpktkm)
     |  kt,    ! molecular diffusion (ergs/cm/deg/sec) (sub cpktkm)
     |  km,    ! molecular viscosity (gm/cm/sec)       (sub cpktkm)
     |  hdt,   ! horizontal diffusion of tn (from sub hdif3, hdif.F)
     |  qji_tn,! joule heating for tn (from sub qjoule_tn, qjoule.F)
     |  cool_imp, ! implicit cooling (newton.F)
     |  cool_exp, ! explicit cooling (newton.F)
     |  w_upd  ! updated vertical velocity (swdot.F)
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) ::
     |  tn_upd,   ! updated tn (output)
     |  tn_nm_upd ! updated tn at time n-1 (output)
!
! Local:
      integer :: k,i,lonbeg,lonend,lat
      integer :: nk,nkm1,nlevs
      complex :: expt,expt2
      real :: rstep
      real :: tnlbc(lon0:lon1,lat0:lat1) ! lower boundary condition
!
! Local at 2d:
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  cptn,       ! cp*exp(-s)*(V.del(T(n))-1/(2*DT))*T(n-1) (k+1/2)
     |  mbar,       ! mean molecular weight
     |  qm,         ! heating due to molecular diffusion
     |  total_heat, ! total heating
     |  dudz,       ! du/dz(k)
     |  dvdz,       ! du/dz(k)
     |  g,          ! g*KT/(p0*H*Ds**2)
     |  f,          ! g*eps/(p0*2*H*Ds)
     |  h,          ! scale height R*T/(M*g) (cm)
     |  rho,        ! density
     |  tni,        ! tn at interfaces
     |  p, q, r,    ! coefficients for tridiagonal solver
     |  rhs,        ! right hand side of trsolv
     |  qpart,      ! part of q coeff
     |  tnlbc_diag  ! tnlbc_diag is redundant in vertical
!
! Local at 3d (tnsmooth needs lat dimension only for sub smooth):
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  tnsmooth,   ! zonal and meridional smoothing of tn_nm
     |  advec_tn    ! horizontal advection (output of sub advec)
!
      nk = lev1-lev0+1
      nkm1 = nk-1
      nlevs = nk
! 
! First latitude scan for dt:
      do lat=lat0,lat1
!
! Lower boundary t_lbc was calculated by sub tuvz_lbc (lbc.F)
        tnlbc(:,lat) = t_lbc(lon0:lon1,lat)
        do k=lev0,lev1
          tnlbc_diag(k,:) = tnlbc(:,lat)
        enddo
!       call addfld('TNLBC1',' ',' ',tnlbc_diag,
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Horizontal advection (pass k vs i slices at full task subdomain 
! longitudes, and the 5 latitudes centered over the current latitude).
!
        call advec(tn(:,:,lat-2:lat+2),advec_tn(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat)

!       call addfld('HADVECTN',' ',' ',advec_tn(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Vertical advection. Sub advecv adds vertical advection to advec_tn.
        call advecv(tn(:,:,lat),tnlbc(:,lat),advec_tn(:,:,lat),
     |    lev0,lev1,lon0,lon1,lat)

!       call addfld('ADVEC_TN',' ',' ',advec_tn(lev0:lev1-1,:,lat),
!    |    'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! End first latitude scan:
      enddo ! lat=lat0,lat1
!
! Shapiro smoother for tn at time n-1:
      call smooth(tn_nm,tnsmooth,lev0,lev1,lon0,lon1,lat0,lat1,0)
!
! Begin second latitude scan:
      do lat=lat0,lat1
!
! Set cptn and mbar (k+1/2):
! (Earlier versions apparently assumed zero periodic points for
!  tnsmooth, since they were not set in smoothing. See sub smooth,
!  where the periodic points are set to zero to avoid NaNS fpe
!  in the following loop)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          cptn(k,i) = .5*(cp(k,i,lat)+cp(k+1,i,lat))*expz(k)*
     |      (advec_tn(k,i,lat)-dtx2inv*tnsmooth(k,i,lat))
!
! 2/28/05 btf: Use k rather than k+1 in mbar calculation
!              (as in recent versions of time-gcm, e.g., timegcm1.2)
          mbar(k,i) = 1./(o2(k,i,lat)*rmassinv_o2 +
     |      o1(k,i,lat)*rmassinv_o1 + n2(k,i,lat)*rmassinv_n2)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('CP'   ,' ',' ',cp(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('TNSMOOTH'   ,' ',' ',tnsmooth(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('ADVEC_TNa',' ',' ',advec_tn(lev0:lev1-1,:,lat),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('CPTN0',' ',' ',cptn,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('MBAR' ,' ',' ',mbar(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Total heat sources are in total_heat (s5).
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! Solar heating from qrj:
          total_heat(k,i) = .5*(qtotal(k,i,lat)+qtotal(k+1,i,lat))
!
! Add heating from 4th order horizontal diffusion (hdt from sub hdif3):
          total_heat(k,i) = total_heat(k,i)+hdt(k,i,lat)
!
! Add heating due to atomic oxygen recombination:
          total_heat(k,i) = total_heat(k,i)+tsurplus*rkm12(k,i,lat)*
     |      (p0*expz(k)*mbar(k,i)/(boltz*tn(k,i,lat))*o1(k,i,lat)*
     |      rmassinv_o1)**2*avo/mbar(k,i)
!
! Add ion joule heating (from sub qjoule_tn, qjoule.F)
! (Note joulefac was applied in sub qji_tn)
!
	  total_heat(k,i) = total_heat(k,i)+qji_tn(k,i,lat)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!     call addfld('HEATING1',' ',' ',total_heat(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!
! Add heating due to molecular diffusion:
! du/dz and dv/dz (s10, s11):
      do i=lon0,lon1
        do k=lev0+1,lev1-2
          dudz(k,i) = (un(k+1,i,lat)-un(k-1,i,lat))/(2.*dz)  ! s10
          dvdz(k,i) = (vn(k+1,i,lat)-vn(k-1,i,lat))/(2.*dz)  ! s11
        enddo ! k=lev0+1,lev1-2
!
! Lower boundary is in ulbc,vlbc: 
        dudz(1,i) = (un(1,i,lat)+1./3.*un(2,i,lat)-4./3.*
     |    ulbc(i,lat))/dz
        dvdz(1,i) = (vn(1,i,lat)+1./3.*vn(2,i,lat)-4./3.*
     |    vlbc(i,lat))/dz
!
! Upper boundary:
        dudz(lev1-1,i) = dudz(lev1-2,i)/3.
        dvdz(lev1-1,i) = dvdz(lev1-2,i)/3.
!
! qm = heating due to molecular diffusion:
! (km = molecular viscosity from sub cpktkm)
        do k=lev0,lev1-1
          qm(k,i) = grav**2*mbar(k,i)*.5*(km(k,i,lat)+km(k+1,i,lat))/
     |      (p0*gask*expz(k)*tn(k,i,lat))*(dudz(k,i)**2+dvdz(k,i)**2)
!
! Add qm to total heating:
          total_heat(k,i) = total_heat(k,i)+qm(k,i)
!
! Complete cptn:
! -cp*exp(-s)*(T(k,n-1)/(2*Dt) - V.del(T(k,n)) +Q/cp)
!
          cptn(k,i) = cptn(k,i)-expz(k)*total_heat(k,i) ! s1
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('CPTN',' ',' ',cptn,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QM'     ,' ',' ',qm(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('HEATING2',' ',' ',total_heat(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

      call mkdiag_HEAT('HEATING',total_heat,lev0,lev1,lon0,lon1,lat)
!
! H = R*T/(M*g)                            (s4)
! rho = p0*exp(-s)*M/(R*T)                 (s5)
! tni = T                                  (s6)
!
! Levels 2 through lev1-1:
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          tni(k,i) = .5*(tn(k-1,i,lat)+tn(k,i,lat))
          h(k,i) = gask*tni(k,i)/barm(k,i,lat)
          rho(k,i) = p0*expzmid_inv*expz(k)/h(k,i)
          h(k,i) = h(k,i)/grav
        enddo ! k=lev0+1,lev1-1
!
! Boundaries:
        tni(lev0,i) = tlbc(i,lat)         ! Lower boundary is in tlbc = tn(itp)
        tni(lev1,i) = tn(lev1-1,i,lat)
        h(lev0,i) = gask*tni(lev0,i)/barm(lev0,i,lat)
        h(lev1,i) = gask*tni(lev1,i)/barm(lev1,i,lat)
        rho(lev0,i) = p0*expzmid_inv*expz(lev0)/h(lev0,i)
        rho(lev1,i) = p0*expzmid_inv*expz(lev1-1)/h(lev1,i)
        h(lev0,i) = h(lev0,i)/grav
        h(lev1,i) = h(lev1,i)/grav
!
! G = g*(kT + H**2*rho*cp*kE)/(p0*H*Ds**2) (s2)
! F = g*(kE*H**3*rho*g/T)/(p0*2*H*Ds)      (s3)
!
        do k=lev0,lev1-1 
          g(k,i) = grav*(kt(k,i,lat)+h(k,i)**2*rho(k,i)*cp(k,i,lat)*
     |      dift(k,iday))/(p0*h(k,i)*dz**2)
          f(k,i)=grav**2*dift(k,iday)*h(k,i)**2*rho(k,i)/(tni(k,i)*
     |      p0*2.*dz)
        enddo ! k=lev0,lev1-1 
      enddo ! i=lon0,lon1

!     call addfld('TNI',' ',' ',tni,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('HG' ,' ',' ',h  ,'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('DEN','Neutral Density','g/cm3',rho,
!    |  'ilev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('G'  ,' ',' ',g(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('F'  ,' ',' ',f(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)

! Save diagnostic density:
      call mkdiag_DEN('DEN',rho,lev0,lev1,lon0,lon1,lat)
!
! Coefficients for trsolv:
! Levels 3/2 through K-3/2
      do i=lon0,lon1
        do k=lev0,lev1-2
          p(k,i) = g(k,i)-f(k,i)
          q(k,i) = -g(k,i)-g(k+1,i) - f(k,i)+f(k+1,i)
          r(k,i) = g(k+1,i) + f(k+1,i)
          rhs(k,i) = cptn(k,i)
        enddo ! k=lev0,lev1-2
! Level k-1/2
        p(lev1-1,i) =  g(lev1-1,i)-f(lev1-1,i)
        q(lev1-1,i) = -g(lev1-1,i)-f(lev1-1,i)
        r(lev1-1,i) = 0.
        rhs(lev1-1,i) = cptn(lev1-1,i)
      enddo ! i=lon0,lon1

!     call addfld('P_COEF0' ,' ',' ',p,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF0' ,' ',' ',q,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('R_COEF0' ,' ',' ',r,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RHS0'    ,' ',' ',rhs,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COOL_IMP',' ',' ',cool_imp(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COOL_EXP',' ',' ',cool_exp(:,lon0:lon1,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! qpart = cp*(1/(2*Dt)+ai/cp+w*R/(cp*M))
      do i=lon0,lon1
        do k=lev0,lev1-1
          qpart(k,i) = 
     |      .5*(cp(k,i,lat)+cp(k+1,i,lat))*(dtx2inv+cool_imp(k,i,lat))+
     |      .5*(w_upd(k,i,lat)+w_upd(k+1,i,lat))*gask/mbar(k,i) 
          rhs(k,i) = rhs(k,i)+cool_exp(k,i,lat)
          q(k,i) = q(k,i)-expz(k)*qpart(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('QPART'  ,' ',' ',qpart(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF1',' ',' ',q    ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RHS1'   ,' ',' ',rhs  ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     do k=lev0,lev1
!       tnlbc_diag(k,:) = tnlbc(:,lat)
!     enddo
!     call addfld('TNLBC2',' ',' ',tnlbc_diag,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Lower boundary:
      do i=lon0,lon1
        q(lev0,i) = q(lev0,i)-p(lev0,i)
!
! Diffs in rhs lbc ??:
        rhs(lev0,i) = rhs(lev0,i)-2.*p(lev0,i)*tnlbc(i,lat)
        p(lev0,i) = 0.
      enddo ! i=lon0,lon1

!     call addfld('P_COEF2',' ',' ',p,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('Q_COEF2',' ',' ',q,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('R_COEF2',' ',' ',r,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('RHS2'   ,' ',' ',rhs,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Solve tridiagonal system for new tn:
!     subroutine trsolv(a,b,c,f,x,lev0,lev1,k1,k2,lon0,lon1,lonmax,lat,
!    |  idebug)
!
      call trsolv(p,q,r,rhs,tn_upd(:,lon0:lon1,lat),lev0,lev1,
     |  lev0,lev1-1,lon0,lon1,nlonp4,lat,0)

!     call addfld('TN_SOLV','Updated TN from trsolv','K',
!    |  tn_upd(:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End second latitude scan:
      enddo ! lat=lat0,lat1
!
! Filter updated tn:
      call filter_tn(tn_upd,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
! Third latitude scan:
      do lat=lat0,lat1
!       call addfld('TN_FILT',' ',' ',tn_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Smooth updated tn:
        do i=lon0,lon1 
          do k=lev0,lev1-1
            tn_nm_upd(k,i,lat) = dtsmooth_div2*(tn_nm(k,i,lat)+
     |        tn_upd(k,i,lat)) + dtsmooth*tn(k,i,lat)
          enddo ! k=lev0,lev1-1
        enddo ! i=lon0,lon1 
!       call addfld('TN_NMOUT',' ',' ',tn_nm_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Put spval in top nlevp1 level:
        tn_upd(lev1,lon0:lon1,lat)    = spval
        tn_nm_upd(lev1,lon0:lon1,lat) = spval
!
! Lower boundary is saved in tlbc (fields.F):
        tlbc_nm(lon0:lon1,lat) = tlbc(lon0:lon1,lat)
        tlbc(lon0:lon1,lat)    = tnlbc(:,lat)
!
! Tn must be at least 100 deg:
        lonbeg = lon0-2
        if (lon0==1) lonbeg = 1
        lonend = lon1+2
        if (lon1==nlonp4) lonend = nlonp4
        do i=lonbeg,lonend
          do k=lev0,lev1
            if (tn_upd(k,i,lat) < 100.) tn_upd(k,i,lat) = 100.
          enddo
        enddo
!       call addfld('TN_FINAL',' ',' ',tn_upd(:,lon0:lon1,lat),
!    |    'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! End third lat scan:
      enddo ! lat=lat0,lat1
      end subroutine dt
!-----------------------------------------------------------------------
      subroutine filter_tn(fout,lev0,lev1,lon0,lon1,lat0,lat1,kut)
!
! Filter updated TN
!
      use params_module,only: nlat,nlonp4,nlon
      use filter_module,only: filter
#ifdef MPI
      use mpi_module,only: mp_gatherlons_f3d,mp_scatterlons_f3d,mytidi
      implicit none
#else
      implicit none
      integer :: mytidi=0
#endif
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1,kut(nlat)
      real,intent(inout) :: fout(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2)
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
! Local:
      integer :: i,j,nlevs,nlons,nlats
      real :: fik(nlonp4,lev0:lev1),fkij(lev0:lev1,nlonp4,lat0:lat1)
      real :: fmin,fmax
!
#ifdef VT
!     code = 131 ; state = 'filter_tn' ; activity='Filtering'
      call vtbegin(131,ier)
#endif
!
      nlevs = lev1-lev0+1
      nlons = lon1-lon0+1
      nlats = lat1-lat0+1
!
! Define lons in fkij from current task:
      fkij = 0.
      do j=lat0,lat1
        do i=lon0,lon1
          fkij(:,i,j) = fout(:,i,j)
        enddo
      enddo ! j=lat0,lat1
!
#ifdef MPI
!
! Gather longitudes into tasks in first longitude column of task table
!   (leftmost of each j-row) for global fft. (i.e., tasks with mytidi==0 
!   gather lons from other tasks in that row). This includes all latitudes.
!
      call mp_gatherlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1,'TN')
#endif
!
! Only leftmost tasks at each j-row of tasks does the global filtering:
      if (mytidi==0) then
!
! Define 2d array with all longitudes for filter at each latitude:
        latscan: do j=lat0,lat1
          if (kut(j) >= nlon/2) cycle latscan
          do i=1,nlonp4
            fik(i,:) = fkij(:,i,j)
          enddo ! i=1,nlonp4
!
! Remove wave numbers > kut(lat):
          call filter(fik,lev0,lev1,kut(j),j)
!
! Return filtered array to fkij:
          do i=1,nlonp4
            fkij(:,i,j) = fik(i,:)
          enddo ! i=1,nlonp4
        enddo latscan ! j=lat0,lat1
      endif ! mytidi==0
#ifdef MPI
!
! Now leftmost task at each j-row must redistribute filtered data
! back to other tasks in the j-row (mytidi>0,mytidj) (includes latitude):
!
      call mp_scatterlons_f3d(fkij,lev0,lev1,lon0,lon1,lat0,lat1,1,'TN')
#endif
!
! Return filtered array to fout at current task longitudes and latitudes:
      do j=lat0,lat1
        do i=lon0,lon1
          fout(:,i,j) = fkij(:,i,j)
        enddo
      enddo
!
#ifdef VT
!     code = 131 ; state = 'filter_tn' ; activity='Filtering'
      call vtend(131,ier)
#endif
      end subroutine filter_tn

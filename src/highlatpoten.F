! BOP
! !IROUTINE: set_zigmar
! !INTERFACE:
!
      subroutine set_zigmar
!       
! !USES:    
      use input_module,only: potential_model
      use params_module,only:      
     |  nmlonp1,	! nmlon+1
     |  nmlat   	! number of geomagnetic grid latitudes
      use cons_module,only: 
     |  pi,             ! Pi
     |  dtr,            ! degrees to radians (pi/180)
     |  crit,        	! critical colatitudes crit(2) (radians)
     |  ylatm,   	! magnetic grid lats (radians)
     |  rtd
      use dynamo_module,only: zigm_r
      use aurora_module,only: 
     |  theta0          ! convection reversal boundary in radians
      use addfld_module,only: addfld
!
! !DESCRIPTION: 
! set conductances Sigma^R for prescribed potential at high latitudes
! add to rhs: - d /d phi_m * [Sigma^R /cos lam_m * d (Phi-Phi^R)/d phi_m]-
!     - d /d lam_m * [Sigma^R *cos lam_m * d (Phi-Phi^R)/d lam_m]
! Sigma^R= min(zig_bs*[exp(|lam_m-crit_lat|*2/3.)-1.] , zig_max) poleward of
!    the convection reversal boundary  otherwise 0.
!
      implicit none
!
! !RETURN VALUE:
! !PARAMETERS: 
!
! !REVISION HISTORY:
! 05.03.10  <Astrid Maute> <include header> 
! 
! EOP 
!      
      real, parameter :: zig_bs= 5.,   ! base value for Sigma^R = 5 S
     |      zig_max = 1.0e4            ! max value for Sigma^R = 1e4 S
!
      integer :: i,j,ihem
      real :: crit_lat(2),	! transition latitude for Sigma^R (south,north) 
     |        rad3inv		! 1/(3.*dtr) 
!
! mod_heelis is only working for symmetric potential pattern
! since only one hemisphere is solved for
      write(6,*) 'Modified Heelispotential: make sure pattern is ',
     |         'hemispherical symmetric'
      if(potential_model /= 'HEELIS') then
          write(6,*) 'subroutine set_zigmar: use hemispherical ',
     |         'symmetric Heelis pattern'
          write(6,*) 'subroutine set_zigmar: stop'
          stop 'subroutine set_zigmar'
      endif
! set latitude of convection reversal boundary (theta0(1/2) -> SH/NH)
! Note: theta0 set in aurora module is in co-latitude NH=SH
! 3 deg transition zone
!    
      crit_lat(2) = -theta0(2)+pi/2.-3.*dtr   
!    
      if( theta0(1) < pi/2 ) then  
         crit_lat(1) = -crit_lat(2)
      else 
         crit_lat(1) = -theta0(1)+pi/2.+3.*dtr    ! south: convert co-lat to lat
      endif
!
! set Sigma^R
!      
      rad3inv = (3.*dtr)**(-1)
      do j = 1,nmlat
        ihem = int(ylatm(j)*2./pi+2.)                  ! 1 - SH, 2 - NH
        if(abs(ylatm(j)) >= abs(crit_lat(ihem))) then
	   zigm_r(:,j) = zig_bs*(exp(abs(ylatm(j)-crit_lat(ihem))*2.*
     |         rad3inv)-1.)
	   zigm_r(:,j) = min(zigm_r(:,j),zig_max)
	else 
	   zigm_r(:,j) = 0.
        endif
      enddo	   
!
      end subroutine set_zigmar
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: add_zigmar
! !INTERFACE:
!
      subroutine add_zigmar
!     
! !DESCRIPTION: add  reference conductance to the dynamo conductances:
!    Sigma^T(0)_phiphi + Sigma^R
!    Sigma^T(0)_lamlam + Sigma^R
! modify Sigma^R to fit to Sigma^T(0)_phiphi and Sigma^T(0)_lamlam and add together
! no special treatment for equator or pole since at the equator Sigma^R are
! zero by definition and at the pole Sigma^R is constant
!         
! !USES:    
      use params_module,only: 
     |  nmlon,     ! number of geomagnetic grid longitudes     
     |  nmlonp1,   ! nmlon+1
     |  nmlat,      ! number of geomagnetic grid latitudes
     |  nmlath
      use cons_module,only: dt0dts,rcos0s,pi,dlonm,dlatm
      use dynamo_module,only: unitvm,
     |  zigm11,    !  +/+Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
     |  zigm22,    !  +/+Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
     |  zigm_r,    ! Sigma^R reference conductance
     |  rim_r      ! rhs^R for reference potential
      use addfld_module,only: addfld
!
      implicit none
!
! !RETURN VALUE:
! !PARAMETERS: 
! External:
      real,external :: sddot ! in util.F 
!
! !REVISION HISTORY:
! 05.03.10  <Astrid Maute> <include header> 
! 
! EOP 
!
      integer :: i,j,jj,jjj
      real :: corfac,cs(nmlat),dfac1n,dfac1s,dfac2n
!  
! real ::  zigm_r(nmlonp1,nmlat),! sigma_r
!     call addfld('ZIGM_R','zigm_r','S',zigm_r,
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
! Transforming PDE from original apex (lam_m) to new apex grid (lam_0)
!     lam_m is irregular spaced in latitude
!     lam_0 is regular spaced in latitude (used for derivatives)
! the whole PDE is divided by d lam_m /d lam_0
! DT1DTS : d lam_0/ d lam_m / |sin I_m|
! RCOS0S : cos(lam_0)/ cos(lam_m)
! corfac : |sin I_m|*d lam_m/d lam_0 * cos(lam_0)/ cos(lam_m)
!
! zigm11 = Sigma_(phi phi)(0)= Sigma_(phi phi)(m) * cos lam_0/cos lam_m * d lam_m/d lam_0
! zigm22 = Sigma_(lam lam)(0)= Sigma_(lam lam)(m) * cos lam_m/cos lam_0 * d lam_0/d lam_m
! zigm11 add: Sigma^R* cos lam_0/cos lam_m * d lam_m/d lam_0
! zigm22 add: Sigma^R* cos lam_m/cos lam_0 * d lam_0/d lam_m

      do j=1,nmlath-1                  ! SH 1:48
        jj = nmlat+1 - j               ! NH 97:50
        corfac = rcos0s(j)/dt0dts(j)
        do i=1,nmlon
	  rim_r(i,j,1) = zigm_r(i,j)*corfac  
	  rim_r(i,j,2) = zigm_r(i,j)/corfac
	  rim_r(i,jj,1)= zigm_r(i,jj)*corfac  
	  rim_r(i,jj,2)= zigm_r(i,jj)/corfac
        enddo ! i,nmlon
      enddo ! j=2,nmlat-1
!     
      rim_r(i,nmlath,:) = 0.
!
! Periodic points:
      do j=1,nmlat
        rim_r(nmlonp1,j,1)= rim_r(1,j,1)
        rim_r(nmlonp1,j,2)= rim_r(1,j,2)
      enddo ! j=1,nmlat 
      
! set magnetic latitude cosine array: cos lam_0
!
      do j = 1,nmlat     ! -pi/2 to pi/2
        cs(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo
!  
! Set up difference coefficients. 
! rim_r1 = Sigma^R(0)/ cos(lam_0)/d lon^2 )
! rim_r2 = Sigma^R(0)*cos(lam_0)/d lam_0^2
!
      do j = 2,nmlat-1    ! 2,96  not value at the poles 
	dfac1n = cs(j)/dlatm**2
	dfac2n = cs(j)*dlonm**2
        do i = 1,nmlonp1
          rim_r(i,j,2) = rim_r(i,j,2)*dfac1n	
          rim_r(i,j,1) = rim_r(i,j,1)/dfac2n	
        enddo
      enddo
!
! Values at the poles (1 and 97)
!
      jj  = nmlat
      jjj = 1
      dfac1n = cs(jj)/dlatm**2
      dfac1s = cs(jjj)/dlatm**2   ! is not necessary cos symmetric
      rim_r(:,jj,2)   =  rim_r(:,jj,2)*dfac1n
      rim_r(:,jjj,2)  =  rim_r(:,jjj,2)*dfac1s
!
! Set rim_r(1) to zero at the magnetic poles (1 and 97) to avoid floating 
!   point exception (values at the poles are not used)
! check this later
!
      rim_r(:,jj,1)  = 0.0
      rim_r(:,jjj,1) = 0.0
!  
! add to zigm11+rim_r1 and zigm22+rim_r2 
! Note it's not necessary for both hemisphere (only north is used afterwards)   
      do j=1,nmlath	! from s-pole to n-pole
        jj = nmlath+j-1   ! 49,97
        jjj = nmlath-j+1  ! 49,1
        zigm11(:,jj) = zigm11(:,jj)+ rim_r(:,jj,1)+rim_r(:,jjj,1)
        zigm22(:,jj) = zigm22(:,jj)+ rim_r(:,jj,2)+rim_r(:,jjj,2)
        zigm11(:,jjj)= zigm11(:,jj)
        zigm22(:,jjj)= zigm22(:,jj)
      enddo
      
! Compute polar values for the conductances, 4th order interpolation:
! 
      zigm11(1,    1) = (4.*sddot(nmlon,unitvm,zigm11(1,      2))-
     |  sddot(nmlon,unitvm,zigm11(1,      3)))/(3.*float(nmlon))
      zigm11(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm11(1,nmlat-1))-
     |  sddot(nmlon,unitvm,zigm11(1,nmlat-2)))/(3.*float(nmlon))
       zigm22(1,    1) = (4.*sddot(nmlon,unitvm,zigm22(1,      2))-
     |  sddot(nmlon,unitvm,zigm22(1,      3)))/(3.*float(nmlon))
      zigm22(1,nmlat) = (4.*sddot(nmlon,unitvm,zigm22(1,nmlat-1))-
     |  sddot(nmlon,unitvm,zigm22(1,nmlat-2)))/(3.*float(nmlon))
     
! Extend over longitude                                        
      do i = 2,nmlon
        zigm11(i,    1)  = zigm11(1,    1)
        zigm11(i,nmlat)  = zigm11(1,nmlat)
        zigm22(i,     1) = zigm22(1,    1)
        zigm22(i, nmlat) = zigm22(1,nmlat)
      enddo ! i = 2,nmlon
    
! rim_r vector (I_1,I_2): average over poles:
      do i = 1,nmlon
        rim_r(i,1,1) = .5*(rim_r(i,2,1)-rim_r(1+mod(i-1+nmlon/2,
     |     nmlon),2,1))
        rim_r(i,nmlat,1) = .5*(rim_r(i,nmlat-1,1)-
     |     rim_r(1+mod(i-1+nmlon/2,nmlon),nmlat-1,1))
        rim_r(i,1,2) = .5*(rim_r(i,2,2)-rim_r(1+mod(i-1+nmlon/2,
     |     nmlon),2,2))
        rim_r(i,nmlat,2) = .5*(rim_r(i,nmlat-1,2)-
     |     rim_r(1+mod(i-1+nmlon/2,nmlon),nmlat-1,2))
      enddo ! i = 1,nmlon
      
! Periodic points:
      do j=1,nmlat
        rim_r(nmlonp1,j,:)  = rim_r(1,j,:)
      enddo ! j=1,nmlat
      
! output to secondary history      
!     call addfld('RIM_R1','rim_r1',' ',rim_r(:,:,1),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('RIM_R2','rim_r2',' ',rim_r(:,:,2),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
      end subroutine add_zigmar
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: diff_rim
! !INTERFACE:
      subroutine diff_rimr(phihm)
!  
! !DESCRIPTION:
! calculate the current due to the reference potential
!   d /d phi_m * [Sigma^R /cos lam_m * d (Phi^R)/d phi_m] +
!   d /d lam_m * [Sigma^R *cos lam_m * d (Phi^R)/d lam_m]
!
! - coefficients calculated in subroutine add_zigmar
! - stencile due to d/d phi_m Sigma^R /cos lam_m d /d phi_m
!                   d/d lam_m * [Sigma^R *cos lam_m d/d lam_m
! - insert reference potential (Heelis) Phi^R
!
! Note: caculated for northern and southern hemisphere seperatly 
!         
! !USES:      
      use params_module,only: nmlonp1,nmlat,nmlath
      use cons_module,only: pi,dlatm
      use dynamo_module,only: nmlon0,nmlat0,
     |  rim_r,      ! coefficient for reference potential
     |  nsrim_r1,   ! difference stencil for rim_r(1)
     |  nsrim_r2,   ! difference stencil for rim_r(2)
     |  rhs_r       ! right-hand side 
      use addfld_module,only: addfld
!     
      implicit none
      
! !PARAMETERS: 
      real,intent(in) :: phihm  ! Heelis potential
!
! !REVISION HISTORY:
! 05.03.10  <Astrid Maute> <include header> 
! 
! EOP    
!
      integer :: i,j,jj,jjj,k
      real :: array(-15:nmlon0+16,nmlat0),sum,cs(nmlat)
!      
! set magnetic latitude cosine array: cos lam_0
!
      do j = 1,nmlat     ! -pi/2 to pi/2
        cs(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo
      
! Clear array for difference stencil over northern and southern hemisphere 
!
      nsrim_r1(:,:,:) = 0.0
      nsrim_r2(:,:,:) = 0.0
!  
! Calculate contribution to stencil for each coefficient
! rim_r1 = Sigma^RT(0)*cos(lam_0)/d lam_0^2
! rim_r2 = Sigma^RT(0)/cos(lam_0)/d lon^2 
! stencil stored in nsrim_r1 and nsrim_r2
!
      call nsstencil(rim_r(1,1,1),nmlon0,nmlat,nsrim_r1,array,1,nmlath)
      call nsstencil(rim_r(1,1,2),nmlon0,nmlat,nsrim_r2,array,4,nmlath)      
!
! Divide stencil by cos lam_0 to be conform with left-hand-side
!
      do j = 2,nmlat-1
	nsrim_r1(:,j,:) = nsrim_r1(:,j,:)/cs(j)
	nsrim_r2(:,j,:) = nsrim_r2(:,j,:)/cs(j)
      enddo     
!
! insert the Heelis-potential  
!
      rhs_r = 0.0
      call insert_pot(nsrim_r1,rhs_r,phihm)

!     call addfld('RHS_R1','rhs_r1',' ',rhs_r(:,:),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     
      call insert_pot(nsrim_r2,rhs_r,phihm)

!     call addfld('RHS_R1R2','rhs_r1r2',' ',rhs_r(:,:),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
      end subroutine diff_rimr
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: insert_pot
! !INTERFACE:
      subroutine insert_pot(nscoef,nscrrt,phim)
!     
! !DESCRIPTION: calculate the current nscrrt from the right hand side
!  by inserting the electric potential phim into the coefficient
!  stencil nscoef. Calculated for both hemispheres seperately.
!     
! !USES: 
      use params_module,only: nmlon,nmlonp1,nmlat,nmlath
      use dynamo_module,only: nmlon0,nmlat0
      implicit none
!     
! !PARAMETERS: 
      real,intent(in) :: 
     |   nscoef(nmlon0,nmlat,10),  ! coefficient stencil
     |   phim(nmlonp1,nmlat)	   ! electric potential
! !RETURN VALUE:
      real,intent(inout) ::
     |   nscrrt(nmlonp1,nmlat) 	   ! current
      real,parameter :: unitv(nmlon)=1.   
!
! !REVISION HISTORY:
! 05.03.11  <Astrid Maute> <include header> 
! 
! EOP 
! 
! External:
      real,external :: sddot ! in util.F
     
      integer :: i,j,jj,jjj
!    for i=1
!
      i = 1
      j = 49	 ! equator set J_mr = 0
      nscrrt(i,j) = 0.0	
!
      do j = 2,nmlath-1	  ! 2,48     no poles
        jj  = nmlath+j-1  ! 50,96
        jjj = nmlath-j+1  ! 48,2
!       			     contribution of stencil
! Northern hemisphere
!
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,1)*phim(i+1,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,2)*phim(i+1,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,3)*phim(i  ,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,4)*phim(nmlon0-1,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,5)*phim(nmlon0-1,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,6)*phim(nmlon0-1,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,7)*phim(i  ,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,8)*phim(i+1,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,9)*phim(i  ,jj)
!        nscrrt(i,jj) = nscrrt(i,jj) - nscoef(i,jj,10)
!           
! Southern hemisphere
!
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,1)*phim(i+1,jjj)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,2)*phim(i+1,jjj-1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,3)*phim(i  ,jjj-1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,4)*phim(nmlon0-1,jjj-1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,5)*phim(nmlon0-1,jjj)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,6)*phim(nmlon0-1,jjj+1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,7)*phim(i  ,jjj+1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,8)*phim(i+1,jjj+1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,9)*phim(i  ,jjj)
!	nscrrt(i,jjj)=nscrrt(i,jjj)-nscoef(i,jjj,10)
      enddo
!
      do i = 2,nmlon0-1         ! 2,80
        j = 49     ! equator
        nscrrt(i,j) = 0.0
        do j = 2,nmlath-1   ! 2,48     no poles  
          jj  = nmlath+j-1  ! 50,96
          jjj = nmlath-j+1  ! 48,2
!				       contribution of stencil
! Northern hemisphere
!
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,1)*phim(i+1,jj)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,2)*phim(i+1,jj+1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,3)*phim(i  ,jj+1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,4)*phim(i-1,jj+1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,5)*phim(i-1,jj)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,6)*phim(i-1,jj-1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,7)*phim(i  ,jj-1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,8)*phim(i+1,jj-1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,9)*phim(i  ,jj)
!	 nscrrt(i,jj) = nscrrt(i,jj) - nscoef(i,jj,10)
!
! Southern hemisphere
!
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,1)*phim(i+1,jjj)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,2)*phim(i+1,jjj-1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,3)*phim(i  ,jjj-1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,4)*phim(i-1,jjj-1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,5)*phim(i-1,jjj)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,6)*phim(i-1,jjj+1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,7)*phim(i  ,jjj+1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,8)*phim(i+1,jjj+1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,9)*phim(i  ,jjj)
!	 nscrrt(i,jjj) = nscrrt(i,jjj) - nscoef(i,jjj,10)
	enddo
      enddo
!
! For i=nmlonp1 (81)
!
      i = nmlonp1
      j = nmlath	 ! equator (49)
      nscrrt(i,j) = 0.0
      do j = 2,nmlath-1	  ! 2,48      no poles 
        jj  = nmlath+j-1  ! 50,96
        jjj = nmlath-j+1  ! 48,2
!       			     contribution of stencil
! Northern hemisphere
!
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,1)*phim(2,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,2)*phim(2,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,3)*phim(i  ,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,4)*phim(i-1,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,5)*phim(i-1,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,6)*phim(i-1,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,7)*phim(i  ,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,8)*phim(2,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,9)*phim(i  ,jj)
!	nscrrt(i,jj) = nscrrt(i,jj) - nscoef(i,jj,10)
! 
! Southern hemisphere
!
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,1)*phim(2,jjj)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,2)*phim(2,jjj-1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,3)*phim(i  ,jjj-1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,4)*phim(i-1,jjj-1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,5)*phim(i-1,jjj)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,6)*phim(i-1,jjj+1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,7)*phim(i  ,jjj+1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,8)*phim(2,jjj+1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,9)*phim(i  ,jjj)
!        nscrrt(i,jjj) = nscrrt(i,jjj) -nscoef(i,jjj,10)
      enddo
!
! Poles
!      
      nscrrt(1,1)    = sddot(nmlon,unitv,nscrrt(1,2))/float(nmlon)
      nscrrt(1,nmlat)= sddot(nmlon,unitv,nscrrt(1,nmlat-1))/float(nmlon)
      nscrrt(:,1)     = nscrrt(1,1)		! extend in longitude
      nscrrt(:,nmlat) = nscrrt(1,nmlat)
!      
      end subroutine insert_pot
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: add_rimr
! !INTERFACE:
      subroutine add_rimr
!  
! !DESCRIPTION:
! add current due to reference potential both hemispheres together and 
! add currentto RHS
!      
! !USES:    
      use params_module,only: nmlonp1,nmlath,nmlatp1
      use dynamo_module,only:  
     |  rhs,	! right-hand side 
     |  rhs_r	! right hand side contribution from zigm_r
!     
      implicit none 
! !RETURN VALUE:
! !PARAMETERS:   
!
! !REVISION HISTORY:
! 05.03.11  <Astrid Maute> <include header> 
! 
! EOP 
!
      integer :: i,j
!      
! add both hemispheres together at equator double value- shouldn't
! matter since Heelis-potential should be zero
! rhs(1:nmlath) equator to N-pole
! rhs_r(1:nmlat) S-pole to N-pole
!
      do j=1,nmlath	! from s-pole to n-pole
        rhs_r(:,nmlatp1-j) = (rhs_r(:,nmlatp1-j)+rhs_r(:,j)) 
	rhs_r(:,j)         = rhs_r(:,nmlatp1-j)
        rhs(:,nmlath+1-j)  = rhs(:,nmlath+1-j) + rhs_r(:,j) 
      enddo ! j=1,nmlath
!
      end subroutine add_rimr
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: calrhs_jrr
! !INTERFACE:
      subroutine calrhs_jrr(phihm)
!      
! !DESCRIPTION:
! J_rR contribution to dynamo - field-aligned current part representing
! only region-1 current, therefore fac = 0 below the convection reversal 
! boundary
! calculation of J_rR is set up such that the total current system, horizontal
! and radial adds up to zero
!
! J_rR = fac*{d / d phi_m [ Sigma_phiphi/ cos(lam_m)* d /d |lam_m| Phi_H +
!        Sigma_philam   d /d phi_m Phi_H ] + d / d |lam_m| [ Sigma_lamphi*
!        d /d phi_m Phi_H + Sigma_lamlam*cos(lam_m)*d /d |lam_m| Phi_H ] }
! 
!
! NOTE: - this part wasn't tested
!       - conductances added together - does not matter for
!         derivatives - only factor of two
!       - not necessary to calculate J_rR for both hemispheres
!         since Phi_H / Hellis potential is symmetric
!         could be changed, insert_pot used by other subroutines too.
!   
! !USES:     
      use params_module,only: nmlon,nmlonp1,nmlatp1,nmlat,nmlath
      use cons_module,only: dlatm,pi,crit,ylatm
      use dynamo_module,only: nmlon0,nmlat0,
                   ! zigms from south to north pole (1:97)
     |  zigm11,    ! +/+Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
     |  zigm22,    ! +/+Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
     |  zigmc,     ! +/-Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
     |  zigm2,     ! +/-Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
     |  rhs        ! right-hand side (1:49) equator to north pole
      use addfld_module,only: addfld
!     
      implicit none
! !RETURN VALUE:
! !PARAMETERS: 
      real,intent(in) :: phihm(nmlonp1,nmlat)  ! Heelis potential     
!
! !REVISION HISTORY:
! 05.03.09  <Astrid Maute> <include header> 
! 
! EOP 
! 
! local   
      integer :: i,j  
      real :: nscoef(nmlon0,nmlat,10),
     |   array(-15:nmlon0+16,nmlat0),
     |   rhs_jrr(nmlonp1,nmlat),cs(nmlat),crit_lat1,f(nmlath)
      real,parameter :: unitv(nmlon)=1.
!
! External:
      real,external :: sddot ! in util.F
!
      rhs_jrr = 0.
!
! Calculate cos lam_0
!
      do j = 1,nmlat    
        cs(j) = cos(-pi/2.+(j-1)*dlatm) ! -pi/2 to pi/2
      enddo
!      
! calculate factor (convection reversal boundar at 75 deg)
! lam_m < lam(convection reversal boundary)  -> fac = 0
! lam_m >= lam(convection reversal boundary) -> fac = 1
!
      crit_lat1 = -crit(1) + pi/2.	
      do j = 1,nmlath
        if(abs(ylatm(j)) >= crit_lat1) then
	   f(j) = 1.
	else 
	   f(j) = 0.
        endif
      enddo	               
!
! Calculate contribution to stencil from each coefficient
! nsstencil: set up stencil (finite differencing)
! divide:    stencil by cos lam_0 note: this has to be done before
! 	the potential is inserted to be conform with the left-hand-side
! insert_pot:insert the Heelis-potential
!
! NOTE: calculated for both hemispheres which is not necessary
!
! "NH" array index 49:97 equator to pole
! zigmc  = Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm2  = Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm22 = Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
! zigm11 = Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
!
! "SH" array index 49:1 equator to pole
! zigmc  = -Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm2  = -Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
! zigm22 = Sigma_(lam lam)^T(0)*cos(lam_0)/d lam_0^2
! zigm11 = Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
!
      nscoef = 0.
      call nsstencil(zigm11,nmlon0,nmlat,nscoef,array,1,nmlath)
      do j = 2,nmlat-1
	nscoef(:,j,:) = nscoef(:,j,:)/cs(j)
      enddo     
      call insert_pot(nscoef,rhs_jrr,phihm)
!      
      nscoef = 0.
      call nsstencil(zigm22,nmlon0,nmlat,nscoef,array,4,nmlath)
      do j = 2,nmlat-1
	nscoef(:,j,:) = nscoef(:,j,:)/cs(j)
      enddo     
      call insert_pot(nscoef,rhs_jrr,phihm)
!      
      nscoef = 0.
      call nsstencil(zigmc,nmlon0,nmlat,nscoef,array,2,nmlath)
      do j = 2,nmlat-1
	nscoef(:,j,:) = nscoef(:,j,:)/cs(j)
      enddo     
      call insert_pot(nscoef,rhs_jrr,phihm)
!      
      nscoef = 0.
      call nsstencil(zigm2,nmlon0,nmlat,nscoef,array,3,nmlath)
      do j = 2,nmlat-1
	nscoef(:,j,:) = nscoef(:,j,:)/cs(j)
      enddo     
      call insert_pot(nscoef,rhs_jrr,phihm)
!      
! add both hemispheres of J_rR (rhs_jrr) together, multiply by factor f
!  J_rR (rhs_jrr)(1:97)  south pole to north pole
! and add to rhs (1:49) equator to pole
!
      do j=2,nmlath	
        rhs_jrr(:,j) = rhs_jrr(:,nmlatp1-j)+rhs_jrr(:,j)
	rhs_jrr(:,j) = f(j)*rhs_jrr(:,j)		
	rhs_jrr(:,nmlatp1-j) = rhs_jrr(:,j)
        rhs(:,nmlath+1-j) = rhs(:,nmlath+1-j) + rhs_jrr(:,j) 
      enddo ! j=1,nmlath
!
! pole values
!     
      rhs_jrr(1,1)    = sddot(nmlon,unitv,rhs_jrr(1,2))/float(nmlon)
      rhs_jrr(1,nmlat)= sddot(nmlon,unitv,rhs_jrr(1,nmlat-1))/
     |   float(nmlon)
      rhs_jrr(:,1)    = rhs_jrr(1,1)			   ! extend in longitude
      rhs_jrr(:,nmlat)= rhs_jrr(1,nmlat)
      rhs_jrr(:,1)    = rhs_jrr(:,nmlatp1-1)+rhs_jrr(:,1)  ! add hemispheres
      rhs(:,nmlath) = rhs(:,nmlath) + rhs_jrr(:,1)         ! add to rhs
      
! output secondary history 2d mag     
!     call addfld('RIM_JRR','rim_jrr',' ',rhs_jrr,
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
      end subroutine calrhs_jrr
!------------------------------------------------------------------------------
! BOP
! !IROUTINE: set_cicr
! !INTERFACE:
      subroutine set_cicr
!      
! !DESCRIPTION:
! Set up equivalent magnetospheric conductances 
! C_r/ Sigma_phiphi^M: equivalent magnetospheric zonal Pedersen conductance
! C_i/ Sigma_H^M : equivalent magnetospheric Hall conductance
! for the TIEGCM grid
! these conductances should generate shielding by region 2 current 
! 
!  Values for equivalent magnet. conductances in 1 deg increasing colatitude
! from [ Peymirat, Richmond 1993: JGR, Vol98, pp.15467, fig. 4 ]   
!        
! !USES:     
      use params_module,only:      
     |  nmlonp1,	! nmlon+1
     |  nmlat,   	! number of geomagnetic grid latitudes
     |  nmlath          ! index to magnetic equator
      use cons_module,only: 
     |  pi,             ! Pi
     |  rtd,            ! radians-to-degrees (180./pi)
     |  dtr,            ! degrees-to-radians (pi/180.)
     |  crit,        	! critical colatitudes crit(2) (radians)
     |  ylatm   	! magnetic grid lats (radians)
      use dynamo_module,only: c_r,c_i,zigmc_r,zigmc_i,
     |  eqMg_crit       ! colatitude where C_r,C_i start
      use aurora_module,only: 
     |  theta0          ! convection reversal boundary in radians
!
      implicit none
!
! !RETURN VALUE:
! !PARAMETERS: 
!
!
! !REVISION HISTORY:
! 05.03.10  <Astrid Maute> <include header> 
! 
! EOP 
      integer :: i,j,i_old
      real :: crit_lat1,lat_crci(15),y,fac
      
! set up latitudes for C_r, C_i conductances
! start at 18 deg colatitude convert to latitude
!      
!      crit_lat1   = -eqMg_crit + pi/2.	! starts at 18 deg colatitude
      crit_lat1   = -theta0(2) + pi/2.	! starts at convection reversal boundary
      lat_crci(1) = crit_lat1
      fac = 1.*dtr			! 1 degree in radians
!      
      do i = 2,15
        lat_crci(i) = lat_crci(i-1) - fac ! increments by -1 degree (fac)
      enddo
!      
! interpolate C_r/ Sigma_phiphi^M [S], C_i/ Sigma_H^M [S] conductances to TIEGCM grid
! only one hemisphere is done and copy to other HS
! change order of interpolated values to equator - pole

! loop over half hemisphere (N-pole to equator)
      fac = 1./(lat_crci(2) -lat_crci(1))   ! factor for interpolation
      i_old = 1 			    ! index to begin search
 
      do j = nmlat,nmlath,-1
        if(abs(ylatm(j)) > lat_crci(1)) then	  ! polar cap
	    zigmc_r(:,j) = 0.	        ! N-hemisphere
	    zigmc_i(:,j) = c_i(1)      
	    zigmc_r(:,nmlat-j+1) = 0.   ! S-hemisphere
	    zigmc_i(:,nmlat-j+1) = c_i(1) 
	elseif(abs(ylatm(j)) < lat_crci(15)) then  ! towards equator
	    zigmc_r(:,j) = 0.           ! N-hemisphere
	    zigmc_i(:,j) = 0. 
	    zigmc_r(:,nmlat-j+1) = 0.   ! S-hemisphere
	    zigmc_i(:,nmlat-j+1) = 0. 
	else	! regime of C_r, C_i
    	   do i = i_old,15	! find interval lat_crci(i) < ylatm(j) < lat_crci(i+1)
	     if(lat_crci(i) >= ylatm(j).and.
     |           ylatm(j) > lat_crci(i+1)) then
                y = c_r(i) + (c_r(i+1)-c_r(i))*
     |                (ylatm(j)-lat_crci(i))*fac
	        zigmc_r(:,j)         = y   ! N-hemisphere
	        zigmc_r(:,nmlat-j+1) = y   ! S-hemisphere
!		
                y = c_i(i) + (c_i(i+1)-c_i(i))*
     |                (ylatm(j)-lat_crci(i))*fac
	        zigmc_i(:,j)         = y   ! N-hemisphere
	        zigmc_i(:,nmlat-j+1) = y   ! S-hemisphere
		
		i_old = i
		goto 10
	     endif
	   enddo   ! end loop over lat_crci(i)
!	   
	endif
  10    continue ! next ylatm(j)
      enddo   ! end loop over ylatm(j) 
!
      end subroutine set_cicr
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: add_cicr
! !INTERFACE:
      subroutine add_cicr
! 
! !DESCRIPTION: add  equivalent magnetospheric conductance
! to the dynamo conductances:
!    Sigma^T_phiphi + Sigma_phiphi^M
!    Sigma^T_philam + Sigma_H^M
!    Sigma^T_lamphi - Sigma_H^M
! Sigma_phiphi^M: equivalent magnetospheric zonal Pedersen conductance
! Sigma_H^M: equivalent magnetospheric Hall conductance
! - acts on below the convection reversal boundary
!
! modify Sigma^M to fit to zigm11, zigm12 and zigm21 from the dynamo 
! and add together
! no special treatment for equator or pole since at the equator zigm_r should
! be zero by definition and constant at the pole
!           
! !USES: 
      use params_module,only: 
     |  nmlon,     ! number of geomagnetic grid longitudes     
     |  nmlonp1,   ! nmlon+1
     |  nmlat,      ! number of geomagnetic grid latitudes
     |  nmlath
      use cons_module,only: dt0dts,rcos0s,pi,dlonm,dlatm
      use dynamo_module,only: 
     |  zigm11,    ! Sigma_(phi phi)^T(0)/ cos(lam_0) / d lon^2 )
     |  zigmc,     ! Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
     |  zigm2,     !Sigma_(phi lam)^T(0)/( 4*d lam_0* d lon )
     |  zigmc_r,   ! C_r on tiegcm grid
     |  zigmc_i    ! C_i on tiegcm grid 
      use addfld_module,only: addfld
      implicit none
!
! !RETURN VALUE:
! !PARAMETERS: 
!
! !REVISION HISTORY:
! 05.03.10  <Astrid Maute> <include header> 
! 
! EOP 
!
      integer :: i,j,jj,jjj
      real :: corfac,cs(nmlat),dfac1n,dfac1s,dfac2n,fac,
     |   zig11(nmlonp1,nmlat),zig12(nmlonp1,nmlat)
!  
!     call addfld('ZIGMC_R','zigmc_r','S',zigmc_r,
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('ZIGMC_I','zigmc_i','S',zigmc_i,
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
! Transforming PDE from original apex (lam_m) to new apex grid (lam_0)
!     lam_m is irregular spaced in latitude
!     lam_0 is regular spaced in latitude (used for derivatives)
! the whole PDE is divided by d lam_m /d lam_0
! DT1DTS : d lam_0/ d lam_m / |sin I_m|
! RCOS0S : cos(lam_0)/ cos(lam_m)
! corfac : |sin I_m|*d lam_m/d lam_0 * cos(lam_0)/ cos(lam_m)
! 
! zig11    = Sigma_phiphi^M|*d lam_m/d lam_0 * cos(lam_0)/ cos(lam_m)
! zig12    = Sigma_(H)^M(0)

      do j=1,nmlath-1                  ! SH 1:48
        jj = nmlat+1 - j               ! NH 97:50
        corfac = rcos0s(j)/dt0dts(j)   ! sym. about equator
        do i=1,nmlon
	  zig11(i,j) = zigmc_r(i,j)*corfac  
	  zig12(i,j) = zigmc_i(i,j)
	  zig11(i,jj)= zigmc_r(i,jj)*corfac  
	  zig12(i,jj)= zigmc_i(i,jj)
        enddo ! i,nmlon
      enddo ! j=2,nmlat-1
!      
      zig11(:,nmlath) = 0.
      zig12(:,nmlath) = 0.
!
! Periodic points:
      do j=1,nmlat
        zig11(nmlonp1,j)= zig11(1,j)
        zig12(nmlonp1,j)= zig12(1,j)
      enddo ! j=1,nmlat 
      
! set magnetic latitude cosine array: cos lam_0
!
      do j = 1,nmlat     ! -pi/2 to pi/2
        cs(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo
!   
! Set up difference coefficients. 
! zig11 = Sigma_(phi phi)^M(0)/ cos(lam_0) / d lon^2 )
! zig12 = zig21 = Sigma_(H)^M(0)/( 4*d lam_0* d lon )
!
      fac = 1./(4.*dlatm*dlonm)
      do j = 2,nmlat-1    ! 2,96  not value at the poles 
	dfac1n = cs(j)*dlonm**2
        do i = 1,nmlonp1
          zig11(i,j) = zig11(i,j)/dfac1n	
          zig12(i,j) = zig12(i,j)*fac
        enddo
      enddo
!
! Values at the poles (1 and 97)
!
      zig12(:,1)     = zig12(:,1)*fac
      zig12(:,nmlat) = zig12(:,nmlat)*fac
!
! Set zig11 to zero at the magnetic poles (1 and 97) to avoid floating 
!   point exception (values at the poles are not used)
! check this later
!
      zig11(:,1)     = 0.0
      zig11(:,nmlat) = 0.0
!  
! add to zigm11+zig11, 
!        zigmc+zig12 and zigm2-zig12 
! Note it's not necessary for both hemisphere (only north is used afterwards)
!   
      do j=1,nmlath    ! from s-pole to n-pole
	jj  = nmlath+j-1  ! 49,97
	jjj = nmlath-j+1  ! 49,1
        zigm11(:,jj) = zigm11(:,jj)+ zig11(:,jj)+zig11(:,jjj)  ! NH
        zigmc(:,jj)  = zigmc(:,jj) + zig12(:,jj)+zig12(:,jjj)
        zigm2(:,jj)  = zigm2(:,jj) - zig12(:,jj)-zig12(:,jjj)
!      
        zigm11(:,jjj)= zigm11(:,jj)			      ! SH
        zigmc(:,jjj) = zigmc(:,jjj) - zig12(:,jj)-zig12(:,jjj)
        zigm2(:,jjj) = zigm2(:,jjj) + zig12(:,jj)+zig12(:,jjj)
      enddo
      
! output to secondary history      
!     call addfld('RIMC_11M','rimc_11M',' ',zig11(:,:),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     call addfld('RIMC_12M','rimc_12M',' ',zig12(:,:),
!    |  'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
      end subroutine add_cicr
!----------------------------------------------------------------------------

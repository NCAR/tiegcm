      module magpres_g
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate the contribution from magnetic pressure and gravity to the 
! RHS of the dynamo equation current (mag.pressure, gravity) calculated 
! at half levels
!
      use fields_module,only: levd0,levd1,lond0,lond1,latd0,latd1 ! subdomain dimensions
      use params_module,only: nlevp1,  
     |  nmlonp1,! nmlon+1  
     |  nlon,   ! nmlon
     |  nlonp4, ! nmlon+4
     |  nlat,   ! number of latitudes
     |  nmlat   ! number of geomagnetic grid latitudes
!      
      implicit none
      real, parameter :: matm = 1.6605E-24 ! mass per atomic weight [g]
!
! Module output: 
      real,dimension(:,:,:),pointer,save :: ! (nlevp1,lon0:lon1,lat0:lat1)
     |  je1oD_geo,  ! J_e1(pg)/D [A/cm^2] 
     |  je2oD_geo   ! J_e2(pg)/D [A/cm^2]
!       
      real,dimension(nmlonp1,nmlat) :: tpg1,tpg2 ! only for testing
!
      contains
!-----------------------------------------------------------------------
      subroutine magpres_grav(z,te,ti,ne,op,nplus, 
     |      n2p,nop,o2p,lev0,lev1,lon0,lon1,lat0,lat1)
!
! Calculate current due to plasma pressure gradient and 
!   gravity which is added to the right hand side of the dynamo
!   equation. The gravity driven current changes the electric field
!   before sunrise and after sunset. The influence of the plasma pressure 
!   gradient driven current on the electric field can be neglected, however
!   if the magnetic perturbations are calculated it will have an influence.
!      
      use cons_module,only:
     |  grav,       ! 870.cm/s^2 accel due to gravity (dep on l.bndry)
     |  boltz,      ! boltzman's constant  cm^2*g/s^2/K    
     |  rmass_o2 , rmass_o1, rmass_n2,  ! 32,16,28
     |  rmass_n4s, rmass_no, rmass_op,   ! 14,30,16
     |  ylatg,  	! geographic latitudes (radians)
     |  ylong,    	! geographic longitudes (radians)
     |  dlatg,dlong,    ! pi/float(nlat) 2.*pi/float(nlon)
     |  re,             ! earth radius (cm)  6.37122e8  
     |  r0,             ! r_e+h0
     |  h0              ! h0 = 90 km b, dvec scaled in apex to 90 km
      use magfield_module,only: 
     |  xb,      ! northward component of magnetic field [Gauss = 10^-4 T]
     |  yb,      ! eastward component of magnetic field [Gauss = 10^-4 T]
     |  zb,      ! downward component of magnetic field [Gauss = 10^-4 T]
     |  bmod,    ! magnitude of magnetic field [Gauss = 10^-4 T]
     |  alatm    ! (nlonp1,0:nlatp1)  
      use getapex_module,only:           ! (see sub apxparm, getapex.F)
     |  dvec,     ! d_1,d_2 referenced to h_0 (nlonp1,nlat,3,2)
     |  dddarr    ! D referenced to h_0 (nlonp1,nlat) 
      use addfld_module,only: addfld
#ifdef MPI
      use mpi_module,only: mp_periodic_f3d
#else
      use params_module,only: nlon,nlonp2
#endif     
      implicit none
!
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Input args: press vs longitude input fields (2d (k,i)):
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     | intent(in) ::
     |    z,    ! geopotential input full level
     |    te, 	! electron temperature K
     |    ti, 	! ion temperature K
     |    ne, 	! electron density 1/cm^3 full level
     |    op, 	! O+ 1/cm^3
     |    nplus,! N+  1/cm^3 
     |    n2p,  ! N2+ 1/cm^3 
     |    nop,  ! NO+  1/cm^3
     |    o2p   ! O2+  1/cm^3
!
! Local:
      integer :: i0,i1,nk,nkm1,i,k,lat,lonbeg,lonend
      real :: rho,sumt,gne,fac,sumtm,sumtp,halfne,halfz,
     |  tk,tkp,gt,fac_horlon,fac_horlat,fac_hgt,oB2,
     |  horlon_ne,horlon_t,horlat_ne,horlat_t,
     |  sinalat,clm2,r0or,hgt1,hgt2,
     |  jmagp_horlon,jmagp_horlat,jmagp_vert,jgrav
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2) ::
     |  j_glon,j_glat, j_gz, plt_gt, plt_gne, plt_ne_horlon,
     |  plt_t_horlon, plt_ne_horlat, plt_t_horlat, plt_grav, plt_rho
      real :: ftmp(lev0:lev1,lon0:lon1,lat0:lat1,2)
!
! set longitude index boundary
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = lon1-2
!      
! Convenience ints for addfsech calls:
      i0 = lon0
      i1 = lon1
      nk = lev1-lev0+1
      nkm1 = nk-1
!      
      fac = 10.*grav*matm ! g*mass per atomic weight factor gauss,kg -> T,g

      j_glon=0. ; j_glat=0. ; j_gz=0.
      plt_gt=0. ; plt_gne=0. ; plt_ne_horlon=0. ; plt_t_horlon=0.
      plt_ne_horlat=0. ; plt_t_horlat=0. ; plt_grav=0. ; plt_rho=0.
!
! Latitude loop: 
      do lat=lat0,lat1
	fac_horlon = 2.*dlong*cos(ylatg(lat))  ! for d()/dphi
	fac_horlon = 1./fac_horlon
	fac_horlat = 2.*dlatg
!	
! Longitude loop: 
        do i= lonbeg,lonend
          if(i <= nlon+1) then  		 ! sin(lam)
	     sinalat	= sin(alatm(i,lat))	   
	  else
	     sinalat	= sin(alatm(i-nlon,lat))	
	  endif
          clm2       = 1. - sinalat*sinalat  ! cos^2(lam) 
!	  
! Levels loop:
          do k=lev0,lev1-1 
	     halfz = 0.5*(z(k,i,lat) + z(k+1,i,lat))    ! Z to half level - k'
	     fac_hgt = (re+halfz)/(re+h0)               ! r/(r_e+h_0) to take weaker 
! ion density                                           ! field g,B into account - k'
             rho = op(k,i,lat)*rmass_o1 + o2p(k,i,lat)*rmass_o2+  ! k'
     |          nplus(k,i,lat)*rmass_n4s+ n2p(k,i,lat)*rmass_n2+
     |          nop(k,i,lat)*rmass_no 
! gravity contribution
             jgrav = -rho*fac/fac_hgt**2  	! k' gravity height variation [(r_e+h_0)/r]^2
!	     		
! plasma pressure
             sumt = te(k,i,lat)  +ti(k,i,lat)           ! T half level - k'
	     if(k == lev0 ) then                        ! T at bottom -1 not defined
	       sumtm= 2.*te(k,i,lat)-te(k+1,i,lat) +   ! extrapolte like in postprocessor
     |                2.*ti(k,i,lat)-ti(k+1,i,lat)
	     else
               sumtm= te(k-1,i,lat)+ti(k-1,i,lat)      ! T half level (k-1)' 
	     endif
	     if(k == lev1-1) then                       ! T at top level not defined
	        sumtp = 2.*sumt - sumtm                ! extrapolte like in postprocessor
	     else
                sumtp= te(k+1,i,lat)+ti(k+1,i,lat)         ! T half level (k+1)' 
	     endif
	     				
	     halfne= 0.5*(ne(k,i,lat) + ne(k+1,i,lat))  ! Ne to half level - k'
! vertical gradient 
! Ne = [Ne(k)-Ne(k+1)]/dz  positive increasing downward to fit to gravity term
             gne = ne(k,i,lat) - ne(k+1,i,lat)            ! ne full level - k'
	     gne = gne/(z(k+1,i,lat) - z(k,i,lat))        ! z full level  - k'
! T = [T(k)-T(k+1)]/dz  positive increasing downward to fit to gravity term
             tk  =  0.5*(sumt+sumtm)	! full level      ! k
             tkp =  0.5*(sumtp+sumt)	! full level      ! k+1
             gt  =  (tk-tkp)/(z(k+1,i,lat) - z(k,i,lat))  ! k'
	     jmagp_vert = gne*sumt+gt*halfne              ! k'
	     jmagp_vert = 10.*boltz*jmagp_vert ! k'
!	     
! horizontal gradient - longitudinal direction
             horlon_ne =  0.5*(ne(k,i+1,lat) + ne(k+1,i+1,lat)) -   ! dne/dphi half levels - k'
     |                    0.5*(ne(k,i-1,lat) + ne(k+1,i-1,lat))	
             horlon_t  =  te(k,i+1,lat)  +ti(k,i+1,lat)  -          ! dt/dphi half level - k'
     |                    (te(k,i-1,lat) +ti(k,i-1,lat))  
             horlon_ne =  horlon_ne*sumt 	! (Te+Ti)*dne/dphi  ! k'	
             horlon_t  =  horlon_t*halfne	! Ne*d(Te+Ti)/dphi  ! k'
	     jmagp_horlon = (horlon_ne+horlon_t)*fac_horlon/        ! k'
     |                (re+halfz)
	     jmagp_horlon = -10.*boltz*jmagp_horlon      ! k'
!	     
! horizontal gradient - latitudinal direction
             if(lat /= 1 .or. lat /= nlat) then
               horlat_ne =  0.5*(ne(k,i,lat+1) + ne(k+1,i,lat+1)) -  ! dne/dlat half levels - k'
     |          	    0.5*(ne(k,i,lat-1) + ne(k+1,i,lat+1)) 
               horlat_t  =  te(k,i,lat+1)  +ti(k,i,lat+1)  -	     ! dt/dlat half levelk'
     |          	   (te(k,i,lat-1)  +ti(k,i,lat-1))  
               horlat_ne =  horlat_ne*sumt	  ! (Te+Ti)*dne/dphi ! k'	  
               horlat_t  =  horlat_t*halfne	  ! Ne*d(Te+Ti)/dphi ! k'
	       jmagp_horlat = (horlat_ne+horlat_t)/(fac_horlat*      ! k'
     |          	(re+halfz))
               jmagp_horlat = -10.*boltz*jmagp_horlat     ! k'
	     else	! set value to zero at the poles other gradients are calculated
	       horlat_ne = 0.
	       horlat_t  = 0.
	       jmagp_horlat = 0.
             endif
!	     
! ()x B / B_mag^2      magnetic field in [Gauss] -> j_* in [A/cm^2]
! xb north, yb east, zb down -> (yb,xb,-zb)
             oB2 = fac_hgt**3/bmod(i-2,lat)**2    ! height variation of B [(r_e+h_0)/r]^3
             j_glon(k,i,lat) = -jmagp_horlat*zb(i-2,lat) - 
     |          (jmagp_vert+jgrav)*xb(i-2,lat)
             j_glat(k,i,lat) = (jmagp_vert+jgrav)*yb(i-2,lat)+ 
     |          jmagp_horlon*zb(i-2,lat)
             j_gz(k,i,lat)   = jmagp_horlon*xb(i-2,lat) - 
     |          jmagp_horlat*yb(i-2,lat)
             j_glon(k,i,lat) = j_glon(k,i,lat)*oB2  ! 1/B^2
             j_glat(k,i,lat) = j_glat(k,i,lat)*oB2
             j_gz(k,i,lat)   = j_gz(k,i,lat)*oB2
!
! d_1/D dot J(pg) -> J_e1_pg
! d_2/D dot J(pg) -> J_e2_pg
!
! am 6/04: d_1,D_2,D are referenced to h_0 (90 km) s. apex.F
! in the following they are scaled from h_0 to the actual altitude
! half levels - k'
! hgt1: d_1(90km) -> [(R_e+h_0)/(R_e+z(half)))]^1.5
! hgt2: d_2(90km) -> [(R_e+h_0)/(R_e+z(half)))]^1.5* sqrt[(4.-3.*cos^2(lam))/(4.-3.*r0or*clm2)]
! je1oD_geo : d_1/D -> d_1/[d_1*d_2] = 1/d_2
! je2oD_geo : d_2/D -> d_2/[d_1*d_2] = 1/d_1
!
             r0or = r0/(r0 + halfz - h0)  ! R_0/R  
             hgt1 = r0or**1.5             ! for d_1
             hgt2 = hgt1*sqrt((4.-3.*clm2)/(4.-3.*r0or*clm2)) ! for d_2	

            if(i <= nlon+1) then	     
! 1/D *d_1 dot J(pg) = [ d_1(1)*j_glon  + d_1(2)*j_glat  + d_1(3)*j_gz ]/hgt2/D 
             je1oD_geo(k,i,lat) = (                              
     |        dvec(i,lat,1,1)*j_glon(k,i,lat) +  		       
     |        dvec(i,lat,2,1)*j_glat(k,i,lat) +
     |        dvec(i,lat,3,1)*j_gz(k,i,lat))/hgt2/dddarr(i,lat)
	     
! 1/D *d_2 dot J(pg) = [ d_2(1)*j_glon  + d_2(2)*j_glat  + d_2(3)*j_gz ]/hgt1/D
             je2oD_geo(k,i,lat) = (                              
     |        dvec(i,lat,1,2)*j_glon(k,i,lat) +		      
     |        dvec(i,lat,2,2)*j_glat(k,i,lat) +
     |        dvec(i,lat,3,2)*j_gz(k,i,lat))/hgt1/dddarr(i,lat) 
            else	     
! 1/D *d_1 dot J(pg) = [ d_1(1)*j_glon  + d_1(2)*j_glat  + d_1(3)*j_gz ]/hgt2/D
             je1oD_geo(k,i,lat) = (                              
     |        dvec(i-nlon,lat,1,1)*j_glon(k,i,lat) +		      
     |        dvec(i-nlon,lat,2,1)*j_glat(k,i,lat) +
     |        dvec(i-nlon,lat,3,1)*j_gz(k,i,lat))/hgt2/
     |        dddarr(i-nlon,lat)	  
	     
! 1/D *d_2 dot J(pg) = [ d_2(1)*j_glon  + d_2(2)*j_glat  + d_2(3)*j_gz ]/hgt1/D
             je2oD_geo(k,i,lat) = (                              
     |        dvec(i-nlon,lat,1,2)*j_glon(k,i,lat) +		     
     |        dvec(i-nlon,lat,2,2)*j_glat(k,i,lat) +
     |        dvec(i-nlon,lat,3,2)*j_gz(k,i,lat))/hgt1/
     |        dddarr(i-nlon,lat)	
	    endif  ! i <= nlon+1
!	     
! for output only  [kg/m^2/s^2]	     
	     plt_grav(k,i,lat)= jgrav
	     plt_rho(k,i,lat) = rho
	     plt_gne(k,i,lat) = -10.*boltz*(gne*sumt)*fac_hgt**3
	     plt_gt(k,i,lat)  = -10.*boltz*(gt*halfne)*fac_hgt**3
	     plt_ne_horlon(k,i,lat)= -10.*boltz*horlon_ne*fac_horlon/   
     |                (re+halfz)*fac_hgt**3
	     plt_t_horlon(k,i,lat) = -10.*boltz*horlon_t*fac_horlon/   
     |                (re+halfz)*fac_hgt**3
	     plt_ne_horlat(k,i,lat)= -10.*boltz*horlat_ne/(2.*dlatg*  
     |                (re+halfz))*fac_hgt**3
	     plt_t_horlat(k,i,lat) = -10.*boltz*horlat_t/(2.*dlatg* 
     |                (re+halfz))*fac_hgt**3
!     
! test for NE field line integrated 
!            je1oD_geo(k,i,lat) = ne(k,i,lat)
!
! End levels loop:
          enddo   ! lev-loop lev0,lev1-1
!	  
! upper boundary
	  k = lev1
	 je1oD_geo(k,i,lat) = 1.5*je1oD_geo(k-1,i,lat) - 0.5*
     |  	 je1oD_geo(k-2,i,lat) ! extrapolation
	 je2oD_geo(k,i,lat) = 1.5*je2oD_geo(k-1,i,lat) - 0.5*
     |  	 je2oD_geo(k-2,i,lat) ! extrapolation
!
! output	  
	  plt_grav(k,i,lat)= 1.5*plt_grav(k-1,i,lat)-0.5*
     |          plt_grav(k-2,i,lat)
	  plt_rho(k,i,lat) = 1.5*plt_rho(k-1,i,lat)-0.5*
     |          plt_rho(k-2,i,lat)
	  plt_gne(k,i,lat) = 1.5*plt_gne(k-1,i,lat)-0.5*
     |          plt_gne(k-2,i,lat)
	  plt_gt(k,i,lat)  = 1.5*plt_gt(k-1,i,lat)-0.5*
     |          plt_gt(k-2,i,lat)
	  plt_ne_horlon(k,i,lat)= 1.5*plt_ne_horlon(k-1,i,lat)-0.5*
     |          plt_ne_horlon(k-2,i,lat)
	  plt_t_horlon(k,i,lat)= 1.5*plt_t_horlon(k-1,i,lat)-0.5*
     |          plt_t_horlon(k-2,i,lat)
	  plt_ne_horlat(k,i,lat)= 1.5*plt_ne_horlat(k-1,i,lat)-0.5*
     |          plt_ne_horlat(k-2,i,lat)
	  plt_t_horlat(k,i,lat)= 1.5*plt_t_horlat(k-1,i,lat)-0.5*
     |          plt_t_horlat(k-2,i,lat)
!	  
! End longitude loop:
        enddo  ! lon-loop
! End latitude loop:
      enddo  ! lat-loop
!
! periodic point longitude
#ifdef MPI
      ftmp(:,:,:,1) = je1oD_geo(:,lon0:lon1,lat0:lat1)
      ftmp(:,:,:,2) = je2oD_geo(:,lon0:lon1,lat0:lat1)

      call mp_periodic_f3d(ftmp,lev0,lev1,lon0,lon1,lat0,lat1,2)

      je1oD_geo(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,1)
      je2oD_geo(:,lon0:lon1,lat0:lat1) = ftmp(:,:,:,2)
#else	
      do lat=lat0,lat1
        do i=1,2
          je1oD_geo(:,i,lat)	    = je1oD_geo(:,nlon+i,lat)
          je1oD_geo(:,nlonp2+i,lat) = je1oD_geo(:,i+2,lat)
          je2oD_geo(:,i,lat)	    = je2oD_geo(:,nlon+i,lat)
          je2oD_geo(:,nlonp2+i,lat) = je2oD_geo(:,i+2,lat)
        enddo
      enddo
#endif
!      
!     do lat=lat0,lat1
!       call addfld('je1oD_pg_geo' ,'je1oD_pg_geo','A/cm^2',je1oD_geo
!    |	  (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('je2oD_pg_geo' ,'je2oD_pg_geo','A/cm^2',je2oD_geo
!    |	  (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('JGRAV' ,'J(G)','kg/m^2/s^2',plt_grav
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('JGNE' ,'J(GNE_Z)','kg/m^2/s^2',plt_gne
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('JGT' ,'J(GT_Z)','kg/m^2/s^2',plt_gt
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('JNEHLON' ,'J(NE_horlon)','kg/m^2/s^2',
!    |     plt_ne_horlon(:,lon0:lon1,lat),'lev',lev0,lev1,
!    |       'lon',lon0,lon1,lat)
!       call addfld('JTHLON' ,'J(T_horlon)','kg/m^2/s^2',plt_t_horlon
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('JNEHLAT' ,'J(NE_horlat)','kg/m^2/s^2',
!    |     plt_ne_horlat(:,lon0:lon1,lat),'lev',lev0,lev1,
!    |     'lon',lon0,lon1,lat)
!       call addfld('JTHLAT' ,'J(T_horlat)','kg/m^2/s^2',plt_t_horlat
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('RHO_M' ,'RHO_M','g/cm^3',plt_rho
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('OP_M' ,'OP_M','1/cm^3',op
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('O2P_M' ,'O2P_M','1/cm^3',o2p
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('NPLUS_M' ,'NPLUS_M','1/cm^3',nplus
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('N2P_M' ,'N2P_M','1/cm^3',n2p
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('NOP_M' ,'NOP_M','1/cm^3',nop
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('NE_M' ,'NE_M','1/cm^3',ne
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('Z_M' ,'Z_M','cm',z
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('TE_M' ,'TE_M','K',te
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!       call addfld('TI_M' ,'TI_M','K',ti
!    |     (:,lon0:lon1,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     enddo
!      
      end subroutine magpres_grav
!-----------------------------------------------------------------------
      subroutine alloc_pg(lon0,lon1,lat0,lat1)
!
! Allocating of arrays for calculation of 
!   plasma pressure and gravity driven current
!
      integer,intent(in) :: lon0,lon1,lat0,lat1
!
! Local:
      integer :: istat
!
      allocate(je1oD_geo(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_pg: error allocating',
     |  ' je1oD_geo: stat=',i3)") istat
      je1oD_geo = 0.
      allocate(je2oD_geo(nlevp1,lon0:lon1,lat0:lat1),stat=istat)
      if (istat /= 0) write(6,"('>>> alloc_pg: error allocating',
     |  ' je2oD_geo: stat=',i3)") istat
      je2oD_geo = 0.
      end subroutine alloc_pg
!-----------------------------------------------------------------------
      end module magpres_g

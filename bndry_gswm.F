! calculate the lower boundary condition for TN (dt_gswm)
! Z (addiag_gswm) and UN,VN (duv_gswm) 
! if no GSWM-input is provided the corresponding tidal
! modes will be used
!
!-----------------------------------------------------------------------
      subroutine lbc_gswm_dt(tnlbc,lon0,lon1,lat0,lat1)
      
      use input_module,only: step
      use init_module,only:  iter,igetgswmdi,igetgswmsdi,
     |   igetgswmnmidi,igetgswmnmisdi
      use bndry_module,only: tb,tb2,tba,bnd,bnd2,bnda,ci
      use cons_module,only:  freq_semidi,tbound
      use gswm_module,only:  tndi_gswm,tnsdi_gswm,     ! GSWM tides
     |   tnnmidi_gswm,tnnmisdi_gswm 
      
      implicit none
      
      
 ! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1      ! dimension for tnlbc    
      real,intent(out) :: tnlbc(lon0:lon1,lat0:lat1) ! lower boundary condition

! Local:
      integer :: i,lat
      real :: rstep
      complex :: expt,expt2,expta
!     
! Calculate exponentials
!      
      rstep = float(step)
      expt  = cexp(ci*freq_semidi*rstep*iter)
      expt2 = cexp(ci*.5*freq_semidi*rstep*iter)
      expta = 1.
!
      if(igetgswmdi == 1.and.igetgswmsdi == 1) then
	do lat=lat0,lat1
	  do i=lon0,lon1
	    tnlbc(i,lat) = tnsdi_gswm(i,lat)+tbound	 ! semidiurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+ tndi_gswm(i,lat)! diurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tba(lat)*bnda(i)*expta)		 ! annual tide
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
!      
      elseif(igetgswmdi == 0.and.igetgswmsdi == 1) then
	do lat=lat0,lat1
	  do i=lon0,lon1
	    tnlbc(i,lat) = tnsdi_gswm(i,lat)+tbound    ! semidiurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tb2(lat)*bnd2(i)*expt2)	       ! diurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tba(lat)*bnda(i)*expta)	       ! annual tide
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
!      
      elseif(igetgswmdi == 1.and.igetgswmsdi == 0) then
	do lat=lat0,lat1
	  do i=lon0,lon1
	    tnlbc(i,lat) = real(tb(lat)*bnd(i)*expt)+tbound  ! semidiurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+ tndi_gswm(i,lat)    ! diurnal tide
	    tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tba(lat)*bnda(i)*expta)		     ! annual tide
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
      else
! 
        do lat=lat0,lat1
          do i=lon0,lon1
            tnlbc(i,lat) = real(tb(lat)*bnd(i)*expt)+tbound ! semidiurnal tide
            tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tb2(lat)*bnd2(i)*expt2) 		    ! diurnal tide
            tnlbc(i,lat) = tnlbc(i,lat)+
     |        real(tba(lat)*bnda(i)*expta) 		    ! annual tide
          enddo   ! i=lon0,lon1
        enddo	  ! j=lat0,lat1
      endif
!
      if(igetgswmnmidi == 1) then  ! nonmigrating diurnal tide
	do lat=lat0,lat1
	  do i=lon0,lon1
	    tnlbc(i,lat) = tnlbc(i,lat)+ tnnmidi_gswm(i,lat)
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
      endif
!
      if(igetgswmnmisdi == 1) then  ! nonmigrating semidiurnal tide
	do lat=lat0,lat1
	  do i=lon0,lon1
	    tnlbc(i,lat) = tnlbc(i,lat)+ tnnmisdi_gswm(i,lat)
	  enddo   ! i=lon0,lon1
	enddo	 ! j=lat0,lat1
      endif
!      
      
      end subroutine lbc_gswm_dt
!
!-----------------------------------------------------------------------
! calculate geopotential height
! z(1) = zb
!
      subroutine lbc_gswm_addiag(z,lon0,lon1,lat0,lat1)
!      
      use input_module,only: step
      use init_module,only:  iter,igetgswmdi,igetgswmsdi,
     |   igetgswmnmidi,igetgswmnmisdi
      use bndry_module,only: zb,zb2,zba,bnd,bnd2,bnda,ci
      use cons_module,only:  freq_semidi,dt
      use gswm_module,only:  zdi_gswm,zsdi_gswm,        ! GSWM tides
     |   znmidi_gswm,znmisdi_gswm
!      
      implicit none
!     
! Args:
      integer,intent(in) :: lon0,lon1,lat0,lat1     ! dimension for tnlbc    
      real,intent(out) :: z(lon0:lon1,lat0:lat1) ! lower boundary condition

! Local:
      integer :: i,j
      real :: rstep
      complex :: expt,expt2,expta
!
! Calculate exponentials
!      
      expt=cexp(ci*freq_semidi*dt*iter)
      expt = cexp(ci*.5*freq_semidi*dt*iter)
      expt = 1.
!
      if(igetgswmdi == 1.and.igetgswmsdi == 1) then
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = zsdi_gswm(i,j)		  ! semidiurnal tide
	    z(i,j) = z(i,j)+ zdi_gswm(i,j)	    ! diurnal tide
	    z(i,j) = z(i,j)+real(zba(j)*bnda(i)*expt) ! annual tide
	  enddo
	enddo
      elseif(igetgswmdi == 0.and.igetgswmsdi == 1) then
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = zsdi_gswm(i,j)		  ! semidiurnal tide
	    z(i,j) = z(i,j)+real(zb2(j)*bnd2(i)*expt)! diurnal tide
	    z(i,j) = z(i,j)+real(zba(j)*bnda(i)*expt) ! annual tide
	  enddo
	enddo
      elseif(igetgswmdi == 1.and.igetgswmsdi == 0) then
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = real(zb(j)*bnd(i)*expt)	      ! semidiurnal tide
	    z(i,j) = z(i,j)+ zdi_gswm(i,j)	      ! diurnal tide
	    z(i,j) = z(i,j)+real(zba(j)*bnda(i)*expt) ! annual tide
	  enddo
	enddo
      else
        do j = lat0,lat1
          do i=lon0,lon1
            z(i,j) = real(zb(j)*bnd(i)*expt)	      ! semidiurnal tide
            z(i,j) = z(i,j)+real(zb2(j)*bnd2(i)*expt) ! diurnal tide
            z(i,j) = z(i,j)+real(zba(j)*bnda(i)*expt) ! annual tide
          enddo
        enddo
      endif
!     
      if(igetgswmnmidi == 1) then	! nonmigrating diurnal tide
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = z(i,j)+ znmidi_gswm(i,j)
	  enddo
	enddo
      endif
!     
      if(igetgswmnmisdi == 1) then	! nonmigrating semidiurnal tide
	do j = lat0,lat1
	  do i=lon0,lon1
	    z(i,j) = z(i,j)+ znmisdi_gswm(i,j)
	  enddo
	enddo
      endif
!      
      end subroutine lbc_gswm_addiag
!-----------------------------------------------------------------------
      subroutine lbc_gswm_duv(unlbc,vnlbc,unlbc_diag,vnlbc_diag,
     |    lon0,lon1,lat,lev0,lev1,expt,expt2,expta)
      
      use init_module,only:  igetgswmdi,igetgswmsdi,
     |   igetgswmnmidi,igetgswmnmisdi
      use bndry_module,only: ub,ub2,uba,vb,vb2,vba,bnd,bnd2,bnda,ci
      use gswm_module,only:  			! GSWM tides
     |   undi_gswm,unsdi_gswm,unnmidi_gswm,unnmisdi_gswm,
     |   vndi_gswm,vnsdi_gswm,vnnmidi_gswm,vnnmisdi_gswm
      
      implicit none
      
      
 ! Args:
      integer,intent(in) :: lon0,lon1,lat,lev0,lev1  ! dimension for lbc    
      complex,intent(in) :: expt,expt2,expta 
      real,intent(out) :: unlbc(lon0:lon1),	! lower boundary condition
     |                    vnlbc(lon0:lon1), 
     |                    unlbc_diag(lev0:lev1,lon0:lon1), ! for output
     |                    vnlbc_diag(lev0:lev1,lon0:lon1)
! Local:
      integer :: i
!     
! 
      if(igetgswmdi == 1.and.igetgswmsdi == 1) then
        do i=lon0,lon1
          unlbc(i) = unsdi_gswm(i,lat)  	       ! semidiurnal tide
          vnlbc(i) = vnsdi_gswm(i,lat)
          unlbc(i) = unlbc(i) + undi_gswm(i,lat)           ! diurnal tide
          vnlbc(i) = vnlbc(i) + vndi_gswm(i,lat)
          unlbc(i) = unlbc(i) + real(uba(lat)*bnda(i)*expta)! annual tide
          vnlbc(i) = vnlbc(i) + real(vba(lat)*bnda(i)*expta)
!
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      elseif(igetgswmdi == 0.and.igetgswmsdi == 1) then
        do i=lon0,lon1
          unlbc(i) = unsdi_gswm(i,lat)  	         ! semidiurnal tide
          vnlbc(i) = vnsdi_gswm(i,lat)
          unlbc(i) = unlbc(i) + real(ub2(lat)*bnd2(i)*expt2) ! diurnal tide
          vnlbc(i) = vnlbc(i) + real(vb2(lat)*bnd2(i)*expt2)
          unlbc(i) = unlbc(i) + real(uba(lat)*bnda(i)*expta) ! annual tide
          vnlbc(i) = vnlbc(i) + real(vba(lat)*bnda(i)*expta)
!
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      elseif(igetgswmdi == 1.and.igetgswmsdi == 0) then
        do i=lon0,lon1
          unlbc(i) =  real(ub(lat)*bnd(i)*expt)         ! semidiurnal tide
          vnlbc(i) =  real(vb(lat)*bnd(i)*expt)
          unlbc(i) = unlbc(i) + undi_gswm(i,lat)            ! diurnal tide
          vnlbc(i) = vnlbc(i) + vndi_gswm(i,lat)
          unlbc(i) = unlbc(i) + real(uba(lat)*bnda(i)*expta)! annual tide
          vnlbc(i) = vnlbc(i) + real(vba(lat)*bnda(i)*expta)
!
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      else
        do i=lon0,lon1
          unlbc(i) = real(ub(lat)*bnd(i)*expt)          ! Semidiurnal tide
          vnlbc(i) = real(vb(lat)*bnd(i)*expt)
          unlbc(i) = unlbc(i) + real(ub2(lat)*bnd2(i)*expt2)! diurnal tide
          vnlbc(i) = vnlbc(i) + real(vb2(lat)*bnd2(i)*expt2)
          unlbc(i) = unlbc(i) + real(uba(lat)*bnda(i)*expta)! annual tide
          vnlbc(i) = vnlbc(i) + real(vba(lat)*bnda(i)*expta)
!
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      endif
!    
      if(igetgswmnmidi == 1) then      ! nonmigrating diurnal tide
        do i=lon0,lon1
          unlbc(i) = unlbc(i) + unnmidi_gswm(i,lat)
          vnlbc(i) = vnlbc(i) + vnnmidi_gswm(i,lat)
!
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      endif
!    
      if(igetgswmnmisdi == 1) then      ! nonmigrating semidiurnal tide
        do i=lon0,lon1
          unlbc(i) = unlbc(i) + unnmisdi_gswm(i,lat)
          vnlbc(i) = vnlbc(i) + vnnmisdi_gswm(i,lat)
!
          unlbc_diag(:,i) = unlbc(i)	! for output
          vnlbc_diag(:,i) = vnlbc(i)
        enddo ! i=lon0,lon1
      endif
!		
      end subroutine lbc_gswm_duv

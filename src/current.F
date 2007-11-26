!
#include <defs.h>
!
! bf 6/02: Modified for tiegcm1 from Astrid Maute tgcm15 mod files 
!          nscoef.F, nscrdens.F, and nscrrt.F.
!
! Subs nosocoef, nosocrrt, and nosocrdens are called from the dynamo
!   module (dynamo.F). These routines are not in a module so as to avoid 
!   circular dependency with the dynamo module.
!
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: nosocoef
! !INTERFACE:
      subroutine nosocoef
!
! !DESCRIPTION:
! am_02/02: calculate the coefficient stencil for both hemisphere which 
!   is used to calculate the height integrated current densities 
!   K_(q,phi)K_(q,lam)
!      
! !USES:
      use params_module,only: nmlon,nmlonp1,nmlat,nmlath
      use cons_module,only: dlonm,dlatm,pi,r0
      use dynamo_module,only: nmlon0,zigm11,zigmc,zigm2,zigm22,rim,
     |  nscoef,unitvm
      use addfld_module,only: addfld
      implicit none
! !PARAMETERS:
! !RETURN VALUE:
!
! !REVISION HISTORY:
! 05.03.15  <Astrid Maute> <include header> 
! 
! EOP 
!
! Calculate coefficients for dynamo pde for both hemisphere
!     
! Local:
      real :: array(0:nmlon0+1,nmlat),cs(nmlat)
      real :: nszigmc(nmlonp1,nmlat),nszigm2(nmlonp1,nmlat)
      real :: nszigm11(nmlonp1,nmlat),nszigm22(nmlonp1,nmlat)
      real :: nsrhs(nmlonp1,nmlat)
      real :: dfac,dfac1n,dfac1s,dfac2n,dfac2s
      integer :: j,je,jj,jjj,i,n
! Externals:
      real,external :: sddot ! in util.F
!
! Clear arrays
      nszigmc(:,:)  = 0.0 
      nszigm2(:,:)  = 0.0 
      nszigm11(:,:) = 0.0 
      nszigm22(:,:) = 0.0 
      nsrhs(:,:)    = 0.0 
! 
! Calculate magnetic latitude cosin array
!
      do j = 1,nmlat     ! -pi/2 to pi/2
        cs(j) = cos(-pi/2.+(j-1)*dlatm)
      enddo
!
! Reverse sign of ZIGMC to be compatible with Cicely's (richmond)
! Calculate difference
!   for zigmc & zigm2 no sign change because values from the
!                     corresponding hemisphere are used
!   zigmc : sum_{phi lam} = +/(-) (sum_H - sum_C) -> C Ridley
!   zigm2 : sum_{lam phi} = -/(+) (sum_H + sum_C) -> D Ridley
!   zigm11: sum_{phi phi}   -> A Ridley
!   zigm22: sum_{lam lam}   -> B Ridley
!   
!   factors from difference quotients and derivatives
!       4.*dlatm*dlonm for mixed terms zigmc and zigm2
!       dlatm**2 or dlonm**2 for zigm22 and zigm11
!   
!   factor cosin (cs and 1/cs) from pde for zigm22 and zigm11
!
      dfac = 4.*dlatm*dlonm
      do j = 2,nmlat-1    ! 2,96  not value at the poles 
	dfac1n = cs(j)/dlatm**2
	dfac2n = cs(j)*dlonm**2
        do i = 1,nmlonp1
	  nszigmc(i,j)  = -zigmc(i,j)
          nszigmc(i,j)  = (nszigmc(i,j)+zigm2(i,j))/dfac   
          nszigm2(i,j)  = nszigmc(i,j)-2.*zigm2(i,j)/dfac
          nszigm22(i,j) = zigm22(i,j)*dfac1n
          nszigm11(i,j) = zigm11(i,j)/dfac2n
        enddo
      enddo
!
! Change sign for values at equator maybe not necessary, but
!   then the sign for coefficient has to be changed too
!
      j = (nmlat+1)/2.0   
      do i = 1,nmlonp1
	nszigmc(i,j)  = -nszigmc(i,j)
	nszigm2(i,j)  = -nszigm2(i,j)
      enddo
!
! Values at the poles (1 and 97)
!
      jj  = nmlat
      jjj = 1
      dfac1n = cs(jj)/dlatm**2
      dfac1s = cs(jjj)/dlatm**2   ! is not necessary cos symmetric
!
      do i = 1,nmlonp1
	nszigmc(i,1)     = -zigmc(i,1)            ! 1
	nszigmc(i,nmlat) = -zigmc(i,nmlat)	  ! 97
        nszigmc(i,jj)    = (nszigmc(i,jj)+zigm2(i,jj))/dfac   
        nszigm2(i,jj)    =  nszigmc(i,jj)-2.*zigm2(i,jj)/dfac
        nszigm22(i,jj)   =  zigm22(i,jj)*dfac1n
        nszigmc(i,jjj)   = (nszigmc(i,jjj)+zigm2(i,jjj))/dfac
        nszigm2(i,jjj)   =  nszigmc(i,jjj)-2.*zigm2(i,jjj)/dfac
        nszigm22(i,jjj)  =  zigm22(i,jjj)*dfac1s
!
! Set zigm11 to zero at the magnetic poles (1 and 97) to avoid floating 
!   point exception (values at the poles are not used)
!
        nszigm11(i,jj)  = 0.0
        nszigm11(i,jjj) = 0.0
      enddo
!
! Clear array for difference stencil over north and south hemisphere 
!
      nscoef(:,:,:) = 0.0
!  
! Calculate contribution to stencil from each pde coefficient
!   one at a time because of smaller working arrays 
!
      call nsstencil(nszigm11,nmlon0,nmlat,nscoef,array,1,nmlath)
      call nsstencil(nszigm22,nmlon0,nmlat,nscoef,array,4,nmlath)
      call nsstencil(nszigmc ,nmlon0,nmlat,nscoef,array,2,nmlath)
      call nsstencil(nszigm2 ,nmlon0,nmlat,nscoef,array,3,nmlath) 
!
! Set boundary conditions at pole 
!   value change from 1.0 to 0.5 for each hemisphere
!
      do i = 1,nmlon0
        do n = 1,8
          nscoef(i,nmlat,n) = 0.
          nscoef(i,1,n)     = 0.
        enddo
        nscoef(i,nmlat,9) = 0.5
        nscoef(i,1,9)     = 0.5
      enddo
!
! Divide stencil by cos(theta)
!
      do j = 2,nmlat-1
        do n = 1,9
	  nscoef(:,j,n) = nscoef(:,j,n)/cs(j)
	enddo
      enddo 
! 
! Calculate right hand side of pde from rim(1) and rim(2)
!
      do j = 2,nmlath-1  ! 2,48  south pole-1 to equator-1
        jj = j+nmlath-1     ! 50,96 equator-1 to north pole-1 
!  
! Differentiate rim(1) w.r.t lambda
!  
        do i = 2,nmlon-1
          nsrhs(i,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim(i+1,j,1)-rim(i-1,j,1))
          nsrhs(i,jj) = 1.0/(dlonm*cs(jj))*
     |                      0.5*(rim(i+1,jj,1)-rim(i-1,jj,1))
	enddo
!  
! Values at the poles
!  
        nsrhs(1,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim(2,j,1)-rim(nmlon,j,1))
        nsrhs(1,jj) = 1.0/(dlonm*cs(jj))*
     |                      0.5*(rim(2,jj,1)-rim(nmlon,jj,1))
        nsrhs(nmlon,j)  = 1.0/(dlonm*cs(j))*
     |                      0.5*(rim(1,j,1)-rim(nmlon-1,j,1))
        nsrhs(nmlon,jj) = 1.0/(dlonm*cs(jj))*
     |                      0.5*(rim(1,jj,1)-rim(nmlon-1,jj,1))
      enddo
!  
! Differentiate rim(2) w.r.t theta0
!  
      do j = 2,nmlath-1  ! 2,48  south pole -1 to equator-1
	jj = j+nmlath-1  ! 50,96 equator+1 to north pole-1 
	do i = 1,nmlon
	 nsrhs(i,j)  = nsrhs(i,j) - 1.0/(dlatm*cs(j))*0.5*
     |  	      (rim(i,j+1,2)*cs(j+1)-rim(i,j-1,2)*cs(j-1))
	 nsrhs(i,jj) = nsrhs(i,jj) + 1.0/(dlatm*cs(jj))*0.5*
     |  	      (rim(i,jj+1,2)*cs(jj+1)-rim(i,jj-1,2)*cs(jj-1))
       enddo
      enddo
!    
! Calculate value at the poles by averaging over i:nmlon
!    
      nsrhs(1,nmlat) = -2./float(nmlon)*
     |                  sddot(nmlon,unitvm,rim(1,nmlat-1,2))/cs(nmlat-1)
      nsrhs(1,1)     = -2./float(nmlon)*
     |                  sddot(nmlon,unitvm,rim(1,2,2))/cs(2)
!
! Extend over longitude
!
      nsrhs(:,nmlat) = nsrhs(1,nmlat)
      nsrhs(:,1)     = nsrhs(1,1)
!
! note: for calculating J_mr values with the stencil at the equator not used
! note: for the test case tstjmrim when both hemisphere are added together
!        get double values at the equator, since no seperate value for south and
!        north of the equator, at the equator jump in rhs, therefore values
!        doubled at equator
! note: nsrhs stencil is the same as rhs in transf.f if average (north & south
!       of equator is taken for derivative in lam direction 
! note: introduced 0.5 to fit coefficient stencil = double nscoef(49)
!       for consistency also for nsrhs: c0(j=49,10) =  double nsrhs(49)
!
      je = nmlath
      i = 1
      nsrhs(i,je) = 0.5/dlonm*(rim(i+1,je,1)-rim(nmlon,je,1))
      nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |  		rim(i,je,2)+ cs(je+1)*rim(i,je+1,2))
      nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |  		rim(i,je,2)+ cs(je-1)*rim(i,je-1,2))
      nsrhs(i,je) = 0.5*nsrhs(i,je)
!
      do i = 2,nmlon-1
        nsrhs(i,je) = 0.5/dlonm*(rim(i+1,je,1)-rim(i-1,je,1))
	nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |                    rim(i,je,2)+ cs(je+1)*rim(i,je+1,2))
	nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |                    rim(i,je,2)+ cs(je-1)*rim(i,je-1,2))
        nsrhs(i,je) = 0.5*nsrhs(i,je)
      enddo
!
      i = nmlon
      nsrhs(i,je) = 0.5/dlonm*(rim(1,je,1)-rim(i-1,je,1))
      nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |  		rim(i,je,2)+ cs(je+1)*rim(i,je+1,2))
      nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)*
     |  		rim(i,je,2)+ cs(je-1)*rim(i,je-1,2))
      nsrhs(i,je) = 0.5*nsrhs(i,je)
!
! Periodic points
!
      nsrhs(nmlonp1,:) = nsrhs(1,:)
!
! Scale rhs by refernce radius (R_E + H0) in meters dfac = r0*1e-2
!
      dfac = r0*1.0e-2
      nsrhs(:,:) = nsrhs(:,:)*dfac
!
! Insert nsrhs into coefficient : from south to north pole
!   and divide by cos(theta) =1.0 not necessary !
!
      nscoef(:,:,10) = nsrhs(:,:)
!
! Set value of solution to 1.0 at poles ! change control later with equilibrium
!
      nscoef(:,nmlat,10) = 0.5
      nscoef(:,1,10)	 = 0.5
!
      end subroutine nosocoef
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: nsstencil
! !INTERFACE:
      subroutine nsstencil(zigm,nmlon0,nmlat,nscoef,array,ncoef,nmlath)
      implicit none
!
! !DESCRIPTION: set up the differencing stencil for coefficient zigm. 
!     the output is the stencil nscoef
!     in contrast to subroutine stencil in dynamo module it's done
!     for both hemispheres
!
! !PARAMETERS: 
      integer,intent(in):: nmlon0,nmlat,ncoef,nmlath
      real,intent(in)   :: zigm(nmlon0,nmlat)
! !RETURN VALUE:
      real,intent(out)  :: 
     |  nscoef(nmlon0,nmlat,10),array(0:nmlon0+1,nmlat)
!
! !REVISION HISTORY:
! 05.03.09  <Astrid Maute> <include header> 
! 
! EOP 
!
! Local:
      integer :: i,j
!
! Copy coefficients (south to north pole) into array (south to north pole)
! 
      do j = 1,nmlat    ! 1,97
        do i = 1,nmlon0
	  array(i,j) = zigm(i,j)
        enddo
      enddo
!
! Extend one additional grid point on both sides of the array for
! differencing
!     
      i = 1
      do j = 1,nmlat
        array(i-1,j)	  = array(nmlon0-i,j)
        array(nmlon0+i,j) = array(1+i,j)
      enddo
!      
! Calculate contribution to stencil for each grid point
!  
      call nscnm(array,nmlon0,nmlat,nscoef,ncoef,nmlath)
!  
      end subroutine nsstencil
!----------------------------------------------------------------------------
! BOP
! !IROUTINE: nscnm
! !INTERFACE:
      subroutine nscnm(array,nmlon0,nmlat,nsc,ncoef,nmlath)
      implicit none
!      
! !DESCRIPTION: 
!      calculate contribution for each zigm
!      the diagonial dominance of the stencil is not checked
!      since the coefficients are not used for solving
!      one reason for the difference between added north-south
!      and seperated north south hemisphere
!
!      stencil for southern hemisphere changed, since
!      also for the southern hemisphere latudinal counts from
!      the equator to the south pole as for the northern hemisphere
!      from the equator to the north pole    
!      nsc(i,j,n) n:
!      northern hemisphere stencil n=(1,2,3,4,5,6,7,8,9)
!      southern hemisphere stencil n=(1,-8,-7,-5,-4,-2,9)
!
!      values at the equator (j=49): not separeted for northern and southern
!      hemisphere- only one stencil and later in nscrrt equally  
!      distributed to northern and southern hemisphere
!
! !PARAMETERS: 
      integer,intent(in) :: nmlon0,nmlat,ncoef,nmlath
      real,intent(in)    :: array(0:nmlon0+1,nmlat)
! !RETURN VALUE:
      real,intent(out)   :: nsc(nmlon0,nmlat,10)
!
!
! !REVISION HISTORY:
! 05.03.09  <Astrid Maute> <include header> 
! 
! EOP  
! 
! Local:
      integer :: i,j,jj,jjj
!
! Check: what's about boundary equator values for nsc, why not used?
! Calculate contribution for zigm11
!
      if(ncoef.eq.1)then	
        do j = 2,nmlath-1   ! 2,48  not value at the pols (jj=97 jjj=1)
          jj  = nmlath+j-1      ! 50,96 northern hemisphere
          jjj = nmlath-j+1      ! 48,2  southern hemisphere
!
          do i = 1,nmlon0
	    nsc(i,jj,1) = nsc(i,jj,1)+.5*(array(i,jj)+array(i+1,jj))
	    nsc(i,jj,5) = nsc(i,jj,5)+.5*(array(i,jj)+array(i-1,jj))
	    nsc(i,jj,9) = nsc(i,jj,9)-.5*(array(i+1,jj)+2.*array(i,jj)+
     |                                 array(i-1,jj))
	    nsc(i,jjj,1) = nsc(i,jjj,1)+.5*(array(i,jjj)+array(i+1,jjj))
	    nsc(i,jjj,5) = nsc(i,jjj,5)+.5*(array(i,jjj)+array(i-1,jjj))
	    nsc(i,jjj,9) = nsc(i,jjj,9)-.5*(array(i+1,jjj)+
     |                                 2.*array(i,jjj)+array(i-1,jjj))
          enddo
	enddo
!
! am 2001-6-27 include boundary condition at equator
!	
        jj  = nmlath      ! 49
        do i = 1,nmlon0
	  nsc(i,jj,1) = nsc(i,jj,1)+.5*(array(i,jj)+array(i+1,jj))
	  nsc(i,jj,5) = nsc(i,jj,5)+.5*(array(i,jj)+array(i-1,jj))
	  nsc(i,jj,9) = nsc(i,jj,9)-.5*(array(i+1,jj)+2.*array(i,jj)+
     |  			     array(i-1,jj))
        enddo
!
! Calculate contribution for zigm12 (=ZIGMC+ZIGM2)
! 
      else if(ncoef.eq.2) then
        do j = 2,nmlath-1       ! 2,48  not value at the pols (jj=97 jjj=1)
          jj  = nmlath+j-1      ! 50,96 northern hemisphere
          jjj = nmlath-j+1      ! 48,2  southern hemisphere
!
          do i = 1,nmlon0 
	    nsc(i,jj,2) = nsc(i,jj,2)+.5*(array(i,jj)+array(i+1,jj))
	    nsc(i,jj,4) = nsc(i,jj,4)-.5*(array(i,jj)+array(i-1,jj))
	    nsc(i,jj,6) = nsc(i,jj,6)+.5*(array(i,jj)+array(i-1,jj))
	    nsc(i,jj,8) = nsc(i,jj,8)-.5*(array(i,jj)+array(i+1,jj))
	    nsc(i,jj,3) = nsc(i,jj,3)+.5*(-array(i-1,jj)+array(i+1,jj))
	    nsc(i,jj,7) = nsc(i,jj,7)-.5*(-array(i-1,jj)+array(i+1,jj))
	    
	    nsc(i,jjj,2)=nsc(i,jjj,2)+.5*(array(i,jjj)+array(i+1,jjj))
	    nsc(i,jjj,4)=nsc(i,jjj,4)-.5*(array(i,jjj)+array(i-1,jjj))
	    nsc(i,jjj,6)=nsc(i,jjj,6)+.5*(array(i,jjj)+array(i-1,jjj))
	    nsc(i,jjj,8)=nsc(i,jjj,8)-.5*(array(i,jjj)+array(i+1,jjj))
	    nsc(i,jjj,3)=nsc(i,jjj,3)+.5*(-array(i-1,jjj)+array(i+1,jjj))
	    nsc(i,jjj,7)=nsc(i,jjj,7)-.5*(-array(i-1,jjj)+array(i+1,jjj))
           enddo
	enddo
!
! Calculate contribution for zigm21 (=ZIGMC-ZIGM2)
! 
      else if(ncoef.eq.3) then
        do j = 3,nmlath-1   ! 3,48  not value at the pols (jj=97 jjj=1)
          jj  = nmlath+j-1      ! 51,96 northern hemisphere
          jjj = nmlath-j+1      ! 47,2  southern hemisphere
!
          do i = 1,nmlon0 
	    nsc(i,jj,2) = nsc(i,jj,2)+.5*(array(i,jj)+array(i,jj+1))
	    nsc(i,jj,4) = nsc(i,jj,4)-.5*(array(i,jj)+array(i,jj+1))
	    nsc(i,jj,6) = nsc(i,jj,6)+.5*(array(i,jj)+array(i,jj-1))
	    nsc(i,jj,8) = nsc(i,jj,8)-.5*(array(i,jj)+array(i,jj-1))
	    nsc(i,jj,1) = nsc(i,jj,1)+.5*(array(i,jj+1)-array(i,jj-1))
	    nsc(i,jj,5) = nsc(i,jj,5)-.5*(array(i,jj+1)-array(i,jj-1))
!
	    nsc(i,jjj,8)=nsc(i,jjj,8)-.5*(array(i,jjj)+array(i,jjj+1))
	    nsc(i,jjj,6)=nsc(i,jjj,6)+.5*(array(i,jjj)+array(i,jjj+1))
	    nsc(i,jjj,4)=nsc(i,jjj,4)-.5*(array(i,jjj)+array(i,jjj-1))
	    nsc(i,jjj,2)=nsc(i,jjj,2)+.5*(array(i,jjj)+array(i,jjj-1))
  	    nsc(i,jjj,1)=nsc(i,jjj,1)-.5*(array(i,jjj+1)-array(i,jjj-1))
	    nsc(i,jjj,5)=nsc(i,jjj,5)+.5*(array(i,jjj+1)-array(i,jjj-1))
          enddo
	enddo
!
	j = 2	! 2 change in sign for equatorial values 
	jj  = nmlath+j-1      ! 50 northern hemisphere
	jjj = nmlath-j+1      ! 48  southern hemisphere
!
       do i = 1,nmlon0   !am2001-7-3
 	 nsc(i,jj,2) = nsc(i,jj,2)+.5*(array(i,jj)+array(i,jj+1))
	 nsc(i,jj,4) = nsc(i,jj,4)-.5*(array(i,jj)+array(i,jj+1))
	 nsc(i,jj,6) = nsc(i,jj,6)+.5*(array(i,jj)-array(i,jj-1))
	 nsc(i,jj,8) = nsc(i,jj,8)-.5*(array(i,jj)-array(i,jj-1))
	 nsc(i,jj,1) = nsc(i,jj,1)+.5*(array(i,jj+1)+array(i,jj-1))
	 nsc(i,jj,5) = nsc(i,jj,5)-.5*(array(i,jj+1)+array(i,jj-1))
!
	 nsc(i,jjj,8)=nsc(i,jjj,8)-.5*(array(i,jjj)-array(i,jjj+1))
	 nsc(i,jjj,6)=nsc(i,jjj,6)+.5*(array(i,jjj)-array(i,jjj+1))
	 nsc(i,jjj,4)=nsc(i,jjj,4)-.5*(array(i,jjj)+array(i,jjj-1))
	 nsc(i,jjj,2)=nsc(i,jjj,2)+.5*(array(i,jjj)+array(i,jjj-1))
	 nsc(i,jjj,1)=nsc(i,jjj,1)-.5*(-array(i,jjj+1)-array(i,jjj-1))
	 nsc(i,jjj,5)=nsc(i,jjj,5)+.5*(-array(i,jjj+1)-array(i,jjj-1))
       enddo
!    
! Low latitude boundary conditions: contribution from zigm21(i,j-1/2)=0
!
	j = 1
        jj  = nmlath+j-1      ! 49
        jjj = nmlath-j+1      ! 49
!
	do i = 1,nmlon0 
	  nsc(i,jj,2) = nsc(i,jj,2)+.25*(-array(i,jj)+array(i,jj+1))
	  nsc(i,jj,4) = nsc(i,jj,4)-.25*(-array(i,jj)+array(i,jj+1))
	  nsc(i,jj,1) = nsc(i,jj,1)+.25*(-array(i,jj)+array(i,jj+1))
	  nsc(i,jj,5) = nsc(i,jj,5)-.25*(-array(i,jj)+array(i,jj+1))
	
	  nsc(i,jj,2) = nsc(i,jj,2)+.25*(-array(i,jjj)+array(i,jjj-1))
	  nsc(i,jj,4) = nsc(i,jj,4)-.25*(-array(i,jjj)+array(i,jjj-1))
	  nsc(i,jj,1) = nsc(i,jj,1)-.25*(array(i,jjj)-array(i,jjj-1))
	  nsc(i,jj,5) = nsc(i,jj,5)+.25*(array(i,jjj)-array(i,jjj-1))
	enddo
!
! Calculate contribution for zigm22 
! 
      else if(ncoef.eq.4) then
        do j = 2,nmlath-1   ! 2,48  not value at the pols (jj=97 jjj=1)
          jj  = nmlath+j-1      ! 50,96 northern hemisphere
          jjj = nmlath-j+1      ! 48,2  southern hemisphere
!
          do i = 1,nmlon0 
	    nsc(i,jj,3) = nsc(i,jj,3)+.5*(array(i,jj)+array(i,jj+1))
	    nsc(i,jj,7) = nsc(i,jj,7)+.5*(array(i,jj)+array(i,jj-1))
	    nsc(i,jj,9) = nsc(i,jj,9)-.5*(array(i,jj-1)+2.0*array(i,jj)+
     |                                 array(i,jj+1))
	    nsc(i,jjj,7)=nsc(i,jjj,7)+.5*(array(i,jjj)+array(i,jjj+1))
	    nsc(i,jjj,3)=nsc(i,jjj,3)+.5*(array(i,jjj)+array(i,jjj-1))
	    nsc(i,jjj,9)=nsc(i,jjj,9)-.5*(array(i,jjj-1)+
     |                                 2.0*array(i,jjj)+ array(i,jjj+1))
           enddo
	enddo
!    
! Low latitude boundary conditions: contribution from zigm22(i,j-1/2)=0
!
	j = 1
        jj  = nmlath+j-1      ! 49
        jjj = nmlath-j+1      ! 49
	do i = 1,nmlon0 
	   nsc(i,jj,3) = nsc(i,jj,3)+.25*(array(i,jj)+array(i,jj+1))
	   nsc(i,jj,9) = nsc(i,jj,9)-.25*(array(i,jj)+array(i,jj+1))
!  am 2001-7-02 otherwise double coefficients
	   nsc(i,jj,3)=nsc(i,jj,3)+.25*(array(i,jjj)+array(i,jjj-1))
	   nsc(i,jj,9)=nsc(i,jj,9)-.25*(array(i,jjj-1)+array(i,jjj))
        enddo
      endif
!      
      end subroutine nscnm
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: nosocrrt
! !INTERFACE:
      subroutine nosocrrt
!
! !DESCRIPTION:
! am_02/02: 
! Calculate current for both hemisphere: 
! [stencil*potential -RHS] = R**2 * J_mr / dt0dts / rcos0s
!
! !USES:
      use params_module,only: nmlon,nmlonp1,nmlat,nmlath,nlev,nlevp1
      use cons_module,only: dt0dts,rcos0s,pi,r0,re,ylatm
      use dynamo_module,only: nscrrt,nscoef,phim,nmlon0,unitvm
      use addfld_module,only: addfld
      implicit none
!
! !PARAMETERS:
! !RETURN VALUE:
!
! !REVISION HISTORY:
! 05.03.15  <Astrid Maute> <include header> 
! 
! EOP 
!
! Local:
      real ::    vtmp,r0sq,fac,facmax,lat,lmin,lmax,pol
      real :: tout(nmlonp1,nmlat,-2:nlevp1)
      integer :: j,jj,jjj,i,k,n,jmod,jmin,jmax
!
! External:
      real,external :: sddot ! in util.F
!  
      nscrrt(:,:) = 0.0
!
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
        nscrrt(i,jj) =  	      nscoef(i,jj,1)*phim(i+1,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,2)*phim(i+1,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,3)*phim(i  ,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,4)*phim(nmlon0-1,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,5)*phim(nmlon0-1,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,6)*phim(nmlon0-1,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,7)*phim(i  ,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,8)*phim(i+1,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,9)*phim(i  ,jj)
        nscrrt(i,jj) = nscrrt(i,jj) - nscoef(i,jj,10)
!           
! Southern hemisphere
!
        nscrrt(i,jjj)=  	    nscoef(i,jjj,1)*phim(i+1,jjj)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,2)*phim(i+1,jjj-1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,3)*phim(i  ,jjj-1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,4)*phim(nmlon0-1,jjj-1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,5)*phim(nmlon0-1,jjj)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,6)*phim(nmlon0-1,jjj+1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,7)*phim(i  ,jjj+1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,8)*phim(i+1,jjj+1)
        nscrrt(i,jjj)=nscrrt(i,jjj)+nscoef(i,jjj,9)*phim(i  ,jjj)
	nscrrt(i,jjj)=nscrrt(i,jjj)-nscoef(i,jjj,10)
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
	 nscrrt(i,jj) =  	      nscoef(i,jj,1)*phim(i+1,jj)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,2)*phim(i+1,jj+1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,3)*phim(i  ,jj+1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,4)*phim(i-1,jj+1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,5)*phim(i-1,jj)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,6)*phim(i-1,jj-1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,7)*phim(i  ,jj-1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,8)*phim(i+1,jj-1)
	 nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,9)*phim(i  ,jj)
	 nscrrt(i,jj) = nscrrt(i,jj) - nscoef(i,jj,10)
!
! Southern hemisphere
!
	 nscrrt(i,jjj) = 		nscoef(i,jjj,1)*phim(i+1,jjj)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,2)*phim(i+1,jjj-1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,3)*phim(i  ,jjj-1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,4)*phim(i-1,jjj-1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,5)*phim(i-1,jjj)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,6)*phim(i-1,jjj+1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,7)*phim(i  ,jjj+1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,8)*phim(i+1,jjj+1)
	 nscrrt(i,jjj) = nscrrt(i,jjj) + nscoef(i,jjj,9)*phim(i  ,jjj)
	 nscrrt(i,jjj) = nscrrt(i,jjj) - nscoef(i,jjj,10)
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
        nscrrt(i,jj) =		    nscoef(i,jj,1)*phim(2,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,2)*phim(2,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,3)*phim(i  ,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,4)*phim(i-1,jj+1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,5)*phim(i-1,jj)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,6)*phim(i-1,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,7)*phim(i  ,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,8)*phim(2,jj-1)
        nscrrt(i,jj) = nscrrt(i,jj) + nscoef(i,jj,9)*phim(i  ,jj)
	nscrrt(i,jj) = nscrrt(i,jj) - nscoef(i,jj,10)
! 
! Southern hemisphere
!
        nscrrt(i,jjj) =  	     nscoef(i,jjj,1)*phim(2,jjj)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,2)*phim(2,jjj-1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,3)*phim(i  ,jjj-1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,4)*phim(i-1,jjj-1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,5)*phim(i-1,jjj)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,6)*phim(i-1,jjj+1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,7)*phim(i  ,jjj+1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,8)*phim(2,jjj+1)
        nscrrt(i,jjj) = nscrrt(i,jjj) +nscoef(i,jjj,9)*phim(i  ,jjj)
        nscrrt(i,jjj) = nscrrt(i,jjj) -nscoef(i,jjj,10)
      enddo
!
! Poles
!      
      nscrrt(1,1) = (9.*sddot(nmlon,unitvm,nscrrt(1,2))-
     |                   sddot(nmlon,unitvm,nscrrt(1,3)))/
     |                   (8.*float(nmlon))
      nscrrt(1,nmlat) = (9.*sddot(nmlon,unitvm,nscrrt(1,nmlat-1))-
     |                   sddot(nmlon,unitvm,nscrrt(1,nmlat-2)))/
     |                   (8.*float(nmlon))
      nscrrt(:,1)     = nscrrt(1,1)		! extend in longitude
      nscrrt(:,nmlat) = nscrrt(1,nmlat)
      do j = 1,nmlat
       lat = ylatm(j)
!       if(lat.le.-pi/18.) then
       if(lat.le.-pi/15.) then
         lmin = lat
	 jmin = j
!       elseif(lat.gt.pi/18.) then
       elseif(lat.gt.pi/15.) then
         lmax = lat
	 jmax = j
	 goto 300
       endif
      enddo   
 300  r0sq = r0*r0*1.e-4 ! in meter   
      facmax = dt0dts(jmax)/r0sq*rcos0s(jmax)
      do j= 1,nmlat
        fac = dt0dts(j)/r0sq*rcos0s(j)
	lat = ylatm(j)
        do i = 1,nmlonp1
!
! Linear interpolation of J_mr between -10 and 10 degrees
!
	  if(j.gt.jmin.and.j.lt.jmax) then
	    pol = nscrrt(i,jmin)-nscrrt(i,jmax)*facmax
	    pol = pol / (lmin-lmax)
	    nscrrt(i,j)= nscrrt(i,jmin) + pol*(lat-lmin)
	  else
            nscrrt(i,j) = nscrrt(i,j)*fac
	  endif	
	enddo 		! endo i-loop
!
! Simple smoothing in longitudinal direction (weighted average)
! (nonsmooth values due to interpolation)
!
        do i = 1,nmlonp1	
	  if (i.eq.1) then 
	     tout(i,j,1)= (nscrrt(i+1,j)+ 
     |  	    3.*nscrrt(i,j)+ nscrrt(nmlonp1-1,j))/5.0
	  elseif (i.eq.nmlonp1) then 
	     tout(i,j,1)= (nscrrt(2,j)+ 
     |  	     3.*nscrrt(i,j)+nscrrt(i-1,j))/5.0
	  else
	     tout(i,j,1)= (nscrrt(i+1,j)+ 
     |  	     3.*nscrrt(i,j)+nscrrt(i-1,j))/5.0
	  endif 
	  nscrrt(i,j) = tout(i,j,1) 
	  tout(i,j,:) = nscrrt(i,j) ! copy for secondary history field
        enddo 		! endo i-loop
      enddo 		! enddo j-loop
!	
!      call addfld('JQR',' j_qr (magnetic)','[A/m^2]',
!    |   nscrrt(:,:),'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!
      end subroutine nosocrrt
!-----------------------------------------------------------------------
! BOP
! !IROUTINE: nosocrdens
! !INTERFACE:
      subroutine nosocrdens
!
! !DESCRIPTION:
! am_02/02: calculate current density J_e1 (half level)
!  K_(q,phi) & K_(q,lam) (full level)
!
      use params_module,only: nmlon,nmlonp1,nmlat,nmlath,nlev,nlevp1,
     |  nmlev
      use cons_module,only: dlonm,re,r0,ylatm
      use dynamo_module,only: ed23d,ed13d,bmodm3d,sinim3d,sigma2m3d,
     |  adotv2m3d,a1a2m3d,sigma1m3d,adotam3d,adotv1m3d,zpotenm3d,
     |  nscrrt,
     |  je13d,je23d,
     |  je1oD_pg3d                     ! J_e1/D (plasma pressure + gravity)
      use magpres_g_module,only: j_pg  ! flag for mag.pressure and gravity 
      use addfld_module,only: addfld
      implicit none
!
! !PARAMETERS: 
! !RETURN VALUE:
!
! !REVISION HISTORY:
! 05.03.15  <Astrid Maute> <include header> 
! 
! EOP 
!
! Local:
      real :: kqphi_int(nmlonp1,nmlat,-2:nlevp1)
      real :: kqphi3d(nmlonp1,nmlat),dkqphi(nmlonp1,nmlat)
      real :: kqlam(nmlonp1,nmlat)
      real :: fac,lamm,sinlamm,sinim,dh,act,actpk,facq,facsin
      real :: coslamq2,sinlamq,r0m,ed1h,ed2h,facqm,facqp
      real :: fsumn(nmlat),afsumn(nmlat),epsn,difflm
      real :: fsums(nmlat),afsums(nmlat),epss
      real :: fmin,fmax
      real :: adotam23d
      integer :: i,j,jj,k,l,ip1f,ip2f,ip3f
!  
! Calculate current density component Je1 (Richmond: Ionospheric
! Electrodynamics using magnetic apex coordinates pp.203 (eq 5.7))
!   at 1/2 level
!   ar half levels: sig_ped,sig_hall, d_1**2/D, d1*d2/D, ue1, ue2, be3
!   at full levels: ed1, ed2
!
!   je1/d = (sig_ped * d_1**2/D * (ed1 + ue2*be3) + 
!           (sig_ped* d1*d2/D - sig_hall) * (ed2 - ue1*be3)
!   je13d = je1/d 
!   je2/d = (sig_ped * d_2**2/D * (ed2 - ue1*be3) + 
!           (sig_ped* d1*d2/D + sig_hall) * (ed2 + ue1*be3)
!   je23d = je2/d 
!   for j_pg: je1/D=(rho*g-grad p)xB/B_mag^2* d_1/D ! gravity and plasma pressure
!         units A/cm^2  factor of 1e4 to convert from A/cm^2 ->A/m^2
!
      adotam23d = 1.0  ! approximation of d_2^2
      do j = 1,nmlat
        do  k = -2,nlev
          do i = 1,nmlonp1
	    fac = sigma1m3d(i,j,k)*adotam3d(i,j)
            ed1h = 0.5*(ed13d(i,j,k)+ed13d(i,j,k+1))	! ed1 at half level
	    je13d(i,j,k) = fac*(ed1h+
     |  		       adotv2m3d(i,j,k)*bmodm3d(i,j))
	    fac = sigma1m3d(i,j,k)*a1a2m3d(i,j) - sigma2m3d(i,j,k)
            ed2h = 0.5*(ed23d(i,j,k)+ed23d(i,j,k+1))	! ed2 at half level
	    je13d(i,j,k) = je13d(i,j,k) + fac*(ed2h-
     |  		       adotv1m3d(i,j,k)*bmodm3d(i,j))
!           
	    fac = sigma1m3d(i,j,k)*adotam23d
	    je23d(i,j,k) = fac*(ed2h-
     |  		       adotv1m3d(i,j,k)*bmodm3d(i,j))
	    fac = sigma1m3d(i,j,k)*a1a2m3d(i,j) + sigma2m3d(i,j,k)
	    je23d(i,j,k) = je23d(i,j,k) + fac*(ed1h+
     |  		       adotv2m3d(i,j,k)*bmodm3d(i,j))
!           
!           je13d(i,j,k) = 0.  ! only contribution due to pressure/gravity
	    if(j_pg) je13d(i,j,k) = je13d(i,j,k)+1e4*je1oD_pg3d(i,j,k) ! convert cm -> m
	  enddo
	  je13d(nmlonp1,j,k) = je13d(1,j,k)
	  je23d(nmlonp1,j,k) = je23d(1,j,k)
	enddo
	je13d(:,j,nlevp1) = je13d(:,j,nlev)
	je23d(:,j,nlevp1) = je23d(:,j,nlev) 
!        call addfld('JE13D',' J_e1/D (magnetic)','A/m^2',
!    |     je13d(:,j,:),'mlon',1,nmlonp1,'mlev',1,nlevp1,j)
      enddo
!      
!     call addfld('bmodm3d','bmodm3d (magnetic)',' ',
!    |  bmodm3d(:,:),'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!  
! Calculate K_(q,phi) (Richmond: Ionospheric at full level
! Electrodynamics using magnetic apex coordinates pp.208 (eq 7.4))
!   K_(q,phi) = int_(h_l)^(h_u) [( [R_0/R]^0.5 * 
!               * je1/D * sin(lam_q)/sin(lam_m) sin(I_m)/sin(I)*D] dh
!   with F = D*sin(lam_m)/sin(lam_q)*sin(I)/sin(I_m)*[R/R_0]^3
!  
      do j = 1,nmlat
        sinlamq	= sin(ylatm(j))				! sin(lam_q)
	if(j.eq.nmlath) sinlamq	= sin(ylatm(j-1))
        coslamq2  = 1. - sinlamq*sinlamq		! cos^2(lam_q)
        do i = 1,nmlonp1
!
! At equator sin lam_q/sin I is set to the average otherwise quotient = 0
!   check this later 010611
!
	  if(j.eq.nmlath) then
	    facsin = sinlamq/sinim3d(i,j-1)	! sin(lam_q)/sin(I)		
	    k = -2
            fac	  = r0/(re + zpotenm3d(i,j,k))	! R_0 / R
	    lamm  = acos(sqrt(coslamq2*fac))	! cos^2(lam_m) = R_0/R*cos^2(lam_q)
	    lamm  = sign(lamm,ylatm(j-1))	! lam_m
	    sinlamm = sin(lamm)			! sin(lam_m)
!
! sin(I_m) = sin(lam_m)/sqrt(1/4+3/4*sin^2(lam_m))
	    sinim = sinlamm/sqrt(.25+.75*sinlamm**2)	
	    facq  = sinim/sinlamm*facsin	! sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
            fac   = fac**0.5			! sqrt(R_0 / R)
!
! sqrt(R_0 / R)*sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
            act   = fac*facq			
	  else
	    facsin = sinlamq/sinim3d(i,j)	! sin(lam_q)/sin(I)	
	    k = -2
            fac	  = r0/(re + zpotenm3d(i,j,k))	! R_0 / R
	    lamm  = acos(sqrt(coslamq2*fac))	! cos^2(lam_m) = R_0/R*cos^2(lam_q)
	    lamm  = sign(lamm,ylatm(j))		! lam_m
	    sinlamm = sin(lamm)			! sin(lam_m)
!
! sin(I_m) = sin(lam_m)/sqrt(1/4+3/4*sin^2(lam_m))
	    sinim = sinlamm/sqrt(.25+.75*sinlamm**2)	
	    facq  = sinim/sinlamm*facsin	! sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
            fac   = fac**0.5			! sqrt(R_0 / R)
!
! sqrt(R_0 / R)*sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
            act   = fac*facq			
	  endif
!
	  kqphi3d(i,j)     = 0.0 
          kqphi_int(i,j,:) = 0.0   
	  do k = -2,nlev  			! height integration from k=-2 to nlev
! (R_e + h)/R_0  : R_0 = R_e + h0 ; h0 = 9.e6 cm
            fac  = r0/(re + zpotenm3d(i,j,k+1))	! R_0 / R
	    lamm = acos(sqrt(coslamq2*fac))	! cos^2(lam_m) = R_0/R*cos^2(lam_q)
	    lamm = sign(lamm,ylatm(j))		! lam_m
	    sinlamm = sin(lamm)			! sin(lam_m)
!
! sin(I_m) = sin(lam_m)/sqrt(1/4+3/4*sin^2(lam_m))
	    sinim = sinlamm/sqrt(.25+.75*sinlamm**2)	
	    facq  = sinim/sinlamm*facsin	! sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
	    dh    = zpotenm3d(i,j,k+1) - zpotenm3d(i,j,k)	! dh
	    dh    = dh*1.e-2			! convertion [cm] to [m]
            fac   = fac**0.5			! sqrt(R_0 / R)
!
! sqrt(R_0 / R)*sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
            actpk = fac*facq			
!
! integration value at 1/2 level
            kqphi3d(i,j) = kqphi3d(i,j) + 0.5*(actpk+act)*	
     |                     je13d(i,j,k)*dh
            kqphi_int(i,j,k+1) = kqphi3d(i,j)
	    act = actpk
          enddo
! lowest level integration value is zero since this is below the 
! ionospheric current layer 110 km it doesn't matter
	  kqphi_int(i,j,0)  = 0.  
        enddo
!        if(j_pg) call addfld('KQPHI_TOT','kqphi_int (u,E,g,p)',
!    |    '[A/m]',kqphi_int(:,j,:),'mlon',1,nmlonp1,'imlev',1,nlevp1,j)
!        if(.not.j_pg) call addfld('KQPHI_UE','kqphi_int (u,E)',
!    |    '[A/m]',kqphi_int(:,j,:),'mlon',1,nmlonp1,'imlev',1,nlevp1,j)
!        call addfld('zpotm3d','zpotm3d','cm',
!    |     zpotenm3d(:,j,:),'mlon',1,nmlonp1,'imlev',1,nlevp1,j)
      enddo
! 
! save on secondary history 
!      call addfld('KQPHI',' kqphi (magnetic)','[A/m]',
!    |    kqphi3d(:,:),'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!                
!  
! Calculate K_(q,lam) (Richmond: Ionospheric at full level
! Electrodynamics using magnetic apex coordinates pp.208 (eq 7.5))
!   K_(q,lam) = -1/cos(lam_q) int_(-pi/2)^(lam_q) [J_mr*R*cos(lam_q) + 
!                d(K_(q,phi))/d(phi_q) ] d(lam_q)
!  
!     d(K_(q,phi))/d(phi_q)
!
      fac = 0.5/dlonm		! 1/(2*d lon_m)
      do j = 1,nmlat 
        do i = 2,nmlonp1-1
!
! central difference
! (kqphi3d(i+1/2)-kqphi3d(i-1/2))/2
!
          dkqphi(i,j) = (kqphi3d(i+1,j)-kqphi3d(i-1,j))*fac 
        enddo				
	dkqphi(1,j)      = (kqphi3d(2,j)-kqphi3d(nmlonp1-1,j))*fac
	dkqphi(nmlonp1,j) = dkqphi(1,j) 
      enddo
      r0m = r0*1.e-2
      do i = 1,nmlonp1
       fsums(1)   = 0.0 
       afsums(1)  = 0.0 
       kqlam(i,1)= -dkqphi(i,1) ! at south pole K_(q,lam) = -d(K_(q,phi))/d(phi)
!
! [J_mr*R*cos(lam_q)+d K_(qphi)]_(j=1)
       act = nscrrt(i,1)*r0m*cos(ylatm(1))+dkqphi(i,1)	
! 
       do j = 2,nmlath
!
! difflm: d | (lam_q(j)-lam_q(j-1))|
! actpk: [J_mr*R*cos(lam_q)+d K_(qphi)]_(j)
! act:  [integrand]_(j-1)+[integrand]_(j)
! act: -[integrand]_(j-1/2)] d lam_q
! fsums: int_(-pi/2)^(lam_q)[integrand] d lam_q
! afsums: int_(-pi/2)^(lam_q) | [integrand] | d lam_q
!
         difflm   = abs(ylatm(j)-ylatm(j-1))	
	 actpk    = nscrrt(i,j)*r0m*cos(ylatm(j))+dkqphi(i,j)	
	 act	  = act +actpk		    
	 act	  = -act/2.0*difflm	    
	 fsums(j) = fsums(j-1)  + act	    
	 afsums(j)= afsums(j-1) + abs(act)  
	 act = actpk
       enddo
! 
! Integrate from the north pole to equator
       j=nmlat
       kqlam(i,nmlat)= dkqphi(i,nmlat) ! at north pole K_(q,lam) = d(K_(q,phi)
       fsumn(j)  = 0.0
       afsumn(j) = 0.0 
!       
! act: [J_mr*R*cos(lam_q)+d K_(qphi)]_(nmlat)
! difflm: d | (lam_q(j)-lam_q(j-1))|
! actpk: [J_mr*R*cos(lam_q)+d K_(qphi)]_(j)
! act: [integrand]_(j+1)+[integrand]_(j)
! act: [integrand]_(j+1/2)] d lam_q
! fsumn: int_(pi/2)^(lam_q)[integrand] d lam_q
! afsumn: int_(pi/2)^(lam_q) | [integrand] | d lam_q 
!       
       act = nscrrt(i,j)*r0m*cos(ylatm(j))+dkqphi(i,j)	
       do j = nmlat-1,nmlath,-1
         difflm    = abs(ylatm(j)-ylatm(j+1))		
	 actpk    = nscrrt(i,j)*r0m*cos(ylatm(j))+dkqphi(i,j)	
	 act	  = act +actpk				
	 act	  = act/2.0*difflm			
	 fsumn(j) = fsumn(j+1)  + act			
	 afsumn(j)= afsumn(j+1) + abs(act)		
	 act = actpk
       enddo
! 
! correction to equal integration from north to equator with
!                     integration from south to equator
! epsn: half of error to south weighted by absolute value
! epss: half of error to north weighted by absolute value
! kqlam_cor = kqlam - err/2/abs(kqlam)_south*abs(kqlam)_(lam_q)
! kqlam_cor /cos(lam_q) 
! 
       epsn = 0.5*(fsums(nmlath)-fsumn(nmlath))/afsumn(nmlath)	
       epss = 0.5*(fsums(nmlath)-fsumn(nmlath))/afsums(nmlath)	
!       
       do j = 2,nmlath			! correct and copy into kqlam
	kqlam(i,j)= 0.0
	kqlam(i,j)= fsums(j) - epss*afsums(j)	
	kqlam(i,j)= kqlam(i,j)/cos(ylatm(j))	
       enddo
!       
! kqlam_cor = kqlam + err/2/abs(kqlam)_north*abs(kqlam)_(lam_q)
! kqlam_cor /cos(lam_q) 
!
       do j = nmlat-1,nmlath,-1
	kqlam(i,j)= 0.0
	kqlam(i,j)= fsumn(j) + epsn*afsumn(j)	
	kqlam(i,j)= kqlam(i,j)/cos(ylatm(j))	
       enddo
      enddo ! end of i-loop
!      
!      call addfld('KQLAM',' kqlam (magnetic)','[A/m]',
!    |      kqlam(:,:),'mlon',1,nmlonp1,'mlat',1,nmlat,0)
!     
      end subroutine nosocrdens
!-----------------------------------------------------------------------

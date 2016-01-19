module current
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate plasma current diagnostics. 
!
! Subroutines in this module is called only when namelist parameter 
! current_kq > 0. Calling these routines will seriously degrade performance
! of the model as a whole.
!
! Subs noso_coef and noso_crrt are serial, called by the root task only.
! Sub noso_crdens is parallel on magnetic subdomains, called by all tasks.
!
! Sub noso_coef is called by stub sub nosocoef (at end of this file),
!   which is called from sub complete_integrals (pdynamo.F) by the
!   root task only. This is necessary to avoid circular dependency.
! Sub noso_crrt is called from advance by the root task after pdynamo.
! Sub noso_crdens is called from advance by all tasks.
!
  use params_module,only: nmlon,nmlonp1,nmlat,nmlevp1,nmlath
!
! Global integrated conductivities are available only on the root task:
!
  use pdynamo_module,only: & 
    zigm11_glb,zigm22_glb,zigmc_glb,zigm2_glb,rhs_glb,rim_glb,nmlon0
  use addfld_module,only: addfld
  use diags_module,only: mkdiag_JE13D, mkdiag_JE23D, mkdiag_KQLAM, &
    mkdiag_KQPHI, mkdiag_JQR
  implicit none
  real,parameter :: unitv(nmlon)=1.
  real,dimension(nmlon0,nmlat,10) :: nscoef ! nmlon0==nmlonp1
  real,dimension(nmlon0,nmlat) :: nscrrt

contains
!-----------------------------------------------------------------------
subroutine noso_coef
!
! am_02/02: calculate the coefficient stencil for both hemisphere which
!   is used to calculate the height integrated current densities 
!   K_(q,phi)K_(q,lam)
!      
  use cons_module,only: dlatm,dlonm,pi,r0
!
! Local:
!
  integer :: i,j,jnpole,jspole,n,jn,js,je
  real :: dfac,dfac1n,dfac1s,dfac2n,dfac2s
  real,dimension(nmlonp1,nmlat) :: nszigmc,nszigm2,nszigm11,nszigm22,nsrhs
  real,dimension(0:nmlon0+1,nmlat) :: array
  real :: cs(nmlat)
!
! Externals:
  real,external :: sddot ! in util.F
!
! Init:
!
  nszigmc  = 0.
  nszigm2  = 0.
  nszigm11 = 0.
  nszigm22 = 0.
  nsrhs    = 0.
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
      nszigmc(i,j)  = -zigmc_glb(i,j)
      nszigmc(i,j)  = (nszigmc(i,j)+zigm2_glb(i,j))/dfac
      nszigm2(i,j)  = nszigmc(i,j)-2.*zigm2_glb(i,j)/dfac
      nszigm22(i,j) = zigm22_glb(i,j)*dfac1n
      nszigm11(i,j) = zigm11_glb(i,j)/dfac2n
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

  jnpole = nmlat
  jspole = 1
  dfac1n = cs(jnpole)/dlatm**2
  dfac1s = cs(jspole)/dlatm**2   ! is not necessary cos symmetric
  do i = 1,nmlonp1
    nszigmc (i,jnpole) = -zigmc_glb(i,jnpole)
    nszigmc (i,jnpole) = (nszigmc(i,jnpole)+zigm2_glb(i,jnpole))/dfac
    nszigm2 (i,jnpole) = nszigmc(i,jnpole)-2.*zigm2_glb(i,jnpole)/dfac
    nszigm22(i,jnpole) = zigm22_glb(i,jnpole)*dfac1n

    nszigmc (i,jspole) = -zigmc_glb(i,jspole)
    nszigmc (i,jspole) = (nszigmc(i,jspole)+zigm2_glb(i,jspole))/dfac
    nszigm2 (i,jspole) = nszigmc(i,jspole)-2.*zigm2_glb(i,jspole)/dfac
    nszigm22(i,jspole) = zigm22_glb(i,jspole)*dfac1s
!
! Set zigm11 to zero at the magnetic poles (1 and 97) to avoid floating
!   point exception (values at the poles are not used)
!
    nszigm11(i,jnpole) = 0.0
    nszigm11(i,jspole) = 0.0
  enddo
!
! Clear array for difference stencil over north and south hemisphere
!
  nscoef = 0.0
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
! Calculate right hand side of pde from rim_glb(1) and rim_glb(2)
!
  do js = 2,nmlath-1  ! 2,48  south pole-1 to equator-1
    jn = js+nmlath-1   ! 50,96 equator-1 to north pole-1
!
! Differentiate rim(1) w.r.t lambda
!
    do i = 2,nmlon-1
      nsrhs(i,js) = 1.0/(dlonm*cs(js))*0.5*(rim_glb(i+1,js,1)-rim_glb(i-1,js,1))
      nsrhs(i,jn) = 1.0/(dlonm*cs(jn))*0.5*(rim_glb(i+1,jn,1)-rim_glb(i-1,jn,1))
    enddo
!
! Values at the poles
!
    nsrhs(1,js) = 1.0/(dlonm*cs(js))*0.5*(rim_glb(2,js,1)-rim_glb(nmlon,js,1))
    nsrhs(1,jn) = 1.0/(dlonm*cs(jn))*0.5*(rim_glb(2,jn,1)-rim_glb(nmlon,jn,1))
    nsrhs(nmlon,js) = 1.0/(dlonm*cs(js))*0.5*(rim_glb(1,js,1)-rim_glb(nmlon-1,js,1))
    nsrhs(nmlon,jn) = 1.0/(dlonm*cs(jn))*0.5*(rim_glb(1,jn,1)-rim_glb(nmlon-1,jn,1))
  enddo
!
! Differentiate rim(2) w.r.t theta0
!
  do js = 2,nmlath-1  ! 2,48  south pole -1 to equator-1
    jn = js+nmlath-1  ! 50,96 equator+1 to north pole-1
    do i = 1,nmlon
      nsrhs(i,js) = nsrhs(i,js) - 1.0/(dlatm*cs(js))*0.5* &
                    (rim_glb(i,js+1,2)*cs(js+1)-rim_glb(i,js-1,2)*cs(js-1))
      nsrhs(i,jn) = nsrhs(i,jn) + 1.0/(dlatm*cs(jn))*0.5* &
                    (rim_glb(i,jn+1,2)*cs(jn+1)-rim_glb(i,jn-1,2)*cs(jn-1))
   enddo
  enddo
!
! Calculate value at the poles by averaging over i:nmlon
!
  nsrhs(1,nmlat) = -2./float(nmlon)*sddot(nmlon,unitv,rim_glb(1,nmlat-1,2))/cs(nmlat-1)
  nsrhs(1,1)     = -2./float(nmlon)*sddot(nmlon,unitv,rim_glb(1,2,2))/cs(2)
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
  nsrhs(i,je) = 0.5/dlonm*(rim_glb(i+1,je,1)-rim_glb(nmlon,je,1))
  nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)* &
                    rim_glb(i,je,2)+ cs(je+1)*rim_glb(i,je+1,2))
  nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)* &
                    rim_glb(i,je,2)+ cs(je-1)*rim_glb(i,je-1,2))
  nsrhs(i,je) = 0.5*nsrhs(i,je)
!
  do i = 2,nmlon-1
    nsrhs(i,je) = 0.5/dlonm*(rim_glb(i+1,je,1)-rim_glb(i-1,je,1))
    nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)* &
                      rim_glb(i,je,2)+ cs(je+1)*rim_glb(i,je+1,2))
    nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)* &
                      rim_glb(i,je,2)+ cs(je-1)*rim_glb(i,je-1,2))
    nsrhs(i,je) = 0.5*nsrhs(i,je)
  enddo
!
  i = nmlon
  nsrhs(i,je) = 0.5/dlonm*(rim_glb(1,je,1)-rim_glb(i-1,je,1))
  nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)* &
                    rim_glb(i,je,2)+ cs(je+1)*rim_glb(i,je+1,2))
  nsrhs(i,je) = nsrhs(i,je) + 1./dlatm*(cs(je)* &
                    rim_glb(i,je,2)+ cs(je-1)*rim_glb(i,je-1,2))
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
  nscoef(:,1,10)     = 0.5

end subroutine noso_coef
!-----------------------------------------------------------------------
subroutine nsstencil(zigm,nmlon0,nmlat,nscoef,array,ncoef,nmlath)
!
! Args:
  integer,intent(in):: nmlon0,nmlat,ncoef,nmlath
  real,intent(in)   :: zigm(nmlon0,nmlat)
  real,intent(out)  :: nscoef(nmlon0,nmlat,10),array(0:nmlon0+1,nmlat)
!
! Local:
  integer :: i,j
!
! Copy zigm (south to north pole) into array (south to north pole)
!   check: is this change necessary?
!
  do j = 1,nmlat    ! 1,97
    do i = 1,nmlon0
      array(i,j) = zigm(i,j)
    enddo
  enddo
!
! Extend over 32 grid spaces to alow a total of five grid levels
!   check: is this necessary because no solving with these coefficients?
!
  i = 1
  do j = 1,nmlat
    array(i-1,j)    = array(nmlon0-i,j)
    array(nmlon0+i,j) = array(1+i,j)
  enddo
!
! Calculate contribution to stencil for each grid point
!
  call nscnm(array,nmlon0,nmlat,nscoef,ncoef,nmlath)

end subroutine nsstencil
!-----------------------------------------------------------------------
subroutine nscnm(array,nmlon0,nmlat,nsc,ncoef,nmlath)
!
!      calculate contribution for each zigm
!      the diagonial dominance of the stencil is not checked
!      since the coefficients are not used for solving
!      one reason for difference between added north-south
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
! Args:
  integer,intent(in) :: nmlon0,nmlat,ncoef,nmlath
  real,intent(in)    :: array(0:nmlon0+1,nmlat)
  real,intent(out)   :: nsc(nmlon0,nmlat,10)
!
! Local:
  integer :: i,j,js,jn 

!
! Check: what's about boundary equator values for nsc, why not used?
! Calculate contribution for zigm11
!
  if (ncoef==1) then
    do j = 2,nmlath-1   ! 2,48  not value at the pols (jn=97 js=1)
      jn  = nmlath+j-1      ! 50,96 northern hemisphere
      js  = nmlath-j+1      ! 48,2  southern hemisphere
!
      do i = 1,nmlon0
        nsc(i,jn,1) = nsc(i,jn,1)+.5*(array(i,jn)+array(i+1,jn))
        nsc(i,jn,5) = nsc(i,jn,5)+.5*(array(i,jn)+array(i-1,jn))
        nsc(i,jn,9) = nsc(i,jn,9)-.5*(array(i+1,jn)+2.*array(i,jn)+array(i-1,jn))

        nsc(i,js,1) = nsc(i,js,1)+.5*(array(i,js)+array(i+1,js))
        nsc(i,js,5) = nsc(i,js,5)+.5*(array(i,js)+array(i-1,js))
        nsc(i,js,9) = nsc(i,js,9)-.5*(array(i+1,js)+2.*array(i,js)+array(i-1,js))
      enddo
    enddo
!
! am 2001-6-27 include boundary condition at equator
!
    do i = 1,nmlon0
      nsc(i,nmlath,1) = nsc(i,nmlath,1)+.5*(array(i,nmlath)+array(i+1,nmlath))
      nsc(i,nmlath,5) = nsc(i,nmlath,5)+.5*(array(i,nmlath)+array(i-1,nmlath))
      nsc(i,nmlath,9) = nsc(i,nmlath,9)-.5*(array(i+1,nmlath)+ &
        2.*array(i,nmlath)+array(i-1,nmlath))
    enddo
!
! Calculate contribution for zigm12 (=ZIGMC+ZIGM2)
!
  elseif (ncoef==2) then
    do j = 2,nmlath-1      ! 2,48  not value at the pols (jn=97 js=1)
      jn  = nmlath+j-1     ! 50,96 northern hemisphere
      js = nmlath-j+1      ! 48,2  southern hemisphere
!
      do i = 1,nmlon0
        nsc(i,jn,2) = nsc(i,jn,2)+.5*(array(i,jn)+array(i+1,jn))
        nsc(i,jn,4) = nsc(i,jn,4)-.5*(array(i,jn)+array(i-1,jn))
        nsc(i,jn,6) = nsc(i,jn,6)+.5*(array(i,jn)+array(i-1,jn))
        nsc(i,jn,8) = nsc(i,jn,8)-.5*(array(i,jn)+array(i+1,jn))
        nsc(i,jn,3) = nsc(i,jn,3)+.5*(-array(i-1,jn)+array(i+1,jn))
        nsc(i,jn,7) = nsc(i,jn,7)-.5*(-array(i-1,jn)+array(i+1,jn))

        nsc(i,js,2)=nsc(i,js,2)+.5*(array(i,js)+array(i+1,js))
        nsc(i,js,4)=nsc(i,js,4)-.5*(array(i,js)+array(i-1,js))
        nsc(i,js,6)=nsc(i,js,6)+.5*(array(i,js)+array(i-1,js))
        nsc(i,js,8)=nsc(i,js,8)-.5*(array(i,js)+array(i+1,js))
        nsc(i,js,3)=nsc(i,js,3)+.5*(-array(i-1,js)+array(i+1,js))
        nsc(i,js,7)=nsc(i,js,7)-.5*(-array(i-1,js)+array(i+1,js))
       enddo
    enddo
  elseif (ncoef==3) then
!
! Calculate contribution for zigm21 (=ZIGMC-ZIGM2)
!
    do j = 3,nmlath-1   ! 3,48  not value at the pols (jn=97 js=1)
      jn  = nmlath+j-1      ! 51,96 northern hemisphere
      js = nmlath-j+1      ! 47,2  southern hemisphere
!
      do i = 1,nmlon0
        nsc(i,jn,2) = nsc(i,jn,2)+.5*(array(i,jn)+array(i,jn+1))
        nsc(i,jn,4) = nsc(i,jn,4)-.5*(array(i,jn)+array(i,jn+1))
        nsc(i,jn,6) = nsc(i,jn,6)+.5*(array(i,jn)+array(i,jn-1))
        nsc(i,jn,8) = nsc(i,jn,8)-.5*(array(i,jn)+array(i,jn-1))
        nsc(i,jn,1) = nsc(i,jn,1)+.5*(array(i,jn+1)-array(i,jn-1))
        nsc(i,jn,5) = nsc(i,jn,5)-.5*(array(i,jn+1)-array(i,jn-1))

        nsc(i,js,8)=nsc(i,js,8)-.5*(array(i,js)+array(i,js+1))
        nsc(i,js,6)=nsc(i,js,6)+.5*(array(i,js)+array(i,js+1))
        nsc(i,js,4)=nsc(i,js,4)-.5*(array(i,js)+array(i,js-1))
        nsc(i,js,2)=nsc(i,js,2)+.5*(array(i,js)+array(i,js-1))
        nsc(i,js,1)=nsc(i,js,1)-.5*(array(i,js+1)-array(i,js-1))
        nsc(i,js,5)=nsc(i,js,5)+.5*(array(i,js+1)-array(i,js-1))
      enddo
    enddo
!
    j = 2   ! 2 change in sign for equatorial values
    jn = nmlath+j-1      ! 50 northern hemisphere
    js = nmlath-j+1      ! 48  southern hemisphere
!
   do i = 1,nmlon0   !am2001-7-3
     nsc(i,jn,2) = nsc(i,jn,2)+.5*(array(i,jn)+array(i,jn+1))
     nsc(i,jn,4) = nsc(i,jn,4)-.5*(array(i,jn)+array(i,jn+1))
     nsc(i,jn,6) = nsc(i,jn,6)+.5*(array(i,jn)-array(i,jn-1))
     nsc(i,jn,8) = nsc(i,jn,8)-.5*(array(i,jn)-array(i,jn-1))
     nsc(i,jn,1) = nsc(i,jn,1)+.5*(array(i,jn+1)+array(i,jn-1))
     nsc(i,jn,5) = nsc(i,jn,5)-.5*(array(i,jn+1)+array(i,jn-1))
!
     nsc(i,js,8)=nsc(i,js,8)-.5*(array(i,js)-array(i,js+1))
     nsc(i,js,6)=nsc(i,js,6)+.5*(array(i,js)-array(i,js+1))
     nsc(i,js,4)=nsc(i,js,4)-.5*(array(i,js)+array(i,js-1))
     nsc(i,js,2)=nsc(i,js,2)+.5*(array(i,js)+array(i,js-1))
     nsc(i,js,1)=nsc(i,js,1)-.5*(-array(i,js+1)-array(i,js-1))
     nsc(i,js,5)=nsc(i,js,5)+.5*(-array(i,js+1)-array(i,js-1))
   enddo
!
! Low latitude boundary conditions: contribution from zigm21(i,j-1/2)=0
!
    j = 1
    jn = nmlath+j-1      ! 49
    js = nmlath-j+1      ! 49
!
    do i = 1,nmlon0
      nsc(i,jn,2) = nsc(i,jn,2)+.25*(-array(i,jn)+array(i,jn+1))
      nsc(i,jn,4) = nsc(i,jn,4)-.25*(-array(i,jn)+array(i,jn+1))
      nsc(i,jn,1) = nsc(i,jn,1)+.25*(-array(i,jn)+array(i,jn+1))
      nsc(i,jn,5) = nsc(i,jn,5)-.25*(-array(i,jn)+array(i,jn+1))

      nsc(i,jn,2) = nsc(i,jn,2)+.25*(-array(i,js)+array(i,js-1))
      nsc(i,jn,4) = nsc(i,jn,4)-.25*(-array(i,js)+array(i,js-1))
      nsc(i,jn,1) = nsc(i,jn,1)-.25*(array(i,js)-array(i,js-1))
      nsc(i,jn,5) = nsc(i,jn,5)+.25*(array(i,js)-array(i,js-1))
    enddo
!
! Calculate contribution for zigm22
!
  elseif (ncoef==4) then
    do j = 2,nmlath-1   ! 2,48  not value at the pols (jn=97 js=1)
      jn  = nmlath+j-1      ! 50,96 northern hemisphere
      js = nmlath-j+1      ! 48,2  southern hemisphere
!
      do i = 1,nmlon0
        nsc(i,jn,3) = nsc(i,jn,3)+.5*(array(i,jn)+array(i,jn+1))
        nsc(i,jn,7) = nsc(i,jn,7)+.5*(array(i,jn)+array(i,jn-1))
        nsc(i,jn,9) = nsc(i,jn,9)-.5*(array(i,jn-1)+2.0*array(i,jn)+array(i,jn+1))
        nsc(i,js,7)=nsc(i,js,7)+.5*(array(i,js)+array(i,js+1))
        nsc(i,js,3)=nsc(i,js,3)+.5*(array(i,js)+array(i,js-1))
        nsc(i,js,9)=nsc(i,js,9)-.5*(array(i,js-1)+2.0*array(i,js)+ array(i,js+1))
       enddo
    enddo
!
! Low latitude boundary conditions: contribution from zigm22(i,j-1/2)=0
!
    j = 1
    jn  = nmlath+j-1      ! 49
    js = nmlath-j+1      ! 49
    do i = 1,nmlon0
       nsc(i,jn,3) = nsc(i,jn,3)+.25*(array(i,jn)+array(i,jn+1))
       nsc(i,jn,9) = nsc(i,jn,9)-.25*(array(i,jn)+array(i,jn+1))
!  am 2001-7-02 otherwise double coefficients
       nsc(i,jn,3)=nsc(i,jn,3)+.25*(array(i,js)+array(i,js-1))
       nsc(i,jn,9)=nsc(i,jn,9)-.25*(array(i,js-1)+array(i,js))
    enddo
  endif ! ncoef == 1,2,3 or 4
end subroutine nscnm
!-----------------------------------------------------------------------
subroutine noso_crrt
  use cons_module,only: dt0dts,rcos0s,pi,r0,re,ylatm
  use mpi_module,only: mlon0,mlon1,mlat0,mlat1,mp_allgather_2d
  use pdynamo_module,only: phim=>phim_glb ! (nmlonp1,nmlat) on root task only
!
! Calculate current for both hemisphere:
! [stencil*potential -RHS] = R**2 * J_mr / dt0dts / rcos0s
!
! Local:
  integer :: j,i,k,n,jmod,jmin,jmax,jn,js
  real :: vtmp,r0sq,fac,facmax,lat,lmin,lmax,pol
  real,dimension(nmlonp1,nmlat,nmlevp1) :: tout
!
! External:
  real,external :: sddot ! in util.F
!
  nscrrt = 0.0

  i = 1
  j = nmlath     ! equator set J_mr = 0
  nscrrt(i,j) = 0.0
!
  do j = 2,nmlath-1  ! 2,48  no poles
    jn = nmlath+j-1  ! 50,96 nhem
    js = nmlath-j+1  ! 48,2  shem
!                                    contribution of stencil
! Northern hemisphere
!
    nscrrt(i,jn) =                nscoef(i,jn,1)*phim(i+1,jn)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,2)*phim(i+1,jn+1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,3)*phim(i  ,jn+1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,4)*phim(nmlon0-1,jn+1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,5)*phim(nmlon0-1,jn)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,6)*phim(nmlon0-1,jn-1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,7)*phim(i  ,jn-1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,8)*phim(i+1,jn-1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,9)*phim(i  ,jn)
    nscrrt(i,jn) = nscrrt(i,jn) - nscoef(i,jn,10)
!
! Southern hemisphere
!
    nscrrt(i,js)=              nscoef(i,js,1)*phim(i+1,js)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,2)*phim(i+1,js-1)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,3)*phim(i  ,js-1)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,4)*phim(nmlon0-1,js-1)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,5)*phim(nmlon0-1,js)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,6)*phim(nmlon0-1,js+1)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,7)*phim(i  ,js+1)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,8)*phim(i+1,js+1)
    nscrrt(i,js)=nscrrt(i,js)+nscoef(i,js,9)*phim(i  ,js)
    nscrrt(i,js)=nscrrt(i,js)-nscoef(i,js,10)
  enddo ! j=2,nmlath-1

  do i = 2,nmlon0-1        ! 2,80
    nscrrt(i,nmlath) = 0.0 ! equator
    do j = 2,nmlath-1      ! 2,48  no poles
      jn = nmlath+j-1      ! 50,96 nhem
      js = nmlath-j+1      ! 48,2  shem
!                                      contribution of stencil
! Northern hemisphere
!
      nscrrt(i,jn) =                nscoef(i,jn,1)*phim(i+1,jn)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,2)*phim(i+1,jn+1)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,3)*phim(i  ,jn+1)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,4)*phim(i-1,jn+1)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,5)*phim(i-1,jn)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,6)*phim(i-1,jn-1)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,7)*phim(i  ,jn-1)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,8)*phim(i+1,jn-1)
      nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,9)*phim(i  ,jn)
      nscrrt(i,jn) = nscrrt(i,jn) - nscoef(i,jn,10)
!
! Southern hemisphere
!
      nscrrt(i,js) =                nscoef(i,js,1)*phim(i+1,js)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,2)*phim(i+1,js-1)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,3)*phim(i  ,js-1)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,4)*phim(i-1,js-1)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,5)*phim(i-1,js)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,6)*phim(i-1,js+1)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,7)*phim(i  ,js+1)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,8)*phim(i+1,js+1)
      nscrrt(i,js) = nscrrt(i,js) + nscoef(i,js,9)*phim(i  ,js)
      nscrrt(i,js) = nscrrt(i,js) - nscoef(i,js,10)
    enddo ! j=2,nmlath-1
  enddo ! i=2,nmlon0-1
!
! For i=nmlonp1 (81)
!
  nscrrt(nmlonp1,nmlath) = 0.0 ! equator
  i = nmlonp1
  do j = 2,nmlath-1   ! 2,48  no poles
    jn  = nmlath+j-1  ! 50,96 nhem
    js = nmlath-j+1   ! 48,2  shem
!                                    contribution of stencil
! Northern hemisphere
!
    nscrrt(i,jn) =                nscoef(i,jn,1)*phim(2,jn)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,2)*phim(2,jn+1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,3)*phim(i  ,jn+1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,4)*phim(i-1,jn+1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,5)*phim(i-1,jn)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,6)*phim(i-1,jn-1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,7)*phim(i  ,jn-1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,8)*phim(2,jn-1)
    nscrrt(i,jn) = nscrrt(i,jn) + nscoef(i,jn,9)*phim(i  ,jn)
    nscrrt(i,jn) = nscrrt(i,jn) - nscoef(i,jn,10)
!
! Southern hemisphere
!
    nscrrt(i,js) =               nscoef(i,js,1)*phim(2,js)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,2)*phim(2,js-1)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,3)*phim(i  ,js-1)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,4)*phim(i-1,js-1)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,5)*phim(i-1,js)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,6)*phim(i-1,js+1)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,7)*phim(i  ,js+1)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,8)*phim(2,js+1)
    nscrrt(i,js) = nscrrt(i,js) +nscoef(i,js,9)*phim(i  ,js)
    nscrrt(i,js) = nscrrt(i,js) -nscoef(i,js,10)
  enddo ! j=2,nmlath-1
!
! Poles
!
  nscrrt(1,1)    = sddot(nmlon,unitv,nscrrt(1,2))/float(nmlon)
  nscrrt(1,nmlat)= sddot(nmlon,unitv,nscrrt(1,nmlat-1))/float(nmlon)
  nscrrt(:,1)    = nscrrt(1,1)             ! extend in longitude
  nscrrt(:,nmlat)= nscrrt(1,nmlat)

  do j = 1,nmlat
    lat = ylatm(j)
    if(lat <= -pi/18.) then
      lmin = lat
      jmin = j
    elseif(lat.gt.pi/18.) then
      lmax = lat
      jmax = j
      exit
    endif
  enddo
  r0sq = r0*r0*1.e-4 ! in meter
  facmax = dt0dts(jmax)/r0sq*rcos0s(jmax)

  do j= 1,nmlat
    fac = dt0dts(j)/r0sq*rcos0s(j)
    lat = ylatm(j)
    do i = 1,nmlonp1
!
! Linear interpolation of J_mr between -10 and 10 degrees
!
      if(j > jmin.and.j < jmax) then
        pol = nscrrt(i,jmin)-nscrrt(i,jmax)*facmax
        pol = pol / (lmin-lmax)
        nscrrt(i,j)= nscrrt(i,jmin) + pol*(lat-lmin)
      else
        nscrrt(i,j) = nscrrt(i,j)*fac
      endif
    enddo           ! endo i-loop
!
! Simple smoothing in longitudinal direction (weighted average)
! (nonsmooth values due to interpolation)
!
    do i = 1,nmlonp1
      if (i == 1) then
         tout(i,j,1)= (nscrrt(i+1,j)+ &
                3.*nscrrt(i,j)+ nscrrt(nmlonp1-1,j))/5.0
      elseif (i == nmlonp1) then
         tout(i,j,1)= (nscrrt(2,j)+ &
                 3.*nscrrt(i,j)+nscrrt(i-1,j))/5.0
      else
         tout(i,j,1)= (nscrrt(i+1,j)+ &
                 3.*nscrrt(i,j)+nscrrt(i-1,j))/5.0
      endif
      nscrrt(i,j) = tout(i,j,1)
      tout(i,j,:) = nscrrt(i,j) ! copy for secondary history field
    enddo           ! endo i-loop
  enddo             ! enddo j-loop
!
! Save JQR: Upward current density (2d)
!
  call mkdiag_JQR('JQR',nscrrt,1,nmlon0,1,nmlat)

end subroutine noso_crrt
!-----------------------------------------------------------------------
subroutine noso_crdens
!
! am_02/02: calculate current density J_e1, K_(q,phi) & K_(q,lam)
! 8/1/13 btf: Parallelize in distributed magnetic grid subdomains.
!
  use mpi_module,only: mlon0,mlon1,mlat0,mlat1,mp_allgather_2d,&
    mp_mag_periodic_f2d,mp_bcast_real
  use pdynamo_module,only: &
    mlev0,mlev1,           & ! -2,nlevp1 at 5-deg, -5,nlevp1 at 2.5-deg
    ped_mag,hal_mag,       & ! 3d (mlon0:mlon1,mlat0:mlat1,mlev0:mlev1)
    ed13d,ed23d,           & ! 3d
    adotv1_mag,adotv2_mag, & ! 3d
    je1pg_mag,             & ! 3d
    zpotm3d                  ! 3d
  use pdynamo_module,only: &
    adota1_mag,            & ! 2d (mlon0:mlon1,mlat0:mlat1)
    be3_mag,               & ! 2d
    a1dta2_mag,            & ! 2d
    sini_mag                 ! 2d
  use input_module,only: current_pg    
  use cons_module,only: ylatm,h0,re,r0,dlonm
!
! Local:
!
  integer :: k,i,j,len
  real :: fac,adotam,ed1h,ed2h,sinlamq,coslamq2,facsin,lamm,sinlamm,act
  real :: sinim,facq,dh,actpk,r0m,difflm,epsn,epss,eps(2)
  real,dimension(nmlat) :: fsums,afsums,fsumn,afsumn
  real,dimension(mlon0:mlon1,mlat0:mlat1,mlev0:mlev1) :: &
    je13d,je23d,  & ! output current density
    kqphi           ! output K_(q,phi)
  real,dimension(mlon0:mlon1,mlat0:mlat1) :: kqphi2d,kqlam
  real,dimension(nmlonp1,nmlat) :: dkqphi_glb,kqphi_glb
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
!   for current_press=1: je1/D=(rho*g-grad p)xB/B_mag^2* d_1/D ! gravity and plasma pressure
!         units A/cm^2  factor of 1e4 to convert from A/cm^2 ->A/m^2
!
! do j=mlat0,mlat1
!   call addfld('JE1PG_MAG',' ',' ', & ! available only when current_pg==1
!     je1pg_mag(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
!   call addfld('PED',' ',' ', &
!     ped_mag(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
!   call addfld('HAL',' ',' ', &
!     hal_mag(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
!   call addfld('ED13D',' ',' ', &
!     ed13d(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
!   call addfld('ED23D',' ',' ', &
!     ed23d(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
!   call addfld('ADOTV1',' ',' ', &
!     adotv1_mag(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
!   call addfld('ADOTV2',' ',' ', &
!     adotv2_mag(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
! enddo

! call addfld('A1DTA2',' ',' ', &
!   a1dta2_mag,'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)
! call addfld('BE3',' ',' ', &
!   be3_mag,   'mlon',mlon0,mlon1,'mlat',mlat0,mlat1,0)

  je13d = 0.
  je23d = 0.
  adotam = 1.0 ! approximation of d_2^2 (local scalar)
  do k=mlev0,mlev1-1
    do j=mlat0,mlat1
      do i=mlon0,mlon1
        fac = ped_mag(i,j,k)*adota1_mag(i,j)
        ed1h = 0.5*(ed13d(i,j,k)+ed13d(i,j,k+1)) ! ed1 at midpoint
        je13d(i,j,k) = fac*(ed1h+adotv2_mag(i,j,k)*be3_mag(i,j))

        fac = ped_mag(i,j,k)*a1dta2_mag(i,j)-hal_mag(i,j,k)
        ed2h = 0.5*(ed23d(i,j,k)+ed23d(i,j,k+1)) ! ed2 at midpoint
        je13d(i,j,k) = je13d(i,j,k) + fac*(ed2h-adotv1_mag(i,j,k)*be3_mag(i,j))

        fac = ped_mag(i,j,k)*adotam
        je23d(i,j,k) = fac*(ed2h-adotv1_mag(i,j,k)*be3_mag(i,j))
        
        fac = ped_mag(i,j,k)*a1dta2_mag(i,j)+hal_mag(i,j,k)
        je23d(i,j,k) = je23d(i,j,k) + fac*(ed1h+adotv2_mag(i,j,k)*be3_mag(i,j))

        if (current_pg > 0) je13d(i,j,k) = je13d(i,j,k) + &
          1.e4*je1pg_mag(i,j,k) ! convert cm -> m
      enddo ! i=mlon0,mlon1
    enddo ! j=mlat0,mlat1
  enddo ! k=mlev0,mlev1-1
!
! Periodic point for je13d is necessary for calculation of kqphi2d below.
  do k=mlev0,mlev1
    call mp_mag_periodic_f2d(je13d(:,:,k),mlon0,mlon1,mlat0,mlat1,1)
  enddo
!
! JE13D: Eastward current density (3d)
! JE23D: Downward current density (3d)
!
  do j=mlat0,mlat1
    call mkdiag_JE13D('JE13D',je13d(:,j,1:mlev1),mlon0,mlon1,1,mlev1,j)
    call mkdiag_JE23D('JE23D',je23d(:,j,1:mlev1),mlon0,mlon1,1,mlev1,j)
  enddo  
!
! Calculate K_(q,phi) (Richmond: Ionospheric
! Electrodynamics using magnetic apex coordinates pp.208 (eq 7.4))
!   K_(q,phi) = int_(h_l)^(h_u) [( [R_0/R]^0.5 *
!               * je1/D * sin(lam_q)/sin(lam_m) sin(I_m)/sin(I)*D] dh
!   with F = D*sin(lam_m)/sin(lam_q)*sin(I)/sin(I_m)*[R/R_0]^3
!
  do j=mlat0,mlat1
    sinlamq = sin(ylatm(j))                   ! sin(lam_q)
    if (j==nmlath) sinlamq = sin(ylatm(j-1))  ! ylatm is global to all tasks
    coslamq2 = 1. - sinlamq*sinlamq           ! cos^2(lam_q)

    do i=mlon0,mlon1
!
! At equator sin lam_q/sin I is set to the average otherwise quotient = 0
!   check this later 010611
!
      k = mlev0
      if (j==nmlath) then
        facsin = sinlamq/sini_mag(i,j-1) ! sin(lam_q)/sin(I)
        fac = r0/(re + max(zpotm3d(i,j-1,k),h0))
        fac = min(fac,1.)
        lamm = acos(sqrt(coslamq2*fac))  ! cos^2(lam_m) = R_0/R*cos^2(lam_q)
        lamm = sign(lamm,ylatm(j-1))     ! lam_m
      else
        facsin = sinlamq/sini_mag(i,j)         ! sin(lam_q)/sin(I)
        fac = r0/(re + max(zpotm3d(i,j,k),h0)) ! sqrt(R_0 / R)
        fac = min(fac,1.)
        lamm = acos(sqrt(coslamq2*fac))  ! cos^2(lam_m) = R_0/R*cos^2(lam_q)
        lamm = sign(lamm,ylatm(j))       ! lam_m
      endif
      sinlamm = sin(lamm)              ! sin(lam_m)
!
! sin(I_m) = sin(lam_m)/sqrt(1/4+3/4*sin^2(lam_m))
      sinim = sinlamm/sqrt(.25+.75*sinlamm**2)
      facq  = sinim/sinlamm*facsin        ! sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
      fac   = fac**0.5                    ! sqrt(R_0 / R)
!
! sqrt(R_0 / R)*sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
      act   = fac*facq

! These appear to be ok:
!     write(6,"('crdens: i=',i4,' j=',i4,' facsin,fac,lamm=',3e12.4)") &
!       i,j,facsin,fac,lamm

      kqphi2d(i,j) = 0.
      kqphi(i,j,:) = 0.
      do k = mlev0,mlev1-1
! (R_e + h)/R_0  : R_0 = R_e + h0 ; h0 = 9.e6 cm
        fac  = r0/(re + max(zpotm3d(i,j,k+1),h0)) ! R_0 / R
        lamm = acos(sqrt(coslamq2*fac))     ! cos^2(lam_m) = R_0/R*cos^2(lam_q)
        lamm = sign(lamm,ylatm(j))          ! lam_m
        sinlamm = sin(lamm)                 ! sin(lam_m)
!
! sin(I_m) = sin(lam_m)/sqrt(1/4+3/4*sin^2(lam_m))
        sinim = sinlamm/sqrt(.25+.75*sinlamm**2)
        facq  = sinim/sinlamm*facsin        ! sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
        dh    = max(zpotm3d(i,j,k+1),h0) - max(zpotm3d(i,j,k),h0)  ! dh
        dh    = dh*1.e-2                    ! convertion [cm] to [m]
        fac   = fac**0.5                    ! sqrt(R_0 / R)
!
! sqrt(R_0 / R)*sin(I_m)/sin(lam_m)*sin(lam_q)/sin(I)
        actpk = fac*facq
!
! integration value at 1/2 level
        kqphi2d(i,j) = kqphi2d(i,j) + 0.5*(actpk+act)*je13d(i,j,k)*dh
        act = actpk
        kqphi(i,j,k) = kqphi2d(i,j)  ! for output
      enddo ! k=mlev0,mlev1-1
    enddo ! i=mlon0,mlon1

!   if (current_pg > 0) call addfld('KQPHI_TOT','kqphi_int (u,E,g,p)', &
!    '[A/m]',kqphi(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)
!   if (current_pg <= 0) call addfld('KQPHI_UE','kqphi_int (u,E)', &
!    '[A/m]',kqphi(:,j,:),'mlon',mlon0,mlon1,'mlev',mlev0,mlev1,j)

  enddo ! mlat0,mlat1
!
! KQPHI: Height-integrated current density (+east)
! Note kqphi2d is dimensioned with halo points.
!
  call mkdiag_KQPHI('KQPHI',kqphi2d(mlon0:mlon1,mlat0:mlat1),mlon0,mlon1,mlat0,mlat1)
!
! Calculate K_(q,lam) (Richmond: Ionospheric
! Electrodynamics using magnetic apex coordinates pp.208 (eq 7.5))
!   K_(q,lam) = -1/cos(lam_q) int_(-pi/2)^(lam_q) [J_mr*R*cos(lam_q) +
!                d(K_(q,phi))/d(phi_q) ] d(lam_q)
!     d(K_(q,phi))/d(phi_q)
!
  call mp_allgather_2d(kqphi2d(mlon0:mlon1,mlat0:mlat1),kqphi_glb, &
    nmlonp1,nmlat,mlon0,mlon1,mlat0,mlat1,1)

  fac = 0.5/dlonm           ! 1/(2*d lon_m)
  dkqphi_glb = 0.
  do j = 1,nmlat
    do i = 2,nmlonp1-1
      dkqphi_glb(i,j) = (kqphi_glb(i+1,j)-kqphi_glb(i-1,j))*fac
    enddo
    dkqphi_glb(1,j)       = (kqphi_glb(2,j)-kqphi_glb(nmlonp1-1,j))*fac
    dkqphi_glb(nmlonp1,j) = dkqphi_glb(1,j)
  enddo ! j=mlon0,mlat1
!
! Scan longitude subdomain:
!
! Note that all tasks now have global ylatm (cons.F),
!   global dkqphi_glb(nmlonp1,nmlat), and global
!   nscrrt(nmlon0,nmlat) (from sub noso_crrt above)
!
! But noso_crrt was called by root task only, so need to broadcast 
! global nscrrt to all tasks:
!
  call mp_bcast_real(nscrrt,nmlonp1*nmlat,0)

  r0m = r0*1.e-2
  do i=mlon0,mlon1
    fsums(1) = 0.  ! south pole
    afsums(1) = 0. 
!
! [J_mr*R*cos(lam_q)+d K_(qphi)]_(j=1)
    act = nscrrt(i,1)*r0m*cos(ylatm(1))+dkqphi_glb(i,1)
!
! difflm: d | (lam_q(j)-lam_q(j-1))|
! actpk: [J_mr*R*cos(lam_q)+d K_(qphi)]_(j)
! act:  [integrand]_(j-1)+[integrand]_(j)
! act: -[integrand]_(j-1/2)] d lam_q
! fsums: int_(-pi/2)^(lam_q)[integrand] d lam_q
! afsums: int_(-pi/2)^(lam_q) | [integrand] | d lam_q
!
! Integrate from the south pole+1 to equator
!
    do j=2,nmlath
      difflm = abs(ylatm(j)-ylatm(j-1)) ! ylatm has global lats
      actpk    = nscrrt(i,j)*r0m*cos(ylatm(j))+dkqphi_glb(i,j)
      act      = act +actpk
      act      = -act/2.0*difflm
      fsums(j) = fsums(j-1)  + act
      afsums(j)= afsums(j-1) + abs(act)
      act = actpk
    enddo ! j=2,nmlath
!
! Integrate from the north pole to equator
!
    fsumn(nmlat)  = 0. ! north pole
    afsumn(nmlat) = 0. 
    act = nscrrt(i,nmlat)*r0m*cos(ylatm(nmlat))+dkqphi_glb(i,nmlat)
    do j=nmlat-1,nmlath,-1
      difflm    = abs(ylatm(j)-ylatm(j+1))
      actpk    = nscrrt(i,j)*r0m*cos(ylatm(j))+dkqphi_glb(i,j)
      act      = act +actpk
      act      = act/2.0*difflm
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
    kqlam(i,:) = 0.
    do j=mlat0,mlat1
      if (j==1) then ! south pole
        kqlam(i,j) = -dkqphi_glb(i,j) ! at south pole K_(q,lam) = -d(K_(q,phi))/d(phi)
      endif
      if (j==nmlat) then ! north pole
        kqlam(i,j) =  dkqphi_glb(i,j) ! at north pole K_(q,lam) = d(K_(q,phi)
      endif
    enddo

    epsn = 0.5*(fsums(nmlath)-fsumn(nmlath))/afsumn(nmlath)
    epss = 0.5*(fsums(nmlath)-fsumn(nmlath))/afsums(nmlath)

    do j=mlat0,mlat1
      if (j > 1 .and. j <= nmlath) then ! southern hem, excluding pole
        kqlam(i,j)= fsums(j) - epss*afsums(j)
        kqlam(i,j)= kqlam(i,j)/cos(ylatm(j))
      endif
    enddo
!
! kqlam_cor = kqlam + err/2/abs(kqlam)_north*abs(kqlam)_(lam_q)
! kqlam_cor /cos(lam_q)
!
    do j=mlat0,mlat1
      if (j < nmlat .and. j >= nmlath) then ! northern hem, excluding pole
        kqlam(i,j)= fsumn(j) + epsn*afsumn(j)
        kqlam(i,j)= kqlam(i,j)/cos(ylatm(j))
      endif
    enddo

  enddo ! i=mlon0,mlon1
!
! KQLAM: Height-integrated current density (+north)
!
  call mkdiag_KQLAM('KQLAM',kqlam,mlon0,mlon1,mlat0,mlat1)

end subroutine noso_crdens
!-----------------------------------------------------------------------
end module current

!-----------------------------------------------------------------------
subroutine nosocoef
!
! This stub is outside the current module to avoid circular dependency
! between current and pdynamo modules. This routine is called from
! sub complete_integrals (pdynamo.F).
!
  use current,only: noso_coef
  call noso_coef
end subroutine nosocoef

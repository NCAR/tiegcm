!
      module heelis_module
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Module used to calculcate the Heelis model potential in both hemispheres
! Byimf, Ctpoten and Power at a minimum using paramaters from aurora_cons
!
      use params_module,only: nmlat,nmlonp1,nmlon,nmlonp1
      implicit none
      integer,parameter :: nmlat0=(nmlat+1)/2
      real,dimension(nmlonp1,nmlat)  :: phihm  ! potential in magnetic
      contains 
!-----------------------------------------------------------------------
      subroutine heelis
!
! Heelis driver, called from sub pdynamo (pdynamo.F).
! These routines return phihm to the pdynamo.
! phihm:  Heelis potential in magnetic coordinates (single level).
!
      implicit none
!
! Calculation of pfrac for dynamo will be done by sub colath, 
! called from advance.
! pfrac:  Fractional presence of dynamo equation given critical 
!           convection colatitudes crit(2).
!
! Calculate  the heelis potential phihm in geomagnetic coordinates:
! (potm calls sub flwv32)
!
      call potm
!
      end subroutine heelis
!-----------------------------------------------------------------------
      subroutine potm
      use magfield_module,only: sunlons 
      use cons_module,only: 
     |  ylonm,ylatm, ! magnetic grid lons, lats
     |  pi           ! pi used in dynamo calculations
      implicit none
!
! Calculate heelis potential in geomagnetic coordinates.
!
! Local:
      integer :: j
      real,dimension(nmlon) :: dlon
!
      dlon(:) = ylonm(1:nmlon)-sunlons(1)
      do j=1,nmlat
!
! flwv32 returns single-level Heelis potential in geomag coords:
!
        if (abs(ylatm(j)) > pi/6.) then
          call flwv32(ylatm(j),dlon,phihm(1:nmlon,j))
        else
          phihm(1:nmlon,j) = 0.
        endif
      enddo ! j=1,nmlat
!
! Periodic points:
      do j=1,nmlat
        phihm(nmlonp1,j) = phihm(1,j)
      enddo ! j=1,nmlat
      end subroutine potm
!-----------------------------------------------------------------------
      subroutine flwv32(dlat,dlon,poten)
!
! Calculate heelis potential at current magnetic latitude mlat.
!
      use aurora_module,only: offc,dskofc,phin,phid,theta0,
     |  psie,psim,pcen,phidm0,phidp0,phinm0,phinp0,rr1
      use cons_module,only: pi
      implicit none
!
! Args:
      real,intent(in) :: dlat
      real,dimension(nmlon),intent(in) :: dlon
      real,dimension(nmlon),intent(out) :: poten
!
! Local:
      integer :: i,n,ihem
      real :: 
     |  sinthr1,psi(8),phirc,sinth0,
     |  ofdc,cosofc,sinofc,aslonc,
     |  phdpmx,phnpmx,phnmmx,phdmmx,
     |  sinlat,coslat,sincolat,
     |  width,fac1,fac2
      real,dimension(nmlon) :: sinlon,coslon,alon,colat,
     |  alonrc,alonrcrot,phifun,phifn2
      real :: phi(nmlon,8)
!
! Set ihem=1,2 for South,North hemisphere:
!
      if (dlat < 0) then
        ihem = 1
      else
        ihem = 2
      endif
!
      ofdc = sqrt(offc(ihem)**2+dskofc(ihem)**2)
      cosofc = cos(ofdc)
      sinofc = sin(ofdc)
      aslonc = asin(dskofc(ihem)/ofdc)
!
      width = phin(ihem)-phid(ihem)
      phdpmx = .5*min(pi,width)
      phnpmx = .5*min(pi,2.*pi-width)
      phnmmx = phdpmx
      phdmmx = phnpmx
!
!     write(6,"('flwv32: cosofc=',e12.4)") cosofc
!     write(6,"('flwv32: sinofc=',e12.4)") sinofc
!     write(6,"('flwv32: aslonc=',e12.4)") aslonc
!     write(6,"('flwv32: phdpmx=',e12.4,' phnpmx=',e12.4)")
!    |  phdpmx,phnpmx
!     write(6,"('flwv32: phnmmx=',e12.4,' phdmmx=',e12.4)")
!    |  phnmmx,phdmmx

      psi(1) = psie(ihem)
      psi(3) = psim(ihem)
      do n=2,4,2
        psi(n) = psi(n-1)
      enddo ! n=2,4,2
      do n=1,4
        psi(n+4) = psi(n)
      enddo ! n=1,4
!
! Transform to auroral circle coordinates:
!
      sinlat = sin(abs(dlat))
      coslat = cos(dlat)
      do i=1,nmlon
        sinlon(i) = sin(dlon(i)+aslonc)
        coslon(i) = cos(dlon(i)+aslonc)
        colat(i) = acos(cosofc*sinlat-sinofc*coslat*coslon(i))
        alon(i) = amod(atan2(sinlon(i)*coslat,
     |    sinlat*sinofc+cosofc*coslat*coslon(i))-
     |    aslonc+3.*pi,2.*pi)-pi
!
! Boundaries for longitudinal function:
!
        width = ((colat(i)-theta0(ihem))/theta0(ihem))**2
        phi(i,4)=phid(ihem)-min(phidm0(ihem)+width*
     |    (.5*pi-phidm0(ihem)),phdmmx)
        phi(i,5)=phid(ihem)+min(phidp0(ihem)+width*
     |    (.5*pi-phidp0(ihem)),phdpmx)
        phi(i,6)=phin(ihem)-min(phinm0(ihem)+width*
     |    (.5*pi-phinm0(ihem)),phnmmx)
        phi(i,7)=phin(ihem)+min(phinp0(ihem)+width*
     |    (.5*pi-phinp0(ihem)),phnpmx)
        phi(i,1)=phi(i,5)-2.*pi
        phi(i,2)=phi(i,6)-2.*pi
        phi(i,3)=phi(i,7)-2.*pi
        phi(i,8)=phi(i,4)+2.*pi
!
! Add ring current rotation to potential (phirc)
!
        phirc = 0.
        alonrc(i) = amod(alon(i)+phirc+5.*pi,2.*pi)-pi
        alonrcrot(i) = amod(alon(i)+phirc+6.*pi,2.*pi)-pi
      enddo ! i=1,nmlon
!
! Longitudinal variation:
!
      phifun(:) = 0.
      phifn2(:) = 0.
      do n=1,7
        do i=1,nmlon
          if ((alonrc(i)-phi(i,n  ))*
     |        (alonrc(i)-phi(i,n+1)) >= 0) then
            fac1 = 0.
          else
            fac1 = .5*(psi(n)+psi(n+1)+(psi(n)-psi(n+1))*
     |        cos(pi*(alonrc(i)-phi(i,n))/(phi(i,n+1)-phi(i,n))))
          endif
          phifun(i) = phifun(i)+fac1
!
          if ((alonrcrot(i)-phi(i,n  ))*
     |        (alonrcrot(i)-phi(i,n+1)) >= 0) then
            fac2 = 0.
          else
            fac2 = .5*(psi(n)+psi(n+1)+(psi(n)-psi(n+1))*
     |        cos(pi*(alonrcrot(i)-phi(i,n))/(phi(i,n+1)-phi(i,n))))
          endif
          phifn2(i) = phifn2(i)+fac2
        enddo
      enddo
!
! Evaluate total potential:
!
      do i=1,nmlon
        if (colat(i) < theta0(ihem)) then
          width = colat(i)/theta0(ihem)
          poten(i) = (2.*(pcen(ihem)-phifun(i))+
     |      (phifun(i)-phifn2(i))*0.75)*width**3 +
     |      (1.5*(phifun(i)+phifn2(i))-3.*pcen(ihem))*width**2 +
     |      0.75*(phifun(i)-phifn2(i))*width + pcen(ihem)
        else
          sincolat = sin(colat(i))
!
          sinth0 = sin(theta0(ihem))
          if (sincolat < sinth0) then
            fac1 = 1
          else
            fac1 = (sincolat/sinth0)**rr1(ihem)
          endif
!
! Average amie results show r1=-2.6 for 11.3 degrees
!   (0.1972 rad) beyond theta0.
!
          sinthr1 = sin(theta0(ihem)+0.1972)
          if (sincolat < sinthr1) then
            fac2 = 1
          else
            fac2 = exp(7.*(1.-sincolat/sinthr1))
          endif
!
          poten(i) = phifun(i)*fac1*fac2
        endif
      enddo

!     write(6,"(/'flwv32: ihem=',i2)") ihem
!     write(6,"('  theta0(ihem)=',e12.4,' pcen(ihem)=',e12.4,
!    |  ' rr1(ihem)=',e12.4)") theta0(ihem),pcen(ihem),rr1(ihem)
!     write(6,"('  sinth0=',e12.4,' sinthr1=',e12.4)") sinth0,sinthr1
!     write(6,"('  colat=',/,(6e12.4))") colat
!     write(6,"('  phifun=',/,(6e12.4))") phifun
!     write(6,"('  phifn2=',/,(6e12.4))") phifn2
!     write(6,"('  poten =',/,(6e12.4))") poten
!     write(6,"(/'flwv32: ihem=',i2,' poten=',/,(6e12.4))") ihem,poten

      end subroutine flwv32
!-----------------------------------------------------------------------
      end module heelis_module

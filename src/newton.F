!
      subroutine newton(tn,o2,o1,n2,no,barm,xnmbar,cp,cool_implicit,
     |  cool_explicit,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate implicit and explicit cooling terms.
!
      use params_module,only: dz,zibot,zpmid
      use cons_module,only: rmass_co2,p0,expz,expzmid_inv,boltz,
     |  avo,rmassinv_o2,rmassinv_o1,rmassinv_n2,rmassinv_no
      use addfld_module,only: addfld
      use diags_module,only: mkdiag_CO2COOL,mkdiag_NOCOOL
      use init_module,only: iday,iyear,secs
      use interp_module,only: interp1d
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  no,      ! nitric oxide (mmr)
     |  cp,      ! specific heat
     |  barm,    ! mean molecular weight
     |  xnmbar   ! p0*e(-z)*mbar/kT
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  cool_implicit, cool_explicit ! output
!
! Local:
      real :: co2u
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  xmco2,o2_cm3,o1_cm3,n2_cm3,no_cm3,
     |  nco2,aco2,bco2,  ! n(co2), a(co2), b(co2) (s10,s11,s12)
     |  ano,             ! a(no) (s11)
     |  co2_cool,        ! co2 contribution to cooling
     |  no_cool          ! no contribution to cooling
      integer :: k,i
      real :: uttime,yfrac
!
! Change O(3P) cooling factor from 0.5 to 1.0 of the Bates expression.
      real,parameter :: 
     |  an(3) = (/0.835E-18,0.6,0.2/),
     |  bn(3) = (/228.,228.,325./)
      real :: xfac(lev0:lev1)
      real,dimension(lev0:lev1,lon0:lon1) :: fac1,fac2,
     |  part1,part2,parta
      real,dimension(29),parameter ::
     |  zpint_5 = (/(zibot+(k-1)*0.5, k=1,29)/),
     |  xfac_5 =
     |    (/0.1000E-01, 0.1000E-01, 0.1000E-01, 0.5000E-01, 0.1000E+00,
     |      0.2000E+00, 0.4000E+00, 0.5500E+00, 0.7000E+00, 0.7500E+00,
     |      (0.8000E+00, k=1,19)/)
!
! set CO2 concentration
! calculate the fractional year value for the model time
      uttime= secs/3600.
      if (mod(iyear,4) == 0) then
         yfrac=iyear+(iday-1+(uttime/24.))/366.
      else
         yfrac=iyear+(iday-1+(uttime/24.))/365.
      endif
      if (yfrac > 1953.0) then 
         co2u=248.0+30.0*exp(0.014*(yfrac-1900.))
      else
         co2u=280.0+0.30*(yfrac-1850.0)
      endif
      if (co2u < 280.0) co2u=280.0
      co2u=co2u*1.e-6
!
      do i=lon0,lon1
        do k=lev0,lev1
          o2_cm3(k,i) = xnmbar(k,i)*o2(k,i)*rmassinv_o2
          o1_cm3(k,i) = xnmbar(k,i)*o1(k,i)*rmassinv_o1
          n2_cm3(k,i) = xnmbar(k,i)*n2(k,i)*rmassinv_n2
          no_cm3(k,i) = xnmbar(k,i)*no(k,i)*rmassinv_no
        enddo
      enddo
!
! Integral(m(co2)*dz/mbar) (s10).
! Low boundary: 
      do i=lon0,lon1
        xmco2(1,i) = .25*rmass_co2*(1./barm(1,i)+1./barm(2,i))*dz
!
! Integrate upwards:
        do k=lev0+1,lev1-1
          xmco2(k,i) = xmco2(k-1,i)+rmass_co2/barm(k,i)*dz
        enddo ! k=lev0+1,lev1-1
!
! n(co2) (s10):
        do k=lev0,lev1-1
          nco2(k,i) = co2u*p0*expz(1)*expzmid_inv/(boltz*tn(k,i))*
     |      exp(-xmco2(k,i))
!
! a(co2) (s11):
          aco2(k,i) = 2.5e-15*(1.+0.03*max(tn(k,i)-200.,0.))
!
! b(co2) (s12):
          if (tn(k,i) < 260.) then
             bco2(k,i) = 1.56e-12
          elseif (tn(k,i) >= 260..and.tn(k,i) <= 300.) then
             bco2(k,i) = (2.6-0.004*tn(k,i))*1.0e-12
          else
             bco2(k,i) = 1.4e-12
          endif
!
! Co2 cooling (s12):
          co2_cool(k,i) = 2.65e-13*nco2(k,i)*exp(-960./tn(k,i))*
     |      avo*((o2_cm3(k,i) + n2_cm3(k,i))*aco2(k,i)+
     |      o1_cm3(k,i)*bco2(k,i))/xnmbar(k,i)
!
! a(no) (s11):
          ano(k,i) = 4.2e-11*o1_cm3(k,i)+2.4e-14*o2_cm3(k,i)
!
! No cooling (s11):
          no_cool(k,i) = 4.956e-12*avo*no_cm3(k,i)/xnmbar(k,i)*
     |      ano(k,i)/(ano(k,i)+13.3)*exp(-2700./tn(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('XMCO2',' ',' ',xmco2(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('CO2_COOL',' ',' ',co2_cool(lev0:lev1-1,lon0:lon1),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('NO_COOL' ,' ',' ',no_cool(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Save as diagnostics on secondary histories (define top level as k-1):
      co2_cool(lev1,:) = co2_cool(lev1-1,:)
      no_cool (lev1,:) = no_cool (lev1-1,:)
      call mkdiag_CO2COOL('CO2_COOL',co2_cool,lev0,lev1,lon0,lon1,lat)
      call mkdiag_NOCOOL ('NO_COOL' ,no_cool ,lev0,lev1,lon0,lon1,lat)
!
! Add O(3P) cooling to implicit and explicit cooling terms.
!
      xfac = interp1d(zpmid,zpint_5,xfac_5)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          part1(k,i) = an(1)*exp(-bn(1)/tn(k,i))
          part2(k,i) = 1.+an(2)*exp(-bn(2)/tn(k,i))+
     |                    an(3)*exp(-bn(3)/tn(k,i))
          parta(k,i) = xfac(k)*avo*o1_cm3(k,i)/xnmbar(k,i)
!
! fac1 = AN(1)*XFAC*N0*PSI2/M2*EXP(-BN(1)/T) (s11)
          fac1(k,i) = part1(k,i)/part2(k,i)*parta(k,i)
!
! fac2 = DF/DT=S11/S12/T**2*(BN(1)+(BN(1)-BN(3))*AN(3)*EXP(-BN(3)/T)) (s12)
          fac2(k,i) = fac1(k,i)/part2(k,i)*
     |      (bn(1)+(bn(1)-bn(3))*an(3)*exp(-bn(3)/tn(k,i)))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! Total cooling terms.
          cool_implicit(k,i) =
     |      (2700.*no_cool(k,i)+960.*co2_cool(k,i)+fac2(k,i))/
     |      (.5*(cp(k,i)+cp(k+1,i))*tn(k,i)**2)
          cool_explicit(k,i) = expz(k)*((1.-2700./tn(k,i))*no_cool(k,i)+
     |      (1.-960./tn(k,i))*co2_cool(k,i)+fac1(k,i)-fac2(k,i)/tn(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('COOL_IMP' ,' ',' ',cool_implicit(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('COOL_EXP' ,' ',' ',cool_explicit(:,lon0:lon1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      end subroutine newton

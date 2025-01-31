!
      subroutine ionvel(z,tn,un,vn,o2,o1,he,n2,w,ti,xnmbar,scht,
     |  ui,vi,wi,Etot,lev0,lev1,lon0,lon1,lat0,lat1)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! tiegcmlicense.txt.
!
! Calculate ion velocities from electric field.
! (electric field ex,ey,ez was calculated in pdynamo module)
!
! vi = (nu*a + omega*axb) / (nu^2 + omega^2)
!   where a = g + bz*g*b + omega*E/B + nu*(u - (b.u)*b)
!
      use params_module,only: nlonp4,nlonp2
      use cons_module,only: re,grav,
     |  rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2
      use pdynamo_module,only: ! (nlevp1,lon0-2:lon1+2,lat0:lat1)
     |  ex_dyn=>ex,ey_dyn=>ey,ez_dyn=>ez
      use input_module,only: colfac
      use magfield_module,only: rjac,bx,by,bz,bmod2
      use addfld_module,only: addfld
      use mpi_module,only: mp_periodic_f3d
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(in) :: z,tn,un,vn,o2,o1,he,n2,w,ti,xnmbar,scht
      real,dimension(lev0:lev1,lon0-2:lon1+2,lat0-2:lat1+2),
     |  intent(out) :: ui,vi,wi,Etot
!
! Local:
      real,parameter :: qe = 1.602e-19, mp = 1.6726e-24, maxvel = 2e5
      integer :: k,i,lat,lonbeg,lonend,ier
      real,dimension(lon0:lon1,lat0:lat1) :: omega
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1) ::
     |  ex,ey,ez,exi,eyi,ezi,ax,ay,az,uii,vii,wii,
     |  wn,tr,o2_cm3,o1_cm3,he_cm3,n2_cm3,nu,bu,denom
      real,dimension(lev0:lev1,lon0:lon1,lat0:lat1,3) :: ftmp
!
! VT vampir tracing:
!
#ifdef VT
#include <VT.inc>
#endif
!
#ifdef VT
!     code = 121 ; state = 'ionvel' ; activity='ModelCode'
      call vtbegin(121,ier)
#endif
!
      lonbeg = lon0
      if (lon0==1) lonbeg = 3
      lonend = lon1
      if (lon1==nlonp4) lonend = nlonp2
!
      do lat=lat0,lat1
        do i=lonbeg,lonend
!
! Rotate ex and ey to geographic orientation using Jacobian
! Divide by distance from center of earth, 1e2 cm -> 1 m
          do k=lev0,lev1
            ex(k,i,lat) = 1e2/(re+z(k,i,lat))*
     |        (rjac(i-2,lat,1,1)*ex_dyn(k,i,lat)+
     |         rjac(i-2,lat,2,1)*ey_dyn(k,i,lat))
            ey(k,i,lat) = 1e2/(re+z(k,i,lat))*
     |        (rjac(i-2,lat,1,2)*ex_dyn(k,i,lat)+
     |         rjac(i-2,lat,2,2)*ey_dyn(k,i,lat))
          enddo ! k=lev0,lev1
!
! Centered difference
          do k=lev0+1,lev1-1
            ez(k,i,lat) = ez_dyn(k,i,lat)*1e2/
     |        (z(k+1,i,lat)-z(k-1,i,lat))
          enddo ! k=lev0+1,lev1-1
!
! One-sided difference for lower and upper boundaries
          ez(lev0,i,lat) = ez_dyn(lev0,i,lat)*1e2/
     |      ((z(lev0+1,i,lat)-z(lev0,i,lat))*2)
          ez(lev1,i,lat) = ez_dyn(lev1,i,lat)*1e2/
     |      ((z(lev1,i,lat)-z(lev1-1,i,lat))*2)
        enddo ! i=lonbeg,lonend
      enddo ! lat=lat0,lat1
!
! Periodic points for electric fields
      ftmp(:,:,:,1) = ex
      ftmp(:,:,:,2) = ey
      ftmp(:,:,:,3) = ez
!
      call mp_periodic_f3d(ftmp,lev0,lev1,lon0,lon1,lat0,lat1,3)
!
      ex = ftmp(:,:,:,1)
      ey = ftmp(:,:,:,2)
      ez = ftmp(:,:,:,3)
!
      do lat=lat0,lat1
        do i=lon0,lon1
          omega(i,lat) = qe*bmod2(i,lat)*rmassinv_o1/(mp*10)
!
          do k=lev0,lev1-1
!
! Interfaces to midpoints
            wn(k,i,lat) = scht(k,i,lat)*
     |        (w(k,i,lat)+w(k+1,i,lat))/2
            exi(k,i,lat) = (ex(k,i,lat)+ex(k+1,i,lat))/2
            eyi(k,i,lat) = (ey(k,i,lat)+ey(k+1,i,lat))/2
            ezi(k,i,lat) = (ez(k,i,lat)+ez(k+1,i,lat))/2
!
            tr(k,i,lat) = (tn(k,i,lat)+ti(k,i,lat))/2
            o2_cm3(k,i,lat) = xnmbar(k,i,lat)*o2(k,i,lat)*rmassinv_o2
            o1_cm3(k,i,lat) = xnmbar(k,i,lat)*o1(k,i,lat)*rmassinv_o1
            he_cm3(k,i,lat) = xnmbar(k,i,lat)*he(k,i,lat)*rmassinv_he
            n2_cm3(k,i,lat) = xnmbar(k,i,lat)*n2(k,i,lat)*rmassinv_n2
            nu(k,i,lat) = 1e-10*
     |        (6.64 *o2_cm3(k,i,lat)+
     |         1.32 *he_cm3(k,i,lat)+
     |         6.82 *n2_cm3(k,i,lat)+
     |         0.367*o1_cm3(k,i,lat)*sqrt(tr(k,i,lat))*
     |         (1-0.064*log10(tr(k,i,lat)))**2*colfac)
            bu(k,i,lat) = bx(i,lat)*un(k,i,lat)+
     |        by(i,lat)*vn(k,i,lat)+bz(i,lat)*wn(k,i,lat)
!
! Auxiliary acceleration
! 1e4 G -> 1 T, 1 m/s2 -> 1e2 cm/s2
            ax(k,i,lat) = bz(i,lat)*grav*bx(i,lat)+
     |        omega(i,lat)*exi(k,i,lat)*1e6/bmod2(i,lat)+
     |        nu(k,i,lat)*(un(k,i,lat)-bu(k,i,lat)*bx(i,lat))
            ay(k,i,lat) = bz(i,lat)*grav*by(i,lat)+
     |        omega(i,lat)*eyi(k,i,lat)*1e6/bmod2(i,lat)+
     |        nu(k,i,lat)*(vn(k,i,lat)-bu(k,i,lat)*by(i,lat))
            az(k,i,lat) = -grav*(1-bz(i,lat)**2)+
     |        omega(i,lat)*ezi(k,i,lat)*1e6/bmod2(i,lat)+
     |        nu(k,i,lat)*(wn(k,i,lat)-bu(k,i,lat)*bz(i,lat))
!
! Ion velocity
            denom(k,i,lat) = nu(k,i,lat)**2+omega(i,lat)**2
            uii(k,i,lat) = 1/denom(k,i,lat)*
     |        (nu(k,i,lat)*ax(k,i,lat)+omega(i,lat)*
     |        (ay(k,i,lat)*bz(i,lat)-az(k,i,lat)*by(i,lat)))
            vii(k,i,lat) = 1/denom(k,i,lat)*
     |        (nu(k,i,lat)*ay(k,i,lat)+omega(i,lat)*
     |        (az(k,i,lat)*bx(i,lat)-ax(k,i,lat)*bz(i,lat)))
            wii(k,i,lat) = 1/denom(k,i,lat)*
     |        (nu(k,i,lat)*az(k,i,lat)+omega(i,lat)*
     |        (ax(k,i,lat)*by(i,lat)-ay(k,i,lat)*bx(i,lat)))
          enddo ! k=lev0,lev1-1
!
! Midpoints to interfaces
          do k=lev0+1,lev1-1
            ui(k,i,lat) = (uii(k-1,i,lat)+uii(k,i,lat))/2
            vi(k,i,lat) = (vii(k-1,i,lat)+vii(k,i,lat))/2
            wi(k,i,lat) = (wii(k-1,i,lat)+wii(k,i,lat))/2
          enddo ! k=lev0+1,lev1-1
!
! Extrapolate lower and upper boundaries
          ui(lev0,i,lat) = 1.5*uii(lev0,i,lat)-0.5*uii(lev0+1,i,lat)
          vi(lev0,i,lat) = 1.5*vii(lev0,i,lat)-0.5*vii(lev0+1,i,lat)
          wi(lev0,i,lat) = 1.5*wii(lev0,i,lat)-0.5*wii(lev0+1,i,lat)
          ui(lev1,i,lat) = 1.5*uii(lev1-1,i,lat)-0.5*uii(lev1-2,i,lat)
          vi(lev1,i,lat) = 1.5*vii(lev1-1,i,lat)-0.5*vii(lev1-2,i,lat)
          wi(lev1,i,lat) = 1.5*wii(lev1-1,i,lat)-0.5*wii(lev1-2,i,lat)
!
          do k=lev0,lev1
!
! MAGE modification
#ifdef GAMERA
            if (ui(k,i,lat) >  maxvel) ui(k,i,lat) =  maxvel
            if (ui(k,i,lat) < -maxvel) ui(k,i,lat) = -maxvel
            if (vi(k,i,lat) >  maxvel) vi(k,i,lat) =  maxvel
            if (vi(k,i,lat) < -maxvel) vi(k,i,lat) = -maxvel
#endif
            Etot(k,i,lat) = sqrt(ex(k,i,lat)**2+
     |        ey(k,i,lat)**2+ez(k,i,lat)**2)
          enddo ! k=lev0,lev1
        enddo ! i=lon0,lon1
      enddo ! lat=lat0,lat1
!
      do lat=lat0,lat1
        call addfld('EX_GEO','Eastward Electric Field','V/m',
     |    ex(:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('EY_GEO','Northward Electric Field','V/m',
     |    ey(:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
        call addfld('EZ_GEO','Upward Electric Field','V/m',
     |    ez(:,:,lat),'lev',lev0,lev1,'lon',lon0,lon1,lat)
      enddo ! lat=lat0,lat1
!
#ifdef VT
!     code = 121 ; state = 'ionvel' ; activity='ModelCode'
      call vtend(121,ier)
#endif
      end subroutine ionvel

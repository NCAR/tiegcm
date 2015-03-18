!
      subroutine lsqdsq(fjp,fj,fjm,fout,lonbeg,lonend,lev0,lev1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Subroutine to calculate L**2DEL**2(PSI)
! On input:
!   lat is current latitude index
!   fjp,fj,fjm are input arrays at lat+1, lat, and lat-1 respectively
! On output:
!   fout is defined with result
!
      use cons_module,only: cs,dphi,dlamda,re
      implicit none
!
! Args:
      integer,intent(in) :: lat,lonbeg,lonend,lev0,lev1
      real,intent(in) :: fjp(lev0:lev1,lonbeg-1:lonend+1),
     |                    fj(lev0:lev1,lonbeg-1:lonend+1),
     |                   fjm(lev0:lev1,lonbeg-1:lonend+1)
      real,intent(out) :: fout(lev0:lev1,lonbeg:lonend)
!
! Local:
      real :: flenl,flenp,con1,con2,dpsi,dpsip,dpsim
      integer :: k,i
!
! Longitude and latitude length functions
! flenl = a*abs(cos(lat))*pi/kut(lat)
! flenp = a*dphi
!
      flenl=re*abs(cs(lat))*dlamda
      flenp=re*dphi
      con1=flenl*flenl/(re*cs(lat)*dlamda)**2
      con2=flenp*flenp/((re*dphi)**2*cs(lat))
!
! fout = L*L/(A*CS(J)*DLAMBDA)**2*(PSI(I+1) - 2*PSI(I) + PSI(I-1))
!
      do i=lonbeg,lonend
        do k=lev0,lev1
          fout(k,i) = ((fj(k,i)*(-2.)+fj(k,i+1))+fj(k,i-1))*con1
!
! DPSIP=(PSI(J+1)-PSI(J))*.5*(CS(J+1)+CS(J))
!
          dpsip = (fjp(k,i)-fj(k,i))*(0.5*(cs(lat+1)+cs(lat)))
!
! DPSIM=(PSI(J)-PSI(J-1))*.5*(CS(J)+CS(J-1))
!
          dpsim = (fj(k,i)-fjm(k,i))*(0.5*(cs(lat)+cs(lat-1)))
!
! DPSI = (DPSIP-DPSIM)*(LP**2/(A*A*CS(J)*DPHI*2)
!
          dpsi = (dpsip-dpsim)*con2
!
! S = S + DPSI = (L**2)(D**2)PSI
!
          fout(k,i) = fout(k,i) + dpsi
        enddo
      enddo
      end subroutine lsqdsq

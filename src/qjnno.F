!
      subroutine qjnno(o2,o1,ne,no,n4s,n2d,xnmbar,
     |  lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculates additions to neutral gas heating and O2 dissociation
! due to N, NO chemistry.
!
      use cons_module,only: evergs,avo,rmassinv_n4s,rmassinv_no,
     |  rmassinv_o2,rmassinv_o1,rmassinv_n2d
      use chemrates_module,only: beta1,beta2,beta3,beta4,beta5,beta6
      use qrj_module,only: qtotal ! qtotal(nlevp1,lon0:lon1,lat0:lat1)
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  o2,    ! molecular oxygen (mmr)
     |  o1,    ! atomic oxygen (mmr)
     |  ne,    ! electron density (cm3)
     |  no,    ! nitric oxide
     |  n4s,   ! N(4S)
     |  n2d,   ! N(2D) (updated from comp_n2d)
     |  xnmbar ! p0*e(-z)*mbar/kT
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) :: deltaq,dq,
     |  o2_cm3,o1_cm3,no_cm3,n4s_cm3,n2d_cm3
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          o2_cm3(k,i) = xnmbar(k,i)*o2(k,i)*rmassinv_o2
          o1_cm3(k,i) = xnmbar(k,i)*o1(k,i)*rmassinv_o1
          no_cm3(k,i) = xnmbar(k,i)*no(k,i)*rmassinv_no
          n4s_cm3(k,i) = xnmbar(k,i)*n4s(k,i)*rmassinv_n4s
          n2d_cm3(k,i) = xnmbar(k,i)*n2d(k,i)*rmassinv_n2d
!
          deltaq(k,i) = 
     |      evergs*avo/xnmbar(k,i)*
     |      (n4s_cm3(k,i)*
     |       (beta1(k,i,lat)*o2_cm3(k,i)*1.4+
     |        beta3(k,i,lat)*no_cm3(k,i)*2.68)+
     |       n2d_cm3(k,i)*
     |       (beta2(k,i,lat)*o2_cm3(k,i)*1.84+
     |        beta4*o1_cm3(k,i)*2.38+
     |        beta5(k,i,lat)*.5*(ne(k,i)+ne(k+1,i))*2.38+
     |        beta6*no_cm3(k,i)*5.63))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('DELTAQ' ,' ',' ',deltaq(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',lon0,lon1,lat)
!     call addfld('QRJ_Q'  ,' ',' ',qtotal(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Add delta(q) to q:
      do i=lon0,lon1
        do k=lev0+1,lev1-1
          dq(k,i) = deltaq(k-1,i)*deltaq(k,i)
          if (dq(k,i) < 1.e-20) dq(k,i) = 1.e-20
          qtotal(k,i,lat) = qtotal(k,i,lat)+sqrt(dq(k,i)) 
        enddo ! k=lev0+1,lev1-1
!
! Lower and upper boundaries:
        qtotal(lev0,i,lat) = qtotal(lev0,i,lat)+
     |    1.5*deltaq(lev0,i)-0.5*deltaq(lev0+1,i)
        qtotal(lev1,i,lat) = qtotal(lev1,i,lat)+
     |    1.5*deltaq(lev1-1,i)-0.5*deltaq(lev1-2,i)
      enddo ! i=lon0,lon1

!     call addfld('QJNNO_Q' ,' ',' ',qtotal(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine qjnno

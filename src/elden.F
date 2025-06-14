!
      subroutine elden(xnmbar,op,op_upd,o2,o1,n2,n2d,no,n4s,xiop2p,
     |  xiop2d,nplus,n2p,nop,o2p,electrons,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Solve for electron density.
! Also calculate ions N+ (nplus), N2+ (n2p), NO+ (nop), and O2+ (o2p).
!
      use cons_module,only: rmassinv_o2,rmassinv_n2d,
     |  rmassinv_o1,rmassinv_n2,rmassinv_no,rmassinv_n4s
      use chemrates_module,only: rk1,rk2,rk3,rk4,rk5,rk6,rk7,rk8,
     |  rk9,rk10,rk16,rk23,rk26,beta9,ra1,ra2,ra3
      use qrj_module,only: qnp,qnop,qo2p,qn2p
      use addfld_module,only: addfld
      use params_module,only: rp
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Input args: press vs longitude input fields (2d (k,i)):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  xnmbar,! for conversion from mmr to cm3
     |  op,    ! O+ ion (current time-step)
     |  op_upd,! O+ ion (updated, from sub oplus)
     |  o2,    ! molecular oxygen (mmr)
     |  o1,    ! atomic oxygen (mmr)
     |  n2,    ! molecular nitrogen (mmr)
     |  n2d,   ! n2d
     |  no,    ! nitric oxide
     |  n4s,   ! n4s
     |  xiop2p,! from oplus
     |  xiop2d ! from oplus
!
! Output args (particles/cm3):
      real(rp),dimension(lev0:lev1,lon0-2:lon1+2),intent(out) ::
     |  nplus,     ! N+  output
     |  n2p,       ! N2+ output
     |  nop,       ! NO+ output
     |  o2p,       ! O2+ output
     |  electrons  ! electron density (output to f4d(ne))
!
! Local:
      integer :: k,i,i0,i1
      real(rp),dimension(lev0:lev1,lon0:lon1) ::
     |  a0,a1,a2,a3,a4, ! coefficients for quartic solver
     |  a,b,c,d,e,fg,h, ! terms for quartic coefficients
     |  root,           ! output from quartic solver
     |  o2_cm3,o1_cm3

      real(rp) :: nelimit = 100._rp

!     write(6,"('enter elden: lat=',i2)") lat
      i0 = lon0 ; i1 = lon1
!
!     call addfld('OP_ELD' ,' ',' ',op    (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('OP1_ELD',' ',' ',op_upd(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O2_ELD' ,' ',' ',o2    (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('O1_ELD' ,' ',' ',o1    (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N2D_ELD',' ',' ',n2d   (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('NO_ELD' ,' ',' ',no    (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('N4S_ELD',' ',' ',n4s   (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('XIOP2P' ,' ',' ',xiop2p(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('XIOP2D' ,' ',' ',xiop2d(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! N+ ion (output): (Haonan Wu updated 2024/01)
! With upward extension, around Z=10 near the winter pole,
! the divisor becomes too small, leading to N+ explosion.
! This is either due to O+ or N2D too high at these altitudes.
! Cap minimum O1 number density here to prevent N+ explosion.
          o2_cm3(k,i) = o2(k,i)*rmassinv_o2*xnmbar(k,i)
          o1_cm3(k,i) = o1(k,i)*rmassinv_o1*xnmbar(k,i)
          if (o1_cm3(k,i) < 1e6_rp) o1_cm3(k,i) = 1e6_rp
          nplus(k,i) = (0.5_rp*(qnp(k,i,lat)+qnp(k+1,i,lat))+
     |      rk10*op(k,i)*n2d(k,i)*xnmbar(k,i)*rmassinv_n2d) /
     |      ((rk6+rk7)*o2_cm3(k,i)+rk8*o1_cm3(k,i))
!
! Set up terms for quartic coefficients:
!
! A = QI(NO+)+K2*N(O+)*N(N2)+K7*N(N+)*N(02)+B9*N(NO) (s10)
!
          a(k,i) = .5_rp*(qnop(k,i,lat)+qnop(k+1,i,lat))+xnmbar(k,i)*
     |      (rk2(k,i,lat)*op_upd(k,i)*n2(k,i)*rmassinv_n2+
     |      rk7*nplus(k,i)*o2(k,i)*rmassinv_o2+
     |      .5_rp*(beta9(k,i,lat)+beta9(k+1,i,lat))*no(k,i)*rmassinv_no)
!
! B = QI(O2+)+K1*N(O+)*N(O2)+K6*N(N+)*N(02) (s9)
! (very small "diamond diffs" with tgcm15 due to op_upd)
!
          b(k,i) = .5_rp*(qo2p(k,i,lat)+qo2p(k+1,i,lat))+xnmbar(k,i)*
     |      (rk1(k,i,lat)*op_upd(k,i)+rk6*nplus(k,i))*o2(k,i)*
     |       rmassinv_o2+rk26*xiop2d(k,i)*o2(k,i)*rmassinv_o2
!
! C = K4*N(N4S)+K5*N(NO) (s8)
!
          c(k,i) = xnmbar(k,i)*(rk4*n4s(k,i)*rmassinv_n4s+
     |       rk5*no(k,i)*rmassinv_no)
!
! D = QI(N2+) (s7)
!
          d(k,i) = .5_rp*(qn2p(k,i,lat)+qn2p(k+1,i,lat))+
     |      (rk16*xiop2p(k,i)+
     |      rk23*xiop2d(k,i))*n2(k,i)*rmassinv_n2
!
! E =K3*N(O)+K9*N(O2) (s6)
!
          e(k,i) = xnmbar(k,i)*(rk3(k,i,lat)*o1(k,i)*rmassinv_o1+
     |       rk9*o2(k,i)*rmassinv_o2)
!
! F+G = N(O+)+N(N+) (s5)
! (very small "diamond diffs" with tgcm15 due to op_upd)
!
          fg(k,i) = op_upd(k,i)+nplus(k,i)
!
! H = K9*N(02)
!
          h(k,i) = xnmbar(k,i)*rk9*o2(k,i)*rmassinv_o2
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Coefficients for quartic solver: a0,a1,a2,a3,a4
!
! A0 = -e*c*(a+b+d) (s15)
!
          a0(k,i) = -e(k,i) * c(k,i) * (a(k,i) + b(k,i) + d(k,i))
!
! A1 = -(ra1*(e*(c*fg+b)+d*(c+h))+ra2*(e*(a+d)-h*d)+ra3*c*(a+b))/4. (s14)
!
          a1(k,i) = -(ra1(k,i,lat)*(e(k,i)*(c(k,i)*fg(k,i)+b(k,i))+
     |      d(k,i)*(c(k,i)+h(k,i)))+ra2(k,i,lat)*(e(k,i)*(a(k,i)+
     |      d(k,i))-h(k,i)*d(k,i))+ra3(k,i,lat)*c(k,i)*(a(k,i)+b(k,i)))/
     |      4._rp
!
! A2 = (ra1*(e*c-(ra2*e+ra3*c)*fg-ra2*d-ra3*b)-ra2*ra3*a)/6. (s13)
!
          a2(k,i) = (ra1(k,i,lat)*(e(k,i)*c(k,i)-(ra2(k,i,lat)*e(k,i)+
     |      ra3(k,i,lat)*c(k,i))*fg(k,i)-ra2(k,i,lat)*d(k,i)-
     |      ra3(k,i,lat)*b(k,i))-ra2(k,i,lat)*ra3(k,i,lat)*a(k,i))/6.
!
! A3 = (ra1*(ra2*e+ra3*c-ra2*ra3*fg))/4. (s12)
!
          a3(k,i) = (ra1(k,i,lat)*(ra2(k,i,lat)*e(k,i)+ra3(k,i,lat)*
     |      c(k,i)-ra2(k,i,lat)*ra3(k,i,lat)*fg(k,i)))/4.
!
! A4 = ra1*ra2*ra3 (s11)
!
          a4(k,i) = ra1(k,i,lat)*ra2(k,i,lat)*ra3(k,i,lat)
!
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1
!
!     call addfld('XNMBAR' ,' ',' ',xnmbar(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('NPLUS'   ,' ',' ',nplus(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('A_COEF'  ,' ',' ',a,'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('B_COEF'  ,' ',' ',b,'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('C_COEF'  ,' ',' ',c,'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('D_COEF'  ,' ',' ',d,'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('E_COEF'  ,' ',' ',e,'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('FG_COEF' ,' ',' ',fg,'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('H_COEF'  ,' ',' ',h,'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('A0'      ,' ',' ',a0(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('A1'      ,' ',' ',a1(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('A2'      ,' ',' ',a2(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('A3'      ,' ',' ',a3(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!     call addfld('A4'      ,' ',' ',a4(lev0:lev1-1,:),
!    |  'lev',lev0,lev1-1,'lon',i0,i1,lat)
!
! Solve quartic. Vquart returns electron density Ne in root:
!
      call vquart(a0,a1,a2,a3,a4,root,lev0,lev1,lon0,lon1)

!     call addfld('ROOT',' ',' ',root,'lev',lev0,lev1,'lon',i0,i1,lat)
!
! 1/24/08 btf, maute: Minimum Ne is replaced by new values for flux 
!                     parameter al in qinite.F
!
! Calculate N2+, O2+, NO+ (cm3)
!
      do i=lon0,lon1      
        do k=lev0,lev1-1
          if (root(k,i) < 1._rp) root(k,i) = 1.0_rp  ! insure positive Ne from solver
	                                       ! in case there is a problem
           n2p(k,i) = d(k,i)/(e(k,i)+ra3(k,i,lat)*root(k,i))    
          o2p(k,i) = (b(k,i)+h(k,i)*d(k,i)/(e(k,i)+ra3(k,i,lat)*
     |      root(k,i)))/(c(k,i)+ra2(k,i,lat)*root(k,i))
!
! nop = (a+d*(e-h)/(e+ra3*root)+c*(b+h*d/(e+ra3*root))/(c+ra2*root))/(ra1*root)
!
          nop(k,i)=(a(k,i)+d(k,i)*(e(k,i)-h(k,i))/(e(k,i)+ra3(k,i,lat)*
     |     root(k,i))+c(k,i)*(b(k,i)+h(k,i)*d(k,i)/(e(k,i)+ra3(k,i,lat)*
     |     root(k,i)))/(c(k,i)+ra2(k,i,lat)*root(k,i)))/(ra1(k,i,lat)*
     |     root(k,i))
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1      
!
!     call addfld('NPLUSb' ,' ',' ',nplus(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('N2P_ELD',' ',' ',n2p  (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('O2P_ELD',' ',' ',o2p  (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',i0,i1,lat)
!     call addfld('NOP_ELD',' ',' ',nop  (:,i0:i1),
!    |  'lev',lev0,lev1,'lon',i0,i1,lat)
!
! Transfer root to electrons output array:
      do i=lon0,lon1
        do k=lev0,lev1-2
          electrons(k+1,i) = sqrt(root(k,i)*root(k+1,i))
          if (electrons(k+1,i) .lt. nelimit)
     |       electrons(k+1,i) = nelimit
        enddo ! k=lev0,lev1-2
!
! Lower and upper boundaries:
        electrons(lev0,i) = sqrt(root(lev0  ,i)**3/root(lev0+1,i))
        electrons(lev1,i) = sqrt(root(lev1-1,i)**3/root(lev1-2,i))
        if (electrons(lev0,i) .lt. nelimit)
     |      electrons(lev0,i) = nelimit
        if (electrons(lev1,i) .lt. nelimit)
     |      electrons(lev1,i) = nelimit
      enddo ! i=lon0,lon1
!     call addfld('NE_ELDEN',' ',' ',electrons(:,i0:i1),
!    |  'lev',lev0,lev1,'lon',i0,i1,lat)

      end subroutine elden
!-----------------------------------------------------------------------
      subroutine vquart(a0,a1,a2,a3,a4,root,lev0,lev1,lon0,lon1)
      use addfld_module,only: addfld
      use params_module,only: rp
      implicit none
!
! Determines five roots of the equation:
!   a4*x**4 + 4.*a3*x**3 + 6.*a2*x**2 + 4.*a1*x + a0 = 0.
!
! Procedure is specificlly designed for real quartics with real roots
! only one of which is positive.
!
! This is called by elden for electron density.
! 
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1
      real,dimension(lev0:lev1,lon0:lon1),intent(in)  :: a0,a1,a2,a3,a4
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: root
!
! Local:
      integer :: k,i,nlevs,i0,i1
      real :: delta
      real,dimension(lev0:lev1,lon0:lon1) :: w1,w2,w3 ! work arrays
!
      nlevs = lev1-lev0+1
      i0 = lon0 ; i1 = lon1
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! w1 = ch
          w1(k,i) = -(a4(k,i)*a0(k,i)-4._rp*a3(k,i)*a1(k,i)+
     |                3._rp*a2(k,i)**2)/12._rp
!
! w2 = cg
          w2(k,i) = (a4(k,i)*(a2(k,i)*a0(k,i)-a1(k,i)**2)-
     |      a3(k,i)*(a3(k,i)*a0(k,i)-a1(k,i)*a2(k,i))+
     |      a2(k,i)*(a3(k,i)*a1(k,i)-a2(k,i)**2))/4._rp
!
! root=rlam=-2.*real((.5*(cg+sqrt(cmplx(cg**2+4.*ch**3))))**(1./3.))
! 
          root(k,i) = -2._rp*real((.5*(w2(k,i)+
     |      sqrt(cmplx(w2(k,i)**2+4._rp*w1(k,i)**3))))**(1./3._rp))
!
! W1=P=SQRT(A(5)*RLAM+A(4)**2-A(5)*A(3))
!
          delta = a4(k,i)*root(k,i)+a3(k,i)**2-a4(k,i)*a2(k,i)
          if (delta <= 0._rp) then
            w1(k,i) = 0._rp
          else
            w1(k,i) = sqrt(delta)
          endif
!
! W2=Q=SQRT((2.*RLAM+A(3))**2-A(5)*A(1))
!
          delta = (2._rp*root(k,i)+a2(k,i))**2-a4(k,i)*a0(k,i)
          if (delta <= 0._rp) then
            w2(k,i) = 0._rp
          else
            w2(k,i) = sqrt(delta)
          endif
!
! W3=PQ=2.*A(4)*RLAM+A(4)*A(3)-A(5)*A(2)
!
          w3(k,i) = 2._rp*a3(k,i)*root(k,i)+a3(k,i)*a2(k,i)-
     |              a4(k,i)*a1(k,i)
!
!  W1=P=SIGN(P,Q*PQ)
!
          w1(k,i) = sign(w1(k,i),w2(k,i)*w3(k,i))
!
! W3=P-A4
!
          w3(k,i) = w1(k,i)-a3(k,i)
!
! Final evaluation of root:
!
          delta = w3(k,i)**2-a4(k,i)*(a2(k,i)+2._rp*root(k,i)-w2(k,i))
          if (delta <= 0._rp) then
            root(k,i) = w3(k,i)/a4(k,i)
          else
            root(k,i) = (w3(k,i)+sqrt(delta))/a4(k,i)
          endif
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

      end subroutine vquart

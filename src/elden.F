!
      subroutine elden(tn,barm,op,op_upd,o2,o1,n2d,no,n4s,xiop2p,xiop2d,
     |  nplus,n2p,nop,o2p,electrons,lev0,lev1,lon0,lon1,lat)
!
! Solve for electron density.
! Also calculate ions N+ (nplus), N2+ (n2p), NO+ (nop), and O2+ (o2p).
!
      use cons_module,only: p0,expz,boltz,rmassinv_o2,rmassinv_n2d,
     |  rmassinv_o1,rmassinv_n2,rmassinv_no,rmassinv_n4s
      use chemrates_module,only: rk1,rk2,rk3,rk4,rk5,rk6,rk7,rk8,
     |  rk9,rk10,rk16,rk23,rk26,beta9,ra1,ra2,ra3
      use qrj_module,only: qnp,qnop,qo2p,qn2p
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! Input args: press vs longitude input fields (2d (k,i)):
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  tn,    ! neutral temperature (deg K)
     |  barm,  ! mean molecular weight
     |  op,    ! O+ ion (current time-step)
     |  op_upd,! O+ ion (updated, from sub oplus)
     |  o2,    ! molecular oxygen (mmr)
     |  o1,    ! atomic oxygen (mmr)
     |  n2d,   ! n2d
     |  no,    ! nitric oxide
     |  n4s,   ! n4s
     |  xiop2p,! from oplus
     |  xiop2d ! from oplus
!
! Output args (particles/cm3):
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(out) :: 
     |  nplus,     ! N+  output
     |  n2p,       ! N2+ output
     |  nop,       ! NO+ output
     |  o2p,       ! O2+ output
     |  electrons  ! electron density (output to f4d(ne))
!
! Local:
      integer :: k,i
      integer ::i0,i1,nk,nkm1 ! for addfsech
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  xnmbarm,        ! for conversion from mmr to cm3
     |  a0,a1,a2,a3,a4, ! coefficients for quartic solver
     |  a,b,c,d,e,fg,h, ! terms for quartic coefficients
     |  xn2,            ! n2 (mmr) (1-o2-o)
     |  root            ! output from quartic solver

!     write(6,"('enter elden: lat=',i2)") lat
!
! Convenience ints for addfsech calls:
      i0 = lon0
      i1 = lon1
      nk = lev1-lev0+1
      nkm1 = nk-1

!
! Bad OP1_ELD
!     call addfsech('TN_ELD' ,' ',' ',tn    (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('BAR_ELD',' ',' ',barm  (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('OP_ELD' ,' ',' ',op    (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('OP1_ELD',' ',' ',op_upd(:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('O2_ELD' ,' ',' ',o2    (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('O1_ELD' ,' ',' ',o1    (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('N2D_ELD',' ',' ',n2d   (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('NO_ELD' ,' ',' ',no    (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('N4S_ELD',' ',' ',n4s   (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('XIOP2P' ,' ',' ',xiop2p(:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('XIOP2D' ,' ',' ',xiop2d(:,i0:i1),i0,i1,nk,nkm1,lat)
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xnmbarm(k,i) = p0*expz(k)*.5*(barm(k,i)+barm(k+1,i))/
     |      (boltz*tn(k,i))
          xn2(k,i) = (1.-o2(k,i)-o1(k,i)) ! n2 (mmr)
!
! N+ ion (output):
          nplus(k,i) = (0.5*(qnp(k,i,lat)+qnp(k+1,i,lat))+
     |      rk10*op(k,i)*n2d(k,i)*xnmbarm(k,i)*rmassinv_n2d) /
     |      (xnmbarm(k,i)*((rk6+rk7)*o2(k,i)*rmassinv_o2+
     |      rk8*o1(k,i)*rmassinv_o1))
!
! Set up terms for quartic coefficients:
!
! A = QI(NO+)+K2*N(O+)*N(N2)+K7*N(N+)*N(02)+B9*N(NO) (s10)
!
          a(k,i) = .5*(qnop(k,i,lat)+qnop(k+1,i,lat))+xnmbarm(k,i)*
     |      (rk2(k,i,lat)*op_upd(k,i)*xn2(k,i)*rmassinv_n2+
     |       rk7*nplus(k,i)*o2(k,i)*rmassinv_o2+
     |       .5*(beta9(k,i,lat)+beta9(k+1,i,lat))*no(k,i)*rmassinv_no)
!
! B = QI(O2+)+K1*N(O+)*N(O2)+K6*N(N+)*N(02) (s9)
! (very small "diamond diffs" with tgcm15 due to op_upd)
!
          b(k,i) = .5*(qo2p(k,i,lat)+qo2p(k+1,i,lat))+xnmbarm(k,i)*
     |      (rk1(k,i,lat)*op_upd(k,i)+rk6*nplus(k,i))*o2(k,i)*
     |       rmassinv_o2+rk26*xiop2d(k,i)*o2(k,i)*rmassinv_o2
!
! C = K4*N(N4S)+K5*N(NO) (s8)
!
          c(k,i) = xnmbarm(k,i)*(rk4*n4s(k,i)*rmassinv_n4s+
     |       rk5*no(k,i)*rmassinv_no)
!
! D = QI(N2+) (s7)
!
          d(k,i) = .5*(qn2p(k,i,lat)+qn2p(k+1,i,lat))+(rk16*xiop2p(k,i)+
     |      rk23*xiop2d(k,i))*xn2(k,i)*rmassinv_n2
!
! E =K3*N(O)+K9*N(O2) (s6)
!
          e(k,i) = xnmbarm(k,i)*(rk3(k,i,lat)*o1(k,i)*rmassinv_o1+
     |       rk9*o2(k,i)*rmassinv_o2)
!
! F+G = N(O+)+N(N+) (s5)
! (very small "diamond diffs" with tgcm15 due to op_upd)
!
          fg(k,i) = op_upd(k,i)+nplus(k,i)
!
! H = K9*N(02)
!
          h(k,i) = xnmbarm(k,i)*rk9*o2(k,i)*rmassinv_o2
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
     |      4.
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
!     call addfsech('XNMBARM' ,' ',' ',xnmbarm,i0,i1,nk,nkm1,lat)
!     call addfsech('NPLUS'   ,' ',' ',nplus(:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('A_COEF'  ,' ',' ',a      ,i0,i1,nk,nkm1,lat)
!     call addfsech('B_COEF'  ,' ',' ',b      ,i0,i1,nk,nkm1,lat)
!     call addfsech('C_COEF'  ,' ',' ',c      ,i0,i1,nk,nkm1,lat)
!     call addfsech('D_COEF'  ,' ',' ',d      ,i0,i1,nk,nkm1,lat)
!     call addfsech('E_COEF'  ,' ',' ',e      ,i0,i1,nk,nkm1,lat)
!     call addfsech('FG_COEF' ,' ',' ',fg     ,i0,i1,nk,nkm1,lat)
!     call addfsech('H_COEF'  ,' ',' ',h      ,i0,i1,nk,nkm1,lat)
!     call addfsech('A0'      ,' ',' ',a0     ,i0,i1,nk,nkm1,lat)
!     call addfsech('A1'      ,' ',' ',a1     ,i0,i1,nk,nkm1,lat)
!     call addfsech('A2'      ,' ',' ',a2     ,i0,i1,nk,nkm1,lat)
!     call addfsech('A3'      ,' ',' ',a3     ,i0,i1,nk,nkm1,lat)
!     call addfsech('A4'      ,' ',' ',a4     ,i0,i1,nk,nkm1,lat)
!
! Solve quartic. Vquart returns electron density Ne in root:
!
      call vquart(a0,a1,a2,a3,a4,root,lev0,lev1,lon0,lon1,lat)

      call addfsech('ROOT',' ',' ',root,i0,i1,nk,nkm1,lat)
!
! Insure positive Ne (at least 3100):
!     where(root < 3.1e3) root = 3.1e3
      do i=lon0,lon1
        do k=lev0,lev1-1
          if (root(k,i) < 3.1e3) root(k,i) = 3.1e3
        enddo 
      enddo
!
! Calculate N2+, O2+, NO+ (cm3)
!
      do i=lon0,lon1      
        do k=lev0,lev1-1
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
! Small diffs in O2P reported as "global diffs" in tgcmproc, but
! not seen in any plots.
!     call addfsech('NPLUSb' ,' ',' ',nplus(:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('N2P_ELD',' ',' ',n2p  (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('O2P_ELD',' ',' ',o2p  (:,i0:i1),i0,i1,nk,nkm1,lat)
!     call addfsech('NOP_ELD',' ',' ',nop  (:,i0:i1),i0,i1,nk,nkm1,lat)
!
! Transfer root to electrons output array:
      do i=lon0,lon1
        do k=lev0,lev1-2
          electrons(k+1,i) = sqrt(root(k,i)*root(k+1,i))
        enddo ! k=lev0,lev1-2
!
! Lower and upper boundaries:
        electrons(lev0,i) = sqrt(root(lev0  ,i)**3/root(lev0+1,i))
        electrons(lev1,i) = sqrt(root(lev1-1,i)**3/root(lev1-2,i))
      enddo ! i=lon0,lon1

!     call addfsech('NE_ELDEN',' ',' ',electrons(:,i0:i1),
!    |  i0,i1,nk,nk,lat)

      end subroutine elden
!-----------------------------------------------------------------------
      subroutine vquart(a0,a1,a2,a3,a4,root,lev0,lev1,lon0,lon1,lat)
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
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0:lon1),intent(in)  :: a0,a1,a2,a3,a4
      real,dimension(lev0:lev1,lon0:lon1),intent(out) :: root
!
! Local:
      integer :: k,i,nlevs
      real,dimension(lev0:lev1,lon0:lon1) :: w1,w2,w3 ! work arrays
      real,parameter :: e=1.e-300 ! largest exponent on ieee is about 307
!
      nlevs = lev1-lev0+1
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! w1 = ch
          w1(k,i) = -(a4(k,i)*a0(k,i)-4.*a3(k,i)*a1(k,i)+3.*a2(k,i)**2)/
     |                12.
!
! w2 = cg
          w2(k,i) = (a4(k,i)*(a2(k,i)*a0(k,i)-a1(k,i)**2)-a3(k,i)*
     |      (a3(k,i)*a0(k,i)-a1(k,i)*a2(k,i))+a2(k,i)*(a3(k,i)*a1(k,i)-
     |      a2(k,i)**2))/4.
!
! root=rlam=-2.*real((.5*(cmplx(cg,0.)+csqrt(cmplx(cg**2+4.
!      *ch**3+e,0.)))+cmplx(e,0.))**(1./3.))
! 
          root(k,i) = -2.*real((.5*(cmplx(w2(k,i),0.)+
     |      csqrt(cmplx(w2(k,i)**2+4.*w1(k,i)**3+e,0.)))+
     |      cmplx(e,0.))**(1./3.))
!
! W1=P=SQRT(A(5)*RLAM+A(4)**2-A(5)*A(3)+E)
!
          w1(k,i) = a4(k,i)*root(k,i)+a3(k,i)**2-a4(k,i)*a2(k,i)+e
          if (w1(k,i) < 0.) w1(k,i) = 0.
          w1(k,i) = sqrt(w1(k,i))
!
! W2=Q=SQRT((2.*RLAM+A(3))**2-A(5)*A(1)+E)
!
          w2(k,i) = sqrt((2.*root(k,i)+a2(k,i))**2-a4(k,i)*a0(k,i)+e)
!
! W3=PQ=2.*A(4)*RLAM+A(4)*A(3)-A(5)*A(2)+E
!
          w3(k,i) = 2.*a3(k,i)*root(k,i)+a3(k,i)*a2(k,i)-a4(k,i)*a1(k,i)
     |      +e
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
          root(k,i) = (w3(k,i)+sqrt(w3(k,i)**2-a4(k,i)*(a2(k,i)+2.*
     |      root(k,i)-w2(k,i))))/a4(k,i)
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfsech('W1'    ,' ',' ',w1,lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('W2'    ,' ',' ',w2,lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('W3'    ,' ',' ',w3,lon0,lon1,nlevs,nlevs-1,lat)
!     call addfsech('VQROOT',' ',' ',root,lon0,lon1,nlevs,nlevs-1,lat)

      end subroutine vquart

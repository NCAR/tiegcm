!
      subroutine qjion(tn,o2,o1,n2,o2p,op,n4s,n2d,no,ne,barm,
     |  n2p,nplus,nop,xiop2p,xiop2d, lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculate ion chemistry contribution to neutral gas heating
!   and O2 dissociation. 
! This routine modifies Q(k,i) from qrj module (all args are input).
!
      use cons_module,only: avo,evergs,rmassinv_o2,rmassinv_o1,
     |  rmassinv_n2,rmassinv_no,rmassinv_n2d,rmassinv_n4s,p0,boltz,
     |  expz,expzmid,expzmid_inv
      use chemrates_module,only: 
     |  rk1,rk2,rk3,rk4,rk5,rk6,rk7,rk8,rk9,rk10,ra1,ra2,ra3,rk16,rk17,
     |  rk18,rk19,rk20,rk21,rk22,rk23,rk24,rk25,rk26,rk27
      use qrj_module,only: ! Q is modified, all others are input.
     |  qtotal,! total heating 
     |  qop2p, ! o+(2p)
     |  qop2d, ! o+(2d)
     |  qo2p,  ! o2+ ionization
     |  qop,   ! o+  ionization
     |  qn2p,  ! n2+ ionization
     |  qnp,   ! n+  ionization
     |  qnop   ! no+ ionization
      use addfld_module,only: addfld
      implicit none
!
! Args (all input):
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) :: 
     |  tn,      ! neutral temperature (deg K)
     |  o2,      ! molecular oxygen (mmr)
     |  o1,      ! atomic oxygen (mmr)
     |  n2,      ! molecular nitrogen (mmr)
     |  o2p,     ! O2+ ion
     |  op,      ! O+ ion
     |  n4s,     ! N(4S)
     |  n2d,     ! N(2D) (updated from comp_n2d)
     |  no,      ! nitric oxide
     |  ne,      ! electron density (cm3)
     |  barm   , ! p0*e(-z)*barm/kT
     |  n2p,     ! N2+ (from elden)
     |  nplus,   ! N+  (from elden)
     |  nop,     ! NO+ (from elden)
     |  xiop2p,  ! from oplus
     |  xiop2d   ! from oplus
!
! Local:
      integer :: k,i
      real,dimension(lev0:lev1,lon0:lon1) :: 
     |  qtot,    ! total ionization rate                 (s1)
     |  qphoto,  ! photo-electron heating of neutral gas (s2)
     |  qic,     ! ion chemistry heating of neutral gas  (s3)
     |  xnmbarm, ! p0*e(-z)*barm/kT ("(K)" -- averaged)        (s12)
     |  xnmbari  ! p0*e(-z)*barm/kT ("(K+1/2)" -- not averaged (s11)
      real,dimension(lev0:lev1) :: aureff
!
      aureff(:) = 0.05 ! this is local (not the aureff from aurora module).
!
      do i=lon0,lon1
        do k=lev0,lev1-1
          xnmbarm(k,i) = p0*expz(k)*.5*(barm(k,i)+barm(k+1,i))/    ! s12 (K+1/2)
     |      (boltz*tn(k,i))
        enddo ! k=lev0,lev1
        do k=lev0+1,lev1-1
          xnmbari(k,i) = p0*expzmid_inv*expz(k)*barm(k,i)/         ! s11 (K)
     |      (boltz*.5*(tn(k,i)+tn(k-1,i)))
        enddo ! k=lev0+1,lev1-1
        xnmbari(lev0,i) = p0*expzmid_inv*expz(lev0)*barm(lev0,i)/  ! s11 bottom
     |    (boltz*.5*(3.*tn(lev0,i)-tn(lev0+1,i)))
        xnmbari(lev1,i) = p0*expzmid*expz(lev1-1)*barm(lev1,i)/    ! s11 top
     |    (boltz*.5*(3.*tn(lev1-1,i)-tn(lev1-2,i)))
      enddo ! i=lon0,lon1
!
! xnmbari = xnmbar at interfaces (averaged in the column == S11 (K))
! xnmbarm = xnmbar at interfaces (NOT averaged in the column == S12 (K+1/2))
!
!     call addfld('XNMBARM' ,' ',' ',xnmbarm,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('XNMBARI' ,' ',' ',xnmbari,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! qtot = total ionization rate = sum(Qxx) = 
! (QO2+) + (QO+) + (QN2+) + (QNO+) + (QN+) + (QO+(2D)) + (QO+(2P))
!
      qtot = 0. ! whole array init
      do i=lon0,lon1
        do k=lev0,lev1
          qtot(k,i) = qtot(k,i)+(qo2p(k,i,lat)+qop(k,i,lat)+
     |      qn2p(k,i,lat)+qnop(k,i,lat)+qnp(k,i,lat)+qop2d(k,i,lat)+
     |      qop2p(k,i,lat))
        enddo ! k=lev0,lev1
      enddo ! i=lon0,lon1

!     call addfld('QTOT',' ',' ',qtot,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      do i=lon0,lon1
        do k=lev0,lev1
          qphoto(k,i) = qtot(k,i)*aureff(k)*35.*avo*evergs/xnmbari(k,i)
        enddo ! k=lev0,lev1
 
        do k=lev0,lev1-1
          qic(k,i) = 
     |      (avo*(o2(k,i)*rmassinv_o2*(rk1(k,i,lat)*op(k,i)*1.555+
     |      (rk6*2.486+rk7*6.699)*nplus(k,i)+
     |      rk9*n2p(k,i)*3.52)+
     |      op(k,i)*(rk2(k,i,lat)*n2(k,i)*rmassinv_n2*1.0888+
     |      rk10*n2d(k,i)*rmassinv_n2d*1.45)+
     |      o1(k,i)*rmassinv_o1*(rk3(k,i,lat)*n2p(k,i)*0.70+
     |      rk8*nplus(k,i)*0.98)+
     |      o2p(k,i)*(rk4*n4s(k,i)*rmassinv_n4s*4.21+
     |      rk5*no(k,i)*rmassinv_no*2.813))+
     |      .5*(ne(k,i)+ne(k+1,i))*(ra1(k,i,lat)*nop(k,i)*0.854+
     |       ra2(k,i,lat)*o2p(k,i)*5.2755+
     |       ra3(k,i,lat)*n2p(k,i)*3.678)/xnmbarm(k,i))*evergs+
!
     |      (avo*(((rk16*3.02+rk17*0.7)*n2(k,i)*rmassinv_n2+
     |       rk18*o1(k,i)*rmassinv_o1*5.0)*xiop2p(k,i)+
     |      (rk23*n2(k,i)*rmassinv_n2*1.33+
     |       rk24*o1(k,i)*rmassinv_o1*3.31+
     |       rk26*4.87*o2(k,i)*rmassinv_o2)*xiop2d(k,i))+
     |      (.5*(ne(k,i)+ne(k+1,i))*((rk19(k,i,lat)*5.0+
     |       rk20(k,i,lat)*1.69)*xiop2p(k,i)+rk25(k,i,lat)*3.31*
     |       xiop2d(k,i))-(rk21*5.02+rk22*1.69)*xiop2p(k,i)-
     |       rk27*3.33*xiop2d(k,i))/xnmbarm(k,i))*evergs
!
! Insure qic > 0:
          if (qic(k,i) < 1.e-30) qic(k,i) = 1.e-30
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     call addfld('QPHOTO' ,' ',' ',qphoto,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!     call addfld('QIC'    ,' ',' ',qic   ,
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)
!
! Add qphoto and qic to Q (reuse qtot = qphoto+qic):
      do i=lon0,lon1
        do k=lev0,lev1-2
          qtot(k+1,i) = qphoto(k+1,i)+sqrt(qic(k,i)*qic(k+1,i)) 
          qtotal(k+1,i,lat) = qtotal(k+1,i,lat)+qtot(k+1,i)
        enddo ! k=lev0,lev1-1
!
! Bottom and top boundaries:
        qtot(lev0,i) = qphoto(lev0,i)+sqrt(qic(lev0,i)**3/qic(lev0+1,i))
        qtot(lev1,i) = qphoto(lev1,i)+sqrt(qic(lev1-1,i)**3/
     |    qic(lev1-2,i))
        qtotal(lev0,i,lat) = qtotal(lev0,i,lat)+qtot(lev0,i)
        qtotal(lev1,i,lat) = qtotal(lev1,i,lat)+qtot(lev1,i)
      enddo ! i=lon0,lon1

!     call addfld('QTOTAL' ,' ',' ',qtotal(:,:,lat),
!    |  'lev',lev0,lev1,'lon',lon0,lon1,lat)

      end subroutine qjion

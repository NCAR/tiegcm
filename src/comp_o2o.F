!
      subroutine comp_o2o(tn,o2,o1,n2,he,barm,op,no,n4s,n2d,o2p,ne,
     |  n2p,nplus,nop,lev0,lev1,lon0,lon1,lat)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! Calculates fs array which gives sources and sinks for comp (o2,o). 
!
      use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,
     |  rmassinv_he,p0,expz,boltz,rmass_no,rmass_n4s,rmass_n2d,rmass_o2,
     |  rmass_o1
      use qrj_module,only: rj,qo2p,qop
      use chemrates_module,only: rk4,rk5,rk6,rk7,rk8,rk9,rk10,beta2,
     |  beta6,ra1,ra2,beta3,beta8,rk1,beta1,rkm12,rk3,fs
      use addfld_module,only: addfld
      implicit none
!
! Args:
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat
!
! All input arrays are at current time step (itp): 
      real,dimension(lev0:lev1,lon0-2:lon1+2),intent(in) ::
     |  tn,    ! neutral temperature (deg K) 
     |  o2,    ! O2 (mmr)
     |  o1,    ! O1 (mmr)
     |  n2,    ! N2 (mmr)
     |  he,    ! He (mmr)
     |  barm,  ! mean molecular weight
     |  op,    ! O+
     |  no,    ! nitric oxide
     |  n4s,   ! N(4s)
     |  n2d,   ! N(2d)
     |  o2p,   ! O2+
     |  ne,    ! electron density
     |  n2p,   ! N(2p)
     |  nplus, ! N+
     |  nop    ! NO+
!
! Local:
      integer :: k,i,nk
      real,dimension(lev0:lev1,lon0:lon1) ::
     |  rji,            ! RJ at interfaces (qrj module)
     |  qo2pi,          ! O2+ ionization at interfaces (qrj module)
     |  qopi,           ! O+ ionization at interfaces
     |  xnmbar,         ! p0*e(-z)*barm
     |  pox1,pox2,      ! OX production terms
     |  lox1,lox2,lox3, ! OX loss terms
     |  po21,po22,po23, ! O2 production terms
     |  lo21,lo22       ! O2 loss terms
!
      nk = lev1-lev0+1

!     call addfld('TN_O2O',' ',' ',tn(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('O2_O2O',' ',' ',o2(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('O1_O2O',' ',' ',o1(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('N2_O2O',' ',' ',n2(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('HE_O2O',' ',' ',he(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('BARM_O2O',' ',' ',barm(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('OP_O2O',' ',' ',op(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('NO_O2O',' ',' ',no(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('N4S_O2O',' ',' ',n4s(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('N2D_O2O',' ',' ',n2d(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('O2P_O2O',' ',' ',o2p(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('NE_O2O',' ',' ',ne(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!
! Diffs with r1014 in N2P, NPLUS, NOP (see elden.F)
!     call addfld('N2P_O2O',' ',' ',n2p(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('NPLUS_O2O',' ',' ',nplus(:,lon0:lon1),'lev',lev0,lev1
!    |  ,'lon',lon0,lon1,lat)
!     call addfld('NOP_O2O',' ',' ',nop(:,lon0:lon1),'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!
! Loop over subdomain at current latitude:
      do i=lon0,lon1
        do k=lev0,lev1-1
!
! Qrj ionization rates at interfaces:
          rji  (k,i) = 0.5*(rj(k,i,lat)+rj(k+1,i,lat))
          qo2pi(k,i) = 0.5*(qo2p(k,i,lat)+qo2p(k+1,i,lat))
          qopi (k,i) = 0.5*(qop(k,i,lat)+qop(k+1,i,lat))
!
! n*mbar (k+1/2):
          xnmbar(k,i) = expz(k)*p0/(boltz*tn(k,i)*(o2(k,i)*rmassinv_o2+
     |     o1(k,i)*rmassinv_o1+he(k,i)*rmassinv_he+n2(k,i)*rmassinv_n2))
!
! OX production:
! s1
          pox1(k,i) = xnmbar(k,i)**2*
     |      (beta3(k,i,lat)*n4s(k,i)/rmass_n4s*no(k,i)/rmass_no+
     |      beta6*n2d(k,i)/rmass_n2d*no(k,i)/rmass_no)+
     |      0.5*(beta8(k,i,lat)+beta8(k+1,i,lat))*
     |      no(k,i)/rmass_no*xnmbar(k,i)+xnmbar(k,i)*
     |      (rk4*o2p(k,i)*n4s(k,i)/rmass_n4s+
     |      rk10*op(k,i)*n2d(k,i)/rmass_n2d)+
     |      (ra1(k,i,lat)*nop(k,i)+2.*ra2(k,i,lat)*o2p(k,i))*
     |      sqrt(ne(k,i)*ne(k+1,i))
! s2
          pox2(k,i) = xnmbar(k,i)*(beta1(k,i,lat)*n4s(k,i)/rmass_n4s+
     |      beta2*n2d(k,i)/rmass_n2d)+rk1(k,i,lat)*op(k,i)+rk7*
     |      nplus(k,i)+2.*rji(k,i)
!
! OX loss:
          lox1(k,i) = 2.*rkm12(k,i,lat)*xnmbar(k,i)/
     |      (.5*(barm(k,i)+barm(k+1,i)))
          lox2(k,i) = rk3(k,i,lat)*n2p(k,i)+rk8*nplus(k,i)
          lox3(k,i) = qopi(k,i)
!
! O2 production:
          po21(k,i) = rkm12(k,i,lat)*xnmbar(k,i)/
     |      (.5*(barm(k,i)+barm(k+1,i)))
!
          po22(k,i) = 0.
          po23(k,i) = rk5*no(k,i)/rmass_no*o2p(k,i)*xnmbar(k,i)
!
! O2 loss:
          lo21(k,i) = xnmbar(k,i)*(beta1(k,i,lat)*n4s(k,i)/rmass_n4s+
     |      beta2*n2d(k,i)/rmass_n2d)+rk1(k,i,lat)*op(k,i)+(rk6+rk7)*
     |      nplus(k,i)+rk9*n2p(k,i)+rji(k,i)
!
          lo22(k,i) = qo2pi(k,i)
!
! Matrix coefficients for O-O2-N2 solution:
          fs(i,k,1,1,lat) = -lo21(k,i)
          fs(i,k,1,2,lat) = xnmbar(k,i)*po21(k,i)*o1(k,i)*rmassinv_o1*
     |      rmass_o2*rmassinv_o1
          fs(i,k,2,1,lat) = pox2(k,i)*rmass_o1*rmassinv_o2
          fs(i,k,2,2,lat) = -lox2(k,i)-lox1(k,i)*o1(k,i)*rmassinv_o1*
     |      xnmbar(k,i)
          fs(i,k,1,0,lat) = (po23(k,i)-lo22(k,i))*rmass_o2/xnmbar(k,i)
          fs(i,k,2,0,lat) = (pox1(k,i)-lox3(k,i))*rmass_o1/xnmbar(k,i)
! Helium prod/loss set to zero for now
          fs(i,k,3,0,lat) = 0.
          fs(i,k,1,3,lat) = 0.
          fs(i,k,2,3,lat) = 0.
          fs(i,k,3,1,lat) = 0.
          fs(i,k,3,2,lat) = 0.
          fs(i,k,3,3,lat) = 0.
        enddo ! k=lev0,lev1-1
      enddo ! i=lon0,lon1

!     do k=lev0,lev1-1
!       write(6,"('comp_o2o: lat=',i3,' k=',i3)") lat,k
!       write(6,"('fs(:,k,1,1,lat)=',/,(6e12.4))") fs(:,k,1,1,lat)
!       write(6,"('fs(:,k,1,2,lat)=',/,(6e12.4))") fs(:,k,1,2,lat)
!       write(6,"('fs(:,k,2,1,lat)=',/,(6e12.4))") fs(:,k,2,1,lat)
!       write(6,"('fs(:,k,2,2,lat)=',/,(6e12.4))") fs(:,k,2,2,lat)
!       write(6,"('fs(:,k,1,0,lat)=',/,(6e12.4))") fs(:,k,1,0,lat)
!       write(6,"('fs(:,k,2,0,lat)=',/,(6e12.4))") fs(:,k,2,0,lat)
!     enddo ! k=lev0,lev1-1

!     call addfld('MBAR_O2O',' ',' ',xnmbar,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('POX1',' ',' ',pox1,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('POX2',' ',' ',pox2,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('LOX1',' ',' ',lox1,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('LOX2',' ',' ',lox2,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('LOX3',' ',' ',lox3,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('PO21',' ',' ',po21,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('PO22',' ',' ',po22,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('PO23',' ',' ',po23,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('LO21',' ',' ',lo21,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)
!     call addfld('LO22',' ',' ',lo22,'lev',lev0,lev1,
!    |  'lon',lon0,lon1,lat)

      end subroutine comp_o2o

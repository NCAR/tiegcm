!#include "dims.h"
c
c     file mud2cr.f
c
C  am 01/03/22 new version mudpack5.0:  
c  am 01/04/16 mudpack modified:-residual calculated with unmodified stencils
c                             converged to same solution as with direct solver
c                            -the number of restriction and prolongation      
c				mgopt(2) = 3
c       			mgopt(3) = 2 
c                             has to be increased in order to get convergence
c                             another possibility would be to use a
c			      preconditioner (so far not necessary and therefore
c			      not implemented)     
c                            -unmodified stencil saved in: 
c                             common/mudmd/cofum(1) 
c 
C
c                       MUDPACK version 5.0 
!
! 5/02 B. Foster:
! Use-associate coefficients cee and cofum from dynamo module in tiegcm1 
! (dynamo.F).
!
      subroutine mudmod(pe,jntl,isolve,ier)
      use dynamo_module,only: cee
      implicit none
      integer jntl,ier  ! output: not converged ier < 0
      integer,intent(in) :: isolve
c
c     set grid size params
c
      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork
      parameter (iixp = 5 , jjyq = 3, iiex = 5, jjey = 5 )
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
!
!     integer imx0, jmx0, imx1, jmx1,imx2,jmx2,imx3,jmx3,imx4,jmx4
!     PARAMETER(IMX0=NNX,JMX0=NNY,IMX1=(IMX0+1)/2,JMX1=(JMX0+1)/2,
!    1  IMX2=(IMX1+1)/2,JMX2=(JMX1+1)/2,IMX3=(IMX2+1)/2,
!    2  JMX3=(JMX2+1)/2,IMX4=(IMX3+1)/2,JMX4=(JMX3+1)/2)
!     integer nc0,nc1,nc2,nc3,nc4,ncee,nny2
!     PARAMETER (NC0=1,NC1=NC0+10*(IMX0*JMX0),NC2=NC1+9*(IMX1*JMX1),
!    1  NC3=NC2+9*(IMX2*JMX2),NC4=NC3+9*(IMX3*JMX3),
!    2  NCEE=NC4+9*(IMX4*JMX4)-1)
!     PARAMETER (NNY2=2*NNY+1)
!#include "ceee.h"
c
c     estimate work space for point relaxation (see mud2cr.d)
c
      parameter (llwork=(7*(nnx+2)*(nny+2)+76*nnx*nny)/3 )
      real phi(nnx,nny),rhs(nnx,nny),work(llwork)
c
c     put integer and floating point argument names in contiguous
c     storage for labelling in vectors iprm,fprm
c
      integer iprm(17),mgopt(4)
      real fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmud2cr/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      real PE(NNX,1)
      integer maxcya
      DATA MAXCYA/50/
      integer mm,nn,jj,jjj
      real pi

!     write(6,"('Enter mudmod: jntl=',i3)") jntl
c
c     set input integer arguments
c
      MM = NNX
      NN = NNY
      PI = 4.*ATAN(1.)
C
C     SET INPUT INTEGER PARAMETERS
C
      INTL = JNTL
c
c     set boundary condition flags
c
      nxa = 0
      nxb = 0
      nyc = 2
      nyd = 1
c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = jjyq
      iex = iiex
      jey = jjey
      nx = nnx
      ny = nny
c
c     set multigrid arguments (w(2,1) cycling with fully weighted
c     residual restriction and cubic prolongation)
c
      mgopt(1) = 2
      mgopt(2) = 3
      mgopt(3) = 2
      mgopt(4) = 3
c
c     set for one cycle
c
      maxcy = maxcya
c
c     set no initial guess forcing full multigrid cycling
c
      iguess = 0
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set line z relaxation
c
      method = 3
c
c     set end points of solution rectangle in (x,y) space
c
      xa = -pi
      xb =  pi
      yc = 0.0
      yd = 0.5*pi
c
c     set error control flag
c
      tolmax = 0.01
c
c     set right hand side in rhs
c     initialize phi to zero
c
      do i=1,nx
        do j=1,ny

!         write(6,"('mudmod: j=',i2,' i=',i2,' cee=',e12.4)")
!    |      j,i,CEE(I+(J-1)*NX+9*NX*NY)

          RHS(I,J) = CEE(I+(J-1)*NX+9*NX*NY)
          phi(i,j) = 0.0
        end do
      end do
c
c     set specified boundaries in phi
c
      DO I=1,NX
        PHI(I,NY) = RHS(I,NY)/CEE(I+(NY-1)*NX+8*NX*NY)
      END DO
      
!     write(*,100)
  100 format(//' mud2cr test ')
!     write (*,101) (iprm(i),i=1,15)
  101 format(/,' integer input arguments ',/,
     |  ' intl =  ',i2,/,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,
     |  ' nyd =   ',i2,/,' ixp = ',i2,' jyq = ',i2,' iex = ',i2,
     |  ' jey =   ',i2,/,' nx =  ',i3,' ny =  ',i3,' iguess = ',i2,
     |  ' maxcy = ',i3,/,' method = ',i2, ' work space estimate = ',i7)
!     write (*,102) (mgopt(i),i=1,4)
  102 format(/' multigrid option arguments ',
     |  /,' kcycle = ',i2,
     |  /,' iprer = ',i2,
     |  /,' ipost = ',i2
     |  /,' intpol = ',i2)
!     write(*,103) xa,xb,yc,yd,tolmax
  103 format(/' floating point input parameters ',
     |  /,' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3,
     |  /,' tolerance (error control) =   ',e10.3)
!     write(6,"('fprm(1-5) (xa,xb,yc,yd,tolmax=',6f8.3)") fprm(1:5)
c
c     intialization call
c
!     write(*,104) intl
  104 format(/' discretization call to mud2cr', ' intl = ', i2)
      call mud2cm(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)
!     write (*,200) ierror,iprm(16)
  200 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     attempt solution
c
      intl = 1
!     write(*,106) intl,method,iguess
  106 format(/' approximation call to mud2cr',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      
      call mud2cm(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)
      ier = ierror ! ier < 0 not converged
      if(ier < 0 )  goto 108
      
!     write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
C
C     COPY PHI TO PE
C
      DO J = 1,NY
	JJ = NY+J-1
	JJJ = NY+1-J
	DO I = 1,NX
	  PE(I,JJ) = PHI(I,J)
	  PE(I,JJJ) = PHI(I,J)
	END DO
      END DO
      
  108 continue   
      end
C-------------------------------------------------------------------
      subroutine mud2cm(iparm,fparm,work,rhs,phi,mgopt,
     +                  ierror,isolve)
      implicit none
      integer,intent(in) :: isolve
      integer iparm,mgopt,ierror
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm,xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,iw,k,kb,nx,ny,ic,itx,ity
      dimension iparm(17),fparm(6),mgopt(4)
      real work(*),phi(*),rhs(*)
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      
      data int / 0 /
      save int
!     write(6,"('Enter mud2cm')")

      ierror = 1
      intl = iparm(1)    ! set and check intl on all calls
      if (intl*(intl-1).ne.0) return
      if (int.eq.0) then
	int = 1
	if (intl.ne.0) return  ! very first call is not intl=0
      end if
      ierror = 0
c
c     set  arguments internally
c     these will not be rechecked if intl=1!
c
      nxa = iparm(2)
      nxb = iparm(3)
      nyc = iparm(4)
      nyd = iparm(5)
      ixp = iparm(6)
      jyq = iparm(7)
      iex = iparm(8)
      jey = iparm(9)
      ngrid = max0(iex,jey)
      nfx = iparm(10)
      nfy = iparm(11)
      iguess = iparm(12)
      maxcy = iparm(13)
      method = iparm(14)
      nwork = iparm(15)
      kcycle = mgopt(1)
      if (kcycle .eq. 0) then
c       set defaults
	kcycle = 2
	iprer = 2
	ipost = 1
	intpol = 3
      else
	iprer = mgopt(2)
	ipost = mgopt(3)
	intpol = mgopt(4)
      end if
      xa = fparm(1)
      xb = fparm(2)
      yc = fparm(3)
      yd = fparm(4)
      tolmax = fparm(5)
      if (intl .eq. 0) then  ! intialization call
c
c     check input arguments
c
	ierror = 2   ! check boundary condition flags
	if (max0(nxa,nxb,nyc,nyd).gt.2) return
	if (min0(nxa,nxb,nyc,nyd).lt.0) return
	if (nxa.eq.0.and.nxb.ne.0) return
	if (nxa.ne.0.and.nxb.eq.0) return
	if (nyc.eq.0.and.nyd.ne.0) return
	if (nyc.ne.0.and.nyd.eq.0) return
	ierror = 3   ! check grid sizes
	if (ixp.lt.2) return
	if (jyq.lt.2) return
	ierror = 4
	ngrid = max0(iex,jey)
	if (iex.lt.1) return
	if (jey.lt.1) return
	if (ngrid.gt.50) return
	ierror = 5
	if (nfx.ne.ixp*2**(iex-1)+1) return
	if (nfy.ne.jyq*2**(jey-1)+1) return
	ierror = 6
	if (iguess*(iguess-1).ne.0) return
	ierror = 7
	if (maxcy.lt.1) return
	ierror = 8
	if (method.lt.0 .or. method.gt.3) return
	ierror = 9
c       compute and test minimum work space
	isx = 0
	if (method.eq.1 .or. method.eq.3) then
	  if (nxa.ne.0) isx = 3
	  if (nxa.eq.0) isx = 5
	end if
	jsy = 0
	if (method.eq.2 .or. method.eq.3) then
	  if (nyc.ne.0) jsy = 3
	  if (nyc.eq.0) jsy = 5
	end if
	kps = 1
	do k=1,ngrid
c       set subgrid sizes
	  nxk(k) = ixp*2**(max0(k+iex-ngrid,1)-1)+1
	  nyk(k) = jyq*2**(max0(k+jey-ngrid,1)-1)+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kps = kps+(nx+2)*(ny+2)+nx*ny*(10+isx+jsy)
	end do
	iparm(16) = kps+(nfx+2)*(nfy+2)   ! exact minimum work space
	lwork = iparm(16)
	if (lwork .gt. nwork) return
	ierror = 10   ! check solution region
	if (xb.le.xa .or. yd.le.yc) return
	ierror = 11
	if (tolmax .lt. 0.0) return
	ierror = 12   ! multigrid parameters
	if (kcycle.lt.0) return
	if (min0(iprer,ipost).lt.1) return
	if ((intpol-1)*(intpol-3).ne.0) return
	if (max0(kcycle,iprer,ipost).gt.2) then
	  ierror = -5   ! inefficient multigrid cycling
	end if
	if (ierror .gt. 0) ierror = 0   ! no fatal errors
c
c     set work space pointers and discretize pde at each grid level
c
	iw = 1
	do kb=1,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k)
	  ny = nyk(k)
	  kpbgn(k) = iw
	  kcbgn(k) = kpbgn(k)+(nx+2)*(ny+2)
	  ktxbgn(k) = kcbgn(k)+10*nx*ny
	  ktybgn(k) = ktxbgn(k)+isx*nx*ny
	  iw = ktybgn(k)+jsy*nx*ny
	  ic = kcbgn(k)
	  itx = ktxbgn(k)
	  ity = ktybgn(k)
	  klevel = k
	  call dismd2cr(nx,ny,work(ic),work(itx),work(ity),
     +                  work,ierror,isolve)
	  end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call mud2c1m(nx,ny,rhs,phi,work)
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) then
	
!	   ierror = -1   ! flag convergenc failure
	   write(6,*) "no convergence with mudmod"
!         
           iguess = 1 
           iparm(12)= iguess                           	   
           call mud2cr1(nx,ny,rhs,phi,work) !  solve with modified stencils
	   
	   fparm(6) = relmax
	   if (relmax.gt.tolmax) then
	     write(6,*) "no convergence with mud"
	     ierror = -1   ! flag convergenc failure
	   end if
	   
	end if
      end if
      
      return
      end
c------------------------------------------------------------------------      
      subroutine mud2c1m(nx,ny,rhsf,phif,wk)
      implicit none
      integer nx,ny
      real phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ic,ir,ipc,irc,icc
      integer ncx,ncy,jj,ij,i,j,iter
      integer iw,itx,ity,ierror
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy

!     write(6,"('Enter mud2c1m')")
      
      nx = nxk(ngrid)
      ny = nyk(ngrid)
      ip = kpbgn(ngrid)
      ic = kcbgn(ngrid)
      ir = ic+9*nx*ny
c
c     set phif,rhsf in wk and adjust right hand side
c
      call swk2(nx,ny,phif,rhsf,wk(ip),wk(ir))
      if (iguess.eq.0) then
c
c     no initial guess at finest grid level!
c
	do kb=2,ngrid
	  k = ngrid-kb+1
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ir = kcbgn(k+1)+9*nx*ny
	  ncx = nxk(k)
	  ncy = nyk(k)
	  ipc = kpbgn(k)
	  icc = kcbgn(k)
	  irc = icc+9*ncx*ncy
c
c     transfer down to all grid levels
c
	  call trsfc2(nx,ny,wk(ip),wk(ir),ncx,ncy,
     +                wk(ipc),wk(irc))
	end do
c
c     adjust right hand side at all grid levels in case
c     rhs or specified b.c. in phi or gbdy changed
c
	do k=1,ngrid
	  nx = nxk(k)
	  ny = nyk(k)
	  ip = kpbgn(k)
	  ic = kcbgn(k)
	  call adjmd2cr(nx,ny,wk(ip),wk(ic))
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcymd2cr(wk)
	  nx = nxk(k+1)
	  ny = nyk(k+1)
	  ip = kpbgn(k+1)
	  ipc = kpbgn(k)
	  ncx = nxk(k)
	  ncy = nyk(k)
c
c     lift or prolong approximation from k to k+1
c
	  call prolon2(ncx,ncy,wk(ipc),nx,ny,wk(ip),nxa,nxb,
     +                 nyc,nyd,intpol)
	end do
      else
c
c     adjust rhs at finest grid level only
c
	nx = nxk(ngrid)
	ny = nyk(ngrid)
	ip = kpbgn(ngrid)
	ic = kcbgn(ngrid)
	call adjmd2cr(nx,ny,wk(ip),wk(ic))
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcym2cm(wk)
	if (tolmax.gt.0.0) then
c
c      error control
c
	  relmax = 0.0
	  phmax = 0.0
	  
	  do j=1,nfy
	    jj = j*(nfx+2)
	    do i=1,nfx
	      ij = jj+i+1
	      phmax = max(phmax,abs(wk(ij)))
	      relmax = max(relmax,abs(wk(ij)-phif(i,j)))
	      
	      phif(i,j) = wk(ij)
	    end do
	  end do
c
c     set maximum relative difference and check for convergence
c
	  if (phmax.gt.0.0) relmax = relmax/phmax
	  if (relmax.le.tolmax) return
	end if
      end do
c
c     set final interate after maxcy cycles in phif
c
      do j=1,nfy
	jj = j*(nfx+2)
	do i=1,nfx
	  ij = jj+i+1
	  phif(i,j) = wk(ij)
	end do
      end do
      return
      end

c------------------------------------------------------------------------
      subroutine kcym2cm(wk)
      use dynamo_module,only: cofum
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      real wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ip,ic,ipc,irc,itx,ity,ncx,ncy,l,nrel
      real xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      integer kount(50)
!     real ::  cofum
!     common/mudmd/cofum(1) 
      
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
c
c     prerelax at current finest grid level
c
      do l=1,iprer
	call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kcur .eq. 1) go to 5
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+9*ncx*ncy
c     call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))

      call resm2cm(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     |             wk(kps),cofum)
!    |             wk(kps),cofum(1))
c
c    set counter for grid levels to zero
c
      do l = 1,kcur
	kount(l) = 0
      end do
c
c    set new grid level and continue k-cycling
c
      klevel = kcur-1
      nrel = iprer
c
c   kcycle control point
c
   10 continue
c
c      post relax when kcur revisited
c
      if (klevel .eq. kcur) go to 5
c
c   count hit at current level
c
      kount(klevel) = kount(klevel)+1
c
c   relax at current level
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      do l=1,nrel
	call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      if (kount(klevel) .eq. kcycle+1) then
c
c     kcycle complete at klevel
c
	ipc = ip
	ip = kpbgn(klevel+1)
	ncx = nxk(klevel)
	ncy = nyk(klevel)
	nx = nxk(klevel+1)
	ny = nyk(klevel+1)
c
c    inject correction to finer grid
c
	call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
c
c    reset counter to zero
c
	kount(klevel) = 0
c
c     ascend to next higher level and set to postrelax there
c
	klevel = klevel+1
	nrel = ipost
	go to 10
      else
	if (klevel .gt. 1) then
c
c    kcycle not complete so descend unless at coarsest grid
c
	  ipc = kpbgn(klevel-1)
	  ncx = nxk(klevel-1)
	  ncy = nyk(klevel-1)
	  irc = kcbgn(klevel-1)+9*ncx*ncy
	  call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     +                wk(kps))
c
c     prerelax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 10
	else
c
c    postrelax at coarsest level
c
	  do l=1,ipost
	    call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
	  end do
	  ipc = ip
	  ip = kpbgn(2)
	  ncx = nxk(1)
	  ncy = nyk(1)
	  nx = nxk(2)
	  ny = nyk(2)
c
c    inject correction to level 2
c
	call cor2(nx,ny,wk(ip),ncx,ncy,wk(ipc),nxa,nxb,nyc,nyd,
     +            intpol,wk(kps))
c
c     set to postrelax at level 2
c
	  nrel = ipost
	  klevel = 2
	  go to 10
	end if
      end if
    5 continue
c
c     post relax at current finest grid level
c
      nx = nxk(kcur)
      ny = nyk(kcur)
      ip = kpbgn(kcur)
      ic = kcbgn(kcur)
      itx = ktxbgn(kcur)
      ity = ktybgn(kcur)
      do l=1,ipost
	call relmd2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end
c----------------------------------------------------------------------      
      subroutine resm2cm(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf,cofum)
c
c     restrict residual from fine to coarse mesh using fully weighted
c     residual restriction
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,ncx,ncy,i,j,ic,jc
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real rhsc(ncx,ncy),resf(nx,ny)
      real phi(0:nx+1,0:ny+1),phic(0:ncx+1,0:ncy+1)
      real cof(nx,ny,10),cofum(nx,ny,9)
c
c     set phic zero
c
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc) = 0.0
	end do
      end do
      
      call bnd2cm(nx,ny,cofum)
c
c     compute residual on fine mesh in resf
c
!$OMP PARALLEL DO SHARED(resf,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny
	do i=1,nx
!	  resf(i,j) = cof(i,j,10)-(
!     +                cof(i,j,1)*phi(i+1,j)+
!     +                cof(i,j,2)*phi(i+1,j+1)+
!     +                cof(i,j,3)*phi(i,j+1)+
!     +                cof(i,j,4)*phi(i-1,j+1)+
!     +                cof(i,j,5)*phi(i-1,j)+
!     +                cof(i,j,6)*phi(i-1,j-1)+
!     +                cof(i,j,7)*phi(i,j-1)+
!     +                cof(i,j,8)*phi(i+1,j-1)+
!     +                cof(i,j,9)*phi(i,j))
	  resf(i,j) = cof(i,j,10)-(
     +                cofum(i,j,1)*phi(i+1,j)+
     +                cofum(i,j,2)*phi(i+1,j+1)+
     +                cofum(i,j,3)*phi(i,j+1)+
     +                cofum(i,j,4)*phi(i-1,j+1)+
     +                cofum(i,j,5)*phi(i-1,j)+
     +                cofum(i,j,6)*phi(i-1,j-1)+
     +                cofum(i,j,7)*phi(i,j-1)+
     +                cofum(i,j,8)*phi(i+1,j-1)+
     +                cofum(i,j,9)*phi(i,j))
	end do

!       do i=1,9
!         write(6,"('resm2cm: j=',i3,' i=',i3,' cofum(:,j,i)=',/,
!    |      (6e12.4))") j,i,cofum(:,j,i)
!       enddo

      end do
c
c     restrict resf to coarse mesh in rhsc
c
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end

c-----------------------------------------------------------------------
      subroutine bnd2cm(nx,ny,cf)
c
c     set stencil & boundary condition for finest stencil
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,i,j,kbdy,l,im1,jm1,ier,jc,nnx,nny
      real cf(nx,ny,*)
      real dlx,dlx2,dlxx,dly,dly2,dlyy,cmin,alfmax,cemax
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
 
c
c     set coefficient for specified boundaries
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  do l=1,9
	    cf(i,j,l) = 0.0
	  end do
	  cf(i,j,9) = 1.0
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  do l=1,9
	    cf(i,j,l) = 0.0
	  end do
	  cf(i,j,9) = 1.0
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  do l=1,9
	    cf(i,j,l) = 0.0
	  end do
	  cf(i,j,9) = 1.0
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  do l=1,9
	    cf(i,j,l) = 0.0
	  end do
	  cf(i,j,9) = 1.0
	end do
      end if
c
      return
      end

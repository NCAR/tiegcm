!#include "dims.h"
c
c     file muh2cr.f
!     am 01/04/16 runs only as direct solver or with modified stencils
!                 it's very slow
! 5/02 B. Foster:
! Use-associate coefficients from dynamo module in tiegcm1 (dynamo.F).
c
c ... required mudpack files
c
c     mudcom.f
c
c
      subroutine muh(pe,jntl)
      use dynamo_module,only: nc,cee
      implicit none
      integer jntl
c
c     set grid size params
c
      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork,iiwork
c      parameter (iixp = 5 , jjyq = 3, iiex = 5, jjey = 5 )
c      parameter (iixp = 10, jjyq = 6,iiex = 4, jjey = 4 )
c      parameter (iixp = 20, jjyq = 12,iiex = 3, jjey = 3 )
c      parameter (iixp = 40, jjyq = 24,iiex = 2, jjey = 2 )
      parameter (iixp = 80 , jjyq = 48,iiex = 1, jjey = 1 )
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
c
c     estimate work space for point relaxation (see muh2cr.d)
c
      parameter (llwork=(5*((nnx+2)*(nny+2)+18*nnx*nny)/3+
     |            (nnx+2)*(nny+2)+ (iixp+1)*(jjyq+1)*(2*iixp+3)))
      parameter (iiwork=(iixp+1)*(jjyq+1))
      real phi(nnx,nny),rhs(nnx,nny),work(llwork)
      integer iwork(iiwork)
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
c      DATA MAXCYA/20/
      DATA MAXCYA/1/
      integer mm,nn,jj,jjj
      real pi
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
      mgopt(2) = 2
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
c
c     set specified boundaries in phi
c
      DO I=1,NX
        PHI(I,NY) = RHS(I,NY)/CEE(I+(NY-1)*NX+8*NX*NY)
      END DO
      
!      write(*,100)
  100 format(//' mud2cr test ')
!      write (*,101) (iprm(i),i=1,15)
  101 format(/,' integer input arguments ',/,
     |  ' intl =  ',i2,/,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,
     |  ' nyd =   ',i2,/,' ixp = ',i2,' jyq = ',i2,' iex = ',i2,
     |  ' jey =   ',i2,/,' nx =  ',i3,' ny =  ',i3,' iguess = ',i2,
     |  ' maxcy = ',i3,/,' method = ',i2, ' work space estimate = ',i7)
!      write (*,102) (mgopt(i),i=1,4)
  102 format(/' multigrid option arguments ',
     |  /,' kcycle = ',i2,
     |  /,' iprer = ',i2,
     |  /,' ipost = ',i2
     |  /,' intpol = ',i2)
!      write(*,103) xa,xb,yc,yd,tolmax
  103 format(/' floating point input parameters ',
     |  /,' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3,
     |  /,' tolerance (error control) =   ',e10.3)
!      write(6,"('fprm(1-5) (xa,xb,yc,yd,tolmax=',6f8.3)") fprm(1:5)
c
c     intialization call
c
!      write(*,104) intl
  104 format(/' discretization call to muh2cr', ' intl = ', i2)
      call muh2cr(iprm,fprm,work,iwork,rhs,phi,mgopt,ierror)
!      write (*,200) ierror,iprm(16)
  200 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     attempt solution
c
      intl = 1
!      write(*,106) intl,method,iguess
  106 format(/' approximation call to muh2cr',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      call muh2cr(iprm,fprm,work,iwork,rhs,phi,mgopt,ierror)
!      write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
C
C     COPY PHI TO PE
C
      DO J = 1,NY
	JJ = NY+J-1
	JJJ = NY+1-J
	DO I = 1,NX
	  PE(I,JJ)  = PHI(I,J)
	  PE(I,JJJ) = PHI(I,J)
	END DO
      END DO
      end
!-----------------------------------------------------------------------      
      subroutine muh2cr(iparm,fparm,wk,iwk,rhs,phi,mgopt,
     +                  ierror)
      implicit none
      integer iparm(17),mgopt(4),ierror,iwk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real fparm(6),xa,xb,yc,yd,tolmax,relmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer int,iw,k,kb,nx,ny,ic,itx,ity
      real wk(*),phi(*),rhs(*)
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      integer ibeta,ialfa,izmat,idmat
      common/mh2cr/ibeta,ialfa,izmat,idmat
      data int / 0 /
      save int
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
c
c     set pointers for direct at coarse grid
c
	nx = ixp+1
	ny = jyq+1
	ibeta = kps+1
	if (nyc .eq. 0) then
	  ialfa = ibeta + nx*nx*(ny-1)
	  izmat = ialfa+nx*nx*(ny-1)
	  idmat = izmat+nx*nx*(ny-2)
	  kps = idmat+nx*nx*(ny-2)
	else
	  ialfa = ibeta + nx*nx*ny
	  kps = ialfa+nx*nx*ny
	end if
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
	  call dismh2cr(nx,ny,wk(ic),wk(itx),wk(ity),
     +                  wk,iwk,ierror)
	  end do
	return
      end if   ! end of intl=0 initialization call block
      nx = nfx
      ny = nfy
      call muh2cr1(nx,ny,rhs,phi,wk,iwk)
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) ierror = -1   ! flag convergenc failure
      end if
      return
      end

      subroutine muh2cr1(nx,ny,rhsf,phif,wk,iwk)
      implicit none
      integer nx,ny,iwk(*)
      real phif(nx,ny),rhsf(nx,ny),wk(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax,phmax
      integer kpbgn,kcbgn,ktxbgn,ktybgn,nxk,nyk,isx,jsy
      integer k,kb,ip,ic,ir,ipc,irc,icc
      integer ncx,ncy,jj,ij,i,j,iter
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
      integer ibeta,ialfa,izmat,idmat
      common/mh2cr/ibeta,ialfa,izmat,idmat
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
	  call adjmh2cr(nx,ny,wk(ip),wk(ic))
	end do
c
c     execute one full multigrid cycle
c
	do k=1,ngrid-1
	  kcur = k
	  call kcymh2cr(wk,iwk)
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
	call adjmh2cr(nx,ny,wk(ip),wk(ic))
      end if
c
c     execute maxcy more multigrid k cycles from finest level
c
      kcur = ngrid
      do iter=1,maxcy
	itero = iter
	call kcymh2cr(wk,iwk)
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

      subroutine kcymh2cr(wk,iwk)
c
c     execute multigrid k cycle from kcur grid level
c     kcycle=1 for v cycles, kcycle=2 for w cycles
c
      implicit none
      integer iwk(*)
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
      integer ibeta,ialfa,izmat,idmat
      common/mh2cr/ibeta,ialfa,izmat,idmat
      integer kount(50)
      klevel = kcur
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      if (kcur .eq. 1) then
c
c     solve at coarse level with direct method and return
c
	if (nyc .ne. 0) then
	  call dir2cr(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),iwk,nxa)
	  return
	else
	  call dir2crp(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +               wk(izmat),wk(idmat),iwk,nxa)
	  return
	end if
      end if
c
c     prerelax at current finest grid level > 1
c
      do l=1,iprer
	call relmh2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
c
c     restrict residual to kcur-1 level
c
      ipc = kpbgn(klevel-1)
      ncx = nxk(klevel-1)
      ncy = nyk(klevel-1)
      irc = kcbgn(klevel-1)+9*ncx*ncy
      call resmh2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     |     wk(kps))
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
c   relax or solve directly at current level
c
      nx = nxk(klevel)
      ny = nyk(klevel)
      ip = kpbgn(klevel)
      ic = kcbgn(klevel)
      itx = ktxbgn(klevel)
      ity = ktybgn(klevel)
      if (klevel.gt.1) then
	do l=1,nrel
	  call relmh2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
	end do
      else
c
c     use direct method at coarsest level
c
	if (nyc .ne. 0) then
	  call dir2cr(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),iwk,nxa)
	else
	  call dir2crp(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +               wk(izmat),wk(idmat),iwk,nxa)
	end if
c
c     insure direct method is not called again at coarse level
c
	kount(1) = kcycle+1
      end if
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
	  call resmh2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),
     +                wk(kps))
c
c     prerelax at next coarser level
c
	  klevel = klevel-1
	  nrel = iprer
	  go to 10
	else
c
c    direct  at coarsest level takes place of postrelax
c
	  ip = kpbgn(1)
	  ic = kcbgn(1)
	  nx = nxk(1)
	  ny = nyk(1)
	  if (nyc .ne. 0) then
	    call dir2cr(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     |                  iwk,nxa)
	  else
	    call dir2crp(nx,ny,wk(ip),wk(ic),wk(ibeta),wk(ialfa),
     +                 wk(izmat),wk(idmat),iwk,nxa)
	  end if
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
	call relmh2cr(nx,ny,wk(ip),wk(ic),wk(itx),wk(ity),wk(kps))
      end do
      return
      end

      subroutine dismh2cr(nx,ny,cf,tx,ty,wk,iwk,ier)
      use dynamo_module,only: nc,cee,ceee
c
c     discretize elliptic pde for muh2cr, set nonfatal errors
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,iwk(*),i,j,kbdy,l,im1,jm1,ier,jc
      real cf(nx,ny,10),tx(nx,ny,*),ty(ny,nx,*)
      real wk(*)
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      integer ibeta,ialfa,izmat,idmat
      common/mh2cr/ibeta,ialfa,izmat,idmat
      integer nnx,nny
c
c     CHECK FOR CONSISTENCYT WRT KLEVEL
c
      NNX = ixp*2**(KLEVEL-1)+1
      NNY = jyq*2**(KLEVEL-1)+1
      IF(NNX.NE.NX.OR.NNY.NE.NY)THEN
	WRITE(6,100)NX,NY,NNX,NNY,ixp,jyq,KLEVEL
  100   FORMAT(' INCONSISTENCY WRT LEVEL. NX,NY,NNX,NNY,ixp,jyq,',
     |    'klevel = ',8I6)
c	STOP
      ENDIF
      CALL CEEE(CEE(NC(6-KLEVEL-4)),NX,NY,CF)  ! sub ceee is in mud.F
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
      if (klevel .eq. 1) then
c
c     set block tri-diagonal coefficient matrix and do lu decomposition
c     for direct method at coarsest grid level
c
	nx = ixp+1
	ny = jyq+1
	if (nyc .ne. 0) then
c     factor non-periodic block matrix
	  call lud2cr(nx,ny,cf,wk(ibeta),wk(ialfa),iwk,nxa)
	  return
	else
c     factor periodic block matrix

	  do j =1,ny-1
	    call setbcr(nx,ny,cf,wk(ibeta),j,nxa)
	    call setacr(nx,ny,cf,wk(ialfa),j,nxa)
	  end do
	  call lud2crp(nx,ny,cf,wk(ibeta),wk(ialfa),wk(izmat),
     +                 wk(idmat),iwk,nxa)
	  return
	end if
      end if
c
c     set and factor tridiagonal matrices for line relaxation(s) if flagged
c
      if (method.eq.1.or.method.eq.3) then
	if (nxa.ne.0) then
c
c    nonperiodic x line relaxation
c
	  do i=1,nx
	    im1 = max0(i-1,1)
	    do j=1,ny
	      tx(im1,j,1) = cf(i,j,5)
	      tx(i,j,2) = cf(i,j,9)
	      tx(i,j,3) = cf(i,j,1)
	    end do
	  end do
	  call factri(ny,nx,tx(1,1,1),tx(1,1,2),tx(1,1,3))
	else
c
c     periodic x line relaxation
c
	  if (nx .gt. 3) then
c
c     set and factor iff nx > 3
c
	    do i=1,nx-1
	      do j=1,ny
		tx(i,j,1) = cf(i,j,5)
		tx(i,j,2) = cf(i,j,9)
		tx(i,j,3) = cf(i,j,1)
	      end do
	    end do
	    call factrp(ny,nx,tx,tx(1,1,2),tx(1,1,3),tx(1,1,4),
     +                  tx(1,1,5),wk(kps))
	  end if
	end if
      end if

      if (method.eq.2.or.method.eq.3) then
	if (nyc.ne.0) then
c
c     nonperiodic y line relaxation
c
	  do j=1,ny
	    jm1 = max0(j-1,1)
	    do i=1,nx
	      ty(jm1,i,1) = cf(i,j,7)
	      ty(j,i,2) = cf(i,j,9)
	      ty(j,i,3) = cf(i,j,3)
	    end do
	  end do
	  call factri(nx,ny,ty(1,1,1),ty(1,1,2),ty(1,1,3))
	else
c
c      periodic y line relaxation
c
	  if (ny .gt. 3) then
c
c     set and factor iff ny > 3
c
	    do j=1,ny-1
	      do i=1,nx
		ty(j,i,1) = cf(i,j,7)
		ty(j,i,2) = cf(i,j,9)
		ty(j,i,3) = cf(i,j,3)
	      end do
	    end do
	    call factrp(nx,ny,ty,ty(1,1,2),ty(1,1,3),ty(1,1,4),
     +                  ty(1,1,5),wk(kps))
	  end if
	end if
      end if
      return
      end

      subroutine lud2cr(nx,ny,cof,beta,alfa,index,nxa)
c
c     decompose nonperiodic block coefficient matrix
c
      implicit none
      integer nx,ny,nxa,index(nx,ny)
      real cof(nx,ny,10),beta(nx,nx,*),alfa(nx,nx,*)
      integer iz,i1,jcur,jm1,l,lm1,lp1,k,i
      real gama,sum
      iz = 0
      i1 = 1
c
c     set and factor umat(1) in beta(1)
c
      jcur = 1
      call setbcr(nx,ny,cof,beta,jcur,nxa)
      call sgfa(beta,nx,nx,index,iz)

      do jcur=2,ny
c
c     solve transpose of lmat(jcur)*beta(jcur-1) = alfa(jcur) in alfa(jcur)
c
	call setacr(nx,ny,cof,alfa,jcur,nxa)
	call transp(nx,alfa(1,1,jcur))
	jm1 = jcur-1
	do l=1,nx
	  call sgsl(beta(1,1,jm1),nx,nx,index(1,jm1),alfa(1,l,jcur),i1)
	end do
	call transp(nx,alfa(1,1,jcur))
	call setbcr(nx,ny,cof,beta,jcur,nxa)
	do i=1,nx
	  do l=1,nx
	    sum = 0.0
	    lm1=max0(1,l-1)
	    lp1=min0(l+1,nx)
	    do k=lm1,lp1
	      if (k .eq. l+1) then
		gama = cof(k,jcur-1,4)
	      else if (k.eq. l) then
		gama = cof(k,jcur-1,3)
	      else if (k .eq. l-1) then
		gama = cof(k,jcur-1,2)
	      else
		gama=0.0
	      end if
	      sum = sum+alfa(i,k,jcur)*gama
	    end do
	    if (nxa.eq.0) then
	      if (l .eq. 2) then
		sum=sum+alfa(i,nx,jcur)*cof(nx,jcur-1,2)
	      end if
	      if (l .eq. nx-1) then
		sum=sum+alfa(i,1,jcur)*cof(1,jcur-1,4)
	      end if
	    end if
	    beta(i,l,jcur) = beta(i,l,jcur)-sum
	  end do
	end do
c
c     factor current beta for next pass
c
	iz = 0
	call sgfa(beta(1,1,jcur),nx,nx,index(1,jcur),iz)
      end do
      return
      end

      subroutine dir2cr(nx,ny,phi,cof,beta,alfa,index,nxa)
c
c     direct solve at coarsest grid
c
      implicit none
      integer nx,ny,index(nx,ny),nxa
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      real beta(nx,nx,*),alfa(nx,nx,*)
c     forward sweep
      call for2cr(nx,ny,phi,cof(1,1,10),alfa)
c     backward sweep
      call bkw2cr(nx,ny,phi,cof,beta,index,nxa)
      return
      end

      subroutine for2cr(nx,ny,phi,frhs,alfa)
c
c     forward sweep
c
      implicit none
      integer nx,ny,i,j,l
      real phi(0:nx+1,0:ny+1),frhs(nx,ny),alfa(nx,nx,*),sum
      do j=1,ny
	do i=1,nx
	  phi(i,j)=frhs(i,j)
	end do
      end do
      do j=2,ny
	do i=1,nx
	  sum=0.0
	  do l=1,nx
	    sum=sum+alfa(i,l,j)*phi(l,j-1)
	  end do
	  phi(i,j)=phi(i,j)-sum
	end do
      end do
      return                                                                    
      end                                                                       

      subroutine bkw2cr(nx,ny,phi,cof,beta,index,nxa)
      implicit none
      integer nx,ny,index(nx,ny),nxa
      real beta(nx,nx,*),sum
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      integer iz,jcur,jb,j,i
      iz = 0
      jcur=ny
      call sgsl(beta(1,1,jcur),nx ,nx ,index(1,jcur),phi(1,jcur),iz)
      do jb=2,ny
	j=ny-jb+1
	jcur=j
	do i=2,nx-1
	  sum=cof(i,j,2)*phi(i+1,j+1)+cof(i,j,3)*phi(i,j+1)+cof(i,j,4)*
     +    phi(i-1,j+1)
	  phi(i,j)=phi(i,j)-sum
	end do
	phi(1,j)=phi(1,j)-(cof(1,j,2)*phi(2,j+1)+cof(1,j,3)*phi(1,j+1))
	phi(nx,j)=phi(nx,j)-(cof(nx,j,3)*phi(nx,j+1)+cof(nx,j,4)*
     +  phi(nx-1,j+1))
	if (nxa .eq.0) then
	  phi(1,j)=phi(1,j)-cof(1,j,4)*phi(nx-1,j+1)
	  phi(nx,j)=phi(nx,j)-cof(nx,j,2)*phi(2,j+1)
	end if
	call sgsl(beta(1,1,jcur),nx ,nx ,index(1,jcur),phi(1,jcur),iz)
      end do
      return                                                                    
      end                                                                       

      subroutine lud2crp(nx,ny,cof,beta,alfa,zmat,dmat,index,nxa)
c
c     decompose periodic block tridiagonal matrix for direct at coarsest grid
c
      implicit none
      integer nx,ny,index(nx,ny),nxa
      real cof(nx,ny,10),alfa(nx,nx,*),beta(nx,nx,*)
      real dmat(nx,nx,*),zmat(nx,nx,*),sum,gama
      integer iz,j,jcur,i,l,jm1,i1,lm1,lp1,k
      jcur = 1
c
c     set dmat(1)=alfa(1)
c
      call setacr(nx,ny,cof,alfa,jcur,nxa)
      do i=1,nx
	do l=1,nx
	  dmat(i,l,1) = alfa(i,l,1)
	end do
      end do
      iz = 0
c
c     factor umat(1) in beta(1)
c
      call setbcr(nx,ny,cof,beta,jcur,nxa)
      call sgfa(beta(1,1,1),nx,nx,index(1,1),iz)
      do jcur=2,ny-2
c
c     solve transpose of lmat(jcur)umat(jcur-1)=alfa(jcur) in alfa(jcur)
c
	call setacr(nx,ny,cof,alfa,jcur,nxa)
	call transp(nx,alfa(1,1,jcur))
	jm1 = jcur-1
	i1 = 1
	do l=1,nx
	 call sgsl(beta(1,1,jm1),nx,nx,index(1,jm1),alfa(1,l,jcur),i1)
	end do
	call transp(nx,alfa(1,1,jcur))
	call setbcr(nx,ny,cof,beta,jcur,nxa)
	do i=1,nx
	  do l=1,nx
	    sum = 0.0
	    lm1=max0(1,l-1)
	    lp1=min0(l+1,nx)
	    do k=lm1,lp1
	      if (k .eq. l+1) then
		gama = cof(k,jcur-1,4)
	      else if (k.eq. l) then
		gama = cof(k,jcur-1,3)
	      else if (k .eq. l-1) then
		gama = cof(k,jcur-1,2)
	      else
		gama=0.0
	      end if
	      sum = sum+alfa(i,k,jcur)*gama
	    end do
	    if (nxa.eq.0) then
	      if (l .eq. 2) then
		sum=sum+alfa(i,nx,jcur)*cof(nx,jcur-1,2)
	      end if
	      if (l .eq. nx-1) then
		sum=sum+alfa(i,1,jcur)*cof(1,jcur-1,4)
	      end if
	    end if
	    beta(i,l,jcur)=beta(i,l,jcur)-sum
	  end do
	end do
c
c     factor current beta(1,1,jcur) for next pass
c
	call sgfa(beta(1,1,jcur),nx ,nx,index(1,jcur),iz)
c
c     set dmat(jcur) = -alfa(jcur)*dmat(jcur-1)
c
	do i=1,nx
	  do j=1,nx
	    dmat(i,j,jcur) = 0.0
	    do l=1,nx
	      dmat(i,j,jcur) = dmat(i,j,jcur)-alfa(i,l,jcur)*
     +                         dmat(l,j,jcur-1)
	    end do
	  end do
	end do
	if (jcur .eq. ny-2) then
c
c     adjust dmat(ny-2) = gama(ny-2)-alfa(ny-2)*dmat(ny-3)
c
	  dmat(1,1,jcur) = cof(1,jcur,3) + dmat(1,1,jcur)
	  dmat(1,2,jcur) = cof(1,jcur,2) + dmat(1,2,jcur)
c
c     adjust for periodic b.c. in x
c
	  if (nxa .eq. 0) then
	    dmat(1,nx-1,jcur) = cof(1,jcur,4) + dmat(1,nx-1,jcur)
	    dmat(nx,2,jcur) = cof(nx,jcur,2) + dmat(nx,2,jcur)
	  end if
c
c     matrix interior
c
	  do i=2,nx-1
	    dmat(i,i,jcur) = cof(i,jcur,3) + dmat(i,i,jcur)
	    dmat(i,i-1,jcur) = cof(i,jcur,4) + dmat(i,i-1,jcur)
	    dmat(i,i+1,jcur) = cof(i,jcur,2) + dmat(i,i+1,jcur)
	  end do
	  dmat(nx,nx,jcur) = cof(nx,jcur,3) + dmat(nx,nx,jcur)
	  dmat(nx,nx-1,jcur) = cof(nx,jcur,4) + dmat(nx,nx-1,jcur)
	end if
      end do
c
c     final phase with periodic factorization
c
c     solve transpose of zmat(1) beta(1) = gama(ny-1)
c
      zmat(1,1,1) = cof(1,ny-1,3)
      zmat(1,2,1) = cof(1,ny-1,2)
      do l=3,nx
	zmat(1,l,1) = 0.0
      end do

      do i=2,nx-1
	do l=1,nx
	  zmat(i,l,1) = 0.0
	end do
	zmat(i,i,1) = cof(i,ny-1,3)
	zmat(i,i+1,1) = cof(i,ny-1,2)
	zmat(i,i-1,1) = cof(i,ny-1,4)
      end do
      zmat(nx,nx-1,1) = cof(nx,ny-1,4)
      zmat(nx,nx,1) = cof(nx,ny-1,3)
      do l=1,nx-2
	zmat(nx,l,1) = 0.0
      end do
c
c     adjust for periodic x b.c.
c
      if (nxa .eq.0) then
	zmat(1,nx-1,1) = cof(1,ny-1,4)
	zmat(nx,2,1) = cof(nx,ny-1,2)
      end if
      call transp(nx,zmat(1,1,1))
      do l=1,nx
	call sgsl(beta(1,1,1),nx,nx,index(1,1),zmat(1,l,1),i1)
      end do
      call transp(nx,zmat(1,1,1))
      do jcur = 2,ny-3
c
c     solve transpose of zmat(jcur) umat(jcur) = -zmat(jcur-1) gama(jcur-1)
c
	do i=1,nx
	  zmat(i,1,jcur) = -(zmat(i,1,jcur-1)*cof(1,jcur-1,3) +
     +                         zmat(i,2,jcur-1)*cof(2,jcur-1,4))
	end do
	do i=1,nx
	  do l=2,nx-1
	    zmat(i,l,jcur) = -(zmat(i,l-1,jcur-1)*cof(l-1,jcur-1,2) +
     +                         zmat(i,l,jcur-1)*cof(l,jcur-1,3) +
     +                     zmat(i,l+1,jcur-1)*cof(l+1,jcur-1,4))
	  end do
	end do
	do i=1,nx
	  zmat(i,nx,jcur) = -(zmat(i,nx-1,jcur-1)*cof(nx-1,jcur-1,2) +
     +                       zmat(i,nx,jcur-1)*cof(nx,jcur-1,3))
	end do
c
c     adjust j=2 and j=nx-1 column if periodic in x
c
	if (nxa .eq. 0) then
	  do i=1,nx
	    zmat(i,2,jcur)=zmat(i,2,jcur)-zmat(i,nx,jcur-1)*
     +                     cof(nx,jcur-1,2)
	    zmat(i,nx-1,jcur)=zmat(i,nx-1,jcur)-zmat(i,1,jcur-1)*
     +                        cof(1,jcur-1,4)
	  end do
	end if
	call transp(nx,zmat(1,1,jcur))
	do l=1,nx
	  call sgsl(beta(1,1,jcur),nx,nx,index(1,jcur),zmat(1,l,jcur),
     |              i1)
	end do
	call transp(nx,zmat(1,1,jcur))
      end do
c
c     solve transpose of zmat(ny-2)umat(ny-2)=alfa(ny-1)-zmat(ny-3)gama(ny-3)
c
      jcur = ny-2
      do i=1,nx
	zmat(i,1,jcur) = -(zmat(i,1,jcur-1)*cof(1,jcur-1,3) +
     +                     zmat(i,2,jcur-1)*cof(2,jcur-1,4))
      end do

      do i=1,nx
	do l=2,nx-1
	  zmat(i,l,jcur) = -(zmat(i,l-1,jcur-1)*cof(l-1,jcur-1,2) +
     +                   zmat(i,l,jcur-1)*cof(l,jcur-1,3) +
     +                   zmat(i,l+1,jcur-1)*cof(l+1,jcur-1,4))
	end do
      end do
      do i=1,nx
	zmat(i,nx,jcur) = -(zmat(i,nx-1,jcur-1)*cof(nx-1,jcur-1,2) +
     +                   zmat(i,nx,jcur-1)*cof(nx,jcur-1,3))
      end do
c
c     adjust j=2 and j=nx-1 column if periodic in x
c
      if (nxa .eq. 0) then
	do i=1,nx
	  zmat(i,2,jcur)=zmat(i,2,jcur)-zmat(i,nx,jcur-1)*
     |                       cof(nx,jcur-1,2)
	  zmat(i,nx-1,jcur)=zmat(i,nx-1,jcur)-zmat(i,1,jcur-1)*
     +                      cof(1,jcur-1,4)
	end do
      end if
      call setacr(nx,ny,cof,alfa,ny-1,nxa)
      do i=1,nx
	do l=1,nx
	  zmat(i,l,ny-2) = alfa(i,l,ny-1) + zmat(i,l,ny-2)
	end do
      end do
      call transp(nx,zmat(1,1,ny-2))
      do l=1,nx
	call sgsl(beta(1,1,ny-2),nx,nx,index(1,ny-2),zmat(1,l,ny-2),i1)
      end do
      call transp(nx,zmat(1,1,ny-2))
c
c     set umat(ny-1) = beta(ny-1)-(zmat(1)*dmat(1)+...+zmat(ny-2)*dmat(ny-2))
c     in beta(ny-1)
c
      call setbcr(nx,ny,cof,beta,ny-1,nxa)
      do i=1,nx
	do j=1,nx
	  sum = 0.0
	  do jcur=1,ny-2
	    do l=1,nx
	      sum = sum + zmat(i,l,jcur)*dmat(l,j,jcur)
	    end do
	  end do
	  beta(i,j,ny-1) = beta(i,j,ny-1) - sum
	end do
      end do
c
c     factor bmat(ny-1) for backward sweep
c
      call sgfa(beta(1,1,ny-1),nx,nx,index(1,ny-1),iz)
c
c     lud is now complete
c
      return
      end

      subroutine dir2crp(nx,ny,phi,cof,beta,alfa,zmat,dmat,index,nxa)
      implicit none
      integer nx,ny,index(nx,ny),nxa
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      real beta(nx,nx,*),alfa(nx,nx,*)
      real zmat(nx,nx,*), dmat(nx,nx,*)
c     forward sweep
      call for2crp(nx,ny,phi,cof(1,1,10),alfa,zmat)
c     backward sweep
      call bkw2crp(nx,ny,phi,cof,beta,dmat,index,nxa)
      return
      end

      subroutine for2crp(nx,ny,phi,frhs,alfa,zmat)
      implicit none
      integer nx,ny,i,j,l,jcur,k
      real frhs(nx,ny)
      real phi(0:nx+1,0:ny+1)

      real alfa(nx,nx,*),zmat(nx,nx,*)
      real sum
      do j=1,ny-1
	do i=1,nx
	  phi(i,j)=frhs(i,j)
	end do
      end do
      do jcur=2,ny-2
	do i=1,nx
	  sum=0.0
	  do l=1,nx
	    sum=sum+alfa(i,l,jcur)*phi(l,jcur-1)
	  end do
	  phi(i,jcur)=phi(i,jcur)-sum
	end do
      end do
c
c     solve:
c     zmat(1)*phi(1)+...+zmat(ny-2)*phi(ny-2) + phi(ny-1) = f(ny-1)
c
      do i=1,nx
	sum = 0.0
	do k=1,ny-2
	  do l=1,nx
	    sum = sum + zmat(i,l,k)*phi(l,k)
	  end do
	end do
	phi(i,ny-1) = phi(i,ny-1) - sum
      end do
      return
      end

      subroutine bkw2crp(nx,ny,phi,cof,beta,dmat,index,nxa)
      implicit none
      integer nx,ny,index(nx,ny),nxa
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      real beta(nx,nx,ny),dmat(nx,nx,*)
      integer iz,i,l,kb,k
       real sum
      iz = 0
      call sgsl(beta(1,1,ny-1),nx,nx,index(1,ny-1),phi(1,ny-1),iz)
c
c     solve beta(ny-2)*phi(ny-2) = phi(ny-2)-dmat(ny-2)*phi(ny-1)
c
      do i=1,nx
	sum = 0.0
	do l=1,nx
	  sum = sum + dmat(i,l,ny-2)*phi(l,ny-1)
	end do
	phi(i,ny-2) = phi(i,ny-2) - sum
      end do
      call sgsl(beta(1,1,ny-2),nx,nx,index(1,ny-2),phi(1,ny-2),iz)
c
c     solve beta(k)*phi(k) = phi(k) - gama(k)*phi(k+1)-dmat(k)*phi(ny-1)
c     k=ny-3,...,1
c
      do kb=4,ny
	k = ny-kb+1
	sum = 0.0
	do l=1,nx
	  sum = sum+dmat(1,l,k)*phi(l,ny-1)
	end do
	phi(1,k) = phi(1,k)-sum - (  cof(1,k,3)*phi(1,k+1) +
     +                               cof(1,k,2)*phi(2,k+1))
	do i=2,nx-1
	  sum = 0.0
	  do  l=1,nx
	    sum = sum+dmat(i,l,k)*phi(l,ny-1)
	  end do
	  phi(i,k) = phi(i,k) - sum - (cof(i,k,4)*phi(i-1,k+1) +
     +                               cof(i,k,3)*phi(i,k+1)   +
     +                               cof(i,k,2)*phi(i+1,k+1))
	end do
	sum = 0.0
	do l=1,nx
	  sum = sum+dmat(nx,l,k)*phi(l,ny-1)
	end do
	phi(nx,k) = phi(nx,k) - sum - (cof(nx,k,4)*phi(nx-1,k+1) +
     +                                 cof(nx,k,3)*phi(nx,k+1))
c
c     adjust for periodic x b.c.
c
	if (nxa .eq. 0) then
	  phi(1,k) = phi(1,k) - cof(1,k,4)*phi(nx-1,k+1)
	  phi(nx,k) = phi(nx,k) - cof(nx,k,2)*phi(2,k+1)
	end if
	call sgsl(beta(1,1,k),nx,nx,index(1,k),phi(1,k),iz)
      end do
c
c     set j=ny by periodicity
c
      do i=1,nx
	phi(i,ny) = phi(i,1)
      end do
      return
      end

      subroutine setbcr(nx,ny,cof,beta,jcur,nxa)
c
c     set diagonal matrix on block
c
      implicit none
      integer nx,ny,jcur,nxa,i,l
      real cof(nx,ny,10),beta(nx,nx,*)
      do i=1,nx
	do l=1,nx
	  beta(i,l,jcur)=0.0
	end do
      end do
      do i=1,nx
	beta(i,i,jcur) = cof(i,jcur,9)
      end do
      do i=2,nx
	beta(i,i-1,jcur) = cof(i,jcur,5)
      end do
      do i=1,nx-1
	beta(i,i+1,jcur) = cof(i,jcur,1)
      end do
      if (nxa.eq.0) then                                                        
	beta(1,nx-1,jcur) = cof(1,jcur,5)
	beta(nx,2,jcur) = cof(nx,jcur,1)
      end if                                                                    
      return                                                                    
      end                                                                       

      subroutine setacr(nx,ny,cof,alfa,jcur,nxa)
      implicit none
      integer nx,ny,jcur,nxa,i,j
      real cof(nx,ny,10),alfa(nx,nx,*)
      do i=1,nx
	do j=1,nx
	  alfa(i,j,jcur)=0.0
	end do
      end do
      do i=2,nx
	alfa(i,i-1,jcur)=cof(i,jcur,6)
      end do
      do i=1,nx
	alfa(i,i,jcur)=cof(i,jcur,7)
      end do
      do i=1,nx-1
	alfa(i,i+1,jcur)=cof(i,jcur,8)
      end do
      if (nxa .eq. 0) then
c     adjust for x periodicity
	alfa(1,nx-1,jcur)=cof(1,jcur,6)
	alfa(nx,2,jcur)=cof(nx,jcur,8)
      end if
      return
      end

      subroutine adjmh2cr(nx,ny,phi,cf)
c
c     adjust righthand side in cf(i,j,10) for boundary conditions
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,i,j,kbdy
      real cf(nx,ny,10),phi(0:nx+1,0:ny+1)
      real dlx,dlx2,dlxx,dly,dly2,dlyy,dlxy,dlxy2,dlxy4,dxoy,dyox
      real x,y,cxx,cxy,cyy,cx,cy,ce,c1,c2,c3,c4,c5
      real c6,c7,c8
      real alfaa,alfab,alfac,alfad,betaa,betab,betac,betad,det
      real gamaa,gamab,gamac,gamad
      real alfim1,alfi,alfip1,betim1,beti,betip1,gamim1,gami,gamip1
      real alfjm1,alfj,alfjp1,betjm1,betj,betjp1,gamjm1,gamj,gamjp1
      real gbdim1,gbdi,gbdip1,gbdj,gbdjm1,gbdjp1
      real gbdya,gbdyb,gbdyc,gbdyd
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      

c
c     set specified boundaries in rhs from phi
c
      if (nxa.eq.1) then
	i = 1
	do j=1,ny
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      if (nxb.eq.1) then
	i = nx
	do j=1,ny
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      if (nyc.eq.1) then
	j = 1
	do i=1,nx
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      if (nyd.eq.1) then
	j = ny
	do i=1,nx
	  cf(i,j,10) = phi(i,j)
	end do
      end if
      return
      end

      subroutine resmh2cr(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
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
      real cof(nx,ny,10)
c
c     set phic zero
c
      do jc=0,ncy+1
	do ic=0,ncx+1
	  phic(ic,jc) = 0.0
	end do
      end do
c
c     compute residual on fine mesh in resf
c
CC!$OMP PARALLEL DO SHARED(resf,cof,phi,nx,ny) PRIVATE(i,j)
      do j=1,ny
	do i=1,nx
	  resf(i,j) = cof(i,j,10)-(
     +                cof(i,j,1)*phi(i+1,j)+
     +                cof(i,j,2)*phi(i+1,j+1)+
     +                cof(i,j,3)*phi(i,j+1)+
     +                cof(i,j,4)*phi(i-1,j+1)+
     +                cof(i,j,5)*phi(i-1,j)+
     +                cof(i,j,6)*phi(i-1,j-1)+
     +                cof(i,j,7)*phi(i,j-1)+
     +                cof(i,j,8)*phi(i+1,j-1)+
     +                cof(i,j,9)*phi(i,j))
	end do
      end do
c
c     restrict resf to coarse mesh in rhsc
c
      call res2(nx,ny,resf,ncx,ncy,rhsc,nxa,nxb,nyc,nyd)
      return
      end

      subroutine relmh2cr(nx,ny,phi,cof,tx,ty,sum)
c
c     relaxation for muh2cr
c
      implicit none
      integer nx,ny
      real phi(*),cof(*),tx(*),ty(*),sum(*)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      if (method.eq.0) then                ! point relaxation
	call relmh2crp(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
	call slxmh2cr(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slymh2cr(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxmh2cr(nx,ny,phi,cof,tx,sum)
	call slymh2cr(nx,ny,phi,cof,ty,sum)
      end if
      return
      end

      subroutine relmh2crp(nx,ny,phi,cof)
c
c     gauss-seidel four color point relaxation
c
      implicit none
      integer nx,ny,i,j,lcolor,i1,i2,i3,i4,it
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10)
      i1 = 1
      i2 = 4
      i3 = 3
      i4 = 2
c
c     sweep four colored grid points
c
      do lcolor=1,4
CC!$OMP PARALLEL DO SHARED(i1,cof,phi,nx,ny) PRIVATE(i,j)
	do j=1,ny,4
	  do i=i1,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
CC!$OMP PARALLEL DO SHARED(i2,cof,phi,nx,ny) PRIVATE(i,j)
	do j=2,ny,4
	  do i=i2,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
CC!$OMP PARALLEL DO SHARED(i3,cof,phi,nx,ny) PRIVATE(i,j)
	do j=3,ny,4
	  do i=i3,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
CC!$OMP PARALLEL DO SHARED(i4,cof,phi,nx,ny) PRIVATE(i,j)
	do j=4,ny,4
	  do i=i4,nx,4
	      phi(i,j) = (cof(i,j,10) - (
     +                    cof(i,j,1)*phi(i+1,j)   +
     +                    cof(i,j,2)*phi(i+1,j+1) +
     +                    cof(i,j,3)*phi(i,j+1)   +
     +                    cof(i,j,4)*phi(i-1,j+1) +
     +                    cof(i,j,5)*phi(i-1,j)   +
     +                    cof(i,j,6)*phi(i-1,j-1) +
     +                    cof(i,j,7)*phi(i,j-1)   +
     +                    cof(i,j,8)*phi(i+1,j-1)))/cof(i,j,9)
	  end do
	end do
c
c     set periodic virtual boundaries as necessary
c
	if (nxa.eq.0) then
	  do j=1,ny
	    phi(0,j) = phi(nx-1,j)
	    phi(nx+1,j) = phi(2,j)
	  end do
	end if
	if (nyc.eq.0) then
	  do i=1,nx
	    phi(i,0) = phi(i,ny-1)
	    phi(i,ny+1) = phi(i,2)
	  end do
	end if
c
c    permute (i1,i2,i3,i4) for next color
c
	it = i4
	i4 = i3
	i3 = i2
	i2 = i1
	i1 = it
      end do
      return
      end

      subroutine slxmh2cr(nx,ny,phi,cof,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,ib,j
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10),tx(nx,ny,*),sum(ny)
c
c     set periodic y virtual boundary if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if

      if (nxa.ne.0) then
CC!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
c
c     x direction not periodic, sweep odd j lines
c
	do j=1,ny,2
	  do i=1,nx
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
c
c     forward sweep
c
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
c
c     backward sweep
c
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	  do ib=2,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	  end do
	end do
c
c     sweep even j lines forward and back
c
CC!$OMP PARALLEL DO SHARED(cof,phi,tx,nx,ny) PRIVATE(i,ib,j)
	do j=2,ny,2
	  do i=1,nx
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
	  do i=2,nx
	    phi(i,j) = phi(i,j)-tx(i-1,j,1)*phi(i-1,j)
	  end do
	  phi(nx,j) = phi(nx,j)/tx(nx,j,2)
	  do ib=2,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j))/tx(i,j,2)
	  end do
	end do
      else
c
c     x direction periodic
c
	do j=1,ny
	  sum(j) = 0.0
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
c
c      sweep odd lines forward and back
c
CC!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=1,ny,2
	  do i=1,nx-1
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
c
c     forward sweep
c
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic and virtual points for j odd
c
	do j=1,ny,2
	  phi(nx,j) = phi(1,j)
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
c
c     sweep even j lines
c
CC!$OMP PARALLEL DO SHARED(sum,cof,phi,tx,nx,ny) PRIVATE(i,j,ib)
	do j=2,ny,2
	  do i=1,nx-1
	    phi(i,j) = cof(i,j,10)-(cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,3)*phi(i,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,7)*phi(i,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
c
c     forward sweep
c
	  do i=2,nx-2
	    phi(i,j) = phi(i,j)-tx(i,j,1)*phi(i-1,j)
	  end do
	  do i=1,nx-2
	    sum(j) = sum(j)+tx(i,j,5)*phi(i,j)
	  end do
	  phi(nx-1,j) = phi(nx-1,j)-sum(j)
c
c     backward sweep
c
	  phi(nx-1,j) = phi(nx-1,j)/tx(nx-1,j,2)
	  phi(nx-2,j) = (phi(nx-2,j)-tx(nx-2,j,4)*phi(nx-1,j))/
     +                   tx(nx-2,j,2)
	  do ib=4,nx
	    i = nx-ib+1
	    phi(i,j) = (phi(i,j)-tx(i,j,3)*phi(i+1,j)-tx(i,j,4)*
     +                 phi(nx-1,j))/tx(i,j,2)
	  end do
	end do
c
c     set periodic and virtual points for j even
c
	do j=2,ny,2
	  phi(nx,j) = phi(1,j)
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
c
c     set periodic y virtual boundaries if necessary
c
      if (nyc.eq.0) then
	do i=1,nx
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
      return
      end

      subroutine slymh2cr(nx,ny,phi,cof,ty,sum)
c
c     line relaxation in the y direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,j,jb
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10),ty(ny,nx,*),sum(nx)
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx,j) = phi(1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if

      if (nyc.ne.0) then
c
c     y direction not periodic
c
CC!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=1,nx,2
	  do j=1,ny
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
c
c     forward sweep thru odd x lines
c
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
c
c     forward sweep even x lines
c
CC!$OMP PARALLEL DO SHARED(cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
	  do j=2,ny
	    phi(i,j) = phi(i,j)-ty(j-1,i,1)*phi(i,j-1)
	  end do
c
c      backward sweep
c
	  phi(i,ny) = phi(i,ny)/ty(ny,i,2)
	  do jb=2,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1))/ty(j,i,2)
	  end do
	end do
      else
c
c     y direction periodic
c
	do i=1,nx
	  sum(i) = 0.0
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
c
c     forward sweep odd x lines
c
CC!$OMP PARALLEL DO SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=1,nx,2
	  do j=1,ny-1
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
c
c       set odd periodic and virtual y boundaries
c
	do i=1,nx,2
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
c
c     forward sweep even x lines
c
CC!$OMP PARALLEL DO SHARED(sum,cof,phi,ty,nx,ny) PRIVATE(i,j,jb)
	do i=2,nx,2
	  do j=1,ny-1
	    phi(i,j) = cof(i,j,10)-(cof(i,j,1)*phi(i+1,j)+
     +                              cof(i,j,2)*phi(i+1,j+1)+
     +                              cof(i,j,4)*phi(i-1,j+1)+
     +                              cof(i,j,5)*phi(i-1,j)+
     +                              cof(i,j,6)*phi(i-1,j-1)+
     +                              cof(i,j,8)*phi(i+1,j-1))
	  end do
	  do j=2,ny-2
	    phi(i,j) = phi(i,j)-ty(j,i,1)*phi(i,j-1)
	  end do
	  do j=1,ny-2
	    sum(i) = sum(i)+ty(j,i,5)*phi(i,j)
	  end do
	  phi(i,ny-1) = phi(i,ny-1)-sum(i)
c
c     backward sweep
c
	  phi(i,ny-1) = phi(i,ny-1)/ty(ny-1,i,2)
	  phi(i,ny-2) = (phi(i,ny-2)-ty(ny-2,i,4)*phi(i,ny-1))/
     +                   ty(ny-2,i,2)
	  do jb=4,ny
	    j = ny-jb+1
	    phi(i,j) = (phi(i,j)-ty(j,i,3)*phi(i,j+1)-ty(j,i,4)*
     +                  phi(i,ny-1))/ty(j,i,2)
	  end do
	end do
c
c       set even periodic and virtual y boundaries
c
	do i=2,nx,2
	  phi(i,0) = phi(i,ny-1)
	  phi(i,ny) = phi(i,1)
	  phi(i,ny+1) = phi(i,2)
	end do
      end if
c
c      set periodic and virtual x boundaries if necessary
c
      if (nxa.eq.0) then
	do j=1,ny
	  phi(0,j) = phi(nx-1,j)
	  phi(nx+1,j) = phi(2,j)
	end do
      end if
      return
      end

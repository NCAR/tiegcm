!#include "dims.h"
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
! 2/99 B. Foster:
! Introducing new mudpack v4.0.
!
c
c *** 2/99 John Adams:
c *** modification of mud(pe,jntl)for Ridley's PDE to use new mud2cr
c *** removal of coef, bndyc to input pde (they are redundant).
c     Coefficients cee(ncee) are passed thru ceee_com common block in 
c     ceee.h  subroutine dismd2cr and adjmd2cr are modified
c     all references to coef and bndyc are removed
c     the common subroutine file "mudcom.f" is attached to the end
!
! 5/02 B. Foster:
! Use-associate coefficients from dynamo module in tiegcm1 (dynamo.F).
!
! BOP
! !IROUTINE: mud
! !INTERFACE:
      subroutine mud(pe,jntl,isolve,ier)
!
! !USES:
      use pdynamo_module,only: nc,ncee,cee
      use timing_module,only: timer,timing
      use init_module,only: istep
!
! !DESCRIPTION: multigrid solver with modified stencil
!     input 
!       int jntl: error check first (0/1)
!       int isolve: type of solver (here 0)
!     output
!        pe: electric potential
!        ier: not converged ier <0
!
      implicit none
! !ARGUMENTS:
      integer jntl         ! error check first
      integer,intent(in) :: isolve  ! type of solver here isolve =0
! !RETURN VALUE:
      integer  ier             ! output: not converged ier < 0
!        
! !REVISION HISTORY:
! 21.02.05  <Astrid Maute> <include header> 
! 
! EOP
c
c     set grid size params
c
      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork
      parameter (iixp = 5 , jjyq = 3, iiex = 5, jjey = 5 )
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
c
c     estimate work space for point relaxation (see mud2cr.d)
c
      parameter (llwork=(7*(nnx+2)*(nny+2)+76*nnx*nny)/3 )
      real phi(nnx,nny),rhs(nnx,nny),work(llwork)
c
c     put integer and floating point argument names in contiguous
c     storage for labelling in vectors iprm,fprm
c
      integer iprm(16),mgopt(4)
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
      DATA MAXCYA/150/
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
      mgopt(3) = 1
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
      if (isolve >= 0) then ! called from dynamo
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
      else                 ! isolve < 0 -> called from dynamo_old
        write(6,"('>>> mud: invalid isolve=',i3)") isolve
        call shutdown('isolve')
c
c     set specified boundaries in phi
c
      endif ! isolve
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
      call mud2cr(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)
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
      call mud2cr(iprm,fprm,work,rhs,phi,mgopt,ierror,isolve)
      ier = ierror ! < 0 not converged
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
!-----------------------------------------------------------------------
c
c     file mud2cr.f  (version 4.0 modified for Cicley 2/99)
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 1998 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      MUDPACK version 4.0                    .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c
c
c ... author and specialist
c
c          John C. Adams (National Center for Atmospheric Research)
c          email: johnad@ucar.edu, phone: 303-497-1213

c ... For MUDPACK 4.0 information, visit the website:
c     (http://www.scd.ucar.edu/css/software/mudpack)
c
c ... purpose (see mud2cr.d for details)
c
c     mud2cr attempts to produce a second order finite difference
c     approximation to the two dimensional nonseparable elliptic
c     partial differential equation with cross derivative
c
c       cxx(x,y)*pxx + cxy(x,y)*pxy + cyy(x,y)*pyy +
c
c       cx(x,y)*px + cy(x,y)*py + ce(x,y)*pe(x,y) = r(x,y)
c
c ... documentation and test files
c
c     see the documentation file "mud2cr.d" for a complete discussion
c     of how to use subroutine mud2cr.  file "tmud2cr.f" is a test/driver
c     sample program illustrating use of mud2cr
c
c ... required MUDPACK files
c
c     mudcom.f
c
c
c
      subroutine mud2cr(iparm,fparm,work,rhs,phi,mgopt,
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
!
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
      call mud2cr1(nx,ny,rhs,phi,work)
      iparm(17) = itero
      if (tolmax.gt.0.0) then   ! check for convergence
	fparm(6) = relmax
	if (relmax.gt.tolmax) then
	   write(6,*) "no convergence with mud"
	   ierror = -1   ! flag convergenc failure
        end if
      end if
      return
      end
!-----------------------------------------------------------------------
      subroutine mud2cr1(nx,ny,rhsf,phif,wk)
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
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
      common/mud2crc/kpbgn(50),kcbgn(50),ktxbgn(50),ktybgn(50),
     +nxk(50),nyk(50),isx,jsy
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
!       write(6,"('mud2cr1: iter=',i3,' maxcy=',i3)") iter,maxcy
	itero = iter
	call kcymd2cr(wk)
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

      subroutine kcymd2cr(wk)
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
      call resmd2cr(nx,ny,wk(ip),ncx,ncy,wk(ipc),wk(irc),wk(ic),wk(kps))
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

      subroutine dismd2cr(nx,ny,cf,tx,ty,wk,ier,isolve)
      use pdynamo_module,only: nc,ncee,cee,ceee
c
c     discretize elliptic pde for mud2cr, set nonfatal errors
c
      implicit none
      integer,intent(in) :: isolve
!     integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
!    +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
!    +             kcycle,iprer,ipost,intpol,kps
      integer nx,ny,i,j,l,im1,jm1,ier,nnx,nny
      real cf(nx,ny,10),tx(nx,ny,*),ty(ny,nx,*)
      real wk(*)
      integer        intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps

      real           xa,xb,yc,yd,tolmax,relmax
      common/fmud2cr/xa,xb,yc,yd,tolmax,relmax
c
c     CHECK FOR CONSISTENCYT WRT KLEVEL
c
      NNX = ixp*2**(KLEVEL-1)+1
      NNY = jyq*2**(KLEVEL-1)+1
      IF(NNX.NE.NX.OR.NNY.NE.NY)THEN
	WRITE(6,100)NX,NY,NNX,NNY,ixp,jyq,KLEVEL
  100   FORMAT(' INCONSISTENCY WRT LEVEL. NX,NY,NNX,NNY,ixp,jyq,',
     |    'klevel = ',7I5)
	call shutdown('dismd2cr')
      ENDIF
      if (isolve >= 0) then
        CALL CEEE(CEE(NC(6-KLEVEL)),NX,NY,CF)
      else
        write(6,"('>>> mud: invalid isolve=',i3)") isolve
        call shutdown('isolve')
!       CALL CEEE(CEE0(nc_old(6-KLEVEL)),NX,NY,CF)
      endif
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

      subroutine adjmd2cr(nx,ny,phi,cf)
c
c     adjust righthand side in cf(i,j,10) for boundary conditions
c
      implicit none
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      real xa,xb,yc,yd,tolmax,relmax
      integer nx,ny,i,j
      real cf(nx,ny,10),phi(0:nx+1,0:ny+1)
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

      subroutine resmd2cr(nx,ny,phi,ncx,ncy,phic,rhsc,cof,resf)
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

      subroutine relmd2cr(nx,ny,phi,cof,tx,ty,sum)
c
c     relaxation for mud2
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
	call relmd2crp(nx,ny,phi,cof)
      else if (method.eq.1) then           ! line x relaxation
	call slxmd2cr(nx,ny,phi,cof,tx,sum)
      else if (method.eq.2) then           ! line y relaxation
	call slymd2cr(nx,ny,phi,cof,ty,sum)
      else if (method.eq.3) then           ! line x&y relaxation
	call slxmd2cr(nx,ny,phi,cof,tx,sum)
	call slymd2cr(nx,ny,phi,cof,ty,sum)
      end if
      return
      end

      subroutine relmd2crp(nx,ny,phi,cof)
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

      subroutine slxmd2cr(nx,ny,phi,cof,tx,sum)
c
c     line relaxation in the x direction (periodic or nonperiodic)
c
      implicit none
      integer nx,ny,i,ib,j,ii
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,iguess,
     +             maxcy,method,nwork,lwork,itero,ngrid,klevel,kcur,
     +             kcycle,iprer,ipost,intpol,kps
      common/imud2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nfx,nfy,
     +               iguess, maxcy,method,nwork,lwork,itero,ngrid,
     +               klevel,kcur,kcycle,iprer,ipost,intpol,kps
      real phi(0:nx+1,0:ny+1),cof(nx,ny,10),tx(nx,ny,*),sum(ny)
c
c     replace line x with point gauss-seidel if
c     x direction is periodic and nx = 3 (coarsest)
c
      if (nxa .eq. 0 .and. nx .eq. 3) then
	call relmd2crp(nx,ny,phi,cof)
	return
      end if
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

      subroutine slymd2cr(nx,ny,phi,cof,ty,sum)
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
c     replace line y with point gauss-seidel if
c     y direction is periodic and ny = 3
c
      if (nyc .eq. 0 .and. ny .eq. 3) then
	call relmd2crp(nx,ny,phi,cof)
	return
      end if
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


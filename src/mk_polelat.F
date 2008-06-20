!
      subroutine mk_polelat(jnew,jold,ixt)
!
! This software is part of the NCAR TIE-GCM.  Use is governed by the 
! Open Source Academic Research License Agreement contained in the file 
! tiegcmlicense.txt.
!
      use params_module,only: nlon,nlevp1
      use fields_module,only: nf4d,f4d
      use mpi_module,only: lon0,lon1
      implicit none
!
! Make boundary latitude jnew, using jold, i.e., 
!   read  f4d(*)%data(:,:,jold,ixt) and 
!   write f4d(*)%data(:,:,jnew,ixt)
! (jnew will be -1, 0, jmax+1, or jmax+2)
! (jold will be  2, 1, jmax  , or jmax-1)
!
! This is called only when not an MPI run (for MPI runs, see
!   mp_bndlats in mpi_module)
!
! Args:
      integer,intent(in) :: jnew,jold,ixt
!
! Locals:
      integer,parameter :: imaxh = nlon/2
      integer :: k,i,n,id,ilim
      real :: sign
!
! Fields loop:
      ilim = nlon/2+2
!     write(6,"(/,'mklatbnd: jnew=',i2,' jold=',i2,' lon0,1=',2i3,
!    |  ' ilim=',i3)") jnew,jold,lon0,lon1,ilim
      field_loop: do n=1,nf4d
        if (f4d(n)%polesign==0.) cycle field_loop
        do k=1,nlevp1
          id = imaxh+1
          if (f4d(n)%polesign==-1.) then
            do i=1,ilim
              f4d(n)%data(k,id+2,jnew,ixt) = 
     |       -f4d(n)%data(k,i+2 ,jold,ixt)
              f4d(n)%data(k,i   ,jnew,ixt) = 
     |       -f4d(n)%data(k,id  ,jold,ixt)
              id = id+1
            enddo
          else
            do i=1,ilim
              f4d(n)%data(k,id+2,jnew,ixt) = 
     |        f4d(n)%data(k,i+2 ,jold,ixt)
              f4d(n)%data(k,i   ,jnew,ixt) = 
     |        f4d(n)%data(k,id  ,jold,ixt)
              id = id+1
            enddo
          endif
        enddo
      enddo field_loop
      end subroutine mk_polelat

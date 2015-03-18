#ifdef CISMAH
!-----------------------------------------------------------------------
!
! Functions for using ad hoc coupling in CISM_LTR
!
!-------------------------------------------------------------------
       subroutine ca_init
!              
! This subroutine initialize ad hoc interface for data exchange 
!
! Open Fortran I/O units for M-I coupling outputs, will be replaced with 
! netcdf calls later
!
!
! Local:
!
      character(len=128) :: filename       
      integer :: npath
      logical isFileOK

      call get_file_path(npath,filename)
!
! Initialize ting status
!     
      filename(npath+1:npath+11)='ting_status'
      open(59,file=filename(1:npath+11))
      write(59,"(a1)")'0'
      close(59)     
!
!     End initialize ad hoc coupling 
!
      write(6,*) 'Ad hoc coupling Init Complete'

      end subroutine ca_init
!-------------------------------------------------------------------
      subroutine ca_close
!              
! This subroutine close cism ad hoc coupling interface for data exchange, and
! close Fortran I/O units 
! 
      write(6,*)'finished closing Fortran Units for coupling'

      end subroutine ca_close
!-------------------------------------------------------------------
      subroutine ca_receive
!              
! This subroutine receive auroral data via file exchanges from the
! M-I coupler 
!
      use params_module,only: nlat,nlonp1 
      use cism_coupling_module,only: gpot,geng,gflx,gpotm,
     |     ting_pot_interp,
     |     ting_eng_interp,
     |     ting_flux_interp

!     
! Local:
!
      character(len=128) :: filename
      integer :: npath
      integer ii,jj
      logical isFileOK

      call get_file_path(npath,filename)
!
! Check if data is available from M-I couplier
!
      close(57)
      filename(npath+1:npath+10)='mhd_status'
c     
c     check to make sure file exists before proceeding
c
      isFileOK = .FALSE.
      do while (isFileOK .NEQV. .TRUE.)
         inquire(file=filename(1:npath+10),exist=isFileOK)
#ifdef AIX
         call sleep_(5)
#else
         call sleep(5)
#endif
      enddo 

!
!     Lock file exists... Wait until mhd status is "new"
!
      mhd_data = 0

      do while(mhd_data.ne.1)
#ifdef AIX
       call sleep_(5)
#else
       call sleep(5)
#endif
       close(57)
       open(57,file=filename(1:npath+10))
       read(57,331)mhd_data
       close(57)
      enddo
!
! open data file 'ting-inputs.txt' and read data
!
      close(58)
       filename(npath+1:npath+15)='ting-inputs.txt'
      open(58,file=filename(1:npath+15))
      read(58,581)nday,nhour,nmin,nsec
581   format(i3,2x,i2,2x,i2,2x,i2)
      write(6,581)nday,nhour,nmin,nsec
      do j=1,nlat
       do i=1,nlonp1
         read(58,582)a,b,
     | ting_pot_interp(j,i),ting_eng_interp(j,i),ting_flux_interp(j,i)
       enddo
      enddo
582   format(5e12.4)
      close(58)  
!
! Change mhd status to old
!
      filename(npath+1:npath+10)='mhd_status'
      open(57,file=filename(1:npath+10))
      mhd_data=0
      write(57,331)mhd_data
      close(57)
331   format(1I2) 

!
      end subroutine ca_receive
!-------------------------------------------------------------------
      subroutine ca_send(cmodeltime)
!              
! This subroutine send TIEGCM data via ad hoc coupling to M-I coupler 
!
      use params_module,only: nlat,nlonp1,dlon,glon,glat
      use cism_coupling_module,only:
     |     ting_ped,
     |     ting_hall,
     |     ting_gnsrhs
!
! Args:
!
      integer,intent(in) :: cmodeltime(4)

!     
! Local:
!
      character(len=128) :: filename       
      integer :: npath
      real :: tie_send(nlat,nlonp1),gglon(nlonp1)
      integer :: jj,ii
      integer ting_data
      logical isFileOK

      call get_file_path(npath, filename)
!
! Don't bother checking lock files.  MIX reads & processes TIEGCM data
! before sending to TIEGCM.  So once TIEGCM receives data, the 
! TIEGCM export/send is already old.
!
!     Open ting-outputs.txt for writing.
!     
      close(56)
      filename(npath+1:npath+16)='ting-outputs.txt'
      open(56,file=filename(1:npath+16))

      write(56,333)cmodeltime(1),cmodeltime(2),
     1             cmodeltime(3),cmodeltime(4)

      do ii=1,nlonp1-1
         gglon(ii)=glon(ii)
      enddo
      gglon(nlonp1)=gglon(nlonp1-1)+dlon

!  FIXME:  Reading data into temporary arrays:
!          ting_ped,ting_hall,ting_gnsrhs
!          ... and then copying to real arrays
!          gzigm1, gzigm2, gnsrhs in cism_coupling.F
!          Could re-code to prevent temporary array slowness.
      do jj=1,nlat
         do ii=1,nlonp1
            write(56,561)glat(jj),gglon(ii),
!     |           gzigm1(ii,jj),gzigm2(ii,jj),gnsrhs(ii,jj) 
     |           ting_ped(jj,ii),ting_hall(jj,ii),ting_gnsrhs(jj,ii)
         enddo
      enddo
 561  format(5e12.4)

 333  format(i3,2x,i2,2x,i2,2x,i2)

      close(56)
      filename(npath+1:npath+11)='ting_status'
      open(56,file=filename(1:npath+11))     
      write(56,"(A1)")'1'
      close(56)     

      end subroutine ca_send
!-------------------------------------------------------------------
      subroutine get_file_path(npath,filename)
! ... Local variables ..................................................
      logical isFileOK
      character(len=60) :: ans
! ... Return variable ..................................................
      integer :: npath
      character(len=128) :: filename
! ... Begin ............................................................      
      inquire(file='path', exist=isFileOK)
      if (isFileOK .NEQV. .TRUE.) then
         write(*,*) "'path' file not found.  Couple in  $PWD"
         ans = "./"
      else      
         open(59,file='path')
         read(59,"(a)")ans
         close(59)
      endif
      npath=len_trim(ans)
      filename(1:npath)=ans(1:npath)
      end subroutine get_file_path
!-----------------------------------------------------------------------
#else
!-----------------------------------------------------------------------
      ! Intel Fortran compiler chokes on empty source files.  
      ! This subroutine is empty so this file will have SOMETHING in it
      subroutine adhoc_null
      end subroutine adhoc_null
!-----------------------------------------------------------------------
#endif

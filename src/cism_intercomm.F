#ifdef INTERCOMM
!-----------------------------------------------------------------------
!
! Functions for using InterComm coupling in CISM_LTR
!
!-----------------------------------------------------------------------
      subroutine ci_init
! ... Shared Module variables ..........................................
      use cism_coupling_module,only: 
     |     nlat,
     |     nlonp1,
     |     mytid,
     |     ntask,
     |     ting_pot_interp,
     |     ting_eng_interp,
     |     ting_flux_interp,
     |     ting_ped, 
     |     ting_hall, 
     |     ting_gnsrhs, 
     |     xjd      ! XJD handle

! ... Local variables ..................................................
      integer status ! InterComm Status flag

      ! IC_Create_bdecomp_desc variables:
      integer iaBlocks2d(1,2,2), tasks(1), iDesc2d

      ! IC_Create_block_region variables:
      integer iaLower2d(2), iaUpper2d(2), iaStride2d(2), iRegion2d

! ... Begin ............................................................      

      write(*,*) "Initializing TIEGCM InterComm Interface"
      call IC_Initialize("TIEGCM", mytid, "cmit.xjd", xjd, status)
      call handle_IC_error("cism_intercomm.F::IC_Initialize", status)

      tasks(1) = mytid
      iaBlocks2d(1,1,1) = 1
      iaBlocks2d(1,2,1) = nlat
      iaBlocks2d(1,1,2) = 1
      iaBlocks2d(1,2,2) = nlonp1

      call IC_Create_bdecomp_desc(2, iaBlocks2d, tasks, 1, iDesc2d, 1)
      !  IC_ROW_MAJOR = 0                                          ^^^
      !  IC_COLUMN_MAJOR  = 1                                      ^^^

      iaLower2d(1) = 1
      iaLower2d(2) = 1
      iaUpper2d(1) = nlat
      iaUpper2d(2) = nlonp1
      iaStride2d(1) = 1
      iaStride2d(2) = 1
      call IC_Create_block_region(2, iaLower2d, iaUpper2d, 
     |                            iaStride2d, iRegion2d)
      

      ! Register the "receive" arrays:
      call IC_Register_region(
     $     xjd, iRegion2d, 1, "pot", iDesc2d, ting_pot_interp, status)
      call handle_IC_error(
     $     "cism_intercomm.F::IC_Register_region(...pot...)", status)
      call IC_Register_region(
     $     xjd, iRegion2d, 1, "eng", iDesc2d, ting_eng_interp, status)
      call handle_IC_error(
     $     "cism_intercomm.F::IC_Register_region(...eng...)", status)
      call IC_Register_region(
     $     xjd, iRegion2d, 1, "flux", iDesc2d, ting_flux_interp, status)
      call handle_IC_error(
     $     "cism_intercomm.F::IC_Register_region(...flux...)", status)

      ! Register the "send" arrays:
      call IC_Register_region(
     $     xjd, iRegion2d, 1, "ped", iDesc2d, ting_ped, status)
      call handle_IC_error(
     $     "cism_intercomm.F::IC_Register_region(...ped...)", status)
      call IC_Register_region(
     $     xjd, iRegion2d, 1, "hall", iDesc2d, ting_hall, status)
      call handle_IC_error(
     $     "cism_intercomm.F::IC_Register_region(...hall...)", status)
      call IC_Register_region(
     $     xjd, iRegion2d, 1, "wfac", iDesc2d, ting_gnsrhs, status)
      call handle_IC_error(
     $     "cism_intercomm.F::IC_Register_region(...wfac...)", status)

      ! Commit all the regions:
      call IC_Commit_region(xjd, status)
      call handle_IC_error("cism_intercomm.F::IC_Commit_region", status)



      end subroutine ci_init
!-----------------------------------------------------------------------
      subroutine ci_close
! ... Shared Module variables ..........................................
      use cism_coupling_module,only: xjd
! ... Local Variables ..................................................
      integer status
! ... Begin ............................................................
      call IC_Finalize(xjd, status)
      call handle_IC_error("cism_intercomm.F::IC_Finalize", status)
      end subroutine ci_close
!-----------------------------------------------------------------------
      subroutine ci_receive
! ... Shared Module variables ..........................................
      use cism_coupling_module,only: 
     |     xjd      ! XJD handle
! ... Local Variables ..................................................
      integer status
! ... Begin ............................................................
      call IC_Import(xjd, "pot", status)
      call handle_IC_error("ITM_demo.F::IC_Import(...pot...)", status)
      call IC_Import(xjd, "eng", status)
      call handle_IC_error("ITM_demo.F::IC_Import(...eng...)", status)
      call IC_Import(xjd, "flux", status)
      call handle_IC_error("ITM_demo.F::IC_Import(...flux...)", status)

      end subroutine ci_receive
!-----------------------------------------------------------------------
      subroutine ci_send(modeltime)
! ... Shared Module variables ..........................................
      use cism_coupling_module,only: xjd
! ... Local Variables ..................................................
      integer status
! ... Parameter variables ..............................................
      integer, intent(in) :: modeltime(4)      
! ... Begin ............................................................
      call IC_Export(xjd, "ped", status)
      call handle_IC_error("ITM_demo.F::IC_Export(...ped...)", status)
      call IC_Export(xjd, "hall", status)
      call handle_IC_error("ITM_demo.F::IC_Export(...hall...)", status)
      call IC_Export(xjd, "wfac", status)
      call handle_IC_error("ITM_demo.F::IC_Export(...wfac...)", status)

      end subroutine ci_send
!-----------------------------------------------------------------------
      subroutine handle_IC_error(message, status)

      character message*(*)
      integer status

      if ( status .ne. 0 ) then
         write(*,"(A, A)") "Error in ", message
         write(*,"(A, I5)") "Status code: ", status
      end if

      end subroutine
!-----------------------------------------------------------------------
#else
!-----------------------------------------------------------------------
      ! Intel Fortran compiler chokes on empty source files.  
      ! This subroutine is empty so this file will have SOMETHING in it
      subroutine intercomm_null
      end subroutine intercomm_null
!-----------------------------------------------------------------------
#endif

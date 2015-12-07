      subroutine inp_read
!
! Namelist read. 
!
! This assumes the model executable has been given the path to
!   the namelist input file via a single program argument, e.g.:
!   mpirun -np $nproc $model $input >&! $output
! The path to the input file is accessed with fortran intrinsic getarg, 
!   then broadcast to all tasks, which open the file on unit 10, do the 
!   namelist read, and close the unit.
! This subroutine is in a separate file, instead of in input.F, to avoid
!   circular dependencies.
!
      use input_module,only: tgcm_input
      use mpi
      implicit none
      character(len=1024) :: inputfile
      integer :: ier
      
      inputfile = ' '
      call getarg(1,inputfile)
      write(6,"(/,'Reading namelist input data from ',a)") 
     |  trim(inputfile)

      call mpi_bcast(inputfile,1024,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

      open(10,file=trim(inputfile)) 
      read(10, nml=tgcm_input)
      close(10)
      write(6,"('Completed successful read of namelist inputs.',/)")
      end subroutine inp_read

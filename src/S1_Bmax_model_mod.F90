module S1_Bmax_model_mod

  use common_model_mod, only: eps,fit_by_MLT

  implicit none

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Arrays for fit coefficients:
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Number of fit coefficients for each mlat:
  integer :: nfit
  ! Number of mlat intervals:
  integer :: nmlat
  ! Fit coefficients (dimension nmlat x nfit):
  real, allocatable, dimension(:,:) :: b ! nmlat x nfit
  ! Interpolated fit coefficients (dimension nfit):
  real, allocatable, dimension(:) :: b_interp
  ! Fitting functions (dimension nfit):
  real, allocatable, dimension(:) :: x
  ! Fit coefficients averaged over the pole (dimension nfit):
  real, allocatable, dimension(:) :: b_90 
  ! for correction of coefficients
  integer, allocatable, dimension(:,:) :: num_reduc_coef
  ! Number of fit functions used by fit_data:
  integer, parameter :: number_of_MLT_coeffs = 1,&
       &number_of_Btrans_coeffs = 1,&
       &number_of_imf_angle_coeffs = 3,&
       &number_of_sinT_coeffs = 1
  ! Indices:
  integer :: imlat,ifit

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Data smoothing:
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  logical :: do_smooth
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! spline:
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter :: islpsw= 3                   ! flag for slope 
  real, parameter :: slp1 = 0., sigma = 1. ! slope at first point, tension factor
  real:: slpn= 0. ! slope at end point

contains
  function find_imf_angle_function (index,imf_angle_in) &
       &result (imf_angle_fxn_out)

    integer :: index
    real :: imf_angle_in,imf_angle_fxn_out

    select case (index)
    case (1)
       imf_angle_fxn_out = 1.
    case (2)
       imf_angle_fxn_out = cos(imf_angle_in)
    case (3)
       imf_angle_fxn_out = sin(imf_angle_in)
    case (4)
       imf_angle_fxn_out = cos(2.*imf_angle_in)
    case (5)
       imf_angle_fxn_out = sin(2.*imf_angle_in)
    end select
  end function find_imf_angle_function

  function find_MLT_function (index,MLT_in) &
       &result (MLT_fxn_out)

    use common_model_mod, only: pi

    integer :: index
    real :: MLT_in, MLT_fxn_out

    select case (index)
    case (1)
       MLT_fxn_out = 1.
    case (2)
       MLT_fxn_out =  cos(pi*MLT_in/12.)
    case (3)
       MLT_fxn_out = sin(pi*MLT_in/12.)
    case (4)
       MLT_fxn_out = cos(pi*MLT_in/6.)
    case (5)
       MLT_fxn_Out = sin(pi*MLT_in/6.)
    end select
  end function find_MLT_function

  function find_sinT_function (index,sinT_in) &
       &result (sinT_fxn_out)

    integer :: index
    real :: sinT_in, sinT_fxn_out

    select case (index)
    case (0)
       sinT_fxn_out = 1.
    case (1) 
       sinT_fxn_out = 1.
    case (2)
       sinT_fxn_out = sinT_in
    case (3)
       sinT_fxn_out = sinT_in**2
    end select
  end function find_sinT_function

end module S1_Bmax_model_mod

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

module S1_Bmax_coeff_mod

  ! This module file created 2007/07/17.
  ! with the following parameters:
  !   for 81252-83047 period 
  !   fit_by_Btrans = 0
  !   fit_by_imf_angle = 1
  !   fit_by_sinT = 1
  !   mlat_resolution = 5

  ! To use this module, include the statement
  !   use coeff_mod
  ! and the subroutine call
  !   call fill_coeff_arrays
  ! in Fortran90 programs.
  ! 
  ! Compilation will be slow; a Makefile is suggested!
  implicit none
  save
  real,  dimension(  9,  3) :: tot_Poynting_Wm2_coeffs
  real,  dimension(  9) :: min_mlat,max_mlat,average_mlat
  integer, parameter :: &
       & fit_by_Btrans = 0, &
       & fit_by_imf_angle = 1, &
       & fit_by_sinT = 0, &
       & mlat_resolution = 5

contains

  subroutine S1_Bmax_fill_coeff_arrays

    ! Fill tot_Poynting_Wm2 arrays:
    tot_Poynting_Wm2_coeffs(  1,  1) =  -0.1129E-04
    tot_Poynting_Wm2_coeffs(  1,  2) =  -0.2730E-06
    tot_Poynting_Wm2_coeffs(  1,  3) =  -0.8407E-05
    tot_Poynting_Wm2_coeffs(  2,  1) =   0.3074E-05
    tot_Poynting_Wm2_coeffs(  2,  2) =  -0.1855E-04
    tot_Poynting_Wm2_coeffs(  2,  3) =  -0.2917E-04
    tot_Poynting_Wm2_coeffs(  3,  1) =   0.1652E-03
    tot_Poynting_Wm2_coeffs(  3,  2) =  -0.2039E-03
    tot_Poynting_Wm2_coeffs(  3,  3) =   0.4913E-04
    tot_Poynting_Wm2_coeffs(  4,  1) =   0.6122E-03
    tot_Poynting_Wm2_coeffs(  4,  2) =  -0.1525E-03
    tot_Poynting_Wm2_coeffs(  4,  3) =  -0.4053E-03
    tot_Poynting_Wm2_coeffs(  5,  1) =   0.2216E-02
    tot_Poynting_Wm2_coeffs(  5,  2) =  -0.2831E-02
    tot_Poynting_Wm2_coeffs(  5,  3) =  -0.1412E-04
    tot_Poynting_Wm2_coeffs(  6,  1) =   0.2869E-02
    tot_Poynting_Wm2_coeffs(  6,  2) =  -0.3530E-02
    tot_Poynting_Wm2_coeffs(  6,  3) =   0.7533E-03
    tot_Poynting_Wm2_coeffs(  7,  1) =   0.1714E-02
    tot_Poynting_Wm2_coeffs(  7,  2) =  -0.3026E-03
    tot_Poynting_Wm2_coeffs(  7,  3) =   0.1705E-02
    tot_Poynting_Wm2_coeffs(  8,  1) =   0.1245E-02
    tot_Poynting_Wm2_coeffs(  8,  2) =   0.7075E-03
    tot_Poynting_Wm2_coeffs(  8,  3) =   0.1013E-03
    tot_Poynting_Wm2_coeffs(  9,  1) =   0.1240E-02
    tot_Poynting_Wm2_coeffs(  9,  2) =   0.2107E-02
    tot_Poynting_Wm2_coeffs(  9,  3) =   0.8019E-03

    ! Fill mlat arrays:
    min_mlat(  1) =     45.00
    max_mlat(  1) =     50.00
    min_mlat(  2) =     50.00
    max_mlat(  2) =     55.00
    min_mlat(  3) =     55.00
    max_mlat(  3) =     60.00
    min_mlat(  4) =     60.00
    max_mlat(  4) =     65.00
    min_mlat(  5) =     65.00
    max_mlat(  5) =     70.00
    min_mlat(  6) =     70.00
    max_mlat(  6) =     75.00
    min_mlat(  7) =     75.00
    max_mlat(  7) =     80.00
    min_mlat(  8) =     80.00
    max_mlat(  8) =     85.00
    min_mlat(  9) =     85.00
    max_mlat(  9) =     90.00
    average_mlat = 0.5*(min_mlat+max_mlat)

  end subroutine S1_Bmax_fill_coeff_arrays

end module S1_Bmax_coeff_mod

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine S1_Bmax_choose_coeff_array

  use S1_Bmax_coeff_mod, only: tot_Poynting_Wm2_coeffs,&
       & S1_Bmax_fill_coeff_arrays
  use S1_Bmax_model_mod, only: nmlat,nfit,b,b_interp,&
       & b_90,x,imlat,ifit,eps,num_reduc_coef
  integer :: binary_fit_by

  !****************************************************************
  ! Fill parameter arrays:
  !****************************************************************
  call S1_Bmax_fill_coeff_arrays

  !****************************************************************
  ! Determine array dimensions.  Arrays for data type 'potential'
  ! are the same size as arrays for other data types.
  !****************************************************************

  nmlat = size(tot_Poynting_Wm2_coeffs,1)
  nfit  = size(tot_Poynting_Wm2_coeffs,2)

  !****************************************************************
  ! Allocate arrays into which coefficients will be transferred:
  !****************************************************************
  allocate(b(nmlat,nfit),b_interp(nfit),x(nfit),b_90(nfit), &
       &  num_reduc_coef(nmlat,nfit))

  !****************************************************************
  ! read in the correction factors
  !****************************************************************

  num_reduc_coef(:,1)  = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,2)  = (/ 0,0,2,0,5,4,0,0,7 /)
  num_reduc_coef(:,3)  = (/ 0,4,0,0,0,0,0,0,0 /)

  !****************************************************************
  ! Transfer values:
  !****************************************************************

  do imlat=1,nmlat
     do ifit=1,nfit
        b(imlat,ifit) = tot_Poynting_Wm2_coeffs(imlat,ifit)* &
             & (1.-eps)**num_reduc_coef(imlat,ifit)
     end do
  end do

end subroutine S1_Bmax_choose_coeff_array

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine S1_Bmax_calc_model_value(mlat_glob,mlt,sinT,imf_angle,&
     &value,value_found)

  ! This version is for public distribution, along with the module
  ! file coeff_mod.f90.
  !
  ! Within the /home/hensel/fortran/fit_data/public subdirectory,
  ! this program is called by call_de2_model.

  use S1_Bmax_coeff_mod
  use S1_Bmax_model_mod

  implicit none

  real, intent(in) :: mlat_glob,mlt,sinT,imf_angle
  real, intent(out) :: value
  logical, intent(out) :: value_found

  integer :: fit_count,ierr
  logical :: mlat_found
  integer :: iMLT,iBtrans,iimf_angle,isinT
  real :: MLT_function,Btrans_function,imf_angle_function,sinT_function,&
       & imf_angle_in, avg_Btrans_array,mlat
  real, dimension(nmlat) :: b_secderiv,temp,b_coeff,lat_array

  ! General:
  logical, parameter :: debug = .false.
  real, parameter :: pi = 3.1415926535898
  real, parameter :: beta_inf = 16.
  real, external :: curv2 

  !****************************************************************
  ! Interpolate to get coefficients for specified mlat:
  !****************************************************************


  mlat_found = .false.
  value_found = .false.

  ! Use absolute value of mlat:
  mlat = abs(mlat_glob)
  lat_array(1:nmlat) = average_mlat(1:nmlat)

  interpolate_coefficients:do ifit = 1,nfit
     b_coeff(1:nmlat) = b(1:nmlat,ifit)

     call curv1(nmlat,lat_array,b_coeff,slp1,slpn, &
          &    islpsw,b_secderiv,temp,sigma,ierr) 
     b_interp(ifit) = curv2(mlat,nmlat,lat_array,&
          & b_coeff,b_secderiv,sigma)

  end do interpolate_coefficients
  mlat_found = .true.


  !****************************************************************
  ! Define which x function to use for each fit coefficient.  As in
  ! the define_function_indices subroutine of fit_data, cycle
  ! through sinT most rapidly; cycle through MLT most slowly.
  !****************************************************************

  inside_mlat_range:if (mlat_found) then
     x = 1.
     fit_count = 1

     MLT_loop2:do iMLT=1,1+(fit_by_MLT*number_of_MLT_coeffs -&
          & fit_by_MLT)
        MLT_function = find_MLT_function(iMLT,MLT)

        Btrans_loop2:do iBtrans=1,1+(fit_by_Btrans&
             &*number_of_Btrans_coeffs - fit_by_Btrans)
           imf_angle_loop2:do iimf_angle=1,1+(fit_by_imf_angle&
                &*number_of_imf_angle_coeffs - fit_by_imf_angle)
              imf_angle_function = find_imf_angle_function(iimf_angle,imf_angle)

              sinT_loop2:do isinT=1,1+(fit_by_sinT&
                   &*number_of_sinT_coeffs - fit_by_sinT) 
                 sinT_function = find_sinT_function(isinT,sinT)

                 x(fit_count) = MLT_function&
                      &*imf_angle_function*sinT_function

                 if (debug) then
                    write(*,'("fit_count = ",I5)') fit_count
                    write(*,'("MLT_function = ",G12.6)') MLT_function
                    write(*,'("Btrans_function = ",G12.6)') Btrans_function
                    write(*,'("imf_angle_function = ",G12.6)') &
                         &imf_angle_function
                    write(*,'("sinT_function = ",G12.6)')&
                         & sinT_function
                    write(*,'("mlat = ",F9.2," MLT = ",F9.2," Btrans = ",&
                         &a," imf_angle = ",F9.2," sinT = ",F9.2,&
                         &" x(",I5,") = ",G13.6)') &
                         &mlat,MLT,'n/a',imf_angle,sinT,&
                         &fit_count,x(fit_count)
                 end if

                 fit_count = fit_count + 1

              end do sinT_loop2
           end do imf_angle_loop2
        end do Btrans_loop2
     end do MLT_loop2

     !****************************************************************
     ! Calculate value with interpolated coefficients:
     !****************************************************************

     value = 0.
     go_through_coefficients:do ifit=1,nfit
        value = value + x(ifit)*b_interp(ifit)
     end do go_through_coefficients
     value_found = .true.

  end if inside_mlat_range

end subroutine S1_Bmax_calc_model_value

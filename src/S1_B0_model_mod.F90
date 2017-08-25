module S1_B0_model_mod

  use common_model_mod, only: fit_by_MLT, eps

  implicit none

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Specify which fit results to use:
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Data is always fitted by MLT:
!!$  integer, parameter :: fit_by_MLT = 1 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Arrays for fit coefficients:
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Number of fit coefficients for each mlat:
  integer :: nfit
  ! Number of mlat intervals:
  integer :: nmlat
  ! Fit coefficients (dimension nmlat x nfit):
  real, allocatable, dimension(:,:) :: b 
  ! Interpolated fit coefficients (dimension nfit):
  real, allocatable, dimension(:) :: b_interp
  ! Fitting functions (dimension nfit):
  real, allocatable, dimension(:) :: x
  ! Fit coefficients averaged over the pole (dimension nfit):
  real, allocatable, dimension(:) :: b_90 
  ! for correction of coefficients
  integer, allocatable, dimension(:,:) :: num_reduc_coef
  ! Number of fit functions used by fit_data:
  integer, parameter :: number_of_MLT_coeffs = 5,&
       &number_of_Btrans_coeffs = 2,&
       &number_of_imf_angle_coeffs = 1,&
       &number_of_sinT_coeffs = 2
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
       sinT_fxn_out = 0.5*(1+(4.*sinT_in+0.5)/SQRT(0.75+(4.*sinT_in+0.5)**2))
    case (2)
       sinT_fxn_out = 0.5*(1+(4.*sinT_in+0.5)/SQRT(0.75+(4.*sinT_in+0.5)**2))
       sinT_fxn_out = sinT_in * sinT_fxn_out
    case (3)
       sinT_fxn_out = sinT_in**2
    end select
  end function find_sinT_function

  function find_Btrans_function (index,Btrans_in) &
       &result (Btrans_fxn_out)

    integer :: index
    real :: Btrans_in, Btrans_fxn_out
    real, parameter :: beta_inf = 16

    select case (index)
    case (1)
       Btrans_fxn_out = 1
    case (2)
       Btrans_fxn_out = Btrans_in /sqrt(1.+(Btrans_in*Btrans_in)/(beta_inf*beta_inf))
    case (3)
       Btrans_fxn_out = (Btrans_in/5.2 - 1.)**2
    end select
  end function find_Btrans_function

end module S1_B0_model_mod

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

module S1_B0_coeff_mod

  ! This module file created 2007/07/24.
  ! with the following parameters:
  !   for 81252-83047 period 
  !   fit_by_Btrans = 1
  !   fit_by_imf_angle = 0
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
  real,  dimension(  9, 20) :: tot_Poynting_Wm2_coeffs
  real,  dimension(  9) :: min_mlat,max_mlat,average_mlat
  integer, parameter :: &
       & fit_by_Btrans = 1, &
       & fit_by_imf_angle = 0, &
       & fit_by_sinT = 1, &
       & mlat_resolution = 5

contains

  subroutine S1_B0_fill_coeff_arrays

    ! Fill tot_Poynting_Wm2 arrays:
    tot_Poynting_Wm2_coeffs(  1,  1) =   0.1923E-05
    tot_Poynting_Wm2_coeffs(  1,  2) =   0.2683E-04
    tot_Poynting_Wm2_coeffs(  1,  3) =  -0.2964E-05
    tot_Poynting_Wm2_coeffs(  1,  4) =  -0.1736E-04
    tot_Poynting_Wm2_coeffs(  1,  5) =  -0.2569E-04
    tot_Poynting_Wm2_coeffs(  1,  6) =   0.2419E-04
    tot_Poynting_Wm2_coeffs(  1,  7) =  -0.4999E-05
    tot_Poynting_Wm2_coeffs(  1,  8) =   0.2534E-04
    tot_Poynting_Wm2_coeffs(  1,  9) =  -0.6518E-05
    tot_Poynting_Wm2_coeffs(  1, 10) =  -0.1529E-04
    tot_Poynting_Wm2_coeffs(  1, 11) =  -0.9450E-05
    tot_Poynting_Wm2_coeffs(  1, 12) =   0.1751E-04
    tot_Poynting_Wm2_coeffs(  1, 13) =  -0.4402E-06
    tot_Poynting_Wm2_coeffs(  1, 14) =   0.7415E-11
    tot_Poynting_Wm2_coeffs(  1, 15) =   0.3702E-05
    tot_Poynting_Wm2_coeffs(  1, 16) =  -0.5246E-11
    tot_Poynting_Wm2_coeffs(  1, 17) =   0.2410E-04
    tot_Poynting_Wm2_coeffs(  1, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  1, 19) =  -0.4875E-05
    tot_Poynting_Wm2_coeffs(  1, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  2,  1) =  -0.8373E-05
    tot_Poynting_Wm2_coeffs(  2,  2) =  -0.1339E-03
    tot_Poynting_Wm2_coeffs(  2,  3) =   0.5384E-05
    tot_Poynting_Wm2_coeffs(  2,  4) =  -0.1719E-04
    tot_Poynting_Wm2_coeffs(  2,  5) =   0.3842E-06
    tot_Poynting_Wm2_coeffs(  2,  6) =   0.1221E-03
    tot_Poynting_Wm2_coeffs(  2,  7) =  -0.1379E-04
    tot_Poynting_Wm2_coeffs(  2,  8) =   0.9573E-05
    tot_Poynting_Wm2_coeffs(  2,  9) =  -0.2105E-04
    tot_Poynting_Wm2_coeffs(  2, 10) =   0.9059E-04
    tot_Poynting_Wm2_coeffs(  2, 11) =  -0.1007E-04
    tot_Poynting_Wm2_coeffs(  2, 12) =  -0.2607E-05
    tot_Poynting_Wm2_coeffs(  2, 13) =  -0.4181E-04
    tot_Poynting_Wm2_coeffs(  2, 14) =  -0.4751E-09
    tot_Poynting_Wm2_coeffs(  2, 15) =   0.6914E-05
    tot_Poynting_Wm2_coeffs(  2, 16) =  -0.7139E-10
    tot_Poynting_Wm2_coeffs(  2, 17) =  -0.1873E-04
    tot_Poynting_Wm2_coeffs(  2, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  2, 19) =   0.1102E-04
    tot_Poynting_Wm2_coeffs(  2, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  3,  1) =   0.7277E-05
    tot_Poynting_Wm2_coeffs(  3,  2) =  -0.2646E-03
    tot_Poynting_Wm2_coeffs(  3,  3) =   0.1549E-05
    tot_Poynting_Wm2_coeffs(  3,  4) =  -0.1693E-05
    tot_Poynting_Wm2_coeffs(  3,  5) =   0.2908E-04
    tot_Poynting_Wm2_coeffs(  3,  6) =   0.7187E-04
    tot_Poynting_Wm2_coeffs(  3,  7) =  -0.2628E-04
    tot_Poynting_Wm2_coeffs(  3,  8) =   0.4439E-04
    tot_Poynting_Wm2_coeffs(  3,  9) =  -0.1147E-04
    tot_Poynting_Wm2_coeffs(  3, 10) =   0.2536E-03
    tot_Poynting_Wm2_coeffs(  3, 11) =  -0.2175E-04
    tot_Poynting_Wm2_coeffs(  3, 12) =   0.2192E-04
    tot_Poynting_Wm2_coeffs(  3, 13) =  -0.2259E-04
    tot_Poynting_Wm2_coeffs(  3, 14) =  -0.4786E-09
    tot_Poynting_Wm2_coeffs(  3, 15) =   0.9214E-05
    tot_Poynting_Wm2_coeffs(  3, 16) =  -0.1188E-09
    tot_Poynting_Wm2_coeffs(  3, 17) =  -0.3691E-04
    tot_Poynting_Wm2_coeffs(  3, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  3, 19) =   0.3797E-05
    tot_Poynting_Wm2_coeffs(  3, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  4,  1) =  -0.2921E-03
    tot_Poynting_Wm2_coeffs(  4,  2) =   0.7416E-03
    tot_Poynting_Wm2_coeffs(  4,  3) =   0.4973E-03
    tot_Poynting_Wm2_coeffs(  4,  4) =  -0.1448E-02
    tot_Poynting_Wm2_coeffs(  4,  5) =   0.2104E-03
    tot_Poynting_Wm2_coeffs(  4,  6) =  -0.1566E-04
    tot_Poynting_Wm2_coeffs(  4,  7) =   0.6384E-04
    tot_Poynting_Wm2_coeffs(  4,  8) =   0.7207E-04
    tot_Poynting_Wm2_coeffs(  4,  9) =   0.4205E-03
    tot_Poynting_Wm2_coeffs(  4, 10) =  -0.1319E-02
    tot_Poynting_Wm2_coeffs(  4, 11) =   0.6828E-04
    tot_Poynting_Wm2_coeffs(  4, 12) =  -0.2248E-03
    tot_Poynting_Wm2_coeffs(  4, 13) =  -0.1979E-04
    tot_Poynting_Wm2_coeffs(  4, 14) =   0.2592E-07
    tot_Poynting_Wm2_coeffs(  4, 15) =  -0.1893E-03
    tot_Poynting_Wm2_coeffs(  4, 16) =  -0.3674E-09
    tot_Poynting_Wm2_coeffs(  4, 17) =  -0.8901E-04
    tot_Poynting_Wm2_coeffs(  4, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  4, 19) =   0.1921E-03
    tot_Poynting_Wm2_coeffs(  4, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  5,  1) =   0.1543E-02
    tot_Poynting_Wm2_coeffs(  5,  2) =  -0.3217E-02
    tot_Poynting_Wm2_coeffs(  5,  3) =   0.5178E-03
    tot_Poynting_Wm2_coeffs(  5,  4) =  -0.1618E-02
    tot_Poynting_Wm2_coeffs(  5,  5) =   0.2300E-02
    tot_Poynting_Wm2_coeffs(  5,  6) =  -0.4980E-02
    tot_Poynting_Wm2_coeffs(  5,  7) =  -0.6277E-04
    tot_Poynting_Wm2_coeffs(  5,  8) =   0.8824E-03
    tot_Poynting_Wm2_coeffs(  5,  9) =   0.1643E-02
    tot_Poynting_Wm2_coeffs(  5, 10) =  -0.4623E-02
    tot_Poynting_Wm2_coeffs(  5, 11) =   0.2627E-03
    tot_Poynting_Wm2_coeffs(  5, 12) =  -0.5532E-03
    tot_Poynting_Wm2_coeffs(  5, 13) =  -0.4375E-04
    tot_Poynting_Wm2_coeffs(  5, 14) =   0.2838E-07
    tot_Poynting_Wm2_coeffs(  5, 15) =  -0.2815E-03
    tot_Poynting_Wm2_coeffs(  5, 16) =   0.5445E-10
    tot_Poynting_Wm2_coeffs(  5, 17) =   0.6141E-03
    tot_Poynting_Wm2_coeffs(  5, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  5, 19) =   0.1871E-04
    tot_Poynting_Wm2_coeffs(  5, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  6,  1) =   0.2832E-02
    tot_Poynting_Wm2_coeffs(  6,  2) =  -0.2433E-02
    tot_Poynting_Wm2_coeffs(  6,  3) =   0.1372E-03
    tot_Poynting_Wm2_coeffs(  6,  4) =  -0.9869E-03
    tot_Poynting_Wm2_coeffs(  6,  5) =  -0.5893E-03
    tot_Poynting_Wm2_coeffs(  6,  6) =  -0.2283E-02
    tot_Poynting_Wm2_coeffs(  6,  7) =   0.1472E-03
    tot_Poynting_Wm2_coeffs(  6,  8) =   0.8081E-03
    tot_Poynting_Wm2_coeffs(  6,  9) =  -0.6501E-04
    tot_Poynting_Wm2_coeffs(  6, 10) =  -0.4668E-02
    tot_Poynting_Wm2_coeffs(  6, 11) =  -0.2429E-03
    tot_Poynting_Wm2_coeffs(  6, 12) =   0.8933E-03
    tot_Poynting_Wm2_coeffs(  6, 13) =  -0.1854E-02
    tot_Poynting_Wm2_coeffs(  6, 14) =   0.2509E-07
    tot_Poynting_Wm2_coeffs(  6, 15) =   0.1698E-03
    tot_Poynting_Wm2_coeffs(  6, 16) =  -0.1354E-08
    tot_Poynting_Wm2_coeffs(  6, 17) =   0.1535E-02
    tot_Poynting_Wm2_coeffs(  6, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  6, 19) =  -0.8631E-04
    tot_Poynting_Wm2_coeffs(  6, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  7,  1) =   0.1517E-02
    tot_Poynting_Wm2_coeffs(  7,  2) =  -0.1098E-02
    tot_Poynting_Wm2_coeffs(  7,  3) =   0.4053E-03
    tot_Poynting_Wm2_coeffs(  7,  4) =  -0.4557E-03
    tot_Poynting_Wm2_coeffs(  7,  5) =  -0.1806E-02
    tot_Poynting_Wm2_coeffs(  7,  6) =   0.2114E-03
    tot_Poynting_Wm2_coeffs(  7,  7) =  -0.1706E-03
    tot_Poynting_Wm2_coeffs(  7,  8) =  -0.6980E-05
    tot_Poynting_Wm2_coeffs(  7,  9) =   0.5251E-03
    tot_Poynting_Wm2_coeffs(  7, 10) =  -0.2451E-02
    tot_Poynting_Wm2_coeffs(  7, 11) =   0.9194E-04
    tot_Poynting_Wm2_coeffs(  7, 12) =   0.5436E-03
    tot_Poynting_Wm2_coeffs(  7, 13) =  -0.4386E-03
    tot_Poynting_Wm2_coeffs(  7, 14) =   0.4865E-08
    tot_Poynting_Wm2_coeffs(  7, 15) =   0.1722E-03
    tot_Poynting_Wm2_coeffs(  7, 16) =  -0.5843E-09
    tot_Poynting_Wm2_coeffs(  7, 17) =   0.6265E-03
    tot_Poynting_Wm2_coeffs(  7, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  7, 19) =  -0.1569E-03
    tot_Poynting_Wm2_coeffs(  7, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  8,  1) =  -0.1021E-02
    tot_Poynting_Wm2_coeffs(  8,  2) =   0.5664E-02
    tot_Poynting_Wm2_coeffs(  8,  3) =   0.1060E-02
    tot_Poynting_Wm2_coeffs(  8,  4) =  -0.9146E-03
    tot_Poynting_Wm2_coeffs(  8,  5) =  -0.1043E-03
    tot_Poynting_Wm2_coeffs(  8,  6) =  -0.7324E-02
    tot_Poynting_Wm2_coeffs(  8,  7) =  -0.3407E-03
    tot_Poynting_Wm2_coeffs(  8,  8) =  -0.7892E-03
    tot_Poynting_Wm2_coeffs(  8,  9) =  -0.7853E-03
    tot_Poynting_Wm2_coeffs(  8, 10) =   0.3407E-02
    tot_Poynting_Wm2_coeffs(  8, 11) =   0.8431E-03
    tot_Poynting_Wm2_coeffs(  8, 12) =  -0.1810E-02
    tot_Poynting_Wm2_coeffs(  8, 13) =   0.1130E-02
    tot_Poynting_Wm2_coeffs(  8, 14) =   0.7713E-08
    tot_Poynting_Wm2_coeffs(  8, 15) =   0.4456E-03
    tot_Poynting_Wm2_coeffs(  8, 16) =  -0.1677E-08
    tot_Poynting_Wm2_coeffs(  8, 17) =   0.1534E-03
    tot_Poynting_Wm2_coeffs(  8, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  8, 19) =  -0.2010E-03
    tot_Poynting_Wm2_coeffs(  8, 20) =    0.000    
    tot_Poynting_Wm2_coeffs(  9,  1) =  -0.1032E-02
    tot_Poynting_Wm2_coeffs(  9,  2) =  -0.1406E-02
    tot_Poynting_Wm2_coeffs(  9,  3) =   0.3905E-03
    tot_Poynting_Wm2_coeffs(  9,  4) =   0.1519E-03
    tot_Poynting_Wm2_coeffs(  9,  5) =   0.7727E-03
    tot_Poynting_Wm2_coeffs(  9,  6) =   0.1404E-02
    tot_Poynting_Wm2_coeffs(  9,  7) =  -0.8205E-03
    tot_Poynting_Wm2_coeffs(  9,  8) =   0.3567E-03
    tot_Poynting_Wm2_coeffs(  9,  9) =  -0.4756E-04
    tot_Poynting_Wm2_coeffs(  9, 10) =   0.1103E-02
    tot_Poynting_Wm2_coeffs(  9, 11) =  -0.9804E-04
    tot_Poynting_Wm2_coeffs(  9, 12) =   0.7322E-04
    tot_Poynting_Wm2_coeffs(  9, 13) =   0.4321E-04
    tot_Poynting_Wm2_coeffs(  9, 14) =  -0.8882E-09
    tot_Poynting_Wm2_coeffs(  9, 15) =   0.3519E-03
    tot_Poynting_Wm2_coeffs(  9, 16) =  -0.9275E-11
    tot_Poynting_Wm2_coeffs(  9, 17) =  -0.3221E-03
    tot_Poynting_Wm2_coeffs(  9, 18) =    0.000    
    tot_Poynting_Wm2_coeffs(  9, 19) =   0.4175E-04
    tot_Poynting_Wm2_coeffs(  9, 20) =    0.000    

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

  end subroutine S1_B0_fill_coeff_arrays

end module S1_B0_coeff_mod

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine S1_B0_choose_coeff_array

  use S1_B0_coeff_mod, only: tot_Poynting_Wm2_coeffs,&
       & S1_B0_fill_coeff_arrays
  use S1_B0_model_mod, only: nmlat,nfit,b,b_interp,&
       & b_90,x,imlat,ifit,eps,num_reduc_coef
  integer :: binary_fit_by

  !****************************************************************
  ! Fill parameter arrays:
  !****************************************************************
  call S1_B0_fill_coeff_arrays

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
  num_reduc_coef(:,2)  = (/ 0,0,5,0,9,6,3,0,7 /)
  num_reduc_coef(:,3)  = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,4)  = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,5)  = (/ 0,0,2,43,5,6,3,7,7 /)
  num_reduc_coef(:,6)  = (/ 0,0,5,0,4,6,0,7,7 /)
  num_reduc_coef(:,7)  = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,8)  = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,9)  = (/ 0,0,0,40,3,0,2,36,0 /)
  num_reduc_coef(:,10) = (/ 0,0,5,8,6,4,1,7,7 /)
  num_reduc_coef(:,11) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,12) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,13) = (/ 0,0,1,0,4,6,3,36,0 /)
  num_reduc_coef(:,14) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,15) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,16) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,17) = (/ 0,0,5,25,8,2,2,11,4 /)
  num_reduc_coef(:,18) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,19) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,20) = (/ 0,0,0,0,0,0,0,0,0 /)

  !****************************************************************
  ! Transfer values:
  !****************************************************************

  do imlat=1,nmlat
     do ifit=1,nfit
        b(imlat,ifit) = tot_Poynting_Wm2_coeffs(imlat,ifit)* &
             & (1.-eps)**num_reduc_coef(imlat,ifit)
     end do
  end do

end subroutine S1_B0_choose_coeff_array

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine S1_B0_calc_model_value(mlat_glob,mlt,sinT,imf_angle,Btrans,&
     &value,value_found)

  ! This version is for public distribution, along with the module
  ! file coeff_mod.f90.
  !
  ! Within the /home/hensel/fortran/fit_data/public subdirectory,
  ! this program is called by call_de2_model.

  use S1_B0_coeff_mod
  use S1_B0_model_mod

  implicit none

  real, intent(in) :: mlat_glob,mlt,sinT,imf_angle,Btrans
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
           Btrans_function = find_Btrans_function(iBtrans,Btrans)

           imf_angle_loop2:do iimf_angle=1,1+(fit_by_imf_angle&
                &*number_of_imf_angle_coeffs - fit_by_imf_angle)
              imf_angle_function = find_imf_angle_function(iimf_angle,imf_angle)

              sinT_loop2:do isinT=1,1+(fit_by_sinT&
                   &*number_of_sinT_coeffs - fit_by_sinT) 
                 sinT_function = find_sinT_function(isinT,sinT)

                 x(fit_count) = MLT_function*Btrans_function&
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
                         &F9.2," imf_angle = ",F9.2," sinT = ",F9.2,&
                         &" x(",I5,") = ",G13.6)') &
                         &mlat,MLT,Btrans,imf_angle,sinT,&
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

end subroutine S1_B0_calc_model_value

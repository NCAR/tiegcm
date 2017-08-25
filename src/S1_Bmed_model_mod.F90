module S1_Bmed_model_mod

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
  integer, parameter :: number_of_MLT_coeffs = 5,&
       &number_of_Btrans_coeffs = 1,&
       &number_of_imf_angle_coeffs = 3,&
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
       sinT_fxn_out = 1.
    case (2)
       sinT_fxn_out = sinT_in
    case (3)
       sinT_fxn_out = sinT_in**2
    end select
  end function find_sinT_function

end module S1_Bmed_model_mod

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

module S1_Bmed_coeff_mod

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
  real,  dimension(  9, 30) :: tot_Poynting_Wm2_coeffs
  real,  dimension(  9) :: min_mlat,max_mlat,average_mlat
  integer, parameter :: &
       & fit_by_Btrans = 0, &
       & fit_by_imf_angle = 1, &
       & fit_by_sinT = 1, &
       & mlat_resolution = 5

contains

  subroutine S1_Bmed_fill_coeff_arrays

    ! Fill tot_Poynting_Wm2 arrays:
    tot_Poynting_Wm2_coeffs(  1,  1) =  -0.3974E-05
    tot_Poynting_Wm2_coeffs(  1,  2) =   0.7041E-05
    tot_Poynting_Wm2_coeffs(  1,  3) =   0.3835E-05
    tot_Poynting_Wm2_coeffs(  1,  4) =  -0.2953E-04
    tot_Poynting_Wm2_coeffs(  1,  5) =  -0.8606E-05
    tot_Poynting_Wm2_coeffs(  1,  6) =   0.6838E-06
    tot_Poynting_Wm2_coeffs(  1,  7) =  -0.8643E-05
    tot_Poynting_Wm2_coeffs(  1,  8) =   0.3111E-04
    tot_Poynting_Wm2_coeffs(  1,  9) =   0.6500E-06
    tot_Poynting_Wm2_coeffs(  1, 10) =  -0.6653E-04
    tot_Poynting_Wm2_coeffs(  1, 11) =  -0.9747E-05
    tot_Poynting_Wm2_coeffs(  1, 12) =   0.1496E-04
    tot_Poynting_Wm2_coeffs(  1, 13) =  -0.1066E-04
    tot_Poynting_Wm2_coeffs(  1, 14) =  -0.2516E-04
    tot_Poynting_Wm2_coeffs(  1, 15) =  -0.1483E-04
    tot_Poynting_Wm2_coeffs(  1, 16) =   0.2030E-04
    tot_Poynting_Wm2_coeffs(  1, 17) =  -0.2489E-05
    tot_Poynting_Wm2_coeffs(  1, 18) =   0.3592E-04
    tot_Poynting_Wm2_coeffs(  1, 19) =  -0.7059E-05
    tot_Poynting_Wm2_coeffs(  1, 20) =   0.1006E-04
    tot_Poynting_Wm2_coeffs(  1, 21) =   0.1091E-04
    tot_Poynting_Wm2_coeffs(  1, 22) =   0.4157E-10
    tot_Poynting_Wm2_coeffs(  1, 23) =  -0.9216E-05
    tot_Poynting_Wm2_coeffs(  1, 24) =  -0.2785E-10
    tot_Poynting_Wm2_coeffs(  1, 25) =  -0.6192E-05
    tot_Poynting_Wm2_coeffs(  1, 26) =   0.1407E-04
    tot_Poynting_Wm2_coeffs(  1, 27) =   0.2603E-05
    tot_Poynting_Wm2_coeffs(  1, 28) =    0.000    
    tot_Poynting_Wm2_coeffs(  1, 29) =   0.3266E-05
    tot_Poynting_Wm2_coeffs(  1, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  2,  1) =  -0.1336E-05
    tot_Poynting_Wm2_coeffs(  2,  2) =   0.2030E-04
    tot_Poynting_Wm2_coeffs(  2,  3) =  -0.1043E-04
    tot_Poynting_Wm2_coeffs(  2,  4) =  -0.9779E-04
    tot_Poynting_Wm2_coeffs(  2,  5) =  -0.1437E-05
    tot_Poynting_Wm2_coeffs(  2,  6) =  -0.5359E-04
    tot_Poynting_Wm2_coeffs(  2,  7) =  -0.3797E-05
    tot_Poynting_Wm2_coeffs(  2,  8) =   0.2427E-04
    tot_Poynting_Wm2_coeffs(  2,  9) =  -0.6620E-05
    tot_Poynting_Wm2_coeffs(  2, 10) =  -0.5491E-04
    tot_Poynting_Wm2_coeffs(  2, 11) =  -0.1468E-05
    tot_Poynting_Wm2_coeffs(  2, 12) =   0.3596E-05
    tot_Poynting_Wm2_coeffs(  2, 13) =  -0.2799E-05
    tot_Poynting_Wm2_coeffs(  2, 14) =  -0.4908E-04
    tot_Poynting_Wm2_coeffs(  2, 15) =  -0.2088E-04
    tot_Poynting_Wm2_coeffs(  2, 16) =   0.9943E-04
    tot_Poynting_Wm2_coeffs(  2, 17) =  -0.9346E-05
    tot_Poynting_Wm2_coeffs(  2, 18) =   0.7715E-04
    tot_Poynting_Wm2_coeffs(  2, 19) =  -0.1891E-04
    tot_Poynting_Wm2_coeffs(  2, 20) =  -0.1310E-04
    tot_Poynting_Wm2_coeffs(  2, 21) =   0.1888E-04
    tot_Poynting_Wm2_coeffs(  2, 22) =  -0.7993E-09
    tot_Poynting_Wm2_coeffs(  2, 23) =  -0.7001E-05
    tot_Poynting_Wm2_coeffs(  2, 24) =  -0.4618E-10
    tot_Poynting_Wm2_coeffs(  2, 25) =  -0.8218E-05
    tot_Poynting_Wm2_coeffs(  2, 26) =  -0.5404E-05
    tot_Poynting_Wm2_coeffs(  2, 27) =  -0.1452E-05
    tot_Poynting_Wm2_coeffs(  2, 28) =    0.000    
    tot_Poynting_Wm2_coeffs(  2, 29) =  -0.4251E-05
    tot_Poynting_Wm2_coeffs(  2, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  3,  1) =   0.1024E-03
    tot_Poynting_Wm2_coeffs(  3,  2) =  -0.2462E-04
    tot_Poynting_Wm2_coeffs(  3,  3) =  -0.9275E-04
    tot_Poynting_Wm2_coeffs(  3,  4) =  -0.1540E-03
    tot_Poynting_Wm2_coeffs(  3,  5) =  -0.6202E-04
    tot_Poynting_Wm2_coeffs(  3,  6) =  -0.2963E-03
    tot_Poynting_Wm2_coeffs(  3,  7) =   0.1845E-04
    tot_Poynting_Wm2_coeffs(  3,  8) =   0.1402E-04
    tot_Poynting_Wm2_coeffs(  3,  9) =   0.1507E-04
    tot_Poynting_Wm2_coeffs(  3, 10) =  -0.3259E-04
    tot_Poynting_Wm2_coeffs(  3, 11) =   0.2933E-04
    tot_Poynting_Wm2_coeffs(  3, 12) =   0.8053E-05
    tot_Poynting_Wm2_coeffs(  3, 13) =  -0.1211E-03
    tot_Poynting_Wm2_coeffs(  3, 14) =  -0.4869E-04
    tot_Poynting_Wm2_coeffs(  3, 15) =   0.1267E-03
    tot_Poynting_Wm2_coeffs(  3, 16) =   0.3015E-04
    tot_Poynting_Wm2_coeffs(  3, 17) =   0.5451E-04
    tot_Poynting_Wm2_coeffs(  3, 18) =   0.2431E-03
    tot_Poynting_Wm2_coeffs(  3, 19) =  -0.1416E-03
    tot_Poynting_Wm2_coeffs(  3, 20) =  -0.4438E-04
    tot_Poynting_Wm2_coeffs(  3, 21) =   0.1957E-03
    tot_Poynting_Wm2_coeffs(  3, 22) =   0.3570E-08
    tot_Poynting_Wm2_coeffs(  3, 23) =   0.4478E-04
    tot_Poynting_Wm2_coeffs(  3, 24) =  -0.4914E-08
    tot_Poynting_Wm2_coeffs(  3, 25) =  -0.1363E-04
    tot_Poynting_Wm2_coeffs(  3, 26) =   0.1551E-03
    tot_Poynting_Wm2_coeffs(  3, 27) =   0.7785E-04
    tot_Poynting_Wm2_coeffs(  3, 28) =    0.000    
    tot_Poynting_Wm2_coeffs(  3, 29) =  -0.7480E-04
    tot_Poynting_Wm2_coeffs(  3, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  4,  1) =   0.4233E-03
    tot_Poynting_Wm2_coeffs(  4,  2) =  -0.3009E-03
    tot_Poynting_Wm2_coeffs(  4,  3) =  -0.2038E-03
    tot_Poynting_Wm2_coeffs(  4,  4) =  -0.3422E-03
    tot_Poynting_Wm2_coeffs(  4,  5) =  -0.3594E-03
    tot_Poynting_Wm2_coeffs(  4,  6) =  -0.8038E-04
    tot_Poynting_Wm2_coeffs(  4,  7) =   0.2550E-03
    tot_Poynting_Wm2_coeffs(  4,  8) =  -0.1677E-03
    tot_Poynting_Wm2_coeffs(  4,  9) =   0.1193E-03
    tot_Poynting_Wm2_coeffs(  4, 10) =   0.4011E-03
    tot_Poynting_Wm2_coeffs(  4, 11) =   0.3983E-04
    tot_Poynting_Wm2_coeffs(  4, 12) =   0.6429E-03
    tot_Poynting_Wm2_coeffs(  4, 13) =  -0.4788E-03
    tot_Poynting_Wm2_coeffs(  4, 14) =  -0.2204E-03
    tot_Poynting_Wm2_coeffs(  4, 15) =   0.4919E-03
    tot_Poynting_Wm2_coeffs(  4, 16) =  -0.2418E-03
    tot_Poynting_Wm2_coeffs(  4, 17) =   0.5177E-03
    tot_Poynting_Wm2_coeffs(  4, 18) =   0.1151E-03
    tot_Poynting_Wm2_coeffs(  4, 19) =  -0.2677E-03
    tot_Poynting_Wm2_coeffs(  4, 20) =  -0.6374E-03
    tot_Poynting_Wm2_coeffs(  4, 21) =   0.4052E-03
    tot_Poynting_Wm2_coeffs(  4, 22) =  -0.5515E-08
    tot_Poynting_Wm2_coeffs(  4, 23) =   0.4063E-03
    tot_Poynting_Wm2_coeffs(  4, 24) =   0.1499E-07
    tot_Poynting_Wm2_coeffs(  4, 25) =  -0.3037E-04
    tot_Poynting_Wm2_coeffs(  4, 26) =   0.3375E-03
    tot_Poynting_Wm2_coeffs(  4, 27) =   0.5479E-04
    tot_Poynting_Wm2_coeffs(  4, 28) =    0.000    
    tot_Poynting_Wm2_coeffs(  4, 29) =  -0.1692E-03
    tot_Poynting_Wm2_coeffs(  4, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  5,  1) =   0.1839E-02
    tot_Poynting_Wm2_coeffs(  5,  2) =  -0.9261E-03
    tot_Poynting_Wm2_coeffs(  5,  3) =  -0.1138E-02
    tot_Poynting_Wm2_coeffs(  5,  4) =  -0.5788E-03
    tot_Poynting_Wm2_coeffs(  5,  5) =  -0.1525E-03
    tot_Poynting_Wm2_coeffs(  5,  6) =  -0.1369E-02
    tot_Poynting_Wm2_coeffs(  5,  7) =   0.9375E-03
    tot_Poynting_Wm2_coeffs(  5,  8) =  -0.6329E-03
    tot_Poynting_Wm2_coeffs(  5,  9) =   0.1722E-03
    tot_Poynting_Wm2_coeffs(  5, 10) =  -0.7318E-03
    tot_Poynting_Wm2_coeffs(  5, 11) =   0.2136E-05
    tot_Poynting_Wm2_coeffs(  5, 12) =   0.2074E-02
    tot_Poynting_Wm2_coeffs(  5, 13) =   0.2172E-04
    tot_Poynting_Wm2_coeffs(  5, 14) =  -0.4351E-04
    tot_Poynting_Wm2_coeffs(  5, 15) =  -0.9445E-04
    tot_Poynting_Wm2_coeffs(  5, 16) =  -0.6352E-03
    tot_Poynting_Wm2_coeffs(  5, 17) =   0.1127E-02
    tot_Poynting_Wm2_coeffs(  5, 18) =  -0.2458E-02
    tot_Poynting_Wm2_coeffs(  5, 19) =  -0.7241E-03
    tot_Poynting_Wm2_coeffs(  5, 20) =  -0.9206E-03
    tot_Poynting_Wm2_coeffs(  5, 21) =   0.9620E-03
    tot_Poynting_Wm2_coeffs(  5, 22) =   0.4102E-07
    tot_Poynting_Wm2_coeffs(  5, 23) =   0.1353E-03
    tot_Poynting_Wm2_coeffs(  5, 24) =  -0.5773E-08
    tot_Poynting_Wm2_coeffs(  5, 25) =   0.9760E-04
    tot_Poynting_Wm2_coeffs(  5, 26) =   0.1635E-02
    tot_Poynting_Wm2_coeffs(  5, 27) =   0.8871E-04
    tot_Poynting_Wm2_coeffs(  5, 28) =    0.000    
    tot_Poynting_Wm2_coeffs(  5, 29) =  -0.1509E-04
    tot_Poynting_Wm2_coeffs(  5, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  6,  1) =   0.2415E-02
    tot_Poynting_Wm2_coeffs(  6,  2) =   0.1082E-02
    tot_Poynting_Wm2_coeffs(  6,  3) =  -0.1259E-02
    tot_Poynting_Wm2_coeffs(  6,  4) =   0.6978E-03
    tot_Poynting_Wm2_coeffs(  6,  5) =  -0.4079E-03
    tot_Poynting_Wm2_coeffs(  6,  6) =  -0.1435E-02
    tot_Poynting_Wm2_coeffs(  6,  7) =  -0.1132E-02
    tot_Poynting_Wm2_coeffs(  6,  8) =  -0.7440E-03
    tot_Poynting_Wm2_coeffs(  6,  9) =   0.1120E-02
    tot_Poynting_Wm2_coeffs(  6, 10) =  -0.1668E-02
    tot_Poynting_Wm2_coeffs(  6, 11) =  -0.1583E-03
    tot_Poynting_Wm2_coeffs(  6, 12) =   0.7790E-03
    tot_Poynting_Wm2_coeffs(  6, 13) =   0.1890E-03
    tot_Poynting_Wm2_coeffs(  6, 14) =  -0.5096E-03
    tot_Poynting_Wm2_coeffs(  6, 15) =  -0.4216E-03
    tot_Poynting_Wm2_coeffs(  6, 16) =  -0.3060E-03
    tot_Poynting_Wm2_coeffs(  6, 17) =  -0.1046E-04
    tot_Poynting_Wm2_coeffs(  6, 18) =  -0.9948E-03
    tot_Poynting_Wm2_coeffs(  6, 19) =  -0.3219E-03
    tot_Poynting_Wm2_coeffs(  6, 20) =   0.2173E-02
    tot_Poynting_Wm2_coeffs(  6, 21) =  -0.1570E-04
    tot_Poynting_Wm2_coeffs(  6, 22) =   0.1608E-07
    tot_Poynting_Wm2_coeffs(  6, 23) =  -0.1288E-02
    tot_Poynting_Wm2_coeffs(  6, 24) =   0.4515E-07
    tot_Poynting_Wm2_coeffs(  6, 25) =   0.6558E-03
    tot_Poynting_Wm2_coeffs(  6, 26) =   0.1187E-02
    tot_Poynting_Wm2_coeffs(  6, 27) =   0.5446E-04
    tot_Poynting_Wm2_coeffs(  6, 28) =  -0.2585E-14
    tot_Poynting_Wm2_coeffs(  6, 29) =   0.1581E-04
    tot_Poynting_Wm2_coeffs(  6, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  7,  1) =   0.1696E-02
    tot_Poynting_Wm2_coeffs(  7,  2) =   0.2360E-02
    tot_Poynting_Wm2_coeffs(  7,  3) =  -0.1760E-03
    tot_Poynting_Wm2_coeffs(  7,  4) =  -0.7032E-04
    tot_Poynting_Wm2_coeffs(  7,  5) =  -0.1985E-04
    tot_Poynting_Wm2_coeffs(  7,  6) =   0.1169E-03
    tot_Poynting_Wm2_coeffs(  7,  7) =  -0.2012E-02
    tot_Poynting_Wm2_coeffs(  7,  8) =  -0.3126E-02
    tot_Poynting_Wm2_coeffs(  7,  9) =  -0.5887E-04
    tot_Poynting_Wm2_coeffs(  7, 10) =  -0.1309E-03
    tot_Poynting_Wm2_coeffs(  7, 11) =   0.1341E-03
    tot_Poynting_Wm2_coeffs(  7, 12) =   0.6342E-03
    tot_Poynting_Wm2_coeffs(  7, 13) =   0.7921E-04
    tot_Poynting_Wm2_coeffs(  7, 14) =  -0.3157E-03
    tot_Poynting_Wm2_coeffs(  7, 15) =  -0.4271E-03
    tot_Poynting_Wm2_coeffs(  7, 16) =   0.8194E-04
    tot_Poynting_Wm2_coeffs(  7, 17) =   0.8875E-03
    tot_Poynting_Wm2_coeffs(  7, 18) =  -0.9680E-03
    tot_Poynting_Wm2_coeffs(  7, 19) =   0.7075E-03
    tot_Poynting_Wm2_coeffs(  7, 20) =   0.2933E-02
    tot_Poynting_Wm2_coeffs(  7, 21) =   0.7861E-03
    tot_Poynting_Wm2_coeffs(  7, 22) =  -0.1833E-07
    tot_Poynting_Wm2_coeffs(  7, 23) =  -0.1965E-03
    tot_Poynting_Wm2_coeffs(  7, 24) =  -0.6069E-07
    tot_Poynting_Wm2_coeffs(  7, 25) =  -0.3574E-04
    tot_Poynting_Wm2_coeffs(  7, 26) =   0.1650E-02
    tot_Poynting_Wm2_coeffs(  7, 27) =   0.4800E-03
    tot_Poynting_Wm2_coeffs(  7, 28) =    0.000    
    tot_Poynting_Wm2_coeffs(  7, 29) =  -0.4759E-03
    tot_Poynting_Wm2_coeffs(  7, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  8,  1) =   0.7035E-03
    tot_Poynting_Wm2_coeffs(  8,  2) =   0.1485E-02
    tot_Poynting_Wm2_coeffs(  8,  3) =   0.4741E-03
    tot_Poynting_Wm2_coeffs(  8,  4) =   0.7274E-03
    tot_Poynting_Wm2_coeffs(  8,  5) =  -0.5028E-03
    tot_Poynting_Wm2_coeffs(  8,  6) =  -0.2536E-03
    tot_Poynting_Wm2_coeffs(  8,  7) =  -0.8216E-03
    tot_Poynting_Wm2_coeffs(  8,  8) =  -0.2373E-02
    tot_Poynting_Wm2_coeffs(  8,  9) =  -0.9793E-04
    tot_Poynting_Wm2_coeffs(  8, 10) =  -0.9764E-03
    tot_Poynting_Wm2_coeffs(  8, 11) =   0.3587E-03
    tot_Poynting_Wm2_coeffs(  8, 12) =   0.2364E-02
    tot_Poynting_Wm2_coeffs(  8, 13) =  -0.9754E-04
    tot_Poynting_Wm2_coeffs(  8, 14) =   0.6408E-03
    tot_Poynting_Wm2_coeffs(  8, 15) =   0.1797E-03
    tot_Poynting_Wm2_coeffs(  8, 16) =  -0.8915E-03
    tot_Poynting_Wm2_coeffs(  8, 17) =   0.6172E-03
    tot_Poynting_Wm2_coeffs(  8, 18) =  -0.7637E-03
    tot_Poynting_Wm2_coeffs(  8, 19) =   0.8731E-03
    tot_Poynting_Wm2_coeffs(  8, 20) =   0.1567E-02
    tot_Poynting_Wm2_coeffs(  8, 21) =   0.1156E-03
    tot_Poynting_Wm2_coeffs(  8, 22) =  -0.5294E-08
    tot_Poynting_Wm2_coeffs(  8, 23) =  -0.3139E-03
    tot_Poynting_Wm2_coeffs(  8, 24) =  -0.1817E-07
    tot_Poynting_Wm2_coeffs(  8, 25) =  -0.4166E-03
    tot_Poynting_Wm2_coeffs(  8, 26) =   0.1080E-02
    tot_Poynting_Wm2_coeffs(  8, 27) =   0.3522E-03
    tot_Poynting_Wm2_coeffs(  8, 28) =    0.000    
    tot_Poynting_Wm2_coeffs(  8, 29) =   0.3808E-03
    tot_Poynting_Wm2_coeffs(  8, 30) =    0.000    
    tot_Poynting_Wm2_coeffs(  9,  1) =   0.2868E-03
    tot_Poynting_Wm2_coeffs(  9,  2) =   0.2850E-03
    tot_Poynting_Wm2_coeffs(  9,  3) =  -0.4226E-03
    tot_Poynting_Wm2_coeffs(  9,  4) =   0.6719E-03
    tot_Poynting_Wm2_coeffs(  9,  5) =  -0.3540E-04
    tot_Poynting_Wm2_coeffs(  9,  6) =   0.8444E-04
    tot_Poynting_Wm2_coeffs(  9,  7) =  -0.5883E-03
    tot_Poynting_Wm2_coeffs(  9,  8) =   0.1311E-04
    tot_Poynting_Wm2_coeffs(  9,  9) =   0.8136E-03
    tot_Poynting_Wm2_coeffs(  9, 10) =  -0.1420E-02
    tot_Poynting_Wm2_coeffs(  9, 11) =  -0.3649E-03
    tot_Poynting_Wm2_coeffs(  9, 12) =   0.2168E-02
    tot_Poynting_Wm2_coeffs(  9, 13) =  -0.1128E-03
    tot_Poynting_Wm2_coeffs(  9, 14) =   0.7001E-03
    tot_Poynting_Wm2_coeffs(  9, 15) =   0.5761E-03
    tot_Poynting_Wm2_coeffs(  9, 16) =  -0.1290E-02
    tot_Poynting_Wm2_coeffs(  9, 17) =   0.1523E-03
    tot_Poynting_Wm2_coeffs(  9, 18) =  -0.1233E-05
    tot_Poynting_Wm2_coeffs(  9, 19) =   0.2396E-03
    tot_Poynting_Wm2_coeffs(  9, 20) =  -0.2943E-03
    tot_Poynting_Wm2_coeffs(  9, 21) =  -0.3917E-04
    tot_Poynting_Wm2_coeffs(  9, 22) =   0.4956E-09
    tot_Poynting_Wm2_coeffs(  9, 23) =  -0.1469E-03
    tot_Poynting_Wm2_coeffs(  9, 24) =  -0.3284E-09
    tot_Poynting_Wm2_coeffs(  9, 25) =  -0.3693E-03
    tot_Poynting_Wm2_coeffs(  9, 26) =   0.7590E-03
    tot_Poynting_Wm2_coeffs(  9, 27) =   0.1286E-03
    tot_Poynting_Wm2_coeffs(  9, 28) =  -0.2136E-16
    tot_Poynting_Wm2_coeffs(  9, 29) =  -0.5438E-04
    tot_Poynting_Wm2_coeffs(  9, 30) =    0.000    

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

  end subroutine S1_Bmed_fill_coeff_arrays

end module S1_Bmed_coeff_mod

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine S1_Bmed_choose_coeff_array

  use S1_Bmed_coeff_mod, only: tot_Poynting_Wm2_coeffs,&
       & S1_Bmed_fill_coeff_arrays
  use S1_Bmed_model_mod, only: nmlat,nfit,b,b_interp,&
       & b_90,x,imlat,ifit,eps,num_reduc_coef
  integer :: binary_fit_by

  !****************************************************************
  ! Fill parameter arrays:
  !****************************************************************
  call S1_Bmed_fill_coeff_arrays

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
  num_reduc_coef(:,2)  = (/ 0,0,2,21,17,12,12,20,1 /)
  num_reduc_coef(:,3)  = (/ 0,0,7,10,6,9,2,11,9 /)
  num_reduc_coef(:,4)  = (/ 9,18,12,9,6,9,5,8,9 /)
  num_reduc_coef(:,5)  = (/ 0,0,13,10,11,3,0,4,0 /)
  num_reduc_coef(:,6)  = (/ 0,7,13,8,12,0,6,5,0 /)
  num_reduc_coef(:,7)  = (/ 0,0,4,20,18,3,7,17,15 /)
  num_reduc_coef(:,8)  = (/ 9,2,0,7,1,9,9,11,0 /)
  num_reduc_coef(:,9)  = (/ 0,0,0,13,6,9,1,0,12 /)
  num_reduc_coef(:,10) = (/ 9,18,5,12,0,9,5,11,13 /)
  num_reduc_coef(:,11) = (/ 0,0,13,8,0,3,6,2,13 /)
  num_reduc_coef(:,12) = (/ 0,0,0,9,12,3,0,7,13 /)
  num_reduc_coef(:,13) = (/ 0,0,10,13,0,9,6,9,0 /)
  num_reduc_coef(:,14) = (/ 0,7,5,16,4,9,7,14,11 /)
  num_reduc_coef(:,15) = (/ 0,0,7,4,6,3,1,8,20 /)
  num_reduc_coef(:,16) = (/ 0,12,2,10,6,6,2,15,20 /)
  num_reduc_coef(:,17) = (/ 0,0,8,11,1,0,6,7,0 /)
  num_reduc_coef(:,18) = (/ 0,7,13,9,12,3,6,7,0 /)
  num_reduc_coef(:,19) = (/ 0,1,18,19,12,12,6,14,0 /)
  num_reduc_coef(:,20) = (/ 0,0,5,22,11,12,10,14,0 /)
  num_reduc_coef(:,21) = (/ 0,1,12,5,6,0,8,0,0 /)
  num_reduc_coef(:,22) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,23) = (/ 0,0,0,4,1,3,2,7,0 /)
  num_reduc_coef(:,24) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,25) = (/ 0,0,0,4,18,6,0,15,7 /)
  num_reduc_coef(:,26) = (/ 0,0,12,12,17,6,14,19,11 /)
  num_reduc_coef(:,27) = (/ 0,0,12,12,6,3,8,7,0 /)
  num_reduc_coef(:,28) = (/ 0,0,0,0,0,0,0,0,0 /)
  num_reduc_coef(:,29) = (/ 0,0,13,10,0,0,6,8,0 /)
  num_reduc_coef(:,30) = (/ 0,0,0,0,0,0,0,0,0 /)

  !****************************************************************
  ! Transfer values:
  !****************************************************************

  do imlat=1,nmlat
     do ifit=1,nfit
        b(imlat,ifit) = tot_Poynting_Wm2_coeffs(imlat,ifit)* &
             & (1.-eps)**num_reduc_coef(imlat,ifit)
     end do
  end do

end subroutine S1_Bmed_choose_coeff_array

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine S1_Bmed_calc_model_value(mlat_glob,mlt,sinT,imf_angle,&
     &value,value_found)

  ! This version is for public distribution, along with the module
  ! file coeff_mod.f90.
  !
  ! Within the /home/hensel/fortran/fit_data/public subdirectory,
  ! this program is called by call_de2_model.

  use S1_Bmed_coeff_mod
  use S1_Bmed_model_mod

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
        MLT_function=find_MLT_function(iMLT,MLT)

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
                    write(*,'("MLT_function = ",G14.6)') MLT_function
                    write(*,'("Btrans_function = ",G14.6)') 1
                    write(*,'("imf_angle_function = ",G14.6)') &
                         &imf_angle_function
                    write(*,'("sinT_function = ",G14.6)')&
                         & sinT_function
!                    write(*,'("mlat = ",F9.2," MLT = ",F9.2,&
!                         &" Btrans = n/a imf_angle = ",F9.2," sinT = ",&
!                         &F9.2,"," x(",I5,") = ",G13.6)') &
!                         &mlat,MLT,imf_angle,sinT,&
!                         &fit_count,x(fit_count)
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

end subroutine S1_Bmed_calc_model_value

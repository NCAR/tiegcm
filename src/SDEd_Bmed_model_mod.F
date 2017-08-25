       module SDEd_Bmed_model_mod
       
         use common_model_mod, only: pi,fit_by_MLT
       
         implicit none   
       
         logical :: first
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Data type
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         integer :: itype
         integer,parameter :: number_of_datatypes=2
         character (len=80), parameter :: 
     |         data_type_SDEd1='SDEd1',
     |         data_type_SDEd2='SDEd2'
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Calculation conditions:
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         real :: mlat,MLT,Btrans,imf_angle,sinT
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Arrays for fit coefficients:
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Number of fit coefficients for each mlat:
         integer :: nfit
         ! Number of mlat intervals:
         integer :: nmlat
         ! Fit coefficients (dimension nmlat x nfit):
         real, allocatable, dimension(:,:) :: b1,b2 ! nmlat x nfit
         ! Interpolated fit coefficients (dimension nfit):
         real, allocatable, dimension(:) :: b_interp1,b_interp2
         ! Fitting functions (dimension nfit):
         real, allocatable, dimension(:) :: x
         ! Fit coefficients averaged over the pole (dimension nfit):
         real, allocatable, dimension(:) :: b_90 
         ! Number of fit functions used by fit_data:
         integer, parameter :: number_of_MLT_coeffs = 5,
     |         number_of_Btrans_coeffs = 1,
     |         number_of_imf_angle_coeffs = 2,
     |         number_of_sinT_coeffs = 2
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
       
         function find_MLT_function (index,MLT_in) result (MLT_fxn_out)
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
              MLT_fxn_out = sin(pi*MLT_in/6.)
           end select
         end function find_MLT_function
       
         function find_Btrans_function (index,Btrans_in) 
     |                  result (Btrans_fxn_out)
           integer :: index
           real :: Btrans_in, Btrans_fxn_out
       
           select case (index)
           case (1)
              Btrans_fxn_out = 1.
           case (2)
              Btrans_fxn_out =  Btrans_in/5.2 - 1.
           case (3)
              Btrans_fxn_out = (Btrans_in/5.2 - 1.)**2
           end select
         end function find_Btrans_function
       
         function find_imf_angle_function (index,imf_angle_in) 
     |         result (imf_angle_fxn_out)
           integer :: index
           real :: imf_angle_in, imf_angle_fxn_out
       
           select case (index)
           case (1)
              imf_angle_fxn_out = 1.-cos(imf_angle_in)
           case (2)
              imf_angle_fxn_out = sin(imf_angle_in)
           case (3)
              imf_angle_fxn_out = sin(imf_angle_in)
           case (4)
              imf_angle_fxn_out = cos(2.*imf_angle_in)
           case (5)
              imf_angle_fxn_out = sin(2.*imf_angle_in)
           end select
         end function find_imf_angle_function
       
         function find_sinT_function(index,sinT_in) result(sinT_fxn_out)
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
       
       
       end module SDEd_Bmed_model_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       module SDEd_Bmed_coeff_mod
       
         ! This module file created 2007/09/12.
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
         real,  dimension(  9, 20) :: Ed1_Var_coeffs
         real,  dimension(  9, 20) :: Ed2_Var_coeffs
         real,  dimension(  9) :: min_mlat,max_mlat,average_mlat
         integer, parameter :: 
     |          fit_by_Btrans = 0, 
     |          fit_by_imf_angle = 1, 
     |          fit_by_sinT = 1, 
     |          mlat_resolution = 5
       
       contains
       
         subroutine SDEd_Bmed_fill_coeff_arrays
       
           ! Fill Ed1_Var arrays:
           Ed1_Var_coeffs(  1,  1) =    1.591    
           Ed1_Var_coeffs(  1,  2) =    11.50    
           Ed1_Var_coeffs(  1,  3) =   -2.128    
           Ed1_Var_coeffs(  1,  4) =   -11.84    
           Ed1_Var_coeffs(  1,  5) =   0.3683    
           Ed1_Var_coeffs(  1,  6) =    14.70    
           Ed1_Var_coeffs(  1,  7) =   -1.381    
           Ed1_Var_coeffs(  1,  8) =   -9.278    
           Ed1_Var_coeffs(  1,  9) =   -6.051    
           Ed1_Var_coeffs(  1, 10) =   -5.829    
           Ed1_Var_coeffs(  1, 11) =    8.579    
           Ed1_Var_coeffs(  1, 12) =    6.372    
           Ed1_Var_coeffs(  1, 13) =   -3.484    
           Ed1_Var_coeffs(  1, 14) =    5.256    
           Ed1_Var_coeffs(  1, 15) =    5.486    
           Ed1_Var_coeffs(  1, 16) =   -2.711    
           Ed1_Var_coeffs(  1, 17) =   -4.383    
           Ed1_Var_coeffs(  1, 18) =   -5.806    
           Ed1_Var_coeffs(  1, 19) =    7.322    
           Ed1_Var_coeffs(  1, 20) =    5.842    
           Ed1_Var_coeffs(  2,  1) =   -1.836    
           Ed1_Var_coeffs(  2,  2) =    16.61    
           Ed1_Var_coeffs(  2,  3) =   0.1880    
           Ed1_Var_coeffs(  2,  4) =   -3.830    
           Ed1_Var_coeffs(  2,  5) =   -3.384    
           Ed1_Var_coeffs(  2,  6) =    14.74    
           Ed1_Var_coeffs(  2,  7) =    1.926    
           Ed1_Var_coeffs(  2,  8) =   -3.799    
           Ed1_Var_coeffs(  2,  9) =   -3.023    
           Ed1_Var_coeffs(  2, 10) =  -0.7545    
           Ed1_Var_coeffs(  2, 11) =    2.276    
           Ed1_Var_coeffs(  2, 12) =    1.842    
           Ed1_Var_coeffs(  2, 13) =   -3.720    
           Ed1_Var_coeffs(  2, 14) =    12.08    
           Ed1_Var_coeffs(  2, 15) =    3.529    
           Ed1_Var_coeffs(  2, 16) =    1.098    
           Ed1_Var_coeffs(  2, 17) =   -3.592    
           Ed1_Var_coeffs(  2, 18) =    5.990    
           Ed1_Var_coeffs(  2, 19) =   0.1472    
           Ed1_Var_coeffs(  2, 20) =    1.040    
           Ed1_Var_coeffs(  3,  1) =    4.491    
           Ed1_Var_coeffs(  3,  2) =    8.602    
           Ed1_Var_coeffs(  3,  3) =    3.015    
           Ed1_Var_coeffs(  3,  4) =   -14.83    
           Ed1_Var_coeffs(  3,  5) =    6.612    
           Ed1_Var_coeffs(  3,  6) =   -6.123    
           Ed1_Var_coeffs(  3,  7) =    6.242    
           Ed1_Var_coeffs(  3,  8) =   -16.59    
           Ed1_Var_coeffs(  3,  9) =    2.387    
           Ed1_Var_coeffs(  3, 10) =   -11.57    
           Ed1_Var_coeffs(  3, 11) =    4.036    
           Ed1_Var_coeffs(  3, 12) =   -2.517    
           Ed1_Var_coeffs(  3, 13) =    1.741    
           Ed1_Var_coeffs(  3, 14) =    6.788    
           Ed1_Var_coeffs(  3, 15) =    4.285    
           Ed1_Var_coeffs(  3, 16) =   -4.778    
           Ed1_Var_coeffs(  3, 17) =   0.5895    
           Ed1_Var_coeffs(  3, 18) =    3.875    
           Ed1_Var_coeffs(  3, 19) =  -0.2139    
           Ed1_Var_coeffs(  3, 20) =   -6.282    
           Ed1_Var_coeffs(  4,  1) =    13.44    
           Ed1_Var_coeffs(  4,  2) =    3.753    
           Ed1_Var_coeffs(  4,  3) =    5.536    
           Ed1_Var_coeffs(  4,  4) =   -20.12    
           Ed1_Var_coeffs(  4,  5) =    12.45    
           Ed1_Var_coeffs(  4,  6) =   -8.855    
           Ed1_Var_coeffs(  4,  7) =    10.72    
           Ed1_Var_coeffs(  4,  8) =   -17.25    
           Ed1_Var_coeffs(  4,  9) =    1.792    
           Ed1_Var_coeffs(  4, 10) =   -8.481    
           Ed1_Var_coeffs(  4, 11) =    6.705    
           Ed1_Var_coeffs(  4, 12) =   -4.262    
           Ed1_Var_coeffs(  4, 13) =    1.502    
           Ed1_Var_coeffs(  4, 14) =    23.39    
           Ed1_Var_coeffs(  4, 15) =    3.736    
           Ed1_Var_coeffs(  4, 16) =    1.493    
           Ed1_Var_coeffs(  4, 17) =    2.478    
           Ed1_Var_coeffs(  4, 18) =    13.38    
           Ed1_Var_coeffs(  4, 19) =    3.188    
           Ed1_Var_coeffs(  4, 20) =   -2.604    
           Ed1_Var_coeffs(  5,  1) =    43.32    
           Ed1_Var_coeffs(  5,  2) =    21.47    
           Ed1_Var_coeffs(  5,  3) =    11.88    
           Ed1_Var_coeffs(  5,  4) =   -41.24    
           Ed1_Var_coeffs(  5,  5) =    14.41    
           Ed1_Var_coeffs(  5,  6) =    25.94    
           Ed1_Var_coeffs(  5,  7) =    23.38    
           Ed1_Var_coeffs(  5,  8) =   -11.21    
           Ed1_Var_coeffs(  5,  9) =   -4.095    
           Ed1_Var_coeffs(  5, 10) =    31.62    
           Ed1_Var_coeffs(  5, 11) =    15.29    
           Ed1_Var_coeffs(  5, 12) =   -12.76    
           Ed1_Var_coeffs(  5, 13) =   -10.77    
           Ed1_Var_coeffs(  5, 14) =    36.15    
           Ed1_Var_coeffs(  5, 15) =    20.23    
           Ed1_Var_coeffs(  5, 16) =    3.517    
           Ed1_Var_coeffs(  5, 17) =    1.997    
           Ed1_Var_coeffs(  5, 18) =    50.18    
           Ed1_Var_coeffs(  5, 19) =    8.885    
           Ed1_Var_coeffs(  5, 20) =    32.14    
           Ed1_Var_coeffs(  6,  1) =    36.57    
           Ed1_Var_coeffs(  6,  2) =   -23.52    
           Ed1_Var_coeffs(  6,  3) =    17.58    
           Ed1_Var_coeffs(  6,  4) =   -44.30    
           Ed1_Var_coeffs(  6,  5) =   -67.02    
           Ed1_Var_coeffs(  6,  6) =    214.7    
           Ed1_Var_coeffs(  6,  7) =    20.88    
           Ed1_Var_coeffs(  6,  8) =   -52.95    
           Ed1_Var_coeffs(  6,  9) =   -4.402    
           Ed1_Var_coeffs(  6, 10) =   -29.77    
           Ed1_Var_coeffs(  6, 11) =    46.52    
           Ed1_Var_coeffs(  6, 12) =   -68.63    
           Ed1_Var_coeffs(  6, 13) =   -21.11    
           Ed1_Var_coeffs(  6, 14) =   -7.363    
           Ed1_Var_coeffs(  6, 15) =   -6.302    
           Ed1_Var_coeffs(  6, 16) =    42.07    
           Ed1_Var_coeffs(  6, 17) =   -29.24    
           Ed1_Var_coeffs(  6, 18) =    142.7    
           Ed1_Var_coeffs(  6, 19) =    2.770    
           Ed1_Var_coeffs(  6, 20) =    12.63    
           Ed1_Var_coeffs(  7,  1) =   -82.36    
           Ed1_Var_coeffs(  7,  2) =    101.1    
           Ed1_Var_coeffs(  7,  3) =    43.68    
           Ed1_Var_coeffs(  7,  4) =   -153.2    
           Ed1_Var_coeffs(  7,  5) =   -41.35    
           Ed1_Var_coeffs(  7,  6) =    252.8    
           Ed1_Var_coeffs(  7,  7) =    68.18    
           Ed1_Var_coeffs(  7,  8) =   -104.5    
           Ed1_Var_coeffs(  7,  9) =   -12.70    
           Ed1_Var_coeffs(  7, 10) =   -21.29    
           Ed1_Var_coeffs(  7, 11) =    53.00    
           Ed1_Var_coeffs(  7, 12) =   -101.7    
           Ed1_Var_coeffs(  7, 13) =    4.534    
           Ed1_Var_coeffs(  7, 14) =   -136.9    
           Ed1_Var_coeffs(  7, 15) =    3.176    
           Ed1_Var_coeffs(  7, 16) =    9.718    
           Ed1_Var_coeffs(  7, 17) =    5.390    
           Ed1_Var_coeffs(  7, 18) =    88.65    
           Ed1_Var_coeffs(  7, 19) =    10.29    
           Ed1_Var_coeffs(  7, 20) =    84.20    
           Ed1_Var_coeffs(  8,  1) =   -193.0    
           Ed1_Var_coeffs(  8,  2) =    242.3    
           Ed1_Var_coeffs(  8,  3) =    28.93    
           Ed1_Var_coeffs(  8,  4) =   -147.2    
           Ed1_Var_coeffs(  8,  5) =    2.774    
           Ed1_Var_coeffs(  8,  6) =    320.3    
           Ed1_Var_coeffs(  8,  7) =    64.60    
           Ed1_Var_coeffs(  8,  8) =   -57.07    
           Ed1_Var_coeffs(  8,  9) =   -8.900    
           Ed1_Var_coeffs(  8, 10) =   -7.707    
           Ed1_Var_coeffs(  8, 11) =    22.00    
           Ed1_Var_coeffs(  8, 12) =   -46.52    
           Ed1_Var_coeffs(  8, 13) =   -103.2    
           Ed1_Var_coeffs(  8, 14) =    27.19    
           Ed1_Var_coeffs(  8, 15) =  -0.5788    
           Ed1_Var_coeffs(  8, 16) =   -52.95    
           Ed1_Var_coeffs(  8, 17) =    64.38    
           Ed1_Var_coeffs(  8, 18) =   -140.7    
           Ed1_Var_coeffs(  8, 19) =    8.194    
           Ed1_Var_coeffs(  8, 20) =   -19.41    
           Ed1_Var_coeffs(  9,  1) =   -151.8    
           Ed1_Var_coeffs(  9,  2) =    228.2    
           Ed1_Var_coeffs(  9,  3) =   -10.74    
           Ed1_Var_coeffs(  9,  4) =   -71.54    
           Ed1_Var_coeffs(  9,  5) =   -16.44    
           Ed1_Var_coeffs(  9,  6) =    276.1    
           Ed1_Var_coeffs(  9,  7) =   -19.71    
           Ed1_Var_coeffs(  9,  8) =    122.3    
           Ed1_Var_coeffs(  9,  9) =    32.70    
           Ed1_Var_coeffs(  9, 10) =   -140.5    
           Ed1_Var_coeffs(  9, 11) =    22.13    
           Ed1_Var_coeffs(  9, 12) =   -98.66    
           Ed1_Var_coeffs(  9, 13) =   -79.11    
           Ed1_Var_coeffs(  9, 14) =    110.7    
           Ed1_Var_coeffs(  9, 15) =   -22.93    
           Ed1_Var_coeffs(  9, 16) =    51.28    
           Ed1_Var_coeffs(  9, 17) =    18.54    
           Ed1_Var_coeffs(  9, 18) =   -27.67    
           Ed1_Var_coeffs(  9, 19) =   -12.30    
           Ed1_Var_coeffs(  9, 20) =    85.26    
           ! Fill Ed2_Var arrays:
           Ed2_Var_coeffs(  1,  1) =    3.249    
           Ed2_Var_coeffs(  1,  2) =    3.652    
           Ed2_Var_coeffs(  1,  3) =  -0.9378    
           Ed2_Var_coeffs(  1,  4) =   -8.083    
           Ed2_Var_coeffs(  1,  5) =    2.569    
           Ed2_Var_coeffs(  1,  6) =   0.2719    
           Ed2_Var_coeffs(  1,  7) =    1.686    
           Ed2_Var_coeffs(  1,  8) =   -8.226    
           Ed2_Var_coeffs(  1,  9) =   -3.842    
           Ed2_Var_coeffs(  1, 10) =   -5.814    
           Ed2_Var_coeffs(  1, 11) =    5.066    
           Ed2_Var_coeffs(  1, 12) =    5.474    
           Ed2_Var_coeffs(  1, 13) =   0.6582    
           Ed2_Var_coeffs(  1, 14) =   -3.307    
           Ed2_Var_coeffs(  1, 15) =    4.185    
           Ed2_Var_coeffs(  1, 16) =  -0.4353    
           Ed2_Var_coeffs(  1, 17) =   0.8775    
           Ed2_Var_coeffs(  1, 18) =   -5.797    
           Ed2_Var_coeffs(  1, 19) =   0.5882    
           Ed2_Var_coeffs(  1, 20) =    8.022    
           Ed2_Var_coeffs(  2,  1) =    10.92    
           Ed2_Var_coeffs(  2,  2) =   -2.365    
           Ed2_Var_coeffs(  2,  3) =   0.4994    
           Ed2_Var_coeffs(  2,  4) =   -19.05    
           Ed2_Var_coeffs(  2,  5) =    12.49    
           Ed2_Var_coeffs(  2,  6) =   -13.38    
           Ed2_Var_coeffs(  2,  7) =    6.056    
           Ed2_Var_coeffs(  2,  8) =   -20.46    
           Ed2_Var_coeffs(  2,  9) =   -7.005    
           Ed2_Var_coeffs(  2, 10) =   -12.87    
           Ed2_Var_coeffs(  2, 11) =    5.092    
           Ed2_Var_coeffs(  2, 12) =    14.21    
           Ed2_Var_coeffs(  2, 13) =    3.082    
           Ed2_Var_coeffs(  2, 14) =   -16.26    
           Ed2_Var_coeffs(  2, 15) =    6.859    
           Ed2_Var_coeffs(  2, 16) =    2.962    
           Ed2_Var_coeffs(  2, 17) =  -0.6709    
           Ed2_Var_coeffs(  2, 18) =   -6.129    
           Ed2_Var_coeffs(  2, 19) =   -1.450    
           Ed2_Var_coeffs(  2, 20) =    4.670    
           Ed2_Var_coeffs(  3,  1) =    39.57    
           Ed2_Var_coeffs(  3,  2) =   -8.510    
           Ed2_Var_coeffs(  3,  3) =   -3.300    
           Ed2_Var_coeffs(  3,  4) =   -41.03    
           Ed2_Var_coeffs(  3,  5) =    52.09    
           Ed2_Var_coeffs(  3,  6) =   -52.29    
           Ed2_Var_coeffs(  3,  7) =    11.88    
           Ed2_Var_coeffs(  3,  8) =   -3.974    
           Ed2_Var_coeffs(  3,  9) =   -39.01    
           Ed2_Var_coeffs(  3, 10) =   -8.030    
           Ed2_Var_coeffs(  3, 11) =    13.02    
           Ed2_Var_coeffs(  3, 12) =    51.48    
           Ed2_Var_coeffs(  3, 13) =    9.016    
           Ed2_Var_coeffs(  3, 14) =   -26.30    
           Ed2_Var_coeffs(  3, 15) =    16.61    
           Ed2_Var_coeffs(  3, 16) =    25.94    
           Ed2_Var_coeffs(  3, 17) =   -27.98    
           Ed2_Var_coeffs(  3, 18) =    40.27    
           Ed2_Var_coeffs(  3, 19) =   -9.981    
           Ed2_Var_coeffs(  3, 20) =   0.6692    
           Ed2_Var_coeffs(  4,  1) =    32.96    
           Ed2_Var_coeffs(  4,  2) =    51.13    
           Ed2_Var_coeffs(  4,  3) =   -7.628    
           Ed2_Var_coeffs(  4,  4) =   -76.49    
           Ed2_Var_coeffs(  4,  5) =    12.57    
           Ed2_Var_coeffs(  4,  6) =   -45.92    
           Ed2_Var_coeffs(  4,  7) =    16.04    
           Ed2_Var_coeffs(  4,  8) =    42.87    
           Ed2_Var_coeffs(  4,  9) =   -36.66    
           Ed2_Var_coeffs(  4, 10) =   -57.20    
           Ed2_Var_coeffs(  4, 11) =    15.54    
           Ed2_Var_coeffs(  4, 12) =    68.17    
           Ed2_Var_coeffs(  4, 13) =  -0.6152    
           Ed2_Var_coeffs(  4, 14) =   -17.83    
           Ed2_Var_coeffs(  4, 15) =    8.626    
           Ed2_Var_coeffs(  4, 16) =    34.49    
           Ed2_Var_coeffs(  4, 17) =    12.42    
           Ed2_Var_coeffs(  4, 18) =    70.15    
           Ed2_Var_coeffs(  4, 19) =   -20.24    
           Ed2_Var_coeffs(  4, 20) =   -55.51    
           Ed2_Var_coeffs(  5,  1) =    133.5    
           Ed2_Var_coeffs(  5,  2) =    23.04    
           Ed2_Var_coeffs(  5,  3) =    16.81    
           Ed2_Var_coeffs(  5,  4) =   -94.37    
           Ed2_Var_coeffs(  5,  5) =    21.05    
           Ed2_Var_coeffs(  5,  6) =    59.35    
           Ed2_Var_coeffs(  5,  7) =   -5.558    
           Ed2_Var_coeffs(  5,  8) =    65.34    
           Ed2_Var_coeffs(  5,  9) =   -56.63    
           Ed2_Var_coeffs(  5, 10) =    111.2    
           Ed2_Var_coeffs(  5, 11) =    75.48    
           Ed2_Var_coeffs(  5, 12) =   -106.0    
           Ed2_Var_coeffs(  5, 13) =   -62.98    
           Ed2_Var_coeffs(  5, 14) =    134.5    
           Ed2_Var_coeffs(  5, 15) =    14.99    
           Ed2_Var_coeffs(  5, 16) =   -16.57    
           Ed2_Var_coeffs(  5, 17) =   -6.674    
           Ed2_Var_coeffs(  5, 18) =    284.4    
           Ed2_Var_coeffs(  5, 19) =    32.59    
           Ed2_Var_coeffs(  5, 20) =    166.7    
           Ed2_Var_coeffs(  6,  1) =    151.9    
           Ed2_Var_coeffs(  6,  2) =   -55.31    
           Ed2_Var_coeffs(  6,  3) =    35.51    
           Ed2_Var_coeffs(  6,  4) =   -132.1    
           Ed2_Var_coeffs(  6,  5) =   -191.4    
           Ed2_Var_coeffs(  6,  6) =    423.6    
           Ed2_Var_coeffs(  6,  7) =    31.97    
           Ed2_Var_coeffs(  6,  8) =    93.90    
           Ed2_Var_coeffs(  6,  9) =   -1.536    
           Ed2_Var_coeffs(  6, 10) =    163.0    
           Ed2_Var_coeffs(  6, 11) =    120.7    
           Ed2_Var_coeffs(  6, 12) =   -136.3    
           Ed2_Var_coeffs(  6, 13) =   -38.49    
           Ed2_Var_coeffs(  6, 14) =    33.62    
           Ed2_Var_coeffs(  6, 15) =   -25.59    
           Ed2_Var_coeffs(  6, 16) =   -66.71    
           Ed2_Var_coeffs(  6, 17) =   -69.17    
           Ed2_Var_coeffs(  6, 18) =    407.4    
           Ed2_Var_coeffs(  6, 19) =    21.07    
           Ed2_Var_coeffs(  6, 20) =   -1.811    
           Ed2_Var_coeffs(  7,  1) =   -10.73    
           Ed2_Var_coeffs(  7,  2) =    95.25    
           Ed2_Var_coeffs(  7,  3) =    46.65    
           Ed2_Var_coeffs(  7,  4) =    8.193    
           Ed2_Var_coeffs(  7,  5) =   -5.722    
           Ed2_Var_coeffs(  7,  6) =   -25.43    
           Ed2_Var_coeffs(  7,  7) =    42.51    
           Ed2_Var_coeffs(  7,  8) =   -161.3    
           Ed2_Var_coeffs(  7,  9) =    88.77    
           Ed2_Var_coeffs(  7, 10) =   -140.4    
           Ed2_Var_coeffs(  7, 11) =    134.4    
           Ed2_Var_coeffs(  7, 12) =   -363.7    
           Ed2_Var_coeffs(  7, 13) =    43.39    
           Ed2_Var_coeffs(  7, 14) =   -104.4    
           Ed2_Var_coeffs(  7, 15) =   -22.28    
           Ed2_Var_coeffs(  7, 16) =    252.4    
           Ed2_Var_coeffs(  7, 17) =   -33.46    
           Ed2_Var_coeffs(  7, 18) =    80.00    
           Ed2_Var_coeffs(  7, 19) =   -12.03    
           Ed2_Var_coeffs(  7, 20) =    262.1    
           Ed2_Var_coeffs(  8,  1) =   -228.6    
           Ed2_Var_coeffs(  8,  2) =    225.1    
           Ed2_Var_coeffs(  8,  3) =  -0.2225    
           Ed2_Var_coeffs(  8,  4) =   -34.32    
           Ed2_Var_coeffs(  8,  5) =    148.5    
           Ed2_Var_coeffs(  8,  6) =    150.4    
           Ed2_Var_coeffs(  8,  7) =    58.66    
           Ed2_Var_coeffs(  8,  8) =   -2.528    
           Ed2_Var_coeffs(  8,  9) =    2.920    
           Ed2_Var_coeffs(  8, 10) =   -153.0    
           Ed2_Var_coeffs(  8, 11) =   -20.77    
           Ed2_Var_coeffs(  8, 12) =   -95.03    
           Ed2_Var_coeffs(  8, 13) =    26.88    
           Ed2_Var_coeffs(  8, 14) =   -320.9    
           Ed2_Var_coeffs(  8, 15) =   -7.202    
           Ed2_Var_coeffs(  8, 16) =   -11.79    
           Ed2_Var_coeffs(  8, 17) =   -3.436    
           Ed2_Var_coeffs(  8, 18) =    252.9    
           Ed2_Var_coeffs(  8, 19) =    20.40    
           Ed2_Var_coeffs(  8, 20) =    58.63    
           Ed2_Var_coeffs(  9,  1) =   -173.2    
           Ed2_Var_coeffs(  9,  2) =    294.7    
           Ed2_Var_coeffs(  9,  3) =    6.188    
           Ed2_Var_coeffs(  9,  4) =   -96.68    
           Ed2_Var_coeffs(  9,  5) =   -9.407    
           Ed2_Var_coeffs(  9,  6) =    353.1    
           Ed2_Var_coeffs(  9,  7) =    16.67    
           Ed2_Var_coeffs(  9,  8) =    68.71    
           Ed2_Var_coeffs(  9,  9) =    1.861    
           Ed2_Var_coeffs(  9, 10) =   -82.29    
           Ed2_Var_coeffs(  9, 11) =    12.17    
           Ed2_Var_coeffs(  9, 12) =    39.78    
           Ed2_Var_coeffs(  9, 13) =    38.71    
           Ed2_Var_coeffs(  9, 14) =   -360.1    
           Ed2_Var_coeffs(  9, 15) =   -14.24    
           Ed2_Var_coeffs(  9, 16) =    46.16    
           Ed2_Var_coeffs(  9, 17) =   -9.628    
           Ed2_Var_coeffs(  9, 18) =    229.1    
           Ed2_Var_coeffs(  9, 19) =   -14.65    
           Ed2_Var_coeffs(  9, 20) =   -18.89    
       
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
       
         end subroutine SDEd_Bmed_fill_coeff_arrays
       
       end module SDEd_Bmed_coeff_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       subroutine SDEd_Bmed_choose_coeff_array
       
         use SDEd_Bmed_coeff_mod, only: Ed1_Var_coeffs,Ed2_Var_coeffs, 
     |         SDEd_Bmed_fill_coeff_arrays
         use SDEd_Bmed_model_mod, only: first,b1,b2,b_interp1,b_interp2,
     |           x,b_90,nmlat,nfit,
     |          itype
         integer :: binary_fit_by
       
         !****************************************************************
         ! Fill parameter arrays:
         !****************************************************************
         call SDEd_Bmed_fill_coeff_arrays
       
         !****************************************************************
         ! Determine array dimensions.  Arrays for data type 'potential'
         ! are the same size as arrays for other data types.
         !****************************************************************
         nmlat = size(Ed1_Var_coeffs,1)
         nfit  = size(Ed1_Var_coeffs,2)
         !****************************************************************
         ! Allocate arrays into which coefficients will be transferred:
         !****************************************************************
         allocate(b1(nmlat,nfit),b2(nmlat,nfit),b_interp1(nfit),
     |            b_interp2(nfit),x(nfit),b_90(nfit))
         first = .false.
       
         !  endif 
       
         !****************************************************************
         ! Transfer values:
         !****************************************************************
       
         do imlat=1,nmlat
            do ifit=1,nfit
               b1(imlat,ifit) = Ed1_Var_coeffs(imlat,ifit)
               b2(imlat,ifit) = Ed2_Var_coeffs(imlat,ifit)
            end do
         end do
       
       end subroutine SDEd_Bmed_choose_coeff_array
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       subroutine SDEd_Bmed_calc_model_value(mlat_glb,mlt,sinT, 
     |        imf_angle,Btrans,value1,value2,value_found)
       
         ! This version is for public distribution, along with the module
         ! file coeff_mod.f90.
         !
         ! Within the /home/hensel/fortran/fit_data/public subdirectory,
         ! this program is called by call_de2_model.
       
         use SDEd_Bmed_coeff_mod,only: average_mlat, fit_by_Btrans, 
     |          fit_by_imf_angle,fit_by_sinT
         use SDEd_Bmed_model_mod, only: x,b1,b2, b_interp1,b_interp2, 
     |           b_90, imlat,ifit,
     |          nfit,nmlat,
     |          do_smooth,
     |          number_of_MLT_coeffs,
     |          number_of_Btrans_coeffs,number_of_imf_angle_coeffs,
     |          number_of_sinT_coeffs,fit_by_mlt,
     |          slp1,slpn,islpsw,sigma,
     |          find_mlt_function,find_btrans_function,
     |           find_imf_angle_function,
     |          find_sinT_function
       
         implicit none
       
         real,intent(in):: mlt,mlat_glb,sinT, 
     |          imf_angle,Btrans
         real,intent(out):: value1,value2
         logical, intent(out) :: value_found
       
         integer :: fit_count,ierr
         logical :: mlat_found
         integer :: iMLT,iBtrans,iimf_angle,isinT
         real :: MLT_function,Btrans_function,imf_angle_function,
     |                  sinT_function,
     |          imf_angle_in, avg_Btrans_array,mlat
       
         ! General:
         logical, parameter :: debug = .false.
         real, parameter :: beta_inf = 16.
         real, parameter :: pi = 3.1415926535898
         real, external :: curv2 
         real, dimension(nmlat) :: b_secderiv1,b_secderiv2,temp,
     |            b_coeff1,b_coeff2,lat_array
       
         !****************************************************************
         ! Interpolate to get coefficients for specified mlat:
         !****************************************************************
       
         mlat_found = .false.
         value_found = .false.
       
         ! Use absolute value of mlat:
         mlat = abs(mlat_glb)
         lat_array(1:nmlat) = average_mlat(1:nmlat)
       
         interpolate_coefficients:do ifit = 1,nfit
            b_coeff1(1:nmlat) = b1(1:nmlat,ifit)
            b_coeff2(1:nmlat) = b2(1:nmlat,ifit)
       
            call curv1(nmlat,lat_array,b_coeff1,slp1,slpn, 
     |                islpsw,b_secderiv1,temp,sigma,ierr) 
            call curv1(nmlat,lat_array,b_coeff2,slp1,slpn, 
     |                islpsw,b_secderiv2,temp,sigma,ierr) 
            b_interp1(ifit) = curv2(mlat,nmlat,lat_array,
     |             b_coeff1,b_secderiv1,sigma)
            b_interp2(ifit) = curv2(mlat,nmlat,lat_array,
     |             b_coeff2,b_secderiv2,sigma)
       
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
       
            MLT_loop2:do iMLT=1,1+(fit_by_MLT*number_of_MLT_coeffs -
     |             fit_by_MLT)
               MLT_function = find_MLT_function(iMLT,MLT)
       
               Btrans_loop2:do iBtrans=1,1+(fit_by_Btrans
     |               *number_of_Btrans_coeffs - fit_by_Btrans)
                  Btrans_function=find_Btrans_function(iBtrans,Btrans)
       
                  imf_angle_loop2:do iimf_angle=1,1+(fit_by_imf_angle
     |                  *number_of_imf_angle_coeffs -fit_by_imf_angle)
                     imf_angle_function = 
     |                  find_imf_angle_function(iimf_angle,imf_angle)
       
                     sinT_loop2:do isinT=1,1+(fit_by_sinT
     |                     *number_of_sinT_coeffs - fit_by_sinT) 
                        sinT_function = find_sinT_function(isinT,sinT)
       
                        x(fit_count) = MLT_function*Btrans_function
     |                        *imf_angle_function*sinT_function
                        if (debug) then
                           write(*,'("fit_count = ",I5)') fit_count
                           write(*,'("MLT_function = ",G12.6)') 
     |                          MLT_function
                           write(*,'("Btrans_function = ",G12.6)') 
     |                          Btrans_function
                           write(*,'("imf_angle_function = ",G12.6)') 
     |                           imf_angle_function
                           write(*,'("sinT_function = ",G12.6)')
     |                            sinT_function
                           write(*,'("mlat = ",F9.2," MLT = ",F9.2,
     |                          " Btrans = ",F9.2," imf_angle = ",
     |                           F9.2," sinT = ",F9.2,
     |                           " x(",I5,") = ",G13.6)') 
     |                           mlat,MLT,Btrans,imf_angle,sinT,
     |                           fit_count,x(fit_count)
                        end if
       
                        fit_count = fit_count + 1
       
                     end do sinT_loop2
                  end do imf_angle_loop2
               end do Btrans_loop2
            end do MLT_loop2
       
            !****************************************************************
            ! Calculate value with interpolated coefficients:
            !****************************************************************
       
            value1 = 0.
            value2 = 0.
            go_through_coefficients:do ifit=1,nfit
               value1 = value1 + x(ifit)*b_interp1(ifit)
               value2 = value2 + x(ifit)*b_interp2(ifit)
            end do go_through_coefficients
            value_found = .true.
         end if inside_mlat_range
       end subroutine SDEd_Bmed_calc_model_value
                                                                                                                                                                 

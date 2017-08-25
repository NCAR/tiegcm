       module SDEd_B0_model_mod
       
         use common_model_mod, only: pi,fit_by_MLT
       
         implicit none   
       
         logical :: first_B0 
       
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
     |         number_of_Btrans_coeffs = 2,
     |         number_of_imf_angle_coeffs = 1,
     |         number_of_sinT_coeffs = 2
         ! Indices:
         integer :: imlat,ifit
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Data smoothing:
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         logical :: do_smooth
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Results
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         real :: value
         logical :: value_found
       
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
     |                 result (Btrans_fxn_out)
           integer :: index
           real :: Btrans_in, Btrans_fxn_out
           real, parameter :: beta_inf = 16.
       
           select case (index)
           case (1)
              Btrans_fxn_out = 1.
           case (2)
              Btrans_fxn_out= Btrans_in /sqrt(1.+(Btrans_in*Btrans_in)/
     |                 (beta_inf*beta_inf))
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
       
       end module SDEd_B0_model_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       module SDEd_B0_coeff_mod
       
         ! This module file created 2007/09/13.
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
         real,  dimension(  9, 20) :: Ed1_Var_coeffs
         real,  dimension(  9, 20) :: Ed2_Var_coeffs
         real,  dimension(  9) :: min_mlat,max_mlat,average_mlat
         integer, parameter :: 
     |          fit_by_Btrans = 1, 
     |          fit_by_imf_angle = 0, 
     |          fit_by_sinT = 1, 
     |          mlat_resolution = 5
       
       contains
       
         subroutine SDEd_B0_fill_coeff_arrays
       
           ! Fill Ed1_Var arrays:
           Ed1_Var_coeffs(  1,  1) =    4.556    
           Ed1_Var_coeffs(  1,  2) =   -6.636    
           Ed1_Var_coeffs(  1,  3) =   0.5487    
           Ed1_Var_coeffs(  1,  4) =   -1.656    
           Ed1_Var_coeffs(  1,  5) =    3.359    
           Ed1_Var_coeffs(  1,  6) =   -6.430    
           Ed1_Var_coeffs(  1,  7) =   0.5088    
           Ed1_Var_coeffs(  1,  8) =    1.946    
           Ed1_Var_coeffs(  1,  9) =  -0.6138    
           Ed1_Var_coeffs(  1, 10) =   -1.927    
           Ed1_Var_coeffs(  1, 11) =   0.5320    
           Ed1_Var_coeffs(  1, 12) =   0.7314    
           Ed1_Var_coeffs(  1, 13) =    2.411    
           Ed1_Var_coeffs(  1, 14) =   -10.75    
           Ed1_Var_coeffs(  1, 15) =  -0.7118    
           Ed1_Var_coeffs(  1, 16) =    1.981    
           Ed1_Var_coeffs(  1, 17) =   0.2438    
           Ed1_Var_coeffs(  1, 18) =    2.938    
           Ed1_Var_coeffs(  1, 19) =  -0.2688    
           Ed1_Var_coeffs(  1, 20) =    2.453    
           Ed1_Var_coeffs(  2,  1) =    7.727    
           Ed1_Var_coeffs(  2,  2) =   -1.947    
           Ed1_Var_coeffs(  2,  3) =  -0.1159    
           Ed1_Var_coeffs(  2,  4) =   -1.552    
           Ed1_Var_coeffs(  2,  5) =    1.965    
           Ed1_Var_coeffs(  2,  6) =   -5.173    
           Ed1_Var_coeffs(  2,  7) =   0.1275    
           Ed1_Var_coeffs(  2,  8) =    2.129    
           Ed1_Var_coeffs(  2,  9) =   -2.388    
           Ed1_Var_coeffs(  2, 10) =    2.474    
           Ed1_Var_coeffs(  2, 11) =   0.9405    
           Ed1_Var_coeffs(  2, 12) =  -0.1718E-01
           Ed1_Var_coeffs(  2, 13) =   -5.385    
           Ed1_Var_coeffs(  2, 14) =    15.63    
           Ed1_Var_coeffs(  2, 15) =   0.8265    
           Ed1_Var_coeffs(  2, 16) =   -2.928    
           Ed1_Var_coeffs(  2, 17) =    1.880    
           Ed1_Var_coeffs(  2, 18) =    3.100    
           Ed1_Var_coeffs(  2, 19) =  -0.4548    
           Ed1_Var_coeffs(  2, 20) =   0.1263    
           Ed1_Var_coeffs(  3,  1) =    12.39    
           Ed1_Var_coeffs(  3,  2) =   -7.218    
           Ed1_Var_coeffs(  3,  3) =  -0.9739    
           Ed1_Var_coeffs(  3,  4) =  -0.3735    
           Ed1_Var_coeffs(  3,  5) =    3.818    
           Ed1_Var_coeffs(  3,  6) =   0.2962E-01
           Ed1_Var_coeffs(  3,  7) =  -0.4840    
           Ed1_Var_coeffs(  3,  8) =   -1.853    
           Ed1_Var_coeffs(  3,  9) =  -0.1573    
           Ed1_Var_coeffs(  3, 10) =    2.710    
           Ed1_Var_coeffs(  3, 11) =  -0.2624    
           Ed1_Var_coeffs(  3, 12) =   0.4606    
           Ed1_Var_coeffs(  3, 13) =   -9.002    
           Ed1_Var_coeffs(  3, 14) =    33.89    
           Ed1_Var_coeffs(  3, 15) =    2.116    
           Ed1_Var_coeffs(  3, 16) =   -8.105    
           Ed1_Var_coeffs(  3, 17) =    5.288    
           Ed1_Var_coeffs(  3, 18) =   -6.815    
           Ed1_Var_coeffs(  3, 19) =   -1.772    
           Ed1_Var_coeffs(  3, 20) =    3.026    
           Ed1_Var_coeffs(  4,  1) =    14.00    
           Ed1_Var_coeffs(  4,  2) =   -5.407    
           Ed1_Var_coeffs(  4,  3) =  -0.6346    
           Ed1_Var_coeffs(  4,  4) =   -2.346    
           Ed1_Var_coeffs(  4,  5) =    3.440    
           Ed1_Var_coeffs(  4,  6) =    12.18    
           Ed1_Var_coeffs(  4,  7) =   0.1217    
           Ed1_Var_coeffs(  4,  8) =   -5.489    
           Ed1_Var_coeffs(  4,  9) =   -2.258    
           Ed1_Var_coeffs(  4, 10) =    3.520    
           Ed1_Var_coeffs(  4, 11) =  -0.5355    
           Ed1_Var_coeffs(  4, 12) =   -1.015    
           Ed1_Var_coeffs(  4, 13) =   -8.185    
           Ed1_Var_coeffs(  4, 14) =    25.98    
           Ed1_Var_coeffs(  4, 15) =    1.339    
           Ed1_Var_coeffs(  4, 16) =   -9.912    
           Ed1_Var_coeffs(  4, 17) =    3.470    
           Ed1_Var_coeffs(  4, 18) =  -0.8594    
           Ed1_Var_coeffs(  4, 19) =   -2.728    
           Ed1_Var_coeffs(  4, 20) =    2.474    
           Ed1_Var_coeffs(  5,  1) =    60.78    
           Ed1_Var_coeffs(  5,  2) =   -102.2    
           Ed1_Var_coeffs(  5,  3) =   -8.991    
           Ed1_Var_coeffs(  5,  4) =    17.02    
           Ed1_Var_coeffs(  5,  5) =    35.51    
           Ed1_Var_coeffs(  5,  6) =   -39.78    
           Ed1_Var_coeffs(  5,  7) =   -1.662    
           Ed1_Var_coeffs(  5,  8) =   -1.980    
           Ed1_Var_coeffs(  5,  9) =   -3.092    
           Ed1_Var_coeffs(  5, 10) =    28.24    
           Ed1_Var_coeffs(  5, 11) =  -0.8015    
           Ed1_Var_coeffs(  5, 12) =    7.101    
           Ed1_Var_coeffs(  5, 13) =   -14.90    
           Ed1_Var_coeffs(  5, 14) =    24.47    
           Ed1_Var_coeffs(  5, 15) =    6.642    
           Ed1_Var_coeffs(  5, 16) =   -24.60    
           Ed1_Var_coeffs(  5, 17) =    9.544    
           Ed1_Var_coeffs(  5, 18) =   -24.25    
           Ed1_Var_coeffs(  5, 19) =   -8.637    
           Ed1_Var_coeffs(  5, 20) =    8.774    
           Ed1_Var_coeffs(  6,  1) =    105.6    
           Ed1_Var_coeffs(  6,  2) =   -144.3    
           Ed1_Var_coeffs(  6,  3) =    17.22    
           Ed1_Var_coeffs(  6,  4) =   -35.96    
           Ed1_Var_coeffs(  6,  5) =    38.89    
           Ed1_Var_coeffs(  6,  6) =   -39.16    
           Ed1_Var_coeffs(  6,  7) =    25.77    
           Ed1_Var_coeffs(  6,  8) =   -63.85    
           Ed1_Var_coeffs(  6,  9) =   -18.98    
           Ed1_Var_coeffs(  6, 10) =    83.52    
           Ed1_Var_coeffs(  6, 11) =    19.54    
           Ed1_Var_coeffs(  6, 12) =   -18.38    
           Ed1_Var_coeffs(  6, 13) =   -10.38    
           Ed1_Var_coeffs(  6, 14) =   -29.68    
           Ed1_Var_coeffs(  6, 15) =    4.386    
           Ed1_Var_coeffs(  6, 16) =   -14.87    
           Ed1_Var_coeffs(  6, 17) =   -12.72    
           Ed1_Var_coeffs(  6, 18) =    17.74    
           Ed1_Var_coeffs(  6, 19) =    16.30    
           Ed1_Var_coeffs(  6, 20) =   -50.11    
           Ed1_Var_coeffs(  7,  1) =    165.3    
           Ed1_Var_coeffs(  7,  2) =   -226.8    
           Ed1_Var_coeffs(  7,  3) =    31.65    
           Ed1_Var_coeffs(  7,  4) =   -58.28    
           Ed1_Var_coeffs(  7,  5) =    28.88    
           Ed1_Var_coeffs(  7,  6) =   -132.8    
           Ed1_Var_coeffs(  7,  7) =    43.94    
           Ed1_Var_coeffs(  7,  8) =   -131.4    
           Ed1_Var_coeffs(  7,  9) =    52.40    
           Ed1_Var_coeffs(  7, 10) =   -79.78    
           Ed1_Var_coeffs(  7, 11) =   -15.05    
           Ed1_Var_coeffs(  7, 12) =    60.02    
           Ed1_Var_coeffs(  7, 13) =    5.639    
           Ed1_Var_coeffs(  7, 14) =    124.8    
           Ed1_Var_coeffs(  7, 15) =    20.17    
           Ed1_Var_coeffs(  7, 16) =   -75.77    
           Ed1_Var_coeffs(  7, 17) =    12.05    
           Ed1_Var_coeffs(  7, 18) =   -4.663    
           Ed1_Var_coeffs(  7, 19) =   -26.07    
           Ed1_Var_coeffs(  7, 20) =    6.933    
           Ed1_Var_coeffs(  8,  1) =    148.9    
           Ed1_Var_coeffs(  8,  2) =   -205.1    
           Ed1_Var_coeffs(  8,  3) =    57.61    
           Ed1_Var_coeffs(  8,  4) =   -28.97    
           Ed1_Var_coeffs(  8,  5) =   -26.16    
           Ed1_Var_coeffs(  8,  6) =   -139.5    
           Ed1_Var_coeffs(  8,  7) =    23.41    
           Ed1_Var_coeffs(  8,  8) =   -202.4    
           Ed1_Var_coeffs(  8,  9) =    58.54    
           Ed1_Var_coeffs(  8, 10) =   -177.8    
           Ed1_Var_coeffs(  8, 11) =   -7.349    
           Ed1_Var_coeffs(  8, 12) =    71.67    
           Ed1_Var_coeffs(  8, 13) =    44.44    
           Ed1_Var_coeffs(  8, 14) =    173.3    
           Ed1_Var_coeffs(  8, 15) =    48.46    
           Ed1_Var_coeffs(  8, 16) =   -82.25    
           Ed1_Var_coeffs(  8, 17) =  -0.4669    
           Ed1_Var_coeffs(  8, 18) =    54.04    
           Ed1_Var_coeffs(  8, 19) =   -38.84    
           Ed1_Var_coeffs(  8, 20) =   -6.120    
           Ed1_Var_coeffs(  9,  1) =    107.1    
           Ed1_Var_coeffs(  9,  2) =   -198.2    
           Ed1_Var_coeffs(  9,  3) =    38.31    
           Ed1_Var_coeffs(  9,  4) =    23.19    
           Ed1_Var_coeffs(  9,  5) =   -51.05    
           Ed1_Var_coeffs(  9,  6) =    116.0    
           Ed1_Var_coeffs(  9,  7) =    51.16    
           Ed1_Var_coeffs(  9,  8) =   -205.8    
           Ed1_Var_coeffs(  9,  9) =   -7.560    
           Ed1_Var_coeffs(  9, 10) =    47.42    
           Ed1_Var_coeffs(  9, 11) =   -13.88    
           Ed1_Var_coeffs(  9, 12) =    70.94    
           Ed1_Var_coeffs(  9, 13) =    15.54    
           Ed1_Var_coeffs(  9, 14) =   -61.12    
           Ed1_Var_coeffs(  9, 15) =    31.05    
           Ed1_Var_coeffs(  9, 16) =   -131.2    
           Ed1_Var_coeffs(  9, 17) =   -19.25    
           Ed1_Var_coeffs(  9, 18) =    48.44    
           Ed1_Var_coeffs(  9, 19) =   -9.100    
           Ed1_Var_coeffs(  9, 20) =    39.75    
           ! Fill Ed2_Var arrays:
           Ed2_Var_coeffs(  1,  1) =    4.617    
           Ed2_Var_coeffs(  1,  2) =   -1.300    
           Ed2_Var_coeffs(  1,  3) =   0.2670    
           Ed2_Var_coeffs(  1,  4) =   -1.146    
           Ed2_Var_coeffs(  1,  5) =    2.263    
           Ed2_Var_coeffs(  1,  6) =   -1.909    
           Ed2_Var_coeffs(  1,  7) =   0.2001E-01
           Ed2_Var_coeffs(  1,  8) =    1.671    
           Ed2_Var_coeffs(  1,  9) =   -2.121    
           Ed2_Var_coeffs(  1, 10) =    3.774    
           Ed2_Var_coeffs(  1, 11) =   0.6875    
           Ed2_Var_coeffs(  1, 12) =   -2.424    
           Ed2_Var_coeffs(  1, 13) =   -2.391    
           Ed2_Var_coeffs(  1, 14) =   0.8198E-01
           Ed2_Var_coeffs(  1, 15) =   0.4508    
           Ed2_Var_coeffs(  1, 16) =   -2.020    
           Ed2_Var_coeffs(  1, 17) =  -0.9139    
           Ed2_Var_coeffs(  1, 18) =  -0.4857    
           Ed2_Var_coeffs(  1, 19) =  -0.4147    
           Ed2_Var_coeffs(  1, 20) =    5.002    
           Ed2_Var_coeffs(  2,  1) =    6.495    
           Ed2_Var_coeffs(  2,  2) =   -3.504    
           Ed2_Var_coeffs(  2,  3) =   0.5132    
           Ed2_Var_coeffs(  2,  4) =   0.1498E-02
           Ed2_Var_coeffs(  2,  5) =    2.004    
           Ed2_Var_coeffs(  2,  6) =   -2.490    
           Ed2_Var_coeffs(  2,  7) =  -0.1507    
           Ed2_Var_coeffs(  2,  8) =   -2.785    
           Ed2_Var_coeffs(  2,  9) =   -6.443    
           Ed2_Var_coeffs(  2, 10) =    4.386    
           Ed2_Var_coeffs(  2, 11) =   0.5665    
           Ed2_Var_coeffs(  2, 12) =   -3.736    
           Ed2_Var_coeffs(  2, 13) =   -4.960    
           Ed2_Var_coeffs(  2, 14) =    4.424    
           Ed2_Var_coeffs(  2, 15) =   0.2701    
           Ed2_Var_coeffs(  2, 16) =   -2.708    
           Ed2_Var_coeffs(  2, 17) =   -1.771    
           Ed2_Var_coeffs(  2, 18) =   0.2872    
           Ed2_Var_coeffs(  2, 19) =   0.6943    
           Ed2_Var_coeffs(  2, 20) =    2.704    
           Ed2_Var_coeffs(  3,  1) =    16.53    
           Ed2_Var_coeffs(  3,  2) =   -3.352    
           Ed2_Var_coeffs(  3,  3) =  -0.7104    
           Ed2_Var_coeffs(  3,  4) =   -7.051    
           Ed2_Var_coeffs(  3,  5) =    9.858    
           Ed2_Var_coeffs(  3,  6) =   -2.970    
           Ed2_Var_coeffs(  3,  7) =   -2.229    
           Ed2_Var_coeffs(  3,  8) =   -6.522    
           Ed2_Var_coeffs(  3,  9) =   -13.48    
           Ed2_Var_coeffs(  3, 10) =   -14.64    
           Ed2_Var_coeffs(  3, 11) =  -0.4067    
           Ed2_Var_coeffs(  3, 12) =    4.832    
           Ed2_Var_coeffs(  3, 13) =   -14.65    
           Ed2_Var_coeffs(  3, 14) =  -0.3277    
           Ed2_Var_coeffs(  3, 15) =   -1.913    
           Ed2_Var_coeffs(  3, 16) =   -7.790    
           Ed2_Var_coeffs(  3, 17) =   -3.169    
           Ed2_Var_coeffs(  3, 18) =   -8.394    
           Ed2_Var_coeffs(  3, 19) =  -0.9566    
           Ed2_Var_coeffs(  3, 20) =    4.816    
           Ed2_Var_coeffs(  4,  1) =    48.49    
           Ed2_Var_coeffs(  4,  2) =    34.44    
           Ed2_Var_coeffs(  4,  3) =  -0.7496    
           Ed2_Var_coeffs(  4,  4) =   -17.40    
           Ed2_Var_coeffs(  4,  5) =    32.31    
           Ed2_Var_coeffs(  4,  6) =    71.59    
           Ed2_Var_coeffs(  4,  7) =    4.260    
           Ed2_Var_coeffs(  4,  8) =   -24.05    
           Ed2_Var_coeffs(  4,  9) =   -39.36    
           Ed2_Var_coeffs(  4, 10) =    29.62    
           Ed2_Var_coeffs(  4, 11) =    3.143    
           Ed2_Var_coeffs(  4, 12) =    2.744    
           Ed2_Var_coeffs(  4, 13) =   -32.85    
           Ed2_Var_coeffs(  4, 14) =    74.87    
           Ed2_Var_coeffs(  4, 15) =    6.246    
           Ed2_Var_coeffs(  4, 16) =   -30.15    
           Ed2_Var_coeffs(  4, 17) =   -13.43    
           Ed2_Var_coeffs(  4, 18) =    44.21    
           Ed2_Var_coeffs(  4, 19) =   -4.612    
           Ed2_Var_coeffs(  4, 20) =   0.9154    
           Ed2_Var_coeffs(  5,  1) =    259.8    
           Ed2_Var_coeffs(  5,  2) =   -482.8    
           Ed2_Var_coeffs(  5,  3) =   -33.98    
           Ed2_Var_coeffs(  5,  4) =    68.22    
           Ed2_Var_coeffs(  5,  5) =    174.8    
           Ed2_Var_coeffs(  5,  6) =   -235.2    
           Ed2_Var_coeffs(  5,  7) =   -9.409    
           Ed2_Var_coeffs(  5,  8) =   -7.739    
           Ed2_Var_coeffs(  5,  9) =   -102.1    
           Ed2_Var_coeffs(  5, 10) =    207.1    
           Ed2_Var_coeffs(  5, 11) =    27.46    
           Ed2_Var_coeffs(  5, 12) =   -54.98    
           Ed2_Var_coeffs(  5, 13) =   -1.156    
           Ed2_Var_coeffs(  5, 14) =   -397.5    
           Ed2_Var_coeffs(  5, 15) =    7.792    
           Ed2_Var_coeffs(  5, 16) =    42.79    
           Ed2_Var_coeffs(  5, 17) =    31.70    
           Ed2_Var_coeffs(  5, 18) =   -41.35    
           Ed2_Var_coeffs(  5, 19) =   -11.93    
           Ed2_Var_coeffs(  5, 20) =    2.204    
           Ed2_Var_coeffs(  6,  1) =    351.4    
           Ed2_Var_coeffs(  6,  2) =   -599.9    
           Ed2_Var_coeffs(  6,  3) =    43.40    
           Ed2_Var_coeffs(  6,  4) =   -173.3    
           Ed2_Var_coeffs(  6,  5) =    55.98    
           Ed2_Var_coeffs(  6,  6) =    19.83    
           Ed2_Var_coeffs(  6,  7) =    109.0    
           Ed2_Var_coeffs(  6,  8) =   -139.2    
           Ed2_Var_coeffs(  6,  9) =   -96.04    
           Ed2_Var_coeffs(  6, 10) =    143.1    
           Ed2_Var_coeffs(  6, 11) =    88.40    
           Ed2_Var_coeffs(  6, 12) =   -154.7    
           Ed2_Var_coeffs(  6, 13) =   -110.2    
           Ed2_Var_coeffs(  6, 14) =    147.6    
           Ed2_Var_coeffs(  6, 15) =    31.61    
           Ed2_Var_coeffs(  6, 16) =   -253.1    
           Ed2_Var_coeffs(  6, 17) =   -85.39    
           Ed2_Var_coeffs(  6, 18) =    466.6    
           Ed2_Var_coeffs(  6, 19) =    123.6    
           Ed2_Var_coeffs(  6, 20) =   -387.0    
           Ed2_Var_coeffs(  7,  1) =    340.2    
           Ed2_Var_coeffs(  7,  2) =   -483.1    
           Ed2_Var_coeffs(  7,  3) =    23.41    
           Ed2_Var_coeffs(  7,  4) =   -41.76    
           Ed2_Var_coeffs(  7,  5) =   -183.0    
           Ed2_Var_coeffs(  7,  6) =    207.0    
           Ed2_Var_coeffs(  7,  7) =    3.173    
           Ed2_Var_coeffs(  7,  8) =   -30.19    
           Ed2_Var_coeffs(  7,  9) =    121.0    
           Ed2_Var_coeffs(  7, 10) =   -253.8    
           Ed2_Var_coeffs(  7, 11) =   -21.83    
           Ed2_Var_coeffs(  7, 12) =    47.32    
           Ed2_Var_coeffs(  7, 13) =   -81.36    
           Ed2_Var_coeffs(  7, 14) =    364.4    
           Ed2_Var_coeffs(  7, 15) =    8.569    
           Ed2_Var_coeffs(  7, 16) =   -91.97    
           Ed2_Var_coeffs(  7, 17) =   -76.95    
           Ed2_Var_coeffs(  7, 18) =    316.8    
           Ed2_Var_coeffs(  7, 19) =    30.85    
           Ed2_Var_coeffs(  7, 20) =   -167.7    
           Ed2_Var_coeffs(  8,  1) =    226.5    
           Ed2_Var_coeffs(  8,  2) =   -178.5    
           Ed2_Var_coeffs(  8,  3) =    51.09    
           Ed2_Var_coeffs(  8,  4) =   -48.03    
           Ed2_Var_coeffs(  8,  5) =   -141.4    
           Ed2_Var_coeffs(  8,  6) =   -185.5    
           Ed2_Var_coeffs(  8,  7) =   -25.75    
           Ed2_Var_coeffs(  8,  8) =   -79.50    
           Ed2_Var_coeffs(  8,  9) =    66.36    
           Ed2_Var_coeffs(  8, 10) =   -2.212    
           Ed2_Var_coeffs(  8, 11) =    52.34    
           Ed2_Var_coeffs(  8, 12) =   -202.1    
           Ed2_Var_coeffs(  8, 13) =   -45.46    
           Ed2_Var_coeffs(  8, 14) =    531.4    
           Ed2_Var_coeffs(  8, 15) =   -12.10    
           Ed2_Var_coeffs(  8, 16) =    84.68    
           Ed2_Var_coeffs(  8, 17) =   -51.99    
           Ed2_Var_coeffs(  8, 18) =   -85.25    
           Ed2_Var_coeffs(  8, 19) =   0.8484    
           Ed2_Var_coeffs(  8, 20) =   -115.6    
           Ed2_Var_coeffs(  9,  1) =    121.0    
           Ed2_Var_coeffs(  9,  2) =   -262.8    
           Ed2_Var_coeffs(  9,  3) =    24.86    
           Ed2_Var_coeffs(  9,  4) =    33.26    
           Ed2_Var_coeffs(  9,  5) =   -36.10    
           Ed2_Var_coeffs(  9,  6) =    22.84    
           Ed2_Var_coeffs(  9,  7) =    26.52    
           Ed2_Var_coeffs(  9,  8) =   -74.75    
           Ed2_Var_coeffs(  9,  9) =   -54.96    
           Ed2_Var_coeffs(  9, 10) =    198.9    
           Ed2_Var_coeffs(  9, 11) =    16.31    
           Ed2_Var_coeffs(  9, 12) =   -86.26    
           Ed2_Var_coeffs(  9, 13) =   -122.6    
           Ed2_Var_coeffs(  9, 14) =    269.7    
           Ed2_Var_coeffs(  9, 15) =    28.05    
           Ed2_Var_coeffs(  9, 16) =    19.44    
           Ed2_Var_coeffs(  9, 17) =    20.48    
           Ed2_Var_coeffs(  9, 18) =    42.77    
           Ed2_Var_coeffs(  9, 19) =    29.52    
           Ed2_Var_coeffs(  9, 20) =   -49.84    
       
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
       
         end subroutine SDEd_B0_fill_coeff_arrays
       
       end module SDEd_B0_coeff_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       subroutine SDEd_B0_choose_coeff_array
       
         use SDEd_B0_coeff_mod, only: Ed1_Var_coeffs,Ed2_Var_coeffs,
     |         SDEd_B0_fill_coeff_arrays
         use SDEd_B0_model_mod, only: first_B0,b1,b2,b_interp1,
     |         b_interp2,x,b_90,nmlat,nfit
       
         integer :: binary_fit_by
       
         !****************************************************************
         ! Fill parameter arrays:
         !****************************************************************
         !  if(first_B0) then
         call SDEd_B0_fill_coeff_arrays
       
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
     |          b_interp2(nfit),
     |         x(nfit),b_90(nfit))
         first_B0 = .false.
       
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
       
       end subroutine SDEd_B0_choose_coeff_array
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       subroutine SDEd_B0_calc_model_value(mlat_glb,mlt,sinT, 
     |        imf_angle,Btrans,value1,value2,value_found)
       
         ! This version is for public distribution, along with the module
         ! file coeff_mod.f90.
         !
         ! Within the /home/hensel/fortran/fit_data/public subdirectory,
         ! this program is called by call_de2_model.
       
         use SDEd_B0_coeff_mod, only: average_mlat, 
     |           fit_by_Btrans,fit_by_imf_angle, 
     |          fit_by_sinT
         use SDEd_B0_model_mod, only: x,b1,b2, b_interp1,
     |           b_interp2, b_90, imlat,ifit,
     |          nmlat, nfit, do_smooth, number_of_MLT_coeffs,
     |          number_of_Btrans_coeffs ,number_of_imf_angle_coeffs ,
     |          number_of_sinT_coeffs, fit_by_mlt,
     |          find_MLT_function, find_Btrans_function,
     |           find_imf_angle_function,
     |          find_sinT_function
         use SDEd_Bmed_model_mod, only:slp1,slpn,islpsw,sigma
       
         implicit none
       
         real,intent(in):: mlt,mlat_glb,sinT, 
     |          imf_angle,Btrans
         real,intent(out):: value1,value2
         logical, intent(out) :: value_found
       
         integer :: fit_count,ierr
         logical :: mlat_found
         integer :: iMLT,iBtrans,iimf_angle,isinT
         real :: MLT_function,Btrans_function,imf_angle_function,
     |           sinT_function,
     |          imf_angle_in, avg_Btrans_array,mlat
       
         ! General:
         logical, parameter :: debug = .false.
         real, parameter :: pi = 3.1415926535898
         real, parameter :: beta_inf = 16.
         real, external :: curv2 
         real, dimension(nmlat) :: b_secderiv1,b_secderiv2,temp,
     |           b_coeff1,b_coeff2,lat_array
       
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
                  Btrans_function = find_Btrans_function(iBtrans,Btrans)
       
                  imf_angle_loop2:do iimf_angle=1,1+(fit_by_imf_angle
     |                  *number_of_imf_angle_coeffs - fit_by_imf_angle)
                     imf_angle_function = 
     |                   find_imf_angle_function(iimf_angle,imf_angle)
       
                     sinT_loop2:do isinT=1,1+(fit_by_sinT
     |                     *number_of_sinT_coeffs - fit_by_sinT) 
                        sinT_function = find_sinT_function(isinT,sinT)
       
                        x(fit_count) = MLT_function*Btrans_function
     |                        *imf_angle_function*sinT_function
       
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
       
       end subroutine SDEd_B0_calc_model_value
                                                                                                                                            

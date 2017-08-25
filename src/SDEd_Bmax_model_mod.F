      module SDEd_Bmax_model_mod
      
        use common_model_mod, only: pi,fit_by_MLT
      
        implicit none   
      
        logical :: first_Bmax 
      
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
     |        number_of_Btrans_coeffs = 1,
     |        number_of_imf_angle_coeffs = 3,
     |        number_of_sinT_coeffs = 1
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
     |              result (Btrans_fxn_out)
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
     |        result (imf_angle_fxn_out)
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
      
        function find_sinT_function (index,sinT_in) result(sinT_fxn_out)
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
      
      end module SDEd_Bmax_model_mod
      
      module SDEd_Bmax_coeff_mod
      
        ! This module file created 2007/08/22.
        ! with the following parameters:
        !   for 81252-83047 period 
        !   fit_by_Btrans = 0
        !   fit_by_imf_angle = 1
        !   fit_by_sinT = 0
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
        real,  dimension(  9, 15) :: Ed1_Var_coeffs
        real,  dimension(  9, 15) :: Ed2_Var_coeffs
        real,  dimension(  9, 15) :: Ed1_Var_DE2_coeffs
        real,  dimension(  9, 15) :: Ed2_Var_DE2_coeffs
        real,  dimension(  9, 15) :: Ed1_Var_mod_coeffs
        real,  dimension(  9, 15) :: Ed2_Var_mod_coeffs
        real,  dimension(  9) :: min_mlat,max_mlat,average_mlat
        integer, parameter :: 
     |         fit_by_Btrans = 0, 
     |         fit_by_imf_angle = 1, 
     |         fit_by_sinT = 0, 
     |         mlat_resolution = 5
      
      contains
      
        subroutine SDEd_Bmax_fill_coeff_arrays
      
          ! Fill Ed1_Var arrays:
          Ed1_Var_coeffs(  1,  1) =    6.786    
          Ed1_Var_coeffs(  1,  2) =   -13.39    
          Ed1_Var_coeffs(  1,  3) =   -3.848    
          Ed1_Var_coeffs(  1,  4) =    2.080    
          Ed1_Var_coeffs(  1,  5) =   -5.577    
          Ed1_Var_coeffs(  1,  6) =   -2.810    
          Ed1_Var_coeffs(  1,  7) =   -11.82    
          Ed1_Var_coeffs(  1,  8) =   -18.53    
          Ed1_Var_coeffs(  1,  9) =    1.781    
          Ed1_Var_coeffs(  1, 10) =   -2.561    
          Ed1_Var_coeffs(  1, 11) =   -16.82    
          Ed1_Var_coeffs(  1, 12) =   -5.613    
          Ed1_Var_coeffs(  1, 13) =   -1.668    
          Ed1_Var_coeffs(  1, 14) =    12.12    
          Ed1_Var_coeffs(  1, 15) =    4.161    
          Ed1_Var_coeffs(  2,  1) =    48.44    
          Ed1_Var_coeffs(  2,  2) =    25.97    
          Ed1_Var_coeffs(  2,  3) =    4.630    
          Ed1_Var_coeffs(  2,  4) =    39.93    
          Ed1_Var_coeffs(  2,  5) =    17.53    
          Ed1_Var_coeffs(  2,  6) =   -3.604    
          Ed1_Var_coeffs(  2,  7) =   -57.25    
          Ed1_Var_coeffs(  2,  8) =   -53.19    
          Ed1_Var_coeffs(  2,  9) =   -37.01    
          Ed1_Var_coeffs(  2, 10) =   -24.95    
          Ed1_Var_coeffs(  2, 11) =   -39.40    
          Ed1_Var_coeffs(  2, 12) =   -40.13    
          Ed1_Var_coeffs(  2, 13) =   -49.83    
          Ed1_Var_coeffs(  2, 14) =   -24.48    
          Ed1_Var_coeffs(  2, 15) =   -38.36    
          Ed1_Var_coeffs(  3,  1) =    46.03    
          Ed1_Var_coeffs(  3,  2) =   -9.839    
          Ed1_Var_coeffs(  3,  3) =   -2.911    
          Ed1_Var_coeffs(  3,  4) =    45.06    
          Ed1_Var_coeffs(  3,  5) =   -3.433    
          Ed1_Var_coeffs(  3,  6) =    9.275    
          Ed1_Var_coeffs(  3,  7) =   -6.640    
          Ed1_Var_coeffs(  3,  8) =    14.92    
          Ed1_Var_coeffs(  3,  9) =   -1.970    
          Ed1_Var_coeffs(  3, 10) =    10.39    
          Ed1_Var_coeffs(  3, 11) =   -21.65    
          Ed1_Var_coeffs(  3, 12) =    1.461    
          Ed1_Var_coeffs(  3, 13) =    10.46    
          Ed1_Var_coeffs(  3, 14) =    54.92    
          Ed1_Var_coeffs(  3, 15) =   -8.033    
          Ed1_Var_coeffs(  4,  1) =    112.6    
          Ed1_Var_coeffs(  4,  2) =   -168.3    
          Ed1_Var_coeffs(  4,  3) =    40.73    
          Ed1_Var_coeffs(  4,  4) =    61.71    
          Ed1_Var_coeffs(  4,  5) =    20.43    
          Ed1_Var_coeffs(  4,  6) =    84.10    
          Ed1_Var_coeffs(  4,  7) =   -19.78    
          Ed1_Var_coeffs(  4,  8) =   -52.41    
          Ed1_Var_coeffs(  4,  9) =    82.06    
          Ed1_Var_coeffs(  4, 10) =   -36.35    
          Ed1_Var_coeffs(  4, 11) =   -40.72    
          Ed1_Var_coeffs(  4, 12) =   -1.303    
          Ed1_Var_coeffs(  4, 13) =    71.04    
          Ed1_Var_coeffs(  4, 14) =    62.12    
          Ed1_Var_coeffs(  4, 15) =    32.98    
          Ed1_Var_coeffs(  5,  1) =    266.4    
          Ed1_Var_coeffs(  5,  2) =   -662.9    
          Ed1_Var_coeffs(  5,  3) =    190.3    
          Ed1_Var_coeffs(  5,  4) =    321.5    
          Ed1_Var_coeffs(  5,  5) =    1141.    
          Ed1_Var_coeffs(  5,  6) =    96.11    
          Ed1_Var_coeffs(  5,  7) =    122.0    
          Ed1_Var_coeffs(  5,  8) =   -95.71    
          Ed1_Var_coeffs(  5,  9) =    69.72    
          Ed1_Var_coeffs(  5, 10) =   -39.64    
          Ed1_Var_coeffs(  5, 11) =   -567.8    
          Ed1_Var_coeffs(  5, 12) =    127.8    
          Ed1_Var_coeffs(  5, 13) =    184.8    
          Ed1_Var_coeffs(  5, 14) =    313.3    
          Ed1_Var_coeffs(  5, 15) =    265.4    
          Ed1_Var_coeffs(  6,  1) =    347.0    
          Ed1_Var_coeffs(  6,  2) =   -89.49    
          Ed1_Var_coeffs(  6,  3) =    15.55    
          Ed1_Var_coeffs(  6,  4) =    24.89    
          Ed1_Var_coeffs(  6,  5) =    247.1    
          Ed1_Var_coeffs(  6,  6) =    153.3    
          Ed1_Var_coeffs(  6,  7) =   -27.99    
          Ed1_Var_coeffs(  6,  8) =   -268.7    
          Ed1_Var_coeffs(  6,  9) =    14.60    
          Ed1_Var_coeffs(  6, 10) =   -49.19    
          Ed1_Var_coeffs(  6, 11) =   -14.69    
          Ed1_Var_coeffs(  6, 12) =   -52.02    
          Ed1_Var_coeffs(  6, 13) =   -70.97    
          Ed1_Var_coeffs(  6, 14) =    43.49    
          Ed1_Var_coeffs(  6, 15) =    27.89    
          Ed1_Var_coeffs(  7,  1) =    649.6    
          Ed1_Var_coeffs(  7,  2) =    377.1    
          Ed1_Var_coeffs(  7,  3) =    39.00    
          Ed1_Var_coeffs(  7,  4) =   -716.2    
          Ed1_Var_coeffs(  7,  5) =   -823.3    
          Ed1_Var_coeffs(  7,  6) =   -85.74    
          Ed1_Var_coeffs(  7,  7) =    368.4    
          Ed1_Var_coeffs(  7,  8) =    292.3    
          Ed1_Var_coeffs(  7,  9) =    130.2    
          Ed1_Var_coeffs(  7, 10) =    564.8    
          Ed1_Var_coeffs(  7, 11) =    680.7    
          Ed1_Var_coeffs(  7, 12) =   -106.4    
          Ed1_Var_coeffs(  7, 13) =   -461.7    
          Ed1_Var_coeffs(  7, 14) =   -630.8    
          Ed1_Var_coeffs(  7, 15) =   -63.38    
          Ed1_Var_coeffs(  8,  1) =    367.6    
          Ed1_Var_coeffs(  8,  2) =    11.94    
          Ed1_Var_coeffs(  8,  3) =   -33.31    
          Ed1_Var_coeffs(  8,  4) =   -100.7    
          Ed1_Var_coeffs(  8,  5) =   -42.09    
          Ed1_Var_coeffs(  8,  6) =   -63.28    
          Ed1_Var_coeffs(  8,  7) =    19.55    
          Ed1_Var_coeffs(  8,  8) =   -78.93    
          Ed1_Var_coeffs(  8,  9) =    39.04    
          Ed1_Var_coeffs(  8, 10) =    281.7    
          Ed1_Var_coeffs(  8, 11) =    56.94    
          Ed1_Var_coeffs(  8, 12) =   -197.5    
          Ed1_Var_coeffs(  8, 13) =   -59.09    
          Ed1_Var_coeffs(  8, 14) =    21.34    
          Ed1_Var_coeffs(  8, 15) =    86.05    
          Ed1_Var_coeffs(  9,  1) =    327.6    
          Ed1_Var_coeffs(  9,  2) =   -47.94    
          Ed1_Var_coeffs(  9,  3) =   -314.5    
          Ed1_Var_coeffs(  9,  4) =    94.05    
          Ed1_Var_coeffs(  9,  5) =    5.883    
          Ed1_Var_coeffs(  9,  6) =   -463.7    
          Ed1_Var_coeffs(  9,  7) =   -86.29    
          Ed1_Var_coeffs(  9,  8) =   -150.6    
          Ed1_Var_coeffs(  9,  9) =    31.17    
          Ed1_Var_coeffs(  9, 10) =    127.3    
          Ed1_Var_coeffs(  9, 11) =   -73.36    
          Ed1_Var_coeffs(  9, 12) =   -538.6    
          Ed1_Var_coeffs(  9, 13) =   -131.7    
          Ed1_Var_coeffs(  9, 14) =   -80.08    
          Ed1_Var_coeffs(  9, 15) =   -34.46    
          ! Fill Ed2_Var arrays:
          Ed2_Var_coeffs(  1,  1) =    78.78    
          Ed2_Var_coeffs(  1,  2) =   -6.681    
          Ed2_Var_coeffs(  1,  3) =    20.31    
          Ed2_Var_coeffs(  1,  4) =    9.724    
          Ed2_Var_coeffs(  1,  5) =   -139.0    
          Ed2_Var_coeffs(  1,  6) =   -2.633    
          Ed2_Var_coeffs(  1,  7) =   -124.0    
          Ed2_Var_coeffs(  1,  8) =   -17.52    
          Ed2_Var_coeffs(  1,  9) =   -71.43    
          Ed2_Var_coeffs(  1, 10) =   -98.30    
          Ed2_Var_coeffs(  1, 11) =   -110.1    
          Ed2_Var_coeffs(  1, 12) =   -68.77    
          Ed2_Var_coeffs(  1, 13) =   -84.93    
          Ed2_Var_coeffs(  1, 14) =    33.75    
          Ed2_Var_coeffs(  1, 15) =   -73.05    
          Ed2_Var_coeffs(  2,  1) =    165.3    
          Ed2_Var_coeffs(  2,  2) =   -100.7    
          Ed2_Var_coeffs(  2,  3) =   -2.699    
          Ed2_Var_coeffs(  2,  4) =    139.7    
          Ed2_Var_coeffs(  2,  5) =   -201.0    
          Ed2_Var_coeffs(  2,  6) =   -39.19    
          Ed2_Var_coeffs(  2,  7) =   -38.17    
          Ed2_Var_coeffs(  2,  8) =    421.4    
          Ed2_Var_coeffs(  2,  9) =   -65.82    
          Ed2_Var_coeffs(  2, 10) =   -56.07    
          Ed2_Var_coeffs(  2, 11) =    9.610    
          Ed2_Var_coeffs(  2, 12) =   -138.4    
          Ed2_Var_coeffs(  2, 13) =   -48.02    
          Ed2_Var_coeffs(  2, 14) =    390.8    
          Ed2_Var_coeffs(  2, 15) =   -80.07    
          Ed2_Var_coeffs(  3,  1) =    138.9    
          Ed2_Var_coeffs(  3,  2) =   -181.5    
          Ed2_Var_coeffs(  3,  3) =   -23.48    
          Ed2_Var_coeffs(  3,  4) =    57.93    
          Ed2_Var_coeffs(  3,  5) =   -322.2    
          Ed2_Var_coeffs(  3,  6) =    63.20    
          Ed2_Var_coeffs(  3,  7) =   -11.21    
          Ed2_Var_coeffs(  3,  8) =    289.4    
          Ed2_Var_coeffs(  3,  9) =    63.41    
          Ed2_Var_coeffs(  3, 10) =   -63.86    
          Ed2_Var_coeffs(  3, 11) =   -221.9    
          Ed2_Var_coeffs(  3, 12) =    76.39    
          Ed2_Var_coeffs(  3, 13) =    63.90    
          Ed2_Var_coeffs(  3, 14) =    402.9    
          Ed2_Var_coeffs(  3, 15) =   -33.12    
          Ed2_Var_coeffs(  4,  1) =    199.8    
          Ed2_Var_coeffs(  4,  2) =   -465.8    
          Ed2_Var_coeffs(  4,  3) =   -16.61    
          Ed2_Var_coeffs(  4,  4) =    149.4    
          Ed2_Var_coeffs(  4,  5) =    167.0    
          Ed2_Var_coeffs(  4,  6) =    130.9    
          Ed2_Var_coeffs(  4,  7) =   -140.6    
          Ed2_Var_coeffs(  4,  8) =   -115.7    
          Ed2_Var_coeffs(  4,  9) =    23.16    
          Ed2_Var_coeffs(  4, 10) =    46.28    
          Ed2_Var_coeffs(  4, 11) =    205.8    
          Ed2_Var_coeffs(  4, 12) =   -3.641    
          Ed2_Var_coeffs(  4, 13) =    69.46    
          Ed2_Var_coeffs(  4, 14) =   -21.97    
          Ed2_Var_coeffs(  4, 15) =    76.34    
          Ed2_Var_coeffs(  5,  1) =    566.4    
          Ed2_Var_coeffs(  5,  2) =   -835.2    
          Ed2_Var_coeffs(  5,  3) =    2.115    
          Ed2_Var_coeffs(  5,  4) =    459.4    
          Ed2_Var_coeffs(  5,  5) =    886.1    
          Ed2_Var_coeffs(  5,  6) =    260.7    
          Ed2_Var_coeffs(  5,  7) =   -117.9    
          Ed2_Var_coeffs(  5,  8) =   -474.8    
          Ed2_Var_coeffs(  5,  9) =    190.3    
          Ed2_Var_coeffs(  5, 10) =    129.5    
          Ed2_Var_coeffs(  5, 11) =    149.0    
          Ed2_Var_coeffs(  5, 12) =    94.19    
          Ed2_Var_coeffs(  5, 13) =    161.2    
          Ed2_Var_coeffs(  5, 14) =    101.5    
          Ed2_Var_coeffs(  5, 15) =    1.250    
          Ed2_Var_coeffs(  6,  1) =    1214.    
          Ed2_Var_coeffs(  6,  2) =    15.23    
          Ed2_Var_coeffs(  6,  3) =   -1.736    
          Ed2_Var_coeffs(  6,  4) =   -306.8    
          Ed2_Var_coeffs(  6,  5) =   -269.5    
          Ed2_Var_coeffs(  6,  6) =    203.2    
          Ed2_Var_coeffs(  6,  7) =    394.2    
          Ed2_Var_coeffs(  6,  8) =   -557.1    
          Ed2_Var_coeffs(  6,  9) =    71.17    
          Ed2_Var_coeffs(  6, 10) =   -245.8    
          Ed2_Var_coeffs(  6, 11) =    473.6    
          Ed2_Var_coeffs(  6, 12) =   -118.3    
          Ed2_Var_coeffs(  6, 13) =   -19.32    
          Ed2_Var_coeffs(  6, 14) =   -309.6    
          Ed2_Var_coeffs(  6, 15) =    40.64    
          Ed2_Var_coeffs(  7,  1) =    1302.    
          Ed2_Var_coeffs(  7,  2) =    561.3    
          Ed2_Var_coeffs(  7,  3) =    86.42    
          Ed2_Var_coeffs(  7,  4) =   -1422.    
          Ed2_Var_coeffs(  7,  5) =   -1454.    
          Ed2_Var_coeffs(  7,  6) =   -14.40    
          Ed2_Var_coeffs(  7,  7) =    199.1    
          Ed2_Var_coeffs(  7,  8) =   -83.23    
          Ed2_Var_coeffs(  7,  9) =    19.13    
          Ed2_Var_coeffs(  7, 10) =    543.4    
          Ed2_Var_coeffs(  7, 11) =    380.1    
          Ed2_Var_coeffs(  7, 12) =   -44.40    
          Ed2_Var_coeffs(  7, 13) =   -56.31    
          Ed2_Var_coeffs(  7, 14) =   -431.5    
          Ed2_Var_coeffs(  7, 15) =    83.93    
          Ed2_Var_coeffs(  8,  1) =    537.0    
          Ed2_Var_coeffs(  8,  2) =    40.26    
          Ed2_Var_coeffs(  8,  3) =    163.4    
          Ed2_Var_coeffs(  8,  4) =   -204.5    
          Ed2_Var_coeffs(  8,  5) =   -402.0    
          Ed2_Var_coeffs(  8,  6) =    117.6    
          Ed2_Var_coeffs(  8,  7) =    208.6    
          Ed2_Var_coeffs(  8,  8) =    144.1    
          Ed2_Var_coeffs(  8,  9) =    80.75    
          Ed2_Var_coeffs(  8, 10) =    90.57    
          Ed2_Var_coeffs(  8, 11) =    186.7    
          Ed2_Var_coeffs(  8, 12) =   -190.8    
          Ed2_Var_coeffs(  8, 13) =    90.50    
          Ed2_Var_coeffs(  8, 14) =   -20.93    
          Ed2_Var_coeffs(  8, 15) =    7.093    
          Ed2_Var_coeffs(  9,  1) =    145.8    
          Ed2_Var_coeffs(  9,  2) =   -204.9    
          Ed2_Var_coeffs(  9,  3) =   -147.5    
          Ed2_Var_coeffs(  9,  4) =   -71.39    
          Ed2_Var_coeffs(  9,  5) =   -217.7    
          Ed2_Var_coeffs(  9,  6) =   -88.63    
          Ed2_Var_coeffs(  9,  7) =    58.56    
          Ed2_Var_coeffs(  9,  8) =   -53.13    
          Ed2_Var_coeffs(  9,  9) =   -67.78    
          Ed2_Var_coeffs(  9, 10) =   -41.76    
          Ed2_Var_coeffs(  9, 11) =   -25.28    
          Ed2_Var_coeffs(  9, 12) =   -15.79    
          Ed2_Var_coeffs(  9, 13) =    23.15    
          Ed2_Var_coeffs(  9, 14) =   -152.2    
          Ed2_Var_coeffs(  9, 15) =   -29.91    
          ! Fill Ed1_Var_DE2 arrays:
          Ed1_Var_DE2_coeffs(  1,  1) =    8.068    
          Ed1_Var_DE2_coeffs(  1,  2) =   -10.59    
          Ed1_Var_DE2_coeffs(  1,  3) =   -3.928    
          Ed1_Var_DE2_coeffs(  1,  4) =   0.2930    
          Ed1_Var_DE2_coeffs(  1,  5) =   -6.077    
          Ed1_Var_DE2_coeffs(  1,  6) =   -4.143    
          Ed1_Var_DE2_coeffs(  1,  7) =   -9.018    
          Ed1_Var_DE2_coeffs(  1,  8) =   -17.58    
          Ed1_Var_DE2_coeffs(  1,  9) =  -0.3879    
          Ed1_Var_DE2_coeffs(  1, 10) =  -0.4107    
          Ed1_Var_DE2_coeffs(  1, 11) =   -12.15    
          Ed1_Var_DE2_coeffs(  1, 12) =   -9.265    
          Ed1_Var_DE2_coeffs(  1, 13) =   -2.101    
          Ed1_Var_DE2_coeffs(  1, 14) =    11.82    
          Ed1_Var_DE2_coeffs(  1, 15) =  -0.2306    
          Ed1_Var_DE2_coeffs(  2,  1) =    51.39    
          Ed1_Var_DE2_coeffs(  2,  2) =    31.72    
          Ed1_Var_DE2_coeffs(  2,  3) =    2.920    
          Ed1_Var_DE2_coeffs(  2,  4) =    38.62    
          Ed1_Var_DE2_coeffs(  2,  5) =    20.58    
          Ed1_Var_DE2_coeffs(  2,  6) =   -6.285    
          Ed1_Var_DE2_coeffs(  2,  7) =   -54.86    
          Ed1_Var_DE2_coeffs(  2,  8) =   -51.15    
          Ed1_Var_DE2_coeffs(  2,  9) =   -41.10    
          Ed1_Var_DE2_coeffs(  2, 10) =   -23.30    
          Ed1_Var_DE2_coeffs(  2, 11) =   -34.78    
          Ed1_Var_DE2_coeffs(  2, 12) =   -47.81    
          Ed1_Var_DE2_coeffs(  2, 13) =   -49.84    
          Ed1_Var_DE2_coeffs(  2, 14) =   -25.82    
          Ed1_Var_DE2_coeffs(  2, 15) =   -44.33    
          Ed1_Var_DE2_coeffs(  3,  1) =    48.25    
          Ed1_Var_DE2_coeffs(  3,  2) =   -14.62    
          Ed1_Var_DE2_coeffs(  3,  3) =   -3.166    
          Ed1_Var_DE2_coeffs(  3,  4) =    41.83    
          Ed1_Var_DE2_coeffs(  3,  5) =   -13.19    
          Ed1_Var_DE2_coeffs(  3,  6) =    11.19    
          Ed1_Var_DE2_coeffs(  3,  7) =   -1.619    
          Ed1_Var_DE2_coeffs(  3,  8) =    22.72    
          Ed1_Var_DE2_coeffs(  3,  9) =   -1.415    
          Ed1_Var_DE2_coeffs(  3, 10) =    11.00    
          Ed1_Var_DE2_coeffs(  3, 11) =   -24.60    
          Ed1_Var_DE2_coeffs(  3, 12) =    3.362    
          Ed1_Var_DE2_coeffs(  3, 13) =    15.91    
          Ed1_Var_DE2_coeffs(  3, 14) =    55.33    
          Ed1_Var_DE2_coeffs(  3, 15) =   -9.681    
          Ed1_Var_DE2_coeffs(  4,  1) =    111.1    
          Ed1_Var_DE2_coeffs(  4,  2) =   -224.6    
          Ed1_Var_DE2_coeffs(  4,  3) =    57.17    
          Ed1_Var_DE2_coeffs(  4,  4) =    74.24    
          Ed1_Var_DE2_coeffs(  4,  5) =   -6.800    
          Ed1_Var_DE2_coeffs(  4,  6) =    91.10    
          Ed1_Var_DE2_coeffs(  4,  7) =   -10.12    
          Ed1_Var_DE2_coeffs(  4,  8) =   -17.88    
          Ed1_Var_DE2_coeffs(  4,  9) =    68.58    
          Ed1_Var_DE2_coeffs(  4, 10) =   -33.24    
          Ed1_Var_DE2_coeffs(  4, 11) =   -76.92    
          Ed1_Var_DE2_coeffs(  4, 12) =    13.61    
          Ed1_Var_DE2_coeffs(  4, 13) =    56.99    
          Ed1_Var_DE2_coeffs(  4, 14) =    34.54    
          Ed1_Var_DE2_coeffs(  4, 15) =    49.95    
          Ed1_Var_DE2_coeffs(  5,  1) =    313.2    
          Ed1_Var_DE2_coeffs(  5,  2) =   -808.5    
          Ed1_Var_DE2_coeffs(  5,  3) =    207.5    
          Ed1_Var_DE2_coeffs(  5,  4) =    320.2    
          Ed1_Var_DE2_coeffs(  5,  5) =    1003.    
          Ed1_Var_DE2_coeffs(  5,  6) =    108.0    
          Ed1_Var_DE2_coeffs(  5,  7) =    94.64    
          Ed1_Var_DE2_coeffs(  5,  8) =   -34.07    
          Ed1_Var_DE2_coeffs(  5,  9) =    39.63    
          Ed1_Var_DE2_coeffs(  5, 10) =   -6.051    
          Ed1_Var_DE2_coeffs(  5, 11) =   -614.7    
          Ed1_Var_DE2_coeffs(  5, 12) =    138.5    
          Ed1_Var_DE2_coeffs(  5, 13) =    221.8    
          Ed1_Var_DE2_coeffs(  5, 14) =    346.2    
          Ed1_Var_DE2_coeffs(  5, 15) =    256.4    
          Ed1_Var_DE2_coeffs(  6,  1) =    524.1    
          Ed1_Var_DE2_coeffs(  6,  2) =   -361.2    
          Ed1_Var_DE2_coeffs(  6,  3) =   -25.67    
          Ed1_Var_DE2_coeffs(  6,  4) =    13.74    
          Ed1_Var_DE2_coeffs(  6,  5) =    306.1    
          Ed1_Var_DE2_coeffs(  6,  6) =    187.2    
          Ed1_Var_DE2_coeffs(  6,  7) =   -93.60    
          Ed1_Var_DE2_coeffs(  6,  8) =   -293.5    
          Ed1_Var_DE2_coeffs(  6,  9) =   -9.537    
          Ed1_Var_DE2_coeffs(  6, 10) =    84.14    
          Ed1_Var_DE2_coeffs(  6, 11) =   -253.2    
          Ed1_Var_DE2_coeffs(  6, 12) =   -83.81    
          Ed1_Var_DE2_coeffs(  6, 13) =    40.80    
          Ed1_Var_DE2_coeffs(  6, 14) =    382.7    
          Ed1_Var_DE2_coeffs(  6, 15) =   -41.65    
          Ed1_Var_DE2_coeffs(  7,  1) =    855.0    
          Ed1_Var_DE2_coeffs(  7,  2) =    6.662    
          Ed1_Var_DE2_coeffs(  7,  3) =   0.1723    
          Ed1_Var_DE2_coeffs(  7,  4) =   -560.5    
          Ed1_Var_DE2_coeffs(  7,  5) =   -418.8    
          Ed1_Var_DE2_coeffs(  7,  6) =   -154.6    
          Ed1_Var_DE2_coeffs(  7,  7) =    99.41    
          Ed1_Var_DE2_coeffs(  7,  8) =   -4.664    
          Ed1_Var_DE2_coeffs(  7,  9) =    122.8    
          Ed1_Var_DE2_coeffs(  7, 10) =    507.5    
          Ed1_Var_DE2_coeffs(  7, 11) =    23.34    
          Ed1_Var_DE2_coeffs(  7, 12) =   -146.8    
          Ed1_Var_DE2_coeffs(  7, 13) =   -509.1    
          Ed1_Var_DE2_coeffs(  7, 14) =   -311.1    
          Ed1_Var_DE2_coeffs(  7, 15) =   -162.6    
          Ed1_Var_DE2_coeffs(  8,  1) =    628.2    
          Ed1_Var_DE2_coeffs(  8,  2) =   -184.3    
          Ed1_Var_DE2_coeffs(  8,  3) =   -148.6    
          Ed1_Var_DE2_coeffs(  8,  4) =   -146.7    
          Ed1_Var_DE2_coeffs(  8,  5) =    28.20    
          Ed1_Var_DE2_coeffs(  8,  6) =   -103.6    
          Ed1_Var_DE2_coeffs(  8,  7) =   -235.9    
          Ed1_Var_DE2_coeffs(  8,  8) =   -274.0    
          Ed1_Var_DE2_coeffs(  8,  9) =    87.38    
          Ed1_Var_DE2_coeffs(  8, 10) =    270.4    
          Ed1_Var_DE2_coeffs(  8, 11) =   -397.9    
          Ed1_Var_DE2_coeffs(  8, 12) =   -277.2    
          Ed1_Var_DE2_coeffs(  8, 13) =   -38.63    
          Ed1_Var_DE2_coeffs(  8, 14) =    347.4    
          Ed1_Var_DE2_coeffs(  8, 15) =   -293.9    
          Ed1_Var_DE2_coeffs(  9,  1) =    683.3    
          Ed1_Var_DE2_coeffs(  9,  2) =   -181.3    
          Ed1_Var_DE2_coeffs(  9,  3) =   -402.3    
          Ed1_Var_DE2_coeffs(  9,  4) =    39.17    
          Ed1_Var_DE2_coeffs(  9,  5) =   -13.36    
          Ed1_Var_DE2_coeffs(  9,  6) =   -537.6    
          Ed1_Var_DE2_coeffs(  9,  7) =   -200.0    
          Ed1_Var_DE2_coeffs(  9,  8) =   -233.0    
          Ed1_Var_DE2_coeffs(  9,  9) =    80.67    
          Ed1_Var_DE2_coeffs(  9, 10) =    248.9    
          Ed1_Var_DE2_coeffs(  9, 11) =   -442.9    
          Ed1_Var_DE2_coeffs(  9, 12) =   -430.2    
          Ed1_Var_DE2_coeffs(  9, 13) =   -101.9    
          Ed1_Var_DE2_coeffs(  9, 14) =    80.37    
          Ed1_Var_DE2_coeffs(  9, 15) =   -220.9    
          ! Fill Ed2_Var_DE2 arrays:
          Ed2_Var_DE2_coeffs(  1,  1) =    75.60    
          Ed2_Var_DE2_coeffs(  1,  2) =   -16.41    
          Ed2_Var_DE2_coeffs(  1,  3) =    21.62    
          Ed2_Var_DE2_coeffs(  1,  4) =    16.01    
          Ed2_Var_DE2_coeffs(  1,  5) =   -141.1    
          Ed2_Var_DE2_coeffs(  1,  6) =   -5.861    
          Ed2_Var_DE2_coeffs(  1,  7) =   -113.0    
          Ed2_Var_DE2_coeffs(  1,  8) =  -0.7517    
          Ed2_Var_DE2_coeffs(  1,  9) =   -73.89    
          Ed2_Var_DE2_coeffs(  1, 10) =   -89.44    
          Ed2_Var_DE2_coeffs(  1, 11) =   -104.7    
          Ed2_Var_DE2_coeffs(  1, 12) =   -72.49    
          Ed2_Var_DE2_coeffs(  1, 13) =   -84.28    
          Ed2_Var_DE2_coeffs(  1, 14) =    39.56    
          Ed2_Var_DE2_coeffs(  1, 15) =   -68.55    
          Ed2_Var_DE2_coeffs(  2,  1) =    159.8    
          Ed2_Var_DE2_coeffs(  2,  2) =   -127.1    
          Ed2_Var_DE2_coeffs(  2,  3) =   -1.339    
          Ed2_Var_DE2_coeffs(  2,  4) =    148.1    
          Ed2_Var_DE2_coeffs(  2,  5) =   -219.0    
          Ed2_Var_DE2_coeffs(  2,  6) =   -44.29    
          Ed2_Var_DE2_coeffs(  2,  7) =   -32.39    
          Ed2_Var_DE2_coeffs(  2,  8) =    441.9    
          Ed2_Var_DE2_coeffs(  2,  9) =   -73.01    
          Ed2_Var_DE2_coeffs(  2, 10) =   -46.28    
          Ed2_Var_DE2_coeffs(  2, 11) =    10.30    
          Ed2_Var_DE2_coeffs(  2, 12) =   -141.0    
          Ed2_Var_DE2_coeffs(  2, 13) =   -53.04    
          Ed2_Var_DE2_coeffs(  2, 14) =    403.0    
          Ed2_Var_DE2_coeffs(  2, 15) =   -77.12    
          Ed2_Var_DE2_coeffs(  3,  1) =    253.7    
          Ed2_Var_DE2_coeffs(  3,  2) =   -188.0    
          Ed2_Var_DE2_coeffs(  3,  3) =   -54.41    
          Ed2_Var_DE2_coeffs(  3,  4) =    59.80    
          Ed2_Var_DE2_coeffs(  3,  5) =   -356.8    
          Ed2_Var_DE2_coeffs(  3,  6) =    99.23    
          Ed2_Var_DE2_coeffs(  3,  7) =   -103.3    
          Ed2_Var_DE2_coeffs(  3,  8) =    452.1    
          Ed2_Var_DE2_coeffs(  3,  9) =    110.0    
          Ed2_Var_DE2_coeffs(  3, 10) =   -171.1    
          Ed2_Var_DE2_coeffs(  3, 11) =   -168.1    
          Ed2_Var_DE2_coeffs(  3, 12) =    115.6    
          Ed2_Var_DE2_coeffs(  3, 13) =    127.0    
          Ed2_Var_DE2_coeffs(  3, 14) =    650.7    
          Ed2_Var_DE2_coeffs(  3, 15) =   -97.98    
          Ed2_Var_DE2_coeffs(  4,  1) =    552.4    
          Ed2_Var_DE2_coeffs(  4,  2) =   -667.7    
          Ed2_Var_DE2_coeffs(  4,  3) =   -14.98    
          Ed2_Var_DE2_coeffs(  4,  4) =    225.0    
          Ed2_Var_DE2_coeffs(  4,  5) =    356.5    
          Ed2_Var_DE2_coeffs(  4,  6) =    184.7    
          Ed2_Var_DE2_coeffs(  4,  7) =   -413.2    
          Ed2_Var_DE2_coeffs(  4,  8) =   -174.9    
          Ed2_Var_DE2_coeffs(  4,  9) =    136.0    
          Ed2_Var_DE2_coeffs(  4, 10) =   -312.7    
          Ed2_Var_DE2_coeffs(  4, 11) =    311.0    
          Ed2_Var_DE2_coeffs(  4, 12) =    14.79    
          Ed2_Var_DE2_coeffs(  4, 13) =    203.3    
          Ed2_Var_DE2_coeffs(  4, 14) =   -218.3    
          Ed2_Var_DE2_coeffs(  4, 15) =    20.32    
          Ed2_Var_DE2_coeffs(  5,  1) =    885.8    
          Ed2_Var_DE2_coeffs(  5,  2) =   -1411.    
          Ed2_Var_DE2_coeffs(  5,  3) =    126.0    
          Ed2_Var_DE2_coeffs(  5,  4) =    476.2    
          Ed2_Var_DE2_coeffs(  5,  5) =    1241.    
          Ed2_Var_DE2_coeffs(  5,  6) =    160.7    
          Ed2_Var_DE2_coeffs(  5,  7) =   -227.1    
          Ed2_Var_DE2_coeffs(  5,  8) =   -859.5    
          Ed2_Var_DE2_coeffs(  5,  9) =    79.25    
          Ed2_Var_DE2_coeffs(  5, 10) =   -143.3    
          Ed2_Var_DE2_coeffs(  5, 11) =    273.2    
          Ed2_Var_DE2_coeffs(  5, 12) =   -14.20    
          Ed2_Var_DE2_coeffs(  5, 13) =    445.6    
          Ed2_Var_DE2_coeffs(  5, 14) =   -360.5    
          Ed2_Var_DE2_coeffs(  5, 15) =    173.2    
          Ed2_Var_DE2_coeffs(  6,  1) =    1674.    
          Ed2_Var_DE2_coeffs(  6,  2) =    302.2    
          Ed2_Var_DE2_coeffs(  6,  3) =    75.52    
          Ed2_Var_DE2_coeffs(  6,  4) =   -564.7    
          Ed2_Var_DE2_coeffs(  6,  5) =   -469.1    
          Ed2_Var_DE2_coeffs(  6,  6) =    184.3    
          Ed2_Var_DE2_coeffs(  6,  7) =    306.2    
          Ed2_Var_DE2_coeffs(  6,  8) =   -930.8    
          Ed2_Var_DE2_coeffs(  6,  9) =   -100.8    
          Ed2_Var_DE2_coeffs(  6, 10) =   -226.3    
          Ed2_Var_DE2_coeffs(  6, 11) =    492.1    
          Ed2_Var_DE2_coeffs(  6, 12) =    67.65    
          Ed2_Var_DE2_coeffs(  6, 13) =    18.33    
          Ed2_Var_DE2_coeffs(  6, 14) =   -294.5    
          Ed2_Var_DE2_coeffs(  6, 15) =    285.4    
          Ed2_Var_DE2_coeffs(  7,  1) =    2206.    
          Ed2_Var_DE2_coeffs(  7,  2) =    1011.    
          Ed2_Var_DE2_coeffs(  7,  3) =    219.8    
          Ed2_Var_DE2_coeffs(  7,  4) =   -2178.    
          Ed2_Var_DE2_coeffs(  7,  5) =   -2366.    
          Ed2_Var_DE2_coeffs(  7,  6) =   -39.31    
          Ed2_Var_DE2_coeffs(  7,  7) =    341.8    
          Ed2_Var_DE2_coeffs(  7,  8) =   -351.8    
          Ed2_Var_DE2_coeffs(  7,  9) =    535.0    
          Ed2_Var_DE2_coeffs(  7, 10) =    934.6    
          Ed2_Var_DE2_coeffs(  7, 11) =    1176.    
          Ed2_Var_DE2_coeffs(  7, 12) =    221.1    
          Ed2_Var_DE2_coeffs(  7, 13) =   -192.6    
          Ed2_Var_DE2_coeffs(  7, 14) =   -820.8    
          Ed2_Var_DE2_coeffs(  7, 15) =    196.3    
          Ed2_Var_DE2_coeffs(  8,  1) =    1085.    
          Ed2_Var_DE2_coeffs(  8,  2) =    10.59    
          Ed2_Var_DE2_coeffs(  8,  3) =    200.0    
          Ed2_Var_DE2_coeffs(  8,  4) =    96.46    
          Ed2_Var_DE2_coeffs(  8,  5) =   -334.9    
          Ed2_Var_DE2_coeffs(  8,  6) =    208.8    
          Ed2_Var_DE2_coeffs(  8,  7) =    190.9    
          Ed2_Var_DE2_coeffs(  8,  8) =   -424.2    
          Ed2_Var_DE2_coeffs(  8,  9) =    446.6    
          Ed2_Var_DE2_coeffs(  8, 10) =    55.33    
          Ed2_Var_DE2_coeffs(  8, 11) =    428.0    
          Ed2_Var_DE2_coeffs(  8, 12) =    189.8    
          Ed2_Var_DE2_coeffs(  8, 13) =    73.65    
          Ed2_Var_DE2_coeffs(  8, 14) =   -435.9    
          Ed2_Var_DE2_coeffs(  8, 15) =    287.3    
          Ed2_Var_DE2_coeffs(  9,  1) =    677.6    
          Ed2_Var_DE2_coeffs(  9,  2) =   -328.8    
          Ed2_Var_DE2_coeffs(  9,  3) =   -79.70    
          Ed2_Var_DE2_coeffs(  9,  4) =   -30.32    
          Ed2_Var_DE2_coeffs(  9,  5) =   -333.5    
          Ed2_Var_DE2_coeffs(  9,  6) =   -113.1    
          Ed2_Var_DE2_coeffs(  9,  7) =   -116.8    
          Ed2_Var_DE2_coeffs(  9,  8) =   -484.1    
          Ed2_Var_DE2_coeffs(  9,  9) =    107.2    
          Ed2_Var_DE2_coeffs(  9, 10) =   -379.2    
          Ed2_Var_DE2_coeffs(  9, 11) =    11.89    
          Ed2_Var_DE2_coeffs(  9, 12) =   -84.65    
          Ed2_Var_DE2_coeffs(  9, 13) =   -217.6    
          Ed2_Var_DE2_coeffs(  9, 14) =   -825.9    
          Ed2_Var_DE2_coeffs(  9, 15) =    317.7    
          ! Fill Ed1_Var_mod arrays:
          Ed1_Var_mod_coeffs(  1,  1) =    1.702    
          Ed1_Var_mod_coeffs(  1,  2) =   0.9159    
          Ed1_Var_mod_coeffs(  1,  3) =  -0.4815E-03
          Ed1_Var_mod_coeffs(  1,  4) =   -1.025    
          Ed1_Var_mod_coeffs(  1,  5) =   -1.610    
          Ed1_Var_mod_coeffs(  1,  6) =   0.4577    
          Ed1_Var_mod_coeffs(  1,  7) =   -1.698    
          Ed1_Var_mod_coeffs(  1,  8) =   -1.790    
          Ed1_Var_mod_coeffs(  1,  9) =   0.9723    
          Ed1_Var_mod_coeffs(  1, 10) =  -0.7679    
          Ed1_Var_mod_coeffs(  1, 11) =   -1.146    
          Ed1_Var_mod_coeffs(  1, 12) =   0.7474    
          Ed1_Var_mod_coeffs(  1, 13) =  -0.8778    
          Ed1_Var_mod_coeffs(  1, 14) =   -1.795    
          Ed1_Var_mod_coeffs(  1, 15) =   0.6407    
          Ed1_Var_mod_coeffs(  2,  1) =    2.480    
          Ed1_Var_mod_coeffs(  2,  2) =    1.608    
          Ed1_Var_mod_coeffs(  2,  3) =  -0.4779    
          Ed1_Var_mod_coeffs(  2,  4) =  -0.9949    
          Ed1_Var_mod_coeffs(  2,  5) =   -1.581    
          Ed1_Var_mod_coeffs(  2,  6) =   0.9520E-01
          Ed1_Var_mod_coeffs(  2,  7) =   -2.884    
          Ed1_Var_mod_coeffs(  2,  8) =   -3.583    
          Ed1_Var_mod_coeffs(  2,  9) =    1.750    
          Ed1_Var_mod_coeffs(  2, 10) =   -2.450    
          Ed1_Var_mod_coeffs(  2, 11) =   -4.001    
          Ed1_Var_mod_coeffs(  2, 12) =    1.633    
          Ed1_Var_mod_coeffs(  2, 13) =   -2.249    
          Ed1_Var_mod_coeffs(  2, 14) =   -3.882    
          Ed1_Var_mod_coeffs(  2, 15) =   0.6752    
          Ed1_Var_mod_coeffs(  3,  1) =    3.497    
          Ed1_Var_mod_coeffs(  3,  2) =  -0.1197    
          Ed1_Var_mod_coeffs(  3,  3) =  -0.5793E-01
          Ed1_Var_mod_coeffs(  3,  4) =  -0.8630    
          Ed1_Var_mod_coeffs(  3,  5) =   -2.012    
          Ed1_Var_mod_coeffs(  3,  6) =  -0.3250    
          Ed1_Var_mod_coeffs(  3,  7) =   -5.107    
          Ed1_Var_mod_coeffs(  3,  8) =   -5.448    
          Ed1_Var_mod_coeffs(  3,  9) =    2.202    
          Ed1_Var_mod_coeffs(  3, 10) =   -4.268    
          Ed1_Var_mod_coeffs(  3, 11) =   -7.569    
          Ed1_Var_mod_coeffs(  3, 12) =    3.254    
          Ed1_Var_mod_coeffs(  3, 13) =   -2.627    
          Ed1_Var_mod_coeffs(  3, 14) =   -6.983    
          Ed1_Var_mod_coeffs(  3, 15) =    1.514    
          Ed1_Var_mod_coeffs(  4,  1) =    8.567    
          Ed1_Var_mod_coeffs(  4,  2) =   -13.53    
          Ed1_Var_mod_coeffs(  4,  3) =    3.608    
          Ed1_Var_mod_coeffs(  4,  4) =    1.006    
          Ed1_Var_mod_coeffs(  4,  5) =   -10.09    
          Ed1_Var_mod_coeffs(  4,  6) =    2.177    
          Ed1_Var_mod_coeffs(  4,  7) =   -7.990    
          Ed1_Var_mod_coeffs(  4,  8) =   -1.834    
          Ed1_Var_mod_coeffs(  4,  9) =  -0.5394    
          Ed1_Var_mod_coeffs(  4, 10) =   -1.474    
          Ed1_Var_mod_coeffs(  4, 11) =   -23.72    
          Ed1_Var_mod_coeffs(  4, 12) =    8.686    
          Ed1_Var_mod_coeffs(  4, 13) =    3.167    
          Ed1_Var_mod_coeffs(  4, 14) =   -3.565    
          Ed1_Var_mod_coeffs(  4, 15) =    1.401    
          Ed1_Var_mod_coeffs(  5,  1) =    59.85    
          Ed1_Var_mod_coeffs(  5,  2) =   -49.87    
          Ed1_Var_mod_coeffs(  5,  3) =    21.59    
          Ed1_Var_mod_coeffs(  5,  4) =    23.30    
          Ed1_Var_mod_coeffs(  5,  5) =    7.593    
          Ed1_Var_mod_coeffs(  5,  6) =    16.58    
          Ed1_Var_mod_coeffs(  5,  7) =   -19.25    
          Ed1_Var_mod_coeffs(  5,  8) =    11.05    
          Ed1_Var_mod_coeffs(  5,  9) =   -6.120    
          Ed1_Var_mod_coeffs(  5, 10) =    30.35    
          Ed1_Var_mod_coeffs(  5, 11) =   -102.5    
          Ed1_Var_mod_coeffs(  5, 12) =    32.13    
          Ed1_Var_mod_coeffs(  5, 13) =    19.70    
          Ed1_Var_mod_coeffs(  5, 14) =    81.33    
          Ed1_Var_mod_coeffs(  5, 15) =   -23.85    
          Ed1_Var_mod_coeffs(  6,  1) =    190.5    
          Ed1_Var_mod_coeffs(  6,  2) =   -121.6    
          Ed1_Var_mod_coeffs(  6,  3) =    57.92    
          Ed1_Var_mod_coeffs(  6,  4) =    39.56    
          Ed1_Var_mod_coeffs(  6,  5) =    66.33    
          Ed1_Var_mod_coeffs(  6,  6) =   -10.98    
          Ed1_Var_mod_coeffs(  6,  7) =   -53.19    
          Ed1_Var_mod_coeffs(  6,  8) =   -14.92    
          Ed1_Var_mod_coeffs(  6,  9) =   -8.494    
          Ed1_Var_mod_coeffs(  6, 10) =    137.5    
          Ed1_Var_mod_coeffs(  6, 11) =   -172.3    
          Ed1_Var_mod_coeffs(  6, 12) =    50.10    
          Ed1_Var_mod_coeffs(  6, 13) =    45.86    
          Ed1_Var_mod_coeffs(  6, 14) =    276.5    
          Ed1_Var_mod_coeffs(  6, 15) =   -92.21    
          Ed1_Var_mod_coeffs(  7,  1) =    273.6    
          Ed1_Var_mod_coeffs(  7,  2) =   -201.7    
          Ed1_Var_mod_coeffs(  7,  3) =    59.61    
          Ed1_Var_mod_coeffs(  7,  4) =    39.63    
          Ed1_Var_mod_coeffs(  7,  5) =    97.72    
          Ed1_Var_mod_coeffs(  7,  6) =   -66.84    
          Ed1_Var_mod_coeffs(  7,  7) =   -57.44    
          Ed1_Var_mod_coeffs(  7,  8) =   -24.25    
          Ed1_Var_mod_coeffs(  7,  9) =    6.242    
          Ed1_Var_mod_coeffs(  7, 10) =    128.9    
          Ed1_Var_mod_coeffs(  7, 11) =   -254.2    
          Ed1_Var_mod_coeffs(  7, 12) =    18.67    
          Ed1_Var_mod_coeffs(  7, 13) =    21.18    
          Ed1_Var_mod_coeffs(  7, 14) =    326.1    
          Ed1_Var_mod_coeffs(  7, 15) =   -191.7    
          Ed1_Var_mod_coeffs(  8,  1) =    327.1    
          Ed1_Var_mod_coeffs(  8,  2) =   -171.5    
          Ed1_Var_mod_coeffs(  8,  3) =    46.44    
          Ed1_Var_mod_coeffs(  8,  4) =    43.99    
          Ed1_Var_mod_coeffs(  8,  5) =   -42.01    
          Ed1_Var_mod_coeffs(  8,  6) =   -40.43    
          Ed1_Var_mod_coeffs(  8,  7) =   -46.20    
          Ed1_Var_mod_coeffs(  8,  8) =   -4.953    
          Ed1_Var_mod_coeffs(  8,  9) =    50.58    
          Ed1_Var_mod_coeffs(  8, 10) =    131.9    
          Ed1_Var_mod_coeffs(  8, 11) =   -232.1    
          Ed1_Var_mod_coeffs(  8, 12) =    36.12    
          Ed1_Var_mod_coeffs(  8, 13) =    19.48    
          Ed1_Var_mod_coeffs(  8, 14) =    373.1    
          Ed1_Var_mod_coeffs(  8, 15) =   -330.6    
          Ed1_Var_mod_coeffs(  9,  1) =    363.9    
          Ed1_Var_mod_coeffs(  9,  2) =   -70.92    
          Ed1_Var_mod_coeffs(  9,  3) =    73.38    
          Ed1_Var_mod_coeffs(  9,  4) =    130.5    
          Ed1_Var_mod_coeffs(  9,  5) =    63.29    
          Ed1_Var_mod_coeffs(  9,  6) =   -42.62    
          Ed1_Var_mod_coeffs(  9,  7) =   -14.99    
          Ed1_Var_mod_coeffs(  9,  8) =    45.12    
          Ed1_Var_mod_coeffs(  9,  9) =    69.58    
          Ed1_Var_mod_coeffs(  9, 10) =    263.8    
          Ed1_Var_mod_coeffs(  9, 11) =   -115.9    
          Ed1_Var_mod_coeffs(  9, 12) =    24.85    
          Ed1_Var_mod_coeffs(  9, 13) =    22.46    
          Ed1_Var_mod_coeffs(  9, 14) =    251.2    
          Ed1_Var_mod_coeffs(  9, 15) =   -258.6    
          ! Fill Ed2_Var_mod arrays:
          Ed2_Var_mod_coeffs(  1,  1) =    1.141    
          Ed2_Var_mod_coeffs(  1,  2) =   0.5483    
          Ed2_Var_mod_coeffs(  1,  3) =   0.2466    
          Ed2_Var_mod_coeffs(  1,  4) =  -0.4721    
          Ed2_Var_mod_coeffs(  1,  5) =  -0.5864E-01
          Ed2_Var_mod_coeffs(  1,  6) =  -0.4945    
          Ed2_Var_mod_coeffs(  1,  7) =  -0.3982    
          Ed2_Var_mod_coeffs(  1,  8) =   0.1374    
          Ed2_Var_mod_coeffs(  1,  9) =   0.1848    
          Ed2_Var_mod_coeffs(  1, 10) =  -0.4744    
          Ed2_Var_mod_coeffs(  1, 11) =   -1.008    
          Ed2_Var_mod_coeffs(  1, 12) =  -0.1914    
          Ed2_Var_mod_coeffs(  1, 13) =  -0.1890E-01
          Ed2_Var_mod_coeffs(  1, 14) =  -0.4914    
          Ed2_Var_mod_coeffs(  1, 15) =   0.5456    
          Ed2_Var_mod_coeffs(  2,  1) =    1.867    
          Ed2_Var_mod_coeffs(  2,  2) =   0.4383    
          Ed2_Var_mod_coeffs(  2,  3) =   0.2357    
          Ed2_Var_mod_coeffs(  2,  4) =  -0.1536    
          Ed2_Var_mod_coeffs(  2,  5) =  -0.8979    
          Ed2_Var_mod_coeffs(  2,  6) =  -0.4724    
          Ed2_Var_mod_coeffs(  2,  7) =  -0.6720    
          Ed2_Var_mod_coeffs(  2,  8) =   0.5050    
          Ed2_Var_mod_coeffs(  2,  9) =   0.5687E-01
          Ed2_Var_mod_coeffs(  2, 10) =  -0.2944    
          Ed2_Var_mod_coeffs(  2, 11) =   -1.135    
          Ed2_Var_mod_coeffs(  2, 12) =  -0.1419    
          Ed2_Var_mod_coeffs(  2, 13) =  -0.1592    
          Ed2_Var_mod_coeffs(  2, 14) =   0.3067    
          Ed2_Var_mod_coeffs(  2, 15) =   0.3532E-01
          Ed2_Var_mod_coeffs(  3,  1) =    67.86    
          Ed2_Var_mod_coeffs(  3,  2) =   -5.001    
          Ed2_Var_mod_coeffs(  3,  3) =    1.169    
          Ed2_Var_mod_coeffs(  3,  4) =    18.72    
          Ed2_Var_mod_coeffs(  3,  5) =   -16.11    
          Ed2_Var_mod_coeffs(  3,  6) =    9.685    
          Ed2_Var_mod_coeffs(  3,  7) =   -53.92    
          Ed2_Var_mod_coeffs(  3,  8) =    47.65    
          Ed2_Var_mod_coeffs(  3,  9) =   -3.892    
          Ed2_Var_mod_coeffs(  3, 10) =   -46.76    
          Ed2_Var_mod_coeffs(  3, 11) =    11.22    
          Ed2_Var_mod_coeffs(  3, 12) =   -8.892    
          Ed2_Var_mod_coeffs(  3, 13) =   0.2049E-01
          Ed2_Var_mod_coeffs(  3, 14) =    66.71    
          Ed2_Var_mod_coeffs(  3, 15) =   -36.64    
          Ed2_Var_mod_coeffs(  4,  1) =    260.1    
          Ed2_Var_mod_coeffs(  4,  2) =   -51.34    
          Ed2_Var_mod_coeffs(  4,  3) =    50.84    
          Ed2_Var_mod_coeffs(  4,  4) =    35.39    
          Ed2_Var_mod_coeffs(  4,  5) =    37.02    
          Ed2_Var_mod_coeffs(  4,  6) =   -10.37    
          Ed2_Var_mod_coeffs(  4,  7) =   -173.8    
          Ed2_Var_mod_coeffs(  4,  8) =   -40.87    
          Ed2_Var_mod_coeffs(  4,  9) =    30.49    
          Ed2_Var_mod_coeffs(  4, 10) =   -269.5    
          Ed2_Var_mod_coeffs(  4, 11) =   0.1814E-01
          Ed2_Var_mod_coeffs(  4, 12) =   -38.75    
          Ed2_Var_mod_coeffs(  4, 13) =    41.16    
          Ed2_Var_mod_coeffs(  4, 14) =   -96.79    
          Ed2_Var_mod_coeffs(  4, 15) =   -3.581    
          Ed2_Var_mod_coeffs(  5,  1) =    377.0    
          Ed2_Var_mod_coeffs(  5,  2) =   -162.6    
          Ed2_Var_mod_coeffs(  5,  3) =    89.29    
          Ed2_Var_mod_coeffs(  5,  4) =    15.65    
          Ed2_Var_mod_coeffs(  5,  5) =    141.5    
          Ed2_Var_mod_coeffs(  5,  6) =   -51.15    
          Ed2_Var_mod_coeffs(  5,  7) =   -141.8    
          Ed2_Var_mod_coeffs(  5,  8) =   -63.84    
          Ed2_Var_mod_coeffs(  5,  9) =   -18.80    
          Ed2_Var_mod_coeffs(  5, 10) =   -308.5    
          Ed2_Var_mod_coeffs(  5, 11) =    23.00    
          Ed2_Var_mod_coeffs(  5, 12) =    5.306    
          Ed2_Var_mod_coeffs(  5, 13) =    232.6    
          Ed2_Var_mod_coeffs(  5, 14) =   -143.1    
          Ed2_Var_mod_coeffs(  5, 15) =    67.59    
          Ed2_Var_mod_coeffs(  6,  1) =    301.5    
          Ed2_Var_mod_coeffs(  6,  2) =    187.4    
          Ed2_Var_mod_coeffs(  6,  3) =    22.07    
          Ed2_Var_mod_coeffs(  6,  4) =   -205.1    
          Ed2_Var_mod_coeffs(  6,  5) =   -183.2    
          Ed2_Var_mod_coeffs(  6,  6) =   -10.71    
          Ed2_Var_mod_coeffs(  6,  7) =    90.66    
          Ed2_Var_mod_coeffs(  6,  8) =    59.91    
          Ed2_Var_mod_coeffs(  6,  9) =   -103.1    
          Ed2_Var_mod_coeffs(  6, 10) =   -37.06    
          Ed2_Var_mod_coeffs(  6, 11) =    1.245    
          Ed2_Var_mod_coeffs(  6, 12) =    72.61    
          Ed2_Var_mod_coeffs(  6, 13) =   -14.17    
          Ed2_Var_mod_coeffs(  6, 14) =   -90.59    
          Ed2_Var_mod_coeffs(  6, 15) =    102.3    
          Ed2_Var_mod_coeffs(  7,  1) =    450.5    
          Ed2_Var_mod_coeffs(  7,  2) =    130.8    
          Ed2_Var_mod_coeffs(  7,  3) =    153.9    
          Ed2_Var_mod_coeffs(  7,  4) =   -340.1    
          Ed2_Var_mod_coeffs(  7,  5) =   -384.4    
          Ed2_Var_mod_coeffs(  7,  6) =   -60.72    
          Ed2_Var_mod_coeffs(  7,  7) =    117.3    
          Ed2_Var_mod_coeffs(  7,  8) =   -245.5    
          Ed2_Var_mod_coeffs(  7,  9) =    461.5    
          Ed2_Var_mod_coeffs(  7, 10) =    58.42    
          Ed2_Var_mod_coeffs(  7, 11) =    270.6    
          Ed2_Var_mod_coeffs(  7, 12) =    27.56    
          Ed2_Var_mod_coeffs(  7, 13) =   -55.75    
          Ed2_Var_mod_coeffs(  7, 14) =   -141.3    
          Ed2_Var_mod_coeffs(  7, 15) =   -91.84    
          Ed2_Var_mod_coeffs(  8,  1) =    402.1    
          Ed2_Var_mod_coeffs(  8,  2) =   -111.8    
          Ed2_Var_mod_coeffs(  8,  3) =    113.2    
          Ed2_Var_mod_coeffs(  8,  4) =    173.8    
          Ed2_Var_mod_coeffs(  8,  5) =    107.4    
          Ed2_Var_mod_coeffs(  8,  6) =    80.07    
          Ed2_Var_mod_coeffs(  8,  7) =    87.02    
          Ed2_Var_mod_coeffs(  8,  8) =   -205.8    
          Ed2_Var_mod_coeffs(  8,  9) =    411.9    
          Ed2_Var_mod_coeffs(  8, 10) =   -63.05    
          Ed2_Var_mod_coeffs(  8, 11) =    359.2    
          Ed2_Var_mod_coeffs(  8, 12) =   -1.692    
          Ed2_Var_mod_coeffs(  8, 13) =   -6.291    
          Ed2_Var_mod_coeffs(  8, 14) =   -257.1    
          Ed2_Var_mod_coeffs(  8, 15) =    259.4    
          Ed2_Var_mod_coeffs(  9,  1) =    468.3    
          Ed2_Var_mod_coeffs(  9,  2) =    41.94    
          Ed2_Var_mod_coeffs(  9,  3) =    137.4    
          Ed2_Var_mod_coeffs(  9,  4) =    203.1    
          Ed2_Var_mod_coeffs(  9,  5) =    179.8    
          Ed2_Var_mod_coeffs(  9,  6) =    56.35    
          Ed2_Var_mod_coeffs(  9,  7) =   -34.17    
          Ed2_Var_mod_coeffs(  9,  8) =   -126.1    
          Ed2_Var_mod_coeffs(  9,  9) =    227.5    
          Ed2_Var_mod_coeffs(  9, 10) =   -188.3    
          Ed2_Var_mod_coeffs(  9, 11) =    165.8    
          Ed2_Var_mod_coeffs(  9, 12) =   -58.34    
          Ed2_Var_mod_coeffs(  9, 13) =    5.571    
          Ed2_Var_mod_coeffs(  9, 14) =   -208.3    
          Ed2_Var_mod_coeffs(  9, 15) =    277.7    
      
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
      
        end subroutine SDEd_Bmax_fill_coeff_arrays
      
      end module SDEd_Bmax_coeff_mod
      
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      
      subroutine SDEd_Bmax_choose_coeff_array
      
        use SDEd_Bmax_coeff_mod, only: Ed1_Var_coeffs,Ed2_Var_coeffs, 
     |        SDEd_Bmax_fill_coeff_arrays
        use SDEd_Bmax_model_mod, only: first_Bmax,b1,b2,b_interp1,
     |              b_interp2,x,b_90,nmlat,nfit
      
        integer :: binary_fit_by
      
        !****************************************************************
        ! Fill parameter arrays:
        !****************************************************************
        call SDEd_Bmax_fill_coeff_arrays
      
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
     |              b_interp2(nfit),x(nfit),b_90(nfit))
        first_Bmax = .false.
      
        !****************************************************************
        ! Transfer values:
        !****************************************************************
      
        do imlat=1,nmlat
           do ifit=1,nfit
              b1(imlat,ifit) = Ed1_Var_coeffs(imlat,ifit)
              b2(imlat,ifit) = Ed2_Var_coeffs(imlat,ifit)
           end do
        end do
      
      end subroutine SDEd_Bmax_choose_coeff_array
      
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      
      subroutine SDEd_Bmax_calc_model_value(mlat_glb,mlt,sinT, 
     |       imf_angle,Btrans,value1,value2,value_found)
      
        ! This version is for public distribution, along with the module
        ! file coeff_mod.f90.
        !
        ! Within the /home/hensel/fortran/fit_data/public subdirectory,
        ! this program is called by call_de2_model.
      
        use SDEd_Bmax_coeff_mod, only: average_mlat, fit_by_Btrans,
     |          fit_by_imf_angle, 
     |         fit_by_sinT
        use SDEd_Bmax_model_mod, only: x,b1,b2, b_interp1,b_interp2, 
     |          b_90, imlat,ifit, 
     |         nmlat, nfit, do_smooth, number_of_MLT_coeffs,
     |         number_of_Btrans_coeffs ,number_of_imf_angle_coeffs ,
     |         number_of_sinT_coeffs , fit_by_mlt,
     |         find_MLT_function,find_Btrans_function,
     |          find_IMF_angle_function,
     |         find_sinT_function
        use SDEd_Bmed_model_mod, only:slp1,slpn,islpsw,sigma
      
        implicit none
      
        real,intent(in):: mlt,mlat_glb,sinT, 
     |         imf_angle,Btrans
        real,intent(out):: value1,value2
        logical, intent(out) :: value_found
      
        integer :: fit_count,ierr
        logical :: mlat_found
        integer :: iMLT,iBtrans,iimf_angle,isinT
        real :: MLT_function,Btrans_function,imf_angle_function,
     |         sinT_function,
     |         imf_angle_in, avg_Btrans_array,mlat
      
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
     |               islpsw,b_secderiv1,temp,sigma,ierr) 
           call curv1(nmlat,lat_array,b_coeff2,slp1,slpn, 
     |               islpsw,b_secderiv2,temp,sigma,ierr) 
           b_interp1(ifit) = curv2(mlat,nmlat,lat_array,
     |            b_coeff1,b_secderiv1,sigma)
           b_interp2(ifit) = curv2(mlat,nmlat,lat_array,
     |            b_coeff2,b_secderiv2,sigma)
      
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
     |            fit_by_MLT)
              MLT_function = find_MLT_function(iMLT,MLT)
      
              Btrans_loop2:do iBtrans=1,1+(fit_by_Btrans
     |              *number_of_Btrans_coeffs - fit_by_Btrans)
                 Btrans_function = find_Btrans_function(iBtrans,Btrans)
      
                 imf_angle_loop2:do iimf_angle=1,1+(fit_by_imf_angle
     |                 *number_of_imf_angle_coeffs - fit_by_imf_angle)
                    imf_angle_function = 
     |                  find_imf_angle_function(iimf_angle,imf_angle)
      
                    sinT_loop2:do isinT=1,1+(fit_by_sinT
     |                    *number_of_sinT_coeffs - fit_by_sinT) 
                       sinT_function = find_sinT_function(isinT,sinT)
      
                       x(fit_count) = MLT_function*Btrans_function
     |                       *imf_angle_function*sinT_function
      
                       if (debug) then
                          write(*,'("fit_count = ",I5)') fit_count
                          write(*,'("MLT_function = ",G12.6)') 
     |                      MLT_function
                          write(*,'("Btrans_function = ",G12.6)') 
     |                       Btrans_function
                          write(*,'("imf_angle_function = ",G12.6)') 
     |                          imf_angle_function
                          write(*,'("sinT_function = ",G12.6)')
     |                           sinT_function
                          write(*,'("mlat = ",F9.2," MLT = ",F9.2,
     |                       " Btrans = ",
     |                       F9.2," imf_angle = ",F9.2," sinT = ",F9.2,
     |                       " x(",I5,") = ",G13.6)') 
     |                       mlat,MLT,Btrans,imf_angle,sinT,
     |                       fit_count,x(fit_count)
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
      
      end subroutine SDEd_Bmax_calc_model_value
                                                                                                

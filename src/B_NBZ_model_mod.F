       module B_NBZ_model_mod
       
         use common_model_mod, only: MMX,MXNMX,ITHMX,NQ,pi,fit_by_MLT
       
         implicit none
       
         logical :: first_time
         logical :: value_dummy
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Arrays for fit coefficients:
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Number of fit coefficients for each mlat:
         integer :: B_NBZ_nfit
         ! Fit coefficients (dimension nmlat x nfit):
         real, allocatable, dimension(:,:) :: B_NBZ_b_Ed1,B_NBZ_b_Ed2 ! nmlat x nfit
         ! Interpolated fit coefficients (dimension nfit):
         real, allocatable, dimension(:) :: b_interp_Ed1,b_interp_Ed2
         ! Fitting functions (dimension nfit):
         real, allocatable, dimension(:) :: x
         ! Fit coefficients averaged over the pole (dimension nfit):
         real, allocatable, dimension(:) :: b_90_Ed1,b_90_Ed2
         ! Number of fit functions used by fit_data:
         integer, parameter :: B_NBZ_number_of_MLT_coeffs = 7,
     |         B_NBZ_number_of_Btrans_coeffs = 1,
     |         B_NBZ_number_of_imf_angle_coeffs = 1,
     |         B_NBZ_number_of_sinT_coeffs = 2
         ! Indices:
         integer :: imlat,ifit
         ! FAC:
         real :: NBZ_N_AO(NQ),NBZ_FAC_AO(-MMX:MMX,MXNMX,
     |       B_NBZ_number_of_Btrans_coeffs,
     |       B_NBZ_number_of_imf_angle_coeffs,
     |       B_NBZ_number_of_sinT_coeffs)
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Data smoothing:
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         logical :: do_smooth
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Results
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         real :: value, value_Ed1, value_Ed2, value_dBr, value_dBm
         logical :: value_found, fit_found
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Qgen variables
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         real ::  B_NBZ_AO(-MMX:MMX,MXNMX,
     |       B_NBZ_number_of_Btrans_coeffs,
     |       B_NBZ_number_of_imf_angle_coeffs,
     |       B_NBZ_number_of_sinT_coeffs)
         real :: B_NBZ_Q(NQ,0:ITHMX),B_NBZ_DQ(NQ,0:ITHMX)
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Coeffmerge
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         integer, parameter :: 
     |         B_NBZ_fit_by_Btrans = 0,
     |         B_NBZ_fit_by_imf_angle = 0,
     |         B_NBZ_fit_by_sinT = 1   
       contains
       
         function find_MLT_function (index,MLT_in) result (MLT_fxn_out)
           integer :: index
           real :: MLT_in, MLT_fxn_out
       
           select case (index)
           case (0) 
              MLT_fxn_out = 1.
           case (1)
              MLT_fxn_out = 1.
           case (2)
              MLT_fxn_out = cos(pi*MLT_in/12.)
           case (3)
              MLT_fxn_out = sin(pi*MLT_in/12.)
           case (4)
              MLT_fxn_out = cos(pi*MLT_in/6.)
           case (5)
              MLT_fxn_out = sin(pi*MLT_in/6.)
           case (6)
              MLT_fxn_out = cos(pi*MLT_in/4.)
           case (7)
              MLT_fxn_out = sin(pi*MLT_in/4.)
           end select
         end function find_MLT_function
       
         function find_Btrans_function (index,Btrans_in) 
     |         result (Btrans_fxn_out)
           integer :: index
           real :: Btrans_in, Btrans_fxn_out
       
           select case (index)
           case (0) 
              Btrans_fxn_out = 1.
           case (1)
              Btrans_fxn_out = 1.
           case (2)
              Btrans_fxn_out = Btrans_in/5.2 - 1.
           case (3)
              Btrans_fxn_out = (Btrans_in/5.2 - 1.)**2
           end select
         end function find_Btrans_function
       
         function find_imf_angle_function (index,imf_angle_in) 
     |         result (imf_angle_fxn_out)
           integer :: index
           real :: imf_angle_in, imf_angle_fxn_out
       
           select case (index)
           case (0) 
              imf_angle_fxn_out = 1.
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
       
         function find_sinT_function(index,sinT_in)result(sinT_fxn_out)
           integer :: index
           real :: sinT_in, sinT_fxn_out
       
           select case (index)
           case (0) 
              sinT_fxn_out = 1.
           case (1)
              sinT_fxn_out = 0.5*(1+(4.*sinT_in+0.5)/SQRT(0.75+ 
     |         (4.*sinT_in+0.5)**2))
           case (2)
              sinT_fxn_out = 0.5*(1+(4.*sinT_in+0.5)/SQRT(0.75+ 
     |         (4.*sinT_in+0.5)**2))
              sinT_fxn_out = sinT_in*sinT_fxn_out
           case (3)
              sinT_fxn_out = sinT_in**2
           end select
         end function find_sinT_function
       
       end module B_NBZ_model_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       module B_NBZ_coeffmerge_mod
       
         ! Merged coefficients from dB_m and dB_r to one set
         ! of Fourier series coefficients
       
         ! This module file created 2006/09/13.
         ! with the following parameters:
         !   for 81249-83047 period 
         !   fit_by_Btrans = 0
         !   fit_by_imf_angle = 0
         !   fit_by_sinT = 1
         !   mlat_resolution = 2
       
         ! To use this module, include the statement
         !   use coeffmerge_mod
         ! and the subroutine call
         !   call fill_coeff_arrays
         ! in Fortran90 programs.
         ! 
         ! Compilation will be slow; a Makefile is suggested!
         implicit none
         save
         real,  dimension( 23, 14) :: dB_d1 
         real,  dimension( 23, 14) :: dB_d2 
         real,  dimension( 23) :: min_mlat,max_mlat!,average_mlat
       
       contains
       
         subroutine B_NBZ_fill_coeff_arrays
       
           ! Fill dB_d1 arrays: 
           dB_d1(  1,  1) =   0.3162    
           dB_d1(  1,  2) =  -0.3415E-01
           dB_d1(  1,  3) =  -0.9423    
           dB_d1(  1,  4) =    2.505    
           dB_d1(  1,  5) =  -0.9385E-01
           dB_d1(  1,  6) =   -11.77    
           dB_d1(  1,  7) =   0.4277    
           dB_d1(  1,  8) =   0.7793    
           dB_d1(  1,  9) =    1.119    
           dB_d1(  1, 10) =   -3.209    
           dB_d1(  1, 11) =   -1.024    
           dB_d1(  1, 12) =  -0.8210    
           dB_d1(  1, 13) =  -0.5211    
           dB_d1(  1, 14) =   0.6550E-01
           dB_d1(  2,  1) =   -3.415    
           dB_d1(  2,  2) =    4.582    
           dB_d1(  2,  3) =   -5.163    
           dB_d1(  2,  4) =   -3.342    
           dB_d1(  2,  5) =   -5.056    
           dB_d1(  2,  6) =   -7.585    
           dB_d1(  2,  7) =  -0.3415    
           dB_d1(  2,  8) =    4.317    
           dB_d1(  2,  9) =    7.921    
           dB_d1(  2, 10) =   -18.45    
           dB_d1(  2, 11) =   0.8495    
           dB_d1(  2, 12) =   0.4625    
           dB_d1(  2, 13) =    3.765    
           dB_d1(  2, 14) =   -5.575    
           dB_d1(  3,  1) =   -7.667    
           dB_d1(  3,  2) =    10.09    
           dB_d1(  3,  3) =   -6.409    
           dB_d1(  3,  4) =   -14.11    
           dB_d1(  3,  5) =   -3.458    
           dB_d1(  3,  6) =   -11.06    
           dB_d1(  3,  7) =    1.095    
           dB_d1(  3,  8) =    11.33    
           dB_d1(  3,  9) =    7.381    
           dB_d1(  3, 10) =   -18.95    
           dB_d1(  3, 11) =   -10.66    
           dB_d1(  3, 12) =    21.54    
           dB_d1(  3, 13) =    2.202    
           dB_d1(  3, 14) =    9.505    
           dB_d1(  4,  1) =   -5.232    
           dB_d1(  4,  2) =   -7.327    
           dB_d1(  4,  3) =   -10.06    
           dB_d1(  4,  4) =   0.6845    
           dB_d1(  4,  5) =   -13.71    
           dB_d1(  4,  6) =  -0.4050    
           dB_d1(  4,  7) =   -1.302    
           dB_d1(  4,  8) =    9.373    
           dB_d1(  4,  9) =    6.187    
           dB_d1(  4, 10) =   -3.594    
           dB_d1(  4, 11) =   0.7140    
           dB_d1(  4, 12) =   0.1075    
           dB_d1(  4, 13) =  -0.5070    
           dB_d1(  4, 14) =    15.13    
           dB_d1(  5,  1) =   -5.642    
           dB_d1(  5,  2) =    2.640    
           dB_d1(  5,  3) =   -14.97    
           dB_d1(  5,  4) =   -3.661    
           dB_d1(  5,  5) =   -22.10    
           dB_d1(  5,  6) =    19.53    
           dB_d1(  5,  7) =    1.661    
           dB_d1(  5,  8) =    14.08    
           dB_d1(  5,  9) =    8.195    
           dB_d1(  5, 10) =  -0.1310    
           dB_d1(  5, 11) =    11.66    
           dB_d1(  5, 12) =   -28.85    
           dB_d1(  5, 13) =    4.809    
           dB_d1(  5, 14) =    14.76    
           dB_d1(  6,  1) =   -5.874    
           dB_d1(  6,  2) =    24.45    
           dB_d1(  6,  3) =   -21.67    
           dB_d1(  6,  4) =   -27.37    
           dB_d1(  6,  5) =   -32.65    
           dB_d1(  6,  6) =    61.11    
           dB_d1(  6,  7) =    8.436    
           dB_d1(  6,  8) =    18.45    
           dB_d1(  6,  9) =    14.86    
           dB_d1(  6, 10) =   -28.19    
           dB_d1(  6, 11) =    14.25    
           dB_d1(  6, 12) =   -55.00    
           dB_d1(  6, 13) =    6.825    
           dB_d1(  6, 14) =    27.25    
           dB_d1(  7,  1) =   -10.33    
           dB_d1(  7,  2) =    36.11    
           dB_d1(  7,  3) =   -17.42    
           dB_d1(  7,  4) =   -39.20    
           dB_d1(  7,  5) =   -35.58    
           dB_d1(  7,  6) =    92.20    
           dB_d1(  7,  7) =    21.40    
           dB_d1(  7,  8) =   -19.32    
           dB_d1(  7,  9) =    18.28    
           dB_d1(  7, 10) =   -46.89    
           dB_d1(  7, 11) =    17.40    
           dB_d1(  7, 12) =   -68.40    
           dB_d1(  7, 13) =    12.03    
           dB_d1(  7, 14) =    46.94    
           dB_d1(  8,  1) =   -18.31    
           dB_d1(  8,  2) =    22.84    
           dB_d1(  8,  3) =   -17.79    
           dB_d1(  8,  4) =   -3.385    
           dB_d1(  8,  5) =   -54.84    
           dB_d1(  8,  6) =    130.3    
           dB_d1(  8,  7) =    30.92    
           dB_d1(  8,  8) =   -69.22    
           dB_d1(  8,  9) =    13.01    
           dB_d1(  8, 10) =   -42.81    
           dB_d1(  8, 11) =    17.68    
           dB_d1(  8, 12) =   -47.86    
           dB_d1(  8, 13) =    17.31    
           dB_d1(  8, 14) =    23.06    
           dB_d1(  9,  1) =   -22.05    
           dB_d1(  9,  2) =    29.70    
           dB_d1(  9,  3) =   -3.631    
           dB_d1(  9,  4) =   -11.50    
           dB_d1(  9,  5) =   -47.29    
           dB_d1(  9,  6) =    102.2    
           dB_d1(  9,  7) =    34.19    
           dB_d1(  9,  8) =   -82.50    
           dB_d1(  9,  9) =   0.4400E-01
           dB_d1(  9, 10) =   -7.245    
           dB_d1(  9, 11) =    14.95    
           dB_d1(  9, 12) =   -21.04    
           dB_d1(  9, 13) =    23.92    
           dB_d1(  9, 14) =    10.51    
           dB_d1( 10,  1) =   -34.90    
           dB_d1( 10,  2) =    31.82    
           dB_d1( 10,  3) =   -29.22    
           dB_d1( 10,  4) =    39.82    
           dB_d1( 10,  5) =   -106.9    
           dB_d1( 10,  6) =    226.6    
           dB_d1( 10,  7) =    26.37    
           dB_d1( 10,  8) =   -64.79    
           dB_d1( 10,  9) =   -13.17    
           dB_d1( 10, 10) =   -20.76    
           dB_d1( 10, 11) =    39.39    
           dB_d1( 10, 12) =   -66.24    
           dB_d1( 10, 13) =    35.33    
           dB_d1( 10, 14) =   -20.85    
           dB_d1( 11,  1) =   -35.61    
           dB_d1( 11,  2) =   -8.350    
           dB_d1( 11,  3) =   -25.54    
           dB_d1( 11,  4) =   -59.83    
           dB_d1( 11,  5) =   -178.4    
           dB_d1( 11,  6) =    358.6    
           dB_d1( 11,  7) =    47.33    
           dB_d1( 11,  8) =   -110.2    
           dB_d1( 11,  9) =   -38.67    
           dB_d1( 11, 10) =    18.81    
           dB_d1( 11, 11) =    60.53    
           dB_d1( 11, 12) =   -126.1    
           dB_d1( 11, 13) =    37.60    
           dB_d1( 11, 14) =   -39.75    
           dB_d1( 12,  1) =   -18.05    
           dB_d1( 12,  2) =   -40.12    
           dB_d1( 12,  3) =    10.07    
           dB_d1( 12,  4) =   -36.52    
           dB_d1( 12,  5) =   -158.4    
           dB_d1( 12,  6) =    305.1    
           dB_d1( 12,  7) =    36.15    
           dB_d1( 12,  8) =   -32.27    
           dB_d1( 12,  9) =   -9.250    
           dB_d1( 12, 10) =   -23.46    
           dB_d1( 12, 11) =    23.94    
           dB_d1( 12, 12) =   -47.02    
           dB_d1( 12, 13) =    36.65    
           dB_d1( 12, 14) =    13.60    
           dB_d1( 13,  1) =   -12.84    
           dB_d1( 13,  2) =   -4.890    
           dB_d1( 13,  3) =   -65.80    
           dB_d1( 13,  4) =    96.89    
           dB_d1( 13,  5) =   -152.6    
           dB_d1( 13,  6) =    241.7    
           dB_d1( 13,  7) =   -5.961    
           dB_d1( 13,  8) =    65.05    
           dB_d1( 13,  9) =    21.80    
           dB_d1( 13, 10) =   -53.83    
           dB_d1( 13, 11) =    21.92    
           dB_d1( 13, 12) =   -71.32    
           dB_d1( 13, 13) =    12.45    
           dB_d1( 13, 14) =    66.18    
           dB_d1( 14,  1) =   -13.37    
           dB_d1( 14,  2) =    56.21    
           dB_d1( 14,  3) =   -83.59    
           dB_d1( 14,  4) =    74.50    
           dB_d1( 14,  5) =   -166.1    
           dB_d1( 14,  6) =    274.8    
           dB_d1( 14,  7) =    19.14    
           dB_d1( 14,  8) =    57.89    
           dB_d1( 14,  9) =    51.39    
           dB_d1( 14, 10) =   -96.66    
           dB_d1( 14, 11) =    44.47    
           dB_d1( 14, 12) =   -156.9    
           dB_d1( 14, 13) =    8.800    
           dB_d1( 14, 14) =    11.20    
           dB_d1( 15,  1) =    17.58    
           dB_d1( 15,  2) =   -46.45    
           dB_d1( 15,  3) =   -101.0    
           dB_d1( 15,  4) =    158.5    
           dB_d1( 15,  5) =   -150.6    
           dB_d1( 15,  6) =    192.3    
           dB_d1( 15,  7) =   -1.550    
           dB_d1( 15,  8) =    142.9    
           dB_d1( 15,  9) =    86.88    
           dB_d1( 15, 10) =   -136.7    
           dB_d1( 15, 11) =    51.25    
           dB_d1( 15, 12) =   -165.3    
           dB_d1( 15, 13) =   -8.200    
           dB_d1( 15, 14) =   -20.40    
           dB_d1( 16,  1) =    12.28    
           dB_d1( 16,  2) =   -76.30    
           dB_d1( 16,  3) =   -63.14    
           dB_d1( 16,  4) =    24.68    
           dB_d1( 16,  5) =    51.24    
           dB_d1( 16,  6) =   -117.9    
           dB_d1( 16,  7) =   -26.47    
           dB_d1( 16,  8) =    112.3    
           dB_d1( 16,  9) =   -5.966    
           dB_d1( 16, 10) =    56.08    
           dB_d1( 16, 11) =    61.13    
           dB_d1( 16, 12) =   -208.1    
           dB_d1( 16, 13) =   -31.43    
           dB_d1( 16, 14) =    1.856    
           dB_d1( 17,  1) =   -8.940    
           dB_d1( 17,  2) =   -81.87    
           dB_d1( 17,  3) =   -18.32    
           dB_d1( 17,  4) =    138.0    
           dB_d1( 17,  5) =    231.5    
           dB_d1( 17,  6) =   -229.0    
           dB_d1( 17,  7) =   -34.34    
           dB_d1( 17,  8) =    58.43    
           dB_d1( 17,  9) =   -48.46    
           dB_d1( 17, 10) =    12.20    
           dB_d1( 17, 11) =    64.10    
           dB_d1( 17, 12) =   -63.47    
           dB_d1( 17, 13) =   -36.41    
           dB_d1( 17, 14) =    50.40    
           dB_d1( 18,  1) =   -10.93    
           dB_d1( 18,  2) =   -182.6    
           dB_d1( 18,  3) =    54.79    
           dB_d1( 18,  4) =    294.0    
           dB_d1( 18,  5) =    198.4    
           dB_d1( 18,  6) =   -97.71    
           dB_d1( 18,  7) =   -78.29    
           dB_d1( 18,  8) =    39.80    
           dB_d1( 18,  9) =    36.15    
           dB_d1( 18, 10) =   -64.17    
           dB_d1( 18, 11) =    15.83    
           dB_d1( 18, 12) =    10.40    
           dB_d1( 18, 13) =   -85.39    
           dB_d1( 18, 14) =    129.3    
           dB_d1( 19,  1) =   -37.42    
           dB_d1( 19,  2) =    87.53    
           dB_d1( 19,  3) =    69.30    
           dB_d1( 19,  4) =   -8.980    
           dB_d1( 19,  5) =    137.3    
           dB_d1( 19,  6) =   -112.6    
           dB_d1( 19,  7) =   -51.06    
           dB_d1( 19,  8) =    179.5    
           dB_d1( 19,  9) =    76.28    
           dB_d1( 19, 10) =   -169.9    
           dB_d1( 19, 11) =   -28.60    
           dB_d1( 19, 12) =   -25.00    
           dB_d1( 19, 13) =   -79.61    
           dB_d1( 19, 14) =    95.03    
           dB_d1( 20,  1) =    12.59    
           dB_d1( 20,  2) =    145.3    
           dB_d1( 20,  3) =    35.94    
           dB_d1( 20,  4) =   -199.8    
           dB_d1( 20,  5) =    148.4    
           dB_d1( 20,  6) =   -499.4    
           dB_d1( 20,  7) =    4.517    
           dB_d1( 20,  8) =    27.95    
           dB_d1( 20,  9) =    96.96    
           dB_d1( 20, 10) =    39.95    
           dB_d1( 20, 11) =   -58.00    
           dB_d1( 20, 12) =    175.4    
           dB_d1( 20, 13) =   -78.18    
           dB_d1( 20, 14) =   -133.7    
           dB_d1( 21,  1) =    16.13    
           dB_d1( 21,  2) =   -153.8    
           dB_d1( 21,  3) =   -54.28    
           dB_d1( 21,  4) =    170.8    
           dB_d1( 21,  5) =    45.28    
           dB_d1( 21,  6) =   -127.5    
           dB_d1( 21,  7) =   -22.36    
           dB_d1( 21,  8) =   -57.03    
           dB_d1( 21,  9) =    103.0    
           dB_d1( 21, 10) =    37.05    
           dB_d1( 21, 11) =    31.98    
           dB_d1( 21, 12) =    120.7    
           dB_d1( 21, 13) =   -73.92    
           dB_d1( 21, 14) =   -19.25    
           dB_d1( 22,  1) =   -14.30    
           dB_d1( 22,  2) =   -6.460    
           dB_d1( 22,  3) =    63.53    
           dB_d1( 22,  4) =   -46.06    
           dB_d1( 22,  5) =   0.8350    
           dB_d1( 22,  6) =    34.76    
           dB_d1( 22,  7) =   -37.37    
           dB_d1( 22,  8) =    37.90    
           dB_d1( 22,  9) =    12.99    
           dB_d1( 22, 10) =   -1.955    
           dB_d1( 22, 11) =    35.20    
           dB_d1( 22, 12) =    109.7    
           dB_d1( 22, 13) =   0.2550    
           dB_d1( 22, 14) =    34.86    
           dB_d1( 23,  1) =   -5.813    
           dB_d1( 23,  2) =  -0.8500E-01
           dB_d1( 23,  3) =    73.84    
           dB_d1( 23,  4) =    7.574    
           dB_d1( 23,  5) =    71.40    
           dB_d1( 23,  6) =    22.77    
           dB_d1( 23,  7) =   -3.433    
           dB_d1( 23,  8) =   -12.05    
           dB_d1( 23,  9) =    73.54    
           dB_d1( 23, 10) =    1.405    
           dB_d1( 23, 11) =   -3.960    
           dB_d1( 23, 12) =    31.30    
           dB_d1( 23, 13) =   -14.52    
           dB_d1( 23, 14) =   -7.385    
       
           ! Fill dB_d2 arrays:
           dB_d2(  1,  1) =  -0.5524    
           dB_d2(  1,  2) =    4.669    
           dB_d2(  1,  3) =   -2.560    
           dB_d2(  1,  4) =    8.872    
           dB_d2(  1,  5) =    1.103    
           dB_d2(  1,  6) =   -2.673    
           dB_d2(  1,  7) =   -1.119    
           dB_d2(  1,  8) =    3.209    
           dB_d2(  1,  9) =   0.4277    
           dB_d2(  1, 10) =   0.7793    
           dB_d2(  1, 11) =   0.5211    
           dB_d2(  1, 12) =  -0.6550E-01
           dB_d2(  1, 13) =   -1.024    
           dB_d2(  1, 14) =  -0.8210    
           dB_d2(  2,  1) =   -9.478    
           dB_d2(  2,  2) =    27.69    
           dB_d2(  2,  3) =   -14.05    
           dB_d2(  2,  4) =    48.53    
           dB_d2(  2,  5) =  -0.8065    
           dB_d2(  2,  6) =    7.512    
           dB_d2(  2,  7) =   -7.921    
           dB_d2(  2,  8) =    18.45    
           dB_d2(  2,  9) =  -0.3415    
           dB_d2(  2, 10) =    4.317    
           dB_d2(  2, 11) =   -3.765    
           dB_d2(  2, 12) =    5.575    
           dB_d2(  2, 13) =   0.8495    
           dB_d2(  2, 14) =   0.4625    
           dB_d2(  3,  1) =   -6.799    
           dB_d2(  3,  2) =    26.32    
           dB_d2(  3,  3) =   -12.11    
           dB_d2(  3,  4) =    69.85    
           dB_d2(  3,  5) =    4.949    
           dB_d2(  3,  6) =   -4.767    
           dB_d2(  3,  7) =   -7.381    
           dB_d2(  3,  8) =    18.95    
           dB_d2(  3,  9) =    1.095    
           dB_d2(  3, 10) =    11.33    
           dB_d2(  3, 11) =   -2.202    
           dB_d2(  3, 12) =   -9.505    
           dB_d2(  3, 13) =   -10.66    
           dB_d2(  3, 14) =    21.54    
           dB_d2(  4,  1) =  -0.9430    
           dB_d2(  4,  2) =   -7.007    
           dB_d2(  4,  3) =    3.333    
           dB_d2(  4,  4) =    64.06    
           dB_d2(  4,  5) =   0.4570    
           dB_d2(  4,  6) =   -10.63    
           dB_d2(  4,  7) =   -6.187    
           dB_d2(  4,  8) =    3.594    
           dB_d2(  4,  9) =   -1.302    
           dB_d2(  4, 10) =    9.373    
           dB_d2(  4, 11) =   0.5070    
           dB_d2(  4, 12) =   -15.13    
           dB_d2(  4, 13) =   0.7140    
           dB_d2(  4, 14) =   0.1075    
           dB_d2(  5,  1) =   -6.245    
           dB_d2(  5,  2) =  -0.5617    
           dB_d2(  5,  3) =   0.1195    
           dB_d2(  5,  4) =    88.31    
           dB_d2(  5,  5) =    1.484    
           dB_d2(  5,  6) =   -11.88    
           dB_d2(  5,  7) =   -8.195    
           dB_d2(  5,  8) =   0.1310    
           dB_d2(  5,  9) =    1.661    
           dB_d2(  5, 10) =    14.08    
           dB_d2(  5, 11) =   -4.809    
           dB_d2(  5, 12) =   -14.76    
           dB_d2(  5, 13) =    11.66    
           dB_d2(  5, 14) =   -28.85    
           dB_d2(  6,  1) =   -5.954    
           dB_d2(  6,  2) =    4.815    
           dB_d2(  6,  3) =   -16.63    
           dB_d2(  6,  4) =    104.2    
           dB_d2(  6,  5) =    3.757    
           dB_d2(  6,  6) =   -35.61    
           dB_d2(  6,  7) =   -14.86    
           dB_d2(  6,  8) =    28.19    
           dB_d2(  6,  9) =    8.436    
           dB_d2(  6, 10) =    18.45    
           dB_d2(  6, 11) =   -6.825    
           dB_d2(  6, 12) =   -27.25    
           dB_d2(  6, 13) =    14.25    
           dB_d2(  6, 14) =   -55.00    
           dB_d2(  7,  1) =   -2.850    
           dB_d2(  7,  2) =    14.78    
           dB_d2(  7,  3) =   -30.83    
           dB_d2(  7,  4) =    93.51    
           dB_d2(  7,  5) =   -2.895    
           dB_d2(  7,  6) =   -11.18    
           dB_d2(  7,  7) =   -18.28    
           dB_d2(  7,  8) =    46.89    
           dB_d2(  7,  9) =    21.40    
           dB_d2(  7, 10) =   -19.32    
           dB_d2(  7, 11) =   -12.03    
           dB_d2(  7, 12) =   -46.94    
           dB_d2(  7, 13) =    17.40    
           dB_d2(  7, 14) =   -68.40    
           dB_d2(  8,  1) =    3.734    
           dB_d2(  8,  2) =    35.10    
           dB_d2(  8,  3) =   -41.06    
           dB_d2(  8,  4) =    55.76    
           dB_d2(  8,  5) =   -18.01    
           dB_d2(  8,  6) =    55.22    
           dB_d2(  8,  7) =   -13.01    
           dB_d2(  8,  8) =    42.81    
           dB_d2(  8,  9) =    30.92    
           dB_d2(  8, 10) =   -69.22    
           dB_d2(  8, 11) =   -17.31    
           dB_d2(  8, 12) =   -23.06    
           dB_d2(  8, 13) =    17.68    
           dB_d2(  8, 14) =   -47.86    
           dB_d2(  9,  1) =    5.899    
           dB_d2(  9,  2) =    30.17    
           dB_d2(  9,  3) =   -42.40    
           dB_d2(  9,  4) =    55.74    
           dB_d2(  9,  5) =   -25.15    
           dB_d2(  9,  6) =    79.76    
           dB_d2(  9,  7) =  -0.4400E-01
           dB_d2(  9,  8) =    7.245    
           dB_d2(  9,  9) =    34.19    
           dB_d2(  9, 10) =   -82.50    
           dB_d2(  9, 11) =   -23.92    
           dB_d2(  9, 12) =   -10.51    
           dB_d2(  9, 13) =    14.95    
           dB_d2(  9, 14) =   -21.04    
           dB_d2( 10,  1) =   -7.321    
           dB_d2( 10,  2) =    87.64    
           dB_d2( 10,  3) =   -51.28    
           dB_d2( 10,  4) =    42.54    
           dB_d2( 10,  5) =   -32.31    
           dB_d2( 10,  6) =    86.18    
           dB_d2( 10,  7) =    13.17    
           dB_d2( 10,  8) =    20.76    
           dB_d2( 10,  9) =    26.37    
           dB_d2( 10, 10) =   -64.79    
           dB_d2( 10, 11) =   -35.33    
           dB_d2( 10, 12) =    20.85    
           dB_d2( 10, 13) =    39.39    
           dB_d2( 10, 14) =   -66.24    
           dB_d2( 11,  1) =   -11.64    
           dB_d2( 11,  2) =    98.88    
           dB_d2( 11,  3) =   -53.37    
           dB_d2( 11,  4) =    20.45    
           dB_d2( 11,  5) =   -41.68    
           dB_d2( 11,  6) =    104.2    
           dB_d2( 11,  7) =    38.67    
           dB_d2( 11,  8) =   -18.81    
           dB_d2( 11,  9) =    47.33    
           dB_d2( 11, 10) =   -110.2    
           dB_d2( 11, 11) =   -37.60    
           dB_d2( 11, 12) =    39.75    
           dB_d2( 11, 13) =    60.53    
           dB_d2( 11, 14) =   -126.1    
           dB_d2( 12,  1) =    3.301    
           dB_d2( 12,  2) =    36.12    
           dB_d2( 12,  3) =   -78.26    
           dB_d2( 12,  4) =    92.70    
           dB_d2( 12,  5) =   -35.81    
           dB_d2( 12,  6) =    44.88    
           dB_d2( 12,  7) =    9.250    
           dB_d2( 12,  8) =    23.46    
           dB_d2( 12,  9) =    36.15    
           dB_d2( 12, 10) =   -32.27    
           dB_d2( 12, 11) =   -36.65    
           dB_d2( 12, 12) =   -13.60    
           dB_d2( 12, 13) =    23.94    
           dB_d2( 12, 14) =   -47.02    
           dB_d2( 13,  1) =    16.03    
           dB_d2( 13,  2) =   -2.125    
           dB_d2( 13,  3) =   -85.08    
           dB_d2( 13,  4) =    83.56    
           dB_d2( 13,  5) =   -26.88    
           dB_d2( 13,  6) =    31.33    
           dB_d2( 13,  7) =   -21.80    
           dB_d2( 13,  8) =    53.83    
           dB_d2( 13,  9) =   -5.961    
           dB_d2( 13, 10) =    65.05    
           dB_d2( 13, 11) =   -12.45    
           dB_d2( 13, 12) =   -66.18    
           dB_d2( 13, 13) =    21.92    
           dB_d2( 13, 14) =   -71.32    
           dB_d2( 14,  1) =    19.20    
           dB_d2( 14,  2) =   -43.34    
           dB_d2( 14,  3) =   -79.74    
           dB_d2( 14,  4) =    90.82    
           dB_d2( 14,  5) =   -14.87    
           dB_d2( 14,  6) =    17.50    
           dB_d2( 14,  7) =   -51.39    
           dB_d2( 14,  8) =    96.66    
           dB_d2( 14,  9) =    19.14    
           dB_d2( 14, 10) =    57.89    
           dB_d2( 14, 11) =   -8.800    
           dB_d2( 14, 12) =   -11.20    
           dB_d2( 14, 13) =    44.47    
           dB_d2( 14, 14) =   -156.9    
           dB_d2( 15,  1) =    15.32    
           dB_d2( 15,  2) =   -48.58    
           dB_d2( 15,  3) =   -107.6    
           dB_d2( 15,  4) =    82.65    
           dB_d2( 15,  5) =    32.15    
           dB_d2( 15,  6) =   -87.00    
           dB_d2( 15,  7) =   -86.88    
           dB_d2( 15,  8) =    136.7    
           dB_d2( 15,  9) =   -1.550    
           dB_d2( 15, 10) =    142.9    
           dB_d2( 15, 11) =    8.200    
           dB_d2( 15, 12) =    20.40    
           dB_d2( 15, 13) =    51.25    
           dB_d2( 15, 14) =   -165.3    
           dB_d2( 16,  1) =   -6.254    
           dB_d2( 16,  2) =    70.61    
           dB_d2( 16,  3) =   -125.2    
           dB_d2( 16,  4) =    132.9    
           dB_d2( 16,  5) =    39.99    
           dB_d2( 16,  6) =   -71.32    
           dB_d2( 16,  7) =    5.966    
           dB_d2( 16,  8) =   -56.08    
           dB_d2( 16,  9) =   -26.47    
           dB_d2( 16, 10) =    112.3    
           dB_d2( 16, 11) =    31.43    
           dB_d2( 16, 12) =   -1.856    
           dB_d2( 16, 13) =    61.13    
           dB_d2( 16, 14) =   -208.1    
           dB_d2( 17,  1) =    1.425    
           dB_d2( 17,  2) =    154.7    
           dB_d2( 17,  3) =   -105.7    
           dB_d2( 17,  4) =   -105.0    
           dB_d2( 17,  5) =    68.47    
           dB_d2( 17,  6) =   -55.82    
           dB_d2( 17,  7) =    48.46    
           dB_d2( 17,  8) =   -12.20    
           dB_d2( 17,  9) =   -34.34    
           dB_d2( 17, 10) =    58.43    
           dB_d2( 17, 11) =    36.41    
           dB_d2( 17, 12) =   -50.40    
           dB_d2( 17, 13) =    64.10    
           dB_d2( 17, 14) =   -63.47    
           dB_d2( 18,  1) =   -15.19    
           dB_d2( 18,  2) =    116.2    
           dB_d2( 18,  3) =   -101.0    
           dB_d2( 18,  4) =    16.59    
           dB_d2( 18,  5) =    37.51    
           dB_d2( 18,  6) =    4.600    
           dB_d2( 18,  7) =   -36.15    
           dB_d2( 18,  8) =    64.17    
           dB_d2( 18,  9) =   -78.29    
           dB_d2( 18, 10) =    39.80    
           dB_d2( 18, 11) =    85.39    
           dB_d2( 18, 12) =   -129.3    
           dB_d2( 18, 13) =    15.83    
           dB_d2( 18, 14) =    10.40    
           dB_d2( 19,  1) =    4.205    
           dB_d2( 19,  2) =    118.7    
           dB_d2( 19,  3) =   -106.7    
           dB_d2( 19,  4) =   -152.6    
           dB_d2( 19,  5) =    68.53    
           dB_d2( 19,  6) =   -146.5    
           dB_d2( 19,  7) =   -76.28    
           dB_d2( 19,  8) =    169.9    
           dB_d2( 19,  9) =   -51.06    
           dB_d2( 19, 10) =    179.5    
           dB_d2( 19, 11) =    79.61    
           dB_d2( 19, 12) =   -95.03    
           dB_d2( 19, 13) =   -28.60    
           dB_d2( 19, 14) =   -25.00    
           dB_d2( 20,  1) =   -13.13    
           dB_d2( 20,  2) =   -127.2    
           dB_d2( 20,  3) =   -33.14    
           dB_d2( 20,  4) =    74.98    
           dB_d2( 20,  5) =    41.62    
           dB_d2( 20,  6) =   -102.9    
           dB_d2( 20,  7) =   -96.96    
           dB_d2( 20,  8) =   -39.95    
           dB_d2( 20,  9) =    4.517    
           dB_d2( 20, 10) =    27.95    
           dB_d2( 20, 11) =    78.18    
           dB_d2( 20, 12) =    133.7    
           dB_d2( 20, 13) =   -58.00    
           dB_d2( 20, 14) =    175.4    
           dB_d2( 21,  1) =    38.72    
           dB_d2( 21,  2) =   -193.9    
           dB_d2( 21,  3) =   -24.78    
           dB_d2( 21,  4) =    107.5    
           dB_d2( 21,  5) =    77.29    
           dB_d2( 21,  6) =    52.17    
           dB_d2( 21,  7) =   -103.0    
           dB_d2( 21,  8) =   -37.05    
           dB_d2( 21,  9) =   -22.36    
           dB_d2( 21, 10) =   -57.03    
           dB_d2( 21, 11) =    73.92    
           dB_d2( 21, 12) =    19.25    
           dB_d2( 21, 13) =    31.98    
           dB_d2( 21, 14) =    120.7    
           dB_d2( 22,  1) =    74.35    
           dB_d2( 22,  2) =   -77.65    
           dB_d2( 22,  3) =   -67.63    
           dB_d2( 22,  4) =    218.9    
           dB_d2( 22,  5) =    62.77    
           dB_d2( 22,  6) =    113.0    
           dB_d2( 22,  7) =   -12.99    
           dB_d2( 22,  8) =    1.955    
           dB_d2( 22,  9) =   -37.37    
           dB_d2( 22, 10) =    37.90    
           dB_d2( 22, 11) =  -0.2550    
           dB_d2( 22, 12) =   -34.86    
           dB_d2( 22, 13) =    35.20    
           dB_d2( 22, 14) =    109.7    
           dB_d2( 23,  1) =    18.32    
           dB_d2( 23,  2) =   -11.86    
           dB_d2( 23,  3) =    6.825    
           dB_d2( 23,  4) =    20.72    
           dB_d2( 23,  5) =    10.56    
           dB_d2( 23,  6) =   -14.29    
           dB_d2( 23,  7) =   -73.54    
           dB_d2( 23,  8) =   -1.405    
           dB_d2( 23,  9) =   -3.433    
           dB_d2( 23, 10) =   -12.05    
           dB_d2( 23, 11) =    14.52    
           dB_d2( 23, 12) =    7.385    
           dB_d2( 23, 13) =   -3.960    
           dB_d2( 23, 14) =    31.30 
       
           ! Fill mlat arrays:
           min_mlat(  1) =     45.00
           max_mlat(  1) =     47.00
           min_mlat(  2) =     47.00
           max_mlat(  2) =     49.00
           min_mlat(  3) =     49.00
           max_mlat(  3) =     51.00
           min_mlat(  4) =     51.00
           max_mlat(  4) =     53.00
           min_mlat(  5) =     53.00
           max_mlat(  5) =     55.00
           min_mlat(  6) =     55.00
           max_mlat(  6) =     57.00
           min_mlat(  7) =     57.00
           max_mlat(  7) =     59.00
           min_mlat(  8) =     59.00
           max_mlat(  8) =     61.00
           min_mlat(  9) =     61.00
           max_mlat(  9) =     63.00
           min_mlat( 10) =     63.00
           max_mlat( 10) =     65.00
           min_mlat( 11) =     65.00
           max_mlat( 11) =     67.00
           min_mlat( 12) =     67.00
           max_mlat( 12) =     69.00
           min_mlat( 13) =     69.00
           max_mlat( 13) =     71.00
           min_mlat( 14) =     71.00
           max_mlat( 14) =     73.00
           min_mlat( 15) =     73.00
           max_mlat( 15) =     75.00
           min_mlat( 16) =     75.00
           max_mlat( 16) =     77.00
           min_mlat( 17) =     77.00
           max_mlat( 17) =     79.00
           min_mlat( 18) =     79.00
           max_mlat( 18) =     81.00
           min_mlat( 19) =     81.00
           max_mlat( 19) =     83.00
           min_mlat( 20) =     83.00
           max_mlat( 20) =     85.00
           min_mlat( 21) =     85.00
           max_mlat( 21) =     87.00
           min_mlat( 22) =     87.00
           max_mlat( 22) =     89.00
           min_mlat( 23) =     89.00
           max_mlat( 23) =     91.00
       
         end subroutine B_NBZ_fill_coeff_arrays
       
       end module B_NBZ_coeffmerge_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       subroutine B_NBZ_calc_model_value(mlat_glob,mlt,sinT,imf_angle,
     |       value_Ed1,value_Ed2,value,value_found,FAC_value)
       
         ! This version is for public distribution, along with the module
         ! file coeff_mod.f90.
         !
         ! Within the /home/hensel/fortran/fit_data/public subdirectory,
         ! this program is called by call_de2_model.
       
         use common_model_mod, only: MMX,MXNMX,ITHMX,mlat_resolution
       
         use B_NBZ_model_mod, only: 
     |      pi,find_Btrans_function,find_imf_angle_function,
     |      find_sint_function,AO =>B_NBZ_AO,Q =>B_NBZ_Q,DQ=>B_NBZ_DQ
       
         use B_NBZ_model_mod, only:  B_NBZ_number_of_MLT_coeffs,
     |      B_NBZ_number_of_Btrans_coeffs, 
     |      B_NBZ_number_of_imf_angle_coeffs,
     |      B_NBZ_number_of_sinT_coeffs,value_dummy,
     |      fit_by_Btrans => B_NBZ_fit_by_Btrans,fit_by_imf_angle=> 
     |      B_NBZ_fit_by_imf_angle,fit_by_sinT => B_NBZ_fit_by_sinT,
     |      FAC_AO => NBZ_FAC_AO
       
         use qgen_mod, only: RI,ST,IBM,IEM,mcoef,dth,ns
       
         implicit none
       
         real, intent(in) :: mlat_glob,mlt,sinT,imf_angle
         real, intent(out) :: value_Ed1,value_Ed2,value
         logical, intent(out) :: value_found
         real, intent(out) :: FAC_value
       
         integer :: ith,m,k,ix,mm,iBtrans,iIMF_angle,isinT,fit_count, 
     |         iMLT
         integer :: m2iMLT(B_NBZ_number_of_MLT_coeffs)
         data m2iMLT /0,-1,1,-2,2,-3,3/ 
         real :: mlat,poten,etheta,ephi,
     |       x,xm1,xp1,xp2,x0,sp,cp,sts,rst,qsoso,qsost,
     |       Btrans_function,imf_angle_function,sint_function,sinIm, 
     |       FAC_poten,get_sinIm
         real :: QS(-MMX:MMX,MXNMX),DQS(-MMX:MMX,MXNMX),F(-MMX:MMX)
         real :: fac_Y,Ymlat
       
         ! General:
         logical, parameter :: debug = .false.
       
         !****************************************************************
         ! Interpolate to get coefficients for specified mlat:
         !****************************************************************
       
         value_found = .false.
         value_dummy = .true.
       
         ! Use absolute value of mlat:
         mlat = abs(mlat_glob)
       
         ! C Construct electric potential and electric field components 
         !C   (etheta equatorward, ephi eastward) 
         !C   on 2x10 degree grid, 0-40 colatitude
         !C Colatitudes (degrees) are in CLATD(0:ITHPLT)
         !C MLTs (hours) are in MLT(0:24)
         poten  = 0.
         etheta = 0.
         ephi   = 0.
         FAC_poten = 0.
         !C Calculate multipliers for cubic interpolation between latitudes where
         !C  Q, DQ are available
         X = (90.-mlat)*pi/180./DTH
         ITH = X
         ITH = MAX0(1,MIN0(ITHMX-2,ITH))
         X   = X - ITH
         XM1 = X*(-2. + X*(3. - X))/6.
         X0  = 1. + X*(-.5 + X*(-1. + .5*X))
         XP1 = X*(1. + X*(.5 - .5*X))
         XP2 = X*(-1. + X*X)/6.
       
         sp = sin(MLT*pi/12.)
         cp = COS(MLT*pi/12.)
         call FCMP_single(MMX,CP,SP,F)
       
         sts =  XM1* ST(ITH-1) + X0*ST(ITH)+ XP1*ST(ITH+1)+XP2*ST(ITH+2)   
       
         rst = RI*STS                                
         do iMLT = 1,B_NBZ_number_of_MLT_coeffs 
            m = m2iMLT(iMLT) 						     
            MM = IABS(M)
       
            Btrans_loop1:do iBtrans=1,1+(fit_by_Btrans
     |            *B_NBZ_number_of_Btrans_coeffs-fit_by_Btrans)	 
               imf_angle_loop1:do iimf_angle=1,1+(fit_by_imf_angle
     |               *B_NBZ_number_of_imf_angle_coeffs-fit_by_imf_angle)
     
                  imf_angle_function = 
     |                     find_imf_angle_function(iimf_angle,imf_angle) 
                  sint_loop1:do isint=1,1+(fit_by_sint
     |                  *B_NBZ_number_of_sint_coeffs - fit_by_sint)
                     sint_function = find_sint_function (isint,sint)   
       
       	            do K=1,MXNMX					     
                     IX = K - 1 					     
                     QS(M,K) = XM1*Q(2*IX+NS(MM),ITH-1)+
     |                  X0*Q(2*IX+NS(MM),ITH)+ 
     |                	XP1*Q(2*IX+NS(MM),ITH+1)+
     |                  XP2*Q(2*IX+NS(MM),ITH+2)
                     DQS(M,K) = XM1*DQ(2*IX+NS(MM),ITH-1)+
     |                  X0*DQ(2*IX+NS(MM),ITH)+ 
     |                	XP1*DQ(2*IX+NS(MM),ITH+1)+
     |                  XP2*DQ(2*IX+NS(MM),ITH+2)   
                     IF (ITH.GT.0) THEN
                     	qsost = QS(M,K)/STS
                     ELSE
                     	qsost = DQS(M,K)
                     ENDIF
                     poten=poten+AO(M,K,iBtrans,iimf_angle,isint)*
     |                  QS(M,K)*F(M)*
     |                  imf_angle_function*sint_function 
                     etheta = etheta + 
     |                  AO(M,K,iBtrans,iimf_angle,isint)*
     |                  DQS(M,K)*F(M)/RI*
     |                  imf_angle_function*sint_function 
                     ephi=ephi- AO(M,K,iBtrans,iimf_angle,isint)*
     |                  qsost*M*F(-M)/RI*
     |                  imf_angle_function*sint_function 
                      FAC_poten  = FAC_poten  +  
     |                  FAC_AO(M,K,iBtrans,iimf_angle,isint)* 
     |                  QS(M,K)*F(M)*
     |                  imf_angle_function*sint_function 
                     enddo
                  end do sint_loop1
               end do imf_angle_loop1
            end do Btrans_loop1
         enddo
       
         sinIm =  get_sinIm (mlat)
         !  (etheta equatorward, ephi eastward)
         value     =  poten
         FAC_value = FAC_poten
         value_Ed1 =  etheta/abs(sin(mlat*pi/180.))  ! dBd1
         value_Ed2 =  ephi/sqrt(0.25+0.75*      ! dBd2
     |         sin(mlat*pi/180.)**2)
         value_found = .true.
         value_dummy = .false.
       
         if (debug) then
            write(*,'(" mlat=",F9.2,", MLT=",F9.2,":  ",G15.6)') 
     |            mlat,MLT*180./24.,value
         end if
       
       
       end subroutine B_NBZ_calc_model_value
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       
       subroutine B_NBZ_set_coeff_array
       
         use common_model_mod, only: nmlat,average_mlat
         use B_NBZ_coeffmerge_mod, only: dB_d1,dB_d2,
     |       B_NBZ_fill_coeff_arrays
         use B_NBZ_model_mod, only: nfit=>B_NBZ_nfit,b_Ed1=>B_NBZ_b_Ed1,
     |      b_Ed2 => B_NBZ_b_Ed2,b_interp_Ed1,
     |       b_interp_Ed2,b_90_Ed1,b_90_Ed2,x,imlat,ifit,
     |      B_NBZ_number_of_MLT_coeffs,
     |      B_NBZ_number_of_Btrans_coeffs,
     |      B_NBZ_number_of_imf_angle_coeffs,
     |      B_NBZ_number_of_sinT_coeffs,
     |      fit_by_MLT,pi,
     |      fit_by_Btrans => B_NBZ_fit_by_Btrans, 
     |      fit_by_imf_angle => B_NBZ_fit_by_imf_angle, 
     |      fit_by_sinT => B_NBZ_fit_by_sinT	   
       
       
         integer :: binary_fit_by,check_dim
         integer :: iMLT,iBtrans,iimf_angle,isinT,fit_count
         real :: fac_ed1, fac_ed2, colat
         real :: fac_Y,Y
       
         !****************************************************************
         ! Fill parameter arrays:
         !****************************************************************
         call B_NBZ_fill_coeff_arrays
       
         !****************************************************************
         ! Determine array dimensions.  Arrays for data type 'potential'
         ! are the same size as arrays for other data types.
         !****************************************************************
       
         nmlat = size(dB_d1,1)
         nfit  = size(dB_d1,2)
       
         ! check dimensions
         check_dim = B_NBZ_number_of_MLT_coeffs*
     |      B_NBZ_number_of_Btrans_coeffs*
     |      B_NBZ_number_of_imf_angle_coeffs*
     |      B_NBZ_number_of_sinT_coeffs
         if(check_dim /= nfit) then
            write(6,'("choose_coeff_array: wrong dimensions ", 2i4)') 
     |             nfit,check_dim
            stop
         endif
         !****************************************************************
         ! Allocate arrays into which coefficients will be transferred:
         !****************************************************************
         allocate(b_Ed1(nmlat,nfit),b_Ed2(nmlat,nfit),
     |         b_interp_Ed1(nfit),
     |         b_interp_Ed2(nfit),x(nfit),
     |         b_90_Ed1(nfit),b_90_Ed2(nfit))
       
         !****************************************************************
         ! Transfer values:
         !****************************************************************
       
         do imlat=1,nmlat
            do ifit=1,nfit
              if(abs(dB_d1(imlat,ifit))>=0.5000E+17)dB_d1(imlat,ifit)=0.
              if(abs(dB_d2(imlat,ifit))>=0.5000E+17)dB_d2(imlat,ifit)=0.
       
              b_Ed1(imlat,ifit) = abs(sin(average_mlat(imlat)*pi/180.))*
     |              dB_d2(imlat,ifit)/2**0.5
              b_Ed2(imlat,ifit) = sqrt(0.25+0.75*
     |              sin(average_mlat(imlat)*pi/180.)**2)*
     |              dB_d1(imlat,ifit)/2**0.5
            end do
         end do
       
       end subroutine B_NBZ_set_coeff_array
                                                                                                                              

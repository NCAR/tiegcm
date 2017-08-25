       module SDEd_NBZ_model_mod
       
         use common_model_mod, only: pi,fit_by_MLT
       
         implicit none   
       
         logical :: first_NBz 
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Data type
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         integer :: itype
         integer,parameter :: number_of_datatypes=6
         character (len=80) :: data_type(number_of_datatypes)
       
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Arrays for fit coefficients:
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! Number of fit coefficients for each mlat:
         integer :: nfit
         ! Number of mlat intervals:
         integer :: nmlat
         ! Fit coefficients (dimension nmlat x nfit):
         real, allocatable, dimension(:,:) :: b1,b2 
         ! Interpolated fit coefficients (dimension nfit):
         real, allocatable, dimension(:) :: b_interp1,b_interp2
         ! Fitting functions (dimension nfit):
         real, allocatable, dimension(:) :: x
         ! Fit coefficients averaged over the pole (dimension nfit):
         real, allocatable, dimension(:) :: b_90 
         ! Number of fit functions used by fit_data:
         integer, parameter :: number_of_MLT_coeffs = 5,
     |         number_of_Btrans_coeffs = 1,
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
         !get_y.f90:
         real :: value
         logical :: value_found
       
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
     |          result (Btrans_fxn_out)
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
       
       end module SDEd_NBZ_model_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       module SDEd_NBZ_coeff_mod
       
         ! This module file created 2007/08/27.
         ! with the following parameters:
         !   for 81252-83047 period 
         !   fit_by_Btrans = 0
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
         real,  dimension(  9, 10) :: Ed1_Var_coeffs
         real,  dimension(  9, 10) :: Ed2_Var_coeffs
         real,  dimension(  9, 10) :: Ed1_Var_DE2_coeffs
         real,  dimension(  9, 10) :: Ed2_Var_DE2_coeffs
         real,  dimension(  9, 10) :: Ed1_Var_mod_coeffs
         real,  dimension(  9, 10) :: Ed2_Var_mod_coeffs
         real,  dimension(  9) :: min_mlat,max_mlat,average_mlat
         integer, parameter :: 
     |          fit_by_Btrans = 0, 
     |          fit_by_imf_angle = 0, 
     |          fit_by_sinT = 1, 
     |          mlat_resolution = 5
       
       contains
       
         subroutine SDEd_NBZ_fill_coeff_arrays
       
           ! Fill Ed1_Var arrays:
           Ed1_Var_coeffs(  1,  1) =    3.592    
           Ed1_Var_coeffs(  1,  2) =   -3.398    
           Ed1_Var_coeffs(  1,  3) =    1.738    
           Ed1_Var_coeffs(  1,  4) =    3.198    
           Ed1_Var_coeffs(  1,  5) =   0.1064    
           Ed1_Var_coeffs(  1,  6) =    2.028    
           Ed1_Var_coeffs(  1,  7) =  -0.6679    
           Ed1_Var_coeffs(  1,  8) =    2.439    
           Ed1_Var_coeffs(  1,  9) =  -0.7194    
           Ed1_Var_coeffs(  1, 10) =    2.767    
           Ed1_Var_coeffs(  2,  1) =    13.48    
           Ed1_Var_coeffs(  2,  2) =   -35.26    
           Ed1_Var_coeffs(  2,  3) =    12.17    
           Ed1_Var_coeffs(  2,  4) =   -29.03    
           Ed1_Var_coeffs(  2,  5) =    1.388    
           Ed1_Var_coeffs(  2,  6) =   0.7410    
           Ed1_Var_coeffs(  2,  7) =    3.378    
           Ed1_Var_coeffs(  2,  8) =   -27.92    
           Ed1_Var_coeffs(  2,  9) =    3.391    
           Ed1_Var_coeffs(  2, 10) =   -9.402    
           Ed1_Var_coeffs(  3,  1) =    13.50    
           Ed1_Var_coeffs(  3,  2) =   -28.25    
           Ed1_Var_coeffs(  3,  3) =    5.939    
           Ed1_Var_coeffs(  3,  4) =   -8.140    
           Ed1_Var_coeffs(  3,  5) =   -3.449    
           Ed1_Var_coeffs(  3,  6) =    12.13    
           Ed1_Var_coeffs(  3,  7) =   -1.363    
           Ed1_Var_coeffs(  3,  8) =   -8.903    
           Ed1_Var_coeffs(  3,  9) =   0.3157E-01
           Ed1_Var_coeffs(  3, 10) =   -5.771    
           Ed1_Var_coeffs(  4,  1) =    21.79    
           Ed1_Var_coeffs(  4,  2) =   -40.64    
           Ed1_Var_coeffs(  4,  3) =    8.159    
           Ed1_Var_coeffs(  4,  4) =   -4.391    
           Ed1_Var_coeffs(  4,  5) =   -9.095    
           Ed1_Var_coeffs(  4,  6) =    21.30    
           Ed1_Var_coeffs(  4,  7) =   -4.443    
           Ed1_Var_coeffs(  4,  8) =   -11.67    
           Ed1_Var_coeffs(  4,  9) =   0.8299    
           Ed1_Var_coeffs(  4, 10) =   -17.51    
           Ed1_Var_coeffs(  5,  1) =    54.79    
           Ed1_Var_coeffs(  5,  2) =   -72.54    
           Ed1_Var_coeffs(  5,  3) =    47.69    
           Ed1_Var_coeffs(  5,  4) =   -44.27    
           Ed1_Var_coeffs(  5,  5) =   -8.085    
           Ed1_Var_coeffs(  5,  6) =    33.59    
           Ed1_Var_coeffs(  5,  7) =    13.06    
           Ed1_Var_coeffs(  5,  8) =   -40.33    
           Ed1_Var_coeffs(  5,  9) =   -5.748    
           Ed1_Var_coeffs(  5, 10) =   -36.73    
           Ed1_Var_coeffs(  6,  1) =    161.9    
           Ed1_Var_coeffs(  6,  2) =   -254.8    
           Ed1_Var_coeffs(  6,  3) =    137.4    
           Ed1_Var_coeffs(  6,  4) =   -230.6    
           Ed1_Var_coeffs(  6,  5) =    26.95    
           Ed1_Var_coeffs(  6,  6) =    21.56    
           Ed1_Var_coeffs(  6,  7) =    24.98    
           Ed1_Var_coeffs(  6,  8) =   -73.65    
           Ed1_Var_coeffs(  6,  9) =    25.81    
           Ed1_Var_coeffs(  6, 10) =   -89.60    
           Ed1_Var_coeffs(  7,  1) =    357.4    
           Ed1_Var_coeffs(  7,  2) =   -542.0    
           Ed1_Var_coeffs(  7,  3) =    97.83    
           Ed1_Var_coeffs(  7,  4) =   -544.6    
           Ed1_Var_coeffs(  7,  5) =    102.8    
           Ed1_Var_coeffs(  7,  6) =   -112.7    
           Ed1_Var_coeffs(  7,  7) =    1.254    
           Ed1_Var_coeffs(  7,  8) =    294.3    
           Ed1_Var_coeffs(  7,  9) =  -0.8758    
           Ed1_Var_coeffs(  7, 10) =   -182.4    
           Ed1_Var_coeffs(  8,  1) =    479.3    
           Ed1_Var_coeffs(  8,  2) =   -590.3    
           Ed1_Var_coeffs(  8,  3) =   -23.22    
           Ed1_Var_coeffs(  8,  4) =   -640.9    
           Ed1_Var_coeffs(  8,  5) =    49.43    
           Ed1_Var_coeffs(  8,  6) =   -132.4    
           Ed1_Var_coeffs(  8,  7) =    213.4    
           Ed1_Var_coeffs(  8,  8) =    15.90    
           Ed1_Var_coeffs(  8,  9) =   -97.99    
           Ed1_Var_coeffs(  8, 10) =    316.1    
           Ed1_Var_coeffs(  9,  1) =    305.3    
           Ed1_Var_coeffs(  9,  2) =   -434.7    
           Ed1_Var_coeffs(  9,  3) =    37.85    
           Ed1_Var_coeffs(  9,  4) =   -461.2    
           Ed1_Var_coeffs(  9,  5) =   -69.25    
           Ed1_Var_coeffs(  9,  6) =    289.5    
           Ed1_Var_coeffs(  9,  7) =    125.7    
           Ed1_Var_coeffs(  9,  8) =   -264.0    
           Ed1_Var_coeffs(  9,  9) =   -23.70    
           Ed1_Var_coeffs(  9, 10) =    77.26    
           ! Fill Ed2_Var arrays:
           Ed2_Var_coeffs(  1,  1) =    6.710    
           Ed2_Var_coeffs(  1,  2) =   -4.679    
           Ed2_Var_coeffs(  1,  3) =    3.303    
           Ed2_Var_coeffs(  1,  4) =    7.527    
           Ed2_Var_coeffs(  1,  5) =   0.7549    
           Ed2_Var_coeffs(  1,  6) =   0.6714    
           Ed2_Var_coeffs(  1,  7) =   -2.622    
           Ed2_Var_coeffs(  1,  8) =    3.205    
           Ed2_Var_coeffs(  1,  9) =   -1.612    
           Ed2_Var_coeffs(  1, 10) =    9.623    
           Ed2_Var_coeffs(  2,  1) =    9.716    
           Ed2_Var_coeffs(  2,  2) =   -11.47    
           Ed2_Var_coeffs(  2,  3) =    2.915    
           Ed2_Var_coeffs(  2,  4) =   0.5536    
           Ed2_Var_coeffs(  2,  5) =   -8.180    
           Ed2_Var_coeffs(  2,  6) =    6.388    
           Ed2_Var_coeffs(  2,  7) =   -6.719    
           Ed2_Var_coeffs(  2,  8) =  -0.2808    
           Ed2_Var_coeffs(  2,  9) =   -1.934    
           Ed2_Var_coeffs(  2, 10) =    6.205    
           Ed2_Var_coeffs(  3,  1) =    25.19    
           Ed2_Var_coeffs(  3,  2) =   -28.63    
           Ed2_Var_coeffs(  3,  3) =    12.28    
           Ed2_Var_coeffs(  3,  4) =   -14.65    
           Ed2_Var_coeffs(  3,  5) =   -31.78    
           Ed2_Var_coeffs(  3,  6) =    15.19    
           Ed2_Var_coeffs(  3,  7) =   -29.56    
           Ed2_Var_coeffs(  3,  8) =   -6.437    
           Ed2_Var_coeffs(  3,  9) =   -12.08    
           Ed2_Var_coeffs(  3, 10) =    10.23    
           Ed2_Var_coeffs(  4,  1) =    82.10    
           Ed2_Var_coeffs(  4,  2) =   -90.94    
           Ed2_Var_coeffs(  4,  3) =    38.37    
           Ed2_Var_coeffs(  4,  4) =   -38.52    
           Ed2_Var_coeffs(  4,  5) =   -93.17    
           Ed2_Var_coeffs(  4,  6) =    216.2    
           Ed2_Var_coeffs(  4,  7) =   -87.78    
           Ed2_Var_coeffs(  4,  8) =    182.9    
           Ed2_Var_coeffs(  4,  9) =   -27.01    
           Ed2_Var_coeffs(  4, 10) =    80.53    
           Ed2_Var_coeffs(  5,  1) =    191.7    
           Ed2_Var_coeffs(  5,  2) =   -185.8    
           Ed2_Var_coeffs(  5,  3) =    198.0    
           Ed2_Var_coeffs(  5,  4) =   -271.5    
           Ed2_Var_coeffs(  5,  5) =   -23.87    
           Ed2_Var_coeffs(  5,  6) =    38.63    
           Ed2_Var_coeffs(  5,  7) =    28.76    
           Ed2_Var_coeffs(  5,  8) =   -31.06    
           Ed2_Var_coeffs(  5,  9) =    6.949    
           Ed2_Var_coeffs(  5, 10) =   -19.91    
           Ed2_Var_coeffs(  6,  1) =    451.3    
           Ed2_Var_coeffs(  6,  2) =   -873.7    
           Ed2_Var_coeffs(  6,  3) =    315.7    
           Ed2_Var_coeffs(  6,  4) =   -427.2    
           Ed2_Var_coeffs(  6,  5) =    163.3    
           Ed2_Var_coeffs(  6,  6) =   -314.1    
           Ed2_Var_coeffs(  6,  7) =    39.45    
           Ed2_Var_coeffs(  6,  8) =   -280.0    
           Ed2_Var_coeffs(  6,  9) =    228.9    
           Ed2_Var_coeffs(  6, 10) =   -478.9    
           Ed2_Var_coeffs(  7,  1) =    586.8    
           Ed2_Var_coeffs(  7,  2) =   -847.8    
           Ed2_Var_coeffs(  7,  3) =   -327.7    
           Ed2_Var_coeffs(  7,  4) =    330.6    
           Ed2_Var_coeffs(  7,  5) =    143.4    
           Ed2_Var_coeffs(  7,  6) =   -363.0    
           Ed2_Var_coeffs(  7,  7) =   -75.97    
           Ed2_Var_coeffs(  7,  8) =    444.8    
           Ed2_Var_coeffs(  7,  9) =   -29.31    
           Ed2_Var_coeffs(  7, 10) =    132.3    
           Ed2_Var_coeffs(  8,  1) =    563.3    
           Ed2_Var_coeffs(  8,  2) =   -544.2    
           Ed2_Var_coeffs(  8,  3) =   -272.2    
           Ed2_Var_coeffs(  8,  4) =   -482.6    
           Ed2_Var_coeffs(  8,  5) =    70.00    
           Ed2_Var_coeffs(  8,  6) =    19.94    
           Ed2_Var_coeffs(  8,  7) =   -91.67    
           Ed2_Var_coeffs(  8,  8) =    739.9    
           Ed2_Var_coeffs(  8,  9) =    10.12    
           Ed2_Var_coeffs(  8, 10) =   -561.3    
           Ed2_Var_coeffs(  9,  1) =    318.5    
           Ed2_Var_coeffs(  9,  2) =   -512.8    
           Ed2_Var_coeffs(  9,  3) =    72.50    
           Ed2_Var_coeffs(  9,  4) =   -443.1    
           Ed2_Var_coeffs(  9,  5) =   -49.64    
           Ed2_Var_coeffs(  9,  6) =    110.0    
           Ed2_Var_coeffs(  9,  7) =   -104.6    
           Ed2_Var_coeffs(  9,  8) =    567.5    
           Ed2_Var_coeffs(  9,  9) =    36.15    
           Ed2_Var_coeffs(  9, 10) =   -345.3    
           ! Fill Ed1_Var_DE2 arrays:
           Ed1_Var_DE2_coeffs(  1,  1) =    3.198    
           Ed1_Var_DE2_coeffs(  1,  2) =   -3.485    
           Ed1_Var_DE2_coeffs(  1,  3) =   0.3380    
           Ed1_Var_DE2_coeffs(  1,  4) =    4.327    
           Ed1_Var_DE2_coeffs(  1,  5) =   0.2082    
           Ed1_Var_DE2_coeffs(  1,  6) =    3.881    
           Ed1_Var_DE2_coeffs(  1,  7) =  -0.9362    
           Ed1_Var_DE2_coeffs(  1,  8) =    4.171    
           Ed1_Var_DE2_coeffs(  1,  9) =  -0.1058    
           Ed1_Var_DE2_coeffs(  1, 10) =    2.514    
           Ed1_Var_DE2_coeffs(  2,  1) =    13.32    
           Ed1_Var_DE2_coeffs(  2,  2) =   -34.15    
           Ed1_Var_DE2_coeffs(  2,  3) =    10.04    
           Ed1_Var_DE2_coeffs(  2,  4) =   -23.75    
           Ed1_Var_DE2_coeffs(  2,  5) =    1.766    
           Ed1_Var_DE2_coeffs(  2,  6) =   0.8761    
           Ed1_Var_DE2_coeffs(  2,  7) =    2.924    
           Ed1_Var_DE2_coeffs(  2,  8) =   -24.94    
           Ed1_Var_DE2_coeffs(  2,  9) =    5.132    
           Ed1_Var_DE2_coeffs(  2, 10) =   -12.97    
           Ed1_Var_DE2_coeffs(  3,  1) =    14.28    
           Ed1_Var_DE2_coeffs(  3,  2) =   -27.20    
           Ed1_Var_DE2_coeffs(  3,  3) =    3.614    
           Ed1_Var_DE2_coeffs(  3,  4) =   -2.302    
           Ed1_Var_DE2_coeffs(  3,  5) =   -3.351    
           Ed1_Var_DE2_coeffs(  3,  6) =    12.01    
           Ed1_Var_DE2_coeffs(  3,  7) =   -1.644    
           Ed1_Var_DE2_coeffs(  3,  8) =   -5.703    
           Ed1_Var_DE2_coeffs(  3,  9) =    1.897    
           Ed1_Var_DE2_coeffs(  3, 10) =   -11.84    
           Ed1_Var_DE2_coeffs(  4,  1) =    23.48    
           Ed1_Var_DE2_coeffs(  4,  2) =   -37.91    
           Ed1_Var_DE2_coeffs(  4,  3) =    7.996    
           Ed1_Var_DE2_coeffs(  4,  4) =   -1.927    
           Ed1_Var_DE2_coeffs(  4,  5) =   -10.45    
           Ed1_Var_DE2_coeffs(  4,  6) =    20.67    
           Ed1_Var_DE2_coeffs(  4,  7) =   -4.760    
           Ed1_Var_DE2_coeffs(  4,  8) =   -10.62    
           Ed1_Var_DE2_coeffs(  4,  9) =    1.160    
           Ed1_Var_DE2_coeffs(  4, 10) =   -20.81    
           Ed1_Var_DE2_coeffs(  5,  1) =    56.71    
           Ed1_Var_DE2_coeffs(  5,  2) =   -71.64    
           Ed1_Var_DE2_coeffs(  5,  3) =    51.56    
           Ed1_Var_DE2_coeffs(  5,  4) =   -61.28    
           Ed1_Var_DE2_coeffs(  5,  5) =   -7.216    
           Ed1_Var_DE2_coeffs(  5,  6) =    19.99    
           Ed1_Var_DE2_coeffs(  5,  7) =    15.59    
           Ed1_Var_DE2_coeffs(  5,  8) =   -60.84    
           Ed1_Var_DE2_coeffs(  5,  9) =   -5.715    
           Ed1_Var_DE2_coeffs(  5, 10) =   -33.31    
           Ed1_Var_DE2_coeffs(  6,  1) =    164.1    
           Ed1_Var_DE2_coeffs(  6,  2) =   -245.6    
           Ed1_Var_DE2_coeffs(  6,  3) =    145.2    
           Ed1_Var_DE2_coeffs(  6,  4) =   -250.6    
           Ed1_Var_DE2_coeffs(  6,  5) =    16.15    
           Ed1_Var_DE2_coeffs(  6,  6) =    12.46    
           Ed1_Var_DE2_coeffs(  6,  7) =    27.94    
           Ed1_Var_DE2_coeffs(  6,  8) =   -85.02    
           Ed1_Var_DE2_coeffs(  6,  9) =    19.35    
           Ed1_Var_DE2_coeffs(  6, 10) =   -73.26    
           Ed1_Var_DE2_coeffs(  7,  1) =    380.4    
           Ed1_Var_DE2_coeffs(  7,  2) =   -549.9    
           Ed1_Var_DE2_coeffs(  7,  3) =    58.95    
           Ed1_Var_DE2_coeffs(  7,  4) =   -495.9    
           Ed1_Var_DE2_coeffs(  7,  5) =    104.9    
           Ed1_Var_DE2_coeffs(  7,  6) =   -68.63    
           Ed1_Var_DE2_coeffs(  7,  7) =   -7.350    
           Ed1_Var_DE2_coeffs(  7,  8) =    331.3    
           Ed1_Var_DE2_coeffs(  7,  9) =   -11.01    
           Ed1_Var_DE2_coeffs(  7, 10) =   -217.3    
           Ed1_Var_DE2_coeffs(  8,  1) =    572.0    
           Ed1_Var_DE2_coeffs(  8,  2) =   -619.4    
           Ed1_Var_DE2_coeffs(  8,  3) =   -136.0    
           Ed1_Var_DE2_coeffs(  8,  4) =   -628.3    
           Ed1_Var_DE2_coeffs(  8,  5) =    42.30    
           Ed1_Var_DE2_coeffs(  8,  6) =   -204.1    
           Ed1_Var_DE2_coeffs(  8,  7) =    239.2    
           Ed1_Var_DE2_coeffs(  8,  8) =    65.71    
           Ed1_Var_DE2_coeffs(  8,  9) =   -107.7    
           Ed1_Var_DE2_coeffs(  8, 10) =    435.5    
           Ed1_Var_DE2_coeffs(  9,  1) =    333.6    
           Ed1_Var_DE2_coeffs(  9,  2) =   -390.8    
           Ed1_Var_DE2_coeffs(  9,  3) =    2.532    
           Ed1_Var_DE2_coeffs(  9,  4) =   -515.4    
           Ed1_Var_DE2_coeffs(  9,  5) =   -88.14    
           Ed1_Var_DE2_coeffs(  9,  6) =    263.0    
           Ed1_Var_DE2_coeffs(  9,  7) =    139.6    
           Ed1_Var_DE2_coeffs(  9,  8) =   -220.3    
           Ed1_Var_DE2_coeffs(  9,  9) =   -18.88    
           Ed1_Var_DE2_coeffs(  9, 10) =    66.78    
           ! Fill Ed2_Var_DE2 arrays:
           Ed2_Var_DE2_coeffs(  1,  1) =    6.695    
           Ed2_Var_DE2_coeffs(  1,  2) =   -3.889    
           Ed2_Var_DE2_coeffs(  1,  3) =    3.860    
           Ed2_Var_DE2_coeffs(  1,  4) =    10.04    
           Ed2_Var_DE2_coeffs(  1,  5) =    2.907    
           Ed2_Var_DE2_coeffs(  1,  6) =   -1.698    
           Ed2_Var_DE2_coeffs(  1,  7) =   -1.622    
           Ed2_Var_DE2_coeffs(  1,  8) =    2.474    
           Ed2_Var_DE2_coeffs(  1,  9) =   -1.393    
           Ed2_Var_DE2_coeffs(  1, 10) =    7.351    
           Ed2_Var_DE2_coeffs(  2,  1) =    9.138    
           Ed2_Var_DE2_coeffs(  2,  2) =   -10.01    
           Ed2_Var_DE2_coeffs(  2,  3) =    3.842    
           Ed2_Var_DE2_coeffs(  2,  4) =    4.749    
           Ed2_Var_DE2_coeffs(  2,  5) =   -4.551    
           Ed2_Var_DE2_coeffs(  2,  6) =    2.376    
           Ed2_Var_DE2_coeffs(  2,  7) =   -4.108    
           Ed2_Var_DE2_coeffs(  2,  8) =   -2.014    
           Ed2_Var_DE2_coeffs(  2,  9) =   -1.647    
           Ed2_Var_DE2_coeffs(  2, 10) =    2.678    
           Ed2_Var_DE2_coeffs(  3,  1) =    29.75    
           Ed2_Var_DE2_coeffs(  3,  2) =   -34.00    
           Ed2_Var_DE2_coeffs(  3,  3) =    14.63    
           Ed2_Var_DE2_coeffs(  3,  4) =   -19.63    
           Ed2_Var_DE2_coeffs(  3,  5) =   -39.84    
           Ed2_Var_DE2_coeffs(  3,  6) =    15.56    
           Ed2_Var_DE2_coeffs(  3,  7) =   -35.49    
           Ed2_Var_DE2_coeffs(  3,  8) =   -13.19    
           Ed2_Var_DE2_coeffs(  3,  9) =   -14.68    
           Ed2_Var_DE2_coeffs(  3, 10) =    19.21    
           Ed2_Var_DE2_coeffs(  4,  1) =    102.0    
           Ed2_Var_DE2_coeffs(  4,  2) =   -115.9    
           Ed2_Var_DE2_coeffs(  4,  3) =    53.14    
           Ed2_Var_DE2_coeffs(  4,  4) =   -74.49    
           Ed2_Var_DE2_coeffs(  4,  5) =   -128.5    
           Ed2_Var_DE2_coeffs(  4,  6) =    265.1    
           Ed2_Var_DE2_coeffs(  4,  7) =   -108.6    
           Ed2_Var_DE2_coeffs(  4,  8) =    201.8    
           Ed2_Var_DE2_coeffs(  4,  9) =   -43.17    
           Ed2_Var_DE2_coeffs(  4, 10) =    128.7    
           Ed2_Var_DE2_coeffs(  5,  1) =    241.9    
           Ed2_Var_DE2_coeffs(  5,  2) =   -221.0    
           Ed2_Var_DE2_coeffs(  5,  3) =    248.5    
           Ed2_Var_DE2_coeffs(  5,  4) =   -369.8    
           Ed2_Var_DE2_coeffs(  5,  5) =   -76.01    
           Ed2_Var_DE2_coeffs(  5,  6) =    1.767    
           Ed2_Var_DE2_coeffs(  5,  7) =    18.85    
           Ed2_Var_DE2_coeffs(  5,  8) =   -127.9    
           Ed2_Var_DE2_coeffs(  5,  9) =   -18.32    
           Ed2_Var_DE2_coeffs(  5, 10) =    33.25    
           Ed2_Var_DE2_coeffs(  6,  1) =    519.8    
           Ed2_Var_DE2_coeffs(  6,  2) =   -882.1    
           Ed2_Var_DE2_coeffs(  6,  3) =    281.5    
           Ed2_Var_DE2_coeffs(  6,  4) =   -379.3    
           Ed2_Var_DE2_coeffs(  6,  5) =    88.98    
           Ed2_Var_DE2_coeffs(  6,  6) =   -331.0    
           Ed2_Var_DE2_coeffs(  6,  7) =   -3.493    
           Ed2_Var_DE2_coeffs(  6,  8) =   -289.8    
           Ed2_Var_DE2_coeffs(  6,  9) =    288.0    
           Ed2_Var_DE2_coeffs(  6, 10) =   -530.1    
           Ed2_Var_DE2_coeffs(  7,  1) =    735.7    
           Ed2_Var_DE2_coeffs(  7,  2) =   -832.0    
           Ed2_Var_DE2_coeffs(  7,  3) =   -444.9    
           Ed2_Var_DE2_coeffs(  7,  4) =    247.0    
           Ed2_Var_DE2_coeffs(  7,  5) =    241.8    
           Ed2_Var_DE2_coeffs(  7,  6) =   -312.3    
           Ed2_Var_DE2_coeffs(  7,  7) =   -121.3    
           Ed2_Var_DE2_coeffs(  7,  8) =    536.4    
           Ed2_Var_DE2_coeffs(  7,  9) =   -55.18    
           Ed2_Var_DE2_coeffs(  7, 10) =    8.500    
           Ed2_Var_DE2_coeffs(  8,  1) =    627.5    
           Ed2_Var_DE2_coeffs(  8,  2) =   -496.7    
           Ed2_Var_DE2_coeffs(  8,  3) =   -278.1    
           Ed2_Var_DE2_coeffs(  8,  4) =   -453.7    
           Ed2_Var_DE2_coeffs(  8,  5) =    54.90    
           Ed2_Var_DE2_coeffs(  8,  6) =    62.67    
           Ed2_Var_DE2_coeffs(  8,  7) =   -150.0    
           Ed2_Var_DE2_coeffs(  8,  8) =    579.1    
           Ed2_Var_DE2_coeffs(  8,  9) =    58.27    
           Ed2_Var_DE2_coeffs(  8, 10) =   -602.8    
           Ed2_Var_DE2_coeffs(  9,  1) =    361.7    
           Ed2_Var_DE2_coeffs(  9,  2) =   -545.6    
           Ed2_Var_DE2_coeffs(  9,  3) =    36.67    
           Ed2_Var_DE2_coeffs(  9,  4) =   -490.9    
           Ed2_Var_DE2_coeffs(  9,  5) =   -41.63    
           Ed2_Var_DE2_coeffs(  9,  6) =    66.74    
           Ed2_Var_DE2_coeffs(  9,  7) =   -93.76    
           Ed2_Var_DE2_coeffs(  9,  8) =    541.7    
           Ed2_Var_DE2_coeffs(  9,  9) =    66.86    
           Ed2_Var_DE2_coeffs(  9, 10) =   -397.6    
           ! Fill Ed1_Var_mod arrays:
           Ed1_Var_mod_coeffs(  1,  1) =   0.3689    
           Ed1_Var_mod_coeffs(  1,  2) =  -0.1464    
           Ed1_Var_mod_coeffs(  1,  3) =   0.2642    
           Ed1_Var_mod_coeffs(  1,  4) =  -0.1071    
           Ed1_Var_mod_coeffs(  1,  5) =  -0.2515    
           Ed1_Var_mod_coeffs(  1,  6) =  -0.2460E-01
           Ed1_Var_mod_coeffs(  1,  7) =   0.1534    
           Ed1_Var_mod_coeffs(  1,  8) =  -0.2606    
           Ed1_Var_mod_coeffs(  1,  9) =  -0.1948    
           Ed1_Var_mod_coeffs(  1, 10) =  -0.1371    
           Ed1_Var_mod_coeffs(  2,  1) =   0.5054    
           Ed1_Var_mod_coeffs(  2,  2) =  -0.1734    
           Ed1_Var_mod_coeffs(  2,  3) =   0.3682    
           Ed1_Var_mod_coeffs(  2,  4) =   0.2055E-01
           Ed1_Var_mod_coeffs(  2,  5) =  -0.3982    
           Ed1_Var_mod_coeffs(  2,  6) =  -0.8546E-01
           Ed1_Var_mod_coeffs(  2,  7) =   0.1320    
           Ed1_Var_mod_coeffs(  2,  8) =  -0.3481    
           Ed1_Var_mod_coeffs(  2,  9) =  -0.3800    
           Ed1_Var_mod_coeffs(  2, 10) =  -0.8100E-01
           Ed1_Var_mod_coeffs(  3,  1) =   0.7624    
           Ed1_Var_mod_coeffs(  3,  2) =  -0.6614E-01
           Ed1_Var_mod_coeffs(  3,  3) =   0.5370    
           Ed1_Var_mod_coeffs(  3,  4) =   0.1576    
           Ed1_Var_mod_coeffs(  3,  5) =  -0.5859    
           Ed1_Var_mod_coeffs(  3,  6) =  -0.1374    
           Ed1_Var_mod_coeffs(  3,  7) =   0.1824    
           Ed1_Var_mod_coeffs(  3,  8) =  -0.4086    
           Ed1_Var_mod_coeffs(  3,  9) =  -0.5924    
           Ed1_Var_mod_coeffs(  3, 10) =  -0.3176    
           Ed1_Var_mod_coeffs(  4,  1) =   0.6714    
           Ed1_Var_mod_coeffs(  4,  2) =   0.7712    
           Ed1_Var_mod_coeffs(  4,  3) =   0.4297    
           Ed1_Var_mod_coeffs(  4,  4) =   0.7992    
           Ed1_Var_mod_coeffs(  4,  5) =  -0.2933    
           Ed1_Var_mod_coeffs(  4,  6) =  -0.9359    
           Ed1_Var_mod_coeffs(  4,  7) =   0.1207E-01
           Ed1_Var_mod_coeffs(  4,  8) =  -0.2629E-01
           Ed1_Var_mod_coeffs(  4,  9) =  -0.2075    
           Ed1_Var_mod_coeffs(  4, 10) =  -0.6793    
           Ed1_Var_mod_coeffs(  5,  1) =    1.991    
           Ed1_Var_mod_coeffs(  5,  2) =   0.4556    
           Ed1_Var_mod_coeffs(  5,  3) =   0.5701    
           Ed1_Var_mod_coeffs(  5,  4) =   0.9651    
           Ed1_Var_mod_coeffs(  5,  5) =   -1.205    
           Ed1_Var_mod_coeffs(  5,  6) =  -0.2708    
           Ed1_Var_mod_coeffs(  5,  7) =  -0.1231    
           Ed1_Var_mod_coeffs(  5,  8) =    1.970    
           Ed1_Var_mod_coeffs(  5,  9) =   0.6183    
           Ed1_Var_mod_coeffs(  5, 10) =  -0.3351    
           Ed1_Var_mod_coeffs(  6,  1) =    11.37    
           Ed1_Var_mod_coeffs(  6,  2) =    3.430    
           Ed1_Var_mod_coeffs(  6,  3) =    4.591    
           Ed1_Var_mod_coeffs(  6,  4) =   0.3142    
           Ed1_Var_mod_coeffs(  6,  5) =   -2.890    
           Ed1_Var_mod_coeffs(  6,  6) =    3.408    
           Ed1_Var_mod_coeffs(  6,  7) =    6.666    
           Ed1_Var_mod_coeffs(  6,  8) =    9.446    
           Ed1_Var_mod_coeffs(  6,  9) =   -2.446    
           Ed1_Var_mod_coeffs(  6, 10) =    3.998    
           Ed1_Var_mod_coeffs(  7,  1) =    26.18    
           Ed1_Var_mod_coeffs(  7,  2) =   -3.227    
           Ed1_Var_mod_coeffs(  7,  3) =   0.3420    
           Ed1_Var_mod_coeffs(  7,  4) =    11.96    
           Ed1_Var_mod_coeffs(  7,  5) =   0.7151E-01
           Ed1_Var_mod_coeffs(  7,  6) =   -2.453    
           Ed1_Var_mod_coeffs(  7,  7) =   -1.087    
           Ed1_Var_mod_coeffs(  7,  8) =    16.00    
           Ed1_Var_mod_coeffs(  7,  9) =   0.3164    
           Ed1_Var_mod_coeffs(  7, 10) =   -16.73    
           Ed1_Var_mod_coeffs(  8,  1) =    34.90    
           Ed1_Var_mod_coeffs(  8,  2) =    2.695    
           Ed1_Var_mod_coeffs(  8,  3) =   -26.66    
           Ed1_Var_mod_coeffs(  8,  4) =    10.12    
           Ed1_Var_mod_coeffs(  8,  5) =   -4.703    
           Ed1_Var_mod_coeffs(  8,  6) =   -21.91    
           Ed1_Var_mod_coeffs(  8,  7) =   -5.001    
           Ed1_Var_mod_coeffs(  8,  8) =   -14.00    
           Ed1_Var_mod_coeffs(  8,  9) =    4.409    
           Ed1_Var_mod_coeffs(  8, 10) =    25.89    
           Ed1_Var_mod_coeffs(  9,  1) =    21.78    
           Ed1_Var_mod_coeffs(  9,  2) =    20.32    
           Ed1_Var_mod_coeffs(  9,  3) =   -12.09    
           Ed1_Var_mod_coeffs(  9,  4) =   -36.59    
           Ed1_Var_mod_coeffs(  9,  5) =   -7.579    
           Ed1_Var_mod_coeffs(  9,  6) =   -43.31    
           Ed1_Var_mod_coeffs(  9,  7) =   -5.271    
           Ed1_Var_mod_coeffs(  9,  8) =   -7.829    
           Ed1_Var_mod_coeffs(  9,  9) =    1.438    
           Ed1_Var_mod_coeffs(  9, 10) =    8.067    
           ! Fill Ed2_Var_mod arrays:
           Ed2_Var_mod_coeffs(  1,  1) =   0.1859    
           Ed2_Var_mod_coeffs(  1,  2) =   0.1158    
           Ed2_Var_mod_coeffs(  1,  3) =  -0.7879E-01
           Ed2_Var_mod_coeffs(  1,  4) =  -0.2283    
           Ed2_Var_mod_coeffs(  1,  5) =  -0.2671    
           Ed2_Var_mod_coeffs(  1,  6) =   0.2002    
           Ed2_Var_mod_coeffs(  1,  7) =  -0.1972    
           Ed2_Var_mod_coeffs(  1,  8) =   0.4213    
           Ed2_Var_mod_coeffs(  1,  9) =   0.8797E-02
           Ed2_Var_mod_coeffs(  1, 10) =   0.2399    
           Ed2_Var_mod_coeffs(  2,  1) =   0.3408    
           Ed2_Var_mod_coeffs(  2,  2) =  -0.2752E-01
           Ed2_Var_mod_coeffs(  2,  3) =  -0.5458E-01
           Ed2_Var_mod_coeffs(  2,  4) =  -0.3105    
           Ed2_Var_mod_coeffs(  2,  5) =  -0.3548    
           Ed2_Var_mod_coeffs(  2,  6) =   0.4124    
           Ed2_Var_mod_coeffs(  2,  7) =  -0.2695    
           Ed2_Var_mod_coeffs(  2,  8) =   0.6864    
           Ed2_Var_mod_coeffs(  2,  9) =   0.5311E-01
           Ed2_Var_mod_coeffs(  2, 10) =   0.1601    
           Ed2_Var_mod_coeffs(  3,  1) =    3.367    
           Ed2_Var_mod_coeffs(  3,  2) =  -0.4027    
           Ed2_Var_mod_coeffs(  3,  3) =    1.518    
           Ed2_Var_mod_coeffs(  3,  4) =    3.844    
           Ed2_Var_mod_coeffs(  3,  5) =   -3.587    
           Ed2_Var_mod_coeffs(  3,  6) =    3.384    
           Ed2_Var_mod_coeffs(  3,  7) =   -2.115    
           Ed2_Var_mod_coeffs(  3,  8) =    4.836    
           Ed2_Var_mod_coeffs(  3,  9) =  -0.9465    
           Ed2_Var_mod_coeffs(  3, 10) =   -1.985    
           Ed2_Var_mod_coeffs(  4,  1) =    14.11    
           Ed2_Var_mod_coeffs(  4,  2) =   -14.22    
           Ed2_Var_mod_coeffs(  4,  3) =    3.760    
           Ed2_Var_mod_coeffs(  4,  4) =   -5.607    
           Ed2_Var_mod_coeffs(  4,  5) =   -18.89    
           Ed2_Var_mod_coeffs(  4,  6) =    32.96    
           Ed2_Var_mod_coeffs(  4,  7) =   -13.36    
           Ed2_Var_mod_coeffs(  4,  8) =    19.98    
           Ed2_Var_mod_coeffs(  4,  9) =   -3.817    
           Ed2_Var_mod_coeffs(  4, 10) =    10.94    
           Ed2_Var_mod_coeffs(  5,  1) =    35.07    
           Ed2_Var_mod_coeffs(  5,  2) =   -6.786    
           Ed2_Var_mod_coeffs(  5,  3) =    9.257    
           Ed2_Var_mod_coeffs(  5,  4) =   -32.83    
           Ed2_Var_mod_coeffs(  5,  5) =   -35.61    
           Ed2_Var_mod_coeffs(  5,  6) =    17.11    
           Ed2_Var_mod_coeffs(  5,  7) =   -22.89    
           Ed2_Var_mod_coeffs(  5,  8) =   -12.63    
           Ed2_Var_mod_coeffs(  5,  9) =   -3.451    
           Ed2_Var_mod_coeffs(  5, 10) =    27.28    
           Ed2_Var_mod_coeffs(  6,  1) =    53.94    
           Ed2_Var_mod_coeffs(  6,  2) =    9.343    
           Ed2_Var_mod_coeffs(  6,  3) =   -25.74    
           Ed2_Var_mod_coeffs(  6,  4) =    12.88    
           Ed2_Var_mod_coeffs(  6,  5) =   -36.61    
           Ed2_Var_mod_coeffs(  6,  6) =   -1.473    
           Ed2_Var_mod_coeffs(  6,  7) =   -40.45    
           Ed2_Var_mod_coeffs(  6,  8) =    13.33    
           Ed2_Var_mod_coeffs(  6,  9) =    21.99    
           Ed2_Var_mod_coeffs(  6, 10) =   -30.78    
           Ed2_Var_mod_coeffs(  7,  1) =    70.94    
           Ed2_Var_mod_coeffs(  7,  2) =   -19.47    
           Ed2_Var_mod_coeffs(  7,  3) =   -25.10    
           Ed2_Var_mod_coeffs(  7,  4) =   -24.09    
           Ed2_Var_mod_coeffs(  7,  5) =    49.63    
           Ed2_Var_mod_coeffs(  7,  6) =   -6.381    
           Ed2_Var_mod_coeffs(  7,  7) =   -32.75    
           Ed2_Var_mod_coeffs(  7,  8) =    52.70    
           Ed2_Var_mod_coeffs(  7,  9) =   -3.273    
           Ed2_Var_mod_coeffs(  7, 10) =   -36.39    
           Ed2_Var_mod_coeffs(  8,  1) =    53.58    
           Ed2_Var_mod_coeffs(  8,  2) =    19.81    
           Ed2_Var_mod_coeffs(  8,  3) =   -7.521    
           Ed2_Var_mod_coeffs(  8,  4) =   -15.55    
           Ed2_Var_mod_coeffs(  8,  5) =    17.08    
           Ed2_Var_mod_coeffs(  8,  6) =   -22.46    
           Ed2_Var_mod_coeffs(  8,  7) =   -48.53    
           Ed2_Var_mod_coeffs(  8,  8) =   -26.41    
           Ed2_Var_mod_coeffs(  8,  9) =    6.781    
           Ed2_Var_mod_coeffs(  8, 10) =   -36.82    
           Ed2_Var_mod_coeffs(  9,  1) =    19.93    
           Ed2_Var_mod_coeffs(  9,  2) =   -6.475    
           Ed2_Var_mod_coeffs(  9,  3) =   -15.08    
           Ed2_Var_mod_coeffs(  9,  4) =    7.115    
           Ed2_Var_mod_coeffs(  9,  5) =   -5.437    
           Ed2_Var_mod_coeffs(  9,  6) =   -9.781    
           Ed2_Var_mod_coeffs(  9,  7) =    6.679    
           Ed2_Var_mod_coeffs(  9,  8) =   -12.50    
           Ed2_Var_mod_coeffs(  9,  9) =    11.44    
           Ed2_Var_mod_coeffs(  9, 10) =   -14.66    
       
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
       
         end subroutine SDEd_NBZ_fill_coeff_arrays
       
       end module SDEd_NBZ_coeff_mod
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       subroutine SDEd_NBZ_choose_coeff_array
       
         use SDEd_NBZ_coeff_mod, only: Ed1_Var_coeffs,Ed2_Var_coeffs, 
     |         SDEd_NBZ_fill_coeff_arrays
         use SDEd_NBZ_model_mod, only: first_NBz,b1,b2,b_interp1,
     |          b_interp2,x,b_90,nmlat,nfit
         integer :: binary_fit_by
       
         !****************************************************************
         ! Fill parameter arrays:
         !****************************************************************
         !  if(first_NBz) then
         call SDEd_NBZ_fill_coeff_arrays
       
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
     |          b_interp2(nfit),x(nfit),b_90(nfit))
         first_NBz = .false.
       
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
       
       end subroutine SDEd_NBZ_choose_coeff_array
       
       !---------------------------------------------------------------------------
       !---------------------------------------------------------------------------
       
       subroutine SDEd_NBZ_calc_model_value(mlat_glb,mlt,sinT, 
     |        imf_angle,Btrans,value1,value2,value_found)
       
         ! This version is for public distribution, along with the module
         ! file coeff_mod.f90.
         !
         ! Within the /home/hensel/fortran/fit_data/public subdirectory,
         ! this program is called by call_de2_model.
       
         use SDEd_NBZ_coeff_mod, only: average_mlat, 
     |          fit_by_Btrans,fit_by_imf_angle, 
     |          fit_by_sinT
         use SDEd_NBZ_model_mod, only: x,b1,b2, b_interp1,
     |          b_interp2, b_90, imlat,ifit, 
     |          nmlat, nfit, do_smooth, number_of_MLT_coeffs,
     |          number_of_Btrans_coeffs ,number_of_imf_angle_coeffs ,
     |          number_of_sinT_coeffs , fit_by_mlt,
     |          find_MLT_function,find_Btrans_function,
     |          find_imf_angle_function,
     |          find_sint_function
         use SDEd_Bmed_model_mod, only:slp1,slpn, islpsw,sigma
       
         implicit none
       
       
         real,intent(in):: mlt,mlat_glb,sinT, 
     |          imf_angle,Btrans
         real,intent(out):: value1,value2
         logical, intent(out) :: value_found
       
         integer :: fit_count,ierr
         logical :: mlat_found
         integer :: iMLT,iBtrans,iimf_angle,isinT
         real :: MLT_function,Btrans_function,
     |          imf_angle_function,sinT_function,
     |          imf_angle_in, avg_Btrans_array,mlat
       
         ! General:
         logical, parameter :: debug = .false.
         real, parameter :: pi = 3.1415926535898
         real, parameter :: beta_inf = 16.
         real, external :: curv2 
         real, dimension(nmlat) :: b_secderiv1,b_secderiv2,temp,
     |          b_coeff1,b_coeff2,lat_array
       
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
     |                    find_imf_angle_function(iimf_angle,imf_angle)
       
                     sinT_loop2:do isinT=1,1+(fit_by_sinT
     |                     *number_of_sinT_coeffs - fit_by_sinT) 
                        sinT_function = find_sinT_function(isinT,sinT)
       
                        x(fit_count) = MLT_function*Btrans_function
     |                        *imf_angle_function*sinT_function
       
                        if (debug) then
                           write(*,'("fit_count = ",I5)') fit_count
                           write(*,'("MLT_function = ",G15.6)') 
     |                            MLT_function
                           write(*,'("Btrans_function = ",G15.6)')
     |                            Btrans_function
                           write(*,'("imf_angle_function = ",G15.6)') 
     |                           imf_angle_function
                           write(*,'("sinT_function = ",G15.6)')
     |                            sinT_function
                           write(*,'("mlat = ",F9.2," MLT = ",F9.2,
     |                           " imf_angle = ",F9.2," sinT = ",F9.2,
     |                           " x(",I5,") = ",G15.6)') 
     |                           mlat,MLT,imf_angle,sinT,
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
       
       end subroutine SDEd_NBZ_calc_model_value
                                                                                                                

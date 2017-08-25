      module E_NBZ_model_mod
      
        use common_model_mod, only: MMX,MXNMX,ITHMX,NQ,pi,fit_by_MLT
      
        implicit none
      
        logical :: first_time
        logical :: value_dummy
      
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Arrays for fit coefficients:
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Number of fit coefficients for each mlat:
        integer :: E_NBZ_nfit
        ! Fit coefficients (dimension nmlat x nfit):
        real, allocatable, dimension(:,:) :: E_NBZ_b_Ed1,E_NBZ_b_Ed2 ! nmlat x nfit
        ! Interpolated fit coefficients (dimension nfit):
        real, allocatable, dimension(:) :: b_interp_Ed1,b_interp_Ed2
        ! Fitting functions (dimension nfit):
        real, allocatable, dimension(:) :: x
        ! Fit coefficients averaged over the pole (dimension nfit):
        real, allocatable, dimension(:) :: b_90_Ed1,b_90_Ed2
        ! Number of fit functions used by fit_data:
        integer, parameter :: E_NBZ_number_of_MLT_coeffs = 7,
     |       E_NBZ_number_of_Btrans_coeffs = 1,
     |        E_NBZ_number_of_imf_angle_coeffs = 1,
     |        E_NBZ_number_of_sinT_coeffs = 2
        ! Indices:
        integer :: imlat,ifit
      
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
        real ::  E_NBZ_AO(-MMX:MMX,MXNMX,E_NBZ_number_of_Btrans_coeffs,
     |    E_NBZ_number_of_imf_angle_coeffs,E_NBZ_number_of_sinT_coeffs)
        real :: E_NBZ_Q(NQ,0:ITHMX),E_NBZ_DQ(NQ,0:ITHMX)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Coeffmerge
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        integer, parameter :: 
     |      E_NBZ_fit_by_Btrans = 0,
     |      E_NBZ_fit_by_imf_angle = 0,
     |      E_NBZ_fit_by_sinT = 1
      
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
      
        function find_Btrans_function(index,Btrans_in) 
     |       result(Btrans_fxn_out)
          integer :: index
          real :: Btrans_in, Btrans_fxn_out
          real, parameter :: beta_inf = 16.
      
          select case (index)
          case (0) 
             Btrans_fxn_out = 1.
          case (1)
             Btrans_fxn_out = 1.
          case (2)
             Btrans_fxn_out = btrans_in /sqrt(1.+ 
     |      (btrans_in*btrans_in)/(beta_inf*beta_inf))
          end select
        end function find_Btrans_function
      
        function find_imf_angle_function (index,imf_angle_in) 
     |       result (imf_angle_fxn_out)
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
      
        function get_sinIm (mlat_in) result (sinIm)
          real :: mlat_in,sinIm
      
          sinIm = sin(mlat_in*pi/180.)/(1-0.75*
     |        cos(mlat_in*pi/180.))**0.5
      
        end function get_sinIm
      
      
      end module E_NBZ_model_mod
      
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      
      module E_NBZ_coeffmerge_mod
      
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
        real,  dimension( 23, 14) :: Ed1_sq2 
        real,  dimension( 23, 14) :: Ed2_sq2
        real,  dimension( 23) :: min_mlat,max_mlat!,average_mlat
      
      contains
      
        subroutine E_NBZ_fill_coeff_arrays
      
          ! Fill Ed1_sq2 arrays:
          Ed1_sq2(  1,  1) =    0.000    
          Ed1_sq2(  1,  2) =    0.000    
          Ed1_sq2(  1,  3) =   0.1000E+18
          Ed1_sq2(  1,  4) =   0.1000E+18
          Ed1_sq2(  1,  5) =    0.000    
          Ed1_sq2(  1,  6) =    0.000    
          Ed1_sq2(  1,  7) =   0.1000E+18
          Ed1_sq2(  1,  8) =   0.1000E+18
          Ed1_sq2(  1,  9) =    0.000    
          Ed1_sq2(  1, 10) =    0.000    
          Ed1_sq2(  1, 11) =   0.1000E+18
          Ed1_sq2(  1, 12) =   0.1000E+18
          Ed1_sq2(  1, 13) =    0.000    
          Ed1_sq2(  1, 14) =    0.000    
          Ed1_sq2(  2,  1) =  -0.3286    
          Ed1_sq2(  2,  2) =   0.4829    
          Ed1_sq2(  2,  3) =   0.4030    
          Ed1_sq2(  2,  4) =  -0.2898    
          Ed1_sq2(  2,  5) =   0.6460E-01
          Ed1_sq2(  2,  6) =   0.5854    
          Ed1_sq2(  2,  7) =  -0.2589    
          Ed1_sq2(  2,  8) =  -0.8461    
          Ed1_sq2(  2,  9) =   0.3520    
          Ed1_sq2(  2, 10) =  -0.4807    
          Ed1_sq2(  2, 11) =  -0.5200E-02
          Ed1_sq2(  2, 12) =   0.8149    
          Ed1_sq2(  2, 13) =   0.2806    
          Ed1_sq2(  2, 14) =  -0.5310    
          Ed1_sq2(  3,  1) =  -0.8238    
          Ed1_sq2(  3,  2) =    1.429    
          Ed1_sq2(  3,  3) =   0.3995    
          Ed1_sq2(  3,  4) =   -1.435    
          Ed1_sq2(  3,  5) =   0.1156    
          Ed1_sq2(  3,  6) =  -0.7950    
          Ed1_sq2(  3,  7) =   0.1127    
          Ed1_sq2(  3,  8) =   -1.668    
          Ed1_sq2(  3,  9) =  -0.3055    
          Ed1_sq2(  3, 10) =  -0.5728    
          Ed1_sq2(  3, 11) =   0.1825    
          Ed1_sq2(  3, 12) =   0.7092    
          Ed1_sq2(  3, 13) =  -0.7650E-02
          Ed1_sq2(  3, 14) =   0.4509    
          Ed1_sq2(  4,  1) =   -1.271    
          Ed1_sq2(  4,  2) =    2.470    
          Ed1_sq2(  4,  3) =  -0.2926    
          Ed1_sq2(  4,  4) =    1.290    
          Ed1_sq2(  4,  5) =   0.2440    
          Ed1_sq2(  4,  6) =   -1.230    
          Ed1_sq2(  4,  7) =  -0.1585    
          Ed1_sq2(  4,  8) =  -0.3055    
          Ed1_sq2(  4,  9) =  -0.6135    
          Ed1_sq2(  4, 10) =  -0.4000E-01
          Ed1_sq2(  4, 11) =  -0.5600E-02
          Ed1_sq2(  4, 12) =    1.339    
          Ed1_sq2(  4, 13) =  -0.2948    
          Ed1_sq2(  4, 14) =    1.508    
          Ed1_sq2(  5,  1) =   -1.343    
          Ed1_sq2(  5,  2) =    2.896    
          Ed1_sq2(  5,  3) =   0.7075E-01
          Ed1_sq2(  5,  4) =  -0.1773    
          Ed1_sq2(  5,  5) =   0.2904    
          Ed1_sq2(  5,  6) =   -2.385    
          Ed1_sq2(  5,  7) =  -0.2340    
          Ed1_sq2(  5,  8) =  -0.4050E-01
          Ed1_sq2(  5,  9) =  -0.6550    
          Ed1_sq2(  5, 10) =  -0.1825    
          Ed1_sq2(  5, 11) =  -0.1451    
          Ed1_sq2(  5, 12) =   0.8284    
          Ed1_sq2(  5, 13) =  -0.1831    
          Ed1_sq2(  5, 14) =    1.060    
          Ed1_sq2(  6,  1) =   -1.263    
          Ed1_sq2(  6,  2) =    2.615    
          Ed1_sq2(  6,  3) =   0.4037    
          Ed1_sq2(  6,  4) =  -0.3821    
          Ed1_sq2(  6,  5) =   0.2103    
          Ed1_sq2(  6,  6) =   -3.027    
          Ed1_sq2(  6,  7) =  -0.1470    
          Ed1_sq2(  6,  8) =  -0.5215    
          Ed1_sq2(  6,  9) =  -0.6810    
          Ed1_sq2(  6, 10) =  -0.3611    
          Ed1_sq2(  6, 11) =  -0.4093    
          Ed1_sq2(  6, 12) =   0.9434    
          Ed1_sq2(  6, 13) =  -0.3034    
          Ed1_sq2(  6, 14) =   0.5560    
          Ed1_sq2(  7,  1) =   -1.166    
          Ed1_sq2(  7,  2) =    2.268    
          Ed1_sq2(  7,  3) =   0.5411    
          Ed1_sq2(  7,  4) =  -0.4909    
          Ed1_sq2(  7,  5) =   0.3160E-01
          Ed1_sq2(  7,  6) =   -3.079    
          Ed1_sq2(  7,  7) =  -0.1888    
          Ed1_sq2(  7,  8) =   -1.173    
          Ed1_sq2(  7,  9) =  -0.6882    
          Ed1_sq2(  7, 10) =  -0.2103    
          Ed1_sq2(  7, 11) =  -0.5001    
          Ed1_sq2(  7, 12) =  -0.1130    
          Ed1_sq2(  7, 13) =  -0.1821    
          Ed1_sq2(  7, 14) =   0.5740    
          Ed1_sq2(  8,  1) =  -0.8435    
          Ed1_sq2(  8,  2) =    2.196    
          Ed1_sq2(  8,  3) =   0.7697    
          Ed1_sq2(  8,  4) =   -1.033    
          Ed1_sq2(  8,  5) =  -0.5880    
          Ed1_sq2(  8,  6) =   -1.545    
          Ed1_sq2(  8,  7) =   0.4065E-01
          Ed1_sq2(  8,  8) =   -1.425    
          Ed1_sq2(  8,  9) =   -1.210    
          Ed1_sq2(  8, 10) =    1.500    
          Ed1_sq2(  8, 11) =  -0.8537    
          Ed1_sq2(  8, 12) =    1.032    
          Ed1_sq2(  8, 13) =  -0.4651    
          Ed1_sq2(  8, 14) =    1.162    
          Ed1_sq2(  9,  1) =  -0.5550    
          Ed1_sq2(  9,  2) =    2.530    
          Ed1_sq2(  9,  3) =   0.8480    
          Ed1_sq2(  9,  4) =  -0.9928    
          Ed1_sq2(  9,  5) =   -1.008    
          Ed1_sq2(  9,  6) =   -1.953    
          Ed1_sq2(  9,  7) =   0.4760    
          Ed1_sq2(  9,  8) =   -1.295    
          Ed1_sq2(  9,  9) =   -1.578    
          Ed1_sq2(  9, 10) =    2.186    
          Ed1_sq2(  9, 11) =  -0.5330    
          Ed1_sq2(  9, 12) =    1.324    
          Ed1_sq2(  9, 13) =  -0.7912    
          Ed1_sq2(  9, 14) =   0.7451    
          Ed1_sq2( 10,  1) =  -0.2857    
          Ed1_sq2( 10,  2) =    3.329    
          Ed1_sq2( 10,  3) =    1.087    
          Ed1_sq2( 10,  4) =  -0.4340    
          Ed1_sq2( 10,  5) =  -0.9108    
          Ed1_sq2( 10,  6) =   -2.910    
          Ed1_sq2( 10,  7) =    1.105    
          Ed1_sq2( 10,  8) =   -1.309    
          Ed1_sq2( 10,  9) =   -1.106    
          Ed1_sq2( 10, 10) =    2.069    
          Ed1_sq2( 10, 11) =   0.5200E-01
          Ed1_sq2( 10, 12) =    1.117    
          Ed1_sq2( 10, 13) =  -0.5524    
          Ed1_sq2( 10, 14) =   0.2066    
          Ed1_sq2( 11,  1) =  -0.2420    
          Ed1_sq2( 11,  2) =    3.230    
          Ed1_sq2( 11,  3) =    1.194    
          Ed1_sq2( 11,  4) =  -0.8658    
          Ed1_sq2( 11,  5) =  -0.9973    
          Ed1_sq2( 11,  6) =   -2.431    
          Ed1_sq2( 11,  7) =    1.864    
          Ed1_sq2( 11,  8) =   -4.136    
          Ed1_sq2( 11,  9) =  -0.8740    
          Ed1_sq2( 11, 10) =   0.1300    
          Ed1_sq2( 11, 11) =    1.010    
          Ed1_sq2( 11, 12) =   -1.077    
          Ed1_sq2( 11, 13) =   -1.279    
          Ed1_sq2( 11, 14) =   -1.082    
          Ed1_sq2( 12,  1) =  -0.3500    
          Ed1_sq2( 12,  2) =    3.775    
          Ed1_sq2( 12,  3) =   0.1705    
          Ed1_sq2( 12,  4) =   -1.453    
          Ed1_sq2( 12,  5) =   -1.472    
          Ed1_sq2( 12,  6) =  -0.4521    
          Ed1_sq2( 12,  7) =    1.437    
          Ed1_sq2( 12,  8) =   -3.680    
          Ed1_sq2( 12,  9) =   -1.507    
          Ed1_sq2( 12, 10) =    1.533    
          Ed1_sq2( 12, 11) =   0.5795    
          Ed1_sq2( 12, 12) =   -3.313    
          Ed1_sq2( 12, 13) =   -1.732    
          Ed1_sq2( 12, 14) =   0.3390E-01
          Ed1_sq2( 13,  1) =   -1.130    
          Ed1_sq2( 13,  2) =    5.840    
          Ed1_sq2( 13,  3) =   -1.680    
          Ed1_sq2( 13,  4) =   0.9511    
          Ed1_sq2( 13,  5) =   -1.406    
          Ed1_sq2( 13,  6) =   0.3850    
          Ed1_sq2( 13,  7) =  -0.7046    
          Ed1_sq2( 13,  8) =   -1.190    
          Ed1_sq2( 13,  9) =   -1.403    
          Ed1_sq2( 13, 10) =    4.472    
          Ed1_sq2( 13, 11) =   0.6750E-01
          Ed1_sq2( 13, 12) =   -1.747    
          Ed1_sq2( 13, 13) =   -1.099    
          Ed1_sq2( 13, 14) =   0.3230    
          Ed1_sq2( 14,  1) =  -0.1160    
          Ed1_sq2( 14,  2) =    4.394    
          Ed1_sq2( 14,  3) =   -1.516    
          Ed1_sq2( 14,  4) =  -0.1285    
          Ed1_sq2( 14,  5) =   0.2321    
          Ed1_sq2( 14,  6) =   -3.149    
          Ed1_sq2( 14,  7) =   -2.918    
          Ed1_sq2( 14,  8) =    2.953    
          Ed1_sq2( 14,  9) =    1.655    
          Ed1_sq2( 14, 10) =   -2.201    
          Ed1_sq2( 14, 11) =  -0.9985    
          Ed1_sq2( 14, 12) =  -0.2875E-01
          Ed1_sq2( 14, 13) =    1.300    
          Ed1_sq2( 14, 14) =   -1.492    
          Ed1_sq2( 15,  1) =   0.5390    
          Ed1_sq2( 15,  2) =   0.2550    
          Ed1_sq2( 15,  3) =   -1.671    
          Ed1_sq2( 15,  4) =   0.6950E-01
          Ed1_sq2( 15,  5) =   0.1725    
          Ed1_sq2( 15,  6) =   -5.611    
          Ed1_sq2( 15,  7) =   -4.070    
          Ed1_sq2( 15,  8) =   -2.603    
          Ed1_sq2( 15,  9) =    1.574    
          Ed1_sq2( 15, 10) =   -2.913    
          Ed1_sq2( 15, 11) =  -0.2685    
          Ed1_sq2( 15, 12) =  -0.1765    
          Ed1_sq2( 15, 13) =    2.362    
          Ed1_sq2( 15, 14) =  -0.8250E-01
          Ed1_sq2( 16,  1) =   0.4599    
          Ed1_sq2( 16,  2) =    3.821    
          Ed1_sq2( 16,  3) =   -4.003    
          Ed1_sq2( 16,  4) =   -1.213    
          Ed1_sq2( 16,  5) =   0.6904    
          Ed1_sq2( 16,  6) =   -4.826    
          Ed1_sq2( 16,  7) =  -0.5175    
          Ed1_sq2( 16,  8) =   -3.527    
          Ed1_sq2( 16,  9) =   0.4255    
          Ed1_sq2( 16, 10) =   -2.434    
          Ed1_sq2( 16, 11) =   -1.202    
          Ed1_sq2( 16, 12) =  -0.4080    
          Ed1_sq2( 16, 13) =    2.569    
          Ed1_sq2( 16, 14) =   -3.134    
          Ed1_sq2( 17,  1) =    2.788    
          Ed1_sq2( 17,  2) =   0.5149    
          Ed1_sq2( 17,  3) =   -3.283    
          Ed1_sq2( 17,  4) =   -8.041    
          Ed1_sq2( 17,  5) =   0.6272    
          Ed1_sq2( 17,  6) =    4.516    
          Ed1_sq2( 17,  7) =    3.415    
          Ed1_sq2( 17,  8) =   0.2550E-01
          Ed1_sq2( 17,  9) =   -2.028    
          Ed1_sq2( 17, 10) =   -3.796    
          Ed1_sq2( 17, 11) =   0.8160    
          Ed1_sq2( 17, 12) =   -7.574    
          Ed1_sq2( 17, 13) =    2.692    
          Ed1_sq2( 17, 14) =    6.291    
          Ed1_sq2( 18,  1) =    1.710    
          Ed1_sq2( 18,  2) =   0.5135    
          Ed1_sq2( 18,  3) =   -1.873    
          Ed1_sq2( 18,  4) =   -7.379    
          Ed1_sq2( 18,  5) =   -1.712    
          Ed1_sq2( 18,  6) =    8.665    
          Ed1_sq2( 18,  7) =   -1.523    
          Ed1_sq2( 18,  8) =    3.629    
          Ed1_sq2( 18,  9) =   -1.830    
          Ed1_sq2( 18, 10) =   -4.333    
          Ed1_sq2( 18, 11) =    5.565    
          Ed1_sq2( 18, 12) =   -8.590    
          Ed1_sq2( 18, 13) =  -0.5000E-01
          Ed1_sq2( 18, 14) =    7.626    
          Ed1_sq2( 19,  1) =   -1.188    
          Ed1_sq2( 19,  2) =  -0.3345    
          Ed1_sq2( 19,  3) =    1.083    
          Ed1_sq2( 19,  4) =   -1.347    
          Ed1_sq2( 19,  5) =   0.6400E-02
          Ed1_sq2( 19,  6) =    1.586    
          Ed1_sq2( 19,  7) =   -9.822    
          Ed1_sq2( 19,  8) =    2.959    
          Ed1_sq2( 19,  9) =    1.032    
          Ed1_sq2( 19, 10) =   -3.166    
          Ed1_sq2( 19, 11) =    11.23    
          Ed1_sq2( 19, 12) =    7.580    
          Ed1_sq2( 19, 13) =   -1.587    
          Ed1_sq2( 19, 14) =    5.287    
          Ed1_sq2( 20,  1) =   -1.394    
          Ed1_sq2( 20,  2) =   -6.108    
          Ed1_sq2( 20,  3) =    2.369    
          Ed1_sq2( 20,  4) =   -7.814    
          Ed1_sq2( 20,  5) =   -1.763    
          Ed1_sq2( 20,  6) =   -15.41    
          Ed1_sq2( 20,  7) =   -9.845    
          Ed1_sq2( 20,  8) =   -1.788    
          Ed1_sq2( 20,  9) =    1.875    
          Ed1_sq2( 20, 10) =   -3.257    
          Ed1_sq2( 20, 11) =    7.332    
          Ed1_sq2( 20, 12) =   -2.372    
          Ed1_sq2( 20, 13) =   -1.451    
          Ed1_sq2( 20, 14) =    3.374    
          Ed1_sq2( 21,  1) =  -0.1220    
          Ed1_sq2( 21,  2) =   -10.43    
          Ed1_sq2( 21,  3) =    2.894    
          Ed1_sq2( 21,  4) =    5.050    
          Ed1_sq2( 21,  5) =    1.722    
          Ed1_sq2( 21,  6) =   -3.618    
          Ed1_sq2( 21,  7) =   -7.803    
          Ed1_sq2( 21,  8) =   -7.703    
          Ed1_sq2( 21,  9) =    2.718    
          Ed1_sq2( 21, 10) =   -8.948    
          Ed1_sq2( 21, 11) =    3.447    
          Ed1_sq2( 21, 12) =    8.184    
          Ed1_sq2( 21, 13) =    1.457    
          Ed1_sq2( 21, 14) =  -0.8225    
          Ed1_sq2( 22,  1) =   -1.507    
          Ed1_sq2( 22,  2) =    10.46    
          Ed1_sq2( 22,  3) =    4.136    
          Ed1_sq2( 22,  4) =   -7.674    
          Ed1_sq2( 22,  5) =   0.7519    
          Ed1_sq2( 22,  6) =    3.399    
          Ed1_sq2( 22,  7) =   -4.582    
          Ed1_sq2( 22,  8) =  -0.4000E-01
          Ed1_sq2( 22,  9) =    1.456    
          Ed1_sq2( 22, 10) =   -2.853    
          Ed1_sq2( 22, 11) =   0.7680    
          Ed1_sq2( 22, 12) =   -1.036    
          Ed1_sq2( 22, 13) =   -2.550    
          Ed1_sq2( 22, 14) =    8.985    
          Ed1_sq2( 23,  1) =    4.859    
          Ed1_sq2( 23,  2) =   -11.08    
          Ed1_sq2( 23,  3) =    6.498    
          Ed1_sq2( 23,  4) =   0.9169    
          Ed1_sq2( 23,  5) =   -8.542    
          Ed1_sq2( 23,  6) =    19.95    
          Ed1_sq2( 23,  7) =   -1.427    
          Ed1_sq2( 23,  8) =   -7.339    
          Ed1_sq2( 23,  9) =  -0.9055    
          Ed1_sq2( 23, 10) =    6.705    
          Ed1_sq2( 23, 11) =   0.2870    
          Ed1_sq2( 23, 12) =    7.553    
          Ed1_sq2( 23, 13) =  -0.1640    
          Ed1_sq2( 23, 14) =    6.526    
      
          ! Fill Ed2_sq2 arrays:
          Ed2_sq2(  1,  1) =   0.1000E+18
          Ed2_sq2(  1,  2) =   0.1000E+18
          Ed2_sq2(  1,  3) =   0.2000E+18
          Ed2_sq2(  1,  4) =   0.2000E+18
          Ed2_sq2(  1,  5) =   0.1000E+18
          Ed2_sq2(  1,  6) =   0.1000E+18
          Ed2_sq2(  1,  7) =    0.000    
          Ed2_sq2(  1,  8) =    0.000    
          Ed2_sq2(  1,  9) =   0.1000E+18
          Ed2_sq2(  1, 10) =   0.1000E+18
          Ed2_sq2(  1, 11) =    0.000    
          Ed2_sq2(  1, 12) =    0.000    
          Ed2_sq2(  1, 13) =   0.1000E+18
          Ed2_sq2(  1, 14) =   0.1000E+18
          Ed2_sq2(  2,  1) =   -1.911    
          Ed2_sq2(  2,  2) =   0.9523    
          Ed2_sq2(  2,  3) =   -1.483    
          Ed2_sq2(  2,  4) =   -1.096    
          Ed2_sq2(  2,  5) =   -1.119    
          Ed2_sq2(  2,  6) =   0.5203    
          Ed2_sq2(  2,  7) =  -0.3520    
          Ed2_sq2(  2,  8) =   0.4807    
          Ed2_sq2(  2,  9) =  -0.2589    
          Ed2_sq2(  2, 10) =  -0.8461    
          Ed2_sq2(  2, 11) =  -0.2806    
          Ed2_sq2(  2, 12) =   0.5310    
          Ed2_sq2(  2, 13) =  -0.5200E-02
          Ed2_sq2(  2, 14) =   0.8149    
          Ed2_sq2(  3,  1) =   -1.758    
          Ed2_sq2(  3,  2) =    1.321    
          Ed2_sq2(  3,  3) =  -0.9631    
          Ed2_sq2(  3,  4) =   -1.547    
          Ed2_sq2(  3,  5) =  -0.8374E-01
          Ed2_sq2(  3,  6) =   0.1711    
          Ed2_sq2(  3,  7) =   0.3055    
          Ed2_sq2(  3,  8) =   0.5728    
          Ed2_sq2(  3,  9) =   0.1127    
          Ed2_sq2(  3, 10) =   -1.668    
          Ed2_sq2(  3, 11) =   0.7650E-02
          Ed2_sq2(  3, 12) =  -0.4509    
          Ed2_sq2(  3, 13) =   0.1825    
          Ed2_sq2(  3, 14) =   0.7092    
          Ed2_sq2(  4,  1) =   -1.762    
          Ed2_sq2(  4,  2) =    1.175    
          Ed2_sq2(  4,  3) =  -0.9038    
          Ed2_sq2(  4,  4) =   -2.253    
          Ed2_sq2(  4,  5) =   0.1336    
          Ed2_sq2(  4,  6) =   0.8620    
          Ed2_sq2(  4,  7) =   0.6135    
          Ed2_sq2(  4,  8) =   0.4000E-01
          Ed2_sq2(  4,  9) =  -0.1585    
          Ed2_sq2(  4, 10) =  -0.3055    
          Ed2_sq2(  4, 11) =   0.2948    
          Ed2_sq2(  4, 12) =   -1.508    
          Ed2_sq2(  4, 13) =  -0.5600E-02
          Ed2_sq2(  4, 14) =    1.339    
          Ed2_sq2(  5,  1) =   -1.882    
          Ed2_sq2(  5,  2) =    1.434    
          Ed2_sq2(  5,  3) =  -0.9953    
          Ed2_sq2(  5,  4) =   -2.317    
          Ed2_sq2(  5,  5) =   0.6148    
          Ed2_sq2(  5,  6) =   0.4279    
          Ed2_sq2(  5,  7) =   0.6550    
          Ed2_sq2(  5,  8) =   0.1825    
          Ed2_sq2(  5,  9) =  -0.2340    
          Ed2_sq2(  5, 10) =  -0.4050E-01
          Ed2_sq2(  5, 11) =   0.1831    
          Ed2_sq2(  5, 12) =   -1.060    
          Ed2_sq2(  5, 13) =  -0.1451    
          Ed2_sq2(  5, 14) =   0.8284    
          Ed2_sq2(  6,  1) =   -1.939    
          Ed2_sq2(  6,  2) =    1.322    
          Ed2_sq2(  6,  3) =   -1.018    
          Ed2_sq2(  6,  4) =   -1.605    
          Ed2_sq2(  6,  5) =    1.269    
          Ed2_sq2(  6,  6) =   0.5990    
          Ed2_sq2(  6,  7) =   0.6810    
          Ed2_sq2(  6,  8) =   0.3611    
          Ed2_sq2(  6,  9) =  -0.1470    
          Ed2_sq2(  6, 10) =  -0.5215    
          Ed2_sq2(  6, 11) =   0.3034    
          Ed2_sq2(  6, 12) =  -0.5560    
          Ed2_sq2(  6, 13) =  -0.4093    
          Ed2_sq2(  6, 14) =   0.9434    
          Ed2_sq2(  7,  1) =   -1.657    
          Ed2_sq2(  7,  2) =   0.2632    
          Ed2_sq2(  7,  3) =  -0.6312    
          Ed2_sq2(  7,  4) =   -2.145    
          Ed2_sq2(  7,  5) =    1.801    
          Ed2_sq2(  7,  6) =   0.3920    
          Ed2_sq2(  7,  7) =   0.6882    
          Ed2_sq2(  7,  8) =   0.2103    
          Ed2_sq2(  7,  9) =  -0.1888    
          Ed2_sq2(  7, 10) =   -1.173    
          Ed2_sq2(  7, 11) =   0.1821    
          Ed2_sq2(  7, 12) =  -0.5740    
          Ed2_sq2(  7, 13) =  -0.5001    
          Ed2_sq2(  7, 14) =  -0.1130    
          Ed2_sq2(  8,  1) =   -1.913    
          Ed2_sq2(  8,  2) =    1.763    
          Ed2_sq2(  8,  3) =  -0.7598    
          Ed2_sq2(  8,  4) =  -0.6953    
          Ed2_sq2(  8,  5) =    3.350    
          Ed2_sq2(  8,  6) =   -2.795    
          Ed2_sq2(  8,  7) =    1.210    
          Ed2_sq2(  8,  8) =   -1.500    
          Ed2_sq2(  8,  9) =   0.4065E-01
          Ed2_sq2(  8, 10) =   -1.425    
          Ed2_sq2(  8, 11) =   0.4651    
          Ed2_sq2(  8, 12) =   -1.162    
          Ed2_sq2(  8, 13) =  -0.8537    
          Ed2_sq2(  8, 14) =    1.032    
          Ed2_sq2(  9,  1) =   -2.161    
          Ed2_sq2(  9,  2) =    2.150    
          Ed2_sq2(  9,  3) =  -0.4980E-01
          Ed2_sq2(  9,  4) =   0.8691    
          Ed2_sq2(  9,  5) =    4.380    
          Ed2_sq2(  9,  6) =   -2.575    
          Ed2_sq2(  9,  7) =    1.578    
          Ed2_sq2(  9,  8) =   -2.186    
          Ed2_sq2(  9,  9) =   0.4760    
          Ed2_sq2(  9, 10) =   -1.295    
          Ed2_sq2(  9, 11) =   0.7912    
          Ed2_sq2(  9, 12) =  -0.7451    
          Ed2_sq2(  9, 13) =  -0.5330    
          Ed2_sq2(  9, 14) =    1.324    
          Ed2_sq2( 10,  1) =   -2.421    
          Ed2_sq2( 10,  2) =    1.112    
          Ed2_sq2( 10,  3) =   0.3859    
          Ed2_sq2( 10,  4) =   0.7716    
          Ed2_sq2( 10,  5) =    5.009    
          Ed2_sq2( 10,  6) =  -0.2110    
          Ed2_sq2( 10,  7) =    1.106    
          Ed2_sq2( 10,  8) =   -2.069    
          Ed2_sq2( 10,  9) =    1.105    
          Ed2_sq2( 10, 10) =   -1.309    
          Ed2_sq2( 10, 11) =   0.5524    
          Ed2_sq2( 10, 12) =  -0.2066    
          Ed2_sq2( 10, 13) =   0.5200E-01
          Ed2_sq2( 10, 14) =    1.117    
          Ed2_sq2( 11,  1) =   -2.989    
          Ed2_sq2( 11,  2) =    1.298    
          Ed2_sq2( 11,  3) =   0.3663    
          Ed2_sq2( 11,  4) =    5.036    
          Ed2_sq2( 11,  5) =    5.825    
          Ed2_sq2( 11,  6) =   0.5822    
          Ed2_sq2( 11,  7) =   0.8740    
          Ed2_sq2( 11,  8) =  -0.1300    
          Ed2_sq2( 11,  9) =    1.864    
          Ed2_sq2( 11, 10) =   -4.136    
          Ed2_sq2( 11, 11) =    1.279    
          Ed2_sq2( 11, 12) =    1.082    
          Ed2_sq2( 11, 13) =    1.010    
          Ed2_sq2( 11, 14) =   -1.077    
          Ed2_sq2( 12,  1) =   -3.357    
          Ed2_sq2( 12,  2) =   0.5270    
          Ed2_sq2( 12,  3) =    1.541    
          Ed2_sq2( 12,  4) =    2.002    
          Ed2_sq2( 12,  5) =    7.355    
          Ed2_sq2( 12,  6) =    2.176    
          Ed2_sq2( 12,  7) =    1.507    
          Ed2_sq2( 12,  8) =   -1.533    
          Ed2_sq2( 12,  9) =    1.437    
          Ed2_sq2( 12, 10) =   -3.680    
          Ed2_sq2( 12, 11) =    1.732    
          Ed2_sq2( 12, 12) =  -0.3390E-01
          Ed2_sq2( 12, 13) =   0.5795    
          Ed2_sq2( 12, 14) =   -3.313    
          Ed2_sq2( 13,  1) =   -2.055    
          Ed2_sq2( 13,  2) =   -2.668    
          Ed2_sq2( 13,  3) =    4.390    
          Ed2_sq2( 13,  4) =   -6.151    
          Ed2_sq2( 13,  5) =    8.965    
          Ed2_sq2( 13,  6) =    2.975    
          Ed2_sq2( 13,  7) =    1.403    
          Ed2_sq2( 13,  8) =   -4.472    
          Ed2_sq2( 13,  9) =  -0.7046    
          Ed2_sq2( 13, 10) =   -1.190    
          Ed2_sq2( 13, 11) =    1.099    
          Ed2_sq2( 13, 12) =  -0.3230    
          Ed2_sq2( 13, 13) =   0.6750E-01
          Ed2_sq2( 13, 14) =   -1.747    
          Ed2_sq2( 14,  1) =   -5.245    
          Ed2_sq2( 14,  2) =    1.836    
          Ed2_sq2( 14,  3) =   0.6143    
          Ed2_sq2( 14,  4) =  -0.6670    
          Ed2_sq2( 14,  5) =    9.800    
          Ed2_sq2( 14,  6) =   -1.013    
          Ed2_sq2( 14,  7) =   -1.655    
          Ed2_sq2( 14,  8) =    2.201    
          Ed2_sq2( 14,  9) =   -2.918    
          Ed2_sq2( 14, 10) =    2.953    
          Ed2_sq2( 14, 11) =   -1.300    
          Ed2_sq2( 14, 12) =    1.492    
          Ed2_sq2( 14, 13) =  -0.9985    
          Ed2_sq2( 14, 14) =  -0.2875E-01
          Ed2_sq2( 15,  1) =   -4.920    
          Ed2_sq2( 15,  2) =    1.495    
          Ed2_sq2( 15,  3) =    2.355    
          Ed2_sq2( 15,  4) =   -6.836    
          Ed2_sq2( 15,  5) =    8.462    
          Ed2_sq2( 15,  6) =   -2.872    
          Ed2_sq2( 15,  7) =   -1.574    
          Ed2_sq2( 15,  8) =    2.913    
          Ed2_sq2( 15,  9) =   -4.070    
          Ed2_sq2( 15, 10) =   -2.603    
          Ed2_sq2( 15, 11) =   -2.362    
          Ed2_sq2( 15, 12) =   0.8250E-01
          Ed2_sq2( 15, 13) =  -0.2685    
          Ed2_sq2( 15, 14) =  -0.1765    
          Ed2_sq2( 16,  1) =   -5.587    
          Ed2_sq2( 16,  2) =  -0.3385    
          Ed2_sq2( 16,  3) =   0.1806    
          Ed2_sq2( 16,  4) =  -0.8163    
          Ed2_sq2( 16,  5) =   -3.221    
          Ed2_sq2( 16,  6) =    13.73    
          Ed2_sq2( 16,  7) =  -0.4255    
          Ed2_sq2( 16,  8) =    2.434    
          Ed2_sq2( 16,  9) =  -0.5175    
          Ed2_sq2( 16, 10) =   -3.527    
          Ed2_sq2( 16, 11) =   -2.569    
          Ed2_sq2( 16, 12) =    3.134    
          Ed2_sq2( 16, 13) =   -1.202    
          Ed2_sq2( 16, 14) =  -0.4080    
          Ed2_sq2( 17,  1) =   -5.378    
          Ed2_sq2( 17,  2) =   -2.421    
          Ed2_sq2( 17,  3) =   -1.328    
          Ed2_sq2( 17,  4) =    9.962    
          Ed2_sq2( 17,  5) =   -16.70    
          Ed2_sq2( 17,  6) =    7.791    
          Ed2_sq2( 17,  7) =    2.028    
          Ed2_sq2( 17,  8) =    3.796    
          Ed2_sq2( 17,  9) =    3.415    
          Ed2_sq2( 17, 10) =   0.2550E-01
          Ed2_sq2( 17, 11) =   -2.692    
          Ed2_sq2( 17, 12) =   -6.291    
          Ed2_sq2( 17, 13) =   0.8160    
          Ed2_sq2( 17, 14) =   -7.574    
          Ed2_sq2( 18,  1) =  -0.2955    
          Ed2_sq2( 18,  2) =  -0.1415    
          Ed2_sq2( 18,  3) =   -2.576    
          Ed2_sq2( 18,  4) =   -3.637    
          Ed2_sq2( 18,  5) =   -17.68    
          Ed2_sq2( 18,  6) =   -11.14    
          Ed2_sq2( 18,  7) =    1.830    
          Ed2_sq2( 18,  8) =    4.333    
          Ed2_sq2( 18,  9) =   -1.523    
          Ed2_sq2( 18, 10) =    3.629    
          Ed2_sq2( 18, 11) =   0.5000E-01
          Ed2_sq2( 18, 12) =   -7.626    
          Ed2_sq2( 18, 13) =    5.565    
          Ed2_sq2( 18, 14) =   -8.590    
          Ed2_sq2( 19,  1) =  -0.7023    
          Ed2_sq2( 19,  2) =    6.200    
          Ed2_sq2( 19,  3) =   -1.563    
          Ed2_sq2( 19,  4) =   -11.79    
          Ed2_sq2( 19,  5) =   -6.885    
          Ed2_sq2( 19,  6) =   -9.368    
          Ed2_sq2( 19,  7) =   -1.032    
          Ed2_sq2( 19,  8) =    3.166    
          Ed2_sq2( 19,  9) =   -9.822    
          Ed2_sq2( 19, 10) =    2.959    
          Ed2_sq2( 19, 11) =    1.587    
          Ed2_sq2( 19, 12) =   -5.287    
          Ed2_sq2( 19, 13) =    11.23    
          Ed2_sq2( 19, 14) =    7.580    
          Ed2_sq2( 20,  1) =   -1.252    
          Ed2_sq2( 20,  2) =   0.1360    
          Ed2_sq2( 20,  3) =    1.789    
          Ed2_sq2( 20,  4) =   -4.063    
          Ed2_sq2( 20,  5) =  -0.1870    
          Ed2_sq2( 20,  6) =    7.701    
          Ed2_sq2( 20,  7) =   -1.875    
          Ed2_sq2( 20,  8) =    3.257    
          Ed2_sq2( 20,  9) =   -9.845    
          Ed2_sq2( 20, 10) =   -1.788    
          Ed2_sq2( 20, 11) =    1.451    
          Ed2_sq2( 20, 12) =   -3.374    
          Ed2_sq2( 20, 13) =    7.332    
          Ed2_sq2( 20, 14) =   -2.372    
          Ed2_sq2( 21,  1) =   -1.479    
          Ed2_sq2( 21,  2) =    3.141    
          Ed2_sq2( 21,  3) =    1.031    
          Ed2_sq2( 21,  4) =   -6.394    
          Ed2_sq2( 21,  5) =    1.042    
          Ed2_sq2( 21,  6) =    8.623    
          Ed2_sq2( 21,  7) =   -2.718    
          Ed2_sq2( 21,  8) =    8.948    
          Ed2_sq2( 21,  9) =   -7.803    
          Ed2_sq2( 21, 10) =   -7.703    
          Ed2_sq2( 21, 11) =   -1.457    
          Ed2_sq2( 21, 12) =   0.8225    
          Ed2_sq2( 21, 13) =    3.447    
          Ed2_sq2( 21, 14) =    8.184    
          Ed2_sq2( 22,  1) =   -3.546    
          Ed2_sq2( 22,  2) =    2.609    
          Ed2_sq2( 22,  3) =   -1.040    
          Ed2_sq2( 22,  4) =    14.43    
          Ed2_sq2( 22,  5) =  -0.9500    
          Ed2_sq2( 22,  6) =    12.73    
          Ed2_sq2( 22,  7) =   -1.456    
          Ed2_sq2( 22,  8) =    2.853    
          Ed2_sq2( 22,  9) =   -4.582    
          Ed2_sq2( 22, 10) =  -0.4000E-01
          Ed2_sq2( 22, 11) =    2.550    
          Ed2_sq2( 22, 12) =   -8.985    
          Ed2_sq2( 22, 13) =   0.7680    
          Ed2_sq2( 22, 14) =   -1.036    
          Ed2_sq2( 23,  1) =   -2.773    
          Ed2_sq2( 23,  2) =    18.02    
          Ed2_sq2( 23,  3) =   0.5180    
          Ed2_sq2( 23,  4) =   -2.186    
          Ed2_sq2( 23,  5) =    11.95    
          Ed2_sq2( 23,  6) =   -13.22    
          Ed2_sq2( 23,  7) =   0.9055    
          Ed2_sq2( 23,  8) =   -6.705    
          Ed2_sq2( 23,  9) =   -1.427    
          Ed2_sq2( 23, 10) =   -7.339    
          Ed2_sq2( 23, 11) =   0.1640    
          Ed2_sq2( 23, 12) =   -6.526    
          Ed2_sq2( 23, 13) =   0.2870    
          Ed2_sq2( 23, 14) =    7.553    
      
      
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
      
        end subroutine E_NBZ_fill_coeff_arrays
      
      end module E_NBZ_coeffmerge_mod
      
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      
      subroutine E_NBZ_calc_model_value(mlat_glob,mlt,sinT,imf_angle,
     |      value_Ed1,value_Ed2,value,value_found)
      
        ! This version is for public distribution, along with the module
        ! file coeff_mod.f90.
        !
        ! Within the /home/hensel/fortran/fit_data/public subdirectory,
        ! this program is called by call_de2_model.
      
        use common_model_mod, only: MMX,MXNMX,ITHMX,mlat_resolution
        use E_NBZ_model_mod, only: 
     |      pi,find_Btrans_function,find_imf_angle_function,
     |      find_sint_function
        ! 
        !   
        use E_NBZ_model_mod, only:  E_NBZ_number_of_MLT_coeffs,
     |     E_NBZ_number_of_Btrans_coeffs,
     |     E_NBZ_number_of_imf_angle_coeffs,
     |     E_NBZ_number_of_sinT_coeffs,value_dummy, AO => E_NBZ_AO,
     |     Q => E_NBZ_Q, DQ => E_NBZ_DQ,
     |     fit_by_Btrans => E_NBZ_fit_by_Btrans,fit_by_imf_angle => 
     |     E_NBZ_fit_by_imf_angle,fit_by_sinT => E_NBZ_fit_by_sinT
        !   
        use qgen_mod, only: RI,ST,IBM,IEM,mcoef,dth,ns
        ! 
        implicit none
        ! 
        real, intent(in) :: mlat_glob,mlt,sinT,imf_angle
        real, intent(out) :: value_Ed1,value_Ed2,value
        logical, intent(out) :: value_found
        ! 
        integer :: ith,m,k,ix,mm,iBtrans,iIMF_angle,isinT,fit_count, 
     |        iMLT
        integer :: m2iMLT(E_NBZ_number_of_MLT_coeffs)
        data m2iMLT /0,-1,1,-2,2,-3,3/ 
        real :: mlat,poten,etheta,ephi,
     |       x,xm1,xp1,xp2,x0,sp,cp,sts,rst,qsoso,qsost,get_sinIm,
     |       Btrans_function,imf_angle_function,sint_function,sinIm
        real :: QS(-MMX:MMX,MXNMX),DQS(-MMX:MMX,MXNMX),F(-MMX:MMX)
        !real :: value
      
        ! General:
        logical, parameter :: debug = .false.
      
        !****************************************************************
        ! Interpolate to get coefficients for specified mlat:
        !****************************************************************
      
        value_found = .false.
      
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
        do iMLT = 1,E_NBZ_number_of_MLT_coeffs 
           m = m2iMLT(iMLT) 						     
           MM = IABS(M)
      
           Btrans_loop1:do iBtrans=1,1+(fit_by_Btrans
     |           *E_NBZ_number_of_Btrans_coeffs - fit_by_Btrans)
      
              imf_angle_loop1:do iimf_angle=1,1+(fit_by_imf_angle*
     |            E_NBZ_number_of_imf_angle_coeffs - fit_by_imf_angle)
                 imf_angle_function = 
     |               find_imf_angle_function(iimf_angle,imf_angle)   
      
      
                 sint_loop1:do isint=1,1+(fit_by_sint
     |                 *E_NBZ_number_of_sint_coeffs - fit_by_sint)
                    sint_function = find_sint_function (isint,sint)   
      
      	      do K=1,	MXNMX					     
                       IX = K - 1						
                       QS(M,K) = XM1*Q(2*IX+NS(MM),ITH-1)+
     |                   X0*Q(2*IX+NS(MM),ITH) 
     |                   +XP1*Q(2*IX+NS(MM),ITH+1)+
     |                   XP2*Q(2*IX+NS(MM),ITH+2)
                       DQS(M,K) = XM1*DQ(2*IX+NS(MM),ITH-1)+
     |                   X0*DQ(2*IX+NS(MM),ITH)
     |                    +XP1*DQ(2*IX+NS(MM),ITH+1)+
     |                   XP2*DQ(2*IX+NS(MM),ITH+2)   
                       IF (ITH.GT.0) THEN
                          qsost = QS(M,K)/STS
                       ELSE
                          qsost = DQS(M,K)
                       ENDIF
                       poten  = poten  + 
     |                   AO(M,K,iBtrans,iimf_angle,isint)*QS(M,K)*F(M)*
     |                     imf_angle_function*sint_function 
                       etheta = etheta+AO(M,K,iBtrans,iimf_angle,isint)*
     |                   DQS(M,K)*F(M)/RI*
     |                   imf_angle_function*sint_function 
                       ephi   = ephi- AO(M,K,iBtrans,iimf_angle,isint)
     |                   *qsost*M*F(-M)/RI*
     |                   imf_angle_function*sint_function 
                    enddo
                 end do sint_loop1
              end do imf_angle_loop1
           end do Btrans_loop1
        enddo
      
        !  (etheta equatorward, ephi eastward)
        sinIm     =  get_sinIm (mlat)
        value     =  poten*0.001              ! conversion from mV to V 
        value_Ed1 =  ephi                    ! mV/m
        value_Ed2 = -etheta/sinIm            ! mV/m  colatitude direction
      
        value_found = .true.
        value_dummy = .false.
      
      
      end subroutine E_NBZ_calc_model_value
      
      !---------------------------------------------------------------------------
      !---------------------------------------------------------------------------
      
      subroutine E_NBZ_set_coeff_array
      
        use common_model_mod, only: nmlat,average_mlat
        use E_NBZ_coeffmerge_mod, only: Ed1_sq2,Ed2_sq2,
     |        E_NBZ_fill_coeff_arrays
        use E_NBZ_model_mod, only: nfit => E_NBZ_nfit,
     |       b_Ed1 => E_NBZ_b_Ed1,b_Ed2 => E_NBZ_b_Ed2,b_interp_Ed1,
     |       b_interp_Ed2,b_90_Ed1,b_90_Ed2,x,imlat,ifit,
     |       E_NBZ_number_of_MLT_coeffs,
     |       E_NBZ_number_of_Btrans_coeffs,
     |       E_NBZ_number_of_imf_angle_coeffs,
     |       E_NBZ_number_of_sinT_coeffs,
     |       fit_by_MLT,pi,
     |        fit_by_Btrans => E_NBZ_fit_by_Btrans, 
     |        fit_by_imf_angle => E_NBZ_fit_by_imf_angle, 
     |        fit_by_sinT => E_NBZ_fit_by_sinT 
      
        integer :: binary_fit_by,check_dim
        integer :: iMLT,iBtrans,iimf_angle,isinT,fit_count
        real :: fac_ed1, fac_ed2, colat
        real :: fac_Y,Y
      
        !****************************************************************
        ! Fill parameter arrays:
        !****************************************************************
        call E_NBZ_fill_coeff_arrays
      
        !****************************************************************
        ! Determine array dimensions.  Arrays for data type 'potential'
        ! are the same size as arrays for other data types.
        !****************************************************************
      
        nmlat = size(Ed1_sq2,1)
        nfit  = size(Ed1_sq2,2)
      
        ! check dimensions
        check_dim = E_NBZ_number_of_MLT_coeffs* 
     |        E_NBZ_number_of_Btrans_coeffs*
     |        E_NBZ_number_of_imf_angle_coeffs* 
     |        E_NBZ_number_of_sinT_coeffs
        if(check_dim /= nfit) then
           write(6,'("NBZ_choose_coeff_array: wrong dimensions ", 2i4)') 
     |            nfit,check_dim
           stop
        endif
        !****************************************************************
        ! Allocate arrays into which coefficients will be transferred:
        !****************************************************************
        allocate(b_Ed1(nmlat,nfit),b_Ed2(nmlat,nfit), 
     |       b_interp_Ed1(nfit),
     |        b_interp_Ed2(nfit),x(nfit),
     |        b_90_Ed1(nfit),b_90_Ed2(nfit))
      
        !****************************************************************
        ! Transfer values:
        !****************************************************************
      
        do imlat=1,nmlat
           do ifit=1,nfit
              if(abs(Ed1_sq2(imlat,ifit)) >= 0.5000E+17) 
     |           Ed1_sq2(imlat,ifit)=0.
              if(abs(Ed2_sq2(imlat,ifit)) >= 0.5000E+17)
     |           Ed2_sq2(imlat,ifit)=0.
      
              b_Ed1(imlat,ifit) = Ed1_sq2(imlat,ifit)/2**0.5
              b_Ed2(imlat,ifit) = -Ed2_sq2(imlat,ifit)/2**0.5
           end do
        end do
      
      
       end subroutine E_NBZ_set_coeff_array
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

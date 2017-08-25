      subroutine calc_pot_S(mlt,mlat,sinT,imf_angle,Btrans,
     |  first_time,value_poynt,value_poten,value_S_background,
     |  value_Ed1,value_Ed2)

      use common_model_mod, only: pi,nq
      use B_bmed_model_mod, only: 
     | B_Bmed_number_of_Btrans_coeffs,
     | B_Bmed_number_of_imf_angle_coeffs,
     | B_Bmed_number_of_sinT_coeffs,B_Bmed_AO,B_Bmed_Q,B_Bmed_DQ,
     | B_Bmed_number_of_MLT_coeffs,B_Bmed_fit_by_Btrans,
     | B_Bmed_fit_by_imf_angle,B_Bmed_fit_by_sinT,
     | B_Bmed_b_Ed1,B_Bmed_b_Ed2,
     | B_Bmed_nfit,Bmed_N_AO,Bmed_FAC_AO
      use B_NBZ_model_mod, only: B_NBZ_number_of_Btrans_coeffs,
     | B_NBZ_number_of_imf_angle_coeffs,B_NBZ_AO,B_NBZ_Q,B_NBZ_DQ,
     | B_NBZ_number_of_MLT_coeffs,B_NBZ_fit_by_Btrans,
     | B_NBZ_fit_by_imf_angle,B_NBZ_fit_by_sinT,
     | B_NBZ_b_Ed1,B_NBZ_b_Ed2,
     | B_NBZ_nfit,NBZ_N_AO,NBZ_FAC_AO,B_NBZ_number_of_sinT_coeffs
      use B_bmax_model_mod, only: B_Bmax_number_of_Btrans_coeffs,
     | B_Bmax_number_of_imf_angle_coeffs,B_Bmax_AO,B_Bmax_Q,B_Bmax_DQ,
     | B_Bmax_number_of_MLT_coeffs,B_Bmax_fit_by_Btrans,
     | B_Bmax_fit_by_imf_angle,B_Bmax_fit_by_sinT,
     | B_Bmax_b_Ed1,B_Bmax_b_Ed2,
     | B_Bmax_nfit,Bmax_N_AO,B_Bmax_number_of_sinT_coeffs,Bmax_FAC_AO
      use B_b0_model_mod, only: B_B0_number_of_Btrans_coeffs,
     | B_B0_number_of_imf_angle_coeffs,B_B0_AO,B_B0_Q,B_B0_DQ,
     | B_B0_number_of_MLT_coeffs,B_B0_fit_by_Btrans,
     | B_B0_fit_by_imf_angle,
     | B_B0_fit_by_sinT,B_B0_b_Ed1,B_B0_b_Ed2,B_B0_nfit,
     | B_B0_nfit,B0_N_AO,
     | B_B0_number_of_sinT_coeffs,B0_FAC_AO
      use E_Bmed_model_mod, only: E_Bmed_number_of_Btrans_coeffs,
     | E_Bmed_number_of_imf_angle_coeffs,
     | E_Bmed_number_of_sinT_coeffs,
     | E_Bmed_number_of_MLT_coeffs,E_Bmed_AO,
     | E_Bmed_Q,E_Bmed_DQ,E_Bmed_fit_by_Btrans,E_Bmed_fit_by_imf_angle,
     | E_Bmed_fit_by_sinT,E_Bmed_b_Ed1,E_Bmed_b_Ed2,E_Bmed_nfit
      use E_NBZ_model_mod, only: E_NBZ_number_of_Btrans_coeffs,
     | E_NBZ_number_of_imf_angle_coeffs,E_NBZ_Q,E_NBZ_DQ,E_NBZ_AO,
     | E_NBZ_number_of_sinT_coeffs,E_NBZ_number_of_MLT_coeffs,
     | E_NBZ_fit_by_Btrans,E_NBZ_fit_by_imf_angle,E_NBZ_fit_by_sinT,
     | E_NBZ_b_Ed1,E_NBZ_b_Ed2,E_NBZ_nfit
      use E_bmax_model_mod, only: E_Bmax_number_of_Btrans_coeffs,
     | E_Bmax_number_of_imf_angle_coeffs,E_Bmax_AO,E_Bmax_Q,E_Bmax_DQ,
     | E_Bmax_number_of_sinT_coeffs,E_Bmax_number_of_MLT_coeffs,
     | E_Bmax_fit_by_Btrans,E_Bmax_fit_by_imf_angle,E_Bmax_fit_by_sinT,
     | E_Bmax_b_Ed1,E_Bmax_b_Ed2,E_Bmax_nfit
      use E_b0_model_mod, only: E_B0_number_of_Btrans_coeffs,
     | E_B0_number_of_imf_angle_coeffs,E_B0_AO,E_B0_Q,E_B0_DQ,
     | E_B0_number_of_MLT_coeffs,E_B0_fit_by_Btrans,
     | E_B0_fit_by_imf_angle,
     | E_B0_fit_by_sinT,E_B0_b_Ed1,E_B0_b_Ed2,E_B0_nfit,
     | E_B0_number_of_sinT_coeffs
      use qgen_mod, only: gen_coeffs
      use calc_ao_mod
      use S1_Bmed_coeff_mod, only: average_mlat

      implicit none
      real :: n_ao_dummy(nq)
      real,intent(in) :: mlt,mlat,Btrans,sinT,imf_angle
      logical,intent(inout) :: first_time
      real,intent(out) :: value_poynt,value_poten,
     |   value_S_background,value_Ed1,value_Ed2
     
      real, parameter :: zero = 0.,ri=6481200 ! R_e(6371 km) + 110 km
      real :: NBZ_value_Ed1,NBZ_value_Ed2,NBZ_value_dB1,NBZ_value_dB2,
     |  Bmax_value_Ed1,Bmax_value_Ed2,Bmax_value_dB1,Bmax_value_dB2,
     |  B0_value_Ed1,B0_value_Ed2,B0_value_dB1,B0_value_dB2,
     |  Bmed_value_Ed1,Bmed_value_Ed2,Bmed_value_dB1,Bmed_value_dB2,
     |  Bmax_MagPot,B0_MagPot,Bmax_ElPot,B0_ElPot
      real :: value_max, value_min,max_lat,min_lat,min_MLT,max_MLT,
     |  Bmed_value_S,Bmax_value_S, B0_value_S,Bmed_ElPot,Bmed_MagPot, 
     |  NBZ_MagPot,NBZ_ElPot,value_MagPot,value_ElPot,
     |  SDEd_Bmed_value1,SDEd_Bmed_value2,SDEd_Bmax_value1,SDEd_Bmax_value2,
     |  SDEd_B0_value1,SDEd_B0_value2,SDEd_NBZ_value1,SDEd_NBZ_value2,
     |  value_SDEd1,value_SDEd2
      logical :: value_found ,NBZ_value_found ,Bmax_value_found ,
     |  B0_value_found, 
     |   S_value_found,SDEd_Bmed_value_found,SDEd_Bmax_value_found,
     |   SDEd_B0_value_found,SDEd_NBZ_value_found

      real,parameter :: E_Bt_max_avg = 14.35
      real,parameter :: E_Bt_NBZ_avg = 6.76
      real,parameter :: E_Bt_med_avg = 4.25
      real,parameter :: E_Bt_0_avg = 4.33
      real,parameter :: B_Bt_max_avg = 14.5
      real,parameter :: B_Bt_NBZ_avg = 6.39
      real,parameter :: B_Bt_med_avg = 4.11
      real,parameter :: B_Bt_0_avg =  3.644
      real,parameter :: Bt_med_avg =  4.186
      real,parameter :: Bt_0_avg   =  2.91
      real,parameter :: Bt_max_avg = 11.47
      real, parameter :: beta_inf = 16.
      real ::beta_max_avg,beta_med_avg,E_beta_max_fac,E_beta_NBZ_fac,
     | E_beta_med_avg,beta_act,beta_NBZ_avg,B_beta_max_fac,beta_NBZ_FAC,
     | B_beta_NBZ_fac,B_beta_med_avg,beta_max_fac
      real :: value,value_S
      real :: value_dB1,value_dB2,get_beta
      real :: Bmed_value_FAC,Bmax_value_FAC,B0_value_FAC,NBZ_value_FAC,
     | value_FAC
      save beta_max_avg,beta_med_avg,E_beta_max_fac,E_beta_NBZ_fac,
     | E_beta_med_avg,beta_act,beta_NBZ_avg,B_beta_max_fac,beta_NBZ_FAC,
     | B_beta_NBZ_fac,B_beta_med_avg,beta_max_fac
!
      if(first_time) then
   	beta_max_avg = get_beta(Bt_max_avg) 
   	beta_max_fac = beta_inf/beta_max_avg 
   	beta_max_avg = get_beta(B_Bt_max_avg) 
   	B_beta_max_fac = beta_inf/beta_max_avg
   	beta_NBZ_avg = get_beta(B_Bt_NBZ_avg) 
   	B_beta_NBZ_fac = beta_inf/beta_NBZ_avg
   	B_beta_med_avg = get_beta(B_Bt_med_avg) 
   	beta_med_avg = get_beta(Bt_med_avg) 
   	beta_max_avg = get_beta(E_Bt_max_avg) 
   	E_beta_max_fac = beta_inf/beta_max_avg
   	beta_NBZ_avg = get_beta(E_Bt_NBZ_avg) 
   	E_beta_NBZ_fac = beta_inf/beta_NBZ_avg
   	E_beta_med_avg = get_beta(E_Bt_med_avg) 
   	beta_max_avg = get_beta(Bt_max_avg) 
   	beta_max_fac = beta_inf/beta_max_avg 

! calculate the Q coefficients
    	call gen_coeffs(E_Bmed_Q,E_Bmed_DQ,n_ao_dummy)
    	call gen_coeffs(E_NBZ_Q,E_NBZ_DQ,n_ao_dummy)
    	call gen_coeffs(E_Bmax_Q,E_Bmax_DQ,n_ao_dummy)
    	call gen_coeffs(E_B0_Q,E_B0_DQ,n_ao_dummy)
    	call gen_coeffs(B_Bmed_Q,B_Bmed_DQ,Bmed_N_AO)
    	call gen_coeffs(B_NBZ_Q,B_NBZ_DQ,NBZ_N_AO)
    	call gen_coeffs(B_Bmax_Q,B_Bmax_DQ,Bmax_N_AO)
    	call gen_coeffs(B_B0_Q,B_B0_DQ,B0_N_AO)

! set E & B coefficient array
    	call E_bmed_set_coeff_array
    	call E_NBZ_set_coeff_array
    	call E_bmax_set_coeff_array
    	call E_b0_set_coeff_array
    	call B_bmed_set_coeff_array
    	call B_NBZ_set_coeff_array
    	call B_bmax_set_coeff_array
    	call B_b0_set_coeff_array
    	call S1_Bmed_choose_coeff_array
    	call S1_Bmax_choose_coeff_array
    	call S1_B0_choose_coeff_array
	
! calculate fitting coefficients ao
        call calc_ao_coeff(E_Bmed_number_of_MLT_coeffs,
     |   E_Bmed_number_of_Btrans_coeffs,
     |   E_Bmed_number_of_imf_angle_coeffs,E_Bmed_number_of_sinT_coeffs,
     |   E_Bmed_fit_by_Btrans,E_Bmed_fit_by_imf_angle,
     |   E_Bmed_fit_by_sinT,E_Bmed_nfit,E_Bmed_b_Ed1,E_Bmed_b_Ed2,
     |   E_Bmed_Q,E_Bmed_DQ,
     |   E_Bmed_AO)
        call calc_ao_coeff(E_NBZ_number_of_MLT_coeffs,
     |   E_NBZ_number_of_Btrans_coeffs,
     |   E_NBZ_number_of_imf_angle_coeffs,E_NBZ_number_of_sinT_coeffs,
     |   E_NBZ_fit_by_Btrans,E_NBZ_fit_by_imf_angle,
     |   E_NBZ_fit_by_sinT,E_NBZ_nfit,E_NBZ_b_Ed1,E_NBZ_b_Ed2,
     |   E_NBZ_Q,E_NBZ_DQ,
     |   E_NBZ_AO)
        call calc_ao_coeff(B_Bmed_number_of_MLT_coeffs,
     |   B_Bmed_number_of_Btrans_coeffs,
     |   B_Bmed_number_of_imf_angle_coeffs,B_Bmed_number_of_sinT_coeffs,
     |   B_Bmed_fit_by_Btrans,B_Bmed_fit_by_imf_angle,
     |   B_Bmed_fit_by_sinT,B_Bmed_nfit,B_Bmed_b_Ed1,B_Bmed_b_Ed2,
     |   B_Bmed_Q,E_Bmed_DQ,
     |   B_Bmed_AO,
     |   Bmed_N_AO,Bmed_FAC_AO)
        call calc_ao_coeff(B_NBZ_number_of_MLT_coeffs,
     |   B_NBZ_number_of_Btrans_coeffs,
     |   B_NBZ_number_of_imf_angle_coeffs,B_NBZ_number_of_sinT_coeffs,
     |   B_NBZ_fit_by_Btrans,B_NBZ_fit_by_imf_angle,
     |   B_NBZ_fit_by_sinT,B_NBZ_nfit,B_NBZ_b_Ed1,B_NBZ_b_Ed2,
     |   B_NBZ_Q,E_NBZ_DQ,
     |   B_NBZ_AO,
     |   NBZ_N_AO,NBZ_FAC_AO)
        call calc_ao_coeff(E_Bmax_number_of_MLT_coeffs,
     |   E_Bmax_number_of_Btrans_coeffs,
     |   E_Bmax_number_of_imf_angle_coeffs,E_Bmax_number_of_sinT_coeffs,
     |   E_Bmax_fit_by_Btrans,E_Bmax_fit_by_imf_angle,
     |   E_Bmax_fit_by_sinT,E_Bmax_nfit,E_Bmax_b_Ed1,E_Bmax_b_Ed2,
     |   E_Bmax_Q,E_Bmax_DQ,
     |   E_Bmax_AO)
        call calc_ao_coeff(E_B0_number_of_MLT_coeffs,
     |   E_B0_number_of_Btrans_coeffs,
     |   E_B0_number_of_imf_angle_coeffs,E_B0_number_of_sinT_coeffs,
     |   E_B0_fit_by_Btrans,E_B0_fit_by_imf_angle,
     |   E_B0_fit_by_sinT,E_B0_nfit,E_B0_b_Ed1,E_B0_b_Ed2,
     |   E_B0_Q,E_B0_DQ,
     |   E_B0_AO)
        call calc_ao_coeff(B_Bmax_number_of_MLT_coeffs,
     |   B_Bmax_number_of_Btrans_coeffs,
     |   B_Bmax_number_of_imf_angle_coeffs,B_Bmax_number_of_sinT_coeffs,
     |   B_Bmax_fit_by_Btrans,B_Bmax_fit_by_imf_angle,
     |   B_Bmax_fit_by_sinT,B_Bmax_nfit,B_Bmax_b_Ed1,B_Bmax_b_Ed2,
     |   B_Bmax_Q,E_Bmax_DQ,
     |   B_Bmax_AO,
     |   Bmax_N_AO,Bmax_FAC_AO)
        call calc_ao_coeff(B_B0_number_of_MLT_coeffs,
     |   B_B0_number_of_Btrans_coeffs,
     |   B_B0_number_of_imf_angle_coeffs,B_B0_number_of_sinT_coeffs,
     |   B_B0_fit_by_Btrans,B_B0_fit_by_imf_angle,
     |   B_B0_fit_by_sinT,B_B0_nfit,B_B0_b_Ed1,B_B0_b_Ed2,
     |   B_B0_Q,E_B0_DQ,
     |   B_B0_AO,
     |   B0_N_AO,B0_FAC_AO)
 
         first_time =.false.
      endif  ! first_time
  
      beta_act = get_beta(Btrans)

      ! B
      !write(6,*) 'IMF_ANGLE=',imf_angle
      call B_bmed_calc_model_value(mlat,mlt,sinT,imf_angle,
     |      Bmed_value_dB1,Bmed_value_dB2,Bmed_MagPot,value_found,
     |      Bmed_value_FAC)
      call B_NBZ_calc_model_value(mlat,mlt,sinT,imf_angle,
     |     NBZ_value_dB1,NBZ_value_dB2,NBZ_MagPot,NBZ_value_found,
     |     NBZ_value_FAC)
      Bmed_value_dB1 = Bmed_value_dB1  + NBZ_value_dB1
      Bmed_value_dB2 = Bmed_value_dB2  + NBZ_value_dB2
      Bmed_MagPot = Bmed_MagPot + NBZ_MagPot
      Bmed_value_FAC = Bmed_value_FAC+NBZ_value_FAC

      if(value_found == .false.) then
         write(6,*) 'value not found ',mlat,MLT,abs(value)
         stop
      endif
      if(NBZ_value_found == .false.) then
         write(6,*) 'NBZ_value not found ',mlat,MLT,abs(NBZ_value_dB1)
         stop
      endif


      if(Btrans.gt.B_Bt_med_avg) then
         ! calculating value for B_max model 
         call B_bmax_calc_model_value(mlat,mlt,sinT,imf_angle,
     |      Bmax_value_dB1,Bmax_value_dB2,Bmax_MagPot,
     |      Bmax_value_found,Bmax_value_FAC)

         Bmax_value_dB1 = Bmax_value_dB1*B_beta_max_fac
         Bmax_value_dB2 = Bmax_value_dB2*B_beta_max_fac
         Bmax_MagPot = Bmax_MagPot*B_beta_max_fac
         Bmax_value_FAC = Bmax_value_FAC*B_beta_max_fac
         if(Bmax_value_found == .false.) then
            write(6,*) 'Bmax_value not found ',mlat,
     |      MLT,abs(Bmax_value_dB1)
            stop
         endif

         NBZ_value_dB1 = NBZ_value_dB1*B_beta_NBZ_fac  ! extrapolating to Bt
         NBZ_value_dB2 = NBZ_value_dB2*B_beta_NBZ_fac  ! extrapolating to Bt
         NBZ_MagPot    = NBZ_MagPot*B_beta_NBZ_fac

         NBZ_value_FAC = NBZ_value_FAC*B_beta_NBZ_fac
         if(NBZ_value_found == .false.) then
            write(6,*)'NBZ_value not found ',mlat,MLT,abs(NBZ_value_dB1)
            stop
         endif

 
         if(imf_angle < pi/2.or. 
     |         imf_angle > 3*pi/2) then
            Bmax_MagPot = Bmax_MagPot*(1-cos(imf_angle)**2)
            Bmax_MagPot = Bmax_MagPot+ NBZ_MagPot*(cos(imf_angle)**2)
            Bmax_value_dB1 = Bmax_value_dB1*(1-cos(imf_angle)**2)
            Bmax_value_dB1 = Bmax_value_dB1+ NBZ_value_dB1*
     |        (cos(imf_angle)**2)
            Bmax_value_dB2 = Bmax_value_dB2*(1-cos(imf_angle)**2)
            Bmax_value_dB2 = Bmax_value_dB2+ NBZ_value_dB2*
     |        (cos(imf_angle)**2)
            Bmax_value_FAC = Bmax_value_FAC*(1-cos(imf_angle)**2)
            Bmax_value_FAC = Bmax_value_FAC+ NBZ_value_FAC*
     |        (cos(imf_angle)**2)
         endif

         ! interpolating between B_max and Bmed model
         value_dB1 = Bmax_value_dB1*(beta_act-B_beta_med_avg)+
     |         Bmed_value_dB1*(beta_inf-beta_act) 
         value_dB1 = value_dB1/(beta_inf-B_beta_med_avg) 
         value_dB2 = Bmax_value_dB2*(beta_act-B_beta_med_avg)+
     |         Bmed_value_dB2*(beta_inf-beta_act) 
         value_dB2 = value_dB2/(beta_inf-B_beta_med_avg) 
         value_MagPot = Bmax_MagPot*(beta_act-B_beta_med_avg)+
     |         Bmed_MagPot*(beta_inf-beta_act)
         value_MagPot = value_MagPot/(beta_inf-B_beta_med_avg)
         value_FAC = Bmax_value_FAC*(beta_act-B_beta_med_avg)+
     |         Bmed_value_FAC*(beta_inf-beta_act)
         value_FAC = value_FAC/(beta_inf-B_beta_med_avg)
      else  
         call B_b0_calc_model_value(mlat,mlt,sinT,imf_angle,
     |         0.,B0_value_dB1,B0_value_dB2,B0_MagPot,B0_value_found,
     |         B0_value_FAC) ! use Btrans = 0 instead of Btrans

         ! interpolating between B_0 and Bmed model
         value_dB1 = Bmed_value_dB1*(beta_act)+B0_value_dB1*
     |         (B_beta_med_avg-beta_act) 
         value_dB1 = value_dB1/(B_beta_med_avg) 
         value_dB2 = Bmed_value_dB2*(beta_act)+B0_value_dB2*
     |         (B_beta_med_avg-beta_act) 
         value_dB2 = value_dB2/(B_beta_med_avg)
         value_MagPot = Bmed_MagPot*(beta_act)+B0_MagPot*
     |         (B_beta_med_avg-beta_act)
         value_MagPot = value_MagPot/(B_beta_med_avg)
         value_FAC = Bmed_value_FAC*(beta_act)+B0_value_FAC*
     |         (B_beta_med_avg-beta_act)
         value_FAC = value_FAC/(B_beta_med_avg)
         if(B0_value_found == .false.) then
            write(6,*) 'B00_value not found ',mlat,MLT,abs(B0_value_dB1)
            stop
         endif

      endif

      !E
      !write(6,*) 'IMF_ANGLE_E=',imf_angle
      call E_bmed_calc_model_value(mlat,mlt,sinT,imf_angle,
     |      Bmed_value_Ed1,Bmed_value_Ed2,Bmed_ElPot,value_found)
      call E_NBZ_calc_model_value(mlat,mlt,sinT,imf_angle,
     |      NBZ_value_Ed1,NBZ_value_Ed2,NBZ_ElPot,
     |      NBZ_value_found)
      Bmed_value_Ed1 = Bmed_value_Ed1  + NBZ_value_Ed1
      Bmed_value_Ed2 = Bmed_value_Ed2  + NBZ_value_Ed2
      !write(6,*) 'Bmed_ElPot=',Bmed_ElPot
      !write(6,*) 'NBZ_ElPot=',NBZ_ElPot
      Bmed_ElPot = Bmed_ElPot + NBZ_ElPot
      
      if(value_found == .false.) then
         write(6,*) 'value not found ',mlat,MLT,abs(value_Ed1)
         stop
      endif
      if(NBZ_value_found == .false.) then
         write(6,*) 'NBZ_value not found ',mlat,MLT,abs(NBZ_value_Ed1)
         stop
      endif

      if(Btrans.gt.E_Bt_med_avg) then
         ! calcualting value for B_max model
         call E_bmax_calc_model_value(mlat,mlt,sinT,imf_angle,
     |         Bmax_value_Ed1,Bmax_value_Ed2,Bmax_ElPot,
     |         Bmax_value_found)
         Bmax_value_Ed1 = Bmax_value_Ed1*E_beta_max_fac
         Bmax_value_Ed2 = Bmax_value_Ed2*E_beta_max_fac
         Bmax_ElPot = Bmax_ElPot*E_beta_max_fac
         !write(6,*) 'Bmax_ElPot=',Bmax_ElPot,E_beta_max_fac
         if(Bmax_value_found == .false.) then
            write(6,*) 'Bmax_value not found ',mlat,MLT,
     |        abs(Bmax_value_Ed1)
            stop
         endif

         NBZ_value_Ed1 = NBZ_value_Ed1*E_beta_NBZ_fac  ! extrapolating to Bmax (beta=16)
         NBZ_value_Ed2 = NBZ_value_Ed2*E_beta_NBZ_fac  ! extrapolating to Bmax (beta=16)
         NBZ_ElPot     = NBZ_ElPot*E_beta_NBZ_fac
         !write(6,*) 'NBZ_ElPot=',NBZ_ElPot*E_beta_NBZ_fac
         if(NBZ_value_found == .false.) then
            write(6,*)'NBZ_value not found ',mlat,MLT,abs(NBZ_value_Ed1)
            stop
         endif

         if(imf_angle < pi/2.or. 
     |        imf_angle > 3*pi/2) then
            Bmax_ElPot     = Bmax_ElPot*(1-cos(imf_angle)**2)
            Bmax_ElPot     = Bmax_ElPot+ NBZ_ElPot*(cos(imf_angle)**2)
            Bmax_value_Ed1 = Bmax_value_Ed1*(1-cos(imf_angle)**2)
            Bmax_value_Ed1 = Bmax_value_Ed1+ NBZ_value_Ed1*
     |        (cos(imf_angle)**2)
            Bmax_value_Ed2 = Bmax_value_Ed2*(1-cos(imf_angle)**2)
            Bmax_value_Ed2 = Bmax_value_Ed2+ NBZ_value_Ed2*
     |        (cos(imf_angle)**2)
           !write(6,*) 'Bmax_ElPot_IMF=',Bmax_ElPot

         endif

         ! interpolating between B_max and Bmed model
         value_Ed1 = Bmax_value_Ed1*(beta_act-E_beta_med_avg)+
     |         Bmed_value_Ed1*(beta_inf-beta_act) 
         value_Ed1 = value_Ed1/(beta_inf-E_beta_med_avg) 
         value_Ed2 = Bmax_value_Ed2*(beta_act-E_beta_med_avg)+
     |         Bmed_value_Ed2*(beta_inf-beta_act) 
         value_Ed2 = value_Ed2/(beta_inf-E_beta_med_avg)
         value_ElPot = Bmax_ElPot*(beta_act-E_beta_med_avg)+
     |         Bmed_ElPot*(beta_inf-beta_act)
         value_ElPot = value_ElPot/(beta_inf-E_beta_med_avg)
          ! write(6,*) 'Inter_ElPot=',value_ElPot
      else
         call E_b0_calc_model_value(mlat,mlt,sinT,imf_angle,
     |         0.,B0_value_Ed1,B0_value_Ed2,B0_ElPot,
     |         B0_value_found) ! use Btrans = 0 instead of Btrans
          ! write(6,*) 'B0_ElPot=',B0_ElPot

         ! interpolating between B_0 and Bmed model
         value_Ed1 = Bmed_value_Ed1*(beta_act)+B0_value_Ed1*
     |         (E_beta_med_avg-beta_act) 
         value_Ed1 = value_Ed1/(E_beta_med_avg) 
         value_Ed2 = Bmed_value_Ed2*(beta_act)+B0_value_Ed2*
     |         (E_beta_med_avg-beta_act) 
         value_Ed2 = value_Ed2/(E_beta_med_avg)
         value_ElPot = Bmed_ElPot*(beta_act)+B0_ElPot*
     |         (E_beta_med_avg-beta_act)
         value_ElPot = value_ElPot/(E_beta_med_avg)
         !  write(6,*) 'Inter_ElPot=',value_ElPot

         if(B0_value_found == .false.) then
            write(6,*) 'B00_value not found ',mlat,MLT,abs(B0_value_Ed1)
            stop
         endif
      endif

      ! Poynting flux
      call S1_Bmed_calc_model_value(mlat,mlt,sinT,imf_angle,
     |      Bmed_value_S,S_value_found)
      if(S_value_found == .false.) then
         write(6,*) 'S_value not found ',mlat,MLT,abs(Bmed_value_S)
         stop
      endif

      if(Btrans.gt.Bt_med_avg) then
         ! calculating value for B_max model
         call S1_Bmax_calc_model_value(mlat,mlt,sinT,imf_angle,
     |         Bmax_value_S,Bmax_value_found)
         Bmax_value_S = Bmax_value_S*beta_max_fac
         if(Bmax_value_found == .false.) then
            write(6,*)'Bmax_value not found ',mlat,MLT,abs(Bmax_value_S)
            stop
         endif

         ! interpolating between B_max and Bmed model
         value_S = Bmax_value_S*(beta_act*beta_act-
     |         beta_med_avg*beta_med_avg)+ 
     |         Bmed_value_S*(beta_inf*beta_inf-
     |         beta_act*beta_act) 
         value_S = value_S/(beta_inf*beta_inf-beta_med_avg*beta_med_avg)
      else
         call S1_b0_calc_model_value(mlat,mlt,sinT,imf_angle,
     |         0.,B0_value_S,B0_value_found) ! use Btrans = 0 instead of Btrans

         ! interpolating between B_0 and Bmed model
         value_S = Bmed_value_S*(beta_act)+B0_value_S*
     |         (beta_med_avg-beta_act) 
         value_S = value_S/(beta_med_avg)

         if(B0_value_found == .false.) then
            write(6,*) 'B00_value not found ',mlat,MLT,abs(B0_value_S)
            stop
         endif
      endif

      !      = (mV/m)*(nT)/(4 pi e-7 T m/A) (W/VA)(V/1000mV)(T/1e9nT)
      !      = 1e-5/4pi W/m^2
      value_S_background = value_Ed1*value_dB2-value_Ed2*value_dB1
      value_S_background = 1.e-5/(4.*pi) *value_S_background
      value_poynt	 = value_S_background + value_S
      
      value_poten = value_ElPot

      end subroutine calc_pot_S

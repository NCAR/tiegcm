      subroutine calc_SDEd(mlt,mlat,sinT,imf_angle,Btrans,
     |    first_time,value_SDEd1,value_SDEd2)

      use common_model_mod, only:pi

      implicit none

      real,intent(in) :: mlt,mlat,Btrans,sinT,imf_angle
      logical,intent(in) :: first_time
      real,intent(out) :: value_SDEd1,value_SDEd2
      real :: SDEd_Bmed_value1,SDEd_Bmed_value2,SDEd_Bmax_value1,
     |	    SDEd_Bmax_value2,SDEd_B0_value1,
     |	    SDEd_B0_value2,SDEd_NBZ_value1,SDEd_NBZ_value2
      logical :: SDEd_Bmed_value_found,SDEd_Bmax_value_found,
     |	      SDEd_B0_value_found,SDEd_NBZ_value_found

      real,parameter :: zero = 0.
      real,parameter :: SDEd_Bt_max_avg = 14.35
      real,parameter :: SDEd_Bt_NBZ_avg = 6.61
      real,parameter :: SDEd_Bt_med_avg = 4.26
      real,parameter :: SDEd_Bt_0_avg	= 1.5
      real,parameter :: beta_inf = 16.
      real,save ::beta_act,SDEd_beta_max_avg,
     |	 SDEd_beta_max_fac,SDEd_beta_NBZ_avg,SDEd_beta_NBZ_fac,
     |	 SDEd_beta_0,
     |	 SDEd_beta_med_avg
      real :: Bmed_value_FAC,Bmax_value_FAC,B0_value_FAC,
     |	 NBZ_value_FAC,
     |	 value_FAC
      real:: get_beta

      if(first_time) then
        SDEd_beta_max_avg = get_beta(SDEd_Bt_max_avg)
        SDEd_beta_max_fac = beta_inf/SDEd_beta_max_avg
        SDEd_beta_NBZ_avg = get_beta(SDEd_Bt_NBZ_avg)
        SDEd_beta_NBZ_fac = beta_inf/SDEd_beta_NBZ_avg
        SDEd_beta_0	  = get_beta(SDEd_Bt_0_avg)
        SDEd_beta_med_avg = get_beta(SDEd_Bt_med_avg)

        call SDEd_Bmed_choose_coeff_array
        call SDEd_NBZ_choose_coeff_array
        call SDEd_Bmax_choose_coeff_array
        call SDEd_B0_choose_coeff_array
      endif  ! first_time
      
      beta_act = get_beta(Btrans)

      call SDEd_Bmed_calc_model_value(mlat,mlt,sinT,imf_angle,Btrans,
     |      SDEd_Bmed_value1,SDEd_Bmed_value2,SDEd_Bmed_value_found)
      call SDEd_NBZ_calc_model_value(mlat,mlt,sinT,Btrans,imf_angle,
     |     SDEd_NBZ_value1,SDEd_NBZ_value2,SDEd_NBZ_value_found)
      if( SDEd_Bmed_value_found == .false. .or. SDEd_NBZ_value_found
     |     == .false.) then
         write(6,*) 'SDEd_value not found ',mlat,MLT,
     |     abs(SDEd_Bmed_value1)
         stop
      endif

      if(SDEd_NBZ_value1 < 0) call fill_neg_value(SDEd_NBZ_value1,
     |      mlat)
      if(SDEd_NBZ_value2 < 0) call fill_neg_value(SDEd_NBZ_value2,
     |      mlat)
      SDEd_Bmed_value1 = SDEd_Bmed_value1 + SDEd_NBZ_value1
      SDEd_Bmed_value2 = SDEd_Bmed_value2 + SDEd_NBZ_value2

      if(SDEd_Bmed_value1 < 0) call fill_neg_value(SDEd_Bmed_value1,
     |      mlat)
      if(SDEd_Bmed_value2 < 0) call fill_neg_value(SDEd_Bmed_value2,
     |      mlat)

      if(Btrans.gt.SDEd_Bt_med_avg) then
         ! calculating value for B_max model
         call SDEd_Bmax_calc_model_value(mlat,mlt,sinT,imf_angle,
     |        Btrans,SDEd_Bmax_value1,SDEd_Bmax_value2,
     |        SDEd_Bmax_value_found)
         if (SDEd_Bmax_value1 < 0) call fill_neg_value(
     |         SDEd_Bmax_value1,mlat)
         if (SDEd_Bmax_value2 < 0) call fill_neg_value(
     |         SDEd_Bmax_value2,mlat)
         SDEd_Bmax_value1 = SDEd_Bmax_value1*SDEd_beta_max_fac
         SDEd_Bmax_value2 = SDEd_Bmax_value2*SDEd_beta_max_fac
         if(SDEd_Bmax_value_found == .false.) then
            write(6,*) 'SDEd_Bmax_value not found ',mlat,MLT,
     |   	 abs(SDEd_Bmax_value1),abs(SDEd_Bmax_value2)
            stop
         endif

         !extrapolating to Bt
         SDEd_NBZ_value1 = SDEd_NBZ_value1*SDEd_beta_NBZ_fac  
         SDEd_NBZ_value2 = SDEd_NBZ_value2*SDEd_beta_NBZ_fac  

         if(imf_angle < pi/2.or.imf_angle > 3*pi/2) then
            SDEd_Bmax_value1 = SDEd_Bmax_value1*(1-cos(imf_angle)**2)
            SDEd_Bmax_value1 = SDEd_Bmax_value1+ SDEd_NBZ_value1*
     |       (cos(imf_angle)**2)
            SDEd_Bmax_value2 = SDEd_Bmax_value2*(1-cos(imf_angle)**2)
            SDEd_Bmax_value2 = SDEd_Bmax_value2+ SDEd_NBZ_value2*
     |       (cos(imf_angle)**2)
         endif

         ! interpolating between Bmax and Bmed model
         value_SDEd1 = SDEd_Bmax_value1*(beta_act-
     |        SDEd_beta_med_avg)+SDEd_Bmed_value1*
     |        (beta_inf-beta_act)
         value_SDEd1 = value_SDEd1/(beta_inf-SDEd_beta_med_avg)
         value_SDEd2 = SDEd_Bmax_value2*(beta_act-
     |         SDEd_beta_med_avg)+SDEd_Bmed_value2*
     |         (beta_inf-beta_act)
         value_SDEd2 = value_SDEd2/(beta_inf-SDEd_beta_med_avg)

      else
         call SDEd_B0_calc_model_value(mlat,mlt,sinT,imf_angle,
     |         0.,SDEd_B0_value1,SDEd_B0_value2,SDEd_B0_value_found)

         if(SDEd_B0_value1 < 0) call fill_neg_value(SDEd_B0_value1,
     |         mlat)
         if(SDEd_B0_value2 < 0) call fill_neg_value(SDEd_B0_value2,
     |         mlat)

         ! interpolating between B_0 and Bmed model
         value_SDEd1 = SDEd_Bmed_value1*(beta_act)+SDEd_B0_value1*
     |         (SDEd_beta_med_avg-beta_act)
         value_SDEd1 = value_SDEd1/(SDEd_beta_med_avg)

         value_SDEd2 = SDEd_Bmed_value2*(beta_act)+SDEd_B0_value2*
     |         (SDEd_beta_med_avg-beta_act)
         value_SDEd2 = value_SDEd2/(SDEd_beta_med_avg)


         if(SDEd_B0_value_found == .false.) then
            write(6,*) 'SDEd_B0_value not found ',mlat,MLT,abs
     |   	 (SDEd_B0_value1),abs(SDEd_B0_value2)
            stop
         endif
      endif
      ! from variance into standard deviation
      if(value_SDEd1 >= 0) then
         value_SDEd1 = sqrt(value_SDEd1)
      endif
      if(value_SDEd2 >= 0) then
         value_SDEd2 = sqrt(value_SDEd2)
      endif

      end subroutine calc_SDEd
!----------------------------------------------------------------------------
      subroutine fill_neg_value(value,mlat)
!
      implicit none
!
      real,intent(in):: mlat
      real,intent(inout):: value
!
      if(mlat .ge. 50. .and. mlat .lt.60) then
    	 value = (mlat - 50.)/10.*1.5+1
      elseif(mlat .ge. 60. .and. mlat .lt.70) then
    	 value = (mlat - 60.)/10.*7.5+2.5
      elseif(mlat .ge. 70.)  then
    	 value = 10.
      elseif(mlat .lt. 50.)  then
    	 value = 1.
      endif
      value = value*value
!
      end subroutine fill_neg_value
!----------------------------------------------------------------------------


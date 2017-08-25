module calc_ao_mod

contains

  subroutine calc_ao_coeff(number_of_MLT_coeffs,&
       &number_of_Btrans_coeffs,number_of_imf_angle_coeffs,&
       &number_of_sinT_coeffs,fit_by_Btrans,&
       &fit_by_imf_angle,fit_by_sinT,nfit,b_Ed1,b_Ed2,Q,DQ,&
       &AO,&
       &N_AO,FAC_AO)

    use common_model_mod, only: MMX,MXNMX,ITHMX,nmlat,average_mlat,&
         &mlat_resolution,NQ,pi

    use qgen_mod,only: inlats,dlatin,ibm,iem,nmx,&
         &mcoef,ri,dth,ns

    implicit none 

    integer, intent(in) :: number_of_MLT_coeffs,&
         &number_of_Btrans_coeffs,number_of_imf_angle_coeffs,&
         &number_of_sinT_coeffs,fit_by_Btrans,&
         &fit_by_imf_angle,fit_by_sinT,nfit

    real, intent(in) :: b_Ed1(1:nmlat,1:nfit),b_Ed2(1:nmlat,1:nfit),&
         &Q(NQ,0:ITHMX),DQ(NQ,0:ITHMX)

    real, intent(out) :: AO(-MMX:MMX,MXNMX,number_of_Btrans_coeffs,&
         &number_of_imf_angle_coeffs,number_of_sinT_coeffs)

    real, intent(in), optional :: N_AO(NQ)

    real, intent(out), optional :: FAC_AO(-MMX:MMX,MXNMX,number_of_Btrans_coeffs,&
         &number_of_imf_angle_coeffs,number_of_sinT_coeffs)
      
    !C ITHMX = number of uniform intervals between pole and equator.
    !C   ITHMX = 45 gives 2 degree latitude grid.
    !C ITHTRNS = index for colatitude where basis functions transition from
    !C   wavy (high latitudes) to decaying (middle and low latitudes).
    !C   This colatitude is (90 degrees)*ITHTRNS/ITHMX.
    !C   ITHTRNS = 18 corresponds to 36 colatitude (54 latitude) when ITHMX=45.
    !C MMX = maximum longitudinal wavenumber of Fourier harmonics 
    !C MXNMX = number of symmetric Q functions for a given longitudinal wavenumber.
    !C NQS = (MMX+1)*MXNMX = total number of symmetric Q functions.
    !C NQ = 2*NQS = total number of Q functions (antisymmetric functions are 
    !C   not used here).
    !C MCOEF = (2*MMX + 1)*MXNMX = total number of fitted coefficients.
    integer, parameter :: LONMX=150,ITHPLT=20

    !C LONMX = number of uniform intervals between in longitude.
    !C ITHPLT = index for outer colatitude for generating fields to plot.
    !C   This colatitude = (90 degrees)*ITHPLT/ITHMX. 
    !C   ITHPLT=20 is 40 degrees colatitude (50 degrees latitude) when ITHMX=45.
    !C

    integer :: i,lon,ith,n,m,mm,k,ix,imx,iBtrans,iIMF_angle,isinT,&
         &fit_count, iMLT
    integer :: m2iMLT(7)
    data m2iMLT /0,-1,1,-2,2,-3,3/ 
    real :: twopi,XLONMX,dlon,ph,dtheta,sth,rst,rdt,rdtst,x,&
         &  xm1,x0,xp1,xp2,rad,curtx,curpx,ed1,ed2
    real :: SL(LONMX),CL(LONMX),QS(-MMX:MMX,MXNMX),DQS(-MMX:MMX,MXNMX)
    real :: F(-MMX:MMX,LONMX), &
         & MI(MCOEF),CURT(LONMX),CURP(LONMX),CURTM(-MMX:MMX), &	      
         & CURPM(-MMX:MMX),MLT(0:LONMX),CLATD(0:ITHPLT)
    !C dlatin is latitude increment (degrees) of input electric fields.
    !C inlats is number of colatitudes where electric fields are non-negligible.

    real :: clatin(0:inlats)
    logical, parameter :: debug = .false.
    !C clatin contains colatitudes (radians) of input electric fields
    !C   (set below to 1./RAD, 3./RAD,...46./RAD);
    !C   Beyond 45 degrees colatitude, electric field values are assumed zero.
    !C             
    RAD = 180./pi                                                            
    TWOPI=2.*PI                                                              
    XLONMX=LONMX                                                              
    DLON = TWOPI/FLOAT(LONMX)                                                 
    do i=0,inlats
       clatin(i) = dlatin*(float(i) - .5)/RAD
    enddo
    MLT(0) = 0.                                                               
    DO 200 LON=1,LONMX                                                        
       MLT(LON) = LON*24/FLOAT(LONMX)                                            
       PH = LON*DLON                                                            
       SL(LON) = SIN(PH)                                                        
200    CL(LON) = COS(PH)                                                        
       DO 201 ITH=0,ITHPLT                                                       
201       CLATD(ITH) = ITH*90/FLOAT(ITHMX)                                         
          IBM(-MMX) = 1                                                             
          I = 0                                                                     
          DO 209 M=-MMX,MMX                                                         
             DO 207 N=1,NMX(IABS(M))                                                  
                I = I + 1                                                               
207             MI(I) = M                                                               
                IEM(M) = IBM(M) - 1 + NMX(IABS(M))                                       
209             IF (M.NE.MMX) IBM(M+1) = NMX(IABS(M)) + IBM(M)                           
                IMX = IEM(MMX)          
                CALL FCMP(MMX,LONMX,CL,SL,F)                                              

                !C COMPUTE COEFFICIENTS FOR MODEL ELECTRIC FIELD, AO                             
                AO = 0.         

                !C Vector-multiply electric field by RI**2 times negative gradient of 
                !C   each basis function, and integrate over hemisphere to get 
                !C   coefficients AO.  
                !C Assume that electric field = 0 at and beyond clat(inlats-1)
                !C LOOP IN COLATITUDE                                                           
                do 370 i=nmlat,1,-1
                   clatin(i) = (90.-average_mlat(i))*pi/180.
                   dtheta = mlat_resolution*pi/180.
                   sth = sin(clatin(i))
                   RST = RI*sth
                   RDT = RI*dtheta
                   RDTST = RDT*sth                                                   

                   !C Calculate multipliers for cubic interpolation between latitudes where
                   !C  Q, DQ are available
                   X = clatin(i)/DTH
                   ITH = X
                   ITH = MAX0(1,MIN0(ITHMX-2,ITH))
                   X = X - ITH
                   XM1 = X*(-2. + X*(3. - X))/6.
                   X0  = 1. + X*(-.5 + X*(-1. + .5*X))
                   XP1 = X*(1. + X*(.5 - .5*X))
                   XP2 = X*(-1. + X*X)/6.

                   DO 305 M=-MMX,MMX                                                        
                      MM = IABS(M)                                                            
                      !C ASSUME THAT THERE ARE EQUAL NUMBERS OF SYMMETRIC AND ANTISYMMETRIC            
                      !C  Q'S FOR A GIVEN VALUE OF M.  NS(MM) GIVES THE FIRST INDEX VALUE              
                      !C  FOR THE FIRST SUBSCRIPT OF Q FOR WAVENUMBER MM.                              
                      DO 303 K=1,  MXNMX                                               
                         IX = K - 1  
                         QS(M,K) = XM1*Q(2*IX+NS(MM),ITH-1)+X0*Q(2*IX+NS(MM),ITH)&
                              &    +XP1*Q(2*IX+NS(MM),ITH+1)+XP2*Q(2*IX+NS(MM),ITH+2)
                         DQS(M,K) = XM1*DQ(2*IX+NS(MM),ITH-1)+X0*DQ(2*IX+NS(MM),ITH)&
                              &    +XP1*DQ(2*IX+NS(MM),ITH+1)+XP2*DQ(2*IX+NS(MM),ITH+2)
303                      continue 
305                      CONTINUE     

                         DO 350 K=1,   MXNMX  
                            fit_count = 1
                            if(present(FAC_AO)) then
                               IX = K - 1
                            end if
                            DO 360 iMLT = 1,number_of_MLT_coeffs 
                               m = m2iMLT(iMLT)
                               !FAC
                               mm = abs(m)
                               Btrans_loop1:do iBtrans=1,1+(fit_by_Btrans&
                                    &*number_of_Btrans_coeffs - fit_by_Btrans)
                                  imf_angle_loop1:do iimf_angle=1,1+(fit_by_imf_angle&
                                       &*number_of_imf_angle_coeffs - fit_by_imf_angle) 
                                     sinT_loop1:do isint=1,1+(fit_by_sint&
                                          &*number_of_sint_coeffs - fit_by_sint) 
                                        AO(M,K,iBtrans,iimf_angle,isint) = AO(M,K,iBtrans,iimf_angle,isint) +  &                               
                                             &    DQS(M,K)*b_Ed2(i,fit_count)*RDTST  
                                        AO(-M,K,iBtrans,iimf_angle,isint) = AO(-M,K,iBtrans,iimf_angle,isint)   &                               
                                             &     -b_Ed1(i,fit_count)*QS(-M,K)*(-M)*RDT

                                        if(present(N_AO) .and. present(FAC_AO)) then

                                           FAC_AO(M,K,iBtrans,iimf_angle,isint) = FAC_AO(M,K,iBtrans,iimf_angle,isint) +  &                               
                                                &      n_ao(2*IX+NS(MM))*(n_ao(2*IX+NS(MM))+1)*DQS(M,K)*b_Ed2(i,fit_count)*RDTST  
                                           FAC_AO(-M,K,iBtrans,iimf_angle,isint) = FAC_AO(-M,K,iBtrans,iimf_angle,isint)   &
                                                &       -n_ao(2*IX+NS(MM))*(n_ao(2*IX+NS(MM))+1)*b_Ed1(i,fit_count)*QS(-M,K)*(-M)*RDT 
                                        end if

                                        fit_count = fit_count + 1 

                                     end do sinT_loop1
                                  end do imf_angle_loop1
                               end do Btrans_loop1
360                            CONTINUE                                              
350                            CONTINUE   

370                            CONTINUE   

380                            FORMAT(7E11.3)                     

                               return
                             end subroutine calc_ao_coeff                                              
                           end module calc_ao_mod

                           !C---------------------------------------------------------------------------                                                                                  
         SUBROUTINE FCMP_single(MMX,CP,SP,F)

           implicit none 

           integer, intent(in) :: mmx
           real, intent(in) ::CP,SP 
           real, intent(out) :: F(-MMX:MMX)			

           integer :: m
           real, parameter :: SQ2=1.4142135623

           F(0) = 1.								 
           F(-1) = SQ2*CP						   
           F(1) = SQ2*SP						   
           do M=2,MMX							      
              F(-M) = CP*F(1-M) - SP*F(M-1)			  
              F(M) = CP*F(M-1) + SP*F(1-M) 
           enddo

           RETURN								     
         END  SUBROUTINE FCMP_single
         !C---------------------------------------------------------------------------  									    
           SUBROUTINE FCMP(MMX,NLONT,CP,SP,F) 

             implicit none             

             integer, intent(in) :: MMX,NLONT
             real, intent(in) ::CP(NLONT),SP(NLONT)
             real, intent(out)::F(-MMX:MMX,NLONT)

             integer :: m,lon                          
             real, parameter :: SQ2=1.4142135623                                              

             DO 60 LON=1,NLONT                                                         
                F(0,LON) = 1.                                                             
                F(-1,LON) = SQ2*CP(LON)                                                   
60              F(1,LON) = SQ2*SP(LON)                                                    
                DO 85 M=2,MMX                                                             
                   DO 85 LON=1,NLONT                                                         
                      F(-M,LON) = CP(LON)*F(1-M,LON) - SP(LON)*F(M-1,LON)                       
85                    F(M,LON) = CP(LON)*F(M-1,LON) + SP(LON)*F(1-M,LON)                        
                      RETURN                                                                    
                    END SUBROUTINE FCMP  

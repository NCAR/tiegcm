module qgen_mod

  use common_model_mod, only: MMX,MXNMX,ITHMX,NQ,pi

  implicit none

  integer, parameter :: ITHTRNS=18,NQS=28,&
       &   MROW=60,KMX=120,MCOEF=49
  integer, parameter :: inlats=23
  integer :: NS(0:MMX),NMX(0:MMX),NSS(0:MMX),IBM(-MMX:MMX),IEM(-MMX:MMX) 

  real :: ST(0:ITHMX),CT(0:ITHMX), &          
       &   VFMPS(NQS,MXNMX),DTH ,ZFMPS(NQS,MXNMX), &                                                       
       &   WNT2(NQ),WNP2(NQ),QNORM(NQ)

  real, parameter :: dlatin=2.,&
       &   Re= 6371200., Ri = 6481200.

contains
!---------------------------------------------------------------------------------
  subroutine gen_coeffs(Q,DQ,N_AO)
    ! ITHMX = number of uniform intervals between pole and equator.
    !   ITHMX = 45 gives 2 degree latitude grid.
    ! ITHTRNS = index for colatitude where basis functions transition from
    !   wavy (high latitudes) to decaying (middle and low latitudes).
    !   This colatitude is (90 degrees)*ITHTRNS/ITHMX.
    !   ITHTRNS = 18 corresponds to 36 colatitude (54 latitude) when ITHMX=45.
    ! MMX = maximum longitudinal wavenumber of Fourier harmonics
    ! MXNMX = number of symmetric Q functions for a given longitudinal wavenumber.
    ! NQS = (MMX+1)*MXNMX = total number of symmetric Q functions.
    ! NQ = 2*NQS = total number of Q functions
    !
    !   NMX(M=0-MMX) = number of Q's for each longitudinal wavenumber M.
    !      (In this version, all NMX's are set to MXNMX.)
    !   MXNMX = max(NMX) (must be > MMX; program fails in Cray double precision
    !	    for MXNMX > 17)
    !   NQS = SUM(NMX)

    real, intent(out) :: Q(NQ,0:ITHMX),DQ(NQ,0:ITHMX)
    real, intent(out), optional :: N_AO(NQ) ! Laplacian used for FAC calculation

    integer :: m,n,i,j,k,kk,ith,nstart,JMAX,js,ja,iter,nx,na

    real :: AOB,AOBN,COA,COATNPO,COAS,UM ,F,TH,TNPO,XPEAK,XNORM,XA,XC,&
         & XB,XN,XNO,QPADQ,DQPADQ,XKPS,XKTS,x,alpha,pnm,dpnmdt
    real :: P(-1:KMX,0:ITHMX),R(0:KMX),DP(KMX,0:ITHMX) &                     
         &  ,FAK(KMX),B(MROW),XJ(MROW) ,FAKZ(KMX)                                                              
    real :: SHT(0:ITHMX),CHT(0:ITHMX),TANHTM(0:ITHMX)                       
    real :: ARRAYS(MXNMX,MXNMX),ARRAYA(MXNMX,MXNMX),AM(MROW,MXNMX)
    logical, parameter :: debug = .false.

    NS(0) = 1                                                                 
    NSS(0) = 1                                                                
    DO 10 M=0,MMX                                                             
       !c       NMX(M) = MXNMX - M 
       NMX(M) = MXNMX 
       IF (M.EQ.MMX) GO TO 10                                                   
       NSS(M+1) = NSS(M) + NMX(M)                                               
       NS(M+1) = NS(M) + 2*NMX(M)                                              
       !C  Zero out entire V,ZFMPS arrays
       DO 5 K=1,NQS	
          ZFMPS(K,M+1) = 0.	
          VFMPS(K,M+1) = 0.   
5      CONTINUE
10  CONTINUE								      
    AOB = Re/Ri 						      
    AOBN = 1.								      
    COA = 6121./6371.							      
    COATNPO = 1./COA							      
    COAS = COA*COA							      
    UM = pi*4.E-7					       
    DO 20 N=1,KMX							      
       AOBN = AOBN*AOB  							
       COATNPO =COATNPO*COAS							
       FAKZ(N) = UM*AOBN*(1. - COATNPO)/FLOAT(2*N+1)/6.371E-3			
20  FAK(N) = UM*AOBN*(1./FLOAT(N) +COATNPO/FLOAT(N+1))/FLOAT(2*N+1)	     
    !C 6.371E-3 IS RADIUS OF EARTH * 1.E-9 TO INCLUDE CONVERSION TO NANOTESLA	     
    DTH = pi/(2.*ITHMX) 					 
    F = FLOAT(ITHMX - ITHTRNS)/FLOAT(ITHMX)				      
    F = F*F								      
    DO 100 ITH=0,ITHMX  						      
       TH = ITH*DTH								
       P(-1,ITH) = 0.								
       P(0,ITH) = 1.								
       CHT(ITH) = COS(.5*ITH*DTH)						
       SHT(ITH) = SIN(.5*ITH*DTH)						
       ST(ITH) = SIN(TH)						       
100 CT(ITH) = COS(TH)							    
    ST(0) = 1.E-20							      
    DO 1000 M=0,MMX							      
       R(M) = 0.								
       DO 110 N=M+1,KMX 							
110  	  R(N) = SQRT((N*N - M*M)/(4.*N*N - 1)) 				  
       DO 160 N=M+1,KMX 							
          TNPO = 2*N + 1							  
          DO 150 ITH=1,ITHMX							  
     	     P(N,ITH) = (CT(ITH)*P(N-1,ITH) - R(N-1)*P(N-2,ITH))/R(N)		    
150  	  DP(N,ITH) = (N*CT(ITH)*P(N,ITH)-TNPO*R(N)*P(N-1,ITH))/ST(ITH) 	 
     	  P(N,0) = (P(N-1,0) - R(N-1)*P(N-2,0))/R(N)				  
     	  DP(N,0) = 0.  							  
     	  IF (M.EQ.1) DP(N,0) = P(N,1)/DTH					  
160    CONTINUE 								
       J = 0									
       NSTART = M								
       IF (M.EQ.0) NSTART = 2							
       DO 170 N=NSTART,KMX,2							
     	     	J = J + 1								
170    XJ(J) = N*(N+1)  						       
       IF (M.EQ.0) THEN 							
          DO 240 ITH=ITHTRNS,ITHMX						  
240          TANHTM(ITH) = 1.							    
       ELSE  								   
             DO 250 ITH=ITHTRNS,ITHMX						     
250     	TANHTM(ITH) = TANHTM(ITH)*SHT(ITH)/CHT(ITH)			       
       ENDIF								      
       DO 300 N=1,NMX(M)						       
          J = NS(M) + 2*(N-1)							 
          DO 270 ITH=ITHTRNS,ITHMX						 
             X = 1./TANHTM(ITH) 						   
             Q(J,ITH)	 = X + TANHTM(ITH)					      
             Q(J+1,ITH)  = X - TANHTM(ITH)					    
             DQ(J,ITH)   = -M*Q(J+1,ITH)/ST(ITH)				     
             DQ(J+1,ITH) = -M*Q(J,ITH)/ST(ITH) 
270       CONTINUE								    
300    CONTINUE 							       
       XPEAK = .75								
       JMAX = NS(M) + 2*NMX(M) -1					       
       if(debug) WRITE (6,"(1X,4X,'M',4X,'J ITER',8X,'XN',7X,'XNO SQRT WNT2 SQRT' &
            & ' WNP2',9X,'X	QNORM')")
       DO 400 J=NS(M),JMAX							
          XPEAK = XPEAK + .5							  
          N = XPEAK								  
          XNORM = Q(J,ITHTRNS)  						  
          IF (M.EQ.0) XNORM = 1.						   
          ALPHA = -DQ(J,ITHTRNS)/XNORM  					   
          IF (M.NE.0) XNORM = ALPHA						  
          XA = 1. - F								  
          XC = M + 2*N - 2*MIN0(M,1)						  
          XB = XC + F/2.							  
          XC = XC*XC + F*M*M							  
          XN = (XB + SQRT(XB*XB - XA*XC))/XA					  
          XNO = XN								  
          XA = .5							       
          DO 310 ITER=1,32							  
             CALL GENPNM(M,XN,CT(ITHTRNS),PNM,DPNMDT)	    				    
             QPADQ = ALPHA*PNM + DPNMDT 			
             IF (ABS(QPADQ).LT..5E-4*ABS(XNORM*PNM)) GO TO 320     
             CALL GENPNM(M,XN-.1,CT(ITHTRNS),PNM,DPNMDT)			    
             DQPADQ = ALPHA*PNM + DPNMDT					    
             CALL GENPNM(M,XN+.1,CT(ITHTRNS),PNM,DPNMDT)  

             DQPADQ = (ALPHA*PNM + DPNMDT - DQPADQ)/.2  			
             X = -QPADQ/DQPADQ  						    
             X =  AMAX1(-XA,AMIN1(XA,X))					    
             XN = XN + X							    
             XA = .75*XA							    
310       CONTINUE								  
320       CONTINUE

          ! Set Laplacian, if present
          if (present(N_AO)) then
             N_AO(J) = XN
          end if

          DO 330 ITH=0,ITHTRNS-1						  
             CALL GENPNM(M,XN,CT(ITH),Q(J,ITH),DQ(J,ITH))   
             IF (M.EQ.0.AND.Q(J,ITHTRNS).LT..5) THEN				    
             	Q(J,ITH) = Q(J,ITH) - PNM			    
             ELSE								    
             	Q(J,ITH) = Q(J,ITH)*Q(J,ITHTRNS)/PNM				      
             	DQ(J,ITH) = DQ(J,ITH)*Q(J,ITHTRNS)/PNM 
             ENDIF
330       CONTINUE								  
          I = J
367       X = .5*Q(I,ITHMX)*Q(I,ITHMX)  					  
          XKPS = M*M*X  							  
          XKTS = .5*DQ(I,ITHMX)*DQ(I,ITHMX)					  
          DO 380 ITH=1,ITHMX-1  						  
             X = X + Q(I,ITH)*Q(I,ITH)*ST(ITH)  					      
             XKTS = XKTS + DQ(I,ITH)*DQ(I,ITH)*ST(ITH)  			
380       XKPS = XKPS + M*M*Q(I,ITH)*Q(I,ITH)/ST(ITH)				 
          XNORM = SQRT(1./(DTH*(XKTS + XKPS)))  				  
          IF (PNM.LT.0..AND.M.NE.0) XNORM = - XNORM				  
          DO 390 ITH=0,ITHMX							  
             Q(I,ITH) = XNORM*Q(I,ITH)  					    
             DQ(I,ITH) = XNORM*DQ(I,ITH)
390       continue					  
          X = .5*Q(I,ITHTRNS)*Q(I,ITHTRNS)*ST(ITHTRNS)  			  

          XKPS = M*M*X/(ST(ITHTRNS)*ST(ITHTRNS))				  

          XKTS = .5*DQ(I,ITHTRNS)*DQ(I,ITHTRNS)*ST(ITHTRNS)			  
          DO 391 ITH=1,ITHTRNS-1						  
             X = X + Q(I,ITH)*Q(I,ITH)*ST(ITH)  				    
             XKTS = XKTS + DQ(I,ITH)*DQ(I,ITH)*ST(ITH)  			    
391       XKPS = XKPS + M*M*Q(I,ITH)*Q(I,ITH)/ST(ITH)				 
          WNT2(I) = XKTS/X							  
          WNP2(I) = XKPS/X							  
          QNORM(I) = SQRT(XN*(XN+1)*(1.-CT(ITHTRNS)))				  
          X = X*DTH*XN*(XN+1)							  
          !C***********************************************************************	   
          if(debug) WRITE(6,392)M,J,ITER,XN,XNO,SQRT(WNT2(I)),SQRT(WNP2(I)),X,QNORM(I)        
392       FORMAT(1X,3I5,6F10.4) 						    
          !C***********************************************************************	   
400    CONTINUE 								
       DO 4019 K=NSS(M),NSS(M)+NMX(M)-1 				       
          DO 4018 I=1,NMX(M)							 
             ZFMPS(K,I) = 0.							   
4018      VFMPS(K,I) = 0.						   
4019   CONTINUE 							       
       DO 405 N=1,NMX(M)							
          NX = NS(M) + 2*(N-1)  						  
          NA = NX + 1								  
          DO 401 K=1,MROW							  
401          AM(K,N) = 0.							    
             !C***********************************************************************        
          DO 402 J=1,NMX(M)							  
             ARRAYS(J,N) = 0.							    
402        ARRAYA(J,N) = 0.							  
           !C***********************************************************************	    
           DO 404 ITH=1,ITHMX							   
              X = DTH*ST(ITH)							     
              IF (ITH.EQ.ITHMX) X = .5*X					     
              !C***********************************************************************        
              !C CALCULATE CROSS-CORRELATIONS OF GRAD(Q)'S				       
              DO 403 J=1,NMX(M) 						    
        	 JS = NS(M) + 2*(J-1)						       
        	 JA = JS + 1							       
        	 ARRAYS(J,N) = ARRAYS(J,N) + X*(DQ(JS,ITH)*DQ(NX,ITH)  &		
        	      &      + M*M*Q(JS,ITH)*Q(NX,ITH)/(ST(ITH)*ST(ITH)))			 
403           ARRAYA(J,N) = ARRAYA(J,N) + X*(DQ(JA,ITH)*DQ(NA,ITH) &		     
        	   &	  + M*M*Q(JA,ITH)*Q(NA,ITH)/(ST(ITH)*ST(ITH)))  		      
              !C***********************************************************************        
              J = 0								     
              DO 4035 K=NSTART,KMX,2						     
        	 J = J + 1							       
        	 AM(J,N) = AM(J,N) + X*(DQ(NX,ITH)*DP(K,ITH) &  			
        	      &      + M*M*Q(NX,ITH)*P(K,ITH)/(ST(ITH)*ST(ITH)))			 
4035          CONTINUE  							     
404        CONTINUE								   
405    CONTINUE 								
       !C***********************************************************************	
       if (debug) WRITE(6,406) M							   
406    FORMAT(1X,'SYMMETRIC INTEGRALS FOR M=',I3)				
       DO 407 J=1,NMX(M)							
          if (debug) WRITE(6,408) (ARRAYS(J,N),N=1,NMX(M))				     
407    CONTINUE 								
408    FORMAT(1X,13F10.7)							
       if (debug) WRITE(6,409) M							   
409    FORMAT(1X,'ANTISYMMETRIC INTEGRALS FOR M=',I3)				
       DO 410 J=1,NMX(M)							
          if (debug) WRITE(6,408) (ARRAYA(J,N),N=1,NMX(M))				     
410    CONTINUE 								
       if (debug) WRITE(6,416) M							   
416    FORMAT(1X,'LEGENDRE COEFFICIENTS (SYMMETRIC) FOR M=',I3) 		
       !C***********************************************************************	
       J = 0									
       DO 420 N=NSTART,KMX,2							
          J = J + 1								  
          DO 417 I=1,NMX(M)							  
417       B(I) = AM(J,I)/XJ(J)  						 
          !C***********************************************************************	   
          if(debug) WRITE(6,408) (B(I),I=1,NMX(M))					    
          !C***********************************************************************	   
          !C CALCULATE VFMPS (CONVERTS EQIV. CUR. TO MAG. POT.) 			   
          DO 419 K=NSS(M),NSS(M)+NMX(M)-1					  
             KK = K - NSS(M) + 1						    
             DO 418 I=1,NMX(M)  						    
        	ZFMPS(K,I) = ZFMPS(K,I) + AM(J,I)*FAKZ(N)*AM(J,KK)		      
418          VFMPS(K,I) = VFMPS(K,I) + AM(J,I)*FAK(N)*AM(J,KK)  		   
419       CONTINUE								  
420    CONTINUE 								
       IF (M.EQ.MMX) GO TO 1000 					       
       X = SQRT(1. + .5/FLOAT(M+1))						
       DO 999 ITH=1,ITHMX							
          P(M+1,ITH) = X*P(M,ITH)*ST(ITH)					  
999    DP(M+1,ITH) = (M+1)*CT(ITH)*P(M+1,ITH)/ST(ITH)			       
       P(M+1,0) = 0.								
       DP(M+1,0) = 0.								
       IF (M.EQ.0) DP(M+1,0) = P(M+1,1)/DTH					
1000   CONTINUE 								 
       RETURN									 

       end subroutine gen_coeffs 
       end module qgen_mod    							       
!C-----------------------------------------------------------------------------                                                                                
      SUBROUTINE GENPNM(MS,N,CT,PNM,DPNMDT)                                     
        !C MS = ORDER                                                                    
        !C N = DEGREE                                                                    
        !C CT = COS(THETA)                                                               
        !C PNM = GENERALIZED ASSOCIATED LEGENDRE FUNC. (ARBITRARILY NORMALIZED)          
        !C DPNMDT = DERIVATIVE OF PNM WITH RESPECT TO THETA                              
        !c      PARAMETER (PRECISE=1.E-16,BIG=1.E60)
        implicit none 

        integer, intent(in):: MS
        real, intent(in) :: ct,N
        real, intent(out) :: PNM,DPNMDT 
        !
        integer :: m,MPREV,i,j,k                                        
        real, PARAMETER :: PRECISE=1.E-7,BIG=1.E20                                      
        real :: CTPREV,OMCTOT,STMM1,RATMAX,&
             &  X,KZERO,KMAX,KPEAK,AO,ak ,ake ,st,stm                                                             
        real ::  PNMD,DPNMDTD,PNME,DPNMDTE 
        SAVE MPREV,CTPREV,STM,ST,OMCTOT                                           
        DATA MPREV,CTPREV/-999999,-99999./  
        logical, parameter :: debug = .false.

        
	M = IABS(MS)                                                              
        IF (N.LT.M) THEN                                                          
           if(debug) WRITE(6,20) M,N                                                          
20         FORMAT (1X,'STOPPED IN GENPNM BECAUSE M =',I6,'N =',F10.2)               
           STOP                                                                     
        ENDIF

!        IF (M.EQ.MPREV.AND.CT.EQ.CTPREV) GO TO 100

        MPREV = M                                                                 
        CTPREV = CT                                                               
        ST = SQRT(1. - CT*CT)                                                     
        OMCTOT = (1. - CT)/2.                                                     
        IF (OMCTOT.GT..5) THEN                                                    
           if(debug) WRITE(6,50) CT                                                           
50         FORMAT (1X,'STOPPED IN GENPNM BECAUSE CT =',F10.5)                       
           STOP                                                                     
        ENDIF
        STM = 1.                                                                  
        STMM1 = 1.                                                             
        IF (M.GT.0) THEN                                                          
           DO 60 I=1,M                                                              
              IF (I.GT.1) STMM1 = STMM1*ST                                            
60         STM = STM*ST                                                            
        ENDIF									 
100     CONTINUE								 
        RATMAX = 1.								  
        IF (M.GT.1) THEN							  
           X = M*(2. + SQRT(2.*(M+1)))/FLOAT(M-1)				     
           RATMAX = (2*M+X)*(X-1.)/(X*(X+M))					     
        ENDIF
        RATMAX = RATMAX*OMCTOT  						  
        X = 0.  								  
        IF (RATMAX.GT.PRECISE)  &						   
             & X = ALOG((1.-RATMAX)*PRECISE)/ALOG(RATMAX)				 
        KZERO = N + 1 - M + .5  						  
        KMAX = N + 1 - M + .5 + X						  
        KMAX = MAX0(int(KMAX),int(KZERO)+1)						    
        X = .5*(OMCTOT + M)/(1.+OMCTOT) 					  
        !C***********************************************************************	 
        !C ALTHOUGH FOLLOWING FORMULA PICKS OUT MAXIMUM AMPLITUDE TERM IN SERIES	 
        !C  FOR PURPOSE OF SUMMING THE SERIES FROM BOTH ENDS, SMALLEST TO		 
        !C  LARGEST, THINGS SEEM TO WORK BETTER IF KPEAK = KZERO.			 
        !C     KPEAK = X + SQRT(X*X + N*(N+1.)*OMCTOT/(OMCTOT+1.)) - M  		 
        !C     KPEAK = MIN0(KPEAK,KZERO)						 
        !C     KPEAK = MAX0(KPEAK,1)							 
        !C***********************************************************************	 
        KPEAK = KZERO								  
        AO = 1. 								  
        AK = AO*(M*(M+1) - N*(N+1.))/FLOAT(M+1) 				  
        PNMD = AK								  
        DPNMDTD = AK								  
        IF (KPEAK.LT.2) GO TO 201						  
        DO 200 K=2,KPEAK							  
           AK = AK*OMCTOT*((K+M-1)*(K+M) - N*(N+1))/FLOAT(K*(K+M))		    
           PNMD = PNMD + AK							    
           DPNMDTD = DPNMDTD + K*AK						    
           IF (AK.GT.BIG) THEN  						    
              AO = AO/AK							      
              PNMD = PNMD/AK							      
              DPNMDTD = DPNMDTD/AK						      
              AK = 1.								      
           ENDIF
200     CONTINUE								 
201     CONTINUE								  
        AKE = 1.								  
        PNME = AKE								  
        DPNMDTE = (KMAX+1)*AKE 					  
        IF (OMCTOT.LT.PRECISE) GO TO 400					  
        DO 300 K=KMAX,KZERO+1,-1						  
           J = K + M + 1							    
           AKE = AKE*J*(J-M)/(OMCTOT*((J-1)*J - N*(N+1)))			    
           PNME = PNME + AKE							     
           DPNMDTE = DPNMDTE + K*AKE						    
           IF (AKE.GT.BIG) THEN 						    
              PNME = PNME/AKE							      
              DPNMDTE = DPNMDTE/AKE						      
              AKE = 1.  							      
           ENDIF
300     CONTINUE								  
        IF (KZERO.EQ.KPEAK) GO TO 350						  
        J = KZERO + 1 + M							  
        X = OMCTOT*((J-1)*J - N*(N-1))/(AKE*(J-M)*J)				  
        PNME = X*PNME + 1.							  
        DPNMDTE = X*DPNMDTE + KZERO						  
        AKE = 1.								  
        IF (KZERO.EQ.KPEAK+1) GO TO 350 					  
        DO 310 K=KZERO-1,KPEAK+1,-1						  
           J = K + M + 1							    
           AKE = AKE*J*(J-M)/(OMCTOT*((J-1)*J - N*(N+1)))			    
           PNME = PNME + AKE							    
           DPNMDTE = DPNMDTE + K*AKE						    
           IF (AKE.GT.BIG) THEN 						    
              PNME = PNME/AKE							      
              DPNMDTE = DPNMDTE/AKE						      
              AKE = 1.  							      
           ENDIF
310     CONTINUE								  
350     K = KPEAK + 1								  
        X = AK*OMCTOT*((K+M-1)*(K+M) - N*(N+1))/(AKE*(K*(K+M))) 		  
        PNMD = PNMD + X*PNME							  
        DPNMDTD = DPNMDTD + X*DPNMDTE						  
400     CONTINUE								  
        PNM = PNMD*OMCTOT + AO  						  
        DPNMDT = DPNMDTD*.5*ST*STM + M*CT*STMM1*PNM				  
        PNM = PNM*STM 

        RETURN   

        end subroutine genpnm	
!---------------------------------------------------------------------------------------------

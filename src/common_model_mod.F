       module common_model_mod
       
         implicit none
       
         integer, parameter :: MMX=3,MXNMX=7,ITHMX=45,NQ=56,ITHTRNS=18
       
         ! Number of mlat intervals:
         integer :: nmlat = 23
       
         integer,parameter :: mlat_resolution = 2
       
         real, dimension( 23) :: average_mlat (1:23) = (/ 46,48,
     |     50,52,54,56,58,
     |     60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90  /) 
       
         real, parameter :: pi = 3.1415926535898
         real, parameter :: magB0_apex= 54312.23 ! magnitude of magnetic field nT
         real, parameter :: hr = 110.   ! reference height 110 km
       
         ! for correction of coefficients
         real, parameter :: eps = 0.05
       
         ! Data is always fitted by MLT:
         integer, parameter :: fit_by_MLT = 1
       
       contains
       
       end module common_model_mod                                                                                                                           

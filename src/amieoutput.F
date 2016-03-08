#include "defs.h"
      subroutine amieoutput
      use params_module,only: nlonp1,nlonp4,nlevp1,nlatp1,nlat,nlev,nlon
      use fields_module,only: ped,hall
      use dyndiag_module,only: 
     |  qwind_sec,qamie_sec,work_sec,wtot_sec,fwindu_sec,
     |  fwindv_sec,famieu_sec,famiev_sec,famie_sec,fwind_sec,
     |  ftot_sec,tec_sec,ped_sec,hall_sec,qjoule_sec
      use cons_module,only: pi, dlamda, dphi, sn
      use magfield_module,only: dipmag
      use addfld_module,only: addfld

      implicit none

      real, dimension(nlonp4,nlat) ::
     |  sigp,    ! height integrated pedersen conductivity
     |  sigh,    ! height integrated hall conductivity
     |  qwind,   ! joule heating associated with neutral winds
     |  qamie,   ! joule heating from amie
     |  wtot,    ! total electromagnetic power
     |  work,    ! mechanical work
     |  fwindu,fwindv, ! horizontal current
     |  famieu,famiev, ! horizontal current from amie
     |  tec      ! total electron content
!
! Local:
      integer :: j1,i,j,k,lonm,lonp
      real :: re
      real, dimension(nlonp4,nlat,nlevp1) :: qjh_tot_sec
!
!CCCCCCCCCCCCCCC calculate the field-aligned currents CCCCCCCCCCCCCCCCCC
        re = 6.481e+6

       fwind_sec(:,:,:) = 0.
       famie_sec(:,:,:) = 0.
!      write(6,"('DIPMAG(23,:) = ',/,(6f12.2))")dipmag(23,:) 
       do j = 2, nlat-1
         j1 = nlat+1-j
          do i=3, nlon+2
!          if (abs(sin(dipmag(i,j))) .gt. 0.25) then
           if (abs(dipmag(i,j)) .gt. 0.55) then
             lonm = i -1
             if (lonm .eq. 2) lonm = nlon+2
             lonp = i + 1
             if (lonp .eq. nlon+3) lonp = 3
C *** calculate height-profiles of field-aligned currents ***
              do k=1,nlevp1
                fwind_sec(i,j,k) =
     |          ((fwindu_sec(lonp,j,k)-fwindu_sec(lonm,j,k))/
     |          dlamda + (sn(j+1)*fwindv_sec(i,j+1,k)-sn(j-1)*
     |          fwindv_sec(i,j-1,k))/dphi) / (2.*re*sn(j)) /
     |          (-sin(dipmag(i,j)))
                famie_sec(i,j,k) =
     |          ((famieu_sec(lonp,j,k)-famieu_sec(lonm,j,k))/
     |          dlamda + (sn(j+1)*famiev_sec(i,j+1,k)-sn(j-1)*
     |          famiev_sec(i,j-1,k))/dphi) / (2.*re*sn(j)) /
     |          (-sin(dipmag(i,j)))
!  Total field-aligned currents
                ftot_sec(i,j,k) = famie_sec(i,j,k)+fwind_sec(i,j,k)
              enddo
C *** end of calculate height-profiles of field-aligned currents ***
            endif

          enddo

       enddo

          do k=1,nlevp1
           fwind_sec(:,1,k)=(fwind_sec(:,2,k)-fwind_sec(:,1,k))/re/dphi
           fwind_sec(:,nlat,k) = (fwind_sec(:,nlat,k)-
     +           fwind_sec(:,nlat-1,k))/re/dphi
           famie_sec(:,1,k)=(famie_sec(:,2,k)-famie_sec(:,1,k))/re/dphi
           famie_sec(:,nlat,k) = (famie_sec(:,nlat,k)-
     +           famie_sec(:,nlat-1,k))/re/dphi
           ftot_sec(:,1,k) = fwind_sec(:,1,k)+famie_sec(:,1,k)
           ftot_sec(:,nlat,k)=fwind_sec(:,nlat,k)+famie_sec(:,nlat,k)
          enddo
CCCCCCCCCCCCCCCCCCC end of FAC calculation CCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCC save secondary histories CCCCCCCCCCCCCCCCCCCCC
!
! Add sigmas to secondary histories:
! ZIMXP = 76, ZIMXP1 = 72+1 = 73, ZKMX = 28, ZKMXP = ZKMX+1 = 29
! From dynamo.h: SIGMA1(ZIMXP1,ZJMX,ZKMX),SIGMA2(ZIMXP1,ZJMX,ZKMX)
!
      do j = 1, nlat
        ped_sec(1,j,:) = ped_sec(nlon+1,j,:) ! 73 <- 1
        ped_sec(2,j,:) = ped_sec(nlon+2,j,:) ! 74 <- 2
        ped_sec(nlonp1+2,j,:) = ped_sec(3,j,:) ! 75 <- 3
        ped_sec(nlonp1+3,j,:) = ped_sec(4,j,:) ! 76 <- 4
        call addfld('PEDERSEN',' ',' ',ped_sec(:,j,:),
     |    'lon',1,nlonp4,'lev',1,nlevp1,j)
        hall_sec(1,j,:) = hall_sec(nlon+1,j,:) ! 73 <- 1
        hall_sec(2,j,:) = hall_sec(nlon+2,j,:) ! 74 <- 2
        hall_sec(nlonp1+2,j,:) = hall_sec(3,j,:) ! 75 <- 3
        hall_sec(nlonp1+3,j,:) = hall_sec(4,j,:) ! 76 <- 4
        call addfld('HALL',' ',' ',hall_sec(1:nlonp4,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)
        qwind_sec(1,j,:) = qwind_sec(nlon+1,j,:) ! 73 <- 1
        qwind_sec(2,j,:) = qwind_sec(nlon+2,j,:) ! 74 <- 2
        qwind_sec(nlonp1+2,j,:) = qwind_sec(3,j,:) ! 75 <- 3
        qwind_sec(nlonp1+3,j,:) = qwind_sec(4,j,:) ! 76 <- 4
        call addfld('QWIND',' ',' ',qwind_sec(:,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)
        qamie_sec(1,j,:) = qamie_sec(nlon+1,j,:) ! 73 <- 1
        qamie_sec(2,j,:) = qamie_sec(nlon+2,j,:) ! 74 <- 2
        qamie_sec(nlonp1+2,j,:) = qamie_sec(3,j,:) ! 75 <- 3
        qamie_sec(nlonp1+3,j,:) = qamie_sec(4,j,:) ! 76 <- 4
        call addfld('QAMIE',' ',' ',qamie_sec(1:nlonp4,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)
        fwind_sec(1,j,:) = fwind_sec(nlon+1,j,:) ! 73 <- 1
        fwind_sec(2,j,:) = fwind_sec(nlon+2,j,:) ! 74 <- 2
        fwind_sec(nlonp1+2,j,:) = fwind_sec(3,j,:) ! 75 <- 3
        fwind_sec(nlonp1+3,j,:) = fwind_sec(4,j,:) ! 76 <- 4
        call addfld('FWIND',' ',' ',fwind_sec(1:nlonp4,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)
        famie_sec(1,j,:) = famie_sec(nlon+1,j,:) ! 73 <- 1
        famie_sec(2,j,:) = famie_sec(nlon+2,j,:) ! 74 <- 2
        famie_sec(nlonp1+2,j,:) = famie_sec(3,j,:) ! 75 <- 3
        famie_sec(nlonp1+3,j,:) = famie_sec(4,j,:) ! 76 <- 4
        call addfld('FAMIE',' ',' ',famie_sec(1:nlonp4,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)
        work_sec(1,j,:) = work_sec(nlon+1,j,:) ! 73 <- 1
        work_sec(2,j,:) = work_sec(nlon+2,j,:) ! 74 <- 2
        work_sec(nlonp1+2,j,:) = work_sec(3,j,:) ! 75 <- 3
        work_sec(nlonp1+3,j,:) = work_sec(4,j,:) ! 76 <- 4
        call addfld('WORK',' ',' ',work_sec(1:nlonp4,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)
        wtot_sec(1,j,:) = wtot_sec(nlon+1,j,:) ! 73 <- 1
        wtot_sec(2,j,:) = wtot_sec(nlon+2,j,:) ! 74 <- 2
        wtot_sec(nlonp1+2,j,:) = wtot_sec(3,j,:) ! 75 <- 3
        wtot_sec(nlonp1+3,j,:) = wtot_sec(4,j,:) ! 76 <- 4
        call addfld('WTOT',' ',' ',wtot_sec(:,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)
        tec_sec(1,j,:) = tec_sec(nlon+1,j,:) ! 73 <- 1
        tec_sec(2,j,:) = tec_sec(nlon+2,j,:) ! 74 <- 2
        tec_sec(nlonp1+2,j,:) = tec_sec(3,j,:) ! 75 <- 3
        tec_sec(nlonp1+3,j,:) = tec_sec(4,j,:) ! 76 <- 4
        call addfld('DEN_TEC',' ',' ',tec_sec(1:nlonp4,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)

        qjoule_sec(1,j,:) = qjoule_sec(nlon+1,j,:) ! 73 <- 1
        qjoule_sec(2,j,:) = qjoule_sec(nlon+2,j,:) ! 74 <- 2
        qjoule_sec(nlonp1+2,j,:) = qjoule_sec(3,j,:) ! 75 <- 3
        qjoule_sec(nlonp1+3,j,:) = qjoule_sec(4,j,:) ! 76 <- 4
        call addfld('QJOULE',' ',' ',qjoule_sec(1:nlonp4,j,:),
     |             'lon',1,nlonp4,'lev',1,nlevp1,j)

        qjh_tot_sec(:,:,:) = qamie_sec(:,:,:) + qwind_sec(:,:,:)
        call addfld('QJH_TOT',' ',' ',
     |     qjh_tot_sec(1:nlonp4,j,:),'lon',1,nlonp4,'lev',1,nlevp1-1,j)
      enddo

!CCCCCCCCCCCCCCCCCC end of save secondary histories CCCCCCCCCCCCCCCCCCCCC
      end subroutine amieoutput
!----------------------------------------------------
      subroutine prep_output(lev0,lev1,lon0,lon1,lat0,lat1)
!
! Prepare geographic-grid fields for input to the dynamo, and gather them
! to the root task. This is executed by all tasks, and is called from
! advance before the dynamo itself (which is executed by master task only).
!
      use cons_module,only: gask,grav
#ifdef MPI
      use dyndiag_module,only: mp_dyndiag_gather
#endif
!
! Args:  
      integer,intent(in) :: lev0,lev1,lon0,lon1,lat0,lat1
!
! Local:
      real :: fmin,fmax
!
#ifdef MPI
!
! Gather dynamo input fields to the root task, defining module data
! above at the global domain on the root task:
!
      call mp_dyndiag_gather
#endif
      end subroutine prep_output


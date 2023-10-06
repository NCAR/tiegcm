      subroutine amieoutput
      use params_module,only: nlevp1,nlat
      use dyndiag_module,only: 
     |  qwind_sec,qamie_sec,work_sec,wtot_sec,fwindu_sec,
     |  fwindv_sec,famieu_sec,famiev_sec,
     |  tec_sec,ped_sec,hall_sec,qjoule_sec
      use cons_module,only: dlamda, dphi, sn
      use magfield_module,only: dipmag
      use addfld_module,only: addfld
      use mpi_module,only: lon0,lon1,lat0,lat1

      implicit none
!
! Local:
      integer :: i,j,k,lonm,lonp
      real :: re
      real, dimension(lon0:lon1,lat0:lat1,nlevp1) ::
     |  qjh_tot_sec,famie_sec,fwind_sec,ftot_sec
!
!CCCCCCCCCCCCCCC calculate the field-aligned currents CCCCCCCCCCCCCCCCCC
      re = 6.481e+6

      fwind_sec = 0.
      famie_sec = 0.
!     write(6,"('DIPMAG(23,:) = ',/,(6f12.2))")dipmag(23,:) 
      do j = lat0,lat1
        if (j/=1 .and. j/=nlat) then
          do i = lon0,lon1
!           if (abs(sin(dipmag(i,j))) > 0.25) then
            if (abs(dipmag(i,j)) > 0.55) then
              lonm = i-1
              lonp = i+1
! *** calculate height-profiles of field-aligned currents ***
              do k = 1,nlevp1
                fwind_sec(i,j,k) =
     |            ((fwindu_sec(lonp,j,k)-fwindu_sec(lonm,j,k))/dlamda+
     |            (sn(j+1)*fwindv_sec(i,j+1,k)-
     |            sn(j-1)*fwindv_sec(i,j-1,k))/dphi)/
     |            (2.*re*sn(j))/(-sin(dipmag(i,j)))
                famie_sec(i,j,k) =
     |            ((famieu_sec(lonp,j,k)-famieu_sec(lonm,j,k))/dlamda+
     |            (sn(j+1)*famiev_sec(i,j+1,k)-
     |            sn(j-1)*famiev_sec(i,j-1,k))/dphi)/
     |            (2.*re*sn(j))/(-sin(dipmag(i,j)))
!  Total field-aligned currents
                ftot_sec(i,j,k) = famie_sec(i,j,k)+fwind_sec(i,j,k)
              enddo
! *** end of calculate height-profiles of field-aligned currents ***
            endif
          enddo
        endif
      enddo

      if (lat0 == 1) then
        do i = lon0,lon1
          do k = 1,nlevp1
            fwind_sec(i,1,k) =
     |        (fwind_sec(i,2,k)-fwind_sec(i,1,k))/re/dphi
            famie_sec(i,1,k) =
     |        (famie_sec(i,2,k)-famie_sec(i,1,k))/re/dphi
            ftot_sec(i,1,k) = fwind_sec(i,1,k)+famie_sec(i,1,k)
          enddo
        enddo
      endif

      if (lat1 == nlat) then
        do i = lon0,lon1
          do k = 1,nlevp1
            fwind_sec(i,nlat,k) =
     |        (fwind_sec(i,nlat,k)-fwind_sec(i,nlat-1,k))/re/dphi
            famie_sec(i,nlat,k) =
     |        (famie_sec(i,nlat,k)-famie_sec(i,nlat-1,k))/re/dphi
            ftot_sec(i,nlat,k) = fwind_sec(i,nlat,k)+famie_sec(i,nlat,k)
          enddo
        enddo
      endif
!CCCCCCCCCCCCCCCCCC end of FAC calculation CCCCCCCCCCCCCCCCCCCCC
!
!CCCCCCCCCCCCCCCCCC save secondary histories CCCCCCCCCCCCCCCCCCCCC
!
! Add sigmas to secondary histories:
! ZIMXP = 76, ZIMXP1 = 72+1 = 73, ZKMX = 28, ZKMXP = ZKMX+1 = 29
! From dynamo.h: SIGMA1(ZIMXP1,ZJMX,ZKMX),SIGMA2(ZIMXP1,ZJMX,ZKMX)
!
      qjh_tot_sec = qamie_sec(lon0:lon1,lat0:lat1,:)
     |  + qwind_sec(lon0:lon1,lat0:lat1,:)
      do j = lat0,lat1
        call addfld('PEDERSEN',' ',' ',ped_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('HALL',' ',' ',hall_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('QWIND',' ',' ',qwind_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('QAMIE',' ',' ',qamie_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('FWIND',' ',' ',fwind_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('FAMIE',' ',' ',famie_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('WORK',' ',' ',work_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('WTOT',' ',' ',wtot_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('DEN_TEC',' ',' ',tec_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('QJOULE',' ',' ',qjoule_sec(lon0:lon1,j,:),
     |     'lon',lon0,lon1,'lev',1,nlevp1,j)
        call addfld('QJH_TOT',' ',' ',qjh_tot_sec(lon0:lon1,j,:),
     |    'lon',lon0,lon1,'lev',1,nlevp1,j)
      enddo

!CCCCCCCCCCCCCCCCCC end of save secondary histories CCCCCCCCCCCCCCCCCCCCC
      end subroutine amieoutput
!----------------------------------------------------
      subroutine prep_output
!
! Prepare geographic-grid fields for input to the dynamo.
! This is executed by all tasks,
! and is called from advance before the dynamo itself.
!
      use dyndiag_module,only: dyndiag_bndlatlon
!
      call dyndiag_bndlatlon
!
      end subroutine prep_output

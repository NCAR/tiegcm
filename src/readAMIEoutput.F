      subroutine readAMIEoutput(iBLK, IsMirror, iError)

        use ModAMIE_Interface
        use ModFiles

        implicit none

        integer, intent(out) :: iError
        logical, intent(in)  :: IsMirror
        integer, intent(in)  :: iBLK

        integer :: iTime, nTimesBig, nTimesTotal
        integer :: nfields
        integer :: ntemp, iyr, imo, ida, ihr, imi
        integer :: i,j, iField, iPot_, iAveE_, iEFlux_
        real*4  :: swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot
        real*8  :: rtime
        integer, dimension(7) :: itime_i

        real*4, allocatable, dimension(:,:,:) :: AllData
        integer, parameter :: nFieldsMax = 100
        character (len=30), dimension(nFieldsMax) :: Fields

        logical :: IsBinary, energyfluxconvert

        real :: dPotential
        integer :: nCellsPad, n

        iError = 0
        open(LunEField_, file=AMIE_FileName, status='old',
     |    form='UNFORMATTED',iostat=iError)
        if (iError.ne.0) then
           write(*,*) "Error opening file:", AMIE_FileName
           stop
        endif
        AMIE_nLats = 0
        IsBinary = .true.

        read(LunEField_,iostat=iError) AMIE_nlats,AMIE_nmlts,
     |    AMIE_ntimes
        if ((iError.ne.0).or.(AMIE_nlats.gt.100)) then
           write(*,*) "Error reading variables AMIE_nlats, AMIE_nmlts, 
     |       AMIE_ntimes"
           IsBinary = .false.
        endif
        close(LunEField_)

        if (IsBinary) then
           open(LunEField_, file=AMIE_FileName, status='old',
     |        form='UNFORMATTED')
           read(LunEField_) AMIE_nlats,AMIE_nmlts,AMIE_ntimes
        else
           open(LunEField_, file=AMIE_FileName, status='old')
           read(LunEField_,*) AMIE_nlats,AMIE_nmlts,AMIE_ntimes
        endif

  !\
  ! We have run into a problem with AMIE during storms.
  ! The potential is not zero at the boundary.  It is sometimes quite
  ! high.  This means that the gradient in the potential will be
  ! large - meaning that very strong flows can exist in this last
  ! cell.  This is not good.
  ! To rectify this, we will pad the AMIE results by 5 grid cells, and
  ! force the potential to go to zero linearly from the last cell to the
  ! new last cell.
  ! Since it is assumed that all of the AMIE quantities are on the same
  ! grid, we have to fill in the eflux and avee also.  We will use the
  ! last cell to fill in those value.
  ! We also have to extend the grid.
  !/

        nCellsPad = 15

        if (allocated(AMIE_Lats)) deallocate(AMIE_Lats)
        allocate(AMIE_Lats(AMIE_nLats+nCellsPad), stat=iError)
        if (iError /= 0) then
           write(*,*) "Error in allocating array AMIE_Lats in "
           stop
        endif

        if (allocated(AMIE_Mlts)) deallocate(AMIE_Mlts)
        allocate(AMIE_Mlts(AMIE_nMlts), stat=iError)
        if (iError /= 0) then
           write(*,*) "Error in allocating array Mlts in "
           stop
        endif

        if (.not.allocated(AMIE_Potential)) then

           allocate(AMIE_Potential(AMIE_nMlts,AMIE_nLats+nCellsPad, 
     |                             AMIE_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array AMIE_Potential in "
              stop
           endif

           allocate(AMIE_EFlux(AMIE_nMlts,AMIE_nLats+nCellsPad, 
     |                         AMIE_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array AMIE_EFlux in "
              stop
           endif

           AMIE_EFlux = 0.0

           allocate(AMIE_AveE(AMIE_nMlts,AMIE_nLats+nCellsPad, 
     |                        AMIE_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array AMIE_AveE in "
              stop
           endif

           AMIE_AveE = 0.0

           allocate(AMIE_Value(AMIE_nMlts,AMIE_nLats+nCellsPad, 
     |                         AMIE_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array AMIE_Value in "
              stop
           endif

           AMIE_Value = 0.0

           allocate(AMIE_Time(AMIE_nTimes,2), stat=iError)
           if (iError /= 0) then
              write(*,*) "Error in allocating array AMIETimes in "
              stop
           endif

        endif

        if (IsBinary) then
           read(LunEField_) (AMIE_Lats(i),i=1,AMIE_nLats)
           read(LunEField_) (AMIE_Mlts(i),i=1,AMIE_nMlts)
           read(LunEField_) nFields
        else
           read(LunEField_,*) (AMIE_Lats(i),i=1,AMIE_nLats)
           read(LunEField_,*) (AMIE_Mlts(i),i=1,AMIE_nMlts)
           read(LunEField_,*) nFields
        endif

        AMIE_Lats = 90.0 - AMIE_Lats
! write(*,*) "AMIE_Lats=",AMIE_Lats
  !\
  ! Extrapolate the latitude grid
  !/

        do i=AMIE_nLats+1,AMIE_nLats+nCellsPad
           AMIE_Lats(i) = AMIE_Lats(i-1) + 
     |           (AMIE_Lats(AMIE_nLats) - AMIE_Lats(AMIE_nLats-1))
        enddo

        if (nFields > nFieldsMax) then
           write(*,*) "Maximum number of fields in AMIE is ",nFieldsMax
           stop
        endif

        allocate(AllData(AMIE_nMlts,AMIE_nLats,nFields), stat=iError)
        if (iError /= 0) then
           write(*,*) "Error in allocating array AllData in "
           stop
        endif

        AMIE_iDebugLevel = 0

        energyfluxconvert = .true.

        do iField=1,nfields
           if (IsBinary) then
              read(LunEField_) Fields(iField)
           else
              read(LunEField_,'(a)') Fields(iField)
           endif

           if (AMIE_iDebugLevel > 1) write(*,*) Fields(iField)

           if ((index(Fields(iField),"Potential") > 0).and. 
     |        (index(Fields(iField),"odel") < 1)) then
              iPot_ = iField
              if (AMIE_iDebugLevel > 1) 
     |          write(*,*) "<--- Potential Found", iPot_
           endif

           if ((index(Fields(iField),"Mean Energy") > 0) .and. 
     |         (index(Fields(iField),"odel") < 1)) then
              iAveE_ = iField
              if (AMIE_iDebugLevel > 1) 
     |          write(*,*) "<--- Mean Energy Found", iAveE_
           endif

      if ((index(Fields(iField),"Energy Flux") > 0) .and. 
     |         (index(Fields(iField),"odel") < 1)) then
        if (index(Fields(iField),"w/m2") > 0) energyfluxconvert = .true.
        iEFlux_ = iField
        if (AMIE_iDebugLevel > 1) 
     |    write(*,*) "<--- Energy Flux Found", iEFlux_
      endif

        enddo

      do iTime=1,AMIE_ntimes

       if (IsBinary) then

        read(LunEField_) ntemp,iyr,imo,ida,ihr,imi
        read(LunEField_) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,pot

        do iField=1,nfields
           read(LunEField_) ((AllData(j,i,iField),j=1,AMIE_nMlts),
     |       i=1,AMIE_nLats)
        enddo

      else

        read(LunEField_,*) ntemp,iyr,imo,ida,ihr,imi
        read(LunEField_,*) swv,bx,by,bz,aei,ae,au,al,dsti,dst,hpi,sjh,
     |    pot

        do iField=1,nfields
           read(LunEField_,*) ((AllData(j,i,iField),j=1,AMIE_nMlts),
     |       i=1,AMIE_nLats)
        enddo

      endif

      itime_i(1) = iyr
      itime_i(2) = imo
      itime_i(3) = ida
      itime_i(4) = ihr
      itime_i(5) = imi
      itime_i(6) = 0
      itime_i(7) = 0
      call time_int_to_real(itime_i,rtime)
      AMIE_Time(iTime,iBLK) = rtime

! We need Potential to be in Volts
!         AveE to be in keV
!         EFlux to be in W/m2

      if (.not.IsMirror) then
        AMIE_Potential(:,1:AMIE_nLats,iTime,iBLK) = 
     |             AllData(:,1:AMIE_nLats,iPot_)
      else
        ! This is typically done for the Southern Hemisphere, when we don't
        ! have Southern Hemisphere Runs.
        do i=1,AMIE_nMlts
           AMIE_Potential(i,1:AMIE_nLats,iTime,iBLK) = 
     |                -AllData(AMIE_nMlts+1 - i,1:AMIE_nLats,iPot_)
        enddo
      endif

      AMIE_AveE(:,1:AMIE_nLats,iTime,iBLK)      = 
     |    AllData(:,1:AMIE_nLats,iAveE_)

! Need to convert from W/m2 to erg/cm2/s
      if (energyfluxconvert) then
        AMIE_EFlux(:,1:AMIE_nLats,iTime,iBLK)     = 
     |   AllData(:,1:AMIE_nLats,iEFlux_) / (1.0e-7 * 100.0 * 100.0)
      else
           AMIE_EFlux(:,1:AMIE_nLats,iTime,iBLK)     =  
     |                AllData(:,1:AMIE_nLats,iEFlux_)
      endif

      do i=1,AMIE_nMlts

        dPotential = AMIE_Potential(i,AMIE_nLats,iTime,iBLK)/nCellsPad

        do j=AMIE_nLats+1, AMIE_nLats+nCellsPad

           AMIE_AveE(i,j,iTime,iBLK)  = 
     |                minval(AMIE_AveE(:,1:AMIE_nLats,iTime,iBLK))
           AMIE_EFlux(i,j,iTime,iBLK) = 
     |           minval(AMIE_EFlux(:,1:AMIE_nLats,iTime,iBLK))

           AMIE_Potential(i,j,iTime,iBLK) = 
     |           AMIE_Potential(i,j-1,iTime,iBLK) - dPotential

        enddo
      enddo

! One problem with the extension to lower latitudes is that the 
! potential can start having large differences in MLT, giving rise
! to large eastward electric fields.  We will try to compensate
! for this by averaging a bit...

      do j=AMIE_nLats+1, AMIE_nLats+nCellsPad
     
        ! We are going to smooth more and more as we go down in latitude

        do n = 1, j-AMIE_nLats

           i = 1
           AMIE_Potential(i,j,iTime,iBLK) = 
     |           (AMIE_Potential(AMIE_nMlts-1,j,iTime,iBLK) + 
     |           2*AMIE_Potential(i,j,iTime,iBLK) + 
     |           AMIE_Potential(i+1,j,iTime,iBLK))/4.0

           do i=2,AMIE_nMlts-1
              AMIE_Potential(i,j,iTime,iBLK) = 
     |              (AMIE_Potential(i-1,j,iTime,iBLK) + 
     |               2*AMIE_Potential(i,j,iTime,iBLK) + 
     |               AMIE_Potential(i+1,j,iTime,iBLK))/4.0
              if (AMIE_EFlux(i,j,iTime,iBLK) > 1000.0) 
     |              nTimesBig = nTimesBig+1
              nTimesTotal = nTimesTotal + 1
           enddo

           i = AMIE_nMlts
           AMIE_Potential(i,j,iTime,iBLK) = 
     |          (AMIE_Potential(i-1,j,iTime,iBLK) + 
     |          2*AMIE_Potential(i,j,iTime,iBLK) + 
     |          AMIE_Potential(2,j,iTime,iBLK))/4.0

        enddo
      enddo

      enddo

      if (nTimesBig > nTimesTotal*0.1) 
     |       AMIE_EFlux(:,:,:,iBLK) = 
     |       AMIE_EFlux(:,:,:,iBLK) * (1.0e-7 * 100.0 * 100.0)

        do iTime=1,AMIE_ntimes
           do j=1, AMIE_nLats+nCellsPad
              do i=1,AMIE_nMlts
                 if (AMIE_EFlux(i,j,iTime,iBLK) > 100.0) 
     |                AMIE_EFlux(i,j,iTime,iBLK) = 100.0
              enddo
           enddo
        enddo

        AMIE_nLats = AMIE_nLats + nCellsPad

        close(LunEField_)

        deallocate(AllData, stat=iError)

      end subroutine readAMIEoutput

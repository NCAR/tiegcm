      subroutine time_int_to_real(itime, timereal)
      
        implicit none

        integer, dimension(1:12) :: dayofmon
        integer, dimension(1:7) :: itime
        real*8 :: timereal
        integer :: nyear, nleap, nmonth, nday, nhour, nmin, nsec, i

        dayofmon(1) = 31
        dayofmon(2) = 28
        dayofmon(3) = 31
        dayofmon(4) = 30
        dayofmon(5) = 31
        dayofmon(6) = 30
        dayofmon(7) = 31
        dayofmon(8) = 31
        dayofmon(9) = 30
        dayofmon(10) = 31
        dayofmon(11) = 30
        dayofmon(12) = 31

        if (mod(itime(1),4).eq.0) dayofmon(2) = 29

        timereal = 0.0;

        if (itime(1) > 1900) then
           nyear = itime(1) - 1965
        else 
           if (itime(1) > 65) then
              nyear = itime(1) - 65 
           else 
              nyear = itime(1) + 100 - 65
           endif
        endif
        nleap = nyear/4

        nmonth = itime(2) - 1

        nday = 0

        do i=1, nmonth
           nday = nday + dayofmon(i)
        enddo

        nday = nday + itime(3) - 1
        nhour = itime(4)
        nmin = itime(5)
        nsec = itime(6)

        timereal = (dble(nsec) * dble(1.0)) +                  
     |  (dble(nmin) * dble(60.0)) +                       
     |  (dble(nhour) * dble(60.0*60.0)) +                 
     |  (dble(nday) * dble(24.0*60.0*60.0)) +             
     |  (dble(nleap) * dble(24.0*60.0*60.0)) +            
     |  (dble(nyear) * dble(365.0*24.0*60.0*60.0)) +      
     |  itime(7)/1000.0

      end subroutine time_int_to_real

!----------------------------------------------------------------------

      subroutine time_real_to_int(timereal, itime)

        implicit none

        integer, dimension(1:12) :: dayofmon
        integer, dimension(1:7) :: itime
        real*8 :: timereal
        integer :: nyear, nleap, nmonth, nday, nhour, nmin, nsec
        real*8 :: speryear = 31536000.0
        real*8 :: sperday = 86400.0
        real*8 :: sperhour = 3600.0
        real*8 :: spermin = 60.0
        real*8 :: timeleft

        dayofmon(1) = 31
        dayofmon(2) = 28
        dayofmon(3) = 31
        dayofmon(4) = 30
        dayofmon(5) = 31
        dayofmon(6) = 30
        dayofmon(7) = 31
        dayofmon(8) = 31
        dayofmon(9) = 30
        dayofmon(10) = 31
        dayofmon(11) = 30
        dayofmon(12) = 31

        nyear = int(timereal/speryear)
        nleap = nyear/4
        nday = int((timereal - (dble(nyear)*speryear))/sperday)

        if (nday.le.nleap) then
           nyear = int((timereal - (dble(nleap)*sperday))/speryear)
           nleap = nyear/4
           nday = int((timereal - (dble(nyear)*speryear))/sperday)
           if (nday.le.nleap) then
              nyear = int((timereal - (dble(nleap)*sperday))/speryear)
              nleap = nyear/4
              nday = int((timereal - (dble(nyear)*speryear))/sperday)
           endif
        endif

        if (mod((nyear+65),4).eq.0) dayofmon(2) = dayofmon(2) + 1

        nday = nday - nleap

        timeleft = timereal - dble(nyear)*speryear
        timeleft = timeleft - dble(nday+nleap)*sperday

        nhour = int(timeleft/sperhour)
        timeleft = timeleft - dble(nhour)*sperhour

        nmin = int(timeleft/spermin)
        timeleft = timeleft - dble(nmin)*spermin

        nsec = int(timeleft)

        nmonth = 1;

        do while (nday.ge.dayofmon(nmonth))
           nday = nday - dayofmon(nmonth)
           nmonth = nmonth + 1
        end do

        itime(1) = nyear + 1965
        itime(2) = nmonth
        itime(3) = nday + 1
        itime(4) = nhour
        itime(5) = nmin
        itime(6) = nsec
        itime(7) = timeleft - nsec

      end subroutine time_real_to_int

!----------------------------------------------------------------------

      subroutine IO_SetnMLTs(nMLTsIn)
        use ModIE_Interface
        implicit none
        integer, intent(in) :: nMLTsIn
        IOi_NeednMLTs = nMLTsIn
      end subroutine IO_SetnMLTs

!----------------------------------------------------------------------

      subroutine IO_SetnLats(nLatsIn)
        use ModIE_Interface
        implicit none
        integer, intent(in) :: nLatsIn
        IOi_NeednLats = nLatsIn
      end subroutine IO_SetnLats

!----------------------------------------------------------------------

      subroutine IO_SetTime(TimeIn)
!  use ModKind
        use ModIE_Interface
        implicit none
!  Integer, parameter :: Real8_ = selected_real_kind(14,200)
        real (kind=Real8_), intent(in) :: TimeIn
        IOd_NeedTime = TimeIn
      end subroutine IO_SetTime

!----------------------------------------------------------------------

      subroutine IO_SetMLTs(MLTsIn, iError)

        use ModErrors
        use ModIE_Interface

        implicit none

        integer, intent(out) :: iError
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(in) :: 
     |    MLTsIn

        integer :: i,j
!        write(*,*)"MLTsIn =",((MLTsIn(i,j),i=1,25),j=1,50)
        iError = 0
        if (allocated(IOr2_NeedMLTs)) deallocate(IOr2_NeedMLTs)
       allocate(IOr2_NeedMLTs(IOi_NeednMLTs,IOi_NeednLats),stat=iError)
        if (iError /= 0) then
           iError = ecAllocationError_
           return
        else
           do i=1,IOi_NeednMLTs
              do j=1,IOi_NeednLats
                 IOr2_NeedMLTs(i,j) = mod((MLTsIn(i,j)+24.0),24.0)
              enddo
           enddo
        endif
!        write(*,*)"NeedMLTs =",((IOr2_NeedMLTs(i,j),i=1,25),j=1,50)
      end subroutine IO_SetMLTs

!----------------------------------------------------------------------

      subroutine IO_SetLats(LatsIn, iError)

        use ModErrors
        use ModIE_Interface

        implicit none

        integer, intent(out) :: iError
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(in) :: 
     |    LatsIn

        integer :: i,j
!        write(*,*)"LatsIn =",((LatsIn(i,j),i=1,25),j=1,50)
        iError = 0
        if (allocated(IOr2_NeedLats)) deallocate(IOr2_NeedLats)
       allocate(IOr2_NeedLats(IOi_NeednMLTs,IOi_NeednLats),stat=iError)
        if (iError /= 0) then
           iError = ecAllocationError_
           return
        else
           IOr2_NeedLats(1:IOi_NeednMLTs,1:IOi_NeednLats) = 
     |          LatsIn(1:IOi_NeednMLTs,1:IOi_NeednLats)
        endif
!        write(*,*)"NeedLats =",((IOr2_NeedLats(i,j),i=1,25),j=1,50)
      end subroutine IO_SetLats

!----------------------------------------------------------------------

      subroutine IO_SetGrid(MLTsIn, LatsIn, iError)
  
        use ModErrors
        use ModIE_Interface

        implicit none

        integer, intent(out) :: iError
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(in) :: 
     |    MLTsIn, LatsIn
        real, dimension(2) :: IOr1_Location
        real, dimension(5) :: IEr1_Location

        integer :: i,j

        iError = 0

        call IO_SetMLTs(MLTsIn, iError)
        if (iError /= 0) return

        call IO_SetLats(LatsIn, iError)
        if (iError /= 0) return

        if (UseGridBasedIE) then

           if (allocated(IOi3_InterpolationIndices)) 
     |      deallocate(IOi3_InterpolationIndices)
      allocate(IOi3_InterpolationIndices(IOi_NeednMLTs,IOi_NeednLats,3), 
     |     stat=iError)
           if (iError /= 0) then
              iError = ecAllocationError_
              return
           endif

           if (allocated(IOr3_InterpolationRatios)) 
     |     deallocate(IOr3_InterpolationRatios)
      allocate(IOr3_InterpolationRatios(IOi_NeednMLTs,IOi_NeednLats,2),
     |     stat=iError)
           if (iError /= 0) then
              iError = ecAllocationError_
              return
           endif
      write(*,*) IOi_NeednLats,IOi_NeednMLTs
         do i=1,IOi_NeednLats
           do j=1,IOi_NeednMLTs

           IOr1_Location(1) = mod((IOr2_NeedMLTs(j,i) + 24.0),24.0)
           IOr1_Location(2) = IOr2_NeedLats(j,i)

           call IE_FindPoint(IOr1_Location, IEr1_Location, iError)

           if (iError == 0) then
              IOi3_InterpolationIndices(j,i,1:3) = IEr1_Location(1:3)
              IOr3_InterpolationRatios(j,i,1:2)  = IEr1_Location(4:5)
           else
              IOi3_InterpolationIndices(j,i,1) = -1
           endif
           
          enddo
        enddo

        endif
!      write(*,*) "IOi3_1", ((IOi3_InterpolationIndices(i,j,1),i=1,24),
!     |  j=1,50)
!      write(*,*) "IOi3_2", ((IOi3_InterpolationIndices(i,j,2),i=1,24),
!     |  j=1,50)
!      write(*,*) "IOi3_3", ((IOi3_InterpolationIndices(i,j,3),i=1,24),
!     |  j=1,50) 

      end subroutine IO_SetGrid

!----------------------------------------------------------------------

      subroutine IO_GetPotential(PotentialOut, iError)
  
        use ModIE_Interface
        use ModErrors

        implicit none

        integer, intent(out) :: iError
        real, dimension(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs) :: 
     |    IE_Value
        real, dimension(IOi_NeednMLTs,IOi_NeednLats)               :: 
     |    ValueOut
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: 
     |    PotentialOut
        real :: Filler = 0.0

        iError = 0

        if (UseGridBasedIE) then

           IE_Value = IEr3_HavePotential

           call IO_GetValue(IE_Value, ValueOut, Filler, iError)

       if (iError /= 0) then
        write(*,*) "Error in routine IO_GetPotential:"
        write(*,*) cErrorCodes(iError)
        stop
       else
        PotentialOut = ValueOut
       endif

      else

       call IO_GetNonGridBasedPotential(ValueOut, iError)

       if (iError == 0) then
        PotentialOut = ValueOut
       else
        PotentialOut = -1.0e32
       endif

      endif

      end subroutine IO_GetPotential

!----------------------------------------------------------------------
      subroutine IO_GetNonGridBasedPotential(PotentialOut, iError)

        use ModIE_Interface
        use ModErrors

        implicit none

        integer, intent(out)                               :: iError
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out) :: 
     |     PotentialOut

        integer, dimension(7) :: itime
        integer :: iYear, iMonth, iDay, iHour, iMinute, iSecond

        integer :: iChange, iHemisphere

        integer :: iMLT, iLat

        real    :: Potential, ETheta, EPhi, MLT, Lat

        iError = 0

        PotentialOut = 0.0

        if (index(IE_NameOfEFieldModel,'zero') > 0) then
           return
        endif

        if (IOd_NeedTime < 0.0) then
           iError = ecTimeNotSet_
           return
        endif

!  if (index(IE_NameOfEFieldModel,'hpi') <= 0) then
!     if (IOr_NeedIMFBz < -1000.0) then
!        iError = ecIMFBzNotSet_
!        return
!     endif
!     if (IOr_NeedIMFBy < -1000.0) then
!        iError = ecIMFByNotSet_
!        return
!     endif
!  endif

!  if (index(IE_NameOfEFieldModel,'weimer96') > 0) then
!     if (IOr_NeedSWV < -1000.0) then
!        iError = ecSWVNotSet_
!        return
!     endif
!     iChange = 1
!  endif

!  if (index(IE_NameOfEFieldModel,'millstone') > 0) then
!     if (IOr_NeedHPI < -1000.0) then
!        iError = ecHPINotSet_
!        return
!     endif
!  endif

!  if (index(IE_NameOfEFieldModel,'hmr89') > 0) then
!     if (IOr_NeedKP < -1000.0) then
!        iError = ecKPNotSet_
!        return
!     endif
!     iChange = 1
!  endif

        call time_real_to_int(IOd_NeedTime, itime)

        iYear   = itime(1)
        iMonth  = itime(2)
        iDay    = itime(3)
        iHour   = itime(4)
        iMinute = itime(5)
        iSecond = itime(6)

        do iMLT = 1, IOi_NeednMLTs
           do iLat = 1, IOi_NeednLats

              MLT = IOr2_NeedMLTs(iMLT,iLat)
              Lat = IOr2_NeedLats(iMLT,iLat)

              if (IOl_IsNorth) then
                 if (Lat > 0.0) then
                    iHemisphere = 1
                 else
                    Lat = abs(Lat)
                    iHemisphere = -1
                 endif
              else
                 iHemisphere = -1
                 if (Lat < 0.0) Lat = abs(Lat)
              endif

!        if (index(IE_NameOfEFieldModel,'weimer96') > 0) then
!           call start_timing("weimer")
!           call WEIEPOT96(iYear, iMonth, iDay, iHour, iMinute, IOr_NeedSWV, &
!                IOr_NeedIMFBy, IOr_NeedIMFBz, iHemisphere, iChange, &
!                Lat, MLT, ETheta, EPhi, Potential)
!           call end_timing("weimer")
!           Potential = Potential * 1000.0
!           iChange = 0
!        endif

!        if (index(IE_NameOfEFieldModel,'millstone_hpi') > 0) then
!           call MHEMODL(Lat,mlt,IOr_NeedHPI,IOr_NeedIMFBy, &
!                IOr_NeedIMFBz,1,ETheta, EPhi, Potential)
!           Potential = Potential * 1000.0
!        endif

!        if (index(IE_NameOfEFieldModel,'millstone_imf') > 0) then
!           call MHEMODL(lat,mlt,IOr_NeedHPI,IOr_NeedIMFBy, &
!                IOr_NeedIMFBz,2,ETheta, EPhi, Potential)
!           Potential = Potential * 1000.0
!        endif

!        if (index(IE_NameOfEFieldModel,'hmr89') > 0) then
!           call HMREPOT(lat,mlt,IOr_NeedIMFBy,IOr_NeedIMFBz, &
!                IOr_NeedKp,iChange,ETheta, EPhi, Potential)
!           Potential = Potential * 1000.0
!           iChange = 0
!        endif

!        if (index(IE_NameOfEFieldModel,'izmem') > 0) then
!           call IZEPOT(iMonth,lat,mlt, iHemisphere * IOr_NeedIMFBy, &
!                IOr_NeedIMFBz,ETheta, EPhi, Potential)
!           Potential = Potential * 1000.0
!        endif

              PotentialOut(iMLT,iLat) = Potential

           enddo
        enddo

      end subroutine IO_GetNonGridBasedPotential

!----------------------------------------------------------------------

      subroutine IO_GetAveE(AveEOut, iError)
  
        use ModIE_Interface
        use ModErrors

        implicit none

        integer, intent(out) :: iError
        real, dimension(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs) :: 
     |    IE_Value
        real, dimension(IOi_NeednMLTs,IOi_NeednLats)               :: 
     |    ValueOut
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: 
     |    AveEOut

        iError = 0

        if (UseGridBasedIE) then

           IE_Value = IEr3_HaveAveE
           call IO_GetValue(IE_Value, ValueOut, 0.1, iError)

           if (iError /= 0) then
              write(*,*) "Error in routine IO_GetAveE:"
              write(*,*) cErrorCodes(iError)
           else
              AveEOut = ValueOut
           endif

         else

           call IO_GetNonGridBasedAveE(ValueOut, iError)
    
           if (iError == 0) then
              AveEOut = ValueOut
           else
              AveEOut = -1.0e32
           endif

         endif

      end subroutine IO_GetAveE

!----------------------------------------------------------------------

      subroutine IO_GetNonGridBasedAveE(AveEOut, iError)

        use ModIE_Interface
        use ModErrors

        implicit none

        integer, intent(out) :: iError
        real, dimension(IOi_NeednMLTs,IOi_NeednLats) :: AveEOut

        integer :: iMLT, iLat
        real :: Lat, MLT, hpi
        real :: ped, hal, avkev, eflx

        iError = 0
        if (IOr_NeedHPINorm < -1000.0) then
           iError = ecHPINotSet_
           return
        endif

        do iMLT = 1, IOi_NeednMLTs
           do iLat = 1, IOi_NeednLats

              MLT = IOr2_NeedMLTs(iMLT,iLat)
              Lat = IOr2_NeedLats(iMLT,iLat)
              hpi = IOr_NeedHPINorm

               avkev = -1.0e32
!        call get_auroral_conductance(Lat, MLT, hpi, &
!             ped, hal, avkev, eflx)

              AveEOut(iMLT, iLat) = avkev

           enddo
        enddo

      end subroutine IO_GetNonGridBasedAveE

!----------------------------------------------------------------------

      subroutine IO_GetNonGridBasedEFlux(EFluxOut, iError)
!  6/12: Emery - U MI AMIE EFluxOut is in erg/cm2-s or mW/m2 (NOT in W/m2!)
        use ModIE_Interface
        use ModErrors

        implicit none

        integer, intent(out) :: iError
        real, dimension(IOi_NeednMLTs,IOi_NeednLats) :: EFluxOut

        integer :: iMLT, iLat
        real :: Lat, MLT, hpi
        real :: ped, hal, avkev, eflx

        iError = 0

        if (IOr_NeedHPINorm < -1000.0) then
           iError = ecHPINotSet_
           return
        endif

        do iMLT = 1, IOi_NeednMLTs
           do iLat = 1, IOi_NeednLats

              MLT = IOr2_NeedMLTs(iMLT,iLat)
              Lat = IOr2_NeedLats(iMLT,iLat)

               eflx = -1.0e32
!        call get_auroral_conductance(Lat, MLT, IOr_NeedHPINorm, &
!             ped, hal, avkev, eflx)

        ! Need to convert from erg/cm2/s to W/m2
!        EFluxOut(iMLT, iLat) = eflx !* 1.0e-7 * 100.0 * 100.0
              EFluxOut(iMLT, iLat) = -1.0e32

           enddo
        enddo

      end subroutine IO_GetNonGridBasedEFlux

!----------------------------------------------------------------------

      subroutine IO_GetEFlux(EFluxOut, iError)
!  6/12: Emery - U MI AMIE EFluxOut is in erg/cm2-s or mW/m2 (NOT in W/m2!)
        use ModIE_Interface
        use ModErrors

        implicit none

        integer, intent(out) :: iError
        real, dimension(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs) :: 
     |    IE_Value
        real, dimension(IOi_NeednMLTs,IOi_NeednLats)               :: 
     |    ValueOut
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out)  :: 
     |    EFluxOut

        iError = 0

        if (UseGridBasedIE) then

           IE_Value = IEr3_HaveEFlux
           call IO_GetValue(IE_Value, ValueOut, 1.0e-10, iError)

       if (iError /= 0) then
        write(*,*) "Error in routine IO_GetEFlux:"
        write(*,*) cErrorCodes(iError)
        stop
       else
        EFluxOut = ValueOut
       endif

      else

!       call IO_GetNonGridBasedEFlux(ValueOut, iError)
	iError = 0
       if (iError == 0) then
        EFluxOut = ValueOut
       else
        EFluxOut = -1.0e32
       endif

      endif

      end subroutine IO_GetEFlux

!----------------------------------------------------------------------

      subroutine IO_GetValue(ValueIn, ValueOut, Filler, iError)
  
        use ModIE_Interface

        implicit none

        integer, intent(out) :: iError
        real, intent(in)     :: Filler
        real, dimension(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), 
     |   intent(in) :: ValueIn       
        real, dimension(IOi_NeednMLTs,IOi_NeednLats), intent(out) :: 
     |    ValueOut

      integer :: iMLT, iLat, iB, iM, iL
      real    :: dM, dL

      iError = 0

      do iMLT = 1, IOi_NeednMLTs
       do iLat = 1, IOi_NeednLats

        iB = IOi3_InterpolationIndices(iMLT,iLat,1)
        iM = IOi3_InterpolationIndices(iMLT,iLat,2)
        iL = IOi3_InterpolationIndices(iMLT,iLat,3)

        if (iB < 0 .or. iM < 0 .or. iL < 0) then
           ValueOut(iMLT, iLat) =  Filler
        else

           dM = IOr3_InterpolationRatios(iMLT,iLat,1)
           dL = IOr3_InterpolationRatios(iMLT,iLat,2)

           ValueOut(iMLT, iLat) =  
     |           (1.0 - dM) * (1.0 - dL) * ValueIn(iM  , IL  , iB) + 
     |          (1.0 - dM) * (      dL) * ValueIn(iM  , IL+1, iB) + 
     |          (      dM) * (      dL) * ValueIn(iM+1, IL+1, iB) + 
     |          (      dM) * (1.0 - dL) * ValueIn(iM+1, IL  , iB)

        endif

       enddo
      enddo

      end subroutine IO_GetValue

!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine AMIE_SetFileName(cFileNameIn)
        use ModAMIE_Interface
        implicit none
        character (len=100), intent(in) :: cFileNameIn
        AMIE_FileName = cFileNameIn
      end subroutine AMIE_SetFileName

!----------------------------------------------------------------------

      subroutine AMIE_GetFileName(cFileNameOut)
        use ModAMIE_Interface
        implicit none
        character (len=100), intent(out) :: cFileNameOut
        cFileNameOut = AMIE_FileName
      end subroutine AMIE_GetFileName

!----------------------------------------------------------------------

      subroutine AMIE_GetnTimes(nTimesOut)
        use ModAMIE_Interface
        implicit none
        integer, intent(out) :: nTimesOut
        nTimesOut = AMIE_nTimes
      end subroutine AMIE_GetnTimes

!----------------------------------------------------------------------

      subroutine AMIE_GetnMLTs(nMLTsOut)
        use ModAMIE_Interface
        implicit none
        integer, intent(out) :: nMLTsOut
        nMLTsOut = AMIE_nMLTs
      end subroutine AMIE_GetnMLTs

!----------------------------------------------------------------------

      subroutine AMIE_GetnLats(nLatsOut)
        use ModAMIE_Interface
        implicit none
        integer, intent(out) :: nLatsOut
        nLatsOut = AMIE_nLats
      end subroutine AMIE_GetnLats

!----------------------------------------------------------------------

      subroutine AMIE_GetLats(IEi_nMLTs, IEi_nLats, IEi_nBLKs, LatsOut)

        use ModAMIE_Interface

        implicit none

        integer, intent(in) :: IEi_nMLTs, IEi_nLats, IEi_nBLKs
        real, dimension(IEi_nMlts,IEi_nLats,IEi_nBLKs), intent(out) :: 
     |     LatsOut
        integer :: i,j

        do i=1,IEi_nMLTs
           do j=IEi_nLats,1,-1
              LatsOut(i,j,AMIE_North_) = AMIE_Lats(IEi_nLats-j+1)
           enddo
           LatsOut(i,1:IEi_nLats,AMIE_South_) = -AMIE_Lats(1:IEi_nLats)
        enddo

      end subroutine AMIE_GetLats

!----------------------------------------------------------------------

      subroutine AMIE_GetMLTs(IEi_nMLTs, IEi_nLats, IEi_nBLKs, MLTsOut)

        use ModAMIE_Interface

        implicit none

        integer, intent(in) :: IEi_nMLTs, IEi_nLats, IEi_nBLKs
        real, dimension(IEi_nMlts,IEi_nLats,IEi_nBLKs), intent(out) :: 
     |    MLTsOut
        integer :: j

        do j=1,IEi_nLats
           MLTsOut(1:IEi_nMLTs,j,1) = AMIE_MLTs(1:IEi_nMLTs)
           MLTsOut(1:IEi_nMLTs,j,2) = AMIE_MLTs(1:IEi_nMLTs)
        enddo

      
      end subroutine AMIE_GetMLTs

!----------------------------------------------------------------------

      subroutine AMIE_GetPotential(TimeIn, Method, 
     |     IEi_nMLTs, IEi_nLats, IEi_nBLKs, PotentialOut, iError)

        use ModAMIE_Interface
        use ModErrors

        implicit none

      real*8, intent(in) :: TimeIn
!       real, intent(in) :: TimeIn
        integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
        real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: 
     |   PotentialOut
        real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: 
     |   ValueOut
        integer, intent(out) :: iError
!	write(*,*) "1) get_AMIE_Potential"
        AMIE_Value = AMIE_Potential
!        write(*,*) "AMIE_Value",AMIE_Value
        call AMIE_GetValue(TimeIn, Method, 
     |    IEi_nMLTs, IEi_nLats, IEi_nBLKs, ValueOut, iError)
!        write(*,*) ValueOut,iError
        if (iError /= 0) then
           write(*,*) "Error in routine AMIE_GetPotential:"
           write(*,*) cErrorCodes(iError)
           stop
        else
           PotentialOut = ValueOut
        endif

      end subroutine AMIE_GetPotential

!----------------------------------------------------------------------

      subroutine AMIE_GetEFlux(TimeIn, Method, 
     |     IEi_nMLTs, IEi_nLats, IEi_nBLKs, EFluxOut, iError)

        use ModAMIE_Interface
        use ModErrors

        implicit none

        real*8, intent(in) :: TimeIn
!       real, intent(in) :: TimeIn
        integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
        real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: 
     |   EFluxOut
        real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs)              :: 
     |   ValueOut
        integer, intent(out) :: iError

        AMIE_Value = AMIE_EFlux

        call AMIE_GetValue(TimeIn, Method, 
     |       IEi_nMLTs, IEi_nLats, IEi_nBLKs, ValueOut, iError)

        if (iError /= 0) then
           write(*,*) "Error in routine AMIE_GetEFlux:"
           write(*,*) cErrorCodes(iError)
           stop
        else
           EFluxOut = ValueOut
        endif

      end subroutine AMIE_GetEFlux

!----------------------------------------------------------------------

      subroutine AMIE_GetAveE(TimeIn, Method, 
     |   IEi_nMLTs, IEi_nLats, IEi_nBLKs, AveEOut, iError)

        use ModAMIE_Interface
        use ModErrors

        implicit none

        real*8, intent(in) :: TimeIn
!       real, intent(in) :: TimeIn
        integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
        real, dimension(IEi_nMLTs,IEi_nLats, IEi_nBLKs), intent(out) :: 
     |   AveEOut
        real, dimension(IEi_nMLTs,IEi_nLats, IEi_nBLKs)              ::
     |   ValueOut
        integer, intent(out) :: iError

        AMIE_Value = AMIE_AveE
        call AMIE_GetValue(TimeIn, Method, 
     |       IEi_nMLTs, IEi_nLats, IEi_nBLKs, ValueOut, iError)

        if (iError /= 0) then
           write(*,*) "Error in routine AMIE_GetAveE:"
           write(*,*) cErrorCodes(iError)
           stop
        else
           AveEOut = ValueOut
        endif

      end subroutine AMIE_GetAveE

!----------------------------------------------------------------------

      subroutine AMIE_GetValue(TimeIn, Method, 
     |     IEi_nMLTs, IEi_nLats, IEi_nBLKs, ValueOut, iError)

        use ModErrors
        use ModAMIE_Interface

        implicit none

        real*8, intent(in) :: TimeIn
!       real, intent(in) :: TimeIn
        integer, intent(in) :: Method, IEi_nMLTs, IEi_nLats, IEi_nBLKs
        real, dimension(IEi_nMLTs,IEi_nLats,IEi_nBLKs), intent(out) :: 
     |   ValueOut
        integer, intent(out) :: iError

        integer :: iTime, i, j, iLat, iBLK
        logical :: IsDone
        real*8  :: dT, VerySmall = 1.0e-6

        iError = 0
!        write(*,*) "AMIE_GetValue TimeIn Method IE_nMLTs,nLats,nBLKs=",
!     |    TimeIn,Method,IE_NMLTs,IE_nLats,nBLKs,AMIE_South_,AMIE_North
        do iBLK = AMIE_South_, AMIE_North_

           IsDone = .false.
           iTime = 1

       do while (.not. IsDone)
        if (TimeIn - AMIE_Time(iTime,iBLK) < VerySmall) IsDone = .true.
        if ((iTime == AMIE_nTimes) .and. (.not.IsDone)) then
           iTime = iTime + 1 
           IsDone = .true.
        endif
        iTime = iTime + 1
       enddo
!       write(*,*) "1) iTime=",iTime
       if (iTime <= AMIE_nTimes+1) then

        iTime = iTime - 1

        if (iTime == 1) then

           ! If we are before the start time, allow users to extrapolate
           ! up to 5 dT.

           dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)
           if (TimeIn + 5*dt < AMIE_Time(1,iBLK)) then
              ValueOut = -1.0e32
              iError = ecBeforeStartTime_
              return
           endif
        endif
       else
        dT = AMIE_Time(2,iBLK) - AMIE_Time(1,iBLK)

        ! If we are after the end time, allow users to extrapolate
        ! up to 5 dT.

        if (TimeIn - 5*dt < AMIE_Time(AMIE_nTimes,iBLK)) then
           iTime = AMIE_nTimes
        else
           ValueOut = -1.0e32
           iError = ecAfterEndTime_
           return
        endif
       endif
!       write(*,*) "2) iTime=",iTime
       if (Method == AMIE_After_) then
        if (iBLK == AMIE_South_) then
           ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) = 
     |        AMIE_Value(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
           ! Reverse the North block of AMIE data for now...
        else
           do iLat = AMIE_nLats,1,-1
              ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  
     |     AMIE_Value(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
           enddo
        endif
      endif
!       write(*,*) "1) AMIE_Value(1,1)=",AMIE_Value(1,1)
      if (Method == AMIE_Closest_) then
        if (iTime > 1) then
           if (abs(TimeIn-AMIE_Time(iTime,iBLK)) > 
     |               abs(TimeIn-AMIE_Time(iTime-1,iBLK))) 
     |         iTime = iTime - 1
        endif
        if (iBLK == AMIE_South_) then
           ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  
     |          AMIE_Value(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)
        else
           ! Reverse the North block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  
     |     AMIE_Value(1:AMIE_nMLTs, AMIE_nLats - iLat + 1,iTime,iBLK)
           enddo
        endif
       endif
!       write(*,*) "2) AMIE_Value(1,1)=",AMIE_Value(1,1)
       if (Method == AMIE_Interpolate_) then
        ! This will do extrapolation if it is before the first time
        if (iTime == 1) iTime = iTime + 1
        ! dT is the percentage of the way away from the current point
        dT = (AMIE_Time(iTime,iBLK) - TimeIn) / 
     |             (AMIE_Time(iTime,iBLK) - AMIE_Time(iTime-1,iBLK))
        ! Use 1-dT for the selected point, since dt = 0 if you are exactly
        ! on the selected point
        if (iBLK == AMIE_South_) then
           ValueOut(1:AMIE_nMLTs, 1:AMIE_nLats,iBLK) =  
     |    (1.0 - dt)*AMIE_Value(1:AMIE_nMLTs, 1:AMIE_nLats,iTime,iBLK)+
     |           dt*AMIE_Value(1:AMIE_nMLTs, 1:AMIE_nLats,iTime-1,iBLK)
        else
           ! Reverse the 2nd block of AMIE data for now...
           do iLat = AMIE_nLats,1,-1
              ValueOut(1:AMIE_nMLTs, iLat,iBLK) =  
     |          (1.0 - dt)*AMIE_Value(1:AMIE_nMLTs,AMIE_nLats-iLat+1,
     |                                   iTime,iBLK) + 
     |      dt*AMIE_Value(1:AMIE_nMLTs, AMIE_nLats-iLat+1,iTime-1,iBLK)
           enddo
        endif
       endif
!       write(*,*) "3) AMIE_Value(1,1)=",AMIE_Value(1,1)
      enddo

      end subroutine AMIE_GetValue
!----------------------------------------------------------------------

      subroutine get_AMIE_values(rtime)

        use ModIE_Interface

       real*8, intent(in) :: rtime
!       real, intent(in) :: rtime
        integer :: iError
!	write(*,*) "1) get_AMIE_values
        call AMIE_GetPotential(rtime, IE_Interpolate_, IEi_HavenMlts,
     |    IEi_HavenLats, IEi_HavenBLKs, IEr3_HavePotential, iError)

!       call AMIE_GetAveE(rtime, IE_Closest_, IEi_HavenMlts,
!    |    IEi_HavenLats, IEi_HavenBLKs, IEr3_HaveAveE, iError)
        call AMIE_GetAveE(rtime, IE_Interpolate_, IEi_HavenMlts,
     |    IEi_HavenLats, IEi_HavenBLKs, IEr3_HaveAveE, iError)

!       call AMIE_GetEFlux(rtime, IE_Closest_, IEi_HavenMlts,
!    |    IEi_HavenLats, IEi_HavenBLKs, IEr3_HaveEFlux, iError)
        call AMIE_GetEFlux(rtime, IE_Interpolate_, IEi_HavenMlts,
     |    IEi_HavenLats, IEi_HavenBLKs, IEr3_HaveEFlux, iError)

      end subroutine get_AMIE_values

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!\
! This routine finds a point on in the Spherical system, given
! a Theta, Phi: 
! LocIn(1) = Phi
! LocIn(2) = Theta
! It returns a 4 element array:
! LocOut(1) = Index of Block
! LocOut(2) = Index of Longitude
! LocOut(3) = Index of Latitude
! LocOut(4) = Multiplication factor for Longitude
! LocOut(5) = Multiplication factor for Latitude
!/

      subroutine IE_FindPoint(LocIn, LocOut, iError)

        use ModErrors
        use ModIE_Interface

        implicit none

        real, dimension(2), intent(in)  :: LocIn
        real, dimension(5), intent(out) :: LocOut
        integer, intent(out) :: iError
        real :: MLTIn, LatIn
        integer :: j,i, iBLK

        logical :: IsFound

        LocOut = -1.0

        iError = 0

  !\
  ! Check to see if the point is even on the grid.
  !/

        MLTIn = mod(LocIn(1),24.0)

        LatIn = LocIn(2)
!        write (*,*) "MLTIn,LatIn=",MLTIn,LatIn
        if (LatIn > 90.0) then
           LatIn = 180.0 - LatIn
           MLTIn = mod(MLTIn+12.0,24.0)
        endif
        if (LatIn < -90.0) then
           LatIn = -180.0 - LatIn
           MLTIn = mod(MLTIn+12.0,24.0)
        endif

        if (MLTIn > 24.0 .or. MLTIn < 0 .or. LatIn > 90.0 
     |    .or. LatIn < -90.0) then
           iError = ecPointOutofRange_
!        write (*,*) "24<MLTIn<0,-90>LatIn>90=",MLTIn,LatIn
           return
        endif

        iBLK = 1
        do while (iBLK <= IEi_HavenBLKs)
           j = 1
           do while (j < IEi_HavenMLTs)
              i = 1
              do while (i < IEi_HavenLats)

           !\
           ! Check to see if the point is within the current cell
           !/

!           if (LatIn <  IEr3_HaveLats(j,i+1,iBLK) .and. &
!               LatIn >= IEr3_HaveLats(j,i,iBLK) .and. &
!               MLTIn <  IEr3_HaveMLTs(j+1,i,iBLK) .and. &
!               MLTIn >= IEr3_HaveMLTs(j,i,iBLK)) then
! 6/12 Emery:  Needed <= for first or +90N called missing (so zero poten, AveE and Eflux)
           if (LatIn <= IEr3_HaveLats(j,i+1,iBLK) .and. 
     |         LatIn >= IEr3_HaveLats(j,i,iBLK) .and. 
     |         MLTIn <= IEr3_HaveMLTs(j+1,i,iBLK) .and. 
     |         MLTIn >= IEr3_HaveMLTs(j,i,iBLK)) then

              !\
              ! If it is, then store the cell number and calculate
              ! the interpolation coefficients.
              !/

              LocOut(1) = iBLK
              LocOut(2) = j
              LocOut(3) = i

              LocOut(4) = (MLTIn            -IEr3_HaveMLTs(j,i,iBLK))/
     |            (IEr3_HaveMLTs(j+1,i,iBLK)-IEr3_HaveMLTs(j,i,iBLK))
              LocOut(5) = (LatIn            -IEr3_HaveLats(j,i,iBLK))/
     |            (IEr3_HaveLats(j,i+1,iBLK)-IEr3_HaveLats(j,i,iBLK))

              iBLK = IEi_HavenBLKs
              j = IEi_HavenMLTs
              i = IEi_HavenLats
      else

           endif

           i = i + 1

        enddo

        j = j + 1

           enddo

           iBLK = iBLK + 1

        enddo

      end subroutine IE_FindPoint
!----------------------------------------------------------------------

      subroutine IE_set_inputs(StringInputLines)

        use ModIE_Interface
        use ModFiles
        use ModExtras

        implicit none

        character (len=100),dimension(*),intent(in) :: StringInputLines
        character (len=100) :: StringLine
        logical :: IsDone
        integer :: iLine, iDebugProc

        iLine = 1

        IsDone = .false.

        do while (.not.IsDone)

           StringLine = StringInputLines(iLine)

           if (StringLine(1:1) == "#") then

        if (index(StringLine,"#BACKGROUND") > 0) then

           call read_in_string(IE_NameOfModelDir)
           call read_in_string(IE_NameOfEFieldModel)
           call read_in_string(IE_NameOfAuroralModel)
           call read_in_string(IE_NameOfSolarModel)

           if (index(IE_NameOfAuroralModel,'IHP') > 0) 
     |                IE_NameOfAuroralModel = 'ihp'
           if (index(IE_NameOfAuroralModel,'PEM') > 0) 
     |                IE_NameOfAuroralModel = 'pem'

           if (index(IE_NameOfEFieldModel,'AMIE') > 0) 
     |                IE_NameOfEFieldModel = 'amie'

           if (index(IE_NameOfEFieldModel,'weimer01') > 0) 
     |          IE_NameOfEFieldModel = 'weimer01'
           if (index(IE_NameOfEFieldModel,'Weimer01') > 0) 
     |          IE_NameOfEFieldModel = 'weimer01'
           if (index(IE_NameOfEFieldModel,'WEIMER01') > 0) 
     |          IE_NameOfEFieldModel = 'weimer01'

           if (index(IE_NameOfEFieldModel,'weimer') > 0 .and. 
     |          index(IE_NameOfEFieldModel,'01') == 0) 
     |          IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'Weimer') > 0 .and. 
     |          index(IE_NameOfEFieldModel,'01') == 0)
     |          IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'WEIMER') > 0 .and. 
     |          index(IE_NameOfEFieldModel,'01') == 0) 
     |          IE_NameOfEFieldModel = 'weimer96'

           if (index(IE_NameOfEFieldModel,'weimer96') > 0) 
     |          IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'Weimer96') > 0) 
     |          IE_NameOfEFieldModel = 'weimer96'
           if (index(IE_NameOfEFieldModel,'WEIMER96') > 0) 
     |          IE_NameOfEFieldModel = 'weimer96'

           if (index(IE_NameOfEFieldModel,'SAMIE') > 0) 
     |          IE_NameOfEFieldModel = 'samie'

           UseGridBasedIE = .false.

        endif

        if (index(StringLine,"#AMIEFILES") > 0) then
           call read_in_string(AMIEFileNorth)
           call read_in_string(AMIEFileSouth)
           IE_NameOfEFieldModel = "amie"
           UseGridBasedIE = .true.
        endif

        if (index(StringLine,"#DEBUG") > 0) then
           call read_in_int(iDebugLevel)
           call read_in_int(iDebugProc)
           if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
              iDebugLevel = -1
           endif
        endif

        if (index(StringLine,"#END") > 0) then
           IsDone = .true.
        endif

        if (iLine >= MaxInputLines) then
           IsDone = .true.
        endif

           else

              iLine = iLine + 1

           endif

        enddo

      contains

        subroutine read_in_int(variable)
          integer, intent(out) :: variable
          iline = iline + 1
          read(StringInputLines(iline),*) variable
        end subroutine read_in_int

        subroutine read_in_logical(variable)
          logical, intent(out) :: variable
          iline = iline + 1
          read(StringInputLines(iline),*) variable
        end subroutine read_in_logical

        subroutine read_in_string(variable)
          character (len=100), intent(out) :: variable
          iline = iline + 1
          variable = StringInputLines(iline)
        end subroutine read_in_string

        subroutine read_in_real(variable)
          real :: variable
          iline = iline + 1
          read(StringInputLines(iline),*) variable
        end subroutine read_in_real

        subroutine read_in_time(variable)
          real*8 :: variable
          integer, dimension(7) :: itime
          iline = iline + 1
          read(StringInputLines(iline),*) itime(1:6)
          itime(7) = 0
          call time_int_to_real(itime, variable)
        end subroutine read_in_time

      end subroutine IE_set_inputs
!----------------------------------------------------------------------
      subroutine IE_Initialize(iOutputError)

        use ModErrors
        use ModIE_Interface
        use ModFiles
        use ModIoUnit, only : UnitTmp_

        implicit none

        integer, intent(out) :: iOutputError
        character (len=100)  :: inFileName
        integer              :: iError

        integer, parameter  :: South_ = 1
        integer, parameter  :: North_ = 2

        logical :: IsFound_EFieldModel

        iError = 0
        iOutputError = 0

        IsFound_EFieldModel = .false.

        call set_error_codes

  !\
  ! --------------------------------------------------------------------
  ! Electric Field Models
  ! --------------------------------------------------------------------
  !/

        if (iDebugLevel > 1) 
     |       write(*,*) "==> Efield Model : ",IE_NameOfEFieldModel

        if (iDebugLevel > 1) 
     |       write(*,*) "==> Model Directory : ",IE_NameOfModelDir

        if (index(IE_NameOfEFieldModel,'zero') > 0) then
           IsFound_EFieldModel = .true.
        endif

        LunEField_      = UnitTmp_
        LunConductance_ = UnitTmp_

!  if (index(IE_NameOfEFieldModel,'weimer96') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, weimer96_file)
!     open(LunEField_,file=weimer96_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',weimer96_file
!        iOutputError = ecFileNotFound_
!     endif
!     call ReadCoef96(LunEField_)
!     close(LunEField_)
!  endif

!  if (index(IE_NameOfEFieldModel,'weimer01') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, weimer01_file)
!     open(LunEField_,file=weimer01_file,status='old',&
!          form='unformatted', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',weimer01_file
!        iOutputError = ecFileNotFound_
!     endif
!     call ReadCoef01(LunEField_)
!     close(LunEField_)
!  endif

!  if (index(IE_NameOfEFieldModel,'samie') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, stat_amie_file)
!     open(LunEField_,file=stat_amie_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',stat_amie_file
!        iOutputError = ecFileNotFound_
!     endif
!     call read_amies(LunEField_)
!     close(LunEField_)
!  endif

!  if (index(IE_NameOfEFieldModel,'millstone_hpi') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, millstone_hill_i_file)
!     open(LunEField_,file=millstone_hill_i_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',millstone_hill_i_file
!        iOutputError = ecFileNotFound_
!     endif
!     call mhinit(1, LunEField_, 1, iDebugLevel)
!     close(LunEField_)
!  endif

!  if (index(IE_NameOfEFieldModel,'millstone_imf') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, millstone_hill_s_file)
!     open(LunEField_,file=millstone_hill_s_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',millstone_hill_s_file
!        iOutputError = ecFileNotFound_
!     endif
!     call mhinit(2, LunEField_, 1, iDebugLevel)
!     close(LunEField_)
!  endif

!  if (index(IE_NameOfEFieldModel,'hmr89') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, hepner_maynard_file)
!     open(LunEField_,file=hepner_maynard_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',hepner_maynard_file
!        iOutputError = ecFileNotFound_
!     endif
!     call gethmr(LunEField_)
!     close(LunEField_)
!  endif

!  if (index(IE_NameOfEFieldModel,'izmem') > 0) then
!     IsFound_EFieldModel = .true.
!     call merge_str(IE_NameOfModelDir, izmem_file)
!     open(LunEField_,file=izmem_file,status='old', iostat = iError)
!     if (iError /= 0) then
!        write(6,*) 'Error opening file :',izmem_file
!        iOutputError = ecFileNotFound_
!     endif
!     call izinit(LunEField_)
!     close(LunEField_)
!  endif

  !\
  ! --------------------------------------------------------------------
  ! Conductance Models
  ! --------------------------------------------------------------------
  !/

        if (iDebugLevel > 1) 
     |    write(*,*) "==> Conductance Model : ",IE_NameOfAuroralModel

!  if (index(IE_NameOfAuroralModel,'ihp') > 0)  &
!       call read_conductance_model(iError)
!  if (index(IE_NameOfAuroralModel,'pem') > 0)  &
!       call read_conductance_model(iError)

        if (iDebugLevel > 4) 
     |    write(*,*) "=====> Back from read conductance"

        if (index(IE_NameOfEFieldModel,'amie') > 0) then

           IsFound_EFieldModel = .true.

           call AMIE_SetFileName(AMIEFileNorth)
           call readAMIEOutput(North_, .false., iError)

      write(*,*) "AMIEFileSouth : ",AMIEFileSouth, 
     |   index(AMIEFileSouth,'mirror')

      if (index(AMIEFileSouth,'mirror') > 0) then
        call AMIE_SetFileName(AMIEFileNorth)
        call readAMIEOutput(South_, .true., iError)
      else
        call AMIE_SetFileName(AMIEFileSouth)
        call readAMIEOutput(South_, .false., iError)
      endif

      call AMIE_GetnLats(IEi_HavenLats)
      call AMIE_GetnMLTs(IEi_HavenMLTs)
      IEi_HavenBLKs = 2

      if (iDebugLevel > 1) then
        write(*,*) "=> IEi_HavenBLKs : ", IEi_HavenBLKs
        write(*,*) "=> IEi_HavenLats : ", IEi_HavenLats
        write(*,*) "=> IEi_HavenMLTs : ", IEi_HavenMLTs
      endif

      allocate(IEr3_HaveLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs)
     |          ,stat=iError)
      if (iError /= 0) then
       write(*,*) "Error in allocating array IEr3_HaveLats in Interface"
        stop
      endif

      allocate(IEr3_HaveMlts(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs), 
     |          stat=iError)
      if (iError /= 0) then
       write(*,*) "Error in allocating array IEr3_HaveMlts in Interface"
        stop
      endif

      allocate(IEr3_HavePotential(IEi_HavenMlts,IEi_HavenLats,
     |  IEi_HavenBLKs), stat=iError)
      if (iError /= 0) then
        write(*,*) "Error in allocating array IEr3_HavePotential ",
     |   "in Interface"
        stop
      endif

      allocate(IEr3_HaveEFlux(IEi_HavenMlts,IEi_HavenLats,
     |  IEi_HavenBLKs), stat=iError)
      if (iError /= 0) then
        write(*,*) "Error in allocating array IEr3_HaveEFlux ",
     |   "in Interface"
        stop
      endif

      allocate(IEr3_HaveAveE(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs),
     |          stat=iError)
      if (iError /= 0) then
       write(*,*)"Error in allocating array IEr3_HaveAveE in Interface"
        stop
      endif

      call AMIE_GetLats(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,
     |          IEr3_HaveLats,iError)

      call AMIE_GetMLTs(IEi_HavenMlts,IEi_HavenLats,IEi_HavenBLKs,
     |         IEr3_HaveMLTs,iError)

      return

      endif

      if (.not.IsFound_EFieldModel) then
       iOutputError = ecEFieldModelNotFound_
      endif

      if (iDebugLevel > 3) write(*,*) "====> Done with IE_Initialize"

      end subroutine IE_Initialize
!----------------------------------------------------------------------

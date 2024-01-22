#ifdef GAMERA
!Various simple routines for managing time/dates
        module dates
            implicit none

            real*8, parameter, private :: mjdScl = 1.0/(60.0*60.0*24.0)

            integer, parameter, private :: i365  =  365, i1582 = -418

            contains

            !Return MJD given t as seconds elapsed from mjd0
            function T2MJD(t,mjd0) result(mjd)
                real*8, intent(in) :: t,mjd0
                real*8 :: mjd

                mjd = mjd0 + t*mjdScl
            end function T2MJD

            subroutine DateTimeStr(dtStr)
                character(len=*), intent(inout) :: dtStr
                call FDATE(dtStr)
            end subroutine DateTimeStr
        !-----------------
        !NOTE: These routines modified from LFM-RCM code by K. Sorathia

            function isLeapYear(iyr) result(isLeap)
                integer, intent(in) :: iyr
                logical :: isLeap

                isLeap = (iyr/4*4 == iyr) .and. (iyr .LT. 1582 .OR.
     |           iyr/100*100 .NE. iyr .OR. iyr/400*400 .EQ. iyr)

            end function isLeapYear
            !------
            !mjd2ut
            !PURPOSE: Convert modified julian date into Yr,Doy,Mon,Day,Hr,Min,Sec
            !INPUTS:
            !   MJD  REAL*8 Modified Julian Date
            !
            ! OUTPUTS:
            !   iYr  INTEGER Year for MJD
            !   iDoy INTEGER Day of Year for MJD
            !   iMon INTEGER Month of Year for MJD
            !   iDay INTEGER Day of Month for MJD
            !   iHr  INTEGER Hour of Day for MJD
            !   iMin INTEGER Minute of Day for MJD
            !   rSec REAL  Second of Minute for MJD

            subroutine mjd2ut(mjd,iyr,idoy,imon,iday,ihr,imin,rsec)
                real*8, intent(inout) :: mjd
                real*8, intent(inout) :: rSec
                integer,  intent(inout) :: iyr,idoy,imon,iday,ihr,imin

                integer :: id
                real*8 :: rdoy,rsod,rsid,jd,jepoch
                character(len=3) :: camon

                camon = 'XXX'
                id = 11
                call julian(id,iyr,rdoy,mjd,jepoch)
                id = 1
                idoy = rdoy
                call date_doy(id,iyr,camon,imon,iday,idoy)

                rsid = 86400.0
                rsod = (mjd-int(mjd))*rsid
                id = 1
                call hms_sod(id,ihr,imin,rsec,rsod)

            end subroutine mjd2ut

            !------
            !ut2mjd
            !PURPOSE: Convert Yr,Doy,Mon,Day,Hr,Min,Sec into modified julian date
            !INPUTS:
            !   iYr  INTEGER Year for MJD
            !   iDoy INTEGER Day of Year for MJD
            !   iMon INTEGER Month of Year for MJD
            !   iDay INTEGER Day of Month for MJD
            !   iHr  INTEGER Hour of Day for MJD
            !   iMin INTEGER Minute of Day for MJD
            !   rSec REAL*8  Second of Minute for MJD
            !
            !OUTPUTS:
            !   MJD  REAL*8 Modified Julian Date

            subroutine ut2mjd(mjd,iyr,idoy,imon,iday,ihr,imin,rsec)
                real*8, intent(inout) :: mjd
                real*8, intent(inout) :: rSec
                integer,  intent(inout) :: iyr,idoy,imon,iday,ihr,imin

                integer :: id
                real*8 :: rdoy,rsod,rsid,jd,jepoch
                character(len=3) :: camon

                camon = 'XXX'
                id=0
                call date_doy(id,iyr,camon,imon,iday,idoy)
                call hms_sod(id,ihr,imin,rsec,rsod)
                rsid = 86400.0
                rdoy = idoy*1.0+rsod/rsid

                id=10
                call julian(id,iyr,rdoy,mjd,jepoch)

            end subroutine ut2mjd

            !------
            !date_doy
            !PURPOSE: Converts dates (iYr,iMon,iDay) to daynumber iDoy and v.v.
            ! INPUTS:
            !   ID  INTEGER     ID=0  date (yr,month,day) ----> doy of year
            !               ID=1  doy of year ----> date (yr,month,day)
            !   iYr INTEGER     year; the year xxxBC should be entered as -xxx+1
            !
            !   ID = 0:
            !    cMon    CHARACTER*3    month   (either cMon or iMon must be specified)
            !    iMon    INTEGER    value 1-12   (see PROCEDURE)
            !    iDay    INTEGER    day of month
            !
            !   ID = 1:
            !    iDoy    INTEGER    day of year
            ! OUTPUTS:
            !   ID = 0:
            !    iDoy    INTEGER    day of year
            !
            !   ID = 1:
            !    cMon    CHARACTER*3    month
            !    iMon    INTEGER    1-12
            !    iDay    INTEGER    day of month
            subroutine date_doy(id,iyr,cmon,imon,iday,idoy)
                integer, intent(in) :: id,iyr
                integer, intent(inout) :: imon,iday,idoy
                character(len=3), intent(inout) :: cmon

                logical :: isLeap

                character(len=3) :: MONTHS(12) = [
     |  'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP',
     |  'OCT','NOV','DEC']
                integer          ::  IDAYS(12) = [
     |  31   ,28   ,31   ,30   ,31   ,30   ,31   ,31   ,30   ,
     |  31   ,30   ,31   ]

                integer :: leap, i
              if (id .NE. 0 .AND. id .NE. 1) then
                  write(*,*) 'DATE_DOY: invalid id'
                  STOP
              endif
              leap = 2

              isLeap = isLeapYear(iyr)
              if (isLeap) leap = 1

              IDAYS(2) = 30-leap            !Days in february

              if (id == 0) then !Date ---> daynumber
                  i = 1
                  do while (i <= 12 .and. cmon /= MONTHS(i))
                      i = i+1
                  enddo
                  if (i <= 12) then
                      imon = i
                  else !If cmon does not make sense, try imon
                      if (imon .LT. 1 .OR. imon .GT. 12) then
                          write(*,*) 'DATE_DOY: invalid month'
                      endif
                      cmon = MONTHS(imon)
                  endif

                  if (iday < 1 .or. iday > IDAYS(imon)) then
                      write(*,*) 'DATE_DOY: invalid iday'
                      stop
                  endif
                  i = imon-1
                  idoy = i*31-MIN(i/2,3)-MAX(0,(i-7)/2)+iday !Day number
                  if (i >= 2) idoy = idoy-leap !Correction for February

                  if (iyr == 1582 .and. idoy > 278) then
                      if (idoy < 288) then
                  write(*,*) 'DATE_DOY: date in gap between calendars'
                        stop
                      endif
                      idoy = idoy-10
                  endif
              else !Daynumber ---> date
                  if (idoy < 1 .or. idoy > 367-leap) then
                      write(*,*) 'DATE_DOY: invalid idoy'
                      stop
                  endif
                  if (iyr == 1582) then
                      if (idoy > 355) then
                          write(*,*) 'DATE_DOY: invalid idoy for 1582'
                          stop
                      endif
                      IDAYS(10) = 21 !Subtract 10 day gap for October
                  endif
                  iday = idoy
                  imon = 1
                  do while ((iday-1)/IDAYS(imon) >= 1)
                      iday = iday-idays(imon)
                      imon = imon+1
                  enddo
                  cmon = MONTHS(imon)
                  if (iyr == 1582 .and. imon == 10 .and. iday >= 5) then
                      iday = iday+10
                  endif
              endif
            end subroutine date_doy

            !------
            !hms_sod
            !PURPOSE: Converts Time (Hr,Min,Sec) into Seconds of Day and v.v.
            ! INPUTS
            !   ID      INTEGER     ID=0 Time (Hr,Min,Sec) ---> Second of Day
            !               ID=1 Second of Day ---> Time (Hr,Min,Sec)
            !
            !   ID=0
            !     ihr   INTEGER     Hour of Day
            !     imin  INTEGER     Minute of Hour
            !     rsec      REAL*8      Second of Minute
            !
            !   ID=1
            !     rsod  REAL*8      Second of Day
            !
            ! OUTPUTS
            !   ID=0
            !     rsod  REAL*8      Seconds of Day
            !
            !   ID=1
            !     ihr       INTEGER         Hour of Day
            !     imin      INTEGER         Minute of Hour
            !     rsec      REAL*8          Second of Minute

            subroutine hms_sod(id,ihr,imin,rsec,rsod)
                integer, intent(in) :: id
                integer, intent(inout) :: ihr,imin
                real*8, intent(inout) :: rsec,rsod

                real*8 :: rSecInHr,rSecInMin,rtmp

                rSecInHr = 3600.0
                rSecInMin = 60.0

                if (id == 0) then
                    rsod = ihr*rSecInHr+imin*rSecInMin+rsec
                endif

                if (id == 1) then
                    ihr = INT(rsod/rSecInHr)
                    imin = INT((rsod-rSecInHr*ihr)/rSecInMin)
                    rtmp = rsod-rsecinhr*ihr-rSecInMin*imin
                    rsec = rtmp
                endif

            end subroutine hms_sod

            !------
            !julian
            !PURPOSE: Of the four quantities, date, Julian day number, Julian and Besselian epoch, one is supplied. The other three are then calculated.
            ! INPUTS:
            !   IDIN    integer     ID=0  date --> JD, Julian epoch, Besselian epoch
            !               =1  JD   --> Date, Julian epoch, Besselian epoch
            !               =2  Julian epoch   --> Date, JD, Besselian epoch
            !       (not active)    =3  Besselian epoch --> Date, JD, Julian epoch
            !               Add 10: use modified Julian days (=JD-2400000.5)
            !               Add 20: use days relative to Jan. 1, 2000
            !                                       (=JD-2451544.5)
            !  INPUTS/OUTPUTS: (depending on value of ID)
            !   iYr integer     year; the year xxxBC should be entered as -xxx+1
            !   Doy real*8      day of year (including fraction for the time of
            !                                            day).
            !   JD  real*8      Julian day
            !   JEpoch  real*8      Julian epoch = time in Julian years
            !   BEpoch  real*8      Besselian epoch

            subroutine julian(idin,iyr,doy,jdio,jepoch)
                integer, intent(in) :: idin
                integer, intent(inout) :: iyr
                real*8, intent(inout) :: doy,jdio,jepoch

                integer :: id,laps,i,j,leap,jd0,idoy
                real*8 :: jd,bepoch,jd2000,frac


                jd2000 = 2451545    ! Jan 1.5, 2000 (noon)
                if (idin >= 10) jd2000 = 51544.5d0
                if (idin >= 20) jd2000 = 0.5d0

                id = MOD(idin,10)   ! ID=0,1,2,3

                if (id == 0) then ! Date ---> Julian day, etc.
                    laps = iyr-2000 ! Whole years from 1 Jan 1.0d 2000 to idem of iYr
                    i = laps
                    if (laps > 0) i = laps-1
                    j = MAX(i1582,i)
                    leap = i/4-j/100+j/400
                    if (laps > 0) leap = leap+1
                    jd = i365*laps+leap+dble(doy-1.5)   ! Rel. to 2000 Jan. 1.5 (noon)
                    if (iyr <= 1582) jd = jd+10
                    jepoch = 2000+jd/365.25
                    bepoch = 1900+(jd+36524.68648)/365.242198781
                else ! Julian day, etc. ---> Date
                    if (id == 2) then ! Julian epoch
                        jd = jepoch-2000
                        bepoch = 2000.0012775+1.000021359*jd
                        jd = jd*365.25
                    else if (id == 3) then
                        jd = bepoch-1900
                        jepoch = 1900.000858+.999978642*jd
                        jd = -36524.68648+jd*365.242198781
                    else           ! Assume ID=1
                        jd = jdio-jd2000
                        jepoch = 2000+jd/365.25
                        bepoch = 1900+(jd+36524.68648)/365.242198781
                    endif


                    jd = jd+0.5D0        ! Shift to previous midnight
                    jd0 = INT(jd)        ! Whole days between previous midnight ...
                    if (jd < jd0) jd0=jd-1 ! .. and midnight of 1 Jan 2000

                    frac = jd-jd0        ! Day fraction from previous midnight
                    laps = 0
                    idoy = jd0
                    do while (abs(idoy) > 365)
                        laps = laps+idoy/366   ! Underestimate number of full years in iDoy
                        i = laps
                        if (laps > 0) i = laps-1
                        j = MAX(i1582,i)
                        leap = i/4-j/100+j/400
                        if (laps > 0) leap = leap+1
                        idoy = jd0-laps*i365-leap
                        if (laps <= i1582) idoy = idoy-10
                    enddo



                    !-----
                    ! Leave do while loop with -365<=iDoy<=365. Doy <= 0 happens for dates earlier
                    ! than 2000 Jan. 1.0; Doy >=0 for dates later than 2000 Jan 1.0.
                    iyr = 2000+laps
                    if (idoy < 0) then
                        iyr = iyr-1
                        if (iyr == 1582) then
                            idoy = idoy+355
                            if (idoy < 0) iyr=iyr-1
                        endif
                    endif

                    leap = 0 ! Suppose it's not a leap year
                    if (isLeapYear(iyr)) leap=1

                    if (idoy < 0) idoy = 365+leap+idoy    ! Now iDoy>=0
                    if (idoy == 365+leap) then
                        idoy = 0
                        iyr = iyr+1
                    endif

                    doy = idoy+1+frac
                    jd = jd-0.5d0
                endif

                jdio = jd2000+jd

            end subroutine julian

        end module dates

#else
!-----------------------------------------------------------------------
      ! Intel Fortran compiler chokes on empty source files.
      ! This subroutine is empty so this file will have SOMETHING in it
      subroutine dates_null
      end subroutine dates_null
!-----------------------------------------------------------------------
#endif

! ====================================================================
!   Subroutine Datediff
!     
!   Purpose: Calculation the difference of two date [d].
! ====================================================================
! =========================Incoming variables=========================
!   date1       The first date [yyyymmdd].
!   date2       The second date [yyyymmdd].
! =========================Outcoming variables========================
!   nd          The different day between date1 and date2 [d].
! ====================================================================
    SUBROUTINE Datediff(date1,date2,nd)
    IMPLICIT NONE
    INTEGER (KIND=4) :: nd
    INTEGER (KIND=4) :: jDATE
    CHARACTER (LEN=8) :: date1,date2

    nd = jDATE(date2)-jDATE(date1)

    RETURN
    END SUBROUTINE Datediff

! ====================================================================
!   Function ifrun
!     
!   Purpose: IF leap year or not.
! ====================================================================
! =========================Incoming variables=========================
!   year        The imput year [yyyymmdd].
! =========================Outcoming variables========================
!   ifrun       1 for leap year; 0 for others.
! ====================================================================
    INTEGER Function ifrun(year)
    IMPLICIT NONE
    INTEGER (KIND=4) :: year
    
    IF (mod(year,4)==0 .and. mod(year,100).ne.0) THEN
        ifrun=1
    ELSEIF (mod(year,400)==0 .and. mod(year,100)==0) THEN
        ifrun=1
    ELSE
        ifrun=0
    ENDIF
    
    RETURN
    END Function ifrun

! ====================================================================
!   Function JDate
!     
!   Purpose: Calculation the Julian day with the beginning of 1960-1-1
!            When the date is 1960-1-1, the result is 1.
! ====================================================================
! =========================Incoming variables=========================
!   date        The imput date [yyyymmdd].
! =========================Outcoming variables========================
!   JDate       The different beween for the date and 1960-1-1.
! ====================================================================
    Function JDate(date)
    IMPLICIT NONE
    CHARACTER (LEN=8) :: date
    INTEGER (KIND=4) :: nY, nM, nD, numrun, nmd, nyd, i
    INTEGER (KIND=4) :: ifrun, JDate

    READ(date(1:4),"(i4)")nY
    READ(date(5:6),"(i2)")nM
    READ(date(7:8),"(i2)")nD
      
    IF (ny<1960) PRINT*,"error in year chioces"  
    nyd=ny-1960 ! Difference of year.
    numrun=0
    
    DO i=1960,nY
        numrun=numrun+ifrun(i)
    ENDDO
      
    SELECT CASE(nm-1)
    CASE(0)
        nmd=0
    CASE(1)
        nmd=31
    CASE(2)
        nmd=59
    CASE(3)
        nmd=90
    CASE(4)
        nmd=120
    CASE(5)
        nmd=151
    CASE(6)
        nmd=181
    CASE(7)
        nmd=212
    CASE(8)
        nmd=243
    CASE(9)
        nmd=273
    CASE(10)
        nmd=304
    CASE(11)
        nmd=334	    
    END SELECT
    
    IF (ifrun(ny).and.nm<3) nmd=nmd-1 ! leap year does not have 2-29
    Jdate=nyd*365+numrun+nmd+nd
    
    END Function JDate

! ====================================================================
!   SUBROUTINE DateAdd
!     
!   Purpose: The new date when a date adds integer days.
! ====================================================================
! =========================Incoming variables=========================
!   date        The imput date [yyyymmdd].
!   interval    The integer that will be added to the date [d].
! =========================Outcoming variables========================
!   date1       the new date which has been added [yyyymmdd].
! ====================================================================
    SUBROUTINE DateAdd(date,interval,date1)
    IMPLICIT NONE

    CHARACTER (LEN=8) :: date
    CHARACTER (LEN=8) :: date1
    INTEGER (kind=4) :: ny, nm, nd, n, interval, i, ifrun
    
    READ(date(1:4),"(i4)")ny
    READ(date(5:6),"(i2)")nm
    READ(date(7:8),"(i2)")nd
    n=interval
    DO i=1,n
        nd=nd+1
        IF (nd>31.and.(nm.eq.1.or.nm.eq.3.or.nm.eq.5.or.nm.eq.7.or.nm.eq.8.or.nm.eq.10.or.nm.eq.12)) THEN
            nd=1
            nm=nm+1
        ELSEIF (nd>30.and.(nm.eq.4.or.nm.eq.6.or.nm.eq.9.or.nm.eq.11)) THEN
            nd=1
            nm=nm+1
        ELSEIF (nm==2.and.ifrun(ny).and.nd>29) THEN
            nd=1
            nm=nm+1
        ELSEIF (nm==2.and.ifrun(ny)==0.and.nd>28) THEN
            nd=1
            nm=nm+1
        ENDIF
!       Month add.
        
        IF (nm>12) THEN
            nm=1
            ny=ny+1
        ENDIF
    ENDDO
    WRITE(date1,"(i4.4,i2.2,i2.2)")ny,nm,nd
    
    END SUBROUTINE DateAdd

! ====================================================================
!   SUBROUTINE Tcontrol
!     
!   Purpose: t+dt.
!            Considering print time, boundary condition time and
!            endding time, self-adaption the dt.
! ====================================================================
    SUBROUTINE Tcontrol
    USE parm
    IMPLICIT NONE
    INTEGER (KIND=KI) :: i
    REAL (KIND=KR) :: dt1

    IF (t+dt>TAtm .and. t<TAtm-1E-8) THEN
        dt1 = TAtm - t
        dt = min(dt,dt1)
    ENDIF

    IF (t+dt>TPrint(PLevel) .and. t<(TPrint(PLevel)-1E-8)) THEN
        dt1 = TPrint(PLevel) - t
        dt = min(dt,dt1)
    ENDIF
    
    IF(t+dt>Tend .and. t<Tend) THEN
        dt1 = Tend-t
        dt = min(dt,dt1)
    ENDIF
!    t = t+dt
    
    END SUBROUTINE Tcontrol


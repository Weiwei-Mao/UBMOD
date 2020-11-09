! ====================================================================
!   Subroutine Selector_In   
!     
!   Purpose: Read information from file selector.in
! ====================================================================
! =========================Incoming variables=========================
!   None
! =========================Outcoming variables========================
!   Hed         Character in the screen.
!   LUnit       Unit of length.
!   TUnit       Unit of time.
!   MUnit       Unit of mass.
!   xConv       The conversion coefficient of length.
!   tConv       The conversion coefficient of time.
!   mConv       The conversion coefficient of mass.
!   AtmBC        If meteorological data to calculate the ET.
!   Bup         Flux upper boundary condition.
!   Bdn         Flux lower boundary condition.
!   lchem       If calculate the solute.
!   parredis    The drainage function. 
!   Dfit        The empirical formula of Diffusion.
!   Nmat        The number of material.[-]
!   NPar        Number of unsaturated hydraulic parameters.
!   par(:,:)	The unsaturated parameter. 
!   thf(:)      The field capacity.[-]
!   thw(:)      The wilting point.[-]
!   sp(4,:)     The discreted function of evaporation.[-]
!   ths(:)      The saturated water content.
!   dt          The time step.[d]
!   ddn         Diffusion related.
!   MPL         Number of time  to print the outcoming.
!   t           The time.[d]
!   date        the time but in a date style.[yyyymmdd]
!   tEnd        the maxinum of time.[d]
!   TPrint(:)	the print time.[d]
!   interval    The Total Calculation Time (Unit Day).
! =========================related files==============================
!   selector.in
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Selector_In
    USE parm
    IMPLICIT NONE

    INTEGER (KIND=KI) :: i,j
    CHARACTER (LEN=100) :: Hed

    WRITE(*,*) 'Reading Basic information' 
    READ(33,*,err=901) 
    READ(33,*,err=901)
    READ(33,'(A100)',err=901) Hed
    WRITE(*,*) 
    WRITE(*,*) Hed
    READ(33,*,err=901) 
    READ(33,*,err=901) LUnit,TUnit,MUnit
    CALL Conversion(LUnit, TUnit, MUnit, xConv, tConv, mConv)
    READ(33,*,err=901)
    READ(33,*,err=901) lCheck,lWat,lChem,lTemp,lCrop,AtmBC
    READ(33,*,err=901)
    READ(33,*,err=901) Bup,Bdn,Drng,Dfit
    WRITE(*,*) 'Reading Material information'
    READ(33,*,err=902)
    READ(33,*,err=902)
    READ(33,*,err=902) NMat,NReg,NPar,CosAlf
    READ(33,*,err=902) 
    DO i=1,NMat
        READ(33,*,err=902) (Par(j,i),j=1,NPar),thF(i),thW(i),(sp(j,i),j=1,4)
        ths(i)=par(2,i)
        par(4,i) = par(4,i)*tConv/xConv
    ENDDO

    WRITE(*,*) 'Reading Time information'
    READ(33,*,err=903)
    READ(33,*,err=903)
    READ(33,*,err=903) dt,ddn,MPL!,MMPL
    dt = dt/tConv
    dtOld = dt
    READ(33,*,err=903)
    MaxAL = 0
    READ(33,*,err=903) date,tinit,tEnd,MaxAL
    t = tinit/tConv+dt
    tEnd = tEnd/tConv
    READ(33,*,err=903)
    READ(33,*,err=903) (TPrint(i),i=1,MPL)
!   READ(33,*)
!   READ(33,*) (TB(i),i=1,MMPL)
!    interval=int(tEnd-t+0.99_KR)! The total simulation period.
    Plevel = 1

!   The solute transport module.
    IF(lchem) THEN
        PAUSE
    ENDIF
!   The temperature module.
    IF(lTemp) THEN
        PAUSE
    ENDIF
    
    CALL Examine1
    CLOSE(33)
    RETURN
    
901 Terr=1
    RETURN
902 Terr=2
    RETURN
903 Terr=3
    RETURN
    
    END SUBROUTINE Selector_In

! ====================================================================
!   Subroutine Profile_In
!
!   Purpose: read information in the uz.in file for 1D model.
! ====================================================================                        
! =========================Incoming variables=========================
!   None
! =========================Outcoming variables========================
!   Nlayer      Number of nodes in every column.
!   zx          The z coordinates of every node descend from up to down.[m]
!   dz          The thickness of every layer.
!   MATuz       The material number.	
!   th(:)       Initial profile moisture.
!   Conc(:)     Initial solute concentration.
! =========================related files==============================
!   uz.in
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Profile_In
    USE parm
    CHARACTER (LEN=20) :: Text
    INTEGER (KIND=KI) :: i, j

    WRITE(*,*) 'Reading information for unsaturated zone'

!   the Nlayer
    READ(32,*,err=901)
    READ(32,*,err=901) Nlayer, NObs
    
    IF (.NOT. ALLOCATED(dz       )) ALLOCATE(dz(Nlayer       ))
    IF (.NOT. ALLOCATED(zx       )) ALLOCATE(zx(Nlayer+1     ))
    IF (.NOT. ALLOCATED(MATuz    )) ALLOCATE(MATuz(Nlayer    ))
    IF (.NOT. ALLOCATED(th       )) ALLOCATE(th(Nlayer       ))
    IF (.NOT. ALLOCATED(Sink1d   )) ALLOCATE(Sink1d(Nlayer   ))
    IF (.NOT. ALLOCATED(Epi      )) ALLOCATE(Epi(Nlayer      ))
    IF (.NOT. ALLOCATED(Tri      )) ALLOCATE(Tri(Nlayer      ))
    IF (.NOT. ALLOCATED(dnn      )) ALLOCATE(dnn(Nlayer      ))
    IF (.NOT. ALLOCATED(Slope    )) ALLOCATE(Slope(Nlayer    ))
    IF (.NOT. ALLOCATED(Intercept)) ALLOCATE(Intercept(Nlayer))
    IF (.NOT. ALLOCATED(par_n    )) ALLOCATE(par_n(Nlayer    ))
    
    zx = 0.0_KR
    dz = 0.0_KR
    MATuz = 0_KI
    th = 0.0_KR

!   the height.
    READ(32,*,err=901)
    READ(32,*,err=901) (zx(j),j=1,Nlayer+1,1) !bottom to surface!
    zx = zx/xConv
    
    DO j=1,Nlayer
        dz(j)=zx(j+1)-zx(j) ! The thickness of each layer
    ENDDO	 
!   the material kind.
    READ(32,*,err=901)
    READ(32,*,err=901) (MATuz(j),j=1,Nlayer,1)
    IF (NReg > 1) THEN
        READ(32,*,err=901)
        READ(32,*,err=901) (REGuz(j), j=1,Nlayer,1)
    ENDIF    
    IF (NObs > 0) THEN
        IF (.NOT. ALLOCATED(Obs)) ALLOCATE(Obs(NObs))
        READ(32,*,err=901)
        READ(32,*,err=901) (Obs(j),j=1,NObs,1)
        Text ='   Theta      '
        WRITE(80,110) (Obs(j),j=1,NObs,1)
        WRITE(80,120) (Text,j=1,NObs,1)
    ENDIF    
!   the initial profile moisture.
    READ(32,*,err=901)
    READ(32,*,err=901) (th(j),j=1,Nlayer,1)
    
    CALL Examine2
    
    CLOSE(32)
110 FORMAT (///4x,5(5x,'Node(',i3,')',5x))
120 FORMAT(/' time ',5(a20))
    RETURN
    
901 Terr=1
    RETURN

    END SUBROUTINE Profile_In

! ====================================================================
!   Subroutine Conversion   
!     
!   Purpose: Conversions unit
! ====================================================================
    SUBROUTINE Conversion(LUnit,TUnit,MUnit,xConv,tConv,mConv)
    USE parm, ONLY : KR, KI
    IMPLICIT NONE

    CHARACTER (len=5) :: LUnit
    CHARACTER (len=5) :: TUnit
    CHARACTER (len=5) :: MUnit
    REAL (kind=KR) :: xConv
    REAL (kind=KR) :: tConv
    REAL (kind=KR) :: mConv
    
    xConv = 1.0_KR
    tConv = 1.0_KR
    mConv = 1.0_KR

    IF (LUnit .eq. "km  ") THEN
        xConv = 0.001_KR
    ELSEIF (LUNIT .eq. "m  ") THEN
        xConv = 1.0_KR
    ELSEIF (LUnit .eq. "dm  ") THEN
        xConv = 10.0_KR
    ELSEIF (LUnit .eq. "cm  ") THEN
        xConv = 100.0_KR
    ELSEIF (LUnit .eq. "mm  ") THEN
        xConv = 1000.0_KR
    ELSE
        WRITE(*,*) 'LUnit is wrong!'
        WRITE(99,*) 'LUnit is wrong!'
        PAUSE
        STOP ! Stop the Program
    ENDIF

!   Base Unit of Time is Set as Day
    IF (TUnit .eq. "s  ") THEN
        tConv = 1.0_KR*60.0_KR*60.0_KR*24.0_KR
    ELSEIF (TUnit .eq. "min ") THEN
        tConv = 1.0_KR*60.0_KR*24.0_KR
    ELSEIF (TUnit .eq. "hours") THEN
        tConv = 1.0_KR*24.0_KR
    ELSEIF (TUnit .eq. "days") THEN
        tConv = 1.0_KR
    ELSEIF (TUnit .eq. "years") THEN
        tConv = 1.0_KR/365.0_KR
    ELSE
        WRITE(*,*) 'TUnit is wrong!'
        WRITE(99,*) 'TUnit is wrong!'
        PAUSE
        STOP ! Stop the Program
    ENDIF
    
    !IF (MUnit .eq. "g  ") THEN
    !    mConv = 1000.0_KR
    !ELSEIF (MUnit .eq. "kg  ") THEN
    !    mConv = 1.0_KR
    !ELSE
    !    WRITE(*,*) 'MUnit is wrong!'
    !    WRITE(99,*) 'MUnit is wrong!'
    !    PAUSE
    !    STOP ! Stop the Program
    !ENDIF

    RETURN
    END SUBROUTINE Conversion

! ====================================================================
!   Subroutine Examine
!     
!   Purpose: Examine the input Parameters
! ====================================================================
    SUBROUTINE Examine1
    USE parm 
    
    INTEGER (kind=4) :: Err
    Err = 0

!   Examine the Number of Hydraulic Parameters
    IF (NPar /= 5) THEN
        WRITE(*,*) 'The numbers of parameters maybe wrong!'
        WRITE(99,*) 'The numbers of parameters maybe wrong!'
        Err = Err + 1
    ENDIF
    
    DO i = 1,NMat
        IF (thF(i) > par(2,i) .or. thW(i) > par(2,i)) THEN
            WRITE(*,*) 'thF or thW is greater than ths!'
            WRITE(99,*) 'thF or thW is greater than ths!'
            Err = Err + 1
        ENDIF
        IF (thF(i) < par(1,i) .or. thW(i) < par(1,i)) THEN
            WRITE(*,*) 'thF or thW is smaller than thr!'
            WRITE(99,*) 'thF or thW is smaller than thr!'
            Err = Err + 1
        ENDIF
    ENDDO

    IF (Err) THEN
        PAUSE
        STOP ! Stop the Program
    ENDIF
    
    IF (dt > 1) THEN
        WRITE(*,*) 'dt > 1d, Caution the time of boundary condition and atmospheric condition'
    ENDIF

    END SUBROUTINE Examine1

    SUBROUTINE Examine2
    USE parm 

    INTEGER (kind=4) :: Err

    Err = 0

    IF (zx(1) /= 0) THEN
        WRITE(*,*) 'The first number of zx does not zero!'
        WRITE(99,*) 'The first number of zx does not zero!'
        Err = Err + 1
    ENDIF
    
    DO i = 1,NMat
        IF (dz(i) <= 0) THEN
            WRITE(*,*) 'The thickness of layer is negative!'
            WRITE(99,*) 'The thickness of layer is negative!'
            Err = Err + 1
        ENDIF
        IF (th(i) > par(2,i) .or. th(i) < par(1,i)) THEN
            WRITE(*,*) 'th is out of bounds!'
            WRITE(99,*) 'th is out of bounds!'
            Err = Err + 1
        ENDIF
    ENDDO

    IF (Err) THEN
        PAUSE
        STOP ! Stop the Program
    ENDIF

    END SUBROUTINE Examine2

! ===================================================================!
!   UBMOD  -- Water balance method of one-dimensional soil water     !
!             movement. Version 1.10.                                !
!                                                                    !
!   Designed by Wei Mao, Yan Zhu and Jinzhong Yang.                  !
!                                                                    !
!   Cited: Mao W, Yang J, Zhu Y, et al. An efficient soil water      !
!          balance model based on hybrid numerical and statistical   !
!          methods[J]. Journal of hydrology, 2018, 559: 721-735,     !
!          doi: 10.1016/j.jhydrol.2018.02.074.                       !
!                                                                    !
!   Feel free to contact us if you have any question.                !
!       Email: weimao@whu.edu.cn, zyan0701@163.com                   !
!                                                                    !
!                            Last modified: Sep, 2019, by Wei Mao.   !
! ===================================================================!
! ====================================================================
!     Input files:
!     1. Essential files.
!         SELECTOR.IN     The basic input information.
!         uz.IN           The discrete information.
!     2. Optional files.
!         cropdat.dat     The simple crop model.
!         01.wea          Meteorological data.
! ====================================================================
! ====================================================================
!     Output files:
!         
! ====================================================================
!   storage Routing Method~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   *************************end of documentation.	
    PROGRAM Main_Program
    USE parm
    IMPLICIT NONE
    REAL (KIND=KR) :: t1, t2, Zero
    INTEGER (KIND=4) :: lengthpath, ierr
    CHARACTER (100) :: filename, datapath

    Zero = 0.0_KR
    Terr = 0_KI
    
    filename = 'LEVEL_01.DIR'
    OPEN(10,file=filename, status='old',err=901)
    READ(10,'(a)',err=904) datapath
    CLOSE(10)
    lengthpath = Len_Trim(datapath)
!---Open the interface files. The relpath is used here.
! ====================================================================
!   input files
    OPEN(33,file=datapath(1:lengthpath)//'/Selector.IN', status='old',err=902)
    OPEN(32,file=datapath(1:lengthpath)//'/Profile.IN' , status='old',err=902)
!   output files
    OPEN(90,file=datapath(1:lengthpath)//'/runtime.OUT'  ,status='unknown',err=903) ! Run time.
    OPEN(80,file=datapath(1:lengthpath)//'/thObs.dat'    ,status='unknown',err=903) ! Node data.
    OPEN(89,file=datapath(1:lengthpath)//'/balance1d.dat',status='unknown',err=903) ! Statistic boundary condition.
    OPEN(81,file=datapath(1:lengthpath)//'/thprofile.dat',status='unknown',err=903) 
    OPEN(99,file=datapath(1:lengthpath)//'/error.txt'    ,status='unknown',err=903) ! Error message.
! ====================================================================
      
!-----Begin of the program.
! ====================================================================
!     subroutine about input information.
!     call for basic information. 
    CALL Selector_In
    IF (Terr.ne.0) GOTO (905) Terr
!     call for node information in 1D.
    CALL Profile_In
! ====================================================================

!-----preparation of the calculation
! ====================================================================
!     CPU time.
    CALL CPU_time (t1)
!     The initial water amount in model.
    CALL Balance_Initial
!     Diffusion model.
    CALL Diffusion_Model
!     Call for reference Evaportranspiration and division of E&T.
    CALL Upper_Boundary(datapath)
! ====================================================================

!-----Begin time loop.
100 CONTINUE

! ====================================================================
!   Set upper boundary condition.
    CALL Set_Input
  
! ====================================================================
!     Four main processes.
!     Firstly, divide infiltration and surface runoff.
!     There is no surface hydrological process by now.
          !if () then
          !    call SCS
          !elseif () then
          !    call Green-Ampt
          !else
          !    qair = Q_infiltration
          !endif

! ====================================================================
!     Secondly, advective movement driven by gravitational potential.
!       "Tipping-bucket" method.
    CALL Water_Redis
    IF (Terr.ne.0) GOTO (930) Terr

! ====================================================================
!     Thirdly, source/sink term.
!     open the Files that stored E&T and the rain, and the writen ETa.
    IF(bup >= Zero) CALL Water_SetET
    IF (Terr.ne.0) GOTO (931) Terr

! ====================================================================
!     Last, Diffusive soil water movement driven by matric potential.
    CALL Water_Diff
    IF (Terr.ne.0) GOTO (932) Terr

! ====================================================================
!     Output control.
!     Output the hydraulic head and soil moisture in 1D model.
    CALL Hthuz_out       
!   Call for the water balance in 1D and 3D model.
    CALL BalanceT   
!   call for new time and time step.
    WRITE(*,*)"t=",sngl(t)
    
!   P-Level information
    IF (abs(TPrint(Plevel)-t) < Tol) THEN
        CALL thOut
        Plevel = Plevel + 1
    ENDIF
    
    IF(abs(t-Tend) <= Tol) THEN
        GOTO 200
    ENDIF
    
    IF (dt < dtOld) THEN
        dt = dtOld
    ENDIF
    dtOld = dt
    CALL Tcontrol
    TLevel = TLevel + 1
    t = t+dt
    
    !IF (bup == 1) THEN
    !    WRITE(150,'(A18,F10.3,2F10.6)')date,t,Epa,Tra
    !    CALL DateAdd(date,1,date)
    !    Epa=Zero
    !    Tra=Zero
    !ENDIF

    GOTO 100
    
200 CALL CPU_time (t2)
    WRITE(90,*)'Real time [sec]',t2-t1
    CLOSE(90)
    CLOSE(80) 
    CLOSE(89) 
    CLOSE(81)
    CLOSE(88)
    CLOSE(98)
    CLOSE(99)
    CLOSE(110) 
    CLOSE(130)
    CLOSE(150) 

    STOP
    
901 ierr=1
    GOTO 999
902 ierr=2
    GOTO 999
903 ierr=3
    GOTO 999
904 ierr=4
    GOTO 999
905 ierr=5
    GOTO 999
930 ierr=30
    GOTO 999
931 ierr=31
    GOTO 999
932 ierr=32
    GOTO 999

999 CALL Error_Out(ierr)
    PAUSE
    STOP
    
    END PROGRAM Main_Program
    
    SUBROUTINE Error_Out(ierr)
    IMPLICIT NONE
    INTEGER (KIND=4) :: ierr
    CHARACTER (LEN=100), DIMENSION(40) :: cErr

    cErr( 1)='Open file error in file LEVEL_01.DIR !'
    cErr( 2)='Open file error for input files !'
    cErr( 3)='Open file error for output files !'
    cErr( 4)='Error when reading from an input file LEVEL_01.DIR !'
    cErr( 5)='Error when reading from an input file Selector.in !'
    cErr(30)='Mass balance error in Water_Redis module !'
    cErr(31)='Mass balance error in Water_SetET module !'
    cErr(32)='Mass balance error in Water_Diff module !'
    
    WRITE(*,*) cErr(ierr)
    WRITE(99,*) cErr(ierr)
    
    END SUBROUTINE Error_Out
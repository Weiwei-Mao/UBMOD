! ====================================================================
!     Subroutine module   
!     
!     Purpose: store the variables that will be used in the model 
!              frequently.
! ==================================================================== 
   
    MODULE parm
    IMPLICIT NONE

    INTEGER, PUBLIC, PARAMETER :: KR = SELECTED_REAL_KIND(r=20, p=20) !Accuracy of Real Number.
    INTEGER, PUBLIC, PARAMETER :: KI = SELECTED_INT_KIND(4)           !Accuracy of Int Number.
    INTEGER, PUBLIC, PARAMETER :: MMPLD=50          ! Print Related.
    INTEGER, PUBLIC, PARAMETER :: NMatD=8          ! Max numbers of soil materials (Vertical)
    INTEGER, PUBLIC, PARAMETER :: NREGD=8
    INTEGER, PUBLIC, PARAMETER :: NPARD=10        ! Max numbers of hydraulic parameters.
    INTEGER, PUBLIC, PARAMETER :: NLAYERD=130      ! Max numbers of soil layers.
    INTEGER, PUBLIC, PARAMETER :: numc=1           ! Max numbers of plant species.

    LOGICAL (KIND=4) :: lCheck
    LOGICAL (KIND=4) :: lWat
    LOGICAL (KIND=4) :: lChem   ! If to calculate the solute.
    LOGICAL (KIND=4) :: lTemp
    LOGICAL (KIND=4) :: lCrop
    LOGICAL (KIND=4) :: AtmBC       ! If to calculate the ET0 with P-M equation.

    INTEGER (KIND=KI) :: Terr
    INTEGER (KIND=KI) :: Nlayer     ! Numbers of real soil layers.
    INTEGER (KIND=KI) :: NObs       ! Numbers of the observation points.
    INTEGER (KIND=KI) :: Nmat       ! Numbers of soil materials.
    INTEGER (KIND=KI) :: NReg
    INTEGER (KIND=KI) :: MPL        ! Print related.
    INTEGER (KIND=KI) :: interval   ! The total calculation time (Unit day).
    INTEGER (KIND=KI) :: Bup        ! Upper boundary condition 0/1/2.
    INTEGER (KIND=KI) :: Bdn        ! Lower boundary condition 0/2.
    INTEGER (KIND=KI) :: Npar       ! Numbers of soil hydraulic parameters.
    INTEGER (KIND=KI) :: Drng       ! Which Drainage Function 1/2/3/4/5/6/7.
    INTEGER (kind=KI) :: Dfit       ! The Empirical Formula of Diffusion -3/-2/-1/1/2/3.
    INTEGER (KIND=KI) :: ddn        ! Diffision process related.
    INTEGER (KIND=KI) :: Nup        ! Numbers of flux upper boundary node.
    INTEGER (KIND=KI) :: Ndn        ! Numbers of flux lower boundary node.
    INTEGER (KIND=KI) :: MaxAL      ! Numbers of atmospheric data-records.
    INTEGER (KIND=KI) :: Plevel
    INTEGER (KIND=KI) :: Tlevel
    INTEGER (KIND=KI) :: Alevel
    INTEGER (KIND=KI), ALLOCATABLE :: Obs(:)    ! The observation points.
    INTEGER (KIND=KI), ALLOCATABLE :: MATuz(:)  ! The Material Serial Number of Each Layer.
    INTEGER (KIND=KI), ALLOCATABLE :: REGuz(:)

    CHARACTER (LEN=8) :: date
    CHARACTER (LEN=5) :: LUnit      ! Length Unit. 
    CHARACTER (LEN=5) :: TUnit      ! Time Unit.
    CHARACTER (LEN=5) :: MUnit      ! Mass Unit.
        
    REAL (KIND=KR) :: Tol=1E-10_KR
    REAL (KIND=KR) :: CosAlf
    REAL (KIND=KR) :: dt            ! Time step.
    REAL (KIND=KR) :: dtOld
    REAL (KIND=KR) :: TotalPerco    ! The Percolation flux of the soil column.
    REAL (KIND=KR) :: t             ! Time.
    REAL (KIND=KR) :: tinit         ! The initial time.
    REAL (KIND=KR) :: tEnd          ! The endding time.
    REAL (KIND=KR) :: tAtm
    REAL (KIND=KR) :: qair          ! Precipitation.
    REAL (KIND=KR) :: qqair 
    REAL (KIND=KR) :: qairt         ! Sum of Precipitation.
    REAL (KIND=KR) :: qbtmt         ! Sum of Bottom Flux.
    REAL (KIND=KR) :: CumE          ! Sum of Evaporation.
    REAL (KIND=KR) :: CumT          ! Sum of Transpiration.
    REAL (KIND=KR) :: sinkt 
    REAL (KIND=KR) :: Qbtm
    REAL (KIND=KR) :: voluz         ! The initial water volumn.
    REAL (KIND=KR) :: voluz1        !
    REAL (KIND=KR) :: Epa	          ! Actual Evaporation.
    REAL (KIND=KR) :: Tra	          ! Actual Transpiration.
    REAL (KIND=KR) :: ptab
    REAL (KIND=KR) :: xConv         ! Length conversion coefficient.
    REAL (KIND=KR) :: tConv         ! Time conversion coefficient.
    REAL (KIND=KR) :: mConv         ! Mass conversion coefficient.
    
    REAL (KINd=KR), ALLOCATABLE :: evatra(:,:),precip(:,:)

    REAL (KIND=KR), DIMENSION(MMPLD) :: TPrint        ! The print time.
    REAL (KIND=KR), DIMENSION(NMATD) :: ths           ! ths for each materials.
    REAL (KIND=KR), DIMENSION(NMATD) :: thF           ! Field capacity.
    REAL (KIND=KR), DIMENSION(NMATD) :: thw           ! Wilting point.
    REAL (KIND=KR), DIMENSION(2,MMPLD) :: up          ! Met.in file, Upper input.
    REAL (KIND=KR), DIMENSION(2,MMPLD) :: dn          ! Met.in file, Lower input.
    REAL (KIND=KR), DIMENSION(4,NMATD) :: sp          ! The Evaporation cumulative distribution function.
    REAL (KIND=KR), DIMENSION(NPARD,NMATD) :: Par     ! The Hydraulic parameters.

!   1D  model.
    REAL (KIND=KR), ALLOCATABLE :: dz(:)     ! Thickness of Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: zx(:)     ! The z Coordinates of Every Node Descend from Up to Down.
    REAL (KIND=KR), ALLOCATABLE :: th(:)     ! Soil Water Content.
    REAL (KIND=KR), ALLOCATABLE :: SInk1d(:) ! The Actural Evapotranspiration for Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: Epi(:)    ! Evaporation for Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: Tri(:)    ! Transpiration for Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: dnn(:)
    REAL (KIND=KR), ALLOCATABLE :: Slope(:)      ! Empirical Formula of Diffusion.
    REAL (KIND=KR), ALLOCATABLE :: Intercept(:)  ! Empirical Formula of Diffusion.
    REAL (KIND=KR), ALLOCATABLE :: par_n(:)      ! Empirical Formula of n.
    END MODULE parm

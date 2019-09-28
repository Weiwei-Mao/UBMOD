! ====================================================================
!   SUBROUTINE BalanceI   
!     
!   PROPOSE: the initial statistics.
! ====================================================================
! =========================Incoming variables=========================
!   par(:,:)    The hydraulic parameters.
!   th(:,:)     The soil water content.
!   dz(:)       The thickness of each layer.
!   matuz(:)    The material serial number.
! =========================Outcoming variables========================
!   voluz       The initial water volumn [m3].
! =========================related files==============================
!   balance1d.dat
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Balance_Initial
    USE parm
    IMPLICIT NONE
    
    voluz = 0.0_KR
    qairt = 0.0_KR
    qbtmt = 0.0_KR
    CumE  = 0.0_KR
    CumT  = 0.0_KR
    sinkt = 0.0_KR
    sink1d = 0.0_KR

    voluz=sum(th(1:Nlayer)*dz(1:Nlayer))  !!!The specific storage is ignored.
    WRITE(89,*)"Variables=t,voluz,Dvoluz,qair,qbtm,CumE,CumT,Error,Error%"

    RETURN
    END SUBROUTINE Balance_Initial

! ====================================================================
!   SUBROUTINE Diffusion_Model   
!     
!   PROPOSE: Diffusion model, by MaoWei
! ====================================================================
! =========================Incoming variables=========================
!   Ks          Saturated soil hydrualic conductivity.
! =========================Outcoming variables========================
!   Slope       Slope.
!   Intercept   Intercept.
! =========================related files==============================
!   None.
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Diffusion_Model
    USE parm
    IMPLICIT NONE
    
    INTEGER (KIND=KI) :: m, i
    
    DO i = 1,Nlayer
        m = MATuz(i)
        ! Curve fitting.
        IF (Dfit == 1) THEN ! quadratic curve, Dfit>0, theta
            Intercept(m) = -2.9602*log10(par(5,m))**2-0.6179*log10(par(5,m))-4.2946
            Slope(m) = 9.4519*log10(par(5,m))**2+6.6871*log10(par(5,m))+11.115
        ELSEIF (Dfit == 2) THEN ! first segmentation fitting
            IF (log10(par(5,m))>=-0.99866) THEN
                Intercept(m) = -0.6245*log10(par(5,m))**2-0.4761*log10(par(5,m))-5.2872
            ELSE
                Intercept(m) = 9.8771*log10(par(5,m))+4.4293
            ENDIF
            IF (log10(par(5,m))>=-0.99866) THEN
                Slope(m) = 1.6683*log10(par(5,m))**2+6.0643*log10(par(5,m))+14.49
            ELSE
                Slope(m) = -27.225*log10(par(5,m))-17.452
            ENDIF                
        ELSEIF (Dfit == 3) THEN ! second segmentation fitting
            IF (log10(par(5,m))>=-0.9896858) THEN
                Intercept(m) = 0.8771*log10(par(5,m))**2+1.0098*log10(par(5,m))-5.0988
            ELSE
                Intercept(m) = 9.8771*log10(par(5,m))+4.4293
            ENDIF
            IF (log10(par(5,m))>=-0.9990038) THEN
                Slope(m) = -3.9735*log10(par(5,m))**2+0.0609*log10(par(5,m))+13.489
            ELSE
                Slope(m) = -27.225*log10(par(5,m))-17.452
            ENDIF
        ELSEIF (Dfit == -1) THEN ! quadratic curve, Dfit<0, S
            Slope(m) = 2.8657*log10(par(5,m))**2+2.5014*log10(par(5,m))+4.2306
            Intercept(m) = -2.4228*log10(par(5,m))**2-0.4907*log10(par(5,m))-3.5813
        ELSEIF (Dfit == -2) THEN ! segmentation fitting, without amendment, left 6, right 6
            IF (log10(par(5,m))>=-1.004802) THEN
                Intercept(m) = -0.9483*log10(par(5,m))**2-0.4217*log10(par(5,m))-4.1855
            ELSE
                Intercept(m) = 7.8942*log10(par(5,m))+3.2129
            ENDIF
            IF (log10(par(5,m))>=-1.016789) THEN
                Slope(m) = 1.6523*log10(par(5,m))**2+2.5041*log10(par(5,m))+4.961
            ELSE
                Slope(m) = -7.2474*log10(par(5,m))-3.511
            ENDIF
        ELSEIF (Dfit == -3) THEN
            IF (log10(par(5,m))>=-0.8989257) THEN
                Intercept(m) = 0.12*log10(par(5,m))**2+0.511*log10(par(5,m))-4.1381
            ELSE
                Intercept(m) = 5.2687*log10(par(5,m))+0.2358
            ENDIF
            IF (log10(par(5,m))>=-0.9171914) THEN
                Slope(m) = -0.4014*log10(par(5,m))**2+0.5861*log10(par(5,m))+4.518
            ELSE
                Slope(m) = -4.8045*log10(par(5,m))-0.7639
            ENDIF
        ENDIF       

        ! par_n(m) = 15.812*thf(m)**2-10.416*thf(m)+2.8891
        ! par_n(m) = 0.1671*log10(par(5,m))**2+0.7103*log10(par(5,m))+1.8925
        par_n(m) = 0.9509*exp(0.72223*log10(par(5,m)))+0.9

    ENDDO
    
    END SUBROUTINE Diffusion_Model

! ====================================================================
!   SUBROUTINE Boundary  
!     
!   PROPOSE: prepare boundary condition
! ====================================================================
! =========================Incoming variables=========================
!   Bup          Up Boundary Condition -1/0/1/2.
!   Bdn          Down Boundary Condition 0/2.
! =========================Outcoming variables========================
!   None.
! =========================related files==============================
!   None.
! =========================related functions==========================
!   None.
! ==================================================================== 
    SUBROUTINE Upper_Boundary(datapath)
    USE parm
    CHARACTER (LEN=8) :: DDate
    CHARACTER (LEN=100) :: datapath
    INTEGER (KIND=4) :: lenpath
    INTEGER (KIND=KI) :: Jd0

    lenpath = Len_Trim(datapath)
    IF (bup == 1) THEN  
        CALL Etp(1,datapath)
        OPEN(110,file=datapath(1:lenpath)//'/'//'01.et0',status='old')
        READ(110,*)
        OPEN(130,file=datapath(1:lenpath)//'/'//'01.eti',status='old')
        OPEN(150,file=datapath(1:lenpath)//'/'//'eta.dat',status='unknown')
        WRITE(150,*)"Variables=DoY,   Ea,   Ta,   Date"
        READ(130,*)
                
        IF (.NOT. ALLOCATED(precip)) ALLOCATE(precip(2,MaxAL))
        IF (.NOT. ALLOCATED(Evatra)) ALLOCATE(Evatra(2*Nlayer,MaxAL))
        
        DO i = 1,MaxAL
            READ(110,*) Nouse,precip(1,i),Nouse,Nouse,Nouse,precip(2,i)
            READ(130,*) Nouse,Nouse,(Evatra(j,i),j=1,2*Nlayer)
        ENDDO
        precip(2,:) = precip(2,:)/1000_KI
        Evatra = Evatra/1000_KI
        
    ELSEIF (bup == 2) THEN
        OPEN(120,file='Rh1D.in/'//trim(iof)//'/'//'met.in',status='old')
        READ(120,*)
        READ(120,*)Nup,Ndn
        READ(120,*)
        READ(120,*)(up(1,i),i=1,Nup)
        READ(120,*)
        READ(120,*)(up(2,i),i=1,Nup)
        READ(120,*)
        READ(120,*)(dn(1,i),i=1,Ndn)
        READ(120,*)
        READ(120,*)(dn(2,i),i=1,Ndn)
    ENDIF
   
    END SUBROUTINE Upper_Boundary
    
! ====================================================================
!   Subroutine Set_Input
!     
!   Purpose: Set the boundary input values
! ====================================================================
    SUBROUTINE SetQ
    USE parm
    INTEGER (kind=KI) :: k, kk
    
    IF (bup == 1) THEN

        CALL FindY_Step(MaxAL,t,qair,precip,k)
        Epi(1:Nlayer) = Evatra(1:NLayer,k)
        Tri(1:Nlayer) = Evatra(NLayer+1:2*Nlayer,k)
        tAtm = k

    ELSEIF (bup == 2) THEN
        CALL FindY_Step(Nup,t,qair,up,k)
    ENDIF
   
    END SUBROUTINE SetQ

! ====================================================================
!   Subroutine FindY_Step 
!     
!   Purpose: interpolation
! ====================================================================
! =========================Incoming variables=========================
!   n           The total numbers of the node.
!   x           The input value.
!   xy(2,:)     The first row is the value corresponding to x,
!               the second row is the value corresponding to y.
! =========================Outcoming variables========================
!   y           The output value.
! ====================================================================
    SUBROUTINE FindY_Step(n,x,y,xy,flag)
    USE parm

    INTEGER (kind=KI) :: n,flag
    REAL (kind=KR) :: x
    REAL (kind=KR) :: y
    REAL (kind=KR), DIMENSION(2,n) :: xy
    
    IF (x <= xy(1,1)) THEN
        Y=xy(2,1)
        flag = 1
        RETURN
    ENDIF
    IF (x >= xy(1,n)) THEN
        y=xy(2,n)
        flag = n
        RETURN
    ENDIF
    DO i=1,n-1
        IF (x > xy(1,i) .and. x <= (xy(1,i+1)+Tol)) THEN
            y=xy(2,i+1)
            flag = i+1
            RETURN
        ENDIF   
    ENDDO

    END SUBROUTINE FindY_Step
    

! ====================================================================
!     subroutine SCS   
!     
!     purpose: division preciptation and runoff 
! ====================================================================
      !subroutine SCS
      !    use parm
      !    
      !    SCS_S = 254*(100/CN-1)
      !    SCS_R0 = (P-0.2*SCS_S)**2/(P+0.8*SCS_S)
      !end

! ====================================================================
!     subroutine Green-Ampt   
!     
!     purpose: division preciptation and runoff 
! ====================================================================
      !subroutine Green-Ampt
      !    use parm
      !    
      !    SCS_S = 254*(100/CN-1)
      !    SCS_R0 = (P-0.2*SCS_S)**2/(P+0.8*SCS_S)
      !end


! ====================================================================
!     Subroutine infil   
!     
!     Purpose: The allocation of the infiltration water.
! ====================================================================
    SUBROUTINE infil
    USE parm
    IMPLICIT NONE
    
    REAL (kind=KR) :: step, qairtemp, thss
    REAL (kind=KR) :: qvoluz      ! The Difference of Saturated and Initial Water Volumn.
    INTEGER (kind=KI) :: i, j, m
          
    qairtemp = qair*dt
    DO j=1,Nlayer
        m=MATuz(j)
        thss=par(2,m)

        IF (qairtemp>(thss-th(j))*dz(j)) THEN	
            qairtemp=qairtemp-(thss-th(j))*dz(j)
            th(j)=thss
        ELSE
            th(j)=th(j)+qairtemp/dz(j)
            qairtemp=0.0_KR
        ENDIF
    ENDDO
    
    END SUBROUTINE infil

! ====================================================================
!     Subroutine Water_Redis   
!     
!     Purpose: calculate the redistribution of the water
! ====================================================================      
    SUBROUTINE Water_Redis 
    USE parm
    IMPLICIT NONE
    
    INTEGER (kind=KI) :: m, i, j
    REAL (kind=KR) :: step,perco,excess,temp,thnew1,sum0,sum1,beta1,delta
    REAL (kind=KR) :: alpha,belta1,A,n,Se,unsatM,unsatK

    step=dt
    TotalPerco=0.0_KR

!    drainage
! ==================================================================== 
         ! do   
    perco=0.0_KR
    temp=0.0_KR
    sum0 = 0.0_KR
    sum1 = 0.0_KR
              
    CALL infil  !allocation of infiltration water
    
    DO i = 1,Nlayer
        sum0 = sum0+th(i)*dz(i)
    ENDDO
              
    DO j=1,Nlayer
        m=MATuz(j)
        th(j)=th(j)+temp

        IF (Drng==1) THEN !SWAT equation linear
            IF (th(j)>thf(m)) THEN
                thnew1=thF(m)+(th(j)-thF(m))*EXP(-par(5,m)*step/dz(j)/(par(2,m)-thf(m)))
            ELSE
                thnew1=th(j)
            ENDIF
        ELSEIF (Drng==2) THEN !exp(-alfa*theta);Kendy model; par(6) 10-30
            IF (th(j)>thf(m)) THEN
                alpha=par(3,m)
                IF (alpha >= 30 .and. alpha <= 10) THEN
                    WRITE(*,*)"Please check alpha"
                    STOP
                ENDIF
                thnew1=ths(M)-(ths(M)-thw(M))/alpha*log(alpha*par(5,m)*step/dz(j)/(ths(M)-thw(M)) &
                      +exp(alpha*(ths(M)-th(j))/(ths(M)-thw(M))))
            ELSE
                thnew1=th(j)
            ENDIF
        ELSEIF (Drng==3) THEN !(theta)beta; Gardner model; beta K=0.002 ks
            IF (th(j)>thf(m)) THEN
                Beta1=-2.655/log10(thf(m)/ths(m))
                thnew1=(th(j)**(1-beta1)+(beta1-1)*par(5,m)*step/dz(j)/ths(M)**beta1)**(1/(1-beta1))
            ELSE
                thnew1=th(j)
            ENDIF
        ELSEIF (Drng==4) THEN !
            IF (th(j)>thf(m)) THEN
                A = -par(5,m)/dz(j)/(exp(par(2,m)-thf(m))-1)
                thnew1=-log(1.0/(1-exp(A*step)*(1.0-1.0/exp(th(j)-thf(m)))))+thf(m)
            ELSE
                thnew1=th(j)
            ENDIF
        ELSEIF (Drng==5) THEN !
            IF (th(j)>thf(m)) THEN
                thnew1=thF(m)+dz(j)*(par(2,m)-thF(m))**2/(par(5,m)*step+dz(j)*(par(2,m)-thF(m))**2/(th(j)-thF(m)))
            ELSE
                thnew1=th(j)
            ENDIF
        ELSEIF (Drng==6) THEN  ! 
            IF (th(j)>thf(m)) THEN
                n=par(4,m)
                thnew1=par(1,m)+(par(5,m)*step*(2/n+2)/dz(j)/((par(2,m)-par(1,m))**(2/n+3))+(th(j)-  &
                           par(1,m))**(-2/n-2))**(1/(-2/n-2))
            ELSE
                thnew1=th(j)
            ENDIF
        ELSEIF (Drng==7) THEN !
            IF (th(j)>thf(m)) THEN
                Se=MIN(1.0, (th(j)-par(1,m))/(par(2,m)-par(1,m)))
                unsatM=1.-1/1.56
                unsatK=par(5,m)*Se**0.5*(1- (1-se**(1/unsatM))**unsatM)**2
                thnew1=thF(M)+(th(j)-thF(M))*EXP(-unsatK*step/dz(j)/(par(2,m)-thf(m)))
            ELSE
                thnew1=th(j)
            ENDIF
        ENDIF
        perco=-(thnew1-th(j))*dz(j) !
        IF (j<Nlayer) THEN
            temp=perco/dz(j+1) !
        ELSEIF (j==Nlayer) THEN
            perco=-(thnew1-th(j))*dz(j)
        ENDIF
        th(j)=thnew1   !
    ENDDO  !

    j=Nlayer
    m=MATuz(j)          
    IF (bdn==0) THEN ! 
!           
        th(j)=th(j)+perco/dz(j)
        perco=0.
    ELSEIF (bdn==2) THEN !
        perco=perco
    ENDIF
    
    TotalPerco=TotalPerco+perco
!
    excess=0.0_KR !

    DO j=Nlayer,1,-1
        m=MATuz(j)
        th(j)=th(j)+excess/dz(j)
        
        IF (th(j)>par(2,m)) THEN	
            excess=(th(j)-par(2,m))*dz(j)
            th(j)=par(2,m)
        ELSE !
            excess=0.0_KR
        ENDIF
    
    ENDDO
    
    DO i = 1,Nlayer
        sum1 = sum1+th(i)*dz(i)
    ENDDO
    delta = sum0-sum1-TotalPerco
    IF (delta > Tol) THEN
        Terr = 1
        RETURN
    ENDIF

    !WRITE(99,"(130F10.7)")th

    END SUBROUTINE Water_Redis
    
! ====================================================================
!   Subroutine Water_SetET
! ====================================================================
!   MATuz(:)    material number
!   thw(:)      wilting point
!   thf(:)      field capacity
! ====================================================================
!   tra(:)      the actual crop transpiration [mm/d]
!   epa(:)      the actual soil evaporation [mm/d]
!   sink1d(:)   the total source/sink term [m/d]
! ====================================================================
    SUBROUTINE Water_SetET
    USE parm
    IMPLICIT NONE
    
    INTEGER (kind=KI) :: m,m1,i,j
    REAL (kind=KR) :: step,p
    REAL (kind=KR) :: sum0,sum1,delta
    REAL (kind=KR) :: thdry
    REAL (kind=KR) :: fEact,fTact
    REAL (kind=KR) :: Eact,Tact
    REAL (kind=KR) :: thcr
    REAL (kind=KR) :: Epy, Tpy

!
! re zero the Tra and epa.
    Tra=0.0_KR
    Epa=0.0_KR
    p=0.0_KR
    sum0 = 0.0_KR
    sum1 = 0.0_KR
!    
    DO i = 1,Nlayer
        sum0 = sum0+th(i)*dz(i)
    ENDDO
    DO j=1,Nlayer
        m=MATuz(j)
        sink1d(j)=0.
        thdry=(Thw(M)+par(1,M))/2. !
        fEact=MAX((Th(j)-Thdry)/(ThF(M)-Thdry),0D0)
!        fEact = 1.0
        Eact=fEact*Epi(j)*dt
        Epa=Epa+Eact !m/d
        p=ptab+0.04*(5-Tri(j))
        thcr=Thw(M)+(1-p)*(ThF(M)-Thw(M)) 
        fTact=MAX(0D0,MIN(1D0,(Th(j)-Thw(M))/(thcr-Thw(M))))
        Tact=fTact*Tri(j)*dt
        Tra=Tra+Tact !m/d
        Sink1d(j)=Tact+Eact !m/d
        !write(99,"(132F10.7)")fEact,fTact,Sink1d
    ENDDO

    DO j=1,Nlayer
        m=MATuz(j)
        Th(j)=Th(j)-sink1d(j)/dz(j)
        IF (Th(j) < par(1,m)) THEN
            IF (j == 1) THEN
                Th(j+1) = Th(j+1)+Th(j)-par(1,m)
                Th(j) = par(1,m)
            ELSEIF (j > 1) THEN
                m1 = MATuz(j-1)
                IF (Th(j-1)+Th(j) <= par(1,m)+par(1,m1)) THEN
                    Th(j+1) = Th(j+1)+Th(j)+Th(j-1)-par(1,m)-par(1,m1)
                    Th(j-1) = par(1,m1)
                    Th(j) = par(1,m)
                ELSEIF (Th(j-1)+Th(j) > par(1,m)+par(1,m1)) THEN
                    Th(j-1) = Th(j-1)+Th(j)-par(1,m)
                    Th(j) = par(1,m)
                ENDIF
            ENDIF
        ENDIF
    ENDDO
    
    DO i = 1,Nlayer
        sum1 = sum1+th(i)*dz(i)
    ENDDO
    delta = sum0-sum1-sum(Sink1d)
    IF (delta > Tol) THEN
        Terr = 1
        RETURN
    ENDIF
   
    !WRITE(99,"(130F10.7)")th

    END SUBROUTINE Water_SetET
      
! ====================================================================
!     Subroutine Water_Diff
!     
!     Purpose: calculate the diffusion process
! ====================================================================  
    SUBROUTINE Water_Diff
    USE parm
    IMPLICIT NONE
    
    INTEGER (kind=KI) :: m,m1,m11,NNN,i,j,k
    REAL (kind=KR) :: deltat,step
    REAL (kind=KR) :: Da,Db
    REAL (kind=KR) :: sum0,sum1,delta
    REAL (kind=KR) :: correction1,correction2
    REAL (kind=KR) :: mq,qku
    REAL (kind=KR) :: v
    REAL (kind=KR), DIMENSION(Nlayer) :: th0,th1,thu,Courant
    REAL (kind=KR), DIMENSION(Nlayer) :: D,S,dern
    REAL (kind=KR), DIMENSION(Nlayer,Nlayer) :: A
    REAL (kind=KR), DIMENSION(Nlayer) :: B
    REAL (kind=KR), DIMENSION(Nlayer) :: dq
    REAL (kind=KR), DIMENSION(Nlayer,2) :: qchange
    
    deltat = dt/ddn
    
    DO j = 1,ddn
        th0 = th  ! 

        DO i = 1,Nlayer
            m = MATuz(i)
            D(i) = 10**(Slope(m)*th0(i)+Intercept(m))
            S(i) = (th0(i)-par(1,m))/(par(2,m)-par(1,m))
            IF (abs(1.0 - S(i)) <= 1e-7) THEN
                S(i) = 1
                dern(i) = 0
            ELSEIF (abs(S(i)) <= 1e-7) THEN
                S(i) = 0
                dern(i) = 0
            ELSE
                dern(i) = (par(2,m)-par(1,m))*S(i)*(log(S(i))/(par_n(m)-1)*par_n(m)**3-(par_n(m)-1)/par_n(m)**2.0*(S(i)**(par_n(m)/(1-par_n(m)))-1)*log(S(i)**(par_n(m)/(1-par_n(m)))-1)*S(i)**(-par_n(m)/(1-par_n(m))))
            ENDIF
        ENDDO

    
        DO i = 2,Nlayer-1
            m   = MATuz(i)
            m1  = MATuz(i-1)
            m11 = MATuz(i+1)
            A(i,i) = dz(i)/deltat+2*sqrt(D(i-1)*D(i))/(dz(i-1)+dz(i))+2*sqrt(D(i+1)*D(i))/(dz(i+1)+dz(i))
            A(i,i+1) = -2*sqrt(D(i+1)*D(i))/(dz(i+1)+dz(i))
            A(i,i-1) = -2*sqrt(D(i-1)*D(i))/(dz(i-1)+dz(i))
            ! correction£¬ 
            correction1 = (S(i)+S(i-1))/2*(par(2,m1)-par(2,m))+(2.0-S(i)-S(i-1))/2*(par(1,m1)-par(1,m))+(dern(i)+dern(i-1))/2*(par_n(m1)-par_n(m))
            correction2 = (S(i)+S(i+1))/2*(par(2,m11)-par(2,m))+(2.0-S(i)-S(i+1))/2*(par(1,m11)-par(1,m))+(dern(i)+dern(i+1))/2*(par_n(m11)-par_n(m))            
            !WRITE(99,*)(S(i)+S(i-1))/2*(par(2,m1)-par(2,m))/correction1,(2.0-S(i)-S(i-1))/2*(par(1,m1)-par(1,m))/correction1,(dern(i)+dern(i-1))/2*(par_n(m1)-par_n(m))/correction1
            !WRITE(99,*)(S(i)+S(i+1))/2*(par(2,m11)-par(2,m))/correction2,(2.0-S(i)-S(i+1))/2*(par(1,m11)-par(1,m))/correction2,(dern(i)+dern(i+1))/2*(par_n(m11)-par_n(m))/correction2
            correction1 = correction1
            correction2 = correction2
            B(i) = th0(i)*dz(i)/deltat-correction1*2*sqrt(D(i-1)*D(i))/(dz(i)+dz(i-1))-correction2*2*sqrt(D(i+1)*D(i))/(dz(i)+dz(i+1))
        END DO

        A(1,1) = dz(1)/deltat+2*sqrt(D(2)*D(1))/(dz(2)+dz(1))
        A(1,2) = -2*sqrt(D(2)*D(1))/(dz(2)+dz(1))
        m  = MATuz(1)
        m1 = MATuz(2)
        correction1 = (S(1)+S(2))/2*(par(2,m1)-par(2,m))+(2-S(1)-S(2))/2*(par(1,m1)-par(1,m))+(dern(1)+dern(2))/2*(par_n(m1)-par_n(m))
        correction1 = correction1
        B(1) = th0(1)*dz(1)/deltat - correction1*2*sqrt(D(2)*D(1))/(dz(1)+dz(2))
        A(Nlayer,Nlayer) = dz(Nlayer)/deltat+2*sqrt(D(Nlayer-1)*D(Nlayer))/(dz(Nlayer-1)+dz(Nlayer))
        A(Nlayer,Nlayer-1) = -2*sqrt(D(Nlayer-1)*D(Nlayer))/(dz(Nlayer-1)+dz(Nlayer))
        m = MATuz(Nlayer)
        m1 = MATuz(Nlayer-1)
        correction1 = (S(Nlayer)+S(Nlayer-1))/2*(par(2,m1)-par(2,m))+(2-S(Nlayer)-S(Nlayer-1))/2*(par(1,m1)-par(1,m))+(dern(Nlayer)+dern(Nlayer-1))/2*(par_n(m1)-par_n(m))
        correction1 = correction1
        B(Nlayer) = th0(Nlayer)*dz(Nlayer)/deltat - correction1*2*sqrt(D(Nlayer-1)*D(Nlayer))/(dz(Nlayer-1)+dz(Nlayer))
        
        !WRITE(99,"(130F10.7)") A,B,th0
        
        CALL chase(A,B,th1,Nlayer)
        
        sum0 = 0.0_KR
        sum1 = 0.0_KR
        DO i = 1,Nlayer
            sum0 = sum0+th0(i)*dz(i)
            sum1 = sum1+th1(i)*dz(i)
        ENDDO
        delta = sum0 - sum1
        IF (delta > Tol) THEN
            Terr = 1
            RETURN
        ENDIF

        mq = 0.0_KR
        DO i = 1,Nlayer
            m = MATuz(i)
            IF (th1(i) > par(2,m)) THEN
                IF (i < Nlayer) THEN
                    mq = th1(i)-par(2,m)
                    th1(i) = par(2,m)
                    th1(i+1) = th1(i+1)+mq
                ELSEIF (i == Nlayer) THEN
                    DO k=Nlayer,1,-1
                        m=MATuz(k)
                        IF (k < Nlayer) THEN
                            th1(k)=th1(k)+mq
                        ENDIF
                        IF (th1(k)>par(2,m)) THEN
                            mq=th1(k)-par(2,m)
                            th1(k)=par(2,m)
                        ELSE 
                            mq=0.0_KR
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
        
        !mq = 0.
        !DO i = Nlayer,1,-1
        !    m = MATuz(i)
        !    IF (th1(i) < par(1,m)) THEN
        !        IF (i > 1) THEN
        !            mq = par(1,m) - th1(i)
        !            th1(i) = par(1,m)
        !            th1(i-1) = th1(i-1) - mq
        !        ELSEIF (i == 1) THEN
        !            DO k=1,Nlayer
              !          m=MATuz(k)
              !          th(k)=th(k)-mq
              !          IF (th(k)<par(1,m)) THEN	
                 !           mq=par(1,m)-th(k)
                 !           th(k)=par(1,m)
        !                ELSE 
                 !           mq=0.
        !                ENDIF
        !            ENDDO
        !        ENDIF
        !    ENDIF
        !ENDDO
        DO i = 1,Nlayer
            m = MATuz(i)
            IF (th1(i) < par(1,m)) THEN
                IF (i==1) THEN
                    th1(i+1) = th1(i+1)+th1(i)-par(1,m)
                    th1(i) = par(1,m)
                ELSEIF (i >1) THEN
                    m1 = MATuz(i-1)
                    IF (th1(i-1)+th1(i) <= par(1,m)+par(1,m1)) THEN
                        th1(i+1) = Th1(i+1)+Th1(i)+Th1(i-1)-par(1,m)-par(1,m1)
                        Th1(i-1) = par(1,m1)
                        Th1(i) = par(1,m)
                    ELSEIF (Th1(i-1)+Th1(i) > par(1,m)+par(1,m1)) THEN
                        Th1(i-1) = Th1(i-1)+Th1(i)-par(1,m)
                        Th1(i) = par(1,m)
                    ENDIF
                ENDIF
            ENDIF
        ENDDO    
        th = th1
        !WRITE(99,"(f7.1,130F10.7,2F10.7)")j,th,correction1,correction2
    ENDDO
    
    !WRITE(99,"(130F10.7)")th
    
    END SUBROUTINE Water_Diff

! ====================================================================
!   Subroutine Chase
!
!   Purpose: Solve the tridiagonal matrix Ax=f.
! ====================================================================                        
! =========================Incoming variables=========================
!   A        The tridiagonal.
!   B        The right item.
!   N        The dimension of the matrix.
! =========================Outcoming variables========================
!   x           The results.
! ====================================================================
    SUBROUTINE chase(A,f,x,N)
    USE parm, ONLY : KI, KR, Nlayer
    
    INTEGER (kind=KI) :: N
    REAL (kind=KR) :: A(Nlayer,Nlayer), f(Nlayer), x(Nlayer)
    REAL (kind=KR), DIMENSION(Nlayer) :: u,b
    REAL (kind=KR), DIMENSION(Nlayer) :: y
    REAL (kind=KR), DIMENSION(Nlayer) :: L,d,c,e
    
    u = 0.0_KR
    b = 0.0_KR
    y = 0.0_KR
    L = 0.0_KR
    d = 0.0_KR
    c = 0.0_KR
    e = 0.0_KR
    
    DO i = 1,N
        b(i) = A(i,i)
    ENDDO
    DO i = 1,N-1
        c(i) = A(i,i+1)
    ENDDO
    DO i = 2,N
        e(i) = a(i,i-1)
    ENDDO
    
    DO i = 1,N-1
        d(i)=c(i)
    ENDDO
    u(1) = b(1)
    DO i = 2,N
        L(i)=e(i)/u(i-1)
        u(i)=b(i)-L(i)*c(i-1)
    ENDDO
!------y

    y(1)=f(1)
    DO i=2,N
        y(i)=f(i)-L(i)*y(i-1)
    ENDDO

!-----x

    x(n)=y(n)/u(n)

    DO i=n-1,1,-1
        x(i)=(y(i)-c(i)*x(i+1))/u(i)
    ENDDO
    END SUBROUTINE chase

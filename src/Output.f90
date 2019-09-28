! ====================================================================
!   Subroutine BalanceT
!
!   Purpose: The mass balance results, including water and solute.
! ====================================================================                        
! =========================Incoming variables=========================
!   voluz       The water volumn at the beginning.
!   voluz1      The water volumn at t time.
!   Svoluz      The solute volumn at the beginning.
!   Smvoluz     The mobile region solute volumn at the beginning.
!   Simvoluz    The immobile region solute volumn at the beginning.
!   dvol        The difference between voluz and voluz1.
!   qairt       Sum of precipitaiton.
!   qirrt       Sum of irrigation.
!   qbtmt       Sum of bottom flux.
!   sinkt       Actual evapotranspiration.
!   CumE        Sum of actual evaporation.
!   CumT        Sum of actual transpiration.
!   wbalt       Absolute water balance error.
!   wbalr       Relative water balance error.
!   sink1d(:)   Actual evaporation and transpiration at time t for
!               each layer.
! =========================Outcoming variables========================
!   None
! =========================Related files==============================
!   89(balance1d.dat)    The whole mass balance results.
! ====================================================================
    SUBROUTINE BalanceT
    USE parm

    REAL (kind=KR) :: dvol,Sdvol
    REAL (kind=KR) :: wbalt,wbalr,Swbalt,Swbalr

!   calculate the total water storage[m]
    voluz1=0.0_KR
    DO j=1,Nlayer
        voluz1=voluz1+th(j)*dz(j)
    !     M=matuz(j)  modified by maowei 2016-10-09 16:34
    !     if(huz(i,j)>0) voluz1=voluz1+SSS(M)*huz(i,j)*dz(j)
    ENDDO
!    calculate the variation of the storage between t=(0,t)
    dvol=voluz1-voluz
    qairt=qairt+qair*dt
    qbtmt=qbtmt+TotalPerco
    DO j=1,Nlayer
        sinkt=sinkt+sink1d(j)
    ENDDO
    CumE=CumE+Epa
    CumT=CumT+Tra
    wbalt=dvol-qairt+qbtmt+sinkt
    wbalr = abs(wbalt*1D2)/MAX(abs(dvol),qair+qbtmt+sinkt)  !2017-04-05
    IF (wbalt < Tol .and. abs(wbalr - 100._KR)<Tol) THEN
        wbalr = 0.0_KR
    ENDIF
    WRITE(89,"(f7.3,10E20.6E3)",err=901) t,voluz1,dvol,qairt,qbtmt,CumE,CumT,wbalt,wbalr
    
    RETURN
    
901 Terr=1
    RETURN
    END SUBROUTINE BalanceT

! ====================================================================
!   Subroutine Result_Out
!
!   Purpose: Print the water content and solute concentration of each
!            layer.
! ====================================================================                        
! =========================Incoming variables=========================
!   None
! =========================Outcoming variables========================
!   None
! =========================Related files==============================
!   80(thObs.dat)       All of the results.
!   81(profile.dat)     The results at print time.    
! ====================================================================
    SUBROUTINE hthuz_out
    USE parm

    WRITE(80,"(1X,F10.4\)",err=901) t,(th(Obs(i)),i=1,NObs)
    WRITE(80,"(1X,F10.4)") !th(Nlayer) 
    !WRITE(81,*)'Zone T="', t  
    !DO j=1,Nlayer
    !    WRITE(81,*)th(j),zx(j)+dz(j)/2,dz(j)
    !ENDDO
    RETURN

901 Terr=1
    RETURN

    END SUBROUTINE hthuz_out
    
    SUBROUTINE thOut
    USE parm
    
    WRITE(81,*,err=901)'Zone T="', t
    WRITE(81,*,err=901)'th,  thickness'
    DO j=1,Nlayer
        WRITE(81,"(2F10.4)",err=901)th(j),dz(j)
    ENDDO
    RETURN

901 Terr=1
    RETURN
    END SUBROUTINE

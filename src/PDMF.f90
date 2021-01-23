SUBROUTINE PDMF(IGIMF1,IGIMF2,IGIMF3)

    INTEGER :: i , t
    REAL :: ALLOCATABLE, DIMENSION(:) :: PDMF1, PDMF2, PDMF3
    REAL :: ALLOCATABLE, DIMENSION(:) :: m, logM, IGIMF1, IGIMF2, IGIMF3, REM_M, M_MAXPA !!!M_max padova
    REAL :: t, dt , gal_age, pad_ini , pad-stp


    CALL input(ncount, mcount, SFR, FeH, M_L, M_tu1, M_tu2, M_U, alpha1, alpha2, ak1, ak2)
    ALLOCATE(m(1:ncount + 1), logM(1:ncount + 1))
    m = 0.0; logM = 0.0
    ALLOCATE(IGIMF1(1:ncount+1),IGIMF2(1:ncount+1),IGIMF3(1:ncount+1))
    IGIMF1 = 0.0; IGIMF2 = 0.0; IGIMF3 = 0.0 
    ALLOCATE(PDMF1(1:ncount+1),PDMF2(1:ncount+1),PDMF3(1:ncount+1), REM_M(1:ncount+1), M_MAXPA(1:ncount+1))
    PDMF1 = 0.0; PDMF2 = 0.0; PDMF3 = 0.0, M_MAXPA = 0.0, REM_M =0.0 
    READ(*,)
    pad_stp=100
    dt = (gal_age - pad_ini)/pad_stp
    t = 0.0
    DO f = 1 , pad_stp 
        t = dt*f 
        DO i = 2, 170
            m(1) = 0.083
            logM(i) = -1.08 + 0.02*(i)  !-1.08+0.02*i
            m(i) = 10.0**(logM(i))
            IF (m(i) <= M_MAXPA(f)) THEN
                PDMF1(f) = IGIMF1(i) + PDMF1(f-1)
            ELSEIF
                IF(m(i) >= 40.0) THEN
			        REM_M(f) = 0.5*m(i)
		        ELSEIF (m(i) >= 8.5 .and. m(i) < 40.0) THEN
			        REM_M(f) = 1.4
		        ELSEIF(m(i) < 8.5) THEN
			        REM_M(f) = (0.077*m(i))+0.48
                ENDIF
		    ENDIF
        ENDDO 
        PDMF_T = PDMF1(f) + REM_M(f) + PDMF_T
    ENDDO
    CLOSE()
    CLOSE()

    RETURN
END     
    
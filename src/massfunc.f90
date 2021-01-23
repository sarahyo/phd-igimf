SUBROUTINE massfunc(Mecl_max, ncount, mcount, SFR, FeH, M_L, M_tu1, M_tu2, M_U, alpha1, alpha2, ak1, ak2, Kecl)

   INTEGER :: i, f
   INTEGER :: ncount, mcount
   REAL, ALLOCATABLE, DIMENSION(:) :: m, dm, logM
   REAL, ALLOCATABLE, DIMENSION(:) :: Mecl, dMecl, Mecl_ave, m_max, m_max1
   REAL, ALLOCATABLE, DIMENSION(:) :: IMF, XI, IGIMF1, IGIMF2, IGIMF3, IX, IMFS, XI2
   REAL :: SFR, FeH, Mecl_max, Nt_IGIMF3, Nt_IMF, Mt_IMF, Mt_IGIMF3, Mt_IGIMF2, Nt_IGIMF2, Nt_IGIMF1, Mt_IGIMF1
   REAL :: M_L, M_tu1, M_tu2, M_U, alpha1, alpha2, K1, K2, K3, K_1, K_2, K_3, alpha3, ak1, ak2, ak3, K11, K22, K33
   REAL :: Beta, delta, m_max2, Kecl

   CALL input(ncount, mcount, SFR, FeH, M_L, M_tu1, M_tu2, M_U, alpha1, alpha2, ak1, ak2)
   ALLOCATE (m(1:ncount + 1), dm(1:ncount + 1), logM(1:ncount + 1))
   m = 0.0; dm = 0.0; logM = 0.0
  ALLOCATE(Mecl(1:ncount+1),dMecl(1:ncount+1),Mecl_ave(1:ncount+1),m_max(1:ncount+1),m_max1(1:ncount+1),IX(1:ncount+1))
   Mecl = 0.0; dMecl = 0.0; Mecl_ave = 0.0; m_max = 0.0; IX = 0.0; m_max1 = 0.0
  ALLOCATE(IMF(1:ncount+1),IMFS(1:ncount+1),XI(1:ncount+1),XI2(1:ncount+1),IGIMF1(1:ncount+1),IGIMF2(1:ncount+1),IGIMF3(1:ncount+1))
   IMF = 0.0; XI = 0.0; XI2 = 0.0; IGIMF1 = 0.0; IGIMF2 = 0.0; IGIMF3 = 0.0; IMFS = 0.0   

   WRITE (500, *) '#', 'M_ecl', '   ', 'm_max'
   WRITE (700,*) 'Nt1' ,'Nt2', 'Nt3' , 'Mt1' , 'Mt2', 'Mt3'
   Nt_IGIMF3 = 0.0 ; Nt_IMF = 0.0
   Mecl(1) = 5.0 ; Mt_IGIMF3 = 0.0 ; Mt_IMF = 0.0
   Mt_IGIMF2 = 0.; Nt_IGIMF1 = 0.; Mt_IGIMF1 = 0.
   logM(1) = -1.08
   m(1) = 10.0**(logM(1))
   Beta = -0.106*log10(SFR) + 2.
   ak3 = 2.3

!==================================================================
   delta = (log10(Mecl_max) - log10(5.))/mcount
   DO f = 1, mcount + 2

      Mecl(f) = 10.0**(delta*f + log10(5.))
      CALL max_mstar(Mecl(f), ncount, SFR, FeH, M_L, M_tu1, M_tu2, M_U, alpha1, alpha2, m_max2)
      m_max(f) = m_max2
   ENDDo
!====================================================================
   Beta = -0.106*log10(SFR) + 2.
   DO i = 2, 170
      m(1) = 0.083
      logM(i) = -1.08 + 0.02*(i)  !-1.08+0.02*i
      m(i) = 10.0**(logM(i))
      dm(i) = m(i) - m(i - 1)
      IMFS(i) = m(i)**(-ak2)    !!!! Saplpeter IMF

      Nt_IMF = IMF(i)*dm(i) + Nt_IMF
      Mt_IMF = IMF(i)*dm(i)*m(i) + Mt_IMF

      Do f = 2, mcount + 2
         Mecl(f) = 10.0**(delta*f + log10(5.))
         IF (Mecl(f) <= Mecl_max) THEN
            Mecl_ave(f) = (Mecl(f) + Mecl(f - 1))/2.0
            IX(f) = 0.99*((0.61*log10(Mecl(f)) + 2.85) - 6.0) - 0.14*FeH
            dMecl(f) = Mecl(f) - Mecl(f - 1)
 
            IF (IX(f) >= -0.87) THEN
               alpha3 = 1.94 - 0.41*(IX(f))
            ELSEIF (IX(f) < -0.87) THEN
               alpha3 = 2.3
            ENDIF
!################################################## IGIMF3
            K11 = 0.0; K22 = 0.0; K33 = 0.; 
            K11 = Mecl(f)/((M_tu1**(-alpha1 + 2.) - M_L**(-alpha1 + 2.))/(-alpha1 + 2.) + (M_tu2**(-alpha2 + 2.) -&
               &M_tu1**(-alpha2 + 2.))*(M_tu1**(alpha2 - alpha1))/(-alpha2 + 2.) + (m_max(f)**(-alpha3 + 2.)&
               &- M_tu2**(-alpha3 + 2.))*(M_tu2**(alpha3 - alpha2))*(M_tu1**(alpha2 - alpha1))/(-alpha3 + 2.))
            K22 = K11*(M_tu1**(alpha2 - alpha1))
            K33 = k22*(M_tu2**(-alpha3 + alpha2))
! ################################################ IGIMF2
            K1=0.0 ; K2 =0.0; K3=0.0;
            K1 = Mecl(f)/((M_tu1**(-ak1 + 2.) - M_L**(-ak1 + 2.))/(-ak1 + 2.) + (M_tu2**(-ak2 + 2.) -&
            &M_tu1**(-ak2 + 2.))*(M_tu1**(ak2 - ak1))/(-ak2 + 2.) + (m_max(f)**(-alpha3 + 2.)&
            &- M_tu2**(-alpha3 + 2.))*(M_tu2**(alpha3 - ak2))*(M_tu1**(ak2 - ak1))/(-alpha3 + 2.))
            K2 = K1*(M_tu1**(ak2 - ak1))
            K3 = K2*(M_tu2**(alpha3 - ak2))
! #################################################IGIMF1
            K_1=0.0; k_2=0.0; K_3=0.0;
            k_2 = (Mecl(f))/(2.*(M_tu1**(-ak1 + 2.) - M_L**(-ak1 + 2.))/(-ak1 + 2.) + &
            &(M_U**(-ak2 + 2.) - M_tu1**(-ak2 + 2.))/(-ak2 + 2.))
            K_1 = k_2*(M_tu1**(ak1 - ak2))
            K_3 = K_2*(M_tu2**(ak2 - ak3))
! -----------------------------------------------------------------------------------------
            IF (m(i) >= 0.08 .AND. m(i) < 0.5) THEN
               IMF(i) = K_1*m(i)**(-ak1)  ! IGIMF1
               XI2(i) = K1*m(i)**(-ak1)   ! IGIMF2
               XI(i) = K11*m(i)**(-alpha1) ! IGIMF3

            ENDIF
            IF (m(i) >= 0.5 .AND. m(i) <= 1.0) THEN
               IMF(i) = k_2*m(i)**(-ak2)
               XI2(i) = K2*m(i)**(-ak2)
               XI(i) = K22*m(i)**(-alpha2)
            ENDIF
            IF (m(i) >= 1.0 .and. m(i) <= 150.0) THEN
               IMF(i) = K_3*m(i)**(-ak3)
               XI2(i) = K3*m(i)**(-alpha3)
               XI(i) = K33*m(i)**(-alpha3)
            ENDIF

            IF (m(i) <= m_max(f)) THEN
               IGIMF1(i) = IMF(i)*((Mecl_ave(f))**(-Beta))*dMecl(f) + IGIMF1(i)
               IGIMF2(i) = XI2(i)*((Mecl_ave(f))**(-Beta))*dMecl(f) + IGIMF2(i)
               IGIMF3(i) = XI(i)*((Mecl_ave(f))**(-Beta))*dMecl(f) + IGIMF3(i)
               ! print*,"here",IGIMF1(i), IGIMF2(i), IGIMF3(i), K_3, K1, K11
            ENDIF


         ENDIF

      ENDDO

      Nt_IGIMF1 = IGIMF1(i)*dm(i) + Nt_IGIMF1
      Mt_IGIMF1 = IGIMF1(i)*dm(i)*m(i) + Mt_IGIMF1

      Nt_IGIMF2 = IGIMF2(i)*dm(i) + Nt_IGIMF2
      Mt_IGIMF2 = IGIMF2(i)*dm(i)*m(i) + Mt_IGIMF2

      Nt_IGIMF3 = IGIMF3(i)*dm(i) + Nt_IGIMF3
      Mt_IGIMF3 = IGIMF3(i)*dm(i)*m(i) + Mt_IGIMF3

   ENDDO
   WRITE (700, *) Nt_IGIMF1, Nt_IGIMF2, Nt_IGIMF3, Mt_IGIMF1, Mt_IGIMF2, Mt_IGIMF3  
   DO i = 2, 170
      IF (m(i) <= 150.0) THEN
         WRITE (200, *) m(i), IGIMF1(i), IGIMF2(i), IGIMF3(i)
         WRITE (240, *) m(i), IMF(i), IMFS(i)
         WRITE (20, *) m(i), IGIMF1(i), IGIMF2(i), IGIMF3(i)
      ENDIF
   ENDDO
   
   CLOSE (200)
   CLOSE (500)
   CLOSE (240)

   RETURN
END


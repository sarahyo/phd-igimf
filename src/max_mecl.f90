SUBROUTINE max_mecl(ncount, SFR, FeH, Mecl_max, Kecl)

   INTEGER :: f
   INTEGER :: ncount
   REAL, ALLOCATABLE, DIMENSION(:) :: Mecl
   REAL :: Beta, Mecl_u, Mecl_min, Kecl, SFR, Mecl_max, FeH, Mtot, M_t

   CALL input(ncount, mcount, SFR, FeH, M_L, M_tu1, M_tu2, M_U, alpha1, alpha2, ak1, ak2)
   ALLOCATE (Mecl(1:ncount + 1))
   Mecl = 0.0

!-----------------------------------------------------
   write(900,*)"f","M_t","Mtot","Mecl(f)","Mecl_max(i)" 
   Beta = -0.106*log10(SFR) + 2.
   Mecl(1) = 5.0
   M_t = 0.
   Mtot = SFR*10.**7.  !Mtot=SFR*dt
   Mecl_u = 10.**9.
   Mecl_min = 5.

   outer: DO f = 2, 900
      M_t = 0.
      Mecl(f) = Mecl_min + 10**(0.01*f)
      IF (Beta == 2.) THEN
         M_t = (log(Mecl(f)) - log(Mecl_min))/(-1.*Mecl(f)**(-1.) + Mecl_u**(-1.))
      ELSE
         M_t = ((1.-Beta)/(2.-Beta))*(Mecl(f)**(2.-Beta) - Mecl_min**(2.-Beta))*((Mecl_u**(1.-Beta) - Mecl(f)**(1.-Beta))**(-1.))
      ENDIF
      
      !   print*,Mtot,M_t, Mecl(f), Beta, SFR
    write(900,*)f,M_t,Mtot,Mecl(f),Mecl_max
      IF (M_t >= Mtot) THEN
         Mecl_max = Mecl(f)
         exit outer
      ENDIF

      !  WRITE(800,*)Mecl_max,Kecl,SFR
    write(900,*)f,M_t,Mtot,Mecl(f),Mecl_max
   ENDDO outer

   IF (M_t < Mtot) THEN
      Mecl_max = 1.e9
   ENDIF

   Kecl = Mtot*(-Beta + 1.)/(Mecl_max**(-Beta + 1.) - Mecl_min**(-Beta + 1.0))
   print *,"inja", Kecl, Mecl_max

   WRITE (800, *) Mecl_max, Kecl, SFR
   CLOSE (800)

   RETURN
END SUBROUTINE


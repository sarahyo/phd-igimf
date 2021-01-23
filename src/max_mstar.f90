SUBROUTINE max_mstar(Mecl1,ncount,SFR,FeH,M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,m_max2)

  INTEGER :: j
  INTEGER :: ncount
  REAL, ALLOCATABLE, DIMENSION(:) :: m,dm,logM
  REAL :: SFR,FeH,M_C
  REAL :: M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,alpha3,k1
  REAL :: m_max2,IX1,Mecl1

  CALL input(ncount,mcount,SFR,FeH,M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,ak1,ak2)  
  ALLOCATE(m(1:ncount+1),dm(1:ncount+1),logM(1:ncount+1)) 
           m=0.0; dm=0.0; logM=0.0
     

  logM(1)=-1.08
  m(1)=10.0**(logM(1))  
             
  IX1=0.99*(0.61*log10(Mecl1)+2.85-6.0)-0.14*FeH
  IF (IX1 >= -0.87 )THEN
     alpha3=(1.0*(1.94-0.41*IX1))
  ELSEIF(IX1 < -0.87)THEN
     alpha3=2.3
  ENDIF

  
  M_C=0.0     
  DO j=2,400

     logM(j)=-1.08+0.01*(j-1)  
     m(j)=10.0**(logM(j))
     dm(j)=m(j)-m(j-1)
     M_C=0.0                  


     IF (m(j) >= 0.08.AND.m(j) < 0.5) THEN
        K1=1./((M_tu1**(1.-alpha1)-m(j)**(1.-alpha1))/(1.-alpha1)+(M_tu2**(1.-alpha2)-&
           &M_tu1**(1.-alpha2))*(M_tu1**(alpha2-alpha1))/(1.-alpha2)+(M_U**(1.-alpha3)-&
           &M_tu2**(1.-alpha3))*(M_tu1**(alpha2-alpha1))*(M_tu2**(alpha3-alpha2))/(1.-alpha3)) 
                   
        M_C=K1*(m(j)**(2.-alpha1)-M_L**(2.-alpha1))/(2.-alpha1)              
     ENDIF 
     IF (m(j) >= 0.5.AND.m(j) < 1.0) THEN
        K1=1./( (M_tu2**(1.-alpha2)-m(j)**(1.-alpha2))*(M_tu1**(alpha2-alpha1))/(1.-alpha2)+&
            &(M_U**(1.-alpha3)-M_tu2**(1.-alpha3))*(M_tu1**(alpha2-alpha1))*&
            &(M_tu2**(alpha3-alpha2))/(1.-alpha3))
                   
        M_C=K1*((M_tu1**(2.-alpha1)-M_L**(2.-alpha1))/(2.-alpha1)+(m(j)**(2.-alpha2)-&
            &M_tu1**(2.-alpha2))*(M_tu1**(alpha2-alpha1))/(2.-alpha2))             
     ENDIF 
     IF (m(j) >= 1..AND.m(j) < 150.) THEN
        K1=(M_tu1**(alpha1-alpha2))*(M_tu2**(alpha2-alpha3))*(1.-alpha3)/&
            &(M_U**(1.-alpha3)-m(j)**(1.-alpha3))
                   
        M_C=K1*((M_tu1**(2.-alpha1)-M_L**(2.-alpha1))/(2.-alpha1)+(M_tu2**(2.-alpha2)-&
            &M_tu1**(2.-alpha2))*(M_tu1**(alpha2-alpha1))/(2.-alpha2)+(M(j)**(2.-alpha3)-&
            &M_tu2**(2.-alpha3))*(M_tu2**(alpha3-alpha2))*(M_tu1**(alpha2-alpha1))/(2.-alpha3))             
     ENDIF 
      
     IF(M_C >= Mecl1) THEN
      m_max2=m(j)
      exit
     ENDIF
   !   print*,m_max1
       
  ENDDO

  IF (M_C < Mecl1)THEN
     m_max2=150.
  ENDIF
!   print*,m_max1

  WRITE(500,*)Mecl1,m_max2 
  RETURN   
END  




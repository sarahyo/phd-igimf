SUBROUTINE mmax_IG1(Mecl1,ncount,SFR,FeH,M_L,M_tu1,M_tu2,M_U,aK1,aK2,mmax_I1)

    INTEGER :: j
    INTEGER :: ncount
    REAL, ALLOCATABLE, DIMENSION(:) :: m,dm,logM
    REAL :: SFR,FeH,M_C
    REAL :: M_L,M_tu1,M_tu2,M_U,ak1,aK2,aK3,k1
    REAL :: mmax_I1,Mecl1
  
    CALL input(ncount,mcount,SFR,FeH,M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,ak1,ak2)  
    ALLOCATE(m(1:ncount+1),dm(1:ncount+1),logM(1:ncount+1)) 
             m=0.0; dm=0.0; logM=0.0
       
  
    logM(1)=-1.08
    m(1)=10.0**(logM(1))  
    aK3=2.3           
    
    M_C=0.0     
    DO j=2,400
  
       logM(j)=-1.08+0.01*(j-1)  
       m(j)=10.0**(logM(j))
       dm(j)=m(j)-m(j-1)
       M_C=0.0                  
  
  
       IF (m(j) >= 0.08.AND.m(j) < 0.5) THEN
          K1=1./((M_tu1**(1.-ak1)-m(j)**(1.-ak1))/(1.-ak1)+(M_tu2**(1.-aK2)-&
             &M_tu1**(1.-ak2))*(M_tu1**(aK2-ak1))/(1.-aK2)+(M_U**(1.-aK3)-&
             &M_tu2**(1.-aK3))*(M_tu1**(aK2-ak1))*(M_tu2**(aK3-aK2))/(1.-aK3)) 
                     
          M_C=K1*(m(j)**(2.-ak1)-M_L**(2.-ak1))/(2.-ak1)              
       ENDIF 
       IF (m(j) >= 0.5.AND.m(j) < 1.0) THEN
          K1=1./( (M_tu2**(1.-aK2)-m(j)**(1.-aK2))*(M_tu1**(aK2-ak1))/(1.-aK2)+&
              &(M_U**(1.-aK3)-M_tu2**(1.-aK3))*(M_tu1**(aK2-ak1))*&
              &(M_tu2**(aK3-aK2))/(1.-aK3))
                     
          M_C=K1*((M_tu1**(2.-ak1)-M_L**(2.-ak1))/(2.-ak1)+(m(j)**(2.-aK2)-&
              &M_tu1**(2.-aK2))*(M_tu1**(aK2-ak1))/(2.-aK2))             
       ENDIF 
       IF (m(j) >= 1..AND.m(j) < 150.) THEN
          K1=(M_tu1**(ak1-aK2))*(M_tu2**(aK2-aK3))*(1.-aK3)/&
              &(M_U**(1.-aK3)-m(j)**(1.-aK3))
                     
          M_C=K1*((M_tu1**(2.-ak1)-M_L**(2.-ak1))/(2.-ak1)+(M_tu2**(2.-aK2)-&
              &M_tu1**(2.-aK2))*(M_tu1**(aK2-ak1))/(2.-aK2)+(M(j)**(2.-aK3)-&
              &M_tu2**(2.-aK3))*(M_tu2**(aK3-aK2))*(M_tu1**(aK2-ak1))/(2.-aK3))             
       ENDIF 
        
       IF(M_C >= Mecl1) THEN
        mmax_I1=m(j)
        exit
       ENDIF
     !   print*,mmax_I1
         
    ENDDO
  
    IF (M_C < Mecl1)THEN
       mmax_I1=150.
    ENDIF
   !  print*,mmax_I1
  
    ! WRITE(500,*)Mecl1,mmax_I1 
    RETURN   
  END  
  
  
  
  

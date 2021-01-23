subroutine input(ncount,mcount,SFR,FeH,M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,ak1,ak2)
  IMPLICIT NONE

     REAL SFR,FeH,M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,ak1,ak2 !!!!!! ak1,2=alpha1,2 korupa
     INTEGER ncount,mcount
     open(10,file='input')
     

     READ(10,*)SFR,FeH
     READ(10,*)ncount,mcount
     READ(10,*)M_L,M_tu1,M_tu2,M_U
     READ(10,*)ak1,ak2
     
    
      
     
     alpha1=1.30+0.5*FeH
     alpha2=2.30+0.5*FeH
     

     CLOSE(10)
RETURN
END


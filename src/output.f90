SUBROUTINE output(FeH,SFR)

  REAL SFR,FeH
  CHARACTER(LEN=45) :: OUTPUT_FILE
  CHARACTER(LEN=7)  :: n2        !SFR
  CHARACTER(LEN=4)  :: n4        !Fe/H 


     WRITE(n2,'(ES7.1)')SFR
     WRITE(n4,'(F4.1)')FeH
     OUTPUT_FILE='SFR'//n2//'FeH'//n4//".txt" 
   


  OPEN(200,FILE=OUTPUT_FILE) 
  OPEN(240,FILE='canonic-IMF.txt')
  OPEN(500,FILE='mmax-Mecl.txt')
  OPEN(700,FILE='Ntot,Mtot.txt')
  open(20,file='IGIMFS')
   
RETURN
END

PROGRAM GWIMF

  IMPLICIT NONE
!   INTEGER :: i,j,f,k
  INTEGER :: ncount,mcount
  REAL, ALLOCATABLE, DIMENSION(:) :: m,dm
  REAL, ALLOCATABLE, DIMENSION(:) :: Mecl,dMecl,Meclave,m_max,IMF,logM
  REAL, ALLOCATABLE, DIMENSION(:) :: phi,IGIMF,N,NN,IX
  REAL :: SFR,Mecl_max,FeH, Kecl
  REAL ::M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,ak1,ak2

  CALL input(ncount,mcount,SFR,FeH,M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,ak1,ak2)
  CALL output(FeH,SFR)

  ALLOCATE(m(1:ncount+1),dm(1:ncount+1),logM(1:ncount+1)) 
     m=0.0; dm=0.0; logM=0.0
  ALLOCATE(Mecl(1:ncount+1),dMecl(1:ncount+1),Meclave(1:ncount+1),m_max(1:ncount+1),IMF(1:ncount+1))
     Mecl=0.0; dMecl=0.0; Meclave=0.0; m_max=0.0; IMF=0.0
  ALLOCATE(phi(1:ncount+1),IGIMF(1:ncount+1),N(1:ncount+1),NN(1:ncount+1),IX(1:ncount+1))
     phi=0.0; IGIMF=0.0; N=0.0; NN=0.0; IX=0.0

  CALL max_mecl(ncount,SFR,FeH,Mecl_max,Kecl)
  CALL massfunc(Mecl_max,ncount,mcount,SFR,FeH,M_L,M_tu1,M_tu2,M_U,alpha1,alpha2,ak1,ak2,Kecl)
  
 
END PROGRAM GWIMF

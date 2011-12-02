CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C       --U S E R   E L E M E N T    S U B R O U T I N E S---          C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  UEL_FEM8.for                                                        C
C  8-NODED ISOPARAMETRIC ELEMENT                                       C
C  UNIVERSIDAD EAFIT                                                   C
C  2011                                                                C
C  MECANICA APLICADA                                                   C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  USER ELEMENT SUBROUTINE-ISOTROPIC ELASTICITY                        C
C                                                                      C
C  STANDARD  ISOPARAMETRIC ELEMENT                                     C
C  FULL GAUSS INTEGRATION                                              C
C  CREATED BY JUAN GOMEZ                                               C
C  STATE VARIABLES DEFINITIONS/GAUSS  POINT                            C
C                                                                      C
C                                                                      C
C     1-4: STRESS VECTOR                       (4)                     C
C     5-8: TOTAL STRAIN VECTOR                 (4)                     C
C                                   TOTAL      72 SVARS                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL_FEMB8(AMATRX,AMASS,NDOFEL,PROPS,NPROPS,COORDS,MCRD,
     1                     NNODE,JELEM)
C     
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.D0,NTENS=4,TWO=2.D0,NGPTS=9,
     1           THREE=3.D0)
C
C     Parameters required by UMAT.f
C
      PARAMETER (NDI=3)
C
C     Parameter arrays from UEL.f
C
      DIMENSION AMATRX(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL),
     1          PROPS(NPROPS),COORDS(MCRD,NNODE)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),
     1          FRST1(NDOFEL,NDOFEL),XX(MCRD,NNODE),
     2          XW(NGPTS),XP(2,NGPTS),AUX1(NTENS,NDOFEL),
     3          FRST3(NDOFEL,NDOFEL)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION DDSDDE(NTENS,NTENS)
C
C     Initializes parameters required by UMAT.f
C
      CALL CLEAR(DDSDDE,NTENS,NTENS)
C
C     Clears RHS vector and Stiffness matrix
C
      CALL CLEAR(AMATRX,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
C**********************************************************************
C     P U R E  D I S P L A C E M E N T  F O R M U L A T I O N
C**********************************************************************
C
      RO=PROPS(3)
C
C     Generates Gauss points and weights.
C
      CALL GPOINTS3X3(XP,XW)
      NGPT=9
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN)
C
C       Assembles B matrix
C
        CALL STDM8(JELEM,NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,SII,
     1             XBAR)
C
C       Assembles material matrix and updates state variables
C
        CALL UMATELA(DDSDDE,NDI,NTENS,PROPS,NPROPS)
C
        CALL MMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL CLEAR(BT,NDOFEL,NTENS)
        CALL MTRAN(B,NTENS,NDOFEL,BT)
C
C       Assembles stiffness matrix and RHS vector contribution
C
        CALL MMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL SMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
C
        CALL AMASS8(FRST3,NDOFEL,RII,SII)
        CALL SMULT(FRST3,NDOFEL,NDOFEL,RO*ALF*DDET*XBAR)
C
C       Updates Stiffness matrix and RHS vector
C
        CALL UPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL UPDMAT(AMASS,NDOFEL,NDOFEL,FRST3)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
        CALL CLEAR(FRST1,NDOFEL,NDOFEL)
        CALL CLEAR(FRST3,NDOFEL,NDOFEL)
C
      END DO
C
      RETURN
C      
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  UEL_FEM9.for 9-Noded Plasticity Isotropic Classical Return Mapp.    C
C                                                                      C
C  9-NODED ISOPARAMETRIC ELEMENT                                       C
C  UNIVERSIDAD EAFIT                                                   C
C  2010                                                                C
C  MECANICA APLICADA                                                   C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  USER ELEMENT SUBROUTINE-ISOTROPIC HARDENING MISES PLASTICITY        C
C  AND ISOTROPIC ELASTICITY                                            C
C                                                                      C
C  STANDARD  ISOPARAMETRIC ELEMENT                                     C
C  FULL GAUSS INTEGRATION                                              C
C  CREATED BY JUAN GOMEZ                                               C
C  STATE VARIABLES DEFINITIONS/GAUSS  POINT                            C
C                                                                      C
C                                                                      C
C     1-4: STRESS VECTOR                       (4)                     C
C     5-8: TOTAL STRAIN VECTOR                 (4)                     C
C                                   TOTAL      72 SVARS                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL_FEM9(RHS,AMATRX,SVARS,NDOFEL,NSVARS,PROPS,NPROPS,
     1                    COORDS,MCRD,NNODE,U,DU,V,A,JELEM)
C     
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.D0,NTENS=4,TWO=2.D0,NGPTS=9,
     1           THREE=3.D0)
C
C     Parameters required by UMAT.f
C
      PARAMETER (NDI=3)
C
C     Parameter arrays from UEL.f
C
      DIMENSION RHS(NDOFEL,1),AMATRX(NDOFEL,NDOFEL),SVARS(NSVARS),
     1          PROPS(NPROPS),COORDS(MCRD,NNODE),U(NDOFEL),
     2          DU(NDOFEL),V(NDOFEL),A(NDOFEL)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),
     1          FRST1(NDOFEL,NDOFEL),FRST2(NDOFEL,1),XX(MCRD,NNODE),
     2          XW(NGPTS),XP(2,NGPTS),AUX1(NTENS,NDOFEL),STRESS(NTENS),
     3          STRAN(NTENS),DSTRAN(NTENS)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION DDSDDE(NTENS,NTENS)
C
C     Initializes parameters required by UMAT.f
C
      CALL CLEARV(STRESS,NTENS)
      CALL CLEARV(STRAN,NTENS)
      CALL CLEARV(DSTRAN,NTENS)
      CALL CLEAR(DDSDDE,NTENS,NTENS)
C
C     Clears RHS vector and Stiffness matrix
C
      CALL CLEARV(RHS,NDOFEL)
      CALL CLEAR(AMATRX,NDOFEL,NDOFEL)
C
C**********************************************************************
C     P U R E  D I S P L A C E M E N T  F O R M U L A T I O N
C**********************************************************************
C
C
      NSVARS_N=0
      NSVARS_F=NSVARS-NSVARS_N
C
C     Generates Gauss points and weights.
C
      CALL GPOINTS3X3(XP,XW)
      NGPT=9
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN)
C
C       Compute State variable index corresponding
C       to current Gauss point and load stress,total strain
C       from state variables as defined in USER ELEMENT. 
C       Different variables are required for different constitutive models.
C
        ISVINT_F=1+(NN-1)*NSVARS_F/NGPT
        JJ=1
        DO II=ISVINT_F,ISVINT_F+3
          STRESS(JJ)=SVARS(II)
          STRAN(JJ)=SVARS(II+4)
          JJ=JJ+1
        END DO
C
C       Starts STATEV array definition as required by UMAT.f
C       Different variables for different constitutive models.
C
C
C       Ends STATEV array definition as required by UMAT.f
C
C       Assembles B matrix
C
        CALL STDM9(JELEM,NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,SII,
     1             XBAR)
        CALL MAVEC(B,NTENS,NDOFEL,DU,DSTRAN)
        CALL UPDVEC(STRAN,NTENS,DSTRAN)
C
C       Assembles material matrix and updates state variables
C
ccc        CALL UMATELA(STRESS,DDSDDE,DSTRAN,NDI,NTENS,PROPS,NPROPS)
C
        CALL MMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL CLEAR(BT,NDOFEL,NTENS)
        CALL MTRAN(B,NTENS,NDOFEL,BT)
C
C       Assembles stiffness matrix and RHS vector contribution
C
        CALL MMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
        CALL MAVEC(BT,NDOFEL,NTENS,STRESS,FRST2)
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL SMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL SMULT(FRST2,NDOFEL,1,-ALF*DDET*XBAR)
C
C       Updates Stiffness matrix and RHS vector
C
        CALL UPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL UPDVEC(RHS,NDOFEL,FRST2)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
        CALL CLEAR(FRST1,NDOFEL,NDOFEL)
C
C       Starts updating of state variables with updated values from UMAT.f
C
        JJ=1
        DO II=ISVINT_F,ISVINT_F+3
          SVARS(II)=STRESS(JJ)
          SVARS(II+4 )=STRAN(JJ)
          JJ=JJ+1
        END DO
C
C       Ends updating of state variables with updated values from UMAT.f
C
      END DO
C
      RETURN
C      
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  UEL_FEM6.for 6-Noded Plasticity Isotropic Classical Return Mapp.    C
C                                                                      C
C  6-NODED ISOPARAMETRIC ELEMENT                                       C
C  UNIVERSIDAD EAFIT                                                   C
C  2005                                                                C
C  MECANICA APLICADA                                                   C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  USER ELEMENT SUBROUTINE-ISOTROPIC HARDENING MISES PLASTICITY        C
C  AND ISOTROPIC ELASTICITY                                            C
C                                                                      C
C  STANDARD  ISOPARAMETRIC ELEMENT                                     C
C  FULL GAUSS INTEGRATION                                              C
C  CREATED BY JUAN GOMEZ                                               C
C  STATE VARIABLES DEFINITIONS/GAUSS  POINT                            C
C                                                                      C
C                                                                      C
C     1-4: STRESS VECTOR                       (4)                     C
C     5-8: TOTAL STRAIN VECTOR                 (4)                     C
C                                   TOTAL      56 SVARS                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL_FEM6(RHS,AMATRX,SVARS,NDOFEL,NSVARS,PROPS,NPROPS,
     1                    COORDS,MCRD,NNODE,U,DU,V,A,JELEM)
C     
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.D0,NTENS=4,TWO=2.D0,NGPTS=7,
     1           THREE=3.D0)
C
C     Parameters required by UMAT.f
C
      PARAMETER (NDI=3)
C
C     Parameter arrays from UEL.f
C
      DIMENSION RHS(NDOFEL,1),AMATRX(NDOFEL,NDOFEL),SVARS(NSVARS),
     1          PROPS(NPROPS),COORDS(MCRD,NNODE),U(NDOFEL),
     2          DU(NDOFEL),V(NDOFEL),A(NDOFEL)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),
     1          FRST1(NDOFEL,NDOFEL),FRST2(NDOFEL,1),XX(MCRD,NNODE),
     2          XW(NGPTS),XP(2,NGPTS),AUX1(NTENS,NDOFEL),STRESS(NTENS),
     3          STRAN(NTENS),DSTRAN(NTENS)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION DDSDDE(NTENS,NTENS)
C
C     Initializes parameters required by UMAT.f
C
      CALL CLEARV(STRESS,NTENS)
      CALL CLEARV(STRAN,NTENS)
      CALL CLEARV(DSTRAN,NTENS)
      CALL CLEAR(DDSDDE,NTENS,NTENS)
C
C     Clears RHS vector and Stiffness matrix
C
      CALL CLEARV(RHS,NDOFEL)
      CALL CLEAR(AMATRX,NDOFEL,NDOFEL)
C
C**********************************************************************
C     P U R E  D I S P L A C E M E N T  F O R M U L A T I O N
C**********************************************************************
C
C
      NSVARS_N=0
      NSVARS_F=NSVARS-NSVARS_N
C
C     Generates Gauss points and weights.
C
      CALL GPOINTS7(XP,XW)
      NGPT=7
C
C     Loops around all Gauss points
C
      DO NN=1,NGPT
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN)
C
C       Compute State variable index corresponding
C       to current Gauss point and load stress,total strain
C       from state variables as defined in USER ELEMENT. 
C       Different variables are required for different constitutive models.
C
        ISVINT_F=1+(NN-1)*NSVARS_F/NGPT
        JJ=1
        DO II=ISVINT_F,ISVINT_F+3
          STRESS(JJ)=SVARS(II)
          STRAN(JJ)=SVARS(II+4)
          JJ=JJ+1
        END DO
C
C       Starts STATEV array definition as required by UMAT.f
C       Different variables for different constitutive models.
C
C
C       Ends STATEV array definition as required by UMAT.f
C
C       Assembles B matrix
C
        CALL STDM6(JELEM,NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,SII,
     1             XBAR)
        CALL MAVEC(B,NTENS,NDOFEL,DU,DSTRAN)
        CALL UPDVEC(STRAN,NTENS,DSTRAN)
C
C       Assembles material matrix and updates state variables
C
ccc        CALL UMATELA(STRESS,DDSDDE,DSTRAN,NDI,NTENS,PROPS,NPROPS)
C
        CALL MMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL CLEAR(BT,NDOFEL,NTENS)
        CALL MTRAN(B,NTENS,NDOFEL,BT)
C
C       Assembles stiffness matrix and RHS vector contribution
C
        CALL MMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
        CALL MAVEC(BT,NDOFEL,NTENS,STRESS,FRST2)
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL SMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL SMULT(FRST2,NDOFEL,1,-ALF*DDET*XBAR)
C
C       Updates Stiffness matrix and RHS vector
C
        CALL UPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL UPDVEC(RHS,NDOFEL,FRST2)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
        CALL CLEAR(FRST1,NDOFEL,NDOFEL)
C
C       Starts updating of state variables with updated values from UMAT.f
C
        JJ=1
        DO II=ISVINT_F,ISVINT_F+3
          SVARS(II)=STRESS(JJ)
          SVARS(II+4 )=STRAN(JJ)
          JJ=JJ+1
        END DO
C
C       Ends updating of state variables with updated values from UMAT.f
C
      END DO
C
      RETURN
C      
      END
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
C  UEL_FEM8.for 8-Noded isoparametric 2d element.                      C
C                                                                      C
C  8-NODED ISOPARAMETRIC ELEMENT                                       C
C  UNIVERSIDAD EAFIT                                                   C
C  2012                                                                C
C  MECANICA APLICADA                                                   C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  USER ELEMENT SUBROUTINE                                             C
C  UEL_FEM8                                                            C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL_FEM8(AMATRX,AMASS,NDOFEL,PROPS,NPROPS,COORDS,MCRD,
     1                    NNODE,JELEM)

C     
      IMPLICIT REAL*8(A-H,O-Z)
C  
      PARAMETER (NTENS=4,NDI=3,NGPTS=9)
C
C     Parameter arrays from UEL.f
C
      DIMENSION AMATRX(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL),
     1          PROPS(NPROPS),COORDS(MCRD,NNODE),FRST2(NDOFEL,NDOFEL)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),FRST1(NDOFEL,NDOFEL),
     1          XX(MCRD,NNODE),XW(NGPTS),XP(2,NGPTS),AUX1(NTENS,NDOFEL)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION DDSDDE(NTENS, NTENS)
C
C     Clears RHS vector and Stiffness matrix
C
      CALL CLEAR(AMATRX,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL, NDOFEL)
C
      RO=PROPS(3)
C
C     Generates Gauss points and weights.
C
      CALL GPOINTS3X3(XP,XW)
C
C     Loops around all Gauss points
C
      DO NN=1,NGPTS
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN)
C
C       Assembles B matrix
C
        CALL STDM8(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,DDET,RII,SII,XBAR)
        CALL RMASS8(FRST2,NDOFEL,RII,SII)
C
C       Assembles material matrix and updates state variables
C
        CALL UMATELA(DDSDDE,PROPS, NPROPS,NTENS,NDI)
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
        CALL SMULT(FRST2,NDOFEL,NDOFEL,ALF*DDET*XBAR*RO)
C
C       Updates Stiffness and mass matrix.
C
        CALL UPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL UPDMAT(AMASS, NDOFEL,NDOFEL,FRST2)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point 
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
        CALL CLEAR(FRST1,NDOFEL,NDOFEL)
        CALL CLEAR(FRST2,NDOFEL,NDOFEL)
C
      END DO
C
      RETURN
C      
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  USER ELEMENT SUBROUTINE                                             C
C  UEL_FEM9                                                            C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL_FEM9(AMATRX,AMASS,NDOFEL,PROPS,NPROPS,COORDS,MCRD,
     1                    NNODE,JELEM)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=4,NDI=3,NGPTS=9)
C
C     Parameter arrays from UEL.f
C
      DIMENSION AMATRX(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL),
     1          PROPS(NPROPS),COORDS(MCRD,NNODE),FRST2(NDOFEL,NDOFEL)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),FRST1(NDOFEL,NDOFEL),
     1          XX(MCRD,NNODE),XW(NGPTS),XP(2,NGPTS),AUX1(NTENS,NDOFEL)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION DDSDDE(NTENS, NTENS)
C
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  UEL_COS8.for 8-Noded isoparametric 2d element.                      C
C                                                                      C
C  8-NODED ISOPARAMETRIC ELEMENT                                       C
C  UNIVERSIDAD EAFIT                                                   C
C  2012                                                                C
C  MECANICA APLICADA                                                   C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C  USER ELEMENT SUBROUTINE                                             C
C  UEL_FEM8                                                            C
C                                                                      C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UEL_COS8(AMATRX,AMASS,NDOFEL,PROPS,NPROPS,COORDS,MCRD,
     1                    NNODE,JELEM)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=7,NDI=3,NGPTS=9)
C
C     Parameter arrays from UEL.f
C
      DIMENSION AMATRX(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL),
     1          PROPS(NPROPS),COORDS(MCRD,NNODE),FRST2(NDOFEL,NDOFEL)
C
C     User defined arrays
C
      DIMENSION B(NTENS,NDOFEL),BT(NDOFEL,NTENS),FRST1(NDOFEL,NDOFEL),
     1          XX(MCRD,NNODE),XW(NGPTS),XP(2,NGPTS),AUX1(NTENS,NDOFEL)
C
C     Arrays to be used in UMAT.f
C
      DIMENSION DDSDDE(NTENS, NTENS)
C
C     Clears RHS vector and Stiffness matrix
C
      CALL CLEAR(AMATRX,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL, NDOFEL)
C
      RO=PROPS(3)
      RI=PROPS(4)
C
C     Generates Gauss points and weights.
C
      CALL GPOINTS3X3(XP,XW)
C
C     Loops around all Gauss points
C
      DO NN=1,NGPTS
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN)
C
C       Assembles B matrix
C
        CALL STDM8C(JELEM,NNODE,NDOFEL,NTENS,COORDS,B,DDET,RII,SII,XBAR)
C
C       Assembles material matrix and updates state variables
C
        CALL UMATCOS(DDSDDE,PROPS, NPROPS,NTENS,NDI)
C
        CALL MMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL CLEAR(BT,NDOFEL,NTENS)
        CALL MTRAN(B,NTENS,NDOFEL,BT)
C
C       Assembles stiffness matrix.
C
        CALL MMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
C
C       Temporary mass matrix.
C
        CALL RMASS8C(FRST2,RO,RI,NDOFEL,RII,SII)
C
C       Considers Gauss weight and Jacobian determinant representing
C       volume differential.
C
        CALL SMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
        CALL SMULT(FRST2,NDOFEL,NDOFEL,ALF*DDET*XBAR)
C
C       Updates Stiffness and mass matrix.
C
        CALL UPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
        CALL UPDMAT(AMASS ,NDOFEL,NDOFEL,FRST2)
C
C       Clears material Jacobian and temporary stiffness matrix
C       array for new Gauss point
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
        CALL CLEAR(FRST1,NDOFEL,NDOFEL)
        CALL CLEAR(FRST2,NDOFEL,NDOFEL)
C
      END DO
C
      RETURN
C
      END
C

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C--C O N S T I T U T I V E  M A T E R I A L S  S U B R O U T I N E S---C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UMATELA                                               C
C     PROPS(1) - E                                                     C
C     PROPS(2) - NU                                                    C
C     PROPS(3) - RO                                                    C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATELA(DDSDDE,PROPS,NPROPS,NTENS,NDI)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0)
C
      DIMENSION DDSDDE(NTENS, NTENS),PROPS(NPROPS)
C
C     ELASTIC PROPERTIES
C
      EMOD=PROPS(1)
      ENU=PROPS(2)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
C
C     ELASTIC STIFFNESS
C
      DO K1=1, NDI
        DO K2=1, NDI
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1, NTENS
        DDSDDE(K1, K1)=EG
      END DO
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UMATCOS                                               C
C                                                                       
C       IN THE INPUT FILE THE PARAMATERS ARE DEFINED AS:               C
C                                                                      C
C       ID  7   BETA    CPOIS   RHO J   XI  GAMMA   MU_C               C
C                                                                      C
C       BETA:   CLASSICAL SHEAR WAVE SPEED  (>0)                       C
C       CPOIS:  POISSON COEFFICIENT  (-1,0.5)                          C
C       RHO:    MASS DENSITY  (>0)                                     C
C       J:      INERTIA MOMENT DENSITY (ROTATIONAL INERTIA)  (>0)      C
C       XI:     DAMPING FACTOR  (>0)                                   C
C       GAMMA:  MICROPOLAR DISTORSION MODULI  (>0)                     C
C       MU_C:   RELATIVE SHEAR MODULI                                  C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: JUAN GOMEZ                                               C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 25 AUGUST 2012 BY NICOLAS GUARIN                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATCOS(DDSDDE,PROPS,NPROPS,NTENS,NDI)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DDSDDE(NTENS, NTENS),PROPS(NPROPS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0,FOUR=4.D0,
     1           SIX=6.0D0,ENUMAX=0.4999D0)
C
C**********************************************************************
C          M I C R O P O L A R     E L A S T I C I T Y
C**********************************************************************
C
C     Real elastic properties
C
      EG2 = TWO*PROPS(1)*PROPS(1)*PROPS(3)
      EG = EG2/TWO
      ENU = PROPS(2)
      EMOD = EG2*(ONE+ENU)
      ELAM = (ENU*EMOD)/((ONE+ENU)*(ONE-TWO*ENU))
      ELPNU = ELAM+EG2
      XI = PROPS(5)
      GAMMA = PROPS(6)
      EGC = PROPS(7)
C
C     Elastic complex stiffness
C
      DO K1=1, 3
        DO K2=1, 3
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=ELPNU
      END DO
      DDSDDE(4, 4) = EG+EGC
      DDSDDE(4, 5) = EG-EGC
      DDSDDE(5, 4) = DDSDDE(4, 5)
      DDSDDE(5, 5) = EG+EGC
      DDSDDE(6, 6) = GAMMA
      DDSDDE(7, 7) = GAMMA
C
      RETURN
C
      END
C

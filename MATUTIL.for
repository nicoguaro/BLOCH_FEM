C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C--C O N S T I T U T I V E  M A T E R I A L S  S U B R O U T I N E S---C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     UNIVERSIDAD EAFIT                                                C
C     2011                                                             C
C     CANNOT BE USED FOR PLANE STRESS                                  C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C                                                                      C
C     PROPS(1) - E                                                     C
C     PROPS(2) - NU                                                    C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATELA(DDSDDE,NDI,NTENS,PROPS,NPROPS)
C
      IMPLICIT REAL*8(A-H,O-Z)
C 
      DIMENSION DDSDDE(NTENS, NTENS),PROPS(NPROPS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0)
C
C**********************************************************************
C     I S O T R O P I C   E L A S T I C I T Y
C     P L A N E  S T R A I N  A N A L Y S I S
C**********************************************************************
C
C     ELASTIC PROPERTIES
C
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
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
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UMATPL                                                C
C     UMAT_PCLI_R.for   Plasticity Classical Isotropic Return Mapping  C
C                       LINEAR ISTROPIC HARDENING                      C
C     UNIVERSIDAD EAFIT                                                C
C     2005                                                             C
C     METODOS NUMERICOS                                                C
C     SEM 01/2005                                                      C
C                                                                      C
C     ISOTROPIC HARDENING PLASYICITY                                   C
C     CANNOT BE USED FOR PLANE STRESS                                  C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C     JANUARY 17/2005                                                  C
C                                                                      C
C     LOCAL ARRAYS                                                     C
C                                                                      C
C     EELAS - ELASTIC STRAINS STATEV(1..NTENS)                         C
C     EPLAS - PLASTIC STRAINS STATEV(NTENS+1...2*NTENS)                C
C     FLOW  - PLASTIC FLOW DIRECTION                                   C
C     HARD  - HARDENING MODULUS                                        C
C                                                                      C
C     PROPS(1) - E                                                     C
C     PROPS(2) - NU                                                    C
C     PROPS(3..) - SYIELD AN HARDENING DATA                            C
C     CALLS KUHARD FOR CURVE OF YIELD STRESS VS. PLASTIC STRAINS       C
C     EQPLAS - EQUIVALENT PLASTIC STRAIN STATEV(2*NTENS+1)             C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATPL(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1                  DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME,
     2                  DTIME, TEMP,DTEMP, PREDEF, DPRED, CMNAME, NDI,
     3                  NSHR, NTENS, NSTATV,PROPS, NPROPS, COORDS, DROT,
     4                  PNEWDT, CELENT, DFGRD0,DFGRD1, NOEL, NPT, LAYER,
     5                  KSPT, KSTEP, KINC)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*80 CMNAME
C 
      DIMENSION STRESS(NTENS), STATEV(NSTATV),DDSDDE(NTENS, NTENS),
     1DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),
     2PREDEF(1), DPRED(1), PROPS(NPROPS),COORDS(3),DROT(3, 3),
     3DFGRD0(3, 3), DFGRD1(3, 3)
C
      DIMENSION EELAS(NTENS), EPLAS(NTENS), FLOW(NTENS),SDEV(NTENS,1),
     1P(NTENS,NTENS)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,MAXITER=10)
C
C**********************************************************************
C     I S O T R O P I C  M I S E S  E L A S T O P L A S T I C I T Y
C     P L A N E  S T R A I N  A N A L Y S I S
C**********************************************************************
C
      CALL CLEAR(SDEV,NTENS,1)
C
C     RECOVER EQUIVALENT PLASTIC STRAIN, ELASTIC STRAINS, AND PLASTIC
C     STRAINS. ALSO INITIALIZE USER-DEFINED DATA SETS     
C
      EQPLAS=STATEV(1+2*NTENS)
      DO K1=1, NTENS
        EELAS(K1)=STATEV(K1)
        EPLAS(K1)=STATEV(K1+NTENS)
        FLOW(K1)=ZERO
      END DO
      SMISES=STATEV(2+2*NTENS)
C
C     ELASTIC PROPERTIES
C
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
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
C     CALCULATE PREDICTOR STRESSES AND ELASTIC STRAIN
C
      DO K1=1, NTENS
        DO K2=1, NTENS
          STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
        EELAS(K1)=EELAS(K1)+DSTRAN(K1)
      END DO
C
C     CALCULATE EQUIVALENT VON MISES STRESS
C
      SHYDRO=(STRESS(1)+STRESS(2)+STRESS(3))/THREE
      CALL PROYECTOR(P)
      CALL MMULT(P,NTENS,NTENS,STRESS,NTENS,1,SDEV)
C
      SMISES=(SDEV(1,1)-SDEV(2,1))**2
     1+(SDEV(2,1)-SDEV(3,1))**2+(SDEV(3,1)-SDEV(1,1))**2
      DO K1=NDI+1, NTENS
        SMISES=SMISES+SIX*SDEV(K1,1)**2
      END DO
      SMISES=SQRT(SMISES/TWO)
      FBAR=DSQRT(TWO/THREE)*SMISES
C
C     GET YIELD STRESS FROM THE SPECIFIED HARDENING CURVE
C
      NVALUE=NPROPS/2-1
      CALL UHARD(SYIEL0,HARD,EQPLAS,NVALUE,PROPS(3))
C
C     DETERMINE IF ACTIVELY YIELDING
C
      SYIELD=SYIEL0
      IF(FBAR.GT.(ONE+TOLER)*SYIEL0) THEN
C
C       Actively yielding==>Compute consistncy parameter and equivalent
C       plastic strain.
C
        DO K1=1, NTENS
          FLOW(K1)=SDEV(K1,1)/FBAR
        END DO
C
        GAM_PAR=(FBAR-SYIEL0)/((HARD/EG3)+1)
        GAM_PAR=GAM_PAR/EG2
        EQPLAS=EQPLAS+DSQRT(TWO/THREE)*GAM_PAR
        CALL UHARD(SYIELD,HARD,EQPLAS,NVALUE,PROPS(3))
C
C       Update plastic strains, equivalent plastic
C       strains, stresses
C
C       Normal components
C
        DO K1=1,NDI
          EPLAS(K1)=EPLAS(K1)+GAM_PAR*FLOW(K1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
          STRESS(K1)=FLOW(K1)*SYIELD+SHYDRO
        END DO
C
C       Shear components
C
        DO K1=NDI+1,NTENS
          EPLAS(K1)=EPLAS(K1)+TWO*GAM_PAR*FLOW(K1)
          EELAS(K1)=EELAS(K1)-EPLAS(K1)
          STRESS(K1)=FLOW(K1)*SYIELD
        END DO
C
C       Formulate the Jacobian (Material tangent)
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
C
C       Calculate effective properties 
C
        BETA1=ONE-EG2*GAM_PAR/FBAR
        C1=ONE+(HARD/EG3)
        BETABAR=(ONE/C1)-(EG2*GAM_PAR/FBAR)
        EFFG=EG*BETA1
        EFFG2=TWO*EFFG
        EFFLAM=ELAM+(EG2/THREE)*(ONE-BETA1)
        EFFHRD=-EG2*BETABAR
C
        DO K1=1, NDI
          DO K2=1, NDI
            DDSDDE(K2, K1)=EFFLAM
          END DO
          DDSDDE(K1, K1)=(EFFG2+EFFLAM)
        END DO
        DO K1=NDI+1, NTENS
          DDSDDE(K1, K1)=EFFG
        END DO
        DO K1=1, NTENS
          DO K2=1, NTENS
            DDSDDE(K2, K1)=DDSDDE(K2, K1)+EFFHRD
     &                     *FLOW(K2)*FLOW(K1)
          END DO
        END DO
      END IF
C
C     STORE ELASTIC STRAINS, (EQUIVALENT) PLASTIC STRAINS
C     IN STATE VARIABLE ARRAY
C
      DO K1=1, NTENS
        STATEV(      K1)=EELAS(K1)
        STATEV(NTENS+K1)=EPLAS(K1)
      END DO
      STATEV(2*NTENS+1)=EQPLAS
      STATEV(2*NTENS+2)=SMISES
C
      RETURN
C
      END
C
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UHARD                                                 C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UHARD(SYIELD,EHARD,EQPLAS,NVALUE,TABLE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION TABLE(2,NVALUE)
C
      PARAMETER(ZERO=0.D0,TWO=2.D0,THREE=3.D0)
C
      SYIEL0=TABLE(1,1)
C
C     Compute hardening modulus
C
      EHARD=(TABLE(1,2)-TABLE(1,1))/TABLE(2,2)
C
C     Compute yield stress corresponding to EQPLAS
C
      SYIELD=DSQRT(TWO/THREE)*(SYIEL0+EHARD*EQPLAS)
C
      RETURN
C
      END
C
C1122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE PROYECTOR                                             C
C     PROJECTS THE STRESS TENSOR INTO THE DEVIATORIC SPACE             C
C     FOR A PLANE STRAIN PROBLEM                                       C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PROYECTOR(P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0)
C
      DIMENSION P(4,4)
C
      CALL CLEAR(P,4,4)
      P(1,1)=TWO/THREE
      P(1,2)=-ONE/THREE
      P(1,3)=-ONE/THREE
      P(2,1)=-ONE/THREE
      P(2,2)=TWO/THREE
      P(2,3)=-ONE/THREE
      P(3,1)=-ONE/THREE
      P(3,2)=-ONE/THREE
      P(3,3)=TWO/THREE
      P(4,4)=ONE
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C     SUBROUTINE UMATCO                                                C
C     UMATCO.for        Isotropic full cosserat                        C
C                                                                      C
C     UNIVERSIDAD EAFIT                                                C
C     OCTOBER 2010                                                     C
C                                                                      C
C     NTENS: LENGTH OF STRESS VECTOR                                   C
C     NDI:   NUMBER OF NORMAL STRESS COMPONENTS                        C
C                                                                      C
C     LOCAL ARRAYS                                                     C
C                                                                      C
C     PROPS(1) - E                                                     C
C     PROPS(2) - NU                                                    C
C     PROPS(3) - SYIELD AN HARDENING DATA                              C
C     PROPS(4) - BENDING MODULUS                                       C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE UMATCO(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1                  DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME,
     2                  DTIME, TEMP,DTEMP, PREDEF, DPRED, CMNAME, NDI,
     3                  NSHR, NTENS, NSTATV,PROPS, NPROPS, COORDS, DROT,
     4                  PNEWDT, CELENT, DFGRD0,DFGRD1, NOEL, NPT, LAYER,
     5                  KSPT, KSTEP, KINC)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*80 CMNAME
C 
      DIMENSION STRESS(NTENS), STATEV(NSTATV),DDSDDE(NTENS, NTENS),
     1          DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2          TIME(2),PREDEF(1), DPRED(1), PROPS(NPROPS),COORDS(3),
     3          DROT(3, 3),DFGRD0(3, 3), DFGRD1(3, 3)
C
      DIMENSION SDEV(NTENS,1)
C
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.0D0,
     1           FOUR=4.D0,ENUMAX=0.4999D0, NEWTON=10, TOLER=1.0D-7,
     2           MAXITER=10)
C
C**********************************************************************
C     ISOTROPIC FULL COSSERAT  ELASTICITY
C     PLANE STRAIN ANALYSIS
C**********************************************************************
C
C     ELASTIC PROPERTIES
C
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2),ENUMAX)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
      ENUC=PROPS(3)
      GAMMA=PROPS(4)
C
C     ELASTIC STIFFNESS-SYMMETRIC PART
C
      DO K1=1, NDI
        DO K2=1, NDI
          DDSDDE(K2, K1)=ELAM
        END DO
        DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1, NTENS-3
        DDSDDE(K1, K1)=EG
      END DO
C
      K1=NTENS-2
      DDSDDE(K1, K1)=ENUC
C
      DO K1=NTENS-1, NTENS
        DDSDDE(K1, K1)=FOUR*GAMMA
      END DO
C
C     UPDATE STRESSES
C
      DO K1=1, NTENS
        DO K2=1, NTENS
          STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
        END DO
      END DO
C
      RETURN
C
      END
C
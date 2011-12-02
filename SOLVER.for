C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE NEWTON                                                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE NEWTON(UH,UG,SVARSEGPT,SVARSGPTH,NUMNP,NUMEL,NUMAT,NNE,
     1                  NSVARS,NMNE,NDOFDIM,MXDOFEL,IELT,IELCON,NDOFN,
     2                  NDOFEL,MATP,NMATP,NMPR,AMATE,COORD,LM,NEQ,SKG,
     3                  IP,MAXA,RHSG,IOUT,KINC,DPT,RR,NINCR,IFLAG,
     4                  ITEROUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,MAXITER=10,TOL=1.0D-7)
C
C     REAL ARRAYS PASSED IN
C
      DIMENSION UH(NEQ,NINCR),UG(NEQ),SVARSEGPT(NUMEL,NSVARS),
     1          AMATE(NMPR,NUMAT),COORD(NDOFDIM,NUMNP),SKG(IP),
     2          RHSG(NEQ,1),DPT(NEQ),RR(NEQ),
     3          SVARSGPTH(NUMEL,NSVARS,NINCR)
C
C     INTEGER ARRAYS PASSED IN
C
      DIMENSION NNE(NUMEL),IELT(NUMEL),IELCON(NMNE,NUMEL),NDOFN(NUMNP),
     1          NDOFEL(NUMEL),MATP(NUMEL),NMATP(NUMAT),
     2          LM(MXDOFEL,NUMEL),MAXA(NEQ+1)
C
C     LOCAL ARRAYS
C
      DIMENSION DU(NEQ),DDU(NEQ)
C
      MNM=NEQ+1
      ITER=1
      IFLAG=0
C
      DO
C
C        SOLVES FOR THE CORRECTION DDU CONSISTENT
C        WITH THE LAST OBTAINED RESIDUAL VECTOR
C
        DO II=1,NEQ
           DDU(II)=RR(II)
        END DO
C
        CALL COLSOL(SKG,DDU,MAXA,NEQ,IP,MNM,1)
        CALL COLSOL(SKG,DDU,MAXA,NEQ,IP,MNM,2)
C
C       UPDATES TO CURRENT:
C       NODAL DISPLACEMENTS UG()
C
        CALL UPDVEC(UG,NEQ,DDU)
C
C       COMPUTES THE ACTUAL INCREMENT
C
        DO K1=1,NEQ
          DU(K1)=UG(K1)-UH(K1,KINC-1)
        END DO
C
C       UPDATES TO CURRENT:
C       STATE VARIABLES SVARSEGPT()
C       INTERNAL NODAL FORCES RHSG()
C
        CALL GSTFASEM(UG,DU,SVARSEGPT,NUMNP,NUMEL,NUMAT,NNE,NSVARS,
     1                NMNE,NDOFDIM,MXDOFEL,IELT,IELCON,NDOFN,NDOFEL,
     2                MATP,NMATP,NMPR,AMATE,COORD,LM,NEQ,SKG,IP,MAXA,
     3                RHSG,IOUT,KINC)
C
C       COMPUTES THE RESIDUAL AND ITS NORM
C
        DO II=1,NEQ
           RR(II)=DPT(II)+RHSG(II,1)
        END DO
C
        CALL SQRNORM(RR,NEQ,E2)
C
C       ACCEPTS OR REJECTS SOLUTION
C       BASED ON TOLERANCE CRITERIA.
C       IF REJECTS STRESS SOLUTION IS
C       RESET TO LAST CONVERGED VALUE.
C
        IF(E2.LT.TOL) THEN
          IFLAG=1
          ITEROUT=ITER
          GOTO 10
C
        ELSE
          IF(ITER.LT.MAXITER) THEN
            ITER=ITER+1
            CALL STRSUPD(SVARSEGPT,SVARSGPTH,MXNO,MXINC,MXEQ,MXEL,
     1                   MXSVARS,KINC,1)
          ELSE
            ITEROUT=MAXITER
            GOTO 10
          END IF
        END IF
      END DO
C
   10 CONTINUE
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SOLVERHS                                              C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SOLVERHS(RFAC,R,U,NEQ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION RFAC(NEQ,NEQ),R(NEQ),U(NEQ)
C
      CALL CLEARV(U,NEQ)
C
cc      CALL DLFSDS(NEQ,RFAC,NEQ,R,U)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE COLSOL(A,V,U,IMIS,IMAX,NEQ,NWK,KKK)                   C
C                 T                                                    C
C     Performs LDL  factorization of the vector form stored stiffness  C
C     matrix A or performs reduction and back substitution of the load C
C     vector V                                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    INPUT PARAMTERS                                                   C
C      A(NWK)      VECTOR FORM OF STIFFNESS MATRIX                     C
C      V(NEQ)      RHS VECTOR                                          C
C      IMIS(NEQ)   ARRAY STORING THE ROW NUMBER OF 1ST NONZERO         C
C                  ELEMENT IN COLUMN J                                 C
C      IMAX(NEQ)   ARRAY STORING THE LOCATION OF EACH DIAGONAL ELEMENT C
C                  OF Kij IN BIG VECTOR A                              C
C      KKK         MODE FLAG                                           C
C                  EQ.1 FACTORIZATION MODE                             C
C                  EQ.2 REDUCTION AND BACKSUBSTITUTION MODE            C
C                                                                      C
C   OUTPUT PARAMTERS                                                   C
C      A(NWK)      VECTOR FORM OF FACTORIZED STIFFNESS MATRIX          C
C      V(NEQ)      REDUCED RHS VECTOR                                  C
C      U(NEQ)      SOLUTION VECTOR                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE COLSOL(A,V,MAXA,NN,NWK,NNM,KKK)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.0D0)
C
      DIMENSION MAXA(NNM),A(NWK),V(NN)
C
C     Perform factorization
C
      IF (KKK-2) 40,150,150
   40 DO 140 N=1,NN
      KN=MAXA(N)
      KL=KN+1
      KU=MAXA(N+1)-1
      KH=KU-KL
      IF (KH) 110,90,50
   50 K=N-KH
      IC=0
      KLT=KU
      DO 80 J=1,KH
      IC=IC+1
      KLT=KLT-1
      KI=MAXA(K)
      ND=MAXA(K+1)-KI-1
      IF (ND) 80,80,60
   60 KK=MIN0(IC,ND)
      C=0.
      DO 70 L=1,KK
   70 C=C+A(KI+L)*A(KLT+L)
      A(KLT)=A(KLT)-C
   80 K=K+1
   90 K=N
      B=0.
      DO 100 KK=KL,KU
      K=K-1
      KI=MAXA(K)
      C=A(KK)/A(KI)
      B=B+C*A(KK)
  100 A(KK)=C
      A(KN)=A(KN)-B
  110 IF (A(KN)) 120,120,140
CCC  120 WRITE (IOUT,2000) N,A(KN)
  120 WRITE (*,*) N,A(KN)
      GOTO 800
C
  140 CONTINUE
      GOTO 900
C
C     Reduce RHS vector.
C
  150 DO 180 N=1,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF (KU-KL) 180,160,160
  160 K=N
      C=0.
      DO 170 KK=KL,KU
      K=K-1
  170 C=C+A(KK)*V(K)
      V(N)=V(N)-C
  180 CONTINUE
C
C     Back-substitute.
C
      DO 200 N=1,NN
      K=MAXA(N)
  200 V(N)=V(N)/A(K)
      IF (NN.EQ.1) GO TO 900
      N=NN
      DO 230 L=2,NN
      KL=MAXA(N)+1
      KU=MAXA(N+1)-1
      IF (KU-KL) 230,210,210
  210 K=N
      DO 220 KK=KL,KU
      K=K-1
  220 V(K)=V(K)-A(KK)*V(N)
  230 N=N-1
      GO TO 900
C
  800 STOP
  900 RETURN
C
 2000 CONTINUE
      END
C
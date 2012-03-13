CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C---------  BLOCH   U  T  I  L  I  T  I  E  S   B  L  O  C  K----------C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE BLOCH                                                 C
C                                                                      C
C     Read the couple of nodes with bloch-periodicity constraints      C
C                                                                      C
C                                                                      C
C     NCOND       :NUMBER OF NODAL CONSTRAINTS                         C
C     IMNODES     :IMAGE NODES                                         C
C     IRNODES     :REFERENCES NODES                                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 05 DECEMBER 2011                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE BLOCH(NDOF,SKG,SMG,NCOR,COORDS,NCOND_WO,NCOND,
     1           IRNODES_WO,IMNODES,IRNODES,NEVALS,AKXMIN,
     2           AKYMIN,AKXMAX,AKYMAX,NKX,NKY,EIGVALS)

      IMPLICIT REAL*8(A-H,O-Z)

      PARAMETER (PI = 3.141592653589793)


      DIMENSION SKG(NDOF,NDOF),SMG(NDOF,NDOF),COORDS(2,NCOR),
     1          IRNODES_WO(NCOND_WO),
     2          IMNODES(NCOND),IRNODES(NCOND)

C     LOCAL VARIABLES

      COMPLEX*16 SKAUX,SMAUX, SKBLOCH, SMBLOCH, FI, FCI, FR, FCR,
     1           XJ, RFACT
      DIMENSION SKAUX(NDOF,NDOF),SMAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     2          SMBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),IANODES(2*NCOND),
     3          EVALS(NEVALS),
     4          EIGVALS(NKX*NKY,NEVALS+2)

      XJ=(0.0,1.0)
      SKBLOCH = 0.0
      SMBLOCH = 0.0
      DO II=1,NCOND
        IANODES(2*II-1) = 2*IMNODES(II)-1
        IANODES(2*II) = 2*IMNODES(II)
      END DO


      DKX = (AKXMAX - AKXMIN)/(NKX-1)  ! X Wave Number step
      DKY = (AKYMAX - AKYMIN)/(NKY-1)  ! Y Wave Number step
      AKX = AKXMIN
      ICONT = 1
      DO IKX=1,NKX

       AKY = AKYMIN
        DO IKY=1,NKY
          SKAUX = SKG
          SMAUX = SMG

          WRITE(*,*)"kx= ", AKX, "ky= ", AKY

C          CALL BLOCH_BC(SKAUX,SMAUX,NDOF,COORDS,NCOND,NCOND_WO,  ! Bloch BC Imposition
C     1                  NCOR,AKX,AKY,IRNODES,IRNODES_WO,
C     2                  IMNODES,IANODES,SKBLOCH,SMBLOCH)


          CALL BLOCH_BC2(SKAUX,SMAUX,NDOF,COORDS,NCOND,NCOND_WO,  ! Bloch BC Imposition
     1                  NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                  IMNODES,IANODES,SKBLOCH,SMBLOCH)

CCCC      SOLUTION  CCCCC

C           CALL GEVALS(SKBLOCH,NDOF-2*NCOND,SMBLOCH,NDOF-2*NCOND,
C     1     NEVALS, EVALS)

          CALL GREVALS(SKBLOCH,NDOF-2*NCOND,SMBLOCH,NDOF-2*NCOND,
     1    NEVALS, EVALS)

C          CALL GEVALSN(SKBLOCH,NDOF-2*NCOND,SMBLOCH,NDOF-2*NCOND,
C     1                 NEVALS, EVALS)  ! Solver with no extras

C          CALL GEVALS(SKAUX,NDOF,SMAUX,NDOF,
C     1    NEVALS, EVALS)

C          CALL GREVALS(SKAUX,NDOF,SMAUX,NDOF,NEVALS,EVALS)

C          CALL GEVALSN(SKAUX,NDOF,SMAUX,NDOF,
C     1                 NEVALS, EVALS)  ! Solver with no extras


          EIGVALS(ICONT,1) = AKX
          EIGVALS(ICONT,2) = AKY
          EIGVALS(ICONT,3:NEVALS+2) = EVALS

          ICONT = ICONT+1

          AKY = AKY + DKY
        END DO

        AKX = AKX + DKX
      END DO

      RETURN

      END SUBROUTINE BLOCH

C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  BLOCH_BC()                                           C
C                                                                      C
C     Pending...                                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 11 MARCH 2012                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE BLOCH_BC(SKAUX,SMAUX,NDOF,COORDS,NCOND,NCOND_WO,
     1                    NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                    IMNODES,IANODES,SKBLOCH,SMBLOCH)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 SKAUX, SMAUX, SKBLOCH, SMBLOCH,FI, FCI, FR,
     1           FCR, XJ, BLOCHAUX

      DIMENSION SKAUX(NDOF,NDOF), SMAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     2          SMBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     3          IANODES(2*NCOND),COORDS(2,NCOR),
     4          IRNODES_WO(NCOND_WO),
     5          IMNODES(NCOND),IRNODES(NCOND),
     6          BLOCHAUX(NDOF-2*NCOND,NDOF-2*NCOND)

      XJ=(0.0,1.0)


      DO II=1,NCOND_WO     ! Assigning the Phase shifts for Reference Nodes
        IR = IRNODES_WO(II)
        XR = COORDS(1,IR); YR = COORDS(2,IR)
        FR = EXP(XJ*AKX*XR)*EXP(XJ*AKY*YR)
        FCR= EXP(-XJ*AKX*XR)*EXP(-XJ*AKY*YR)
        CALL CROWMULT(SKAUX,NDOF,SKAUX,2*IR-1,FCR) ! Stiffness
        CALL CROWMULT(SKAUX,NDOF,SKAUX,2*IR,FCR)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,2*IR-1,FR)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,2*IR,FR)

        CALL CROWMULT(SMAUX,NDOF,SMAUX,2*IR-1,FCR) ! Mass
        CALL CROWMULT(SMAUX,NDOF,SMAUX,2*IR,FCR)
        CALL CCOLMULT(SMAUX,NDOF,SMAUX,2*IR-1,FR)
        CALL CCOLMULT(SMAUX,NDOF,SMAUX,2*IR,FR)
      END DO

      DO II=1,NCOND     ! Assigning the Phase shifts for Image Nodes
        IM = IMNODES(II)
        XI = COORDS(1,IM); YI = COORDS(2,IM)
        FI = EXP(XJ*AKX*XI)*EXP(XJ*AKY*YI)
        FCI= EXP(-XJ*AKX*XI)*EXP(-XJ*AKY*YI)
        CALL CROWMULT(SKAUX,NDOF,SKAUX,2*IM-1,FCI) ! Stiffness
        CALL CROWMULT(SKAUX,NDOF,SKAUX,2*IM,FCI)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,2*IM-1,FI)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,2*IM,FI)

        CALL CROWMULT(SMAUX,NDOF,SMAUX,2*IM-1,FCI) ! Mass
        CALL CROWMULT(SMAUX,NDOF,SMAUX,2*IM,FCI)
        CALL CCOLMULT(SMAUX,NDOF,SMAUX,2*IM-1,FI)
        CALL CCOLMULT(SMAUX,NDOF,SMAUX,2*IM,FI)
      END DO


      DO II=1,NCOND            ! Summing rows and columns
        IR = IRNODES(II); IM = IMNODES(II)
        CALL CROWADD(SKAUX,NDOF,SKAUX,2*IM-1,2*IR-1) ! Stiffness
        CALL CROWADD(SKAUX,NDOF,SKAUX,2*IM,2*IR)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,2*IM-1,2*IR-1)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,2*IM,2*IR)

        CALL CROWADD(SMAUX,NDOF,SMAUX,2*IM-1,2*IR-1) ! Mass
        CALL CROWADD(SMAUX,NDOF,SMAUX,2*IM,2*IR)
        CALL CCOLADD(SMAUX,NDOF,SMAUX,2*IM-1,2*IR-1)
        CALL CCOLADD(SMAUX,NDOF,SMAUX,2*IM,2*IR)
      END DO

      CALL COLROWDEL(SKAUX,NDOF,NCOND,SKBLOCH,IANODES) ! Deleting redundant equations
      CALL COLROWDEL(SMAUX,NDOF,NCOND,SMBLOCH,IANODES)
     
C      The matrices SMBLOCH and SKBLOCH are slightly non-Hermitian
C      so, they are averaged with their conjugated matrices.
      CALL MHERMIT(SKBLOCH,NDOF-2*NCOND,NDOF-2*NCOND,
     1             BLOCHAUX)
      SKBLOCH = 1.0/2.0*(SKBLOCH+BLOCHAUX)
      CALL MHERMIT(SMBLOCH,NDOF-2*NCOND,NDOF-2*NCOND,
     1                 BLOCHAUX)
      SMBLOCH = 1.0/2.0*(SMBLOCH+BLOCHAUX)

      RETURN

      END SUBROUTINE BLOCH_BC
C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  BLOCH_BC2()                                          C
C                                                                      C
C     Pending...                                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 11 MARCH 2012                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE BLOCH_BC2(SKAUX,SMAUX,NDOF,COORDS,NCOND,NCOND_WO,
     1                    NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                    IMNODES,IANODES,SKBLOCH,SMBLOCH)


      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER CONT

      COMPLEX*16 SKAUX, SMAUX, SKBLOCH, SMBLOCH,FAC,
     1           FCR, XJ, BLOCHAUX,TMAT,THMAT, MAUX

      DIMENSION SKAUX(NDOF,NDOF), SMAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     2          SMBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     3          IANODES(2*NCOND),COORDS(2,NCOR),
     4          IRNODES_WO(NCOND_WO),
     5          IMNODES(NCOND),IRNODES(NCOND),
     6          BLOCHAUX(NDOF-2*NCOND,NDOF-2*NCOND),
     7          TMAT(NDOF,NDOF-2*NCOND),
     7          THMAT(NDOF-2*NCOND,NDOF),
     8          INDEXVEC(NDOF),INDEXREF(NDOF),
     9          MAUX(NDOF,NDOF-2*NCOND)

      XJ=(0.0,1.0)


      INDEXVEC = 0
      DO I=1,NCOND
        INDEXVEC(2*IMNODES(I)-1) = -(2*IRNODES(I)-1)
        INDEXVEC(2*IMNODES(I)) = -2*IRNODES(I)
      END DO
      INDEXREF = - INDEXVEC

      DO I=1,NDOF
        CONT = 0
        DO J=1,NCOND
          IF (INDEXVEC(I)<0)THEN
            IF( 2*IMNODES(J)<ABS(INDEXVEC(I)) )THEN
              CONT = CONT + 1
            END IF
          ELSEIF (2*IMNODES(J)<I)THEN
            CONT = CONT + 1
          END IF
        END DO
        IF(INDEXVEC(I)<0)THEN
          INDEXVEC(I) = INDEXVEC(I) + 2*CONT
        ELSE
          INDEXVEC(I) = I - 2*CONT
        END IF
      END DO

      TMAT = 0
      DO I=1,NDOF
        IF(INDEXVEC(I)>0)THEN
          J = INDEXVEC(I)
          TMAT(I,J) = 1
        ELSE
          J = ABS(INDEXVEC(I))
          IR = (INDEXREF(I)+1)/2; IM = (I+1)/2
          XI = COORDS(1,IM); YI = COORDS(2,IM)
          XR = COORDS(1,IR); YR = COORDS(2,IR)
          FAC = EXP(XJ*AKX*(XI-XR))*EXP(XJ*AKY*(YI-YR))
          TMAT(I,J) = FAC
        END IF
      END DO

      CALL MHERMIT(TMAT,NDOF,NDOF-2*NCOND,THMAT)
      CALL CMMULT(SKAUX,NDOF,NDOF,TMAT,NDOF,NDOF-2*NCOND,MAUX)
      CALL CMMULT(THMAT,NDOF-2*NCOND,NDOF,MAUX,NDOF,NDOF-2*NCOND,
     1            SKBLOCH)
      CALL CMMULT(SMAUX,NDOF,NDOF,TMAT,NDOF,NDOF-2*NCOND,MAUX)
      CALL CMMULT(THMAT,NDOF-2*NCOND,NDOF,MAUX,NDOF,NDOF-2*NCOND,
     1            SMBLOCH)


C      The matrices SMBLOCH and SKBLOCH are slightly non-Hermitian
C      so, they are averaged with their conjugated matrices.
      CALL MHERMIT(SKBLOCH,NDOF-2*NCOND,NDOF-2*NCOND,
     1             BLOCHAUX)
      SKBLOCH = 1.0/2.0*(SKBLOCH+BLOCHAUX)
      CALL MHERMIT(SMBLOCH,NDOF-2*NCOND,NDOF-2*NCOND,
     1                 BLOCHAUX)
      SMBLOCH = 1.0/2.0*(SMBLOCH+BLOCHAUX)
      
      RETURN

      END SUBROUTINE BLOCH_BC2
C   
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  GEVALS(A,NRA,B,NRB,NEVALS, EVALS)                    C
C                                                                      C
C     Solves the Complex Generalized EigenValue Problem                C
C       A*x = lambda*B*x                                               C
C     making a Cholesky factorization for the B matrix, and then       C
C     converting the original problem to a Standard Eigenvalue         C
C     Problem.                                                         C
C                                                                      C
C     Both matrices A and B should be Hermitian and additionally       C
C     Positive-Definiteness is needed for the B matrix.                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 05 DECEMBER 2011                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE GEVALS(A,NRA,B,NRB,NEVALS, EVALS)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 A, B, RFACT, RINV, RHINV, C,AUX

      DIMENSION A(NRA,NRA), B(NRB,NRB), RFACT(NRB,NRB), RINV(NRB,NRB),
     1          RHINV(NRB,NRB), C(NRB,NRB), AUX(NRB,NRB), EVALS(NEVALS)


      CALL DLFTDH(NRB,B,NRB,RFACT,NRB)  ! Cholesky Factorization of Mass Matrix

      CALL DLINCT(NRB,RFACT,NRB,2,RINV,NRB) ! Matrix R inversion

      CALL MHERMIT(RINV,NRA,NRA,RHINV)

      CALL CMMULT(A,NRA,NRA,RINV,NRA,NRA,AUX)

      CALL CMMULT(RHINV,NRA,NRA,AUX,NRA,NRA,C)

      CALL MHERMIT(C,NRA,NRA,AUX)

      C = 1.0/2.0*(C + AUX)

      CALL DEVAHF (NRB, NEVALS, C, NRB, .TRUE., EVALS)

      RETURN

      END SUBROUTINE GEVALS

C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  GEVALSN(A,NRA,B,NRB,NEVALS, EVALS)                   C
C                                                                      C
C     Solves the Complex Generalized EigenValue Problem                C
C       A*x = lambda*B*x                                               C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 07 MARCH 2012                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE GEVALSN(A,NRA,B,NRB,NEVALS, EVALS)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 A, B, ALPHA, BETA

      DIMENSION A(NRA,NRA), B(NRB,NRB), ALPHA(NRA), BETA(NRA),
     1          EVALS(NEVALS)


      CALL DGVLCG(NRA, A, NRA, B, NRB, ALPHA, BETA)

      DO I=1,NRA
        ALPHA(I) = ALPHA(I)/BETA(I)
      END DO

      CALL DSVRGN (NRA, ALPHA, ALPHA)  ! Sort the vector

      DO I=1,NEVALS
        EVALS(I) = ALPHA(I)
      END DO

      RETURN

      END SUBROUTINE GEVALSN



C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  GREVALS(A,NRA,B,NRB,NEVALS, EVALS)                   C
C                                                                      C
C     Solves the Complex Generalized EigenValue Problem                C
C       A*x = lambda*B*x                                               C
C     rewritting the problem as a Real Equivalent one                  C
C     doubling the number of DOF.                                      C
C                                                                      C
C     Both matrices A and B should be Hermitian and additionally       C
C     Positive-Definiteness is needed for the B matrix. The resulting  C
C     Real System have Symmetric matrices                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 05 DECEMBER 2011                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE GREVALS(A,NRA,B,NRB,NEVALS, EVALS)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 A, B, ALPHA, BETA

      DIMENSION A(NRA,NRA),B(NRB,NRB),RA(2*NRA,2*NRA),RB(2*NRB,2*NRB),
     1          EVALS(NEVALS),EVAUX(2*NEVALS), ALLEV(2*NRA),
     2          ALPHA(NRA), BETA(NRA)


C     Rewriting the complex A matrix as a real matrix
      RA(1:NRA,1:NRA) = REAL(A)
      RA(1:NRA,NRA+1:2*NRA) = -AIMAG(A)
      RA(NRA+1:2*NRA,1:NRA) = AIMAG(A)
      RA(NRA+1:2*NRA,NRA+1:2*NRA) = REAL(A)

C     Rewriting the complex B matrix as a real matrix
      RB(1:NRB,1:NRB) = REAL(B)
      RB(1:NRB,NRB+1:2*NRB) = -AIMAG(B)
      RB(NRB+1:2*NRB,1:NRB) = AIMAG(B)
      RB(NRB+1:2*NRB,NRB+1:2*NRB) = REAL(B)

      
      CALL DGVLSP(2*NRA, RA, 2*NRA, RB, 2*NRB, ALLEV)

      CALL DSVRGN (2*NRA, ALLEV, ALLEV)  ! Sort the vector

      DO I=1,NEVALS
        EVALS(I) = ALLEV(2*I-1)
      END DO



      RETURN 

      END SUBROUTINE GREVALS



C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE CMMULT                                                  C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE CMMULT(A,NRA,NCA,B,NRB,NCB,C)
C
      COMPLEX*16 A, B, C, DUM

      DIMENSION A(NRA,NCA),B(NRB,NCB),C(NRA,NCB)
C
      CALL CLEAR(C,NRA,NCB)
      DUM=0.0D0
      DO I=1,NRA
        DO J=1,NCB
          DO K=1,NCA
            DUM=DUM+A(I,K)*B(K,J)
          END DO
          C(I,J)=DUM
          DUM=0.0D0
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE CMMULT
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE MHERMIT                                                 C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE MHERMIT(A,NRA,NCA,B)
C
      COMPLEX*16 A, B
C
      DIMENSION A(NRA,NCA),B(NCA,NRA)
C
      CALL CLEAR(B,NCA,NRA)
      DO I=1,NRA
        DO J=1,NCA
          B(J,I)=CONJG(A(I,J))
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE MHERMIT
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE COLROWDEL(A,NA,NRED,B,COLS)                             C
C                                                                      C
C    Delete NRED rows and columns from the matrix A and returns the    C
C    submatrix B. COLS are the columns/rows to be deleted.             C
C    NA is the dimension of A.                                         C
C                                                                      C
C    Both matrices, A and B are complex.
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE COLROWDEL(A,NA,NRED,B,COLS)
C
      COMPLEX*16 A, B
      INTEGER COLS, COLN, CONT, CONTII, CONTJJ
      DIMENSION A(NA,NA),B(NA-NRED,NA-NRED), COLS(NRED), COLN(NA)

      COLN = 0

      DO II=1,NRED
        COLN(COLS(II)) = 1
      END DO

      CONTII = 1
      DO II=1,NA
        IF(COLN(II).EQ.0)THEN
          CONTJJ = 1
          DO JJ=1,NA
            IF(COLN(JJ).EQ.0)THEN
              B(CONTII,CONTJJ) = A(II,JJ)
              CONTJJ = CONTJJ + 1
            END IF
          END DO
         CONTII = CONTII +1
        END IF

      END DO

      RETURN

C
      END SUBROUTINE COLROWDEL
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE CCOLMULT                                                C
C                                                                      C
C    Multiply the COL column of the complex matrix A by a factor FACT, C
C    the result is stored in the complex matrix B.                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE CCOLMULT(A,N,B,COL,FACT)
C
      COMPLEX*16 A, B, FACT
      INTEGER COL
      DIMENSION A(N,N),B(N,N)

       DO II=1,N
          DO JJ=1,N
            IF (JJ.EQ.COL) THEN
              B(II,JJ) = A(II,JJ)*FACT
            ELSE
              B(II,JJ) = A(II,JJ)
            END IF
          END DO
      END DO

      RETURN

      END SUBROUTINE CCOLMULT
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE CROWMULT                                                C
C                                                                      C
C    Multiply the ROW row of the complex matrix A by a factor FACT,    C
C    the result is stored in the complex matrix B.                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE CROWMULT(A,N,B,ROW,FACT)
C
      COMPLEX*16 A, B, FACT
      INTEGER ROW
      DIMENSION A(N,N),B(N,N)

       DO II=1,N
          DO JJ=1,N
            IF (II.EQ.ROW) THEN
              B(II,JJ) = A(II,JJ)*FACT
            ELSE
              B(II,JJ) = A(II,JJ)
            END IF
          END DO
      END DO

      RETURN

      END SUBROUTINE CROWMULT
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE CROWADD                                                 C
C                                                                      C
C    PENDING...                                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE CROWADD(A,N,B,ROW1,ROW2)

      COMPLEX*16 A, B
      INTEGER ROW1, ROW2
      DIMENSION A(N,N),B(N,N)
 
      DO II=1,N
          DO JJ=1,N
            IF (II.EQ.ROW2) THEN
              B(II,JJ) = A(II,JJ) + A(ROW1,JJ)

            ELSE
              B(II,JJ) = A(II,JJ)
            END IF
          END DO
      END DO

      RETURN

      END SUBROUTINE CROWADD
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE CCOLADD                                                 C
C                                                                      C
C    PENDING...                                                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE CCOLADD(A,N,B,COL1,COL2)

      COMPLEX*16 A, B
      INTEGER COL1, COL2
      DIMENSION A(N,N),B(N,N)
 
      DO II=1,N
          DO JJ=1,N
            IF (JJ.EQ.COL2) THEN
              B(II,JJ) = A(II,JJ) + A(II,COL1)

            ELSE
              B(II,JJ) = A(II,JJ)
            END IF
          END DO
      END DO

      RETURN

      END SUBROUTINE CCOLADD
C
C
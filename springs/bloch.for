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
C     LAST MOD: 08 AUGUST 2012                                         C
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
     2          IMNODES(NCOND),IRNODES(NCOND),
     3          EIGVALS(NKX*NKY,NEVALS+2)

C     LOCAL VARIABLES
      COMPLEX*16 SAUX, SKBLOCH, SMBLOCH, XJ, WORK
      DIMENSION SAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF - NCOND,NDOF - NCOND),
     2          SMBLOCH(NDOF - NCOND,NDOF - NCOND),IANODES(NCOND),
     3          EVALS(NDOF - NCOND),
     4          RWORK(3*(NDOF - NCOND)-2), WORK(2*(NDOF - NCOND)-1)

      N = NDOF - NCOND

      XJ=(0.0,1.0)
      SKBLOCH = 0.0
      DO II=1,NCOND
        IANODES(II) = IMNODES(II)
      END DO


      DKX = (AKXMAX - AKXMIN)/(NKX-1)  ! X Wave Number step
      DKY = (AKYMAX - AKYMIN)/(NKY-1)  ! Y Wave Number step
      AKX = AKXMIN
      ICONT = 1
      DO IKX=1,NKX

       AKY = AKYMIN
        DO IKY=1,NKY

          SAUX = DCMPLX(SKG)
          CALL BLOCH_BC(SAUX,NDOF,COORDS,NCOND,NCOND_WO,  ! Bloch BC Imposition
     1                  NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                  IMNODES,IANODES,SKBLOCH)

          SAUX = DCMPLX(SMG)
          CALL BLOCH_BC(SAUX,NDOF,COORDS,NCOND,NCOND_WO,  ! Bloch BC Imposition
     1                  NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                  IMNODES,IANODES,SMBLOCH)

        CALL ZHEGV(1, 'N', 'U', N, SKBLOCH, N, SMBLOCH, N, EVALS, WORK,
     1             2*N-1, RWORK, INFO)

          EIGVALS(ICONT,1) = AKX
          EIGVALS(ICONT,2) = AKY
          EIGVALS(ICONT,3:NEVALS+2) = EVALS(1:NEVALS)

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
      SUBROUTINE BLOCH_BC(SKAUX,NDOF,COORDS,NCOND,NCOND_WO,
     1                    NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                    IMNODES,IANODES,SKBLOCH)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 SKAUX,SKBLOCH,FI, FCI, FR,
     1           FCR, XJ

      DIMENSION SKAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF - NCOND,NDOF - NCOND),
     2          IANODES(NCOND),COORDS(2,NCOR),
     3          IRNODES_WO(NCOND_WO),
     4          IMNODES(NCOND),IRNODES(NCOND)

      N = NDOF - NCOND

      XJ=(0.0,1.0)

      DO II=1,NCOND_WO     ! Assigning the Phase shifts for Reference Nodes
        IR = IRNODES_WO(II)
        XR = COORDS(1,IR); YR = COORDS(2,IR)
        FR = EXP(XJ*AKX*XR)*EXP(XJ*AKY*YR)
        FCR= EXP(-XJ*AKX*XR)*EXP(-XJ*AKY*YR)
        CALL CROWMULT(SKAUX,NDOF,SKAUX,IR,FCR)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,IR,FR)
        
      END DO

      DO II=1,NCOND     ! Assigning the Phase shifts for Image Nodes
        IM = IMNODES(II)
        XI = COORDS(1,IM); YI = COORDS(2,IM)
        FI = EXP(XJ*AKX*XI)*EXP(XJ*AKY*YI)
        FCI= EXP(-XJ*AKX*XI)*EXP(-XJ*AKY*YI)
        CALL CROWMULT(SKAUX,NDOF,SKAUX,IM,FCI)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,IM,FI)
      END DO

      DO II=1,NCOND            ! Summing rows and columns
        IR = IRNODES(II); IM = IMNODES(II)
        CALL CROWADD(SKAUX,NDOF,SKAUX,IM,IR)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,IM,IR)
      END DO

      CALL COLROWDEL(SKAUX,NDOF,NCOND,SKBLOCH,IANODES) ! Deleting redundant equations

      RETURN

      END SUBROUTINE BLOCH_BC



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
      INTEGER COLS, COLN, CONTII, CONTJJ
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

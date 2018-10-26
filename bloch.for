CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE BLOCH                                                 C
C                                                                      C
C     Read the couple of DOF with bloch-periodicity constraints        C
C                                                                      C
C                                                                      C
C     NCOND       :NUMBER OF CONSTRAINTS BETWEEN DOF                   C
C     IMDOF       :IMAGE DOF                                           C
C     IRDOF       :REFERENCES DOF                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 20 AUGUST 2012                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE BLOCH(NDOF,SKG,SMG,NCOR,COORDS,IDOFCOOR,NCON_WO,NCON,
     1           IRDOF_WO,IMDOF,IRDOF,NEVALS,AKXMIN,
     2           AKYMIN,AKXMAX,AKYMAX,NKX,NKY,EIGVALS)

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (PI = 3.141592653589793)

      DIMENSION SKG(NDOF,NDOF),SMG(NDOF,NDOF),COORDS(2,NCOR),
     1          IRDOF_WO(NCON_WO),
     2          IMDOF(NCON),IRDOF(NCON),
     3          EIGVALS(NKX*NKY,NEVALS+2),
     4          IDOFCOOR(NDOF)

C     LOCAL VARIABLES
      COMPLEX*16 SAUX, SKBLOCH, SMBLOCH, XJ, WORK
      DIMENSION SAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF - NCON,NDOF - NCON),
     2          SMBLOCH(NDOF - NCON,NDOF - NCON),
     3          EVALS(NDOF - NCON),
     4          RWORK(3*(NDOF - NCON)-2), WORK(2*(NDOF - NCON)-1)

      N = NDOF - NCON

      XJ=(0.0,1.0)
      SKBLOCH = 0.0

      DKX = (AKXMAX - AKXMIN)/(NKX-1)  ! X Wave Number step
      DKY = (AKYMAX - AKYMIN)/(NKY-1)  ! Y Wave Number step
      AKX = AKXMIN
      ICONT = 1
      DO IKX=1,NKX

        AKY = AKYMIN
        DO IKY=1,NKY


          SAUX = DCMPLX(SKG)
          CALL BLOCH_BC(SAUX,NDOF,COORDS,IDOFCOOR,NCON,NCON_WO,  ! Bloch BC Imposition
     1                  NCOR,AKX,AKY,IRDOF,IRDOF_WO,
     2                  IMDOF,SKBLOCH)

          SAUX = DCMPLX(SMG)
          CALL BLOCH_BC(SAUX,NDOF,COORDS,IDOFCOOR,NCON,NCON_WO,  ! Bloch BC Imposition
     1                  NCOR,AKX,AKY,IRDOF,IRDOF_WO,
     2                  IMDOF,SMBLOCH)

          CALL ZHEGV(1, 'N', 'L', N, SKBLOCH, N, SMBLOCH, N, EVALS,
     1               WORK, 2*N-1, RWORK, INFO)

          IF(INFO.EQ.0)THEN
            WRITE(*,777) ICONT, NKX*NKY
          ELSE
            WRITE(*,666) ICONT, INFO
          END IF

          EIGVALS(ICONT,1) = AKX
          EIGVALS(ICONT,2) = AKY
          EIGVALS(ICONT,3:NEVALS+2) = EVALS(1:NEVALS)

          ICONT = ICONT+1

          AKY = AKY + DKY
        END DO

        AKX = AKX + DKX
      END DO

      RETURN

 777  FORMAT(/,'WAVE NUMBER: ',I5,'/',I5)
 666  FORMAT(/,'ERROR AT WAVE NUMBER: ',I5,'  SOLVER RETURNED: ',I5)

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
      SUBROUTINE BLOCH_BC(SKAUX,NDOF,COORDS,IDOFCOOR,NCOND,NCOND_WO,
     1                    NCOR,AKX,AKY,IRDOF,IRDOF_WO,
     2                    IMDOF,SKBLOCH)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 SKAUX,SKBLOCH,FI, FCI, FR,
     1           FCR, XJ

      DIMENSION SKAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF - NCOND,NDOF - NCOND),
     2          COORDS(2,NCOR),
     3          IRDOF_WO(NCOND_WO),
     4          IMDOF(NCOND),IRDOF(NCOND),
     5          IDOFCOOR(NDOF)

      XJ=(0.0,1.0)

      DO II=1,NCOND_WO    ! Assigning the Phase shifts for Reference DOF
        IR = IRDOF_WO(II)
        IDCOR = IDOFCOOR(IR)
        XR = COORDS(1,IDCOR); YR = COORDS(2,IDCOR)
        FR = EXP(XJ*AKX*XR)*EXP(XJ*AKY*YR)
        FCR= EXP(-XJ*AKX*XR)*EXP(-XJ*AKY*YR)
        CALL CROWMULT(SKAUX,NDOF,SKAUX,IR,FCR)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,IR,FR)
      END DO

      DO II=1,NCOND    ! Assigning the Phase shifts for Image DOF
        IM = IMDOF(II)
        IDCOR = IDOFCOOR(IM)
        XI = COORDS(1,IDCOR); YI = COORDS(2,IDCOR)
        FI = EXP(XJ*AKX*XI)*EXP(XJ*AKY*YI)
        FCI= EXP(-XJ*AKX*XI)*EXP(-XJ*AKY*YI)
        CALL CROWMULT(SKAUX,NDOF,SKAUX,IM,FCI)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,IM,FI)
      END DO

      DO II=1,NCOND    ! Summing rows and columns
        IR = IRDOF(II); IM = IMDOF(II)
        CALL CROWADD(SKAUX,NDOF,SKAUX,IM,IR)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,IM,IR)
      END DO

      CALL COLROWDEL(SKAUX,NDOF,NCOND,SKBLOCH,IMDOF)    ! Deleting redundant equations

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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE NOD2DOF                                               C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE NOD2DOF(NUMN,MXDOFDIM,NEQ,ID,NCOND,NCOND_WO,NCON_DOF,
     1                  NCON_DOF_WO,IMNODES,IRNODES,IRNODES_WO,IMDOF,
     2                  IRDOF,IRDOF_WO,IDOFCOOR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ID(MXDOFDIM,NUMN),IMNODES(NCOND),IRNODES(NCOND),
     1          IRNODES_WO(NCOND_WO), IMDOF(NCON_DOF),IRDOF(NCON_DOF),
     2          IRDOF_WO(NCON_DOF_WO), IDOFCOOR(NEQ)

      ICONT = 1
      JCONT = 1
      DO J=1,NCOND
        DO  I=1,MXDOFDIM
            IF (ID(I,IMNODES(J)).NE.0)THEN
              IMDOF(ICONT) = ID(I,IMNODES(J))
              ICONT = ICONT + 1
            END IF
            IF (ID(I,IRNODES(J)).NE.0)THEN
              IRDOF(JCONT) = ID(I,IRNODES(J))
              JCONT = JCONT + 1
            END IF
        END DO
      END DO

      ICONT = 1
      DO J=1,NCOND_WO
        DO  I=1,MXDOFDIM
            IF (ID(I,IRNODES_WO(J)).NE.0)THEN
              IRDOF_WO(ICONT) = ID(I,IRNODES_WO(J))
              ICONT = ICONT + 1
            END IF
        END DO
      END DO

      DO J=1,NUMN
        DO  I=1,MXDOFDIM
              IDOFCOOR(ID(I,J)) = J
        END DO
      END DO

C
      RETURN
C
      END SUBROUTINE NOD2DOF
C

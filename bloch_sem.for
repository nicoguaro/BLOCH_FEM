CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C---------  BLOCH   U  T  I  L  I  T  I  E  S   B  L  O  C  K----------C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE BLOCH_SEM                                             C
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
      SUBROUTINE BLOCH_SEM(NDOF,SKG,SMG,NCOR,COORDS,NCOND_WO,NCOND,
     1           IRNODES_WO,IMNODES,IRNODES,NEVALS,AKXMIN,
     2           AKYMIN,AKXMAX,AKYMAX,NKX,NKY,EIGVALS)

      IMPLICIT REAL*8(A-H,O-Z)


      PARAMETER (PI = 3.141592653589793)
C      N = NDOF - 2*NCOND
      N = NDOF - NCOND


      DIMENSION SKG(NDOF,NDOF),SMG(NDOF,NDOF),COORDS(2,NCOR),
     1          IRNODES_WO(NCOND_WO),
     2          IMNODES(NCOND),IRNODES(NCOND)

C     LOCAL VARIABLES

      COMPLEX*16 SKAUX,SKBLOCH, FI, FCI, FR, FCR,
     1           XJ, RFACT, SKBLOCH2
      DIMENSION SKAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),IANODES(2*NCOND),
     2          EVALS(NDOF),
     3          EIGVALS(NKX*NKY,NEVALS+2),
     4          SKBLOCH2(NDOF-2*NCOND,NDOF-2*NCOND),
     5          RWORK(3*N-2), WORK(2*N-1)

      DO I=1,NDOF
          DO J=1,NDOF
              SKG(I,J) = SKG(I,J)/(SQRT(SMG(I,I))*SQRT(SMG(J,J)))
          END DO
      END DO

C      DEALLOCATE(SMG)

      XJ=(0.0,1.0)
      SKBLOCH = 0.0
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
          SKAUX = DCMPLX(SKG)

C          WRITE(*,100) ICONT, NKX*NKY, AKX, AKY
C100       FORMAT(I5,'/',I5,' kx= ',F6.4,' ky= ',F6.4)


          CALL BLOCH_SEM_BC(SKAUX,NDOF,COORDS,NCOND,NCOND_WO,  ! Bloch BC Imposition
     1                  NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                  IMNODES,IANODES,SKBLOCH)


C          CALL BLOCH_SEM_BC2(SKAUX,NDOF,COORDS,NCOND,NCOND_WO,  ! Bloch BC Imposition
C     1                  NCOR,AKX,AKY,IRNODES,IRNODES_WO,
C     2                  IMNODES,IANODES,SKBLOCH)


C          CALL BLOCH_SEM_BC3(SKAUX,NDOF,COORDS,NCOND,NCOR,AKX,AKY,  ! Bloch BC Imposition
C     1                     IRNODES,IMNODES,IANODES,SKBLOCH)



CCCC      SOLUTION  CCCCC

C           CALL DEVAHF(NDOF-2*NCOND,NEVALS,SKBLOCH,NDOF-2*NCOND,.TRUE.,
C     1                EVALS)

           CALL ZHEEV( 'N', 'U', N, SKBLOCH, N, EVALS, WORK, -1, RWORK,
     1                  INFO )



          EIGVALS(ICONT,1) = AKX
          EIGVALS(ICONT,2) = AKY
          EIGVALS(ICONT,3:NEVALS+2) = EVALS(1:NEVALS)

          ICONT = ICONT+1

          AKY = AKY + DKY
        END DO

        AKX = AKX + DKX
      END DO

      RETURN

      END SUBROUTINE BLOCH_SEM

C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  BLOCH_SEM_BC()                                       C
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
      SUBROUTINE BLOCH_SEM_BC(SKAUX,NDOF,COORDS,NCOND,NCOND_WO,
     1                    NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                    IMNODES,IANODES,SKBLOCH)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 SKAUX,SKBLOCH,FI, FCI, FR,
     1           FCR, XJ, BLOCHAUX

      DIMENSION SKAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     2          IANODES(2*NCOND),COORDS(2,NCOR),
     3          IRNODES_WO(NCOND_WO),
     4          IMNODES(NCOND),IRNODES(NCOND),
     5          BLOCHAUX(NDOF-2*NCOND,NDOF-2*NCOND)

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
      END DO


      DO II=1,NCOND            ! Summing rows and columns
        IR = IRNODES(II); IM = IMNODES(II)
        CALL CROWADD(SKAUX,NDOF,SKAUX,2*IM-1,2*IR-1) ! Stiffness
        CALL CROWADD(SKAUX,NDOF,SKAUX,2*IM,2*IR)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,2*IM-1,2*IR-1)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,2*IM,2*IR)

      END DO

      CALL COLROWDEL(SKAUX,NDOF,2*NCOND,SKBLOCH,IANODES) ! Deleting redundant equations
     
C      The matrices SMBLOCH and SKBLOCH are slightly non-Hermitian
C      so, they are averaged with their conjugated matrices.
      CALL MHERMIT(SKBLOCH,NDOF-2*NCOND,NDOF-2*NCOND,
     1             BLOCHAUX)
      SKBLOCH = 1.0/2.0*(SKBLOCH+BLOCHAUX)

      RETURN

      END SUBROUTINE BLOCH_SEM_BC
C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  BLOCH_SEM_BC2()                                      C
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
      SUBROUTINE BLOCH_SEM_BC2(SKAUX,NDOF,COORDS,NCOND,NCOND_WO,
     1                    NCOR,AKX,AKY,IRNODES,IRNODES_WO,
     2                    IMNODES,IANODES,SKBLOCH)


      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER CONT

      COMPLEX*16 SKAUX, SKBLOCH, FAC,
     1           FCR, XJ, BLOCHAUX,TMAT,THMAT, AMAUX

      DIMENSION SKAUX(NDOF,NDOF), SMAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     2          IANODES(2*NCOND),COORDS(2,NCOR),
     3          IRNODES_WO(NCOND_WO),
     4          IMNODES(NCOND),IRNODES(NCOND),
     5          BLOCHAUX(NDOF-2*NCOND,NDOF-2*NCOND),
     6          TMAT(NDOF,NDOF-2*NCOND),
     7          THMAT(NDOF-2*NCOND,NDOF),
     8          INDEXVEC(NDOF),INDEXREF(NDOF),
     9          AMAUX(NDOF,NDOF-2*NCOND)

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
      CALL CMMULT(SKAUX,NDOF,NDOF,TMAT,NDOF,NDOF-2*NCOND,AMAUX)
      CALL CMMULT(THMAT,NDOF-2*NCOND,NDOF,AMAUX,NDOF,NDOF-2*NCOND,
     1            SKBLOCH)


C      The matrices SMBLOCH and SKBLOCH are slightly non-Hermitian
C      so, they are averaged with their conjugated matrices.
      CALL MHERMIT(SKBLOCH,NDOF-2*NCOND,NDOF-2*NCOND,
     1             BLOCHAUX)
      SKBLOCH = 1.0/2.0*(SKBLOCH+BLOCHAUX)
      
      RETURN

      END SUBROUTINE BLOCH_SEM_BC2
C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE  BLOCH_SEM_BC3()                                      C
C                                                                      C
C     Pending...                                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA- UNIVERSIDAD EAFIT                    C
C     LAST MOD: 02 APRIL 2012                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE BLOCH_SEM_BC3(SKAUX,NDOF,COORDS,NCOND,NCOR,AKX,
     1                     AKY,IRNODES,IMNODES,IANODES,SKBLOCH)


      IMPLICIT REAL*8(A-H,O-Z)

      COMPLEX*16 SKAUX, SKBLOCH,FI, FCI, FR, XJ

      DIMENSION SKAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND,NDOF-2*NCOND),
     2          IANODES(2*NCOND),COORDS(2,NCOR),
     3          IMNODES(NCOND),IRNODES(NCOND)

      XJ=(0.0,1.0)

      DO II=1,NCOND     ! Assigning the Phase shifts
        IM = IMNODES(II)
        XI = COORDS(1,IM); YI = COORDS(2,IM)
        IR = IRNODES(II)
        XR = COORDS(1,IR); YR = COORDS(2,IR)
        FI = EXP(XJ*AKX*(XI-XR))*EXP(XJ*AKY*(YI-YR))
        FCI= EXP(-XJ*AKX*XI)*EXP(-XJ*AKY*YI)
        CALL CROWMULT(SKAUX,NDOF,SKAUX,2*IM-1,FCI) ! Stiffness
        CALL CROWMULT(SKAUX,NDOF,SKAUX,2*IM,FCI)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,2*IM-1,FI)
        CALL CCOLMULT(SKAUX,NDOF,SKAUX,2*IM,FI)

      END DO


      DO II=1,NCOND            ! Summing rows and columns
        IR = IRNODES(II); IM = IMNODES(II)
        CALL CROWADD(SKAUX,NDOF,SKAUX,2*IM-1,2*IR-1) ! Stiffness
        CALL CROWADD(SKAUX,NDOF,SKAUX,2*IM,2*IR)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,2*IM-1,2*IR-1)
        CALL CCOLADD(SKAUX,NDOF,SKAUX,2*IM,2*IR)
      END DO

      CALL COLROWDEL(SKAUX,NDOF,2*NCOND,SKBLOCH,IANODES) ! Deleting redundant equations
      RETURN

      END SUBROUTINE BLOCH_SEM_BC3

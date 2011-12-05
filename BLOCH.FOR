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
C     LAST MOD: 21 SEPTEMBER 2011                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE BLOCH(NDOF,SKG,SMG,NCOR,COORDS,NCOND_WO,NCOND,
     1           IMNODES_WO,IRNODES_WO,IMNODES,IRNODES,NEVALS,AKXMIN,
     2           AKYMIN,AKXMAX,AKYMAX,NKX,NKY,EIGVALS)

     IMPLICIT REAL*8(A-H,O-Z)

     PARAMETER (PI = 3.141592653589793)


     DIMENSION SKG(NDOF,NDOF),SMG(NDOF,NDOF),COORDS(2,NCOR),
     1          IMNODES_WO(NCOND_WO),IRNODES_WO(NCOND_WO),
     2          IMNODES(NCOND),IRNODES(NCOND)

C     LOCAL VARIABLES

      COMPLEX*16 SKAUX,SMAUX, SKBLOCH, SMBLOCH, FI, FCI, FR, FCR,CALPHA,
     1           XJ, RFACT, BLOCHAUX
      DIMENSION SKAUX(NDOF,NDOF),SMAUX(NDOF,NDOF),
     1          SKBLOCH(NDOF-2*NCOND_WO,NDOF-2*NCOND_WO),
     2          SMBLOCH(NDOF-2*NCOND_WO,NDOF-2*NCOND_WO),IANODES(NDOF),
     3          RBETA(NEVALS),CALPHA(NDOF-2*NCOND_WO),
     4          EIGVALS(NKX*NKY,NEVALS+2),
     5          RFACT(NDOF-2*NCOND_WO,NDOF-2*NCOND_WO),
     6          BLOCHAUX(NDOF-2*NCOND_WO,NDOF-2*NCOND_WO)

      XJ=(0.0,1.0)
      IANODES = 0
      SKBLOCH = 0.
      SMBLOCH = 0.
      DO II=1,NCOND
        IANODES(2*IMNODES(II)-1) = 1
        IANODES(2*IMNODES(II)) = 1
      END DO


      DKX = (AKXMAX - AXKMIN)/(NKX-1)  ! X Wave Number step
      DKY = (AKYMAX - AXYMIN)/(NKY-1)  ! Y Wave Number step
      AKX = AKXMIN
      ICONT = 1
      DO IKX=1,NKX

       AKY = AKYMIN
        DO IKY=1,NKY
          SKAUX = SKG
          SMAUX = SMG

          WRITE(*,*)"kx= ", AKX, "ky= ", AKY

          DO II=1,NCOND_WO         ! Assigning the Phase shifts
            XI = COORDS(1,IMNODES_WO(II)); YI = COORDS(2,IMNODES_WO(II))
            XR = COORDS(1,IRNODES_WO(II)); YR = COORDS(2,IRNODES_WO(II))
            FI = EXP(XJ*AKX*XI)*EXP(XJ*AKY*YI)
            FCI= EXP(-XJ*AKX*XI)*EXP(-XJ*AKY*YI)
            FR = EXP(XJ*AKX*XR)*EXP(XJ*AKY*YR)
            FCR= EXP(-XJ*AKX*XR)*EXP(-XJ*AKY*YR)
            SKAUX(2*IMNODES_WO(II)-1,:)= SKAUX(2*IMNODES_WO(II)-1,:)*FCI ! Stiffness
            SKAUX(2*IMNODES_WO(II),:)  = SKAUX(2*IMNODES_WO(II),:)*FCI
            SKAUX(:,2*IMNODES_WO(II)-1)= SKAUX(:,2*IMNODES_WO(II)-1)*FI
            SKAUX(:,2*IMNODES_WO(II))  = SKAUX(:,2*IMNODES_WO(II))*FI
            SKAUX(2*IRNODES_WO(II)-1,:)= SKAUX(2*IRNODES_WO(II)-1,:)*FCR
            SKAUX(2*IRNODES_WO(II),:)  = SKAUX(2*IRNODES_WO(II),:)*FCR
            SKAUX(:,2*IRNODES_WO(II)-1)= SKAUX(:,2*IRNODES_WO(II)-1)*FR
            SKAUX(:,2*IRNODES_WO(II))  = SKAUX(:,2*IRNODES_WO(II))*FR
      
            SMAUX(2*IMNODES_WO(II)-1,:)= SMAUX(2*IMNODES_WO(II)-1,:)*FCI  ! Mass
            SMAUX(2*IMNODES_WO(II),:)  = SMAUX(2*IMNODES_WO(II),:)*FCI
            SMAUX(:,2*IMNODES_WO(II)-1)= SMAUX(:,2*IMNODES_WO(II)-1)*FI
            SMAUX(:,2*IMNODES_WO(II))  = SMAUX(:,2*IMNODES_WO(II))*FI
            SMAUX(2*IRNODES_WO(II)-1,:)= SMAUX(2*IRNODES_WO(II)-1,:)*FCR
            SMAUX(2*IRNODES_WO(II),:)  = SMAUX(2*IRNODES_WO(II),:)*FCR
            SMAUX(:,2*IRNODES_WO(II)-1)= SMAUX(:,2*IRNODES_WO(II)-1)*FR
            SMAUX(:,2*IRNODES_WO(II))  = SMAUX(:,2*IRNODES_WO(II))*FR
          END DO


          DO II=1,NCOND            ! Summing rows and columns
            SKAUX(2*IRNODES(II)-1,:)=SKAUX(2*IRNODES(II)-1,:)  ! Stiffness
     1                                + SKAUX(2*IMNODES(II)-1,:)
            SKAUX(2*IRNODES(II),:)  =SKAUX(2*IRNODES(II),:)
     1                                + SKAUX(2*IMNODES(II),:)
            SKAUX(:,2*IRNODES(II)-1)=SKAUX(:,2*IRNODES(II)-1)
     1                                + SKAUX(:,2*IMNODES(II)-1)
            SKAUX(:,2*IRNODES(II))  =SKAUX(:,2*IRNODES(II))
     1                                + SKAUX(:,2*IMNODES(II))
     
            SMAUX(2*IRNODES(II)-1,:)=SMAUX(2*IRNODES(II)-1,:)  ! Mass
     1                                + SMAUX(2*IMNODES(II)-1,:)
            SMAUX(2*IRNODES(II),:)  =SMAUX(2*IRNODES(II),:)
     1                                + SMAUX(2*IMNODES(II),:)
            SMAUX(:,2*IRNODES(II)-1)=SMAUX(:,2*IRNODES(II)-1)
     1                                + SMAUX(:,2*IMNODES(II)-1)
            SMAUX(:,2*IRNODES(II))  =SMAUX(:,2*IRNODES(II))
     1                                + SMAUX(:,2*IMNODES(II))
          END DO

          IB = 1      ! Rearranging the Bloch matrices
          DO II=1,NDOF
            JB = 1
            IF(IANODES(II).NE.1)THEN
             DO JJ=1,NDOF
                IF(IANODES(JJ).NE.1)THEN
                  SKBLOCH(IB,JB) = SKAUX(II,JJ)
                  SMBLOCH(IB,JB) = SMAUX(II,JJ)
                  JB = JB+1
                END IF
              END DO
              IB = IB+1
            END IF
          END DO

     
C         The matrices SMBLOCH and SKBLOCH are slightly non-Hermitian
C         so, they are averaged with their conjugated matrices.
          CALL MHERMIT(SKBLOCH,NDOF-2*NCOND+2,NDOF-2*NCOND+2,BLOCHAUX) 
          SKBLOCH = 1.0/2.0*(SKBLOCH+BLOCHAUX)
          CALL MHERMIT(SMBLOCH,NDOF-2*NCOND+2,NDOF-2*NCOND+2,BLOCHAUX) 
          SMBLOCH = 1.0/2.0*(SMBLOCH+BLOCHAUX)


CCCC      SOLUTION  CCCCC

C          CALL GEVALS(SKBLOCH,NDOF-2*NCOND+2,SMBLOCH,NDOF-2*NCOND+2,
C    1     NEVALS, RBETA)


          CALL GREVALS(SKBLOCH,NDOF-2*NCOND+2,SMBLOCH,NDOF-2*NCOND+2,
     1    NEVALS, RBETA)

          EIGVALS(ICONT,1) = AKX;  EIGVALS(ICONT,2) = AKY
          DO II=1, NEVALS
            EIGVALS(ICONT,II+2) = RBETA(II)
          END DO


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
C     SUBROUTINE  GEVALS(A,NRA,B,NRB,NEVALS, EVALS)                    C
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
C     SUBROUTINE  GREVALS(A,NRA,B,NRB,NEVALS, EVALS)                   C
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

      COMPLEX*16 A, B

      DIMENSION A(NRA,NRA), B(NRB,NRB), RA(2*NRA,2*NRA), RB(2*NRB,2*NRB),
     1          EVALS(NEVALS),EVAUX(2*NEVALS), ALLEV(2*NRA)


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
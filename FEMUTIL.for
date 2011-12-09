CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C---------F  E  M   U  T  I  L  I  T  I  E  S   B  L  O  C  K----------C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE GSTFASEM                                              C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE GSTFASEM(NUMNP,NUMEL,NUMAT,NNE,NMNE,MXDOFDIM,MXDOFEL,
     1                    IELT,IELCON,NDOFN,NDOFEL,MATP,NMATP,NMPR,
     2                    AMATE,COORD,LM,NEQ,SKG,SMG,IOUT,KINC)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0)
C
C     INTEGER ARRAYS
C
      ALLOCATABLE SK(:,:),SM(:,:)
C
      DIMENSION NNE(NUMEL),IELT(NUMEL),IELCON(NMNE,NUMEL),MATP(NUMEL),
     1          LM(MXDOFEL,NUMEL),LML(MXDOFEL),NDOFN(NUMNP),
     2          NDOFEL(NUMEL),NMATP(NUMAT)
C
C     MODEL ARRAYS
C
      DIMENSION AMATE(NMPR,NUMAT),COORD(MXDOFDIM,NUMNP),
     1          ELCOOR(MXDOFDIM,NMNE),SKG(NEQ,NEQ),SMG(NEQ,NEQ)
C
C     LOCAL ARRAYS
C
      DIMENSION PROPS(NMPR)
C
      CALL CLEAR(SKG,NEQ,NEQ)
      CALL CLEAR(SMG,NEQ,NEQ)
      CALL CLEARIV(LML,MXDOFEL)
C
C     RETRIEVES ELEMENT MATERIAL PROPERTIES
C
      DO I=1,NUMEL
C
        NPROPS=NMATP(MATP(I))
        DO K1=1,NPROPS
          PROPS(K1)=AMATE(K1,MATP(I))
        END DO
        ND=NDOFEL(I)
        ALLOCATE(SK(ND,ND),SM(ND,ND))
        CALL CLEAR(SK,ND,ND)
        CALL CLEAR(SM,ND,ND)
C
C       RETRIEVES ELEMENT NODAL COORDINATES
C       AND LOCAL DME OPERATOR
C
        K3=0
        DO J=1,NNE(I)
          ELCOOR(1,J)=COORD(1,IELCON(J,I))
          ELCOOR(2,J)=COORD(2,IELCON(J,I))
          K1=NDOFN(IELCON(J,I))
          DO K2=1,K1
            LML(K3+K2)=LM(K3+K2,I)
          END DO
          K3=K3+K1
        END DO
C
C       CALLS EMBEDDED USER ELEMENT SUBROUTINE UEL ACCORDING TO IELT(I)
C               IELT()=1 8 NODED QUAD.
C               IELT()=2 9 NODED QUAD.
C               IELT()=3 6 NODED TRIA.
C               IELT()=4 4 NODED QUAD.
C               IELT()=5 3 NODED TRIA.
C
        SELECT CASE(IELT(I))
C
          CASE(1)
          CALL UEL_FEMB8(SK,SM,ND,PROPS,NPROPS,ELCOOR,MXDOFDIM,NNE(I),I)
          IFLAG=1
cc          CASE(2)
cc          CALL UEL_FEM9(RHS,S,SVARS,ND,NSVARS,PROPS,NPROPS,
cc     1                  ELCOOR,MXDOFDIM,NNE(I),UL,DUL,V,A,I)
cc          CASE(3)
cc          CALL UEL_FEM6(RHS,S,SVARS,ND,NSVARS,PROPS,NPROPS,
cc     1                  ELCOOR,MXDOFDIM,NNE(I),UL,DUL,V,A,I)
C
        END SELECT
C
C       ASSEMBLE GLOBAL STIFFNESS MATRIX AND GLOBAL RHS VECTOR
C
        DO II=1,ND
          KK=LML(II)
          IF(KK.NE.0) THEN
             DO JJ=1,ND
               LL=LML(JJ)
               IF(LL.NE.0) THEN
                  SKG(KK,LL)=SKG(KK,LL)+SK(II,JJ)
                  SMG(KK,LL)=SMG(KK,LL)+SM(II,JJ)
               END IF
             END DO
          END IF
        END DO
        DEALLOCATE(SK,SM)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE ASSEMLIS                                              C
C                                                                      C
C     NUMNP      :NUMBER OF NODAL POINTS                               C
C     NUMEL      :NUMBER OF ELEMENTS                                   C
C     MXNE       :MAXIMUM NUMBER OF NODES IN A GIVEN ELEMENT           C
C     NMDOFEL    :MAXIMUM NUMBER OF DOF ASSIGNED TO A GIVEN ELEMENT    C
C     MXDOFDIM   :PROBLEM DIMENSION (2D, 3D)                           C
C     NNE(:)     :NUMBER OF NODES ASSIGNED TO EACH ELEMENT             C
C     NDOFN(:)   :NODAL DOF PER NODE                                   C
C     IELCON(:,:):ELEMENT CONNECTIVITIES                               C
C     LM(:,:)    :DME OPERATOR                                         C
C     ID(:,:)    :EQUATION NUMBERS ASSIGNED TO EACH NODE              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE ASSEMLIS(NUMNP,NUMEL,MXNE,MXDOFEL,MXDOFDIM,NNE,NDOFN,
     1                    IELCON,LM,ID,IIN,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION NNE(NUMEL),IELCON(MXNE,NUMEL),LM(MXDOFDIM*MXNE,NUMEL),
     1          ID(MXDOFDIM,NUMNP),NDOFN(NUMNP)
C
      DO I=1,NUMEL
        K3=0
        DO J=1,NNE(I)
          K2I=NDOFN(IELCON(J,I))
          DO K2=1,K2I
            LM(K3+K2,I)=ID(K2,IELCON(J,I))
          END DO
          K3=K3+K2I
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE FINDMIS                                               C
C     FINDS THE ROW NUMBER OF THE 1ST NON ZERO ELEMENT IN EACH COLUMN  C
C     AND THE COLUMN HEIGHTS                                           C
C                                                                      C
C     LM(:,:)    :DME OPERATOR                                         C
C     NUMEL      :NUMBER OF ELEMENTS                                   C
C     MXNE       :MAXIMUM NUMBER OF NODES IN A GIVEN ELEMENT           C
C     IMIS(:)         :                                                C
C     ICHS(:)    :COLUMN HEIGHTS                                       C
C     MAXA(:)    :ADDRESS OF DIAGONAL ELEMENTS INDEX                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE FINDMIS(LM,NUMEL,MXNE,NEQ,IMIS,ICHS,MK,MAXA)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION LM(2*MXNE,NUMEL),IMIS(NEQ),ICHS(NEQ),LIST(2*MXNE),
     1          ISL(2*MXNE),MAXA(NEQ+1)
C
      DO I=1,NEQ
        CALL CLEARIV(LIST,2*MXNE)
        CALL CLEARIV(ISL,2*MXNE)
        KK=0
        DO J=1,NUMEL
          DO K=1,2*MXNE
             LIST(K)=LM(K,J)
          END DO
          CALL SEARCH(LIST,2*MXNE,I,IFLAG)
          IF (IFLAG.EQ.1) THEN
            KK=KK+1
            CALL FINDSMALL(LIST,2*MXNE,IS)
            ISL(KK)=IS
          END IF
        END DO
        CALL FINDSMALL(ISL,2*MXNE,IS)
        IMIS(I)=IS
        ICHS(I)=I-IMIS(I)
      END DO
C
      CALL FINDLARGE(ICHS,NEQ,MK)
C
      MAXA(1)=1
      MAXA(2)=2
      DO I=3,NEQ
        ICONT=0
        DO K=1,I-1
          ICONT=ICONT+ICHS(K)
        END DO
        MAXA(I)=ICONT+I
      END DO
      MAXA(NEQ+1)=MAXA(NEQ)+ICHS(NEQ)+1
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SEARCH                                                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SEARCH(LIST,NDAT,NVAL,IFLAG)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION LIST(NDAT)
C
      IFLAG=0
      DO I=1,NDAT
         IF(LIST(I).EQ.NVAL) THEN
            IFLAG=1
         END IF
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE FINDSMALL                                             C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE FINDSMALL(LIST,NDAT,IS)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION LIST(NDAT)
C
      IFLAG=0
      I=1
      DO WHILE (IFLAG.EQ.0.OR.I.LT.NDAT)
         IF(LIST(I).NE.0) THEN
           IS=LIST(I)
           IFLAG=1
         END IF
         I=I+1
      END DO
C
      DO I=1,NDAT
         IF(LIST(I).GT.0.AND.LIST(I).LT.IS) THEN
            IS=LIST(I)
         END IF
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE FINDLARGE                                             C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE FINDLARGE(ICHS,NEQ,MK)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ICHS(NEQ)
C
      MK=ICHS(1)
C
      DO I=2,NEQ
         IF(ICHS(I).GT.MK) THEN
            MK=ICHS(I)
         END IF
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C         --G E N E R A L  F E M  S U B R O U T I N E S--              C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C STRAIN DISPLACEMENT INTERPOLATORS                                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM9:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 9-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM9(IDEL,NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER9(IELT,NDOFEL,NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM8:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 8-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM8(IDEL,NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER8(IELT,NDOFEL,NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM6:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 6-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM6(IDEL,NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER6(IELT,NDOFEL,NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM4:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 4-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM4(IDEL,NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER4(IELT,NDOFEL,NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END SUBROUTINE STDM4
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM3:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 3-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM3(IDEL,NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER3(IELT,NDOFEL,NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END SUBROUTINE STDM3
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SHAPE FUNCTIONS AND THEIR DERIVATIVES                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER9:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER9(IELT,NDOFEL,NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     9-NODED ELEMENT
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C
      P(1,9)=-TWO*R*SMS
      P(1,8)=-HALF*SMS-HALF*P(1,9)
      P(1,7)=-R*SP-HALF*P(1,9)
      P(1,6)= HALF*SMS-HALF*P(1,9)
      P(1,5)=-R*SM-HALF*P(1,9)
      P(1,4)=-QUARTER*SP-HALF*P(1,7)-HALF*P(1,8)-QUARTER*P(1,9)
      P(1,3)= QUARTER*SP-HALF*P(1,7)-HALF*P(1,6)-QUARTER*P(1,9)
      P(1,2)= QUARTER*SM-HALF*P(1,5)-HALF*P(1,6)-QUARTER*P(1,9)
      P(1,1)=-QUARTER*SM-HALF*P(1,8)-HALF*P(1,5)-QUARTER*P(1,9)
C
C        w.r.t.s
C
      P(2,9)=-TWO*S*RMS
      P(2,8)=-S*RM-HALF*P(2,9)
      P(2,7)= HALF*RMS-HALF*P(2,9)
      P(2,6)=-S*RP-HALF*P(2,9)
      P(2,5)=-HALF*RMS-HALF*P(2,9)
      P(2,4)= QUARTER*RM-HALF*P(2,7)-HALF*P(2,8)-QUARTER*P(2,9)
      P(2,3)= QUARTER*RP-HALF*P(2,7)-HALF*P(2,6)-QUARTER*P(2,9)
      P(2,2)=-QUARTER*RP-HALF*P(2,5)-HALF*P(2,6)-QUARTER*P(2,9)
      P(2,1)=-QUARTER*RM-HALF*P(2,8)-HALF*P(2,5)-QUARTER*P(2,9)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER8:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER8(IELT,NDOFEL,NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     8-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C      
      P(1,6)= HALF*SMS
      P(1,5)=-R*SM
      P(1,8)=-HALF*SMS
      P(1,7)=-R*SP
      P(1,2)= QUARTER*SM-HALF*P(1,5)-HALF*P(1,6)
      P(1,1)=-QUARTER*SM-HALF*P(1,8)-HALF*P(1,5)
      P(1,4)=-QUARTER*SP-HALF*P(1,7)-HALF*P(1,8)
      P(1,3)= QUARTER*SP-HALF*P(1,7)-HALF*P(1,6)
C
C     w.r.t.s
C
      P(2,6)=-S*RP
      P(2,5)=-HALF*RMS
      P(2,8)=-S*RM
      P(2,7)= HALF*RMS
      P(2,2)=-QUARTER*RP-HALF*P(2,5)-HALF*P(2,6)
      P(2,1)=-QUARTER*RM-HALF*P(2,8)-HALF*P(2,5)
      P(2,4)= QUARTER*RM-HALF*P(2,7)-HALF*P(2,8)
      P(2,3)= QUARTER*RP-HALF*P(2,7)-HALF*P(2,6)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER6:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER6(IELT,NDOFEL,NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0,THREE=3.0D0, HALF=0.5D0,
     1          QUARTER=0.25D0,FOUR=4.D0,EIGHT=8.0D0, ZERO=0.0D0)
C
      DIMENSION P(2,NNE)
C
C     6-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C
      P(1,6)=-FOUR*S
      P(1,5)= FOUR*S
      P(1,4)= FOUR-EIGHT*R-FOUR*S
      P(1,3)= ZERO
      P(1,2)= FOUR*R-ONE
      P(1,1)=-THREE+FOUR*R+FOUR*S
C
C     w.r.t.s
C
      P(2,6)= FOUR-FOUR*R-EIGHT*S
      P(2,5)= FOUR*R
      P(2,4)=-FOUR*R
      P(2,3)= FOUR*S-ONE
      P(2,2)= ZERO
      P(2,1)=-THREE+FOUR*R+FOUR*S
C
      RETURN
C
      END SUBROUTINE SFDER6
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER4:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER4(IELT,NDOFEL,NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     4-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C      
      P(1,1)=-QUARTER*SM
      P(1,2)= QUARTER*SM
      P(1,3)= QUARTER*SP
      P(1,4)=-QUARTER*SP
C
C     w.r.t.s
C

      P(2,1)=-QUARTER*RM
      P(2,2)=-QUARTER*RP
      P(2,3)= QUARTER*RP
      P(2,4)= QUARTER*RM
C
      RETURN
C
      END SUBROUTINE SFDER4
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER3:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER3(IELT,NDOFEL,NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.0D0,ZERO=0.0D0)
C
      DIMENSION P(2,NNE)
C
C
C     3-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C      
      P(1,1)= -ONE
      P(1,2)= ONE
      P(1,3)= ZERO
C
C     w.r.t.s
C
      P(2,1)= -ONE
      P(2,2)= ZERO
      P(2,3)= ONE
C
      RETURN
C
      END SUBROUTINE SFDER3
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF9(SN,R,S)                                      C
C      9-noded element shape function at point R,S                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF9(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0,QUART=0.25D0)
C
      DIMENSION SN(9)
C
      CALL CLEARV(SN,9)
C
      RP =ONE+R
      RM =ONE-R
      RMS=ONE-R*R
      SP =ONE+S
      SM =ONE-S
      SMS=ONE-S*S
C
      SN(9)=RMS*SMS
      SN(8)=HALF*SMS*RM-HALF*SN(9)
      SN(7)=HALF*RMS*SP-HALF*SN(9)
      SN(6)=HALF*SMS*RP-HALF*SN(9)
      SN(5)=HALF*RMS*SM-HALF*SN(9)
      SN(1)=QUART*RM*SM-HALF*SN(8)-HALF*SN(5)-QUART*SN(9)
      SN(2)=QUART*RP*SM-HALF*SN(6)-HALF*SN(5)-QUART*SN(9)
      SN(3)=QUART*RP*SP-HALF*SN(6)-HALF*SN(7)-QUART*SN(9)
      SN(4)=QUART*RM*SP-HALF*SN(8)-HALF*SN(7)-QUART*SN(9)
C
      RETURN
C
      END SUBROUTINE SHAPEF9
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF8(SN,R,S)                                      C
C      8-noded element shape function at point R,S                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF8(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0,QUART=0.25D0)
C
      DIMENSION SN(8)
C
      CALL CLEARV(SN,8)
C
      RP =ONE+R
      RM =ONE-R
      RMS=ONE-R*R
      SP =ONE+S
      SM =ONE-S
      SMS=ONE-S*S
C
      SN(8)=HALF*SMS*RM
      SN(7)=HALF*RMS*SP
      SN(6)=HALF*SMS*RP
      SN(5)=HALF*RMS*SM
      SN(1)=QUART*RM*SM-HALF*SN(8)-HALF*SN(5)
      SN(2)=QUART*RP*SM-HALF*SN(6)-HALF*SN(5)
      SN(3)=QUART*RP*SP-HALF*SN(6)-HALF*SN(7)
      SN(4)=QUART*RM*SP-HALF*SN(8)-HALF*SN(7)
C
      RETURN
C
      END SUBROUTINE SHAPEF8
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF6(SN,R,S)                                      C
C      6-noded element shape function at point R,S                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF6(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0,FOUR=4.0D0)
C
      DIMENSION SN(6)
C
      CALL CLEARV(SN,6)

      MRS = ONE - R -S
C
      SN(4) = FOUR*R*MRS
      SN(5) = FOUR*R*S
      SN(6) = 4*S*MRS
      SN(1) = MRS - HALF*SN(4) - HALF*SN(6)
      SN(2) = R - HALF*SN(4) - HALF*SN(5)
      SN(3) = S - HALF*SN(5) - HALF*SN(6)
C
      RETURN
C
      END SUBROUTINE SHAPEF6
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF4(SN,R,S)                                      C
C      4-noded element shape function at point R,S                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF4(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0,QUART=0.25D0)
C
      DIMENSION SN(4)
C
      CALL CLEARV(SN,4)
C
      RP =ONE+R
      RM =ONE-R
      RMS=ONE-R*R
      SP =ONE+S
      SM =ONE-S
      SMS=ONE-S*S
C
      SN(1)=QUART*RM*SM
      SN(2)=QUART*RP*SM
      SN(3)=QUART*RP*SP
      SN(4)=QUART*RM*SP
C
      RETURN
C
      END SUBROUTINE SHAPEF4
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF3(SN,R,S)                                      C
C      3-noded element shape function at point R,S                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF3(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0)
C
      DIMENSION SN(3)
C
      CALL CLEARV(SN,3)
C
      SN(1)=ONE - R - S
      SN(2)=R
      SN(3)=S
C
      RETURN
C
      END SUBROUTINE SHAPEF3
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KJACOPER                                                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE JACOPER(NNE,XJA,XCORD,RDER)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
C
      DIMENSION XJA(2,2),XCORD(2,NNE),RDER(2,NNE)
C
      CALL CLEAR(XJA,2,2)
C
      DUM=ZERO
      DO K=1,2
        DO J=1,2
          DO I=1,NNE
            DUM=DUM+RDER(K,I)*XCORD(J,I)
          END DO
          XJA(K,J)=DUM
          DUM=ZERO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE JACOPER
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      2X2 Jacobian inverse                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE JACINVE(XJA,DD,XJAI)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION XJA(2,2),XJAI(2,2),XADJ(2,2)
C
      XADJ(1,1)=XJA(2,2)
      XADJ(1,2)=XJA(1,2)
      XADJ(2,1)=XJA(2,1)
      XADJ(2,2)=XJA(1,1)
      DO J=1,2
        DO K=1,2
          COFA=((-1)**(J+K))*XADJ(J,K)
          XJAI(J,K)=COFA/DD
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE JACINVE
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE AMASS8(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 8-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE AMASS8(AMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(8),
     1          AUX1(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,8)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF8(SN,RII,SII)
C
      DO I=1,8
        KK=2*I-1
        SF(1,  KK)=SN(I)
        SF(2,KK+1)=SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,AUX1)
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMASS(I,J)=AUX1(I,J)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE AMASS8
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE AMASS9(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 9-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE AMASS9(AMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(9),
     1          AUX1(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,9)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF9(SN,RII,SII)
C
      DO I=1,9
        KK=2*I-1
        SF(1,  KK)=SN(I)
        SF(2,KK+1)=SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,AUX1)
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMASS(I,J)=AUX1(I,J)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE AMASS9
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE AMASS4(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 4-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE AMASS4(AMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(4),
     1          AUX1(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,4)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF4(SN,RII,SII)
C
      DO I=1,4
        KK=2*I-1
        SF(1,  KK)=SN(I)
        SF(2,KK+1)=SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,AUX1)
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMASS(I,J)=AUX1(I,J)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE AMASS4
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE AMASS3(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 3-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE AMASS3(AMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(3),
     1          AUX1(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,3)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF3(SN,RII,SII)
C
      DO I=1,3
        KK=2*I-1
        SF(1,  KK)=SN(I)
        SF(2,KK+1)=SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,AUX1)
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMASS(I,J)=AUX1(I,J)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE AMASS3
C
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
      SUBROUTINE GSTFASEM(NUMNP,NUMEL,NUMAT,NNE,NMNE,MXDOFEL,IELT,
     1                    IELCON,NDOFN,NDOFEL,MATP,NMATP,NMPR,
     2                    MCRD,AMATE,COORD,LM,NEQ,SKG,SMG,IATYPE,IOUT)

C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      ALLOCATABLE SK(:,:),SM(:,:),PROPS(:)
C
      PARAMETER (ZERO=0.D0)
C
C     INTEGER ARRAYS
C
      DIMENSION NNE(NUMEL),IELT(NUMEL),IELCON(NMNE,NUMEL),MATP(NUMEL),
     1          LM(NUMEL,MXDOFEL),LML(MXDOFEL),NDOFN(NUMNP),
     2          NDOFEL(NUMEL),NMATP(NUMAT),AMATE(NMPR,NUMAT),
     3          COORD(MCRD,NUMNP),ELCOOR(MCRD,NMNE),SKG(NEQ,NEQ),
     4          SMG(NEQ,NEQ)
C
      CALL CLEAR(SKG,NEQ,NEQ)
      CALL CLEAR(SMG,NEQ,NEQ)
      CALL CLEARIV(LML,MXDOFEL)
C
C     RETRIEVES ELEMENT MATERIAL PROPERTIES
C
      DO I=1,NUMEL
C
        ND=NDOFEL(I)
        NPROPS=NMATP(MATP(I))
        ALLOCATE(PROPS(NPROPS),SK(ND,ND),SM(ND,ND))
        CALL CLEAR(SK,ND,ND)
        CALL CLEAR(SM,ND,ND)
C
        DO K1=1,NPROPS
          PROPS(K1)=AMATE(K1,MATP(I))
        END DO
C
C       RETRIEVES ELEMENT NODAL COORDINATES
C       ASSEMBLY LIST
C
        K3=0
        DO J=1,NNE(I)
          K1=NDOFN(IELCON(J,I))
          DO JJ=1,MCRD
            ELCOOR(JJ,J)=COORD(JJ,IELCON(J,I))
          END DO
          DO K2=1,K1
            LML(K3+K2)=LM(I,K3+K2)
          END DO
          K3=K3+K1
        END DO

C
C       CALLS EMBEDED USER ELEMENT SUBROUTINE UEL ACCORDING TO IELT(I)
C               IELT()=1 8 NODED QUAD.
C
        SELECT CASE (IELT(I))
C
          CASE (1)
          CALL UEL_FEM8(SK,SM,ND,PROPS,NPROPS,ELCOOR,MCRD,NNE(I),I)
          CASE (2)
          CALL UEL_FEM9(SK,SM,ND,PROPS,NPROPS,ELCOOR,MCRD,NNE(I),I)
          CASE (3)
          CALL UEL_COS8(SK,SM,ND,PROPS,NPROPS,ELCOOR,MCRD,NNE(I),I)
C
        END SELECT
C
C       ASSEMBLE GLOBAL STIFFNESS AND MASS MATRICES
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
        DEALLOCATE(PROPS,SK,SM)
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
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE ASSEMLIS(NUMNP,NUMEL,MXNE,MXDOFEL,MXDOFDIM,NNE,NDOFN,
     1                    IELCON,LM,ID,IIN,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION NNE(NUMEL),IELCON(MXNE,NUMEL),LM(NUMEL,MXDOFEL),
     1          ID(MXDOFDIM,NUMNP),NDOFN(NUMNP)
C
      DO I=1,NUMEL
        K3=0
        DO J=1,NNE(I)
          K2I=NDOFN(IELCON(J,I))
          DO K2=1,K2I
            LM(I,K3+K2)=ID(K2,IELCON(J,I))
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
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM9(IDEL,NNE,NDOFEL,NTENS,XX,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),XJ(2,2),
     2          XJI(2,2)
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
      SUBROUTINE STDM8(IDEL,NNE,NDOFEL,NTENS,XX,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),XJ(2,2),
     2          XJI(2,2)
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
      SUBROUTINE STDM6(IDEL,NNE,NDOFEL,NTENS,XX,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),XJ(2,2),
     2          XJI(2,2)
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM8C:GENERATES THE STRAIN-DISPLACEMENT MATRIX B   C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 8-NODED 2D ELEMENT-PLANE STRAIN FULL COSSERAT ELEMENT     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM8C(IDEL,NNE,NDOFEL,NTENS,XX,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(2,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),XJ(2,2),
     1          XJI(2,2),SN(NNE)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions and its derivatives w.r.t natural coordinates
C
      CALL SHAPEF8(SN,R,S)
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
C     Full Cosserat element.
C
      DO I=1,NNE
        II=3*(I-1)+1
C
C       Normal components
C
        B(1,II  )= AUX1(1,I)
        B(2,II+1)= AUX1(2,I)
C
C       Shear components
C
        B(4,II+1)= AUX1(1,I)
        B(4,II+2)= SN(I)
        B(5,II  )= AUX1(2,I)
        B(5,II+2)= -SN(I)
C
C       Curvatures
C
        B(6,II+2)= AUX1(1,I)
        B(7,II+2)= AUX1(2,I)
C
      END DO
C
      RETURN
C
      END SUBROUTINE STDM8C
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
      END
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
      END
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
      END
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
      END
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
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                     MASS- MATRICES                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE CMASS8(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 8-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE RMASS8(RMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(8),
     1          RMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,8)
      CALL CLEAR(RMASS,NDOFEL,NDOFEL)
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
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,RMASS)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE RMASS8C(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 8-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE RMASS8C(RMASS,RO,RI,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=3)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(8),
     1          RMASS(NDOFEL,NDOFEL),AUX1(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,8)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(RMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF8(SN,RII,SII)
C
      SQRRO=DSQRT(RO)
      SQRRI=DSQRT(RI)
C
      DO I=1,8
        KK=3*I-2
        SF(1,  KK)=SQRRO*SN(I)
        SF(2,KK+1)=SQRRO*SN(I)
        SF(3,KK+2)=SQRRI*SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,RMASS)
C
      RETURN
C
      END SUBROUTINE RMASS8C
C

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE EL1D                                                  C
C                                                                      C
C     Returns the mass and stiffness matrices for a 1D element         C
C     with interpolation polynomials of degree N and constant          C
C     parameters RHO and E.                                            C
C                                                                      C
C                                                                      C
C     PARAMETERS       N:   Number of nodes                            C
C                      RHO: Mass density                               C 
C                      E:   Stiffness density                          C 
C                      SML: Lumped mass matrix                         C
C                      SKL: Stiffness matrix                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 28 MAY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012      
      SUBROUTINE EL1D(N,RHO,E,SML,SKL)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION X(N), W(N), XNODES(N)
        DIMENSION H(N), DIFH(N)
        DIMENSION SML(N,N),SKL(N,N)
        
          CALL LAG_NODES(N,XNODES,X,W)
          SML = 0.0
          RHO = 1.0
          DO K=1,N
            XK = X(K) ! Integration points
            CALL LAG_POL(N,XNODES,XK,H,DIFH)
            DO I=1,N
                DO J=1,N
                    SML(I,J) = SML(I,J) + H(I)*H(J)*RHO*W(K)
                    SKL(I,J) = SKL(I,J) + DIFH(I)*DIFH(J)*E*W(K)
                END DO
            END DO
          END DO
 
      END SUBROUTINE EL1D

CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE UEL_LAG                                               C
C                                                                      C
C     Returns the mass and stiffness matrices for a 2D element         C
C     with interpolation polynomials of degree NX in the X             C
C     direction and degree NY in the Y direction. For constant         C
C     parameters RHO and E. It uses complete Lagrange polynomials      C
C     for equispaced nodes.                                            C
C                                                                      C
C                                                                      C
C     PARAMETERS       N:   Number of nodes                            C
C                      RHO: Mass density                               C 
C                      E:   Stiffness density                          C 
C                      SML: Mass matrix                                C
C                      SKL: Stiffness matrix                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 3 JULY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
        SUBROUTINE UEL_LAG(SKL,SML,NDOFEL,PROPS,NPROPS,NX,NY,COORDS,
     1                         MCRD,JELEM)
            IMPLICIT REAL*8(A-H,O-Z)
            DIMENSION PROPS(NPROPS)
            DIMENSION C(3,3)
            DIMENSION COORDS(MCRD,NX*NY)
            DIMENSION X_LAG(NX), Y_LAG(NY)
            DIMENSION WX_LAG(NX), WY_LAG(NY)
            DIMENSION XNODES(NX), YNODES(NY)
            DIMENSION AJINV(2,2)
            DIMENSION H2D(NX*NY), DH2D(NX*NY,2)
            DIMENSION R(NDOFEL,NDOFEL)
            DIMENSION SML(NDOFEL,NDOFEL),SKL(NDOFEL,NDOFEL)
            
            
            E = PROPS(1)
            CPOIS = PROPS(2)
            RHO = PROPS(3)
            CALL LAG_NODES(NX,XNODES,X_LAG,WX_LAG)
            CALL LAG_NODES(NY,YNODES,Y_LAG,WY_LAG)
            
C           Plane Strain
            C(1,1) = (1.0-CPOIS)*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(1,2) = CPOIS*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(1,3) = 0.0
            C(2,1) = CPOIS*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(2,2) = (1.0-CPOIS)*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(2,3) = 0.0
            C(3,1) = 0.0
            C(3,2) = 0.0
            C(3,3) = E/(CPOIS+1.0)/2.0

C           Plane Stress
C            t = 1.0
C            C(1,1) = t*E/(1-CPOIS**2)
C            C(1,2) = CPOIS*t*E/(1-CPOIS**2)
C            C(1,3) = 0
C            C(2,1) = CPOIS*t*E/(1-CPOIS**2)
C            C(2,2) = t*E/(1-CPOIS**2)
C            C(2,3) = 0
C            C(3,1) = 0
C            C(3,2) = 0
C            C(3,3) = (1-CPOIS)*t*E/(1-CPOIS**2)/2.0E

            SML = 0.0
            SKL = 0.0
            DO IX=1,NX
                X = X_LAG(IX)
                WX = WX_LAG(IX)
                DO IY=1,NY
                    Y = Y_LAG(IY)
                    WY = WY_LAG(IY)
                    CALL LAG_POL2D(NX,NY,XNODES,YNODES,X,Y,H2D,DH2D)
                    CALL JACO(NX,NY,COORDS,DH2D,AJINV,DETJ)
                    CALL INTEGRAND(NX,NY,DH2D,AJINV,C,R)
                    DO II=1,NX*NY  !  Mass matrix
                        SML(2*II-1,2*II-1) = SML(2*II-1,2*II-1) + 
     1                                         DETJ*RHO*H2D(II)*WX*WY
                        SML(2*II,2*II) = SML(2*II,2*II) + 
     1                                         DETJ*RHO*H2D(II)*WX*WY
                    END DO
                    DO II=1,2*NX*NY  !  Stiffness matrix
                        DO JJ=1,2*NX*NY
                            SKL(II,JJ) = SKL(II,JJ) + 
     1                                      DETJ*R(II,JJ)*WX*WY 
                        END DO
                   END DO
                END DO
            END DO
 
        END SUBROUTINE UEL_LAG

CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE UEL_LAG_L                                             C
C                                                                      C
C     Returns the mass and stiffness matrices for a 2D element         C
C     with interpolation polynomials of degree NX in the X             C
C     direction and degree NY in the Y direction. For constant         C
C     parameters RHO and E. It uses complete Lagrange polynomials      C
C     for equispaced nodes.                                            C
C                                                                      C
C                                                                      C
C     PARAMETERS       N:   Number of nodes                            C
C                      RHO: Mass density                               C 
C                      E:   Stiffness density                          C 
C                      SML: Mass matrix                                C
C                      SKL: Stiffness matrix                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 3 JULY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
        SUBROUTINE UEL_LAG_L(SKL,SML,NDOFEL,PROPS,NPROPS,NX,NY,COORDS,
     1                         MCRD,JELEM)
            IMPLICIT REAL*8(A-H,O-Z)
            DIMENSION PROPS(NPROPS)
            DIMENSION C(3,3)
            DIMENSION COORDS(MCRD,NX*NY)
            DIMENSION X_LAG(NX), Y_LAG(NY)
            DIMENSION WX_LAG(NX), WY_LAG(NY)
            DIMENSION XNODES(NX), YNODES(NY)
            DIMENSION AJINV(2,2)
            DIMENSION H2D(NX*NY), DH2D(NX*NY,2)
            DIMENSION R(NDOFEL,NDOFEL)
            DIMENSION SML(NDOFEL,NDOFEL),SKL(NDOFEL,NDOFEL)
            
            
            E = PROPS(1)
            CPOIS = PROPS(2)
            RHO = PROPS(3)
            CALL LAG_NODES(NX,XNODES,X_LAG,WX_LAG)
            CALL LAG_NODES(NY,YNODES,Y_LAG,WY_LAG)
            
C           Plane Strain
            C(1,1) = (1.0-CPOIS)*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(1,2) = CPOIS*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(1,3) = 0.0
            C(2,1) = CPOIS*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(2,2) = (1.0-CPOIS)*E/((1.0-2.0*CPOIS)*(CPOIS+1.0))
            C(2,3) = 0.0
            C(3,1) = 0.0
            C(3,2) = 0.0
            C(3,3) = E/(CPOIS+1.0)/2.0

C           Plane Stress
C            t = 1.0
C            C(1,1) = t*E/(1-CPOIS**2)
C            C(1,2) = CPOIS*t*E/(1-CPOIS**2)
C            C(1,3) = 0
C            C(2,1) = CPOIS*t*E/(1-CPOIS**2)
C            C(2,2) = t*E/(1-CPOIS**2)
C            C(2,3) = 0
C            C(3,1) = 0
C            C(3,2) = 0
C            C(3,3) = (1-CPOIS)*t*E/(1-CPOIS**2)/2.0E

            SML = 0.0
            SKL = 0.0
            DO IX=1,NX
                X = X_LAG(IX)
                WX = WX_LAG(IX)
                DO IY=1,NY
                    Y = Y_LAG(IY)
                    WY = WY_LAG(IY)
                    CALL LAG_POL2D(NX,NY,XNODES,YNODES,X,Y,H2D,DH2D)
                    CALL JACO(NX,NY,COORDS,DH2D,AJINV,DETJ)
                    CALL INTEGRAND(NX,NY,DH2D,AJINV,C,R)
                    DO II=1,NX*NY  !  Mass matrix
                        SML(2*II-1,2*II-1) = SML(2*II-1,2*II-1) + 
     1                                         DETJ*RHO*H2D(II)*WX*WY
                        SML(2*II,2*II) = SML(2*II,2*II) + 
     1                                         DETJ*RHO*H2D(II)*WX*WY
                    END DO
                    DO II=1,2*NX*NY  !  Stiffness matrix
                        DO JJ=1,2*NX*NY
                            SKL(II,JJ) = SKL(II,JJ) + 
     1                                      DETJ*R(II,JJ)*WX*WY 
                        END DO
                   END DO
                END DO
            END DO


            DO I=1,16
                DO J =1,16
                    TMASS = TMASS + SML(I,J)
                END DO
                TRMASS = TRMASS + SML(I,I)
            END DO

            DO I=1,16
               AUX = SML(I,I)
               DO J=1,16
                   SML(I,J) = 0.0
                   SML(J,I) = 0.0
               END DO
               SML(I,I) = (AUX/TRMASS)*TMASS
            END DO
 
        END SUBROUTINE UEL_LAG_L

CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LAG_POL                                               C
C                                                                      C
C     Returns the value for the interpolant polynomials                C
C     evaluated at a point X, provided the interpolation               C
C     points for the interpolant Lagrange polynomial of order N.       C
C                                                                      C
C                                                                      C
C     PARAMETERS       N:   Number of nodes                            C
C                      XNODES: Interpolation points                     C 
C                      X:     Evaluation point                         C
C                      H:     Evaluated interpolator array             C
C                      DIFH:  Evaluated interpolator derivatives       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 3 JULY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012      
        SUBROUTINE LAG_POL(N,XNODES,X,H,DIFH)
            IMPLICIT REAL*8(A-H,O-Z)
            REAL*8 PROD1, PROD2
            DIMENSION XNODES(N),H(N),DIFH(N)

            DO I = 1,N
                PROD1 = 1.0
                PROD2 = 1.0
                XI = XNODES(I)
                DO J = 1,N  ! Interpolators
                    IF(I.NE.J) THEN
                        PROD1 = PROD1*(X - XNODES(J))
                        PROD2 = PROD2*(XI - XNODES(J))
                    END IF
                END DO        
                H(I) = PROD1/PROD2
         
                DIFH(I)=0.0
                DO J = 1,N ! Interpolator derivatives
                    IF(J.NE.I) THEN
                        PROD1 = 1.0
                        DO K = 1,N
                            IF(K.NE.I .and. K.NE.J) THEN
                                PROD1 = PROD1*(X - XNODES(K))
                            END IF
                        END DO
                        DIFH(I) = DIFH(I) + PROD1
                    END IF
                END DO
                DIFH(I) = DIFH(I)/PROD2
            END DO
        END SUBROUTINE LAG_POL
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LAG_POL2D                                             C C                                                                      C
C     Returns the value for the 2D interpolant polynomials and         C
C     its derivaties evaluated at a point (X,Y), provided              C
C     the points for the interpolant Lagrange polynomial of            C
C     order N-1.                                                       C
C                                                                      C
C                                                                      C
C     PARAMETERS       NX:   Number of nodes in X direction            C
C                      NY:   Number of nodes in Y direction            C 
C                      XNODES: Interpolation points in X direction     C 
C                      YNODES: Interpolation points in Y direction     C 
C                      X:     Evaluation point 1st component           C
C                      Y:     Evaluation point 2nd component           C
C                      H2D:     Evaluated interpolator array           C
C                      DH2D:  Evaluated interpolator derivatives       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 3 JULY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012      
        SUBROUTINE LAG_POL2D(NX,NY,XNODES,YNODES,X,Y,H2D,DH2D)
            IMPLICIT REAL*8(A-H,O-Z)
            REAL*8 PROD1, PROD2
            DIMENSION XNODES(NX), YNODES(NY)
            DIMENSION HX(NX), HY(NY), DIFHX(NX), DIFHY(NY)
            DIMENSION H2D(NX*NY),DH2D(NX*NY,2)
            
            CALL LAG_POL(NX,XNODES,X,HX,DIFHX)
            CALL LAG_POL(NY,YNODES,Y,HY,DIFHY)
                 
            DO IX=1,NX
                DO IY=1,NY
                    I = IX + (IY-1)*NX
                    H2D(I) = HX(IX)*HY(IY)
                    DH2D(I,1) = DIFHX(IX)*HY(IY)
                    DH2D(I,2) = HX(IX)*DIFHY(IY)
                END DO
            END DO
     
        END SUBROUTINE LAG_POL2D
CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LAG_NODES                                             C
C                                                                      C
C     Returns the Gauss-Legendre nodes(X) and weights(W)               C
C     for an interpolant Lagrange polynomial of order N-1.             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 28 MAY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012      
        SUBROUTINE LAG_NODES(N,XNODES,X,W)
            IMPLICIT REAL*8(A-H,O-Z)
            DIMENSION X(N), W(N), XNODES(N)

            XNODES(1) = -1.0; XNODES(N) = 1.0
            DO I=2,N-1
                XNODES(I) = -1.0 + (I-1)*2.0/REAL(N-1.0)
            END DO
            
            SELECT CASE(N)
              
              CASE(2)
                X(1) = -5.7735026918962573E-1
                X(2) = 5.7735026918962573E-1
                
                W(1) = 1.0E+0
                W(2) = 1.0E+0
                
                
              CASE(3)
                X(1) = -7.745966692414834E-1
                X(2) = 0.0E+0
                X(3) = 7.7459666924148329E-1
                
                W(1) = 5.5555555555555536E-1
                W(2) = 8.8888888888888884E-1
                W(3) = 5.5555555555555569E-1

              CASE(4)
                X(1) = -8.6113631159405257E-1
                X(2) = -3.3998104358485626E-1
                X(3) = 3.3998104358485626E-1
                X(4) = 8.6113631159405257E-1
                
                W(1) = 3.4785484513745396E-1
                W(2) = 6.5214515486254621E-1
                W(3) = 6.5214515486254621E-1
                W(4) = 3.4785484513745396E-1

              CASE(5)
                X(1) = -9.0617984593866585E-1
                X(2) = -5.3846931010568078E-1
                X(3) = 0.0E+0
                X(4) = 5.3846931010568311E-1
                X(5) = 9.0617984593866363E-1
                        
                W(1) = 2.3692688505618448E-1
                W(2) = 4.7862867049936808E-1
                W(3) = 5.6888888888888889E-1
                W(4) = 4.7862867049936625E-1
                W(5) = 2.3692688505618997E-1

              CASE(6)
                X(1) = -9.3246951420315194E-1
                X(2) = -6.6120938646626459E-1
                X(3) = -2.386191860831969E-1
                X(4) = 2.386191860831969E-1
                X(5) = 6.6120938646626493E-1
                X(6) = 9.3246951420315183E-1
                
                W(1) = 1.7132449237917094E-1
                W(2) = 3.60761573048138E-1
                W(3) = 4.6791393457269093E-1
                W(4) = 4.6791393457269093E-1
                W(5) = 3.6076157304813866E-1
                W(6) = 1.7132449237917066E-1
              
              CASE(7)
                X(1) = -9.4910791234275849E-1
                X(2) = -7.4153118559939479E-1
                X(3) = -4.0584515137739685E-1
                X(4) = 0.0E+0
                X(5) = 4.0584515137739718E-1
                X(6) = 7.4153118559939413E-1
                X(7) = 9.4910791234275882E-1
                        
                W(1) = 1.2948496616886973E-1
                W(2) = 2.7970539148927687E-1
                W(3) = 3.8183005050511948E-1
                W(4) = 4.179591836734694E-1
                W(5) = 3.8183005050511881E-1
                W(6) = 2.7970539148927676E-1
                W(7) = 1.2948496616886887E-1
              
              CASE(8)
                X(1) = -9.6028985649753584E-1
                X(2) = -7.9666647741362739E-1
                X(3) = -5.2553240991632888E-1
                X(4) = -1.8343464249564978E-1
                X(5) = 1.8343464249564978E-1
                X(6) = 5.2553240991632899E-1
                X(7) = 7.9666647741362651E-1
                X(8) = 9.6028985649753629E-1
                        
                W(1) = 1.0122853629037777E-1
                W(2) = 2.2238103445337584E-1
                W(3) = 3.1370664587788782E-1
                W(4) = 3.6268378337836182E-1
                W(5) = 3.6268378337836182E-1
                W(6) = 3.1370664587788749E-1
                W(7) = 2.223810344533737E-1
                W(8) = 1.0122853629037602E-1
              
              CASE(9)
                X(1) = -9.6816023950763419E-1
                X(2) = -8.3603110732662922E-1
                X(3) = -6.1337143270058692E-1
                X(4) = -3.2425342340380875E-1
                X(5) = 0.0E+0
                X(6) = 3.2425342340380892E-1
                X(7) = 6.1337143270058669E-1
                X(8) = 8.3603110732663133E-1
                X(9) = 9.6816023950763197E-1
                        
                W(1) = 8.1274388361553818E-2
                W(2) = 1.8064816069486569E-1
                W(3) = 2.6061069640293655E-1
                W(4) = 3.1234707704000314E-1
                W(5) = 3.3023935500125978E-1
                W(6) = 3.1234707704000292E-1
                W(7) = 2.6061069640293805E-1
                W(8) = 1.8064816069486445E-1
                W(9) = 8.1274388361559383E-2
              
              
              CASE(10)
                X(1) = -9.7390652851716797E-1
                X(2) = -8.6506336668899053E-1
                X(3) = -6.794095682990221E-1
                X(4) = -4.3339539412924727E-1
                X(5) = -1.4887433898163122E-1
                X(6) = 1.4887433898163122E-1
                X(7) = 4.3339539412924738E-1
                X(8) = 6.7940956829902277E-1
                X(9) = 8.6506336668898831E-1
                X(10) = 9.739065285171693E-1
                        
                W(1) = 6.6671344308694785E-2
                W(2) = 1.4945134915057012E-1
                W(3) = 2.1908636251598035E-1
                W(4) = 2.6926671930999635E-1
                W(5) = 2.9552422471475293E-1
                W(6) = 2.9552422471475293E-1
                W(7) = 2.6926671930999457E-1
                W(8) = 2.1908636251598526E-1
                W(9) = 1.4945134915057251E-1
                W(10) = 6.6671344308693217E-2
              
              CASE(11)
                X(1) = -9.7822865814604265E-1
                X(2) = -8.870625997681203E-1
                X(3) = -7.3015200557403748E-1
                X(4) = -5.1909612920681303E-1
                X(5) = -2.6954315595234468E-1
                X(6) = 0.0E+0
                X(7) = 2.6954315595234501E-1
                X(8) = 5.1909612920681225E-1
                X(9) = 7.3015200557404658E-1
                X(10) = 8.8706259976809987E-1
                X(11) = 9.7822865814605453E-1
                        
                W(1) = 5.5668567116209156E-2
                W(2) = 1.2558036946487536E-1
                W(3) = 1.8629021092773104E-1
                W(4) = 2.3319376459199076E-1
                W(5) = 2.628045445102466E-1
                W(6) = 2.7292508677790062E-1
                W(7) = 2.6280454451024687E-1
                W(8) = 2.3319376459198801E-1
                W(9) = 1.8629021092773856E-1
                W(10) = 1.255803694649068E-1
                W(11) = 5.5668567116179533E-2
                
              CASE DEFAULT
                WRITE(*,*) "N should be in [2,11]."
                STOP
                
            END SELECT        
            
        END SUBROUTINE LAG_NODES
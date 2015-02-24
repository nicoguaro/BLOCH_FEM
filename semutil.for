CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SPEL1D                                                C
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
      SUBROUTINE SPEL1D(N,RHO,E,SML,SKL)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION X_LOB(N), W_LOB(N)
        DIMENSION H(N), DIFH(N)
        DIMENSION SML(N,N),SKL(N,N)
        
          CALL LOB_NODES(N,X_LOB,W_LOB)
          SML = 0.0
          RHO = 1.0
          DO K=1,N
            XK = X_LOB(K) ! Integration points - GLL
            CALL LOB_POL(N,X_LOB,XK,H,DIFH)
            DO I=1,N
                DO J=1,N
                    SML(I,J) = SML(I,J) + H(I)*H(J)*RHO*W_LOB(K)
                    SKL(I,J) = SKL(I,J) + DIFH(I)*DIFH(J)*E*W_LOB(K)
                END DO
            END DO
          END DO
 
      END SUBROUTINE SPEL1D
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE UEL_SEM_LOB                                           C
C                                                                      C
C     Returns the mass and stiffness matrices for a 2D element         C
C     with interpolation polynomials of degree NX in the X             C
C     direction and degree NY in the Y direction. For constant         C
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
C     LAST MOD: 3 JUNE 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
        SUBROUTINE UEL_SEM_LOB(SKL,SML,NDOFEL,PROPS,NPROPS,NX,NY,COORDS,
     1                         MCRD,JELEM)
            IMPLICIT REAL*8(A-H,O-Z)
            DIMENSION PROPS(NPROPS)
            DIMENSION C(3,3)
            DIMENSION COORDS(MCRD,NX*NY)
            DIMENSION X_LOB(NX), Y_LOB(NY)
            DIMENSION WX_LOB(NX), WY_LOB(NY)
            DIMENSION AJINV(2,2)
            DIMENSION H2D(NX*NY), DH2D(NX*NY,2)
            DIMENSION R(NDOFEL,NDOFEL)
            DIMENSION SML(NDOFEL,NDOFEL),SKL(NDOFEL,NDOFEL)
            
            
            E = PROPS(1)
            CPOIS = PROPS(2)
            RHO = PROPS(3)
            CALL LOB_NODES(NX,X_LOB,WX_LOB)
            CALL LOB_NODES(NY,Y_LOB,WY_LOB)
            
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
                X = X_LOB(IX)
                WX = WX_LOB(IX)
                DO IY=1,NY
                    Y = Y_LOB(IY)
                    WY = WY_LOB(IY)
                    CALL LOB_POL2D(NX,NY,X_LOB,Y_LOB,X,Y,H2D,DH2D)
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
 
        END SUBROUTINE UEL_SEM_LOB
        
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE JACO                                                  C
C                                                                      C
C     Computes the Jacobian matrix for a NX*NY GLL element.            C
C                                                                      C
C                                                                      C
C     PARAMETERS       NX:   Number of nodes                           C
C                      NY: Number of nodes in Y                        C 
C                      COORDS: Array with coordinates for the          C
C                              for the undeformed element.             C
C                      E:   Stiffness density                          C 
C                      SML: Lumped mass matrix                         C
C                      SKL: Stiffness matrix                           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 2 June 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012      
      SUBROUTINE JACO(NX,NY,COORDS,DH2D,AJINV,DETJ)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION COORDS(2,NX*NY), DH2D(NX*NY,2)
        DIMENSION AJAC(2,2), AJINV(2,2)
 
            AJAC = 0.0
            DO I=1,NX*NY
                AJAC(1,1) = AJAC(1,1) + COORDS(1,I)*DH2D(I,1)
                AJAC(1,2) = AJAC(1,2) + COORDS(2,I)*DH2D(I,1)
                AJAC(2,1) = AJAC(2,1) + COORDS(1,I)*DH2D(I,2)
                AJAC(2,2) = AJAC(2,2) + COORDS(2,I)*DH2D(I,2)
            END DO
            
            DETJ = AJAC(1,1)*AJAC(2,2) - AJAC(1,2)*AJAC(2,1)
            
            AJINV(1,1) =  AJAC(2,2)/DETJ
            AJINV(1,2) = -AJAC(1,2)/DETJ
            AJINV(2,1) = -AJAC(2,1)/DETJ
            AJINV(2,2) =  AJAC(1,1)/DETJ
      END SUBROUTINE JACO
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE INTEGRAND                                             C
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
      SUBROUTINE INTEGRAND(NX,NY,DH2D,AJINV,C,R)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION DH2D(NX*NY,2),DHDX(NX*NY,2)
        DIMENSION C(3,3), AJINV(2,2), B(3,2*NX*NY)
        DIMENSION R(2*NX*NY,2*NX*NY)
        
        DO I=1,2  ! Derivatives interpolator in physical coordinates
            DO J=1,NX*NY
                DHDX(J,I) = DH2D(J,2)*AJINV(I,2)
     1                           + DH2D(J,1)*AJINV(I,1)
            END DO
        END DO
        
        B = 0.0  ! Displacement to Strains Matrix as Bathe (1995)
        DO I=1,NX*NY 
            B(1,2*I-1) = DHDX(I,1)
            B(2,2*I)   = DHDX(I,2)
            B(3,2*I-1) = DHDX(I,2)
            B(3,2*I)   = DHDX(I,1)
        END DO
        
        
        DO I=1,2*NX*NY  ! Integrand
            DO J=1,2*NX*NY
                R(I,J) = C(3,3)*B(3,I)*B(3,J) + C(2,3)*B(2,I)*B(3,J)
     1              + C(1,3)*B(1,I)*B(3,J) + B(2,J)*C(3,2)*B(3,I)
     2              + B(1,J)*C(3,1)*B(3,I) + C(2,2)*B(2,I)*B(2,J)
     3              + C(1,2)*B(1,I)*B(2,J) + B(1,J)*C(2,1)*B(2,I)
     4              + C(1,1)*B(1,I)*B(1,J)
            END DO
        END DO
 
      END SUBROUTINE INTEGRAND
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LOB_POL                                               C
C                                                                      C
C     Returns the value for the interpolant polynomials                C
C     evaluated at a point X, provided the GLL for the                 C
C     interpolant Lagrange polynomial of order N.                      C
C                                                                      C
C                                                                      C
C     PARAMETERS       N:   Number of nodes                            C
C                      X_LOB: GLL points                               C 
C                      X:     Evaluation point                         C
C                      H:     Evaluated interpolator array             C
C                      DIFH:  Evaluated interpolator derivatives       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 28 MAY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012      
        SUBROUTINE LOB_POL(N,X_LOB,X,H,DIFH)
            IMPLICIT REAL*8(A-H,O-Z)
            REAL*8 PROD1, PROD2
            DIMENSION X_LOB(N),H(N),DIFH(N)

            DO I = 1,N
                PROD1 = 1.0
                PROD2 = 1.0
                XI = X_LOB(I)
                DO J = 1,N  ! Interpolators
                    IF(I.NE.J) THEN
                        PROD1 = PROD1*(X - X_LOB(J))
                        PROD2 = PROD2*(XI - X_LOB(J))
                    END IF
                END DO        
                H(I) = PROD1/PROD2
         
                DIFH(I)=0.0
                DO J = 1,N ! Interpolator derivatives
                    IF(J.NE.I) THEN
                        PROD1 = 1.0
                        DO K = 1,N
                            IF(K.NE.I .and. K.NE.J) THEN
                                PROD1 = PROD1*(X - X_LOB(K))
                            END IF
                        END DO
                        DIFH(I) = DIFH(I) + PROD1
                    END IF
                END DO
                DIFH(I) = DIFH(I)/PROD2
            END DO
        END SUBROUTINE LOB_POL
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LOB_POL2D                                             C C                                                                      C
C     Returns the value for the 2D interpolant polynomials and         C
C     its derivaties evaluated at a point (X,Y), provided              C
C     the GLL for the interpolant Lagrange polynomial of               C
C     order N.                                                         C
C                                                                      C
C                                                                      C
C     PARAMETERS       NX:   Number of nodes in X direction            C
C                      NY:   Number of nodes in Y direction            C 
C                      X_LOB: GLL points in X direction                C 
C                      Y_LOB: GLL points in Y direction                C 
C                      X:     Evaluation point 1st component           C
C                      Y:     Evaluation point 2nd component           C
C                      H2D:     Evaluated interpolator array           C
C                      DH2D:  Evaluated interpolator derivatives       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     AUTHOR: NICOLAS GUARIN Z.                                        C
C     GRUPO DE MECANICA APLICADA - UNIVERSIDAD EAFIT                   C
C     LAST MOD: 29 MAY 2012                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012      
        SUBROUTINE LOB_POL2D(NX,NY,X_LOB,Y_LOB,X,Y,H2D,DH2D)
            IMPLICIT REAL*8(A-H,O-Z)
            REAL*8 PROD1, PROD2
            DIMENSION X_LOB(NX), Y_LOB(NY)
            DIMENSION HX(NX), HY(NY), DIFHX(NX), DIFHY(NY)
            DIMENSION H2D(NX*NY),DH2D(NX*NY,2)
            
            CALL LOB_POL(NX,X_LOB,X,HX,DIFHX)
            CALL LOB_POL(NY,Y_LOB,Y,HY,DIFHY)
                 
            DO IX=1,NX
                DO IY=1,NY
                    I = IX + (IY-1)*NX
                    H2D(I) = HX(IX)*HY(IY)
                    DH2D(I,1) = DIFHX(IX)*HY(IY)
                    DH2D(I,2) = HX(IX)*DIFHY(IY)
                END DO               
            END DO
     
        END SUBROUTINE LOB_POL2D
        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LOB_NODES                                             C
C                                                                      C
C     Returns the Gauss-Lobatto-Lagrange (GLL) nodes(X) and weights(W) C
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
        SUBROUTINE LOB_NODES(N,X,W)
            IMPLICIT REAL*8(A-H,O-Z)
            DIMENSION X(N), W(N)      
            
            SELECT CASE(N)
              
              CASE(2)
                X(1) = -1.0D+0
                X(2) = 1.0D+0
                
                W(1) = 1.0D+0
                W(2) = 1.0D+0
                
              CASE(3)
                X(1) = -1.0D+0
                X(2) = 0.0D+0
                X(3) = 1.0D+0
                
                W(1) = 3.3333333333333331D-1
                W(2) = 1.3333333333333333D+0
                W(3) = 3.3333333333333331D-1

              CASE(4)
                X(1) = -1.0D+0
                X(2) = -4.4721359549995793D-1
                X(3) = 4.4721359549995793D-1
                X(4) = 1.0D+0
                
                W(1) = 1.6666666666666666D-1
                W(2) = 8.3333333333333337D-1
                W(3) = 8.3333333333333326D-1
                W(4) = 1.6666666666666666D-1

              CASE(5)
                X(1) = -1.0D+0
                X(2) = -6.5465367070797709D-1 
                X(3) = 0.0D+0
                X(4) = 6.5465367070797709D-1
                X(5) = 1.0D+0
                
                W(1) = 1.0000000000000001D-1
                W(2) = 5.444444444444444D-1
                W(3) = 7.1111111111111114D-1
                W(4) = 5.4444444444444462D-1
                W(5) = 1.0000000000000001D-1

              CASE(6)
                X(1) = -1.0D+0
                X(2) = -7.6505532392946474D-1 
                X(3) = -2.852315164806451D-1
                X(4) = 2.852315164806451D-1
                X(5) = 7.6505532392946474D-1
                X(6) = 1.0D+0
                
                W(1) = 6.6666666666666666D-2
                W(2) = 3.7847495629784672D-1
                W(3) = 5.5485837703548657D-1
                W(4) = 5.5485837703548657D-1
                W(5) = 3.7847495629784694D-1
                W(6) = 6.6666666666666666D-2
              
              CASE(7)
                X(1) = -1.0D+0
                X(2) = -8.3022389627856696D-1 
                X(3) = -4.6884879347071429D-1
                X(4) = 0.0D+0
                X(5) = 4.6884879347071429D-1
                X(6) = 8.3022389627856696D-1
                X(7) = 1.0D+0
                
                W(1) = 4.7619047619047616D-2
                W(2) = 2.7682604736156574D-1
                W(3) = 4.3174538120986261D-1
                W(4) = 4.8761904761904762D-1
                W(5) = 4.3174538120986261D-1
                W(6) = 2.7682604736156607D-1
                W(7) = 4.7619047619047616D-2
              
              CASE(8)
                X(1) = -1.0D+0
                X(2) = -8.7174014850960491D-1
                X(3) = -5.9170018143314462D-1
                X(4) = -2.0929921790247885D-1
                X(5) = 2.0929921790247885D-1
                X(6) = 5.9170018143314451D-1
                X(7) = 8.7174014850960513D-1
                X(8) = 1.0D+0
                
                W(1) = 3.5714285714285712D-2
                W(2) = 2.1070422714350615D-1
                W(3) = 3.4112269248350408D-1
                W(4) = 4.1245879465870428D-1
                W(5) = 4.1245879465870428D-1
                W(6) = 3.4112269248350435D-1
                W(7) = 2.1070422714350615D-1
                W(8) = 3.5714285714285712D-2
              
              CASE(9)
                X(1) = -1.0D+0
                X(2) = -8.9975799541145984D-1
                X(3) = -6.7718627951073818D-1
                X(4) = -3.631174638261781D-1
                X(5) = 0.0D+0
                X(6) = 3.6311746382617821D-1
                X(7) = 6.7718627951073807D-1
                X(8) = 8.9975799541145995D-1
                X(9) = 1.0D+0
                
                W(1) = 2.7777777777777776D-2
                W(2) = 1.6549536156080569D-1
                W(3) = 2.7453871250016099D-1
                W(4) = 3.4642851097304617D-1
                W(5) = 3.7151927437641724D-1
                W(6) = 3.4642851097304606D-1
                W(7) = 2.7453871250016132D-1
                W(8) = 1.654953615608056D-1
                W(9) = 2.7777777777777776D-2
              
              
              CASE(10)
                X(1) = -1.0D+0
                X(2) = -9.195339081664603D-1
                X(3) = -7.3877386510550336D-1
                X(4) = -4.7792494981044442D-1
                X(5) = -1.6527895766638703D-1
                X(6) = 1.6527895766638703D-1
                X(7) = 4.7792494981044453D-1
                X(8) = 7.3877386510550502D-1
                X(9) = 9.1953390816645864D-1
                X(10) = 1.0D+0
                
                W(1) = 2.2222222222222223D-2
                W(2) = 1.3330599085106998D-1
                W(3) = 2.2488934206312614D-1
                W(4) = 2.9204268367968328D-1
                W(5) = 3.2753976118389755D-1
                W(6) = 3.2753976118389755D-1
                W(7) = 2.9204268367968406D-1
                W(8) = 2.2488934206312591D-1
                W(9) = 1.3330599085107009D-1
                W(10) = 2.2222222222222223D-2
              
              CASE(11)
                X(1) = -1.0D+0
                X(2) = -9.3400143040805805D-1
                X(3) = -7.8448347366314708D-1
                X(4) = -5.6523532699620305D-1
                X(5) = -2.9575813558693925D-1
                X(6) = 0.0D+0
                X(7) = 2.9575813558693936D-1
                X(8) = 5.6523532699620271D-1
                X(9) = 7.8448347366314808D-1
                X(10) = 9.3400143040805716D-1
                X(11) = 1.0D+0
                
                W(1) = 1.8181818181818181D-2
                W(2) = 1.096122732669947D-1
                W(3) = 1.8716988178030528D-1
                W(4) = 2.4804810426402857D-1
                W(5) = 2.8687912477900762D-1
                W(6) = 3.0021759545569077D-1
                W(7) = 2.8687912477900812D-1
                W(8) = 2.4804810426402779D-1
                W(9) = 1.8716988178030508D-1
                W(10) = 1.0961227326699495D-1
                W(11) = 1.8181818181818181D-2
                
              CASE DEFAULT
                WRITE(*,*) "N should be in [2,11]."
                STOP
                
            END SELECT        
            
        END SUBROUTINE LOB_NODES
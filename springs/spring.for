      PROGRAM SPRING
        IMPLICIT REAL*8(A-H,O-Z)
        REAL*8 SK(16,16), SM(16,16), COORDS(2,16)
        INTEGER IRNODES_WO(6), IRNODES(10), IMNODES(10)
        REAL*8 EIGVALS(41*41,8)
        PARAMETER (PI = 3.141592653589793)

        OPEN(8,FILE ="spring.fre")

        IRNODES_WO(1) = 1; IRNODES_WO(4) = 10
        IRNODES_WO(2) = 2; IRNODES_WO(5) = 15
        IRNODES_WO(3) = 9; IRNODES_WO(6) = 16


        IRNODES(1) = 1; IRNODES(6) = 2
        IRNODES(2) = 2; IRNODES(7) = 15
        IRNODES(3) = 9; IRNODES(8) = 16
        IRNODES(4) = 10; IRNODES(9) = 1
        IRNODES(5) = 1; IRNODES(10) = 2

        IMNODES(1) = 7; IMNODES(6) = 6
        IMNODES(2) = 8; IMNODES(7) = 11
        IMNODES(3) = 13; IMNODES(8) = 12
        IMNODES(4) = 14; IMNODES(9) = 3
        IMNODES(5) = 5; IMNODES(10) = 4

        CALL MAT(16,SK,SM)

        COORDS(1,1) = -1.0; COORDS(2,1) = -1.0
        COORDS(1,2) = -1.0; COORDS(2,2) = -1.0
        COORDS(1,3) = 1.0; COORDS(2,3) = -1.0
        COORDS(1,4) = 1.0; COORDS(2,4) = -1.0
        COORDS(1,5) = 1.0; COORDS(2,5) = 1.0
        COORDS(1,6) = 1.0; COORDS(2,6) = 1.0
        COORDS(1,7) = -1.0; COORDS(2,7) = 1.0
        COORDS(1,8) = -1.0; COORDS(2,8) = 1.0
        COORDS(1,9) = 0.0; COORDS(2,9) = -1.0
        COORDS(1,10) = 0.0; COORDS(2,10) = -1.0
        COORDS(1,11) = 1.0; COORDS(2,11) = 0.0
        COORDS(1,12) = 1.0; COORDS(2,12) = 0.0
        COORDS(1,13) = 0.0; COORDS(2,13) = 1.0
        COORDS(1,14) = 0.0; COORDS(2,14) = 1.0
        COORDS(1,15) = -1.0; COORDS(2,15) = 0.0
        COORDS(1,16) = -1.0; COORDS(2,16) = 0.0



        NCOR = 16
        NDOF = 16
        NCOND = 10
        NCOND_WO = 6

        AKXMIN = -PI/2.0
        AKYMIN = -PI/2.0
        AKXMAX = PI/2.0
        AKYMAX = PI/2.0

        NEVALS = 6

        NKX = 41
        NKY = 41


        CALL BLOCH(NDOF,SK,SM,NCOR,COORDS,NCOND_WO,NCOND,
     1           IRNODES_WO,IMNODES,IRNODES,NEVALS,AKXMIN,
     2           AKYMIN,AKXMAX,AKYMAX,NKX,NKY,EIGVALS)

        DO I=1,NKX*NKY
            DO J=1,NEVALS+2
                WRITE(8,*) EIGVALS(I,J)
            END DO
        END DO

      END PROGRAM SPRING

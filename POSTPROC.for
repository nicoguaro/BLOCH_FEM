CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SVARSUPD                                              C
C                                                                      C
C     IFLAG=0 UPDATES                                                  C
C     IFLAG=1 RESETS                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE SVARSUPD(UH,RH,UG,RF,SVARSEGPT,SVARSGPTH,MXNO,MXINC,
     1                    MXEQ,MXEL,MXSVARS,KINC,IFLAG)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION RH(MXEQ,MXINC),UH(MXEQ,MXINC),UG(MXEQ),RF(MXEQ),
     1          SVARSEGPT(MXEL,MXSVARS),SVARSGPTH(MXEL,MXSVARS,MXINC)
C
      IF(IFLAG.EQ.0) THEN
        CALL VECHAN2D3D(SVARSGPTH,MXEL,MXSVARS,MXINC,KINC,SVARSEGPT,
     1                  MXEL,MXSVARS,IFLAG)
        DO K1=1,MXEQ
          UH(K1,KINC)=UG(K1)
          RH(K1,KINC)=RF(K1)
        END DO
      END IF
C
      IF(IFLAG.EQ.1) THEN
        IF(KINC.EQ.1) THEN
          CALL CLEAR(SVARSEGPT,MXEL,MXSVARS)
        ELSE
          CALL VECHAN2D3D(SVARSGPTH,MXEL,MXSVARS,MXINC,KINC-1,SVARSEGPT,
     1                    MXEL,MXSVARS,IFLAG)
        END IF
      END IF
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE RESOUT                                                C
C                                                                      C
C     SVARSGPTH(:,:,:):SVARS HISTORIES                                 C
C     IDOEL(:)        :ELEMENTS FOR OUTPUT                             C
C     IDGP(:)         :GAUSS POINTS FOR OUTPUT                         C
C     IDSV(:)         :SVARS FOR OUTPUT                                C
C     NOUTEL          :NUMBER OF ELEMENTS FOR OUTPUT                   C
C     NGPRES          :NUMBER OF GAUSS POINTS FOR OUTPUT               C
C     NSVRES          :NUMBER OF SVARS FOR OUTPUT                      C
C     NSVARSE         :NUMBER OF DECLARED SVARS PER ELEMENT            C
C     NSVARSP         :NUMBER OF DECLARED SVARS PER GAUSS POINT        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE RESOUT(SVARSGPTH,IDOEL,IDGP,IDSV,NOUTEL,NGPRES,NSVRES,
     1                  NSVARSE,NSVARSP,NUMEL,NINC,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IDGP(NGPRES),IDSV(NSVRES),SVARSGPTH(NUMEL,NSVARSE,NINC),
     1          IDOEL(NOUTEL)
C
      DO LL=1,NOUTEL
         LI=IDOEL(LL)
         WRITE(IOUT,*) 'ELE ID   ',LI
         DO II=1,NSVRES
            I=IDSV(II)
            WRITE(IOUT,*) 'SDV   ',I
            DO JJ=1,NGPRES
              J=IDGP(JJ)
              WRITE(IOUT,*) 'G POINT-',J
              K=NSVARSP*J-NSVARSP+I
              DO KK=1,NINC
                 WRITE(IOUT,6030) SVARSGPTH(LI,K,KK)
              END DO
            END DO
         END DO
      END DO
C
 6030 FORMAT(F12.3)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE WRSVARS                                               C
C                                                                      C
C     IPROB=1 GENERAL PROBLEM                                          C
C     IPROB=2 STRESS-STRAIN PROBLEM                                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE WRSVARS(NUMNP,NUMEL,NDOFDIM,ID,IELT,SVARSGPTH,UH,RH,
     1                   NSVARS,NEQ,NINC,INC,IOUT,IPROB)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ID(NDOFDIM,NUMNP),D(3),R(3),UH(NEQ,NINC),RH(NEQ,NINC),
     1          SVARSGPTH(NUMEL,NSVARS,NINC),SIGGPT(NUMEL,4,9),
     2          STRGPT(NUMEL,4,9),STRGPE(NUMEL,4,9),STRGPP(NUMEL,4,9),
     3          EQPGPT(NUMEL,9),IELT(NUMEL)
C
C     G E N E R A L  P R O B L E M-IPROB=1
C
      IF(IPROB.EQ.1) THEN
         WRITE(IOUT,2090) INC
ccc         WRITE(IOUT,3000)
ccc         WRITE(IOUT,6000)
ccc         WRITE(IOUT,6010)
ccc         WRITE(IOUT,6020)
ccc         DO K1=1,NUMEL
ccc           DO K3=1,4,4
ccc             WRITE(IOUT,6030) K1,1,SVARSGPTH(K1,K3,INC),
ccc     1       SVARSGPTH(K1,K3+1,INC),SVARSGPTH(K1,K3+2,INC),
ccc     2       SVARSGPTH(K1,K3+3,INC)
ccc           END DO
ccc         END DO
      END IF
C
C     S T R E S S E S A N D  S T R A I N S-IPROB=2
C     (Convert svars to stress, strains, etc)
C
      IF(IPROB.EQ.2) THEN
         DO K1=1,NUMEL
           if(ielt(k1).eq.3.or.ielt(k1).eq.4) then
              ngpt=4
           else
              ngpt=9
           end if
           DO NN=1,ngpt
             ISVINT=1+(NN-1)*NSVARS/ngpt
             JJ=1
             DO II=ISVINT,ISVINT+3
               SIGGPT(K1,JJ,NN)=SVARSGPTH(K1,II,INC)
               STRGPT(K1,JJ,NN)=SVARSGPTH(K1,II+4,INC)
               STRGPE(K1,JJ,NN)=SVARSGPTH(K1,II+8,INC)
               STRGPP(K1,JJ,NN)=SVARSGPTH(K1,II+12,INC)
               JJ=JJ+1
             END DO
             EQPGPT(K1,NN)=SVARSGPTH(K1,ISVINT+16,INC)
           END DO
         END DO
C
         WRITE(IOUT,2090) INC
         WRITE(IOUT,3000)
C
         DO K1=1,NUMEL
           if(ielt(k1).eq.3.or.ielt(k1).eq.4) then
              ngpt=4
           else
              ngpt=9
           end if
C
C          Stresses at the Gauss points
C
           WRITE(IOUT,5000)
           WRITE(IOUT,5020)
           WRITE(IOUT,6020)
           DO K3=1,ngpt
             WRITE(IOUT,6030) K1,K3,SIGGPT(K1,1,K3),SIGGPT(K1,2,K3),
     1       SIGGPT(K1,3,K3),SIGGPT(K1,4,K3)
           END DO
C
C          Strains at the Gauss points
C
           WRITE(IOUT,7000)
           WRITE(IOUT,7010)
           WRITE(IOUT,6020)
           DO K3=1,ngpt
             WRITE(IOUT,7030)K1,K3,STRGPT(K1,1,K3),STRGPT(K1,2,K3),
     1       STRGPT(K1,3,K3),STRGPT(K1,4,K3)
           END DO
C
C          Elastic Strains at the Gauss points
C
           WRITE(IOUT,8000)
           WRITE(IOUT,7010)
           WRITE(IOUT,6020)
           DO K3=1,ngpt
             WRITE(IOUT,7030)K1,K3,STRGPE(K1,1,K3),STRGPE(K1,2,K3),
     1       STRGPE(K1,3,K3),STRGPE(K1,4,K3)
           END DO
C
C          Plastic Strains at the Gauss points
C
           WRITE(IOUT,9000)
           WRITE(IOUT,7010)
           WRITE(IOUT,6020)
           DO K3=1,ngpt
             WRITE(IOUT,7030)K1,K3,STRGPP(K1,1,K3),STRGPP(K1,2,K3),
     1       STRGPP(K1,3,K3),STRGPP(K1,4,K3)
           END DO
C
           WRITE(IOUT,*)
         END DO
      END IF
C
C     Displacements and reactions
C
      WRITE (IOUT,2000)
      DO II=1,NUMNP
        DO I=1,3
          D(I)=0.0D0
          R(I)=0.0D0
        END DO
        D(1)=UH(3*II-2,INC)
        D(2)=UH(3*II-1,INC)
        D(3)=UH(3*II,  INC)
        R(1)=RH(3*II-2,INC)
        R(2)=RH(3*II-1,INC)
        R(3)=RH(3*II  ,INC)
        WRITE(IOUT,2010) II,D,R
      END DO
C
      RETURN
C
C     SOLUTION OUTPUT FORMAT
C
 2000 FORMAT (///,15X,'N  O  D  A  L   S  O  L  U  T  I  O  N',
     1       //,'  NODE',18X,'U1',15X,'U2',15X,'U3',15X,'RF-1',15X,
     2       'RF-2',15X,'RF-3')
C
 2010 FORMAT (1X,I3,8X,6E18.6)
c
 2090 FORMAT (///,15X,'INCREMENT',I5,2X,'SUMMARY',/)
C
 3000 FORMAT (//,20X'E  L  E  M  E  N  T    S  O  L  U  T  I  O  N',
     1       /,30X,'(VALUES AT THE GAUSS POINTS)')
C
 6000 FORMAT(/,'  SVARS     ',/)
 6010 FORMAT('ELEMENT  GPT      SDV-1         SDV-2        SDV-3       
     1SDV-4')
 6020 FORMAT('  ID      ID')
C
 5000 FORMAT(/,' STRESS     ',/)
 7000 FORMAT(/,' TOTAL STRAINS',/)
 8000 FORMAT(/,' ELASTIC STRAINS',/)
 9000 FORMAT(/,' PLASTIC STRAINS',/)
C
 5020 FORMAT('ELEMENT GPT       SIG-X        SIG-Y        SIG-Z    
     1   TAO-XY')
 7010 FORMAT('ELEMENT GPT       EPS-X        EPS-Y        EPS-Z     
     1GAMMA-XY')
C
 6030 FORMAT(I3,5X,I3,X,4(X,F12.3))
 7030 FORMAT(I3,5X,I3,X,4(X,F12.8))
C
      END SUBROUTINE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE WRMSG                                                 C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE WRMSG(ICALL)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      IF(ICALL.EQ.1) WRITE(*,*) 'READING MODEL DATA'
      IF(ICALL.EQ.2) WRITE(*,*) 'GENERATING ASSEMBLY LIST'
      IF(ICALL.EQ.3) WRITE(*,*) 'STARTING INCREMNTAL SOLUTION'
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SVARSUPD                                              C
C                                                                      C
C     IFLAG=0 UPDATES                                                  C
C     IFLAG=1 RESETS                                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE STRSUPD(SVARSEGPT,SVARSGPTH,MXNO,MXINC,MXEQ,MXEL,
     1                   MXSVARS,KINC,IFLAG)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION SVARSEGPT(MXEL,MXSVARS),SVARSGPTH(MXEL,MXSVARS,MXINC)
C
      IF(IFLAG.EQ.0) THEN
        CALL VECHAN2D3D(SVARSGPTH,MXEL,MXSVARS,MXINC,KINC,SVARSEGPT,
     1                  MXEL,MXSVARS,IFLAG)
      END IF
C
      IF(IFLAG.EQ.1) THEN
        IF(KINC.EQ.1) THEN
          CALL CLEAR(SVARSEGPT,MXEL,MXSVARS)
        ELSE
          CALL VECHAN2D3D(SVARSGPTH,MXEL,MXSVARS,MXINC,KINC-1,SVARSEGPT,
     1                    MXEL,MXSVARS,IFLAG)
        END IF
      END IF
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE DISPUPD(UH,RH,UG,RF,MXINC,MXEQ,KINC)                  C
C                                                                      C
C     UH(:,:)         :NODAL DISPLACEMENTS HISTORY                     C
C     RH(:,:)         :NODAL FORCES HISTORY                            C
C     MXEQ            :TOTAL NUMBER OF EQUATIONS                       C
C     MXINC           :NUMBER OF INCREMENTS                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE DISPUPD(UH,RH,UG,RF,MXINC,MXEQ,KINC)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION RH(MXEQ,MXINC),UH(MXEQ,MXINC),UG(MXEQ),RF(MXEQ)
C
      DO K1=1,MXEQ
        UH(K1,KINC)=UG(K1)
        RH(K1,KINC)=RF(K1)
      END DO
C
      RETURN
C
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE WRSTRS(NUMNP,NUMEL,IELT,SVARSGPTH,NSVARS,NINC,IOUT)   C
C                                                                      C
C     NUMNP           :NUMBER OF NODAL POINTS                          C
C     NUMEL           :NUMBER OF ELEMENTS                              C
C     IELT(:)         :ELEMENT TYPE ID'S                               C
C     SVARSGPTH(:,:,:):SVARS HISTORIES                                 C
C     NSVARS          :MAXIMUM NUMBER OF SVARS IN A GIVEN ELEMENT      C
C     NINC            :NUMBER OF INCREMENTS                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE WRSTRS(NUMNP,NUMEL,IELT,SVARSGPTH,NSVARS,NINC,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION SVARSGPTH(NUMEL,NSVARS,NINC),IELT(NUMEL)
C
      DO I=1,NINC
        WRITE(IOUT,500) I
        DO J=1,NUMEL
          WRITE(IOUT,510) J
C
          SELECT CASE(IELT(J))
C
            CASE(1)
            NTSVARS=72
            NGPTS=9
            NSAL=NTSVARS/NGPTS
            DO KK=1,NGPTS
              IGP=1+(KK-1)*NSAL
              WRITE(IOUT,1000) KK, (SVARSGPTH(J,IGP+JJ-1,I),JJ=1,8)
            END DO
C
            CASE(2)
            NTSVARS=72
            NGPTS=9
            NSAL=NTSVARS/NGPTS
            DO KK=1,NGPTS
              IGP=1+(KK-1)*NSAL
              WRITE(IOUT,1000) (SVARSGPTH(J,IGP+JJ,I),JJ=1,7)
            END DO
            CASE(3)
            NTSVARS=56
            NGPTS=7
            NSAL=NTSVARS/NGPTS
            DO KK=1,NGPTS
              IGP=1+(KK-1)*NSAL
              WRITE(IOUT,1000) (SVARSGPTH(J,IGP+JJ,I),JJ=1,7)
            END DO
C
          END SELECT
        END DO
      END DO

C
      RETURN
C
C     SOLUTION OUTPUT FORMAT
C
  500 FORMAT (//,15X,'INCREMENT',I5,/)
  510 FORMAT (/,15X,'ELEMENT',I5,/)
 1000 FORMAT (1X,'GP',I,8(X,F12.3))
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE WRDISP(NUMNP,NDOFDIM,ID,UH,RH,NEQ,NINC,IOUT)          C
C                                                                      C
C     NUMNP           :NUMBER OF NODAL POINTS                          C
C     NDOFDIM         :PROBLEM DIMENSION (2D, 3D)                      C
C     ID              :EQUATION NUMEBERS ASSIGNED TO EACH NODE         C
C     UH(:,:)         :NODAL DISPLACEMENTS HISTORY                     C
C     RH(:,:)         :NODAL FORCES HISTORY                            C
C     NEQ             :TOTAL NUMBER OF EQUATIONS                       C
C     NINC            :NUMBER OF INCREMENTS                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE WRDISP(NUMNP,NDOFDIM,ID,UH,RH,NEQ,NINC,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ID(NDOFDIM,NUMNP),D(2),R(2),UH(NEQ,NINC),RH(NEQ,NINC)
C
      WRITE (IOUT,500)
      DO INC=1,NINC
        WRITE(IOUT,1000) INC
        DO I=1,NUMNP
          DO J=1,NDOFDIM
            K=ID(J,I)
            IF (K.EQ.0) THEN
               D(J)=0.D0
               R(J)=1.D0
            ELSE
               D(J)=UH(K,INC)
               R(J)=RH(K,INC)
            END IF
          END DO
          WRITE(IOUT,2010) I,D
        END DO
      END DO
C
      RETURN
C
C     SOLUTION OUTPUT FORMAT
C
  500 FORMAT (///,15X,'N  O  D  A  L   S  O  L  U  T  I  O  N')
C
 1000 FORMAT (//,15X,'INCREMENT',I5,/)
 1010 FORMAT (/,15X,'NODE',I5,/)
 2010 FORMAT (1X,I3,8X,6E18.6)
C
      END SUBROUTINE
C
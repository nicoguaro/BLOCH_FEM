C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C SUBROUTINE NODINP                                                    C
C                                                                      C
C      NUMNP    :NUMBER OF NODAL POINTS                                C
C      ID       :EQUATION NUMEBERS ASSIGNED TO EACH NODE               C
C      MXDOFDIM :PROBLEM DIMENSION (2D, 3D)                            C
C      COORD    :NODAL COORDINATES                                     C
C      NDOFN    :NODAL DOF PER NODE                                    C
C      NEQ      :TOTAL NUMBER OF EQUATIONS                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE NODINP(NUMNP,ID,MXDOFDIM,COORD,NDOFN,NEQ,IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ID(MXDOFDIM,NUMNP),COORD(MXDOFDIM,NUMNP),NDOFN(NUMNP)
C
      CALL CLEAR(COORD,MXDOFDIM,NUMNP)
      CALL CLEARIM(ID,MXDOFDIM,NUMNP)
C
C     READ AND GENERATE NODAL POINT DATA
C
      WRITE(IOUT,1000)
      WRITE(IOUT,1005)
C
      DO I=1,NUMNP
        READ(IIN,     *) N,NDOFN(N),(ID(II,I),II=1,NDOFN(N)),
     1  (COORD(JJ,I),JJ=1,MXDOFDIM)
        WRITE(IOUT,1010)N,NDOFN(N),(ID(II,I),II=1, MXDOFDIM),
     1  (COORD(JJ,I),JJ=1,MXDOFDIM)
      END DO
C
C     ASSIGN EQUATION NUMBERS TO ACTIVE DOF
C
      ICOUNT=1
      DO K1=1,NUMNP
        NDOF=NDOFN(K1)
        DO K2=1,NDOFN(K1)
          IF(ID(K2,K1).EQ.0) THEN
            ID(K2,K1)=ICOUNT
            ICOUNT=ICOUNT+1
          ELSE
            ID(K2,K1)=0
          END IF
        END DO
      END DO
      NEQ=ICOUNT-1
C
C     WRITE EQUATION NUMBERS
C
      WRITE(IOUT,1020)
      WRITE(IOUT,1025)
      WRITE(IOUT,1030) (N,(ID(I,N),I=1,MXDOFDIM),N=1,NUMNP)
C
 1000 FORMAT(//,6X,' N O D A L  D A T A',//)
 1005 FORMAT('ID',2X,'NDOF',2X,'BC-X',2X,'BC-Y',2X,'BC-R',
     1        6X,'COORD-X',6X,'COORD-Y',/)
 1010 FORMAT(I3,X,I2,4X,I2,4X,I2,4X,I2,6X,F5.2,6X,F5.2)
 1020 FORMAT(/,'E Q U A T I O N  N U M B E R S',/)
 1025 FORMAT(/,'NODE',6X,'DEGREES OF FREEDOM',/,
     1      'NUMBER',8X,'U',9X,'V'/)
 1030 FORMAT(I5,5X,I5,5X,I5)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C SUBROUTINE MATINP                                                    C
C                                                                      C
C     NUMAT       :NUMBER OF MATERIAL PROFILES                         C
C     NMPR        :MAXIMUM NUMBER OF MATERIAL PROPERTIES               C
C     AMATE       :MATERIAL PROFILES ARRAY                             C
C     NMATP       :NUMBER OF MATERIAL PROPERTIES IN EACH PROFILE       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE MATINP(NUMAT,NMPR,AMATE,NMATP,IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AMATE(NMPR,NUMAT),NMATP(NUMAT)
C
      CALL CLEAR(AMATE,NMPR,NUMAT)
C
C     READS AND GENERATES MATERIAL DATA
C
      WRITE (IOUT,1000)
      WRITE (IOUT,3000)
      DO I=1,NUMAT
        READ  (IIN, *) N,NMATP(N),(AMATE(JJ,N),JJ=1,NMATP(N))
        WRITE (IOUT,3020) N,NMATP(N),(AMATE(JJ,N),JJ=1,NMATP(N))
      END DO
C
 1000 FORMAT(//,4X,'M A T E R I A L  D A T A',/)
 3000 FORMAT(' MAT.ID',10X,'YOUNG.S MODULUS',10X,'POIS. RATIO',10X,
     1'YIEL.1',10X,'EP.1',12X,'YIEL.2',8X,'EP.2',//)
 3020 FORMAT(I5,3X,I5,12X,6(5X,F13.3))
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE ELEINP                                                C
C                                                                      C
C     NUMNP           :NUMBER OF NODAL POINTS                          C
C     NUMEL           :NUMBER OF ELEMENTS                              C
C     NNE(:)          :NUMBER OF NODES ASSIGNED TO EACH ELEMENT        C
C     IELT(:)         :ELEMENT TYPE ID'S                               C
C     NDOFEL(:)       :NUMBER OF DOF AT EACH ELEMENT                   C
C     NMNE            :MAXIMUM NUMBER OF NODES AT A GIVEN ELEMENT      C
C     MATP(:)         :MATERIAL PROFILE ASSIGNED TO EACH ELEMENT       C
C     IELCON(:,:)     :ELEMENT CONNECTIVITIES                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE ELEINP(NUMNP,NUMEL,NNE,IELT,NDOFEL,NMNE,MATP,IELCON,
     1                  IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION NNE(NUMEL),IELT(NUMEL),MATP(NUMEL),IELCON(NMNE,NUMEL),
     1          NDOFEL(NUMEL)
C
C     Creates IELCON()
C
      WRITE(IOUT,1000)
      WRITE(IOUT,1001)
      WRITE(IOUT,1002)
      DO I=1,NUMEL
        READ(IIN,*) M,IELT(M),NDOFEL(M),MATP(M),NNE(M),
     1  (IELCON(J,M),J=1,NNE(M))
        WRITE(IOUT,1010) M,IELT(M),NDOFEL(M),MATP(M),NNE(M),
     1  (IELCON(J,M),J=1,NNE(M))
      END DO
C
 1000 FORMAT(///,8X,'E L E M E N T  I N F O R M A T I O N',//)
 1001 FORMAT('ID  TYPE NDOF MAT NNEL NODE NODE NODE NODE NODE NODE
     1 NODE NODE')
 1002 FORMAT(25X'1    2    3    4    5    6   7    8',/)
 1010 FORMAT(5(2X,I2),X,8(2X,I3))
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LOADV                                                 C
C                                                                      C
C     NUMNP           :NUMBER OF NODAL POINTS                          C
C     NDOFDIM         :PROBLEM DIMENSION (2D, 3D)                      C
C     ID(:)           :EQUATION NUMEBERS ASSIGNED TO EACH NODE         C
C     PP(:)           :NODAL LOADS VECTOR                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE LOADV(NUMNP,NDOFDIM,ID,PP,NEQ,IIN,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION ID(NDOFDIM,NUMNP),PP(NEQ)
C
      CALL CLEARV(PP,NEQ)
C
      READ (IIN,    *) NPLOAD
      WRITE(IOUT,2010) NPLOAD
C
C     Assembles Point Loads Vector
C
      IF(NPLOAD.NE.0) CALL PLOADS(NUMNP,NDOFDIM,NPLOAD,ID,PP,NEQ,IIN,
     1                            IOUT)
      WRITE(IOUT,2030)
C
      DO I=1,NEQ
        WRITE(IOUT,2040) I,PP(I)
      END DO
C
 1010 FORMAT (3I5)
 2010 FORMAT (///,'   NUMBER OF CONCENTRATED LOADS .= ',I5)
 2030 FORMAT (//,'----LOAD VECTORS------',//,
     1        'DOF-#        P-LOAD',/)
 2040 FORMAT (I5,4X,F10.5)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C SUBROUTINE PLOADS                                                    C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PLOADS(NUMNP,NDOFDIM,NPLOAD,ID,PR,NEQ,IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ID(NDOFDIM,NUMNP),PR(NEQ),NOD(NPLOAD),IDIRN(NPLOAD),
     1          FLOAD(NPLOAD)
C
      CALL CLEARIV (NOD,NPLOAD)
      CALL CLEARV (FLOAD, NPLOAD)
      CALL CLEARIV(IDIRN,NPLOAD)
      CALL CLEARV (PR,NEQ)
C
      WRITE (IOUT,2000)
      DO K1=1,NPLOAD
        READ(IIN,*) NOD(K1),IDIRN(K1),FLOAD(K1)
      END DO
C
      DO K1=1,NPLOAD
         WRITE (IOUT,2010) NOD(K1),IDIRN(K1),FLOAD(K1)
      END DO
C
      DO L=1,NPLOAD
        LN=NOD(L)
        LI=IDIRN(L)
        II=ID(LI,LN)
        IF (II.NE.0) PR(II)=PR(II)+FLOAD(L)
      END DO
C
 2000 FORMAT(//,'NODE        DIRECTION      LOAD',/,
     1        ' NUMBER',19X,'MAGNITUD')
 2010 FORMAT(' ',I6,9X,I4,7X,E12.5)
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE IDRESP                                                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE IDRESP(NGPRES,NSVRES,NOUTEL,IDOEL,IDGP,IDSV,IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION IDGP(NGPRES),IDSV(NSVRES),IDOEL(NOUTEL)
C
C     READ AND GENERATE NODAL POINT DATA
C
      DO L=1,NOUTEL
        READ (IIN,     *) IDOEL(L)
        WRITE(IOUT,1010 ) IDOEL(L)
      END DO
C
      DO I=1,NSVRES
        READ (IIN,     *) IDSV(I)
        WRITE(IOUT,1010 ) IDSV(I)
      END DO
C
      DO J=1,NGPRES
        READ (IIN,     *) IDGP(J)
        WRITE(IOUT,1010 ) IDGP(J)
      END DO
C
 1010 FORMAT(I5)
C
      RETURN
C
      END
C
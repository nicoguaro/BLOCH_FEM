CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C---------P R O G R A M   D A M I A N  P E R I O D I C---------------- C
C                     Finite Element Method                            C
C                                                                      C
C                                                                      C
C VARIABLE NUMBER OF DOF PER NODE                                      C
C DYNAMIC MEMORY ALLOCATION                                            C
C                                                                      C
C UNIVERSIDAD EAFIT                                                    C
C GRUPO DE INVESTIGACION EN MECANICA APLICADA                          C
C AUGUST 2012                                                          C
C                                                                      C
C UNIX VERSION                                                         C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C**********************************************************************C
C                                                                      C
C                      MAIN PROGRAM STARTS                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      PROGRAM PERIODIC
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      ALLOCATABLE SKG(:,:),SMG(:,:),UG(:),COORD(:,:),AMATE(:,:)
C
      ALLOCATABLE ID(:,:),IELT(:),IELCON(:,:),NNE(:),MATP(:),NMATP(:),
     1            LM(:,:),IBC(:,:),NDOFN(:),NDOFEL(:)
      ALLOCATABLE IMNODES(:), IRNODES_WO(:), IRNODES(:), IMDOF(:),
     1            IRDOF(:), IRDOF_WO(:), IDOFCOOR(:), EIGVALS(:,:)

C
      PARAMETER (ZERO=0.D0,MXLC=1)
C
      CHARACTER*80 HED
      CHARACTER*10 FILENAME
C
C     *****************************************************************C
C     ***                P R O B L E M  F I L E S                    **C
C     *****************************************************************C
C
      IIN =7
      IOUT=8
      IMSG=9
      IFRQ =10
      WRITE(*,*) 'INPUT THE JOB NAME(max 10 characters): '
      READ(*,*) FILENAME
      LST=LEN_TRIM(FILENAME)
      OPEN(UNIT=IIN,FILE=FILENAME(1:LST)//".inp",FORM='FORMATTED')
      OPEN(UNIT=IOUT,FILE=FILENAME(1:LST)//".dat",FORM='FORMATTED')
      OPEN(UNIT=IMSG,FILE=FILENAME(1:LST)//".msg",FORM='FORMATTED')
      OPEN(UNIT=IFRQ,FILE=FILENAME(1:LST)//".fre",FORM='FORMATTED')
C
C     *****************************************************************C
C     ***                 I N P U T   P H A S E                      **C
C     *****************************************************************C
C
C     PROBLEM DEFINITION PARAMETERS
C
      READ(IIN,*,IOSTAT=INOUTSTATUS) HED
      IF(INOUTSTATUS.LT.0) STOP "***COULD NOT OPEN FILE"
      READ(IIN,     *) NUMNP,NUMEL,NUMAT,NDOFDIM,NMNE,NMDOFEL,NMPR,MCRD
      WRITE(IOUT,1900) HED
      WRITE(IOUT,2000) NUMNP,NUMEL,NUMAT,NDOFDIM,NMNE,NMDOFEL,NMPR,MCRD
C
      MXEL=NUMEL
      MXNO=NUMNP
      MXMA=NUMAT
      MXDOFDIM=NDOFDIM
      MXNE=NMNE
      MXDOFEL=NMDOFEL
C
C     INTEGER ARRAYS
C
      ALLOCATE(ID(MXDOFDIM,NUMNP),IELT(NUMEL),IELCON(MXNE,NUMEL),
     1         NNE(NUMEL),MATP(NUMEL),NMATP(NUMAT),LM(NUMEL,MXDOFEL),
     2         IBC(MXDOFDIM,NUMNP),NDOFN(NUMNP),NDOFEL(NUMEL),
     3         COORD(MCRD,NUMNP))
C
      CALL CLEARIM(ID, MXDOFDIM,NUMNP)
      CALL CLEARIM(IBC,MXDOFDIM,NUMNP)
      CALL CLEARIV(NDOFN,NUMNP)
      CALL CLEAR(COORD,MCRD,NUMNP)
C
      CALL NODINP(NUMNP,ID,IBC,NDOFDIM,MCRD,COORD,NDOFN,NEQ,IIN,IOUT)
      MXEQ=NEQ
C
C     SOLUTION ARRAYS
C
      ALLOCATE(SKG(MXEQ,MXEQ),SMG(MXEQ,MXEQ),UG(MXEQ),AMATE(NMPR,NUMAT))
C
C     CLEARS MODEL STORAGE
C
      CALL CLEARIM(IELCON,NMNE,NUMEL)
      CALL CLEARIM(LM,NUMEL,NMDOFEL)
      CALL CLEARIV(IELT,NUMEL)
      CALL CLEARIV(NNE,NUMEL)
      CALL CLEARIV(NDOFEL,NUMEL)
      CALL CLEAR(AMATE,NMPR,NUMAT)
      CALL CLEARIV(NMATP,NUMAT)
C
C     CLEARS SOLUTION STORAGE
C
      CALL CLEAR(SKG,MXEQ,MXEQ)
      CALL CLEAR(SMG,MXEQ,MXEQ)
      CALL CLEARV(UG,MXEQ)
C
C     READS MODEL
C
      CALL MATINP(NUMAT,NMPR,AMATE,NMATP,IIN,IOUT)
      CALL ELEINP(NUMNP,NUMEL,NNE,IELT,NDOFEL,NMNE,MATP,IELCON,IIN,IOUT)
C
C     READS BLOCH INFO
C
      READ(IIN,*) AKXMIN,AKYMIN,AKXMAX,AKYMAX,NKX,NKY,NCOND_WO,
     1            NCOND,NEVALS

      ALLOCATE(IRNODES_WO(NCOND_WO),IMNODES(NCOND),IRNODES(NCOND))

      READ(IIN,*) (IMNODES(I),I=1,NCOND)
      READ(IIN,*) (IRNODES(I),I=1,NCOND)
      READ(IIN,*) (IRNODES_WO(I),I=1,NCOND_WO)

C
C     *****************************************************************C
C     ***          A S  S E M B L Y   P H A S E                      **C
C     *****************************************************************C
C
      CALL ASSEMLIS(NUMNP,NUMEL,NMNE,NMDOFEL,NDOFDIM,NNE,NDOFN,IELCON,
     1              LM,ID,IIN,IOUT)

      CLOSE(IIN)
C
      CALL GSTFASEM(NUMNP,NUMEL,NUMAT,NNE,NMNE,NMDOFEL,IELT,IELCON,
     1              NDOFN,NDOFEL,MATP,NMATP,NMPR,MCRD,AMATE,
     2              COORD,LM,NEQ,SKG,SMG,IATYPE,IOUT)

      WRITE(IMSG, *) SKG


      IACU = 0
      DO I=1,MXDOFDIM
        DO  J=1,NCOND
            IF (ID(I,IMNODES(J)).NE.0)THEN
              IACU = IACU + 1
            END IF
        END DO
      END DO
      NCON_DOF = IACU

      IACU = 0
      DO I=1,MXDOFDIM
        DO  J=1,NCOND_WO
            IF (ID(I,IRNODES_WO(J)).NE.0)THEN
              IACU = IACU + 1
            END IF
        END DO
      END DO
      NCON_DOF_WO = IACU

      ALLOCATE(IRDOF_WO(NCON_DOF_WO),IMDOF(NCON_DOF),IRDOF(NCON_DOF),
     1         IDOFCOOR(NEQ),EIGVALS(NKX*NKY,NEVALS+2))

      CALL NOD2DOF(NUMNP,MXDOFDIM,NEQ,ID,NCOND,NCOND_WO,NCON_DOF,
     1             NCON_DOF_WO,IMNODES,IRNODES,IRNODES_WO,IMDOF,
     2             IRDOF,IRDOF_WO,IDOFCOOR)

C
C     *****************************************************************C
C     ***                 S L N  P H A S E                           **C
C     *****************************************************************C

      CALL BLOCH(NEQ,SKG,SMG,NUMNP,COORD,IDOFCOOR,NCON_DOF_WO,NCON_DOF,
     1           IRDOF_WO,IMDOF,IRDOF,NEVALS,AKXMIN,
     2           AKYMIN,AKXMAX,AKYMAX,NKX,NKY,EIGVALS)

      DO II=1,NKX*NKY
        DO JJ=1,NEVALS+2
          WRITE(IFRQ,*), EIGVALS(II,JJ)
        END DO
      END DO

C
C      WRITE(IMSG,3050)
      STOP

 1900 FORMAT(//,10X,'P R O B L E M   N A M E',10X,A80,//)
 2000 FORMAT (///
     1    ' C O N T R O L  I N F O R M A T I O N',//,
     2    '   NUMBER OF NODAL POINTS',10(' .'),' (NUMNP)=',I5,//,
     3    '   NUMBER OF ELEMENTS',12(' .'),'(NUMEL)=',I5,//,
     4    '   NUMBER OF MATERIALS       ',9(' .'),'(NUMAT)=',I5,//,
     5    '   DEGREE OF FREEDOM DIMENSION',6(' .'),'(NDOFDIM)=',I5,//,
     6    '   MAX.NUMBER OF NODES IN AN ELE.',6(' .'),'(NMNE)=',I5,//,
     7    '   MAX.NUMBER OF DOF IN AN ELE.',6(' .'),'(NMDOFEL)=',I5,//,
     8    '   NUMBER MAX OF MATE,R PROP ',6(' .'),'(NMPR)=',I5,//,
     9    '   PROBLEM DIMENSION ',6(' .'),'(MCRD)=',I5)
C
 3050 FORMAT(/,'PERIODIC JOB=',A8,1X,'COMPLETED')
C
      END PROGRAM

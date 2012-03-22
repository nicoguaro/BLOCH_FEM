CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   -----------------P R O G R A M   M Y F E A P--------------- -----  C
C                       Finite Element Method                          C
C                                                                      C
C GENERALIZED INCREMENTAL PLANE STRAIN ANALYSIS                        C
C                                                                      C
C VARIABLE NUMBER OF DOF PER NODE                                      C
C DYNAMIC MEMORY ALLOCATION                                            C
C OUTPUT IN TERMS OF SVARS                                             C
C SKYLINE SOLVER                                                       C
C                                                                      C
C UNIVERSIDAD EAFIT                                                    C
C APPLIED MECHANICS RESEARCH GROUP                                     C
C SEPTREMBER 2011                                                      C
C                                                                      C
C (Dynamic memory allocation+general element)                          C
C  PC VERSION                                                          C
C                                                                      C
C                                                                      C
C         PROBLEM PARAMETERS                                           C
C         ------------------                                           C
C                                                                      C
C     NUMNP      :NUMBER OF NODAL POINTS                               C
C     NUMEL      :NUMBER OF ELEMENTS                                   C
C     NUMAT      :NUMBER OF MATERIAL PROPERTIES                        C
C     NINCR      :NUMBER OF INCREMENTS (TIME STEPS)                    C
C     NDOFDIM    :PROBLEM DIMENSION (2D, 3D)                           C
C     NMNE       :MAXIMUM NUMBER OF NODES IN A GIVEN ELEMENT           C
C     NMDOFEL    :MAXIMUM NUMBER OF DOF ASSIGNED TO A GIVEN ELEMENT    C
C     NMPR       :MAXIMUM NUMBER OF MATERIAL PROPERTIES IN A PROFILE   C
C     NSVARS     :MAXIMUM NUMBER OF SVARS IN A GIVEN ELEMENT           C
C     IPROB      :OUTPUT MODE                                          C
C     NSVRES     :SVARS FOR OUTPUT                                     C
C     NGPRES     :GAUSS POINTS FOR OUTPUT                              C
C     NOUTEL     :ELEMENTS FOR OUTPUT                                  C
C                                                                      C
C                                                                      C
C         PROBLEM REAL ARRAYS                                          C
C         -------------------                                          C
C                                                                      C
C     SKG(:)          :GLOBAL STIFFNESS MATRIX IN VECTOR STORAGE       C
C     RHSG(:,:)       :CURRENT RIGHT HAND SIDE VECTROR                 C
C     RH(:,:)         :NODAL FORCES HISTORY                            C
C     UH(:,:)         :NODAL DISPLACEMENTS HISTORY                     C
C     UG(:)           :CURRENT NODAL DISPLACEMENTS                     C
C     DU(:)           :NODAL INCREMENTAL DISPLACEMENTS                 C
C     DP(:)           :NODAL INCREMENTAL NODAL FORCES                  C
C     DPT(:)          :                                                C
C     SVARSEGPT(:,:)  :CURRENT SVARS VALUES                            C
C     SVARSGPTH(:,:,:):SVARS HISTORIES                                 C
C     COORD(:,:)      :NODAL COORDINATES                               C
C     AMATE(:,:)      :MATERIAL PROPERTIES                             C
C                                                                      C
C         PROBLEM INTEGER ARRAYS                                       C
C         ----------------------                                       C
C                                                                      C
C     ID(:,:)         :EQUATION ASSIGNED TO EACH NODE                  C
C     IELT(:)         :ELEMENT TYPE ID'S                               C
C     IELCON(:,:)     :ELEMENT CONNECTIVITIES                          C
C     NNE(:)          :NUMBER OF NODES ASSIGNED TO EACH ELEMENT        C
C     MATP(:)         :MATERIAL PROFILE ASSIGNED TO EACH ELEMENT       C
C     NMATP(:)        :NUMBER OF MATERIAL PARAMTERS ASSIGNED TO EACH   C
C                      PROFILE                                         C
C     LM(:,:)         :DME OPERATOR                                    C
C     NDOFN(:)        :NUMBER OF DOF ASSIGNED TO EACH NODE             C
C     NDOFEL(:)       :NUMBER OF DOF ASSIGNED TO EACH ELEMENT          C
C     IDGP(:)         :ID OF GAUSS POINTS FOR OUTPUT                   C
C     IDSV(:)         :ID OF SVARS FOR OUTPUT                          C
C     IDOEL(:)        :ID OF ELEMENTS FOR OUTPUT                       C
C     IMIS(:)         :                                                C
C     ICHS(:)         :COLUMN HEIGHTS                                  C
C     MAXA(:)         :ADDRESS OF DIAGONAL ELEMENTS INDEX              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
***********************************************************************C
C                                                                      C
C                      MAIN PROGRAM STARTS                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      PROGRAM MY_FEAP_BLOCH
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C      INCLUDE 'link_fnl_static_hpc.h'
C
      ALLOCATABLE SKG(:,:),SMG(:,:),COORD(:,:),AMATE(:,:)
C
      ALLOCATABLE ID(:,:),IELT(:),IELCON(:,:),NNE(:),MATP(:),NMATP(:),
     1            LM(:,:),NDOFN(:),NDOFEL(:)
C
      COMPLEX*16 CALPHA
      ALLOCATABLE CALPHA(:),RBETA(:), EIGVALS(:,:)
      ALLOCATABLE IRNODES_WO(:), IMNODES(:), IRNODES(:)
C
      PARAMETER (ZERO=0.D0,MAXITER=10,TOL=1.0D-7, PI=3.141592653589793)
C
      CHARACTER*80 HED
      CHARACTER*10 FILENAME
C
C     *****************************************************************C
C     ***          P R O B L E M  F I L E S                          **C
C     *****************************************************************C
C
      IIN=5
      IOUT=6
      IMSG=7
      IST=8
      IFRQ=9
      WRITE(*,*) 'INPUT THE JOB NAME(max 10 characters):'
      READ(*,*) FILENAME
      LST=LEN_TRIM(FILENAME)
      OPEN(UNIT=IIN,FILE =FILENAME(1:LST)//".inp",FORM='FORMATTED')
      OPEN(UNIT=IOUT,FILE=FILENAME(1:LST)//".dat",FORM='FORMATTED')
      OPEN(UNIT=IMSG,FILE=FILENAME(1:LST)//".msg",FORM='FORMATTED')
      OPEN(UNIT=IST ,FILE=FILENAME(1:LST)//".sta",FORM='FORMATTED')


      OPEN(UNIT=IFRQ,FILE =FILENAME(1:LST)//".freq")
C
C     *****************************************************************C
C     ***                 I N P U T   P H A S E                      **C
C     *****************************************************************C
C
C     PROBLEM DEFINITION PARAMETERS
C
      READ(IIN,*,IOSTAT=INOUTSTATUS) HED
      IF(INOUTSTATUS.LT.0) STOP "***COULD NOT OPEN FILE"
      READ(IIN,     *) NUMNP,NUMEL,NUMAT,NINCR,NDOFDIM,NMNE,NMDOFEL,
     1                 NMPR
      WRITE(IMSG,1900) HED
      WRITE(IMSG,2000) NUMNP,NUMEL,NUMAT,NINCR,NDOFDIM,NMNE,
     1                 NMDOFEL
      MXEL=NUMEL
      MXNO=NUMNP
      MXMA=NUMAT
      MXINC=NINCR
      MXDOFDIM=NDOFDIM
      MXNE=NMNE
      MXDOFEL=NMDOFEL



C
C     INTEGER ARRAYS
C
      ALLOCATE(ID(NDOFDIM,NUMNP),IELT(NUMEL),IELCON(MXNE,NUMEL),
     1         NNE(NUMEL),MATP(NUMEL),NMATP(NUMAT),LM(MXDOFEL,NUMEL),
     2         NDOFN(NUMNP),NDOFEL(NUMEL))
C
      ALLOCATE(COORD(NDOFDIM,NUMNP))
C
      CALL CLEARIM(ID,NDOFDIM,NUMNP)
      CALL CLEARIV(NDOFN,NUMNP)
      CALL CLEAR(COORD,NDOFDIM,NUMNP)
C
      CALL NODINP(NUMNP,ID,NDOFDIM,COORD,NDOFN,NEQ,IIN,IMSG)
      MXEQ=NEQ
      MNM=NEQ+1
C

C
C     STIFFNESS ASSEMBLY ARRAYS
C
      ALLOCATE(SKG(NEQ,NEQ),SMG(NEQ,NEQ))
C
C     SOLUTION ARRAYS
C
      ALLOCATE(AMATE(NMPR,NUMAT))
C
C     CLEARS MODEL STORAGE
C
      CALL CLEARIM(IELCON,NMNE,NUMEL)
      CALL CLEARIM(LM,NMDOFEL,NUMEL)
      CALL CLEARIV(IELT,NUMEL)
      CALL CLEARIV(NNE,NUMEL)
      CALL CLEARIV(NDOFEL,NUMEL)
      CALL CLEAR(AMATE,NMPR,NUMAT)
      CALL CLEARIV(NMATP,NUMAT)
      CALL CLEAR(SKG,NEQ,NEQ)
      CALL CLEAR(SMG,NEQ,NEQ)
C
C     READS IN MODEL
C
      CALL MATINP(NUMAT,NMPR,AMATE,NMATP,IIN,IMSG)
      CALL ELEINP(NUMNP,NUMEL,NNE,IELT,NDOFEL,NMNE,MATP,IELCON,IIN,IMSG)
C
C     CREATES DME OPERATOR LM(:) FOR
C     USE IN SKYLINE STORAGE FORM.
C
      CALL ASSEMLIS(NUMNP,NUMEL,NMNE,NMDOFEL,NDOFDIM,NNE,NDOFN,IELCON,
     1              LM,ID,IIN,IMSG)
C
C     READS EIGENANALYSIS PARAMETERS AND BLOCH CONDITIONS
C


      READ(IIN,*) AKXMIN,AKYMIN,AKXMAX,AKYMAX,NKX,NKY,NCOND_WO, 
     1            NCOND,NEVALS

      ALLOCATE(CALPHA(NEVALS),RBETA(NEVALS),EIGVALS(NKX*NKY,NEVALS+2),
     1        IRNODES_WO(NCOND_WO),IMNODES(NCOND),IRNODES(NCOND))

      READ(IIN,*) (IMNODES(I),I=1,NCOND)
      READ(IIN,*) (IRNODES(I),I=1,NCOND)
      READ(IIN,*) (IRNODES_WO(I),I=1,NCOND_WO)


C
C     *****************************************************************C
C     ***             S L N    P H A S E                             **C
C     *****************************************************************C
C
      WRITE(IST,4000)
C
C     CLEARS SOLUTION STORAGE
C
      CALL GSTFASEM(NUMNP,NUMEL,NUMAT,NNE,NMNE,MXDOFDIM,NMDOFEL,
     1              IELT,IELCON,NDOFN,NDOFEL,MATP,NMATP,NMPR,AMATE,
     2              COORD,LM,NEQ,SKG,SMG,IMSG,JJ)



      CALL BLOCH(NEQ,SKG,SMG,NUMNP,COORD,NCOND_WO,NCOND,
     1           IRNODES_WO,IMNODES,IRNODES,NEVALS,AKXMIN,AKYMIN,
     2           AKXMAX,AKYMAX,NKX,NKY,EIGVALS)

      DO II=1,NKX*NKY
        DO JJ=1,NEVALS+2
          WRITE(IFRQ,*), EIGVALS(II,JJ)
        END DO
      END DO

C
      STOP
C
C     *****************************************************************C
C     ***         S L N  O U T P U T  P H A S E                      **C
C     *****************************************************************C
C
 1900 FORMAT(//,10X,'P R O B L E M   N A M E',10X,A80,//)
 2000 FORMAT (///,
     1    ' C O N T R O L  I N F O R M A T I O N',//,
     2    '   NUMBER OF NODAL POINTS',10(' .'),' (NUMNP)=',I5,//,
     3    '   NUMBER OF MATERIALS',12(' .'),'(NUMEL)=',I5,//,
     4    '   NUMBER OF ELEMENTS       ',9(' .'),'(NUMAT)=',I5,//,
     5    '   NUMBER OF LOAD INCREMENTS',9(' .'),'(NINCR)=',I5,//,
     6    '   DEGREE OF FREEDOM DIMENSION',6(' .'),'(NDOFDIM)=',I,//,
     7    '   MAX.NUMBER OF NODES IN AN ELE.',6(' .'),'(NMNE)=',I2,//,
     8    '   MAX.NUMBER OF DOF IN AN ELE.',6(' .'),'(NMDOFEL)=',I2)
C
 3016 FORMAT(1X,'TIME STEP COMPLETED=',2X,F5.2,2X,'TOT.TIME COMPLETED
     1       =',2X,F5.2,/)
 3017 FORMAT(1X,'TOTAL ITERATIONS=',I5,/)
 3020 FORMAT(5X,I5,10X,I5)
 4000 FORMAT(2X,'INC',6X,'TOT-TIME',6X,'INC-TIME')
 4010 FORMAT(2X,I3,3X,I3,7x,F6.3,9X,F6.3)
 3050 FORMAT(/,'MYABAQUS JOB=',A8,1X,'COMPLETED')
C
      END PROGRAM MY_FEAP_BLOCH
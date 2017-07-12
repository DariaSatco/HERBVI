      SUBROUTINE HERWIGINIT(IBEAM,PROCIN,ENERGY)
c     HERWIG/JIMMY COMMONS. HERWIG COMMON INCLUDES AN IMPLICIT NONE.
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'jimmy.inc'
      DOUBLE PRECISION ENERGY
      INTEGER PROCIN,IBEAM
      EXTERNAL HWUDAT
c     Set beam energies
      IF(IBEAM==0) THEN
         PART1='E+ '
         PART2='E- '
         ZMXISR=0.
      ELSE
         PART1='P  '
         PART2='P  '
      ENDIF
c     SUSY 2->2 with S.U.E suppressed
      IPROC=PROCIN

      PBEAM1=0.5*ENERGY
      PBEAM2=0.5*ENERGY

C     HERWIG initialization
      CALL HWIGIN
C************************inserted from test.f**************************
      LWSUD=77                                                           
      LRSUD=0                                                          
C     NOWGT=.FALSE.                                                     
      NOWGT=.true.                                                      
      MAXPR=0                    
      PRVTX = .FALSE.
C---TOP MASS                                                            
      RMASS(6)=172.5             
C---INPUT SUSY PARTICLE (AND TOP QUARK) DATA
c      CALL HWISSP                                                  
      NRN(1)= 431269609                                                 
      NRN(2)=2127819028                                               
C---HIGGS MASS                                                          
      RMASS(201)=125.      
C***********************************************************************

C     JIMMY initialization
      CALL JIMMIN

      MODPDF(1)=260000
      AUTPDF(1)='HWLHAPDF'
      MODPDF(2)=260000
      AUTPDF(2)='HWLHAPDF'

c     Minimum PT of the hardest scatter
      PTMIN=10.


c     Minimum PT of secondary scatters
C  IBEAM=0
      IF(IBEAM.NE.0) THEN
C         PTJIM=4.9
C         JMUEO = 1
C         JMRAD(73) = 1.8
C         PRSOF=0.
C         IPRINT = 2
         JMUEO=1
         PTJIM=2.8
C         PRRAD=1.8
         IPRINT=2

C---  Turn MI on(1) or off(0)
         MSFLAG=1
      ELSE
         MSFLAG=0
      ENDIF
C      IBEAM=1
      CALL HWUINC
      CALL HWEINI            
      CALL HVINIT

C******changed number of events MAXEV (orig.10000) and MAXER (orig. MAXEV/5)****
      MAXEV=10000
      MAXER=MAXEV/10

C     Initialise your histograms here.
      CALL HWABEG
      
      RETURN

 999    WRITE (6,*)
      WRITE (6,*) 'SUSY input file did not open correctly.'
      WRITE (6,*) 'Please check that it is in the right place.'
      WRITE (6,*) 'Examples can be obtained from the ISAWIG web page.'
      WRITE (6,*)

      STOP

      END

      SUBROUTINE HERWIGGEN
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'jimmy.inc'
      LOGICAL ABORT
C---INITIALISE EVENT
 10      CALL HWUINE
C---GENERATE HARD SUBPROCESS
      CALL HWEPRO
C---GENERATE PARTON CASCADES
      CALL HWBGEN
C---  GENERATE MULTIPARTON INTERACTIONS
      IF (MSFLAG.EQ.1) THEN
         CALL HWMSCT(ABORT)
         IF (ABORT) GOTO 10
      ENDIF
C---DO HEAVY OBJECT DECAYS
      CALL HWDHOB
C---DO CLUSTER FORMATION
      CALL HWCFOR
C---DO CLUSTER DECAYS
      CALL HWCDEC
C---DO UNSTABLE PARTICLE DECAYS
      CALL HWDHAD
C---DO HEAVY FLAVOUR HADRON DECAYS
      CALL HWDHVY
C---ADD SOFT UNDERLYING EVENT IF NEEDED
      CALL HWMEVT
C---FINISH EVENT
      CALL HWUFNE
      IF(IERROR.NE.0) GOTO 10
      CALL HWANAL
      END
      SUBROUTINE HERWIGEND
      INCLUDE 'HERWIG65.INC'
C---TERMINATE ELEMENTARY PROCESS
      CALL HWEFIN
C     Finish off JIMMY
      CALL JMEFIN
C---USER'S TERMINAL CALCULATIONS
      CALL HWAEND
      END
C----------------------------------------------------------------------
      SUBROUTINE HWABEG
C     USER'S ROUTINE FOR INITIALIZATION
C----------------------------------------------------------------------
      END
C----------------------------------------------------------------------
      SUBROUTINE HWAEND
C     USER'S ROUTINE FOR TERMINAL CALCULATIONS, HISTOGRAM OUTPUT, ETC
C----------------------------------------------------------------------
      END
C----------------------------------------------------------------------
      SUBROUTINE HWANAL
C     USER'S ROUTINE TO ANALYSE DATA FROM EVENT
C Analysis of each BVI event 
      INCLUDE 'HERWIG65.INC'
      INTEGER I,CHARGE
      DOUBLE PRECISION EPS,MOM(4)
      PARAMETER(EPS=1D-2)
      CHARGE = 0
      DO I=1,4
        MOM(I) = ZERO
      ENDDO
C--test charge and momentum conservation
      DO I=1,NHEP
        IF(ISTHEP(I).EQ.1) THEN
          CHARGE = CHARGE+ICHRG(IDHW(I))
          CALL HWVSUM(4,PHEP(1,I),MOM,MOM)
        ENDIF
      ENDDO             
C--test charge conservation  
      IF(CHARGE-ICHRG(IDHW(1))-ICHRG(IDHW(2)).NE.0) 
     &   print *,'warning violates charge conservation',nevhep                
C--test momentum conservation
      DO I=1,2
        CALL HWVDIF(4,MOM,PHEP(1,I),MOM)
      ENDDO
      DO I=1,4
        IF(ABS(MOM(I)).GT.EPS) 
     &    print *,'warning violates momentum conservation'
     &            ,nevhep,i,mom(i)
      ENDDO
C----------------------------------------------------------------------
      END
C----------------------------------------------------------------------

      SUBROUTINE CROSSSECTION(OUT)
      INCLUDE 'HERWIG65.INC'
      DOUBLE PRECISION OUT
      OUT=AVWGT
      END

C-----------------------------------------------------------------------
CDECK  ID>,  UPINIT.
*CMZ :-        -16/07/02  10.30.48  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
C      SUBROUTINE UPINIT
C-----------------------------------------------------------------------
C     DUMMY UPINIT ROUTINE DELETE AND REPLACE IF USING LES HOUCHES
C     INTERFACE
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C      WRITE (6,10)
C 10    FORMAT(/10X,'UPINIT CALLED BUT NOT LINKED')
C      STOP
C      END
C-----------------------------------------------------------------------
CDECK  ID>,  UPEVNT.
*CMZ :-        -16/07/02  10.30.48  by  Peter Richardson
*-- Author :    Peter Richardson
C-----------------------------------------------------------------------
C      SUBROUTINE UPEVNT
C-----------------------------------------------------------------------
C     DUMMY UPEVNT ROUTINE DELETE AND REPLACE IF USING LES HOUCHES
C     INTERFACE
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C      WRITE (6,10)
C 10    FORMAT(/10X,'UPEVNT CALLED BUT NOT LINKED')
C      STOP
C      END

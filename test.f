                                 
C                                                                       
      PROGRAM HWIGPR                                                    
                                                                        
C                                                                       
C An example test cradle for HERBVI use in HERWIG                       
C                                                                       
      INCLUDE 'HERWIG65.INC'   
      INCLUDE 'HERBVI.INC'
      LOGICAL HVCA70,HVCA71                                             
      COMMON /HVCETD/ HVCA70,HVCA71                                     
      EXTERNAL HWUDAT       
      LOGICAL ERRFLG                                                    
      COMMON /ERRCOM/ ERRFLG                                            
      INTEGER N,A,B,C                                                   
C---Error count parameter                                               
      A=0                                                               
C---Counts of event types                                               
      B=0                                                               
      C=0                                                               
C---MAX NUMBER OF EVENTS THIS RUN                                       
      MAXEV=1
C---PROCESS (SEE TABLE)                                                 
      IPROC=7200                                                       
      PBEAM1=20000.                                                     
      PBEAM2=20000.                                                     
      PART1='P   '                                                      
      PART2='P   '                                                      
C---INITIALISE OTHER COMMON BLOCKS 
      CALL HWIGIN                                                       
C---USER CAN RESET PARAMETERS AT                                        
C   THIS POINT, OTHERWISE VALUES                                        
C   SET IN HWIGIN WILL BE USED.                                         
      MAXER=MAXEV/5                                                          
      LWSUD=77                                                           
      LRSUD=0                                                          
C     NOWGT=.FALSE.                                                     
      NOWGT=.true.                                                      
      MAXPR=1                    
      PRVTX = .FALSE.
C---TOP MASS                                                            
      RMASS(6)=175.             
C---INPUT SUSY PARTICLE (AND TOP QUARK) DATA
c      CALL HWISSP                                                  
      NRN(1)= 431269609                                                 
      NRN(2)=2127819028                                               
C---HIGGS MASS                                                          
      RMASS(201)=300.                                                   
C---COMPUTE PARAMETER-DEPENDENT CONSTANTS  
      CALL HWUINC                                                       
C---CALL HWUSTA TO MAKE ANY PARTICLE STABLE
c      CALL HWUSTA('PI0 ')                                              
C---Call HERBVI initialisation                                          
      CALL HVINIT                                                       
CC---Set not 1/shat generation                                          
C      HVCONT(4)=.FALSE.                                                
CC---Set LOME boson number estimate                                     
C      HVCONT(7)=.FALSE.                                                
CC---Set limits for Shat generation                                     
C      SLMIN=2.D0*LOG(5.D0/17.0D0)                                      
C---USER'S INITIAL CALCULATIONS                                         
      CALL HWABEG                                                       
C---INITIALISE ELEMENTARY PROCESS                                       
      CALL HWEINI                                                       
C---LOOP OVER EVENTS                                                    
      DO 100 N=1,MAXEV                                                  
C---INITIALISE EVENT                                                    
      CALL HWUINE                                                       
C---GENERATE HARD SUBPROCESS                                            
      CALL HWEPRO                                                       
C---Add event types                                                     
      IF (HVCA70) B=B+1                                                 
      IF (HVCA71) C=C+1                                                 
C---GENERATE PARTON CASCADES                                            
      CALL HWBGEN                                                       
C---Decay heavy objects                                            
      CALL HWDHOB                                                       
C---DO CLUSTER HADRONIZATION                                            
      CALL HWCFOR                                                       
C---DO CLUSTER DECAY                                                    
      CALL HWCDEC                                                       
C---DO UNSTABLE PARTICLE DECAYS                                         
      CALL HWDHAD                                                       
C---DO HEAVY FLAVOUR DECAYS                                             
      CALL HWDHVY                                                       
C---ADD SOFT UNDERLYING EVENT IF NEEDED                                 
      CALL HWMEVT                                                       
C---End event                                                           
      CALL HWUFNE                                                       
C---USER'S EVENT ANALYSIS                                               
      CALL HWANAL                                                       
 100   CONTINUE                                                          
C---TERMINATE ELEMENTARY PROCESS                                        
      CALL HWEFIN                                                       
C---USER'S TERMINAL CALCULATIONS                                        
      PRINT *,'  '                                                      
      PRINT *,' Number of HVCBVI errors ',NUMERU                        
      PRINT *,' Number of events        ',MAXEV                         
      PRINT *,' Percntage of C BVI rej. ',DBLE(100*NUMERU)/DBLE(MAXEV)  
      PRINT *,' Number of BVI events ',B                                
      PRINT *,' Number of m-w events ',C                                
      PRINT *,'  '                                                      
      CALL HWAEND                                                       
 999   CONTINUE                                                          
      STOP                                                              
      END                                                               
C-----------------------------------------------------------------------
      SUBROUTINE HWABEG                                                 
                                                                        
C Initial calculations before the run starts                            
                                                                        
      END                                                               
C-----------------------------------------------------------------------
      SUBROUTINE HWAEND                                                 
                                                                        
C Final analysis after all events generated                             
                                                                        
      END                                                               
C-----------------------------------------------------------------------
      SUBROUTINE HWANAL                                                 
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
      END


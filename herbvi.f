C   HERBVI.  HERBVI, A PROGRAM FOR SIMULATION OF BARYON AND LEPTON      
C   NUMBER VIOLATING PROCESSES.  M.J. GIBBS, B.R. WEBBER.               
C     REF. IN COMP. PHYS. COMMUN. 90 (1995) 369                         
C                                                                       
C                                                                       
C                                                                       
C This is the code for saddle10.for                                     
C                                                                       
C                                                                       
C *******************************************************************   
C                                                                       
C                                                                       
C *******************************************************************   
C                                                                       
C                                                                       
C *******************************************************************   
C                                                                       
C                                                                       
C 17/10/92 Saddle point routines for BVI phase space calcs.             
C                                                                       
C Routines included:                                                    
C                                                                       
C BVPABR S  Asymptotic expansion of modified bessel fn. ratios          
C BVPAIL F  Logarithim of I1(x) modified bessel function                
C BVPAIR S  Ratio of I0/I1 and first three derivatives                  
C BVPAI0 F  Modified bessel function I0                                 
C BVPAI1 F  Modified bessel function I1                                 
C BVPAKL F  Logarithim of K1(x) modified bessel function                
C BVPAKR S  Ratio of K0/K1 and first three derivatives                  
C BVPAKQ F  Pick member of quark doublet using ps ratio                 
C BVPAK0 F  Modified bessel function K0                                 
C BVPAK1 F  Modified bessel function K1                                 
C BVPALF F  Log of Gamma Function                                       
C BVPCAL S  Matrix-Element contrib. for Event given mom. distrib.       
C BVPDIN S  Initialize physical parameters for boson MC                 
C BVPINI S  Initialization routine                                      
C BVPINS F  Instanton contribution to cross section                     
C BVPINT S  Phase space integral                                        
C BVPMAM S  MAMBO phase-space generating subroutine                     
C BVPMBF FH Search in Boson number for limits                           
C BVPMBG SH Calculate boson number contributions                        
C BVPMBI S  Initialize arrays for boson/fermion MC routines             
C BVPMBP S  Find peak boson number for given energy                     
C BVPMBS SH Calaulate boson contributions within limits                 
C BVPMBT S  Total up cross section contributions                        
C BVPMCI S  Initialize MC routine (MAMBO) and momentum array            
C BVPMFC SH Initialize configuration arrays for given event structure   
C BVPMFI S  Initialize CONFIG array for BVI event                       
C BVPMFM S  MC calculation of BVI cs for given number of particles      
C BVPMML S  Monte-Carlo iteration routine (uses MAMBO ps-generator)     
C BVPNBA F  Evaluate binomial expression                                
C BVPNBS S  Select total number of bosons for BVI process               
C BVPNGA S  Gaussian approximation for number of gammas                 
C BVPSBR S  Modified bessel function ratios                             
C BVPSCB S  Initialize all common blocks for given beta value           
C BVPSNR S  Newton-Rhapson iteration for saddle point                   
C BVPSPD S  N'th derivative of integrand                                
C BVPSPS S  Saddle point integration for phase space integral           
C BVPSSC S  Higher order corrections to saddle point integral           
C BVPSZD S  Integrand for saddle point integral                         
C CGEVPB F  Conversion factor for Gev-2 to pb                           
C                                                                       
C ********************************************************************* 
C                                                                       
C Common blocks in SADDLE10.INC:                                        
C                                                                       
C BVPCFG Particle configuration                                         
C BVPCBD Internal bessel function data                                  
C BVPCSD Saddle point data                                              
C BVPCDR Derivatives and Integral results                               
C BVPCBS Internal variables for Boson MC routines BVPMB*                
C BVPCBP Physical paramaters for Boson MC. Set by BVPDIN                
C                                                                       
C ********************************************************************* 
C                                                                       
C Modifications:                                                        
C 03/11/92 Bessel function routines from MAMBO program added            
C 13/02/93 All dependence on NAG now removed                            
C 16/02/93 Boson MC routines added                                      
C 07/04/93 Boson-Fermion MC routines updated                            
C 12/04/93 New param common block added; routines streamlined           
C 11/10/93 INCLUDE files modified for UNIX preprocessing                
C 26/03/94 BVPAKQ routine added                                         
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPABR(X,R0,R1,R2,R3)                                  
C                                                                       
C ASYMPTOTIC EXPANSION OF I0(X)/I1(X) AND ITS FIRST 3 DERIVATIVES       
C                                                                       
C Taken from Kleiss/Stirling MAMBO program                              
C Calculates ratio for X values >10 and corresponding K                 
C ratios using K0(x)/K1(x)=I0(-x)/I1(-x)                                
C                                                                       
C Modified by M Gibbs 15/6/92                                           
C                                                                       
      DOUBLE PRECISION A1(10),A2(10),A3(10),A0(10),                     
     &                 X,Y,R0,R1,R2,R3                                  
      DATA A0/.5000000000D+00,.3750000000D+00,                          
     .        .3750000000D+00,.4921875000D+00,.8437500000D+00,          
     .        .1854492188D+01,.5062500000D+01,.1658578491D+02,          
     .        .6333398438D+02,.2756161079D+03/                          
      DATA A1/.5000000000D+00,.7500000000D+00,                          
     .        .1125000000D+01,.1968750000D+01,.4218750000D+01,          
     .        .1112695313D+02,.3543750000D+02,.1326862793D+03,          
     .        .5700058594D+03,.2756161079D+04/                          
      DATA A2/.1000000000D+01,.2250000000D+01,                          
     .        .4500000000D+01,.9843750000D+01,.2531250000D+02,          
     .        .7788867188D+02,.2835000000D+03,.1194176514D+04,          
     .        .5700058594D+04,.3031777187D+05/                          
      DATA A3/.3000000000D+01,.9000000000D+01,                          
     .        .2250000000D+02,.5906250000D+02,.1771875000D+03,          
     .        .6231093750D+03,.2551500000D+04,.1194176514D+05,          
     .        .6270064453D+05,.3638132625D+06/                          
      IF(DABS(X).LT.10D0) THEN                                          
        WRITE(*,*) 'ROUTINE BVPABR: MAY BE INACCURATE FOR SMALL X=',X   
      ENDIF                                                             
      Y=1./X                                                            
      R0=1.+Y*                                                          
     .     (A0(1)+Y*(A0(2)+Y*(A0(3)+Y*(A0(4)+Y*(A0(5)+Y*(               
     .      A0(6)+Y*(A0(7)+Y*(A0(8)+Y*(A0(9)+Y*A0(10))))))))))          
      R1=-Y*Y*                                                          
     .     (A1(1)+Y*(A1(2)+Y*(A1(3)+Y*(A1(4)+Y*(A1(5)+Y*(               
     .      A1(6)+Y*(A1(7)+Y*(A1(8)+Y*(A1(9)+Y*A1(10))))))))))          
      R2=Y*Y*Y*                                                         
     .     (A2(1)+Y*(A2(2)+Y*(A2(3)+Y*(A2(4)+Y*(A2(5)+Y*(               
     .      A2(6)+Y*(A2(7)+Y*(A2(8)+Y*(A2(9)+Y*A2(10))))))))))          
      R3=-Y*Y*Y*Y*                                                      
     .     (A3(1)+Y*(A3(2)+Y*(A3(3)+Y*(A3(4)+Y*(A3(5)+Y*(               
     .      A3(6)+Y*(A3(7)+Y*(A3(8)+Y*(A3(9)+Y*A3(10))))))))))          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPAIL(X)                                                
C                                                                       
C Modified bessel function of the first kind Log(I1(x))                 
C                                                                       
C Taken from MAMBO by M Gibbs 03/11/92                                  
C                                                                       
      DOUBLE PRECISION BVPAIL,X,T,BETEMP                                
      IF(X.LT.0D0) THEN                                                 
        STOP 'BVPAIL: X<0'                                              
      ELSEIF(X.LE.3.75D0) THEN                                          
        T=(X/3.75D0)**2                                                 
        BETEMP=0.5D0+T*(0.87890594+T*(0.51498869+T*(0.15084934+         
     .  T*(0.02658733+T*(0.00301532+T*0.00032411)))))                   
        BVPAIL=DLOG(X)+DLOG(BETEMP)                                     
      ELSE                                                              
        T=3.75D0/X                                                      
        BETEMP=0.39894228+T*(-0.03988024+T*(-0.00362018+T*(             
     .         0.00163801+T*(-0.01031555+T*(0.02282967+T*(              
     .        -0.02895312+T*(0.01787654-T*0.00420059)))))))             
        BVPAIL=X-DLOG(X)/2.+DLOG(BETEMP)                                
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPAIR(X,FI0,FI1,FI2,FI3)                              
C                                                                       
C I0(x)/I1(x) and the first three derivatives                           
C                                                                       
C Routines used:BVPAI0    Modified bessel function I0                   
C               BVPAI1    Modified bessel function I1                   
C               BVPABR    Modified bessel fn. asymptotic expansion      
C                                                                       
C Taken from MAMBO program by M Gibbs 03/11/92                          
C                                                                       
      DOUBLE PRECISION X,FI0,FI1,FI2,FI3,R1,R2,R3,R4,X2,                
     &             BVPAI1,BVPAI0                                        
      IF(X.LE.0D0) THEN                                                 
        WRITE(*,*) ' ROUTINE BVPAIR: X=',X,' NOT ALLOWED'               
        STOP                                                            
      ENDIF                                                             
      IF(X.LT.10D0) THEN                                                
        R1=BVPAI0(X)/BVPAI1(X)                                          
        R2=R1*R1                                                        
        R3=R2*R1                                                        
        R4=R3*R1                                                        
        X2=X*X                                                          
        FI0=R1                                                          
        FI1=( - R2*X + R1 + X )/X                                       
        FI2=( 2.*R3*X - 3.*R2 - 2.*R1*X + 1. )/X                        
        FI3=( 6.*R4*X2 - 12.*R3*X + R2*(3.-8.*X2) +                     
     .        8.*R1*X + 2.*X2 + 1. )/(-X2)                              
      ELSE                                                              
        CALL BVPABR(X,FI0,FI1,FI2,FI3)                                  
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPAI0(X)                                                
C                                                                       
C Modified bessel function of the first kind I0(x)                      
C                                                                       
C Taken from MAMBO by M Gibbs 03/11/92                                  
C                                                                       
      DOUBLE PRECISION X,BVPAI0,T,BETEMP                                
      IF(X.LT.0D0) THEN                                                 
        STOP 'BVPAI0: X<0'                                              
      ELSEIF(X.LE.3.75D0) THEN                                          
        T=(X/3.75D0)**2                                                 
        BVPAI0=1D0+T*(3.5156229+T*(3.0899424+T*(1.2067492+              
     .  T*(0.2659732+T*(0.0360768+T*0.0045813)))))                      
      ELSE                                                              
        T=3.75D0/X                                                      
        BETEMP=0.39894228+T*(0.01328592+T*(0.00225319+T*(               
     .  -0.00157565+T*(0.00916281+T*(-0.02057706+T*(                    
     .   0.02635537+T*(-0.01647633+T*0.00392377)))))))                  
        BVPAI0=DEXP(X)/DSQRT(X)*BETEMP                                  
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPAI1(X)                                                
C                                                                       
C Modified bessel function of the first kind I1(x)                      
C                                                                       
C Taken from MAMBO by M Gibbs 03/11/92                                  
C                                                                       
      DOUBLE PRECISION X,BVPAI1,T,BETEMP                                
      IF(X.LT.0D0) THEN                                                 
        STOP 'BVPAI1: X<0'                                              
      ELSEIF(X.LE.3.75D0) THEN                                          
        T=(X/3.75D0)**2                                                 
        BETEMP=0.5D0+T*(0.87890594+T*(0.51498869+T*(0.15084934+         
     .  T*(0.02658733+T*(0.00301532+T*0.00032411)))))                   
        BVPAI1=X*BETEMP                                                 
      ELSE                                                              
        T=3.75D0/X                                                      
        BETEMP=0.39894228+T*(-0.03988024+T*(-0.00362018+T*(             
     .         0.00163801+T*(-0.01031555+T*(0.02282967+T*(              
     .        -0.02895312+T*(0.01787654-T*0.00420059)))))))             
        BVPAI1=DEXP(X)/DSQRT(X)*BETEMP                                  
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPAKL(X)                                                
C                                                                       
C Modified bessel function of the first kind Log(K1(x))                 
C                                                                       
C Taken from MAMBO by M Gibbs 03/11/92                                  
C                                                                       
      DOUBLE PRECISION BVPAKL,X,T,BETEMP,BVPAI1,BE                      
      IF(X.LE.0D0) THEN                                                 
        STOP 'BVPAKL: X<=0'                                             
      ELSEIF(X.LE.2D0) THEN                                             
        T=(X/2D0)**2                                                    
        BETEMP=X*DLOG(X/2D0)*BVPAI1(X) +                                
     .  1D0+T*(0.15443144+T*(-0.67278579+T*(-0.18156897+T*(             
     .  -0.01919402+T*(-0.00110404-T*0.00004686)))))                    
        BE=BETEMP/X                                                     
        BVPAKL=DLOG(BE)                                                 
      ELSE                                                              
        T=2D0/X                                                         
        BETEMP=1.25331414+T*(0.23498619+T*(-0.03655620+T*(              
     .  0.01504268+T*(-0.00780353+T*(0.00325614-T*0.00068245)))))       
        BVPAKL=(-X)-0.5D0*DLOG(X)+DLOG(BETEMP)                          
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPAKR(X,FK0,FK1,FK2,FK3)                              
C                                                                       
C K0(x)/K1(x) and the first three derivatives                           
C                                                                       
C Routines used:BVPAK0    Modified bessel function I0                   
C               BVPAK1    Modified bessel function I1                   
C               BVPABR    Modified bessel fn. asymptotic expansion      
C                                                                       
C Taken from MAMBO program by M Gibbs 03/11/92                          
C                                                                       
      DOUBLE PRECISION X,FK0,FK1,FK2,FK3,R1,R2,R3,R4,X2,                
     &         BVPAK0,BVPAK1                                            
      IF(X.LE.0D0) THEN                                                 
        WRITE(*,*) ' ROUTINE BVPAKR: X=',X,' NOT ALLOWED'               
        STOP                                                            
      ENDIF                                                             
      IF(X.LT.10D0) THEN                                                
        R1=BVPAK0(X)/BVPAK1(X)                                          
        R2=R1*R1                                                        
        R3=R2*R1                                                        
        R4=R3*R1                                                        
        X2=X*X                                                          
        FK0=R1                                                          
        FK1=( R2*X + R1 - X )/X                                         
        FK2=( 2.*R3*X + 3.*R2 - 2.*R1*X - 1.)/X                         
        FK3=( 6.*R4*X2 + 12.*R3*X + R2*(3.-8.*X2) -                     
     .        8*R1*X + 2.*X2 + 1. )/X2                                  
      ELSE                                                              
        CALL BVPABR(-X,FK0,FK1,FK2,FK3)                                 
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPAKQ(M1,M2,BETA)                                       
C                                                                       
C Function to calculate ratio of bessel functions given M1,M2,BETA      
C                                                                       
C Returns ps-weighted prob. for m1 particle as opposed to m2 (M8/23)    
C                                                                       
C Only tests for massless M1                                            
C                                                                       
C Written by M Gibbs 31/01/94                                           
C                                                                       
C 26/03/94 Moved hvdevl -> saddle                                       
C                                                                       
      IMPLICIT NONE                                                     
      DOUBLE PRECISION BVPAKQ,M1,M2,BETA,                               
     &            LMCOFF,BVPAK1,TEMP,T1,T2                              
      EXTERNAL BVPAK1                                                   
C Lower mass cut off                                                    
      LMCOFF = 1.0D-4                                                   
C Test for masslessness                                                 
      IF (M1.LT.LMCOFF) THEN                                            
C Massless case                                                         
        TEMP = M2*BETA                                                  
        T2 = TEMP * BVPAK1(TEMP)                                        
      ELSE                                                              
C Both massive                                                          
C        PRINT *,M1,M2,BETA                                             
        TEMP = M1*BETA                                                  
        T1 = M1 * BVPAK1(TEMP)                                          
        T2 = M2 * BETA                                                  
        TEMP = M2 * BVPAK1(T2)                                          
        T2 = TEMP/T1                                                    
      ENDIF                                                             
C      PRINT *,TEMP,T1,T2                                               
C Calc probability of M1                                                
      TEMP = 1.0D0 + T2                                                 
      BVPAKQ = 1.0D0 / TEMP                                             
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPAK0(X)                                                
C                                                                       
C Modified bessel function of the second kind K0(x)                     
C                                                                       
C Taken from MAMBO by M Gibbs 03/11/92                                  
C                                                                       
      DOUBLE PRECISION X,BVPAK0,T,BETEMP,BVPAI0                         
      IF(X.LT.0D0) THEN                                                 
        STOP 'BVPAK0: X<0'                                              
      ELSEIF(X.LE.2D0) THEN                                             
        T=(X/2D0)**2                                                    
        BVPAK0=-DLOG(X/2D0)*BVPAI0(X)                                   
     .  -0.57721566+T*(0.42278420+T*(0.23069756+T*(0.03488590+T*(       
     .  0.00262698+T*(0.00010750+T*0.00000740)))))                      
      ELSE                                                              
        T=2D0/X                                                         
        BETEMP=1.25331414+T*(-0.07832358+T*(+0.02189568+T*(             
     .  -0.01062446+T*(+0.00587872+T*(-0.00251540+T*0.00053208)))))     
        BVPAK0=DEXP(-X)/DSQRT(X)*BETEMP                                 
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPAK1(X)                                                
C                                                                       
C Modified bessel function of the second kind K1(x)                     
C                                                                       
C Taken from MAMBO by M Gibbs 03/11/92                                  
C                                                                       
      DOUBLE PRECISION X,BVPAK1,T,BETEMP,BVPAI1                         
      IF(X.LE.0D0) THEN                                                 
        STOP 'BVPAK1: X<=0'                                             
      ELSEIF(X.LE.2D0) THEN                                             
        T=(X/2D0)**2                                                    
        BETEMP=X*DLOG(X/2D0)*BVPAI1(X) +                                
     .  1D0+T*(0.15443144+T*(-0.67278579+T*(-0.18156897+T*(             
     .  -0.01919402+T*(-0.00110404-T*0.00004686)))))                    
        BVPAK1=BETEMP/X                                                 
      ELSE                                                              
        T=2D0/X                                                         
        BETEMP=1.25331414+T*(0.23498619+T*(-0.03655620+T*(              
     .  0.01504268+T*(-0.00780353+T*(0.00325614-T*0.00068245)))))       
        BVPAK1=DEXP(-X)/DSQRT(X)*BETEMP                                 
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPALF(ZINPUT)                                           
C                                                                       
C Calculation of log of gamma function                                  
C                                                                       
C 01/02/93 Taken from HWSGAM by M Gibbs                                 
C                                                                       
      DOUBLE PRECISION BVPALF                                           
      INTEGER I                                                         
      DOUBLE PRECISION ZINPUT,HLNTPI,Z,SHIFT,G,T,RECZSQ                 
C                                                                       
C   Gamma function computed by eq. 6.1.40, Abramowitz.                  
C   B(M) = B2m/(2m *(2m-1)) where B2m is the 2m'th Bernoulli number.    
C   HLNTPI = .5*LOG(2.*PI)                                              
C                                                                       
      DOUBLE PRECISION B(10)                                            
      DATA B/                                                           
     1     0.83333333333333333333D-01,   -0.27777777777777777778D-02, 
     1     0.79365079365079365079D-03,   -0.59523809523809523810D-03, 
     1     0.84175084175084175084D-03,   -0.19175269175269175269D-02, 
     1     0.64102564102564102564D-02,   -0.29550653594771241830D-01, 
     1     0.17964437236883057316D0  ,    -1.3924322169059011164D0  / 
      DATA HLNTPI/0.91893853320467274178D0/                             
C                                                                       
C   Shift argument to large value ( > 20 )                              
C                                                                       
      Z=ZINPUT                                                          
      SHIFT=1.                                                          
 10    IF (Z.LT.20.D0) THEN                                              
         SHIFT = SHIFT*Z                                                
         Z = Z + 1.D0                                                   
         GO TO 10                                                       
      ENDIF                                                             
C                                                                       
C   Compute asymptotic formula                                          
C                                                                       
      G = (Z-.5D0)*LOG(Z) - Z + HLNTPI                                  
      T = 1.D0/Z                                                        
      RECZSQ = T**2                                                     
      DO 20 I = 1,10                                                    
         G = G + B(I)*T                                                 
         T = T*RECZSQ                                                   
 20       CONTINUE                                                          
      BVPALF = G-LOG(SHIFT)                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPCAL(P,WTFACT)                                       
C                                                                       
C Subroutine to calculate phase space contribution for a given          
C   momentum distribution,using BVI routines.                           
C                                                                       
C Includes: Boson coefficient           BVPPID=1                        
C           Mod Pi for light particles  BVPPID=2                        
C           Boson Coeff. for gamma's    BVPPID=3                        
C           Boson Coeff. for Z0's       BVPPID=4                        
C           No coefficient for Higgs    BVPPID=5                        
C                                                                       
C P         Momentum of particles(in CofM frame)                        
C WTFACT    Calculated log of wt. factor for this distbn.               
C                                                                       
C Routines used:BVPBFA                                                  
C                                                                       
C Written by M Gibbs 26/3/92                                            
C                                                                       
C 19/10/92 Modified to use SADDLE10.INC include file                    
C 29/03/93 Fermion matrix element changed to E_i                        
C 08/04/93 Taken from saddle.for to allow inclusion of gamma's          
C          Currently has very botched approach to masses!               
C 12/04/93 Now takes parameters from common blocks                      
C 16/04/93 Returned to saddle.for                                       
C          Higgs code (PID=5) added                                     
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER I,J,K,PID                                                 
      DOUBLE PRECISION WTFACT,P(5,200),C1,C2,WMASS2                     
      C1 = LOG(CABB)                                                    
      C2 = LOG(1.D0 - CABB)                                             
      WMASS2 = WMASS**2                                                 
      WTFACT=0                                                          
      K=0                                                               
      DO 10 I=1,NTYPE                                                   
        IF (BVPNUM(I).GT.0) THEN                                        
          PID=BVPPID(I)                                                 
          IF (PID.EQ.1) THEN                                            
C---Boson factors                                                       
            DO 20 J=1,BVPNUM(I)                                         
            WTFACT=WTFACT+LOG(((4.D0*P(4,J+K)*P(4,J+K)                  
     &            /(P(5,J+K)*P(5,J+K)))-1.D0)*2./3.)                    
 20              CONTINUE                                                      
          ELSEIF (PID.EQ.2) THEN                                        
C---Light particle factors |Pi|                                         
            DO 30 J=1,BVPNUM(I)                                         
              WTFACT=WTFACT+(0.5*LOG((P(4,J+K)*P(4,J+K))                
C     &                    - (P(5,J+K)*P(5,J+K))                        
     &    ))                                                            
 30                  CONTINUE                                                    
          ELSEIF (PID.EQ.3) THEN                                        
            DO 40 J=1,BVPNUM(I)                                         
              WTFACT=WTFACT+LOG(((4.D0*P(4,J+K)*P(4,J+K)                
     &              /(WMASS2)))*2./3.)+C1                               
 40                  CONTINUE                                                    
          ELSEIF (PID.EQ.4) THEN                                        
            DO 50 J=1,BVPNUM(I)                                         
              WTFACT=WTFACT+LOG(((4.D0*P(4,J+K)*P(4,J+K)                
     &              /(P(5,J+K)**2))-1.D0)*2./3.)+C2                     
 50                  CONTINUE                                                    
C---Do nothing for the Higgs code                                       
          ELSEIF (PID.EQ.5) THEN                                        
            CONTINUE                                                    
          ENDIF                                                         
          K=K+BVPNUM(I)                                                 
        ENDIF                                                           
 10      CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPDIN                                                 
C                                                                       
C Routine to initialise common block with physical params               
C                                                                       
C Written by M Gibbs 12/04/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
C---Cabbibo angle for Z0/gamma probabilities                            
      CABB = 0.23                                                       
C---W boson mass/GeV                                                    
      WMASS = 81.                                                       
C---Mass cut-off for saddle-point calculations                          
      MCOFF = 1.D-5                                                     
C---Set BVPMML to not reject on MAMBO weights                           
      REJCON = .FALSE.                                                  
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPINI                                                 
C                                                                       
C Initialisation routine for saddle point integrals                     
C                                                                       
C 1. Initialize counter BVPTOT                                          
C 2. Zero result common blocks                                          
C                                                                       
C Written by M Gibbs 17/10/92                                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER I                                                         
      BVPTOT=0                                                          
      DO 10 I=1,NTYPE                                                   
         BVPTOT=BVPTOT+BVPNUM(I)                                        
 10       CONTINUE                                                          
      DO 20 I=1,5                                                       
         BVPDER(I)=0.D0                                                 
 20       CONTINUE                                                          
      BVPPSI=0.D0                                                       
      BVPPCO=0.D0                                                       
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPINS(NBOSON,NHIGGS)                                    
C                                                                       
C Instanton matrix element for NBOSON bosons                            
C                                                                       
C NBOSON    Integer number of bosons                                    
C BVPINS    Matrix element contribution                                 
C           (Ringwald Nuc.Phys.B330,(1990),1)                           
C                                                                       
C Written by M Gibbs 24/3/92                                            
C                                                                       
C Routines used:S14ABF (NAG) log n! approx.                             
C                                                                       
C 24/04/92 Modified to use NAG approx. for log n!                       
C 19/10/92 Higgs dependence included also                               
C 13/02/93 NAG dependence removed - now calls BVPALF                    
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER NBOSON,NHIGGS                                             
      DOUBLE PRECISION BVPINS,NB,BVPALF,NH                              
      NB=DBLE(NBOSON)                                                   
      NH=DBLE(NHIGGS)                                                   
C Log G from Ringwald paper.                                            
      BVPINS=((-232.091090763153D0)*2.0D0)                              
     &      +((NB+NH)*LOG(2.0D0))-(2.D0*(NB+NH)*LOG(246.0D0))           
      BVPINS=BVPINS+(NB*LOG(3.0D0))                                     
     &      +(2.0D0*(BVPALF(NB+NH+(103.D0/12.0D0))))                    
     &      -19.44608026D0                                              
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPINT(ENERGY,ACCUR,INT,COR)                           
C                                                                       
C Phase space integral for energy EN and NR accuracy ACCU               
C                                                                       
C Returned variables:                                                   
C INT   Log. of integrand, to base e                                    
C COR   Correction to integral                                          
C                                                                       
C Written by M Gibbs 17/10/92                                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION ENERGY,ACCUR,INT,COR                             
      CALL BVPINI                                                       
      EN=ENERGY                                                         
      ACCU=ACCUR                                                        
      CALL BVPSPS                                                       
      INT=BVPPSI                                                        
      COR=BVPPCO                                                        
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPMAM(N,RS,P,WXI)                                     
C                                                                       
C Routine to generate phase-space configurations for large numbers      
C of massive particles.                                                 
C                                                                       
C AUTHORS: R. KLEISS AND W.J.STIRLING                                   
C VERSION OF DECEMBER 5, 1991                                           
C                                                                       
C Modified 11/2/92 by M Gibbs to allow 5-component P vectors (5=Mass)   
C Now requires HWRGEN of HERWIG for random number generation            
C                                                                       
C N      = NUMBER OF PARTICLES (>2, <NMAX)                              
C RS     = TOTAL ENERGY                                                 
C P(5,N) = MASSES (MAY BE 0 BUT NOT NEGATIVE)                           
C P(4,N) = MOMENTA                                                      
C WXI    = Weight (relative) of config. returned                        
C                                                                       
C Common Blocks required for MAMBO :                                    
C      COMMON / TPARS / TTRY,TACC,XTRY,XACC                             
C      COMMON / BVIMAM / INIT,PRTEN                                     
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      LOGICAL PRTEN                                                     
                                                                        
* PARAMETERS                                                            
      PARAMETER(NMAX=200)                                               
      PARAMETER(TWOPI=6.28318530717960D0)                               
                                                                        
* VARIABLE-DIMENSION ARRAYS (ARGUMENTS)                                 
      DIMENSION P(5,NMAX)                                               
                                                                        
* FIXED-DIMENSION ARRAYS (INTERNAL VARIABLES)                           
      DIMENSION RM2(NMAX),AL(NMAX),UM(NMAX),VM(NMAX)                    
      DIMENSION QK(NMAX,4),Q(NMAX,4)                                    
      DIMENSION QKTOT(4)                                                
      DIMENSION EN(NMAX),Q2(NMAX)                                       
                                                                        
* AUXILIARY VARIABLES FOR EFFICIENCY ANALYSIS                           
      COMMON / TPARS / TTRY,TACC,XTRY,XACC                              
                                                                        
* INITIALIZE                                                            
      COMMON / BVIMAM / INIT,PRTEN                                      
                                                                        
      IF (INIT.EQ.0) THEN                                               
        IF (PRTEN) THEN                                                 
           WRITE(*,*) ' Starting MAMBO initialization'                  
        ENDIF                                                           
* CHECK ON ADMISSIBILITY OF NUMBER OF PARTICLES                         
        IF(N.LE.2) THEN                                                 
          WRITE(*,*) 'TOO FEW PARTICLES: N=',N                          
          STOP                                                          
        ENDIF                                                           
        IF(N.GT.NMAX) THEN                                              
          WRITE(*,*) 'TOO MANY PARTICLES: N=',N                         
          WRITE(*,*) 'INCREASE PARAMETER NMAX'                          
          STOP                                                          
        ENDIF                                                           
                                                                        
* COMPUTE TOTAL MASS AND CHECK ON ADMISSIBILITY                         
        RMTOT=0.                                                        
        RM2TOT=0.                                                       
        DO 1 K=1,N                                                      
          IF(P(5,K).LT.0D0) THEN                                        
            WRITE(*,*) 'MASS ',K,' = ',P(5,K),' IS NEGATIVE'            
            STOP                                                        
          ENDIF                                                         
          RMTOT=RMTOT+P(5,K)                                            
          RM2(K)=P(5,K)**2                                              
          RM2TOT=RM2TOT+RM2(K)                                          
 1           CONTINUE                                                        
        IF(RMTOT.GE.RS) THEN                                            
          WRITE(*,*) 'TOTAL MASS   =',RMTOT                             
          WRITE(*,*) 'TOTAL ENERGY =',RS                                
          WRITE(*,*) 'TOTAL MASS IS TOO HIGH'                           
          STOP                                                          
        ENDIF                                                           
                                                                        
        S=RS*RS                                                         
                                                                        
* THE MASS CONFIGURATION IS NOW ADMISSIBLE: REPRODUCE                   
        IF (PRTEN) THEN                                                 
          WRITE(*,*) 'NUMBER OF PARTICLES  =',N                         
          WRITE(*,*) 'TOTAL INVARIANT MASS =',RS                        
          WRITE(*,*) 'MASSES:'                                          
          WRITE(*,99) (P(5,K),K=1,N)                                    
 99          FORMAT(7F10.4)                                                  
* NOW WE DETERMINE THE OPTIMAL VALUE OF W                               
                                                                        
* THE LOWER AND UPPER BOUNDS FOR THE VALUE OF W                         
          WRITE(*,*) 'STARTING THE DETERMINATION OF W'                  
        ENDIF                                                           
        WB=DSQRT((N*S-RMTOT*RMTOT)/(N*N*(N-1.)))-RMTOT/N                
        WMIN=1./2.*WB                                                   
        WMAX=2./3.*WB                                                   
        IF (PRTEN) THEN                                                 
          WRITE(*,*) 'LOWER BOUND ON W =',WMIN                          
          WRITE(*,*) 'UPPER BOUND ON W =',WMAX                          
                                                                        
* ITERATE NEWTON-RAPHSON TO FIND W                                      
          WRITE(*,*) 'START NEWTON-RAPHSON ITERATIONS'                  
          WRITE(*,8) 'W_OLD','U0(W)-S','U1(W)','W','W/W_OLD-1'          
        ENDIF                                                           
 8         FORMAT(5A15)                                                    
        WOLD=WMAX                                                       
        ITERW=0                                                         
 101       CONTINUE                                                        
          ITERW=ITERW+1                                                 
          SF=0.                                                         
          SF1=0.                                                        
          SFF1=0.                                                       
          SM2F2=0.                                                      
          DO 102 K=1,N                                                  
            RN=P(5,K)/WOLD                                              
            IF(RN.EQ.0D0) THEN                                          
              F=2.*WOLD                                                 
              F1=2.                                                     
            ELSE                                                        
              CALL BVPAKR(RN,FK0,FK1,FK2,FK3)                           
              F=WOLD*(2.+RN*FK0)                                        
              F1=2.-RN*RN*FK1                                           
            ENDIF                                                       
            SF=SF+F                                                     
            SF1=SF1+F1                                                  
            SFF1=SFF1+F*F1                                              
            SM2F2=SM2F2-F*F                                             
 102             CONTINUE                                                      
          U0=SF*SF+SM2F2+RM2TOT                                         
          U1=2.*(SF*SF1-SFF1)                                           
          W=WOLD-(U0-S)/U1                                              
          IF (PRTEN) THEN                                               
            WRITE(*,9) WOLD,U0-S,U1,W,W/WOLD-1.                         
          ENDIF                                                         
 9             FORMAT(5D15.8)                                                
          IF(DABS(W/WOLD-1.).GT.1D-12) THEN                             
            WOLD=W                                                      
            GOTO 101                                                    
          ENDIF                                                         
        CONTINUE                                                        
                                                                        
* W IS NOW FOUND                                                        
        IF (PRTEN) THEN                                                 
          WRITE(*,*) 'VALUE OF W OBTAINED: W=',W,                       
     .             ' AFTER ',ITERW,' ITERATIONS'                        
          WRITE(*,*) 'NORMALIZED VALUE: W/WB=',W/WB                     
        ENDIF                                                           
* COMPUTE AUXILIARY VARIABLES FOR KINDERMAN-MONAHAN ALGORITHM           
* AND COMPUTE EXPECTED EFFICIENCY FOR THIS STEP                         
        SINVEF=0.                                                       
        DO 2 K=1,N                                                      
          B=2.*P(5,K)/W                                                 
          XU=(1.-B+DSQRT(1.+B*B))/2.                                    
          XV=(3.-B+DSQRT(9.+B*(4.+B)))/2.                               
          UM(K)=DEXP(DLOG(XU*(XU+B))/4.-XU/2.)                          
          VM(K)=XV*DEXP(DLOG(XV*(XV+B))/4.-XV/2.)                       
          AL(K)=B                                                       
          IF(P(5,K).EQ.0D0) THEN                                        
            EFNUM=1.                                                    
          ELSE                                                          
            EFNUM=P(5,K)/W*DEXP(P(5,K)/W+BVPAKL(P(5,K)/W))              
          ENDIF                                                         
          SINVEF=SINVEF+2.*UM(K)*VM(K)/EFNUM                            
 2           CONTINUE                                                        
        IF (PRTEN) THEN                                                 
          WRITE(*,*) 'RATIO-OF-UNIFORMS PARAMETERS COMPUTED'            
          WRITE(*,*) 'ESTIMATED EFFICENCY OF ',                         
     .    'ENERGY GENERATION IS',100.*N/SINVEF,' %'                     
        ENDIF                                                           
* END OF INITIALIZATION STEP                                            
                                                                        
        TTRY=0.                                                         
        TACC=0.                                                         
      ENDIF                                                             
                                                                        
* START GENERATION STEP: GENERATE K MOMENTA, COMPUTE X AND REJECT/ACCEPT
                                                                        
* STARTING POINT OF W(XI) REJECTION/ACCEPTANCE LOOP                     
 3000  CONTINUE                                                          
                                                                        
* STARTING POINT OF X REJECTION/ACCEPTANCE LOOP                         
 3001   CONTINUE                                                          
      XTRY=XTRY+1.                                                      
                                                                        
* LOOP OVER THE K MOMENTA                                               
      DO 11 K=1,N                                                       
                                                                        
* STARTING POINT FOR KI(0) REJECTION/ACCEPTANCE STEP                    
 3002       CONTINUE                                                        
          TTRY=TTRY+1.                                                  
C Fix to stop divide by zero                                            
 5801          U=HWRGEN(1)*UM(K)                                             
          IF (U.GT.0.) GOTO 5802                                        
          GOTO 5801                                                     
 5802          V=HWRGEN(2)*VM(K)                                             
          X=V/U                                                         
           IF(X.GT.175D0) GOTO 3002                                     
           IF(U.GT.(X+AL(K)*0.5D0)) GOTO 3002                           
          IF((U*U).GT.DSQRT(X*(X+AL(K)))*DEXP(-X)) GOTO 3002            
        CONTINUE                                                        
        QK(K,4)=P(5,K)+W*X                                              
        TACC=TACC+1.                                                    
                                                                        
* THE ENERGY KI(0) IS NOW ACCEPTED: GENERATE SPACELIKE PART             
        CK=2.*HWRGEN(3)-1.                                              
        FK=TWOPI*HWRGEN(4)                                              
        SK=DSQRT(1.-CK*CK)                                              
        QKV=W*DSQRT(X*(X+AL(K)))                                        
        QK(K,1)=QKV*SK*DSIN(FK)                                         
        QK(K,2)=QKV*SK*DCOS(FK)                                         
        QK(K,3)=QKV*CK                                                  
 11      CONTINUE                                                          
                                                                        
* THE K MOMENTA HAVE NOW BEEN GENERATED                                 
                                                                        
* COMPUTE TOTAL INVARIANT MASS AND X                                    
      DO 22 I=1,4                                                       
        QKTOT(I)=0.                                                     
        DO 21 K=1,N                                                     
          QKTOT(I)=QKTOT(I)+QK(K,I)                                     
 21          CONTINUE                                                        
 22           CONTINUE                                                          
      QKTSQ=QKTOT(4)**2-QKTOT(1)**2-QKTOT(2)**2-QKTOT(3)**2             
      X2=S/QKTSQ                                                        
      X=DSQRT(X2)                                                       
                                                                        
* REJECT ON TOO HIGH X VALUES                                           
      IF(X2.GT.1D0) GOTO 3001                                           
                                                                        
* X HAS NOW BEEN ACCEPTED                                               
      XACC=XACC+1.                                                      
      WK=RS/X                                                           
                                                                        
* BOOST AND SCALE THE K MOMENTA INTO THE Q MOMENTA                      
      DO 33 K=1,N                                                       
        Q(K,4)=(QK(K,4)*QKTOT(4)-QK(K,1)*QKTOT(1)                       
     .         -QK(K,2)*QKTOT(2)-QK(K,3)*QKTOT(3))/WK                   
        T=(Q(K,4)+QK(K,4))/(QKTOT(4)+WK)                                
        DO 31 I=1,3                                                     
          Q(K,I)=QK(K,I)-T*QKTOT(I)                                     
 31          CONTINUE                                                        
        DO 32 I=1,4                                                     
          Q(K,I)=X*Q(K,I)                                               
 32          CONTINUE                                                        
 33           CONTINUE                                                          
                                                                        
* NOW THE STEP TO FIND XI                                               
                                                                        
* SQUARED LENGTHS OF THE Q THREE-MOMENTA                                
      DO 41 K=1,N                                                       
        Q2(K)=Q(K,4)*Q(K,4)-X2*RM2(K)                                   
 41      CONTINUE                                                          
                                                                        
* START NEWTON-RAPHSON ITERATIONS FOR XI                                
      XIOLD=1.                                                          
 50    CONTINUE                                                          
        F0=-RS                                                          
        F1=0.                                                           
        XI2=XIOLD*XIOLD                                                 
        DO 51 K=1,N                                                     
          EN(K)=DSQRT(XI2*Q2(K)+RM2(K))                                 
          F0=F0+EN(K)                                                   
          F1=F1+Q2(K)/EN(K)                                             
 51          CONTINUE                                                        
        XI=XIOLD-F0/(XIOLD*F1)                                          
        IF(DABS(XI/XIOLD-1.).GT.1D-12) THEN                             
          XIOLD=XI                                                      
          GOTO 50                                                       
        ENDIF                                                           
      CONTINUE                                                          
                                                                        
* XI FOUND: THE EN() ARE ALSO THE P ENERGIES                            
      DO 62 K=1,N                                                       
        P(4,K)=EN(K)                                                    
        DO 61 I=1,3                                                     
          P(I,K)=XI*Q(K,I)                                              
 61          CONTINUE                                                        
 62           CONTINUE                                                          
                                                                        
* COMPUTE THE W(XI) WEIGHT                                              
      T1=1.                                                             
      T2=0.                                                             
      T3=0.                                                             
      DO 71 K=1,N                                                       
        T1=T1*Q(K,4)/P(4,K)                                             
        T2=T2+RM2(K)/Q(K,4)                                             
        T3=T3+RM2(K)/P(4,K)                                             
 71      CONTINUE                                                          
      WXI=XI**(3.*N-3.)*T1*(RS-X2*T2)/(RS-T3)                           
      END                                                               
C                                                                       
C********************************************************************   
C                                                                       
      FUNCTION BVPMBF(START,ISTEP,BOUND)                                
C                                                                       
C Routine to search in boson number given limit rel. to max of BOUND    
C Assumes value for max has already been estimated                      
C                                                                       
C Routines used:BVPMFC                                                  
C               BVPMFM                                                  
C               BVPNBS                                                  
C                                                                       
C Written by M Gibbs 13/02/93                                           
C                                                                       
C 08/04/93 Modified to use BVPNBS Boson no. selection routine           
C                                                                       
      INCLUDE 'HERWIG65.INC'                                         
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION AN,BOUND,LNREF,MASST
      INTEGER TP,START,STEP,LNN,BVPMBF,NP,ISTEP,SUCC
      LNREF=RES(NBMAX)                                                  
      LNN=NBMAX                                                         
      TP=START                                                          
      STEP=ISTEP                                                        
 10     TP=TP+STEP                                                        
      IF (TP.LE.20) GOTO 11                                             
      CALL BVPNBS(TP)                                                   
      IF (RES(TP).LT.(-5.0D2)) THEN                                     
         CALL BVPMFC(1,NP,MASST)                                        
         CALL BVPMFM(PHEP,EN,NITS,AN,SUCC)                              
         RES(TP)=AN                                                     
      ELSE                                                              
         AN=RES(TP)                                                     
      ENDIF                                                             
      IF ((LNREF-AN).LE.(0.0D0)) THEN                                   
         LNREF=AN                                                       
         LNN=TP                                                         
      ENDIF                                                             
      RES(TP)=AN                                                        
      STEP=(1+INT(BOUND-LNREF+AN))*(ISTEP/ABS(ISTEP))                   
      IF ((LNREF-AN).LE.BOUND) GOTO 10                                  
 11    BVPMBF=TP                                                         
      NBMAX=LNN                                                         
      END                                                               
C                                                                       
C***********************************************************************
C                                                                       
      SUBROUTINE BVPMBG                                                 
C                                                                       
C Routine to fill gaps between BMIN and BMAX for boson calc             
C                                                                       
C Routines used:BVPMFC                                                  
C               BVPMFM                                                  
C               BVPNBS                                                  
C                                                                       
C Written by M Gibbs 13/02/93                                           
C                                                                       
C 08/04/93 Modified to use BVPNBS boson number sel. routine             
C                                                                       
      INCLUDE 'HERWIG65.INC'                                           
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION MASST
      INTEGER NP,SUCC,I                                          
      DO 10 I=BMIN+1,BMAX-1                                             
        IF (RES(I).LT.(-5.0D2)) THEN                                    
          CALL BVPNBS(I)                                                
          CALL BVPMFC(1,NP,MASST)                                       
          CALL BVPMFM(PHEP,EN,NITS,RES(I),SUCC)                         
          IF (RES(I).GE.RES(NBMAX)) NBMAX = I                           
        ENDIF                                                           
 10      CONTINUE                                                          
      END                                                               
C                                                                       
C***********************************************************************
C                                                                       
      SUBROUTINE BVPMBI                                                 
C                                                                       
C Routine to initialize storage arrays for Boson distribution calc.     
C                                                                       
C Written by M Gibbs 14/02/93                                           
C                                                                       
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER I                                                         
      NBMAX = 0                                                         
      DO 10 I=1,200                                                     
         RES(I)=-1.0D3                                                  
 10       CONTINUE                                                          
      END                                                               
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE BVPMBP                                                 
C                                                                       
C Routine to find peak of distribution                                  
C                                                                       
C Written by M Gibbs 14/02/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER I                                                         
      DOUBLE PRECISION TEMP                                             
      TEMP=RES(NBMAX)                                                   
      DO 10 I=BMIN+1,BMAX                                               
        IF(RES(I).GT.TEMP) THEN                                         
          NBMAX=I                                                       
          TEMP=RES(I)                                                   
        ENDIF                                                           
 10      CONTINUE                                                          
      END                                                               
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE BVPMBS(BOUND)                                          
C                                                                       
C Routine to calc. boson contribs within bound BOUND of maximum         
C If no maximum value available (NBMAX=0) then guesses                  
C a value of EN/(2*MASSZ)                                               
C                                                                       
C Routines used:BVPMBF                                                  
C               BVPMFM                                                  
C               BVPMFC                                                  
C               BVPNBS                                                  
C                                                                       
C Written by M Gibbs 13/02/93                                           
C                                                                       
C 05/03/93 Modified to set up CONFIG array for boson numbers            
C          HERWIG56.INC file now required                               
C 08/04/93 Modified to use BVPNBS boson number selection routine        
C                                                                       
      INCLUDE 'HERWIG65.INC'                                        
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION BOUND,MASST
      INTEGER GU,NP,BVPMBF,SUCC 
      IF (NBMAX.EQ.0) THEN                                              
C Estimate using M_Z0 from HW56 include file                            
        GU = EN/(2*RMASS(200))                                          
        CALL BVPNBS(GU)                                                 
        CALL BVPMFC(1,NP,MASST)                                         
        CALL BVPMFM(PHEP,EN,NITS,RES(GU),SUCC)                          
        NBMAX=GU                                                        
      ENDIF                                                             
      BMIN = BVPMBF(NBMAX,-5,BOUND)                                     
      BMAX = BVPMBF(NBMAX,5,BOUND)                                      
      CALL BVPMBG                                                       
      END                                                               
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE BVPMBT                                                 
C                                                                       
C Routine to sum up boson contribs at a partic. en.                     
C The Log. of the sum of all contributions and                          
C <n_{bos}**2>are returned via common block                             
C                                                                       
C Written by M Gibbs 14/02/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION REF,LNREF
      INTEGER I                                                  
      NBOSSQ = 0.D0                                                     
      REF=RES(NBMAX)                                                    
      TSUM=0.D0                                                         
      DO 10 I=BMIN,BMAX                                                 
         TSUM = TSUM+DEXP(RES(I)-REF)                                   
         NBOSSQ = NBOSSQ + (DBLE(I*I)*DEXP(RES(I)-REF))                 
 10       CONTINUE                                                          
C      PRINT *,' Total sum without ref value is ',TSUM                  
      NBOSSQ = NBOSSQ/TSUM                                              
      LNREF=LOG(TSUM)                                                   
      TSUM=LNREF+REF                                                    
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPMCI(P)                                              
C                                                                       
C Subroutine to initialize monte-carlo simulation and P-vector          
C                                                                       
C P(5,200)  En-Momentum vector                                          
C                                                                       
C Routines used:MAMBO    Monte-carlo event generator                    
C                                                                       
C Written by M Gibbs 18/06/92                                           
C                                                                       
C 19/10/92 SADDLE10.INC include file inserted                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION P(5,200),WXI                                     
      INTEGER INIT,I,J,K                                                
      LOGICAL PRTEN                                                     
      COMMON /BVIMAM/ INIT,PRTEN                                        
      PRTEN=.FALSE.                                                     
      K=0                                                               
      DO 10 I=1,NTYPE                                                   
        DO 20 J=1,BVPNUM(I)                                             
          P(5,K+J)=BVPMAS(I)                                            
 20          CONTINUE                                                        
        K=K+BVPNUM(I)                                                   
 10      CONTINUE                                                          
      INIT=0                                                            
      CALL BVPMAM(BVPTOT,EN,P,WXI)                                      
      INIT=1                                                            
      END                                                               
C                                                                       
C********************************************************************   
C                                                                       
      SUBROUTINE BVPMFC(START,NPS,MASST)                                
C                                                                       
C Routine to set up configuration for fermion part of BVI element       
C The data is stored at position START in HEPEVT common block           
C of HERWIG5.6                                                          
C                                                                       
C CONFIG contains numbers of each particle type                         
C        order of storage d,u,s,c,b,t,e+,neb,mu+,nmub,tau+,nmut         
C        (all quarks are qbar in above statement)                       
C                                                                       
C NP (returned) is the number of particles placed into the event.       
C                                                                       
C Written by M Gibbs 14/02/93                                           
C                                                                       
C 04/03/93 W+ W- Z0 codes added as 13,14,15th elements of CONFIG        
C 08/04/93 Gamma code added as 16th element of config                   
C          Updated to allow new BVPPID codes for Z0/Gamma               
C 09/04/93 Bug fix added - remainder of arrays zeroed at end of routine 
C          Routine transferred into saddev.for for Gamma updates        
C 16/04/93 Routine replaced into saddle.for                             
C          Higgs coefficient added                                      
C                                                                       
      INCLUDE 'HERWIG65.INC'                                         
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION MASST                                            
      INTEGER I,J,K,ID,START,NPS,POINT,CODES(17),PIDS(17)               
      INTEGER ZGCODE(2),ZIS,ZS                                          
      COMMON /BVPCQW/ ZS,ZIS,ZGCODE                                     
      DATA CODES /7,8,9,10,11,12,127,128,129,                           
     &             130,131,132,198,199,200,59,201/                      
      DATA PIDS /2,2,2,2,2,2,2,2,2,2,2,2,1,1,4,3,5/                     
      ZGCODE(1) = CODES(15)                                             
      ZGCODE(2) = CODES(16)                                             
      POINT = START                                                     
      NPS = 0                                                           
      K=1                                                               
      MASST = 0.0D0                                                     
      DO 10 I=1,NCONF                                                   
        IF ((CONFIG(I).GT.0) .OR. (I.EQ.16)) THEN                       
          BVPPID(K)=PIDS(I)                                             
C---Trap the Z0 case                                                    
          IF (I.EQ.15) THEN                                             
C---Now set variables for BVPNBG                                        
            ZIS = K                                                     
            ZS = POINT                                                  
          ENDIF                                                         
C---Set the internal and HEPEVT variables                               
          BVPNUM(K)=CONFIG(I)                                           
          BVPMAS(K)=RMASS(CODES(I))                                     
C---Enforce the mass cut-off on massless particles                      
          IF (BVPMAS(K).LT.MCOFF) BVPMAS(K)=MCOFF                       
C---Enter data into HW common block                                     
          IF (CONFIG(I).GT.0) THEN                                      
            MASST=MASST+(BVPMAS(K)*CONFIG(I))                           
            NPS=NPS+CONFIG(I)                                           
            DO 20 J=1,CONFIG(I)                                         
              ID=CODES(I)                                               
              IDHW(POINT)=ID                                            
              IDHEP(POINT)=IDPDG(ID)                                    
              PHEP(5,POINT)=BVPMAS(K)                                   
              POINT=POINT+1                                             
 20                  CONTINUE                                                    
          ENDIF                                                         
          K=K+1                                                         
        ENDIF                                                           
 10      CONTINUE                                                          
C---Now zero remainder of array                                         
      IF (K .LE. NTYPE) THEN                                            
        DO 30 I=K,NTYPE                                                 
          BVPNUM(I) = 0                                                 
          BVPPID(I) = 0                                                 
 30          CONTINUE                                                        
      ENDIF                                                             
C This next line must go!                                               
      NHEP = NPS                                                        
      BVPTOT = NPS                                                      
      END                                                               
C                                                                       
C********************************************************************   
C                                                                       
      SUBROUTINE BVPMFI(FLAG)                                           
C                                                                       
C Routine to set up CONFIG array for use in BVPMFC                      
C                                                                       
C FLAG(4) Logical variable containing structure description             
C                                                                       
C FLAG(1) .true.  first light quark dbar                                
C         .false. first light quark ubar                                
C FLAG(2) .true.  second light quark dbar                               
C         .false. second light quark ubar                               
C FLAG(3) .true.  cbar sbar sbar nmubar                                 
C         .false. cbar cbar sbar muplus                                 
C FLAG(4) .true.  tbar bbar bbar ntaubar                                
C         .false. tbar tbar tbar tauplus                                
C                                                                       
C FLAG(1) and FLAG(2) combinations:                                     
C                                                                       
C FLAG(1) FLAG(2) Result                                                
C   T        F       ubar eplus                                         
C   F        T       dbar nuebar                                        
C                                                                       
C Written by M Gibbs 16/02/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      LOGICAL FLAG(4)                                                   
      IF (FLAG(4)) THEN                                                 
        CONFIG(5)=2                                                     
        CONFIG(6)=1                                                     
        CONFIG(11)=0                                                    
        CONFIG(12)=1                                                    
      ELSE                                                              
        CONFIG(5)=1                                                     
        CONFIG(6)=2                                                     
        CONFIG(11)=1                                                    
        CONFIG(12)=0                                                    
      ENDIF                                                             
      IF (FLAG(3)) THEN                                                 
        CONFIG(3)=2                                                     
        CONFIG(4)=1                                                     
        CONFIG(9)=0                                                     
        CONFIG(10)=1                                                    
      ELSE                                                              
        CONFIG(3)=1                                                     
        CONFIG(4)=2                                                     
        CONFIG(9)=1                                                     
        CONFIG(10)=0                                                    
      ENDIF                                                             
      IF (FLAG(1)) THEN                                                 
        CONFIG(1)=0                                                     
        CONFIG(2)=1                                                     
      ELSE                                                              
        CONFIG(1)=1                                                     
        CONFIG(2)=0                                                     
      ENDIF                                                             
      IF (FLAG(2)) THEN                                                 
        CONFIG(7)=0                                                     
        CONFIG(8)=1                                                     
      ELSE                                                              
        CONFIG(7)=1                                                     
        CONFIG(8)=0                                                     
      ENDIF                                                             
      END                                                               
C                                                                       
C********************************************************************   
C                                                                       
      SUBROUTINE BVPMFM(P,REN,ITS,RESUL,SUCC)                           
C                                                                       
C Routine to MC configuration up                                        
C Uses BVPINT to calculate PS contribution                              
C                                                                       
C Written by M Gibbs 14/02/93                                           
C                                                                       
C 04/03/93 Modified to include bosons in calculation                    
C 08/04/93 NB changed to add in Gamma contribution also                 
C 12/04/93 In saddev.for for gamma selection improvement                
C 16/04/93 Replaced into saddle.for                                     
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER I,SUCC,ITS,NB,INIT,NWUN,TNITS,TS,SNIT                   
     &    ,ITT,SUT,GMAXP,NPARS                                          
      LOGICAL PRTEN                                                     
      COMMON /BVIMAM/ INIT,PRTEN                                        
      DOUBLE PRECISION P(5,200),REN,RESUL,D,COR,RACCU,PSWT,             
     &  CGEVPB,STOR(100),RWGT,BVPNBA,CM,REF,MASST,GMAXSF,               
     &  BVPINS,BVPALF                                                   
      COMMON /PERFZ/ ITT,SUT                                            
      NB = CONFIG(13)+CONFIG(14)+CONFIG(15)+CONFIG(16)                  
      NWUN = CONFIG(15) + CONFIG(16)                                    
      EN=REN                                                            
C      PRINT *,'PMFM entry, NITS',NITS                                  
C      PRINT *,'PMFM entry,  ITS',ITS                                   
      TNITS = 0                                                         
C      PRINT *,' DBLE(ITS)*0.5 is ',DBLE(ITS)*0.5D0                     
C      PRINT *,' INT of this is   ',INT(DBLE(ITS)*0.5D0)                
      SUCC = 0                                                          
      SNIT = ITS                                                        
C---Evaluate the instanton contribution part                            
      CM = -((3.D0*DBLE(BVPTOT)                                         
     &     -4.D0)*1.8378770664093D0)                                    
     &     +BVPINS(NB,CONFIG(17))+CGEVPB(0)                             
     &     -BVPALF(DBLE(NB+1))                                          
C---Loop over all possible Z/gamma splits                               
      GMAXP = 0                                                         
      GMAXSF = -1.D10                                                   
      DO 10 I=0,NWUN                                                    
        STOR(I+1) = -1.D10                                              
C---Set up the configuration                                            
        CONFIG(16) = I                                                  
        CONFIG(15) = NWUN - I                                           
        CALL BVPMFC(1,NPARS,MASST)                                      
C---Number of iterations:                                               
        RWGT = EXP(BVPNBA(NWUN,I,CABB))                                 
        NITS = INT(RWGT*DBLE(SNIT))                                     
C          PRINT *,'For ',I,' gammas from ',NWUN,' bosons.'             
C          PRINT *,'Total number of its is ',NITS                       
C          PRINT *,'Relative weight was ',RWGT                          
C          PRINT *,'Overall no. of its  ',ITS                           
        IF (NITS.GT.0) THEN                                             
C---Do the MC calculation for this configuration                        
          COR=0.D0                                                      
          RACCU=1.D-7                                                   
          CALL BVPINI                                                   
          INIT = 0                                                      
          PRTEN = .FALSE.                                               
          CALL BVPMAM(BVPTOT,EN,P(1,1),D)                               
          INIT = 1                                                      
          CALL BVPMML(NITS,P,RESUL,TS)                                  
          SUCC = SUCC + TS                                              
          ITT = ITT+NITS                                                
          SUT = SUT+TS                                                  
          CALL BVPINT(EN,RACCU,PSWT,COR)                                
          RESUL=RESUL+PSWT                                              
C          PRINT *,'Result for ',I,' gammas is ',RESUL                  
          STOR(I+1) = RESUL                                             
          IF (RESUL.GT.GMAXSF) THEN                                     
            GMAXSF = RESUL                                              
            GMAXP = I+1                                                 
C            PRINT *,' New peak - for ',I,' gammas'                     
          ENDIF                                                         
        ENDIF                                                           
 10      CONTINUE                                                          
C---Now sum result and add the instanton part                           
      REF = STOR(GMAXP)                                                 
      RESUL = 0.0D0                                                     
C      PRINT *,'Reference value is',REF                                 
      DO 20 I=1,NWUN                                                    
C        PRINT *,I,STOR(I),STOR(I)-REF                                  
        IF (STOR(I).GT.(-5.D9))                                         
     &    RESUL = RESUL + EXP(STOR(I)-REF)                              
C        PRINT *,RESUL                                                  
 20      CONTINUE                                                          
      RESUL = LOG(RESUL) + REF + CM                                     
C      PRINT *,'Result for ',NB,' bosons.'                              
C      PRINT *,RESUL                                                    
      ITS = SNIT                                                        
      NITS = SNIT                                                       
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPMML(ITER,P,WT,SUCC)                                 
C                                                                       
C Subroutine to perform main loop of monte-carlo simulation             
C Assumes MAMBO already initialized                                     
C                                                                       
C P(5,200)  En-mom vector                                               
C WT        Average wt. returned                                        
C ITER      Required number of iterations                               
C SUCC      Number of successful tries                                  
C                                                                       
C Routines used:MAMBO   Phase space generator                           
C               BVPCAL  Wt. calculating routine                         
C               HWRGEN  Random number generator                         
C                                                                       
C Written by M Gibbs 18/06/92                                           
C                                                                       
C 19/10/92 SADDLE10.INC file added                                      
C 29/03/93 MC factor of (sucesses/tries) added                          
C 12/04/93 MAMBO rejection step removed - now controlled by REJCON      
C 14/04/93 Weight calculation corrected                                 
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION P(5,200),WT,WXI,HWRGEN,WREF                                   
      INTEGER ITER,I,SUCC                                               
      SUCC=0                                                            
      WT=0.D0                                                           
      WREF=0.D0                                                         
      DO 10 I=1,ITER                                                    
 20          WXI=0.D0                                                        
        CALL BVPMAM(BVPTOT,EN,P,WXI)                                    
C---Reject event according to WXI if REJCON set                         
        IF (WXI.LT.HWRGEN(1)) THEN                                   
C---If Rejcon, then loop until we get an accepted event                 
          IF (REJCON) GO TO 20                                          
          CONTINUE                                                      
        ELSE                                                            
          SUCC=SUCC+1                                                   
          CALL BVPCAL(P,WXI)                                            
          IF (SUCC.EQ.1) THEN                                           
            WREF=WXI                                                    
            WXI=0.D0                                                    
          ELSE                                                          
            WXI=WXI-WREF                                                
          ENDIF                                                         
            WXI=EXP(WXI)                                                
            WT=WT+WXI                                                   
        ENDIF                                                           
 10      CONTINUE                                                          
C---Check to see that we have a result                                  
      IF (SUCC.GT.0) THEN                                               
C---Divide by the number of iterations to get result                    
C   (/succ as we have unweighted distrib. wrt the MAMBO part            
        WXI = WT/DBLE(SUCC)                                             
        WT = WREF + LOG(WXI)                                            
C---Otherwise return a v.low result                                     
      ELSE                                                              
        WT = -1.D10                                                     
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION BVPNBA(NZERO,NG,PROB)                                    
C                                                                       
C Function to evaluate binomial expression                              
C                                                                       
C NG     Number of gammas                                               
C NZERO  Total number of neutral particles                              
C PROB   Rel. prob. of a gammma (=sin^2 \theta_W)                       
C                                                                       
C Written by M Gibbs 09/04/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      DOUBLE PRECISION A,B,C,PROB,BVPALF,BVPNBA                         
      INTEGER NG,NZERO                                                  
      EXTERNAL BVPALF                                                   
      A = DBLE(NZERO)                                                   
      C = DBLE(NG)                                                      
      B = A-C                                                           
      BVPNBA = BVPALF(A+1.) - BVPALF(B+1.) - BVPALF(C+1.)               
     &      + C*LOG(PROB) + B*LOG(1.-PROB)                              
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPNBS(NBOSON)                                         
C                                                                       
C Routine to select number of bosons for BVI process                    
C                                                                       
C NBOSON  Total number of available bosons                              
C                                                                       
C Written by M Gibbs 08/04/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER NBOSON,NW,NZ                         
      DOUBLE PRECISION A1,A2
C---Original selection code                                             
      NZ = NBOSON / 3                                                   
      NW = NBOSON - NZ                                                  
C---Make sure equal numbers of W+/W- bosons                             
      CONFIG(13) = INT(NW/2)                                            
      CONFIG(14) = CONFIG(13)                                           
      NZ = NBOSON - (CONFIG(13)*2)                                      
C---Now select how many gammas to have                                  
      A1 = DBLE(NZ)*CABB                                                
      A2 = A1 -1.D0                                                     
      CALL BVPNGA(A1,A2,CONFIG(16))                                     
C---Make sure we don't have too many gammas                             
      IF (CONFIG(16).GT.NZ) CONFIG(16) = NZ                             
      CONFIG(15) = NZ - CONFIG(16)                                      
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPNGA(PEAK,SIGMA,ARES)                                
C                                                                       
C Routine to output gaussian approx. to number of gammas                
C                                                                       
C Written by M Gibbs 08/04/93                                           
C                                                                       
      INCLUDE 'HERWIG65.INC'                                        
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER ARES                                                      
      DOUBLE PRECISION PEAK,Z,HWRGAU,SIGMA                        
      Z = HWRGAU(0,PEAK,SIGMA)                                          
      ARES = INT(Z)                                                     
      IF (ARES.LT.0) ARES = 0                                           
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPSBR                                                 
C                                                                       
C Subroutine to evaluate modified bessel function ratios                
C                                                                       
C Derivatives are returned as derivatives of x,i.e. BETA*(MASS_var)     
C Form of derivatives taken from Kleiss and Stirling,                   
C CERN-TH.6293 Oct 1991                                                 
C                                                                       
C Routines used:S18CEF I0 Nag function                                  
C               S18CFF I1 Nag function                                  
C               S18CCF K0 Nag routine                                   
C               S18CDF K1 Nag routine                                   
C                                                                       
C Written by M Gibbs 13/4/92                                            
C                                                                       
C 15/06/92 Modified from BVPSGB to give ratios in common blocks         
C 17/10/92 Modified to use SADDLE10.INC include file                    
C 03/11/92 Modified to use SADDLE10.FOR Modified bessel fn. rtns        
C                                                                       
C ** Currently has no error trapping                                    
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION TEMP
      INTEGER N,I,IFAIL                                                 
      IFAIL=0                                                           
C---Modified bessel functions I0/I1 at energy EN                        
      TEMP=BETA*EN                                                      
      CALL BVPAIR(TEMP,BVPIR(1),BVPIR(2),BVPIR(3),BVPIR(4))             
C      IF ((TEMP).GT.(10.D0)) THEN                                      
C          CALL BVPABR(TEMP,BVPIR(1),BVPIR(2),BVPIR(3),BVPIR(4))        
C      ELSE                                                             
C             D0=S18CEF(TEMP,IFAIL)                                     
C             D1=S18CFF(TEMP,IFAIL)                                     
C             BVPIR(1)=D0/D1                                            
C             D3=BVPIR(1)*BVPIR(1)                                      
C             BVPIR(2)=1.D0+((BVPIR(1)/TEMP)-D3)                        
C             BVPIR(3)=((1.D0-(3.D0*D3))/TEMP)                          
C     &            +(2.D0*BVPIR(1)*D3)                                  
C             D2=D3*BVPIR(1)                                            
C             BVPIR(4)=(-2.D0)+((-1.D0/(TEMP*TEMP))*(1.D0+              
C     &           (8.D0*BVPIR(1)*TEMP)-(12.D0*TEMP*D2)+                 
C     &           ((3.D0-(8.D0*TEMP*TEMP))*D3)+                         
C     &           (6.D0*TEMP*TEMP*D3*D3)))                              
C      ENDIF                                                            
C---Loop over particle types for K0/K1 ratios                           
      DO 100 N=1,NTYPE                                                  
         IF (BVPNUM(N).GT.0) THEN                                       
             TEMP=BETA*BVPMAS(N)                                        
             CALL BVPAKR(TEMP,BVPKR(N,1),BVPKR(N,2),                    
     &        BVPKR(N,3),BVPKR(N,4))                                    
C             IF ((TEMP).GT.(10D0)) THEN                                
C                CALL BVPABR((-1.D0*TEMP),BVPKR(N,1),                   
C     &                 BVPKR(N,2),BVPKR(N,3),BVPKR(N,4))               
C             ELSE                                                      
C                D0=S18CCF(TEMP,IFAIL)                                  
C                D1=S18CDF(TEMP,IFAIL)                                  
C                BVPKR(N,1)=D0/D1                                       
C                D2=BVPKR(N,1)*BVPKR(N,1)                               
C                D3=D2*BVPKR(N,1)                                       
C                BVPKR(N,2)=(D2-1.D0)+(BVPKR(N,1)/TEMP)                 
C                BVPKR(N,3)=(2.D0*BVPKR(N,1)*(D2-1.D0))                 
C     &               +(3.D0*D2/TEMP)-1.D0                              
C                BVPKR(N,4)=2.D0+((1.D0/(TEMP*TEMP))*(1.D0              
C     &                +((3.D0-(8.D0)*TEMP*TEMP)*D2)                    
C     &                +(12.D0*TEMP*D3)+(6.D0*TEMP                      
C     &                *TEMP*D2*D2)-(8.D0*TEMP*BVPKR(N,1))))            
C             ENDIF                                                     
         ELSE                                                           
             DO 101 I=1,4                                               
                BVPKR(N,I)=0.D0                                         
 101                    CONTINUE                                                   
         ENDIF                                                          
 100      CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPSCB                                                 
C                                                                       
C Routine to initialize all common blocks for given BETA value          
C                                                                       
C Written by M Gibbs 17/10/92                                           
C                                                                       
      CALL BVPSBR                                                       
      CALL BVPSPD                                                       
      CALL BVPSZD                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPSNR                                                 
C                                                                       
C Routine to find saddle point using Newton-Rhapson iteration           
C **On exit,bessel fn. common blocks need to be reset.                  
C                                                                       
C OLDBET   Old value of beta in iteration process                       
C STYLE    .TRUE. indicates non-trivial mass config. for light particles
C                                                                       
C Routines used:BVPSBR Obtain bessel function values                    
C               BVPSPD Derivatives of integrand                         
C                                                                       
C Written by M Gibbs 13/4/92                                            
C                                                                       
C ** Contains no error checking                                         
C                                                                       
C 15/06/92 Modified to use better Bessel fn. calculating routine        
C 08/07/92 Now allows a non-trivial light particle config.to be used    
C 17/10/92 Modified to use SADDLE10.INC include file -STYLE removed     
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION TEMP,FIRST,RACCU,                                
     &      OLDBET,BETAS                                                
C      LOGICAL STYLE                                                    
      INTEGER I                                                         
      I=0                                                               
      RACCU=ACCU                                                        
C---Initial value of BETA                                               
      OLDBET=2.D0*DBLE(BVPTOT)/EN                                       
      TEMP=1200.*OLDBET                                                 
      BETAS=OLDBET                                                      
C---Main iteration loop for NR                                          
 10    IF (BETAS.GT.TEMP) THEN                                           
          BETAS=BETAS-1.D0                                              
          GOTO 10                                                       
      ENDIF                                                             
      IF (I.GT.25) THEN                                                 
          OLDBET=BETAS                                                  
          GOTO 101                                                      
      ENDIF                                                             
      IF (BETAS.LT.(1E-10)) THEN                                        
          BETAS=BETAS+1.D0                                              
          GOTO 10                                                       
      ENDIF                                                             
      OLDBET=BETAS                                                      
      BETA=OLDBET                                                       
      CALL BVPSCB                                                       
      FIRST=BVPDER(2)                                                   
      IF (ABS(FIRST).LT.RACCU) THEN                                     
         GOTO 101                                                       
      ENDIF                                                             
      BETAS=OLDBET-(FIRST/BVPDER(3))                                    
      I=I+1                                                             
      OLDBET=BETAS                                                      
      GOTO 10                                                           
C---Have now iterated to required accuracy                              
 101   CONTINUE                                                          
C      PRINT 876,NUMBER(2),OLDBET                                       
C      PRINT 877,I,ACCU                                                 
 877    FORMAT(' Number of Iterations:',I4,' ACCU value:',E25.10)         
 876     FORMAT(' Number of bosons:',I4,' Beta value:',E25.10)             
      BETA=OLDBET                                                       
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPSPD                                                 
C                                                                       
C Routine to evaluate derivatives of integrand at BETA                  
C                                                                       
C STYLE     .TRUE. indicates non-trivial mass config.                   
C                                                                       
C Routines used:S18CCF (NAG) Bessel fn. calculation I0                  
C               S18CDF (NAG) Bessel fn. calculation I1                  
C                                                                       
C Written by M Gibbs 15/6/92                                            
C                                                                       
C 07/08/92   Modified to allow use of a non-trivial mass config.        
C 17/10/92   Now uses SADDLE10.INC include file - STYLE option removed  
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION FACT(4)                                          
      INTEGER J,N                                                       
      DATA FACT/-1.D0,1.D0,-2.D0,6.D0/                                  
      DO 10 N=1,4                                                       
C---First the CofM terms                                                
        BVPDER(N+1)=(((2.D0*DBLE(BVPTOT))-1.D0)/(BETA**DBLE(N)))        
     &       *FACT(N)                                                   
C---Now Energy terms                                                    
        BVPDER(N+1)=BVPDER(N+1)+((EN**DBLE(N))*BVPIR(N))                
C---Now particle terms                                                  
        DO 20 J=1,NTYPE                                                 
           BVPDER(N+1)=BVPDER(N+1)-((BVPMAS(J)**DBLE(N))                
     &           *BVPKR(J,N)*DBLE(BVPNUM(J)))                           
 20           CONTINUE                                                        
C---New code for non-trivial mass configuration                         
C      IF (STYLE) THEN                                                  
C         DO 30 J=1,NUMBER(1)                                           
C            TMASS=P(5,J)                                               
C            TEMP=BETA*TMASS                                            
C            WRITE(*,*) 'BETA*MASS',TEMP,' Mass is',MASS                
C            IF (TEMP.LT.ALACC) TEMP=ALACC                              
C            IF ((TEMP).GT.(10D0)) THEN                                 
C               CALL BVPABR((-1.D0*TEMP),K(1,1),                        
C     &                 K(1,2),K(1,3),K(1,4))                           
C            ELSE                                                       
C                D0=S18CCF(TEMP,IFAIL)                                  
C                D1=S18CDF(TEMP,IFAIL)                                  
C                K(1,1)=D0/D1                                           
C                D2=K(1,1)*K(1,1)                                       
C                D3=D2*K(1,1)                                           
C                K(1,2)=(D2-1.D0)+(K(1,1)/TEMP)                         
C                K(1,3)=(2.D0*K(1,1)*(D2-1.D0))                         
C     &               +(3.D0*D2/TEMP)-1.D0                              
C                K(1,4)=2.D0+((1.D0/(TEMP*TEMP))*(1.D0                  
C     &                +((3.D0-(8.D0)*TEMP*TEMP)*D2)                    
C     &                +(12.D0*TEMP*D3)+(6.D0*TEMP                      
C     &                *TEMP*D2*D2)-(8.D0*TEMP*K(1,1))))                
C            ENDIF                                                      
C         BVPSPD=BVPSPD-((TMASS**DBLE(N))                               
C     &         *K(1,N))                                                
C                                                                       
C   30    CONTINUE                                                      
C      ELSE                                                             
C         BVPSPD=BVPSPD-((MASS(1)**DBLE(N))                             
C     &         *K(1,N)*DBLE(NUMBER(1)))                                
C      ENDIF                                                            
 10            CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPSPS                                                 
C                                                                       
C Routine to perform saddle point integration for BVI phase space       
C Returns log of integral multiplied by 1/4.pi.pi.En                    
C Multiply answer by (2pi)**(4-3N) to get req'd phase space             
C                                                                       
C STYLE     .TRUE. indicates mass distribution is to be taked from P    
C                                                                       
C Routines used:BVPSNR Saddle point location                            
C               BVPSPD Derivs. of integrand                             
C               BVPSZD Value of integrand                               
C               BVPSBR Calculate Bessel function values                 
C               BVPSSC Higher order corrections to saddle point int.    
C                                                                       
C Written by M Gibbs 13/4/92                                            
C                                                                       
C 05/05/92 Modified to use higher order corrections                     
C 15/06/92 Modified from BVPSPI to calc. phase space only.              
C 17/06/92 Extra BETA**N term added to correct calculation              
C 09/07/92 Extra flag included to allow use of light particles of       
C          varying mass                                                 
C 17/10/92 Modified to use SADDLE10.INC include file -STYLE flag removed
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION SEC                                                 
C---First find saddle point                                             
C      WRITE(*,*) 'now doing NR for ps calc'                            
      CALL BVPSNR                                                       
C---Now set bessel function common block up                             
      CALL BVPSCB                                                       
C---Now return value for integral                                       
      SEC=BVPDER(3)                                                     
      BVPPSI=BVPDER(1)-(2.75681559961402D0                              
     &      +(0.5D0*LOG(SEC))                                           
     &      +LOG(EN))                                                   
C---Now the beta term from saddle point integral                        
      BVPPSI=BVPPSI+(DBLE(BVPTOT)*(1.83787706640934D0                   
     &      -LOG(BETA)))                                                
C---Now saddle-point higher order corrections                           
C---This part NOT YET UPDATED                                           
      CALL BVPSSC                                                       
C      WRITE(*,*)'Correction of:',BVPPCO,' to value:',BVPPSI            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPSSC                                                 
C                                                                       
C Function to supply next-order correction to saddle point integral     
C Assumes common block of bessel fn. ratios is set up                   
C                                                                       
C Taken from MAMBO preprint (Kleiss and Stirling)                       
C CERN-TH.6293/91                                                       
C                                                                       
C SECOND    Second derivative of integrand                              
C                                                                       
C Written by M Gibbs 30/4/92                                            
C 15/06/92 Modified to use BVPSPD for phase space integral              
C 17/10/92 Now uses SADDLE10.INC include file                           
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION SECOND,TEMP,TEMP2                                
      SECOND=BVPDER(3)                                                  
      TEMP=SECOND*SECOND*8.D0                                           
      BVPPCO=1.0D0+(BVPDER(5)/TEMP)                                     
      TEMP=TEMP*SECOND*3.D0/5.D0                                        
      TEMP2=BVPDER(4)                                                   
      BVPPCO=BVPPCO-(TEMP2*TEMP2/TEMP)                                  
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE BVPSZD                                                 
C                                                                       
C Saddle point integrand for phase space integral                       
C                                                                       
C TEMP    Temporary store                                               
C STYLE   .TRUE. indicates P mass config. to be used for light particles
C                                                                       
C Routines iused:S18CFF (NAG) I1 modified bessel function (scaled)      
C                S18CDF (NAG) K1 modified bessel function (scaled)      
C                                                                       
C Written by M Gibbs 13/4/92                                            
C                                                                       
C 27/04/92 Modified to use own calls to NAG routines                    
C 15/06/92 Changed to calculate ps integral only                        
C  8/07/92 Variable light particle mass configs. added                  
C 17/10/92 SADDLE10.INC include file added - STYLE removed              
C  4/11/92 NAG routine dependence removed - now uses BVPA**             
C                                                                       
      IMPLICIT NONE                                                     
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION BVPAIL,BVPAKL
      INTEGER IFAIL,J                                                   
      IFAIL=0                                                           
C---First the CofM terms                                                
      BVPDER(1)=2.0D0*LOG(BETA)                                         
C NAG needed for next two lines:                                        
C      I1=S18CFF(BETA*EN,IFAIL)                                         
C      BVPDER(1)=BVPDER(1)+LOG(I1)+(BETA*EN)                            
      BVPDER(1)=BVPDER(1)+BVPAIL(BETA*EN)                               
C---Particle terms                                                      
      DO 20 J=1,NTYPE                                                   
         IF (BVPNUM(J).GT.0) THEN                                       
            BVPDER(1)=BVPDER(1)+(DBLE(BVPNUM(J))*                       
C NAG next three lines, own BVPA** three following                      
C     &             (LOG(S18CDF(BETA*BVPMAS(J),IFAIL))                  
C     &             +LOG(BVPMAS(J))                                     
C     &             -(BETA*BVPMAS(J))))                                 
     &             (LOG(BVPMAS(J))+                                     
     &             BVPAKL(BETA*BVPMAS(J))))                             
         ELSE                                                           
            CONTINUE                                                    
         ENDIF                                                          
 20       CONTINUE                                                          
C---Extra code for non-trivial mass config for light particles          
C      WRITE(*,*) 'BVPSZD eval. of particle functions'                  
C      IF (STYLE) THEN                                                  
C          DO 30 J=1,NUMBER(1)                                          
C            TMASS=P(5,J)                                               
C            BVPSZD=BVPSZD                                              
C     &             +(LOG(S18CDF(BETA*TMASS,IFAIL))                     
C     &             +LOG(TMASS)                                         
C     &             -(BETA*TMASS))                                      
C   30     CONTINUE                                                     
C      ELSE                                                             
C---Code as above for J=1 case                                          
C            BVPSZD=BVPSZD+(DBLE(NUMBER(1))*                            
C     &             (LOG(S18CDF(BETA*MASS(1),IFAIL))                    
C     &             +LOG(MASS(1))                                       
C     &             -(BETA*MASS(1))))                                   
C      ENDIF                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION CGEVPB(VOID)                                             
C                                                                       
C Function to supply log conversion factor for Gev-2 to pb              
C                                                                       
C Written by M Gibbs 24/3/92                                            
C                                                                       
      INTEGER VOID                                                      
      DOUBLE PRECISION CGEVPB                                           
      CGEVPB=+19.780065D0                                               
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
C This is the file herbvi10.for                                         
C                                                                       
C ********************************************************************* 
C                                                                       
C                                                                       
C ********************************************************************* 
C                                                                       
C                                                                       
C ********************************************************************* 
C                                                                       
C                                                                       
C ********************************************************************* 
C                                                                       
C                                                                       
C HERBVI10.FOR  -  BNV MC package for HERWIG5.7 onwards                 
C                                                                       
C M Gibbs 26/03/94                                                      
C                                                                       
C Requirements: HERWIG5.7                                               
C               SADDLE1.0                                               
C                                                                       
C Files needed: herbvi10.inc                                            
C               herwig57.for                                            
C               herwig57.inc                                            
C               saddle10.for                                            
C               saddle10.for                                            
C                                                                       
C Note that include statements disregard the version numbers            
C                                                                       
C                                                                       
C ********************************************************************* 
C                                                                       
C Routines included:                                                    
C                                                                       
C User-supplied routines (optional)                                     
C                                                                       
C HURBOS SD User-supplied boson distribution generation                 
C HUSGEN SD User-supplied energy distribution generation                
C                                                                       
C Analysis routines                                                     
C                                                                       
C HVAFIP S  Find isolated particles                                     
C HVANAL S  Initial analysis routine                                    
C HVASPN S  Sum particle numbers                                        
C                                                                       
C HERWIG interface routines                                             
C                                                                       
C HVCBVI S  Colour-connect quarks in HEPEVT record                      
C HVCBVT F  Find BVI-connected quarks in HEPEVT                         
C                                                                       
C Hard process generation                                               
C                                                                       
C HVETWO S  Set up 2->2+n hard subprocess                               
C                                                                       
C Hard subprocesses                                                     
C                                                                       
C HVHBVI S  BVI interface for HERWIG hard-subprocess                    
C HVHEVT S  Generate configuration for BVI subprocess                   
C HVHGEN S  Generate event data                                         
C HVHINI S  Initialise common block for BVI subprocess                  
C HVHMWE S  Phase-space generation for Multi-W configuration            
C HVHMWG S  Generate event data for Multi-W subprocess                  
C HVHPQM S  Project quarks through CKM matrix                           
C                                                                       
C Initialisation                                                        
C                                                                       
C HVINIT S  Initialise parameters                                       
C                                                                       
C Random number generation of distributions                             
C                                                                       
C HVRBOS S  Generate boson number distribution                          
C HVRBIN F  Generate binomial distribution                              
C                                                                       
C Structure functions and energy distributions                          
C                                                                       
C HVSFUN S  Generate x distribution                                     
C HVSGEN S  Generate \shat distribution                                 
C                                                                       
C Utility routines                                                      
C                                                                       
C HVUSPT F  Set array pointers                                          
C                                                                       
C Other routines                                                        
C                                                                       
C HWUPRA F  Calculate pseudo-rapidity                                   
C                                                                       
C Key: S Subroutine                                                     
C      F Function                                                       
C      D Dummy routine (replacement optional)                           
C      R Replacement (delete dummy version in HERWIG to use)            
C                                                                       
C ********************************************************************* 
C                                                                       
C Modification notes from HVDEVL code:                                  
C                                                                       
C 11/04/93 BVPMFC removed to allow use of saddev.for version            
C          Version in this file renamed QVPMFC                          
C  4/06/93 Routine HVPSXV removed, incorporated into HVINIT             
C 18/06/93 HVHPQM projection routine added                              
C 30/06/93 Considerable mods made to doc.                               
C          Routine HVASPN added                                         
C 21/07/93 Routine HVAFIP added                                         
C 12/10/93 Sun makefile version created                                 
C 13/10/93 HWUPRA replaced as bug fix                                   
C 14/10/93 Include file changes made for HW version number changes      
C 15/10/93 Flag HVFCEN added for BVI routines                           
C          Multi-W routines added                                       
C 18/10/93 Routine HVRBOS to gen boson distrib. added                   
C 22/03/94 Charge difference in boson distrib. parameterised            
C 26/03/94 Major modifications to strucutre for HERBVI package          
C          M8/60 for details. Routines removed: BVPPCS                  
C          Moved to saddle: BVPAKQ                                      
C          Include file renamed HERBVI10.INC / SADDLE10.INC             
C 02/04/95 Total parton level cross section parameter added             
C          Unix version generated.                                      
C          Uses makefile to produce object code.                        
C          Include statements in UNIX format                            
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HURBOS(CHDIFF)                                         
C                                                                       
C User-generated boson distribution generation routine                  
C                                                                       
C Called if HVCONT(8) is set .FALSE. by user                            
C                                                                       
C Written by M Gibbs 18/10/93                                           
C                                                                       
C 22/03/94 CHDIFF boson charge difference added                         
C 23/03/94 CHCONT control flag added                                    
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER CHDIFF                                                 
      CALL HWWARN('HURBOS',500,*999)                                    
 999    CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HUSGEN(WT)                                             
C                                                                       
C User-generated sf generation routine                                  
C                                                                       
C Called if HVCONT(6) is set .FALSE. by user                            
C                                                                       
C Written by M Gibbs 15/10/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      DOUBLE PRECISION WT                                               
      WT = 0.0D0                                                        
      CALL HWWARN('HUSGEN',500,*999)                                    
 999    CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVAFIP(IDAL,R,ENIN,ENAR)                               
C                                                                       
C Routine to find isolated particles of type IDAL in calorimeter        
C                                                                       
C Assumes calorimeter already loaded by HVANAL                          
C                                                                       
C Note that the calorimeter is bigger than our accepted particle        
C limits - but have cuts to make sure anyway                            
C Also, this routine sums a SQUARE around the particle                  
C                                                                       
C Written by M Gibbs 21/07/93                                           
C                                                                       
C On entry:                                                             
C IDAL   IDHEP of particles to be considered                            
C ENIN   Max. amount of energy in calorimeter other than particle       
C ENAR   Max. amount of energy around particle                          
C XD     Distance to sum in calorim X direction                         
C YD     Distance to sum in calorim Y direction                         
C                                                                       
C ISCIS() .true. for particle meeting isolation criterion               
C         This array must be zeroed elsewhere before use of HVAFIP      
C                                                                       
C 15/08/93 Bug fixed, now dependent on R parameter                      
C                                                                       
      INCLUDE 'HERWIG65.INC'                                         
      INCLUDE 'HERBVI.INC'                                                   
C      INCLUDE 'HERWIG57.INC'                                           
C      INCLUDE 'HERBVI10.INC'                                           
C      INCLUDE 'GETJET.INC'                                             
      INTEGER IDAL                
      DOUBLE PRECISION ENIN,ENAR,R
CC---Rapidity limit                                                     
CC      RLIM=3.0D0                                                      
CC---Calculate XD and YD for this R value                               
C      RSQ = R**2                                                       
CC---Angular granularity                                                
C      YD = INT(R/DELPHI)                                               
CC---Rapidity distribution                                              
C      XD = INT(R/DELY)                                                 
CC---Loop over all final state particles                                
C      DO 100 IPTR=1,FSPPTR                                             
C        IHEP=IDHNUM(IPTR)                                              
CC---Only consider correct particles                                    
C        IF (IDAL.NE.IDHEP(IHEP)) GOTO 100                              
CC---Only include accepted ones                                         
C        IF (.NOT.ISACC(IPTR)) GOTO 100                                 
CC---Initialise variables                                               
C        ENSUR=0.0D0                                                    
C        ENRIN=0.0D0                                                    
CC---Get calorimeter position                                           
C        CALX=CALPOS(1,IPTR)                                            
C        CALY=CALPOS(2,IPTR)                                            
CC---Reject particle if too close to rapidity limit                     
CC        RP=ABS(ATAN2(PHEP(2,IHEP),PHEP(1,IHEP)))                      
CC        IF (RP.                                                       
CC---Energy in calorimeter not due to this particle                     
C        ENRIN = ET(CALX,CALY) - PHEP(4,IHEP)                           
CC---Cut if this is too large                                           
CC        IF (ENRIN.GT.ENIN) GOTO 100                                   
CC---Sum surrounding energy                                             
C        ENSUR = -ET(CALX,CALY)                                         
CC---Loop in circular direction                                         
C        DO 10 YP=CALY-YD,CALY+YD                                       
C          YU=YP                                                        
CC---Enforce a circular calorimeter                                     
C          IF (YU.LT.1) YU=YU+NCPHI                                     
C          IF (YU.GT.NCPHI) YU=YU-NCPHI                                 
CC---Now loop over rapidity direction                                   
C          DO 20 XP=CALX-XD,CALX+XD                                     
CC---Check we are in the circle of radius R                             
CC      PRINT *,XP-CALX,YU-CALY                                         
C            T = ((DELY*DBLE(XP-CALX))**2)                              
C     &        + ((DELPHI*DBLE(YP-CALY))**2)                            
C            IF (T.LE.RSQ) ENSUR=ENSUR+ET(XP,YU)                        
CC      PRINT *,T,RSQ                                                   
C   20     CONTINUE                                                     
C   10   CONTINUE                                                       
CC---Surrounding energy cut                                             
C        IF (ENSUR.GT.ENAR) GOTO 100                                    
CC---Debug print                                                        
CC      PRINT *,IPTR,ENSUR,XD,YD                                        
CC---Particle has survived cuts, so is isolated                         
C        ISCIS(IPTR)=.TRUE.                                             
CC---Loop over all particles                                            
 100   CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVANAL                                                 
C                                                                       
C Routine to perform initial analysis on final state particles after BVI
C                                                                       
C This routine should only be called if IERROR is zero                  
C                                                                       
C Error messages (reported by a call to HWWARN in HERWIG):              
C                                                                       
C    1 Bad momentum sum for final state                                 
C    2 Bad charge sum for final state                                   
C  101 Number of final state particles exceeds size of arrays           
C                                                                       
C Written by M Gibbs 19/04/93                                           
C                                                                       
C 17/06/93 ISACC flag added                                             
C 21/06/93 ISLFZ flag added                                             
C          Identification of high energy leptons added                  
C 30/06/93 High energy tagging removed                                  
C          Added storage of Z0 if parent of lepton                      
C                                                                       
      INCLUDE 'HERWIG65.INC'                                         
      INCLUDE 'HERBVI.INC'                                                   
      COMMON/HWBVIC/NBV,IBV(18)                                         
      DOUBLE PRECISION HWVDOT,HWURAP,PSUM(4),HWUPRA             
      INTEGER NBV,IBV,JBV,ICHSUM,ICHINI,IHEP,JHEP,KHEP,ID             
C---Initialise pointers                                                 
      FSPPTR=1                                                          
      CALL HWVSUM(4,PHEP(1,1),PHEP(1,2),PSUM)                           
      CALL HWVSCA(4,-1D0,PSUM,PSUM)                                     
      ICHSUM=0                                                          
      ICHINI=ICHRG(IDHW(1))+ICHRG(IDHW(2))                              
C---Loop over all the particle types                                    
      DO 100 IHEP=1,NHEP                                                
        IF (ISTHEP(IHEP).EQ.1) THEN                                     
C---Add up the total momenta and charge of the outgoing particles       
          CALL HWVSUM(4,PHEP(1,IHEP),PSUM,PSUM)                         
          ICHSUM=ICHSUM+ICHRG(IDHW(IHEP))                               
          ID=IDHEP(IHEP)                                                
C---Calculate rapidities                                                
          PSRFSP(FSPPTR)=HWUPRA(PHEP(1,IHEP))                           
          RAPFSP(FSPPTR)=HWURAP(PHEP(1,IHEP))                           
C---Calculate transverse momenta                                        
          PTFSP(FSPPTR)=SQRT(PHEP(1,IHEP)**2+PHEP(2,IHEP)**2)           
C---Enter particle id's into the event record                           
          IDHFSP(FSPPTR) = ID                                           
          IDHNUM(FSPPTR) = IHEP                                         
          ISFSP(FSPPTR)  =.FALSE.                                       
          ISLFZ(FSPPTR)  =.FALSE.                                       
          ISACC(FSPPTR)  =.TRUE.                                        
          IDZPAR(FSPPTR) = 0                                            
C---Now search out particles from the BVI process                       
          IF (NBV.NE.0) THEN                                            
C---Is it a gamma?                                                      
            IF (IDHW(IHEP).EQ.59)  THEN                                 
C---Is it a primary gamma?                                              
              JHEP = IHEP                                               
 10                    IF ((IDHEP(JHEP).EQ.22).OR.(IDHEP(JHEP).EQ.94)            
     &             .OR.(IDHW(JHEP).EQ.59)) THEN                         
                IF (JMOHEP(1,JHEP).EQ.6) THEN                           
C---This is a primary gamma                                             
                  ISFSP(FSPPTR)=.TRUE.                                  
                ELSE                                                    
                  JHEP = JMOHEP(1,JHEP)                                 
                  GOTO 10                                               
                ENDIF                                                   
              ENDIF                                                     
C---Particle processed, continue on to next one                         
              GOTO 90                                                   
            ENDIF                                                       
C---Is this particle a Baryon?                                          
            IF (ID/1000.NE.0) THEN                                      
C---SEE IF IT CAME FROM A CREATED QUARK                                 
              JHEP=IHEP                                                 
 20                    JHEP=JMOHEP(1,JHEP)                                       
              IF (ISTHEP(JHEP).LT.163.OR.ISTHEP(JHEP).GT.186) GO TO 20  
              JHEP=JMOHEP(1,JHEP)                                       
              KHEP=JMOHEP(2,JHEP)                                       
              DO 30 JBV=1,NBV                                           
                IF (JHEP.EQ.IBV(JBV).OR.KHEP.EQ.IBV(JBV)) THEN          
C---This baryon is formed from the BVI process                          
                  ISFSP(FSPPTR)=.TRUE.                                  
C---Particle processed, go on to next one                               
                  GO TO 90                                              
                ENDIF                                                   
 30                      CONTINUE                                                  
C---Is the particle a lepton?                                           
            ELSEIF (ABS(ID/10).EQ.1) THEN                               
C---SEE IF IT CAME FROM PRIMARY PROCESS                                 
              JHEP=JMOHEP(1,JMOHEP(1,IHEP))                             
              IF (ISTHEP(JHEP).EQ.120) THEN                             
                ISFSP(FSPPTR)=.TRUE.                                    
              ENDIF                                                     
C---See if it is from a Z0                                              
              JHEP=IHEP                                                 
 31                    JHEP=JMOHEP(1,JHEP)                                       
              IF (IDHW(JHEP).EQ.200) THEN                               
                ISLFZ(FSPPTR) =.TRUE.                                   
                IDZPAR(FSPPTR)= JHEP                                    
              ENDIF                                                     
              IF (ABS(IDHEP(JHEP)/10).EQ.1) GOTO 31                     
C---Particle processed, move on to the next one                         
              GO TO 90                                                  
            ENDIF                                                       
          ENDIF                                                         
C---Loop over all particles                                             
 90            FSPPTR = FSPPTR + 1                                           
          IF (FSPPTR.GT.NPMAX) THEN                                     
            CALL HWWARN('HVANAL',101,*999)                               
          ENDIF                                                         
        ENDIF                                                           
 100       CONTINUE                                                        
      FSPPTR=FSPPTR-1                                                   
C---Checks on results - ensure momentum and charge is conserved         
      IF (HWVDOT(3,PSUM,PSUM).GT.1.E-4*PHEP(4,1)**2) THEN               
        CALL HWWARN('HVANAL',1,*101)                                    
 101       PRINT 102, PSUM                                                 
 102          FORMAT(' Bad Momentum Sum = ',4F9.3/)                           
        IERROR=1                                                        
      ENDIF                                                             
      IF (ICHSUM.NE.ICHINI) THEN                                        
        CALL HWWARN('HVANAL',2,*201)                                    
 201       PRINT 202, ICHSUM                                               
 202          FORMAT(' Bad Charge Sum = ',I5/)                                
        IERROR=2                                                        
      ENDIF                                                             
 999   CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVASPN(EM,TOT,HIGH,ID)                                 
C                                                                       
C Function to sum number of particles of type ID in output list         
C                                                                       
C Only includes particles with ISACC() .true.                           
C                                                                       
C Written by M Gibbs 30/06/93                                           
C                                                                       
      INCLUDE 'HERWIG65.INC'                                         
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER TOT,ID,I,HIGH                                             
      DOUBLE PRECISION EM,EL                                            
      TOT = 0                                                           
      EM = 0.0D0                                                        
C---Loop over all particles                                             
      DO 10 I=1,FSPPTR                                                  
        IF ((IDHFSP(I).EQ.ID).AND.ISACC(I)) THEN                        
          TOT = TOT + 1                                                 
          EL = PHEP(4,IDHNUM(I))                                        
          IF (EL.GT.EM) THEN                                            
            EM = EL                                                     
            HIGH = I                                                    
          ENDIF                                                         
        ENDIF                                                           
 10      CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVCBVI                                                 
C                                                                       
C FINDS UNPAIRED PARTONS AFTER BARYON-NUMBER VIOLATION                  
C                                                                       
C 15/10/93 upgraded for HW57 version                                    
C          Now calls HWDHQK itself and resets calling flag              
C                                                                       
      INCLUDE 'HERWIG65.INC'                                          
      COMMON/HWBVIC/NBV,IBV(18)                                         
      DOUBLE PRECISION HWRGEN,PDQ(5)                                    
      INTEGER NBV,IBV,JBV,KBV,LBV,IHEP,IP1,IP2,IP3,JP1,JP2,JP3,         
     & HVCBVT,ID,NBR,MBV,IQ1,IQ2,IQ3,ID1,ID2,IDQ,IDIQK(3,3)             
      LOGICAL SPLIT,TRIPL,DUNBV(18)                                     
      DATA IDIQK/111,110,113,110,109,112,113,112,114/                   
C---TRIPL IS .TRUE. FOR A QUARK OR ANTIDIQUARK                          
      TRIPL(ID)=ID.LT.7.OR.(ID.GT.114.AND.ID.LT.121)                    
C---Check for errors                                                    
      IF (IERROR.NE.0)  RETURN                                          
C      PRINT *,' processing event ',NEVHEP                              
C      PRINT *,' hvcbvi called, now calling HWDHQK'                     
C---Decay heavy quarks                                                  
C      CALL HWDHQK                                                      
C      PRINT *,' HWDHQK called, now calling HWCGSP'                     
C---Gluon splitting                                                     
      CALL HWCGSP                                                       
C---Reset bvi clustering flag                                           
      HVFCEN = .FALSE.                                                  
C---LIST PARTONS WITH WRONG COLOUR PARTNERS                             
 5     NBV=0                                                             
      DO 10 IHEP=1,NHEP                                                 
      IF (ISTHEP(IHEP).GT.149.AND.ISTHEP(IHEP).LT.155) THEN             
        IF (TRIPL(IDHW(IHEP))) THEN                                     
          IF (.NOT.TRIPL(IDHW(JMOHEP(2,IHEP)))) GO TO 10                
        ELSE                                                            
C---Extra check for Gamma's                                             
          IF (IDHW(IHEP).EQ.59) GO TO 10                                
C---End of bug fix.                                                     
          IF (TRIPL(IDHW(JDAHEP(2,IHEP)))) GO TO 10                     
        ENDIF                                                           
        NBV=NBV+1                                                       
        IF (NBV.GT.18) CALL HWWARN('HVCBVI',100,*999)                   
        IBV(NBV)=IHEP                                                   
        DUNBV(NBV)=.FALSE.                                              
      ENDIF                                                             
 10    CONTINUE                                                          
      IF (NBV.EQ.0) RETURN                                              
      IF (MOD(NBV,3).NE.0) CALL HWWARN('HVCBVI',101,*999)               
C      PRINT *,' HVCBVI has found ',NBV,' partons to process'           
C---PROCESS FOUND PARTONS, STARTING AT RANDOM POINT IN LIST             
      NBR=NBV*HWRGEN(0)                                                 
      DO 100 MBV=1,NBV                                                  
      JBV=MBV+NBR                                                       
      IF (JBV.GT.NBV) JBV=JBV-NBV                                       
      IF (.NOT.DUNBV(JBV)) THEN                                         
        DUNBV(JBV)=.TRUE.                                               
        IP1=IBV(JBV)                                                    
        JP1=HVCBVT(IP1)                                                 
C---FIND ASSOCIATED PARTONS                                             
        DO 20 KBV=1,NBV                                                 
        IF (.NOT.DUNBV(KBV)) THEN                                       
          IP2=IBV(KBV)                                                  
          JP2=HVCBVT(IP2)                                               
          IF (JP2.EQ.JMOHEP(2,JP1)) THEN                                
            DUNBV(KBV)=.TRUE.                                           
            DO 15 LBV=1,NBV                                             
            IF (.NOT.DUNBV(LBV)) THEN                                   
              IP3=IBV(LBV)                                              
              JP3=HVCBVT(IP3)                                           
              IF (JP3.EQ.JMOHEP(2,JP2)) THEN                            
                DUNBV(LBV)=.TRUE.                                       
                GO TO 25                                                
              ENDIF                                                     
            ENDIF                                                       
 15                CONTINUE                                                    
          ENDIF                                                         
        ENDIF                                                           
 20        CONTINUE                                                        
        print *,' Cant find partners for ',IBV(JBV)                     
        print *,' IBV follows'                                          
        print *,(IBV(KBV),KBV=1,NBV)                                    
        print *,(DUNBV(KBV),KBV=1,NBV)                                  
        CALL HWWARN('HVCBVI',102,*999)                                  
 25        IQ1=0                                                           
C        print *,' All partons have been found partners'                
C---LOOK FOR DIQUARK                                                    
        IF (ABS(IDHEP(IP1)).GT.100) THEN                                
          IQ1=IP1                                                       
          IQ2=IP2                                                       
          IQ3=IP3                                                       
        ELSEIF (ABS(IDHEP(IP2)).GT.100) THEN                            
          IQ1=IP2                                                       
          IQ2=IP3                                                       
          IQ3=IP1                                                       
        ELSEIF (ABS(IDHEP(IP3)).GT.100) THEN                            
          IQ1=IP3                                                       
          IQ2=IP1                                                       
          IQ3=IP2                                                       
        ENDIF                                                           
        IF (IQ1.EQ.0) THEN                                              
C---NO DIQUARKS: COMBINE TWO (ANTI)QUARKS                               
          IF (ABS(IDHEP(IP1)).GT.3) THEN                                
            IQ1=IP2                                                     
            IQ2=IP3                                                     
            IQ3=IP1                                                     
          ELSEIF (ABS(IDHEP(IP2)).GT.3) THEN                            
            IQ1=IP3                                                     
            IQ2=IP1                                                     
            IQ3=IP2                                                     
          ELSE                                                          
            IQ1=IP1                                                     
            IQ2=IP2                                                     
            IQ3=IP3                                                     
          ENDIF                                                         
          ID1=IDHEP(IQ1)                                                
          ID2=IDHEP(IQ2)                                                
C---CHECK FLAVOURS                                                      
          IF (ID1.GT.0.AND.ID1.LT.4.AND.                                
     &        ID2.GT.0.AND.ID2.LT.4) THEN                               
            IDQ=IDIQK(ID1,ID2)                                          
          ELSEIF (ID1.LT.0.AND.ID1.GT.-4.AND.                           
     &            ID1.LT.0.AND.ID2.GT.-4) THEN                          
            IDQ=IDIQK(-ID1,-ID2)+6                                      
          ELSE                                                          
C---CANT MAKE DIQUARKS WITH HEAVY QUARKS: TRY CLUSTER SPLITTING         
C            print *,' attmpting cluster splitting'                     
            CALL HWVSUM(4,PHEP(1,IQ1),PHEP(1,IQ2),PDQ)                  
            CALL HWUMAS(PDQ)                                            
C            print *,' calling HWCCUT for hvcbvi now'                   
            CALL HWCCUT(IQ1,IQ2,PDQ,.FALSE.,SPLIT)                      
C            print *,' hwccut finished, result is ',SPLIT               
            IF (SPLIT) GO TO 5                                          
C---Unable to form cluster; dispose of event                            
            CALL HWWARN('HVCBVI',-3,*999)                               
          ENDIF                                                         
C---OVERWRITE FIRST AND CANCEL SECOND                                   
          IDHW(IQ1)=IDQ                                                 
          IDHEP(IQ1)=IDPDG(IDQ)                                         
          CALL HWVSUM(4,PHEP(1,IQ1),PHEP(1,IQ2),PHEP(1,IQ1))            
          CALL HWUMAS(PHEP(1,IQ1))                                      
          ISTHEP(IQ2)=0                                                 
C---REMAKE COLOUR CONNECTIONS                                           
          IF (TRIPL(IDQ)) THEN                                          
            JMOHEP(2,IQ1)=IQ3                                           
            JDAHEP(2,IQ3)=IQ1                                           
          ELSE                                                          
            JDAHEP(2,IQ1)=IQ3                                           
            JMOHEP(2,IQ3)=IQ1                                           
          ENDIF                                                         
        ELSE                                                            
C---SPLIT A DIQUARK                                                     
          NHEP=NHEP+1                                                   
          CALL HWVSCA(5,HALF,PHEP(1,IQ1),PHEP(1,IQ1))                   
          CALL HWVEQU(5,PHEP(1,IQ1),PHEP(1,NHEP))                       
          ISTHEP(NHEP)=150                                              
          JMOHEP(1,NHEP)=JMOHEP(1,IQ1)                                  
          JDAHEP(1,NHEP)=0                                              
C---FIND FLAVOURS                                                       
          IDQ=IDHW(IQ1)                                                 
          DO 30 ID2=1,3                                                 
          DO 30 ID1=1,3                                                 
          IF (IDIQK(ID1,ID2).EQ.IDQ) THEN                               
            IDHW(IQ1)=ID1                                               
            IDHW(NHEP)=ID2                                              
C---REMAKE COLOUR CONNECTIONS (DIQUARK)                                 
            JMOHEP(2,IQ1)=IQ2                                           
            JMOHEP(2,IQ2)=NHEP                                          
            JMOHEP(2,IQ3)=IQ1                                           
            JMOHEP(2,NHEP)=IQ3                                          
            JDAHEP(2,IQ1)=IQ3                                           
            JDAHEP(2,IQ2)=IQ1                                           
            JDAHEP(2,IQ3)=NHEP                                          
            JDAHEP(2,NHEP)=IQ2                                          
            GO TO 35                                                    
          ELSEIF (IDIQK(ID1,ID2).EQ.IDQ+6) THEN                         
            IDHW(IQ1)=ID1+6                                             
            IDHW(NHEP)=ID2+6                                            
C---REMAKE COLOUR CONNECTIONS (ANTIDIQUARK)                             
            JMOHEP(2,IQ1)=IQ3                                           
            JMOHEP(2,IQ2)=IQ1                                           
            JMOHEP(2,IQ3)=NHEP                                          
            JMOHEP(2,NHEP)=IQ2                                          
            JDAHEP(2,IQ1)=IQ2                                           
            JDAHEP(2,IQ2)=NHEP                                          
            JDAHEP(2,IQ3)=IQ1                                           
            JDAHEP(2,NHEP)=IQ3                                          
            GO TO 35                                                    
          ENDIF                                                         
 30            CONTINUE                                                      
          CALL HWWARN('HVCBVI',104,*999)                                
 35            IDHEP(IQ1)=IDPDG(IDHW(IQ1))                                   
          IDHEP(NHEP)=IDPDG(IDHW(NHEP))                                 
        ENDIF                                                           
      ENDIF                                                             
 100   CONTINUE                                                          
 999    END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION HVCBVT(IP)                                               
C                                                                       
C     FINDS THE HARD PARTON FROM WHICH PARTON IP ORIGINATED,            
C     SEARCHING BACK THROUGH WEAK DECAYS IF NECESSARY                   
C                                                                       
      INCLUDE 'HERWIG65.INC'                                       
      INTEGER HVCBVT,IP,JP,KP                                           
      JP=IP                                                             
 10    KP=JMOHEP(1,JP)                                                   
      IF (ISTHEP(KP).EQ.120) THEN                                       
        HVCBVT=JP                                                       
        RETURN                                                          
      ELSE                                                              
        JP=KP                                                           
        GO TO 10                                                        
      ENDIF                                                             
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVETWO                                                 
C                                                                       
C Set up hard subprocess 2->2+n particles                               
C                                                                       
C Taken from HW listing 08/07/93 by M Gibbs.                            
C                                                                       
C Documented 20/07/93                                                   
C                                                                       
C Loads particle types according to IDN(4) contents                     
C                                                                       
C On exit:Outgoing quark momenta unset                                  
C         NHEP   Points to second outgoing quark                        
C                                                                       
C 20/07/93 No longer generates PS for outgoing particles                
C          Leaves CofM momenta in ICMF, still sets up colour conn.      
C 15/10/93 Added to hvdevl.F on 1 version                               
C                                                                       
      INCLUDE 'HERWIG65.INC'                                          
      INTEGER ICMF,IBM,I,J,K,IHEP                                       
      DOUBLE PRECISION PA
C---INCOMING LINES                                                      
      ICMF=NHEP+3                                                       
      DO 15 I=1,2                                                       
        IBM=I                                                           
C---FIND BEAM AND TARGET                                                
        IF (JDAHEP(1,I).NE.0) IBM=JDAHEP(1,I)                           
        IHEP=NHEP+I                                                     
        IDHW(IHEP)=IDN(I)                                               
        IDHEP(IHEP)=IDPDG(IDN(I))                                       
        ISTHEP(IHEP)=110+I                                              
        JMOHEP(1,IHEP)=ICMF                                             
        JMOHEP(I,ICMF)=IHEP                                             
        JDAHEP(1,IHEP)=ICMF                                             
        PHEP(1,IHEP)=0.                                                 
        PHEP(2,IHEP)=0.                                                 
        PHEP(5,IHEP)=RMASS(IDN(I))                                      
        PA=XX(I)*(PHEP(4,IBM)+ABS(PHEP(3,IBM)))                         
        PHEP(4,IHEP)=0.5*(PA+PHEP(5,IHEP)**2/PA)                        
 15        PHEP(3,IHEP)=PA-PHEP(4,IHEP)                                    
      PHEP(3,NHEP+2)=-PHEP(3,NHEP+2)                                    
C---HARD CENTRE OF MASS                                                 
      IDHW(ICMF)=IDCMF                                                  
      IDHEP(ICMF)=IDPDG(IDCMF)                                          
      ISTHEP(ICMF)=110                                                  
C---Sum incoming momentum, place into ICMF                              
      CALL HWVSUM(4,PHEP(1,NHEP+1),PHEP(1,NHEP+2),PHEP(1,ICMF))         
C---Calculate invariant mass for ICMF                                   
      CALL HWUMAS(PHEP(1,ICMF))                                         
C---OUTGOING LINES                                                      
      DO 20 I=3,4                                                       
        IHEP=NHEP+I+1                                                   
        IDHW(IHEP)=IDN(I)                                               
        IDHEP(IHEP)=IDPDG(IDN(I))                                       
        ISTHEP(IHEP)=110+I                                              
        JMOHEP(1,IHEP)=ICMF                                             
        JDAHEP(I-2,ICMF)=IHEP                                           
 20        PHEP(5,IHEP)=RMASS(IDN(I))                                      
C---SET UP COLOUR STRUCTURE LABELS                                      
      DO 30 I=1,4                                                       
        J=I                                                             
        IF (J.GT.2) J=J+1                                               
        K=ICO(I)                                                        
        IF (K.GT.2) K=K+1                                               
        JMOHEP(2,NHEP+J)=NHEP+K                                         
 30        JDAHEP(2,NHEP+K)=NHEP+J                                         
      NHEP=NHEP+5                                                       
 999   END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVHBVI                                                 
C                                                                       
C Baryon-number violating processes interface routine for HERWIG        
C                                                                       
C GENEV  .true.  indicates event previously generated to be loaded      
C FIRST  .true.  first call - variables to be initialised               
C        .false. phase space generation to reply with event weight      
C                                                                       
C Routines used:HVHGEN Generate event                                   
C               HVHINI Initialise all variables                         
C               HVHEVT Generate phase space and event weight            
C                                                                       
C Taken from original HWHBVI by M Gibbs 21/01/93                        
C                                                                       
C IPRO calls : 70 BVI event generation using HVHGEN HVHEVT              
C              71 Multi-W events             HVHMWG HVHMWE              
C                                                                       
C 15/10/93 name changed to hVhbvi for compatibility                     
C          calls to multi-w process added                               
C                                                                       
      INCLUDE 'HERWIG65.INC'                                           
      LOGICAL FIRST,HVCA70,HVCA71,HWRLOG                                
      COMMON /HVCETD/ HVCA70,HVCA71                                     
      DATA FIRST /.TRUE./                                               
      EXTERNAL HWRLOG                                                   
C---Main calling structure                                              
      IF (GENEV) THEN                                                   
C---Call phase space generators                                         
        IF (HVCA70) THEN                                                
C---BVI process                                                         
          CALL HVHGEN                                                   
        ELSEIF (HVCA71) THEN                                            
C---Multi-W process                                                     
          CALL HVHMWG                                                   
        ENDIF                                                           
      ELSE                                                              
        IF (FIRST) THEN                                                 
C---Initialise variables if not done so already                         
          CALL HVHINI                                                   
          FIRST = .FALSE.                                               
        ENDIF                                                           
C---Select event type                                                   
        IF (IPRO.EQ.70) HVCA70=.TRUE.                                   
        IF (IPRO.EQ.71) HVCA70=.FALSE.                                  
        IF (IPRO.EQ.72) HVCA70=HWRLOG(HALF)                             
        HVCA71 = (.NOT.HVCA70)                                          
C---Select event type                                                   
        IF (HVCA70) THEN                                                
C---BVI record insertion                                                
          CALL HVHEVT                                                   
        ELSEIF (HVCA71) THEN                                            
C---Multi-w 2->2+n process                                              
          CALL HVHMWE                                                   
        ENDIF                                                           
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVHEVT                                                 
C                                                                       
C Routine to generate ps configuration and return event weight          
C                                                                       
C Routines used:MAMBO  Phase space generator                            
C                                                                       
C Taken from HWHBVI by M Gibbs 21/01/93                                 
C                                                                       
C Error codes:100 Non-BVI event generation attempted                    
C                                                                       
C 22/01/93 HVDEVL.INC include file now required                         
C          Now uses parameters HVPNWP,HVPNWM,HVPNZO                     
C          Modified to use BVPMAM phase space generator(in SADDLE)      
C 16/03/93 Now generates x values using bvi distribution                
C 29/03/93 STARTP to parametrize event location added                   
C          HVNTOT is total number of particles in BVI subprocess        
C 11/04/93 Gamma particles now included                                 
C          Uses BVPNGA to generate Gamma/Z distribution                 
C 12/04/93 Energy Common block COFMEN added                             
C 14/04/93 Nboson estimated using a gaussian approx.                    
C 15/04/93 FLAG for event type now set at random                        
C 26/04/93 BVRBIN generator for binomial distrib. introduced            
C 24/05/93 MAMBO event generation moved to HVHGEN                       
C  4/06/93 COFMEN moved into HVDEVL.INC file                            
C 18/06/93 Only contrib. to weight now SF part.                         
C 22/06/93 MAMBO event gen. call controlled by HVCONT(1)                
C          HVCONT(2) .true. for a BVI event                             
C 29/06/93 MAMBO calls improved by preventing reinitialisation          
C 21/07/93 HVCONT(2) .false. now generates error message                
C 04/08/93 Higgs generation added using a rough approx.                 
C 15/10/93 Boson distribution now 1/alpha_w approx                      
C 18/10/93 HVRBOS now used to generate boson distributions              
C 22/03/94 CHDIFF (charge diff. in bosons) now calculated               
C 23/03/94 Added estimate of BETA from PS integral for fermion distrib. 
C 02/04/95 Total cross section constant added                           
C                                                                       
      INCLUDE 'HERWIG65.INC'                                         
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION WT,LEVWGT,AXWT,PSFACT,PSBETA,SHAT,BVPAKQ                   
      DOUBLE PRECISION MASST,CNB,M1,M2,RELPR
      LOGICAL HWRLOG,FLAG(4),PRTEN                                      
      INTEGER J,I,MINIT,NPS,CHDIFF,N1                            
      EXTERNAL BVPAKQ                                                   
      COMMON /BVIMAM/ MINIT,PRTEN                                     
      MINIT = 0                                                         
C---This is a BVI event                                                 
      HVCONT(2)=.TRUE.                                                  
C---Energy for colliding beams                                          
      COFMEN = PHEP(5,3)                                                
C---Generate the SHAT distribution                                      
      CALL HVSGEN(AXWT)                                                 
      SHAT = EMSCA*EMSCA                                                
C---Generate incoming quark types                                       
      DO 29 I=1,2                                                       
        FLAG(I)=HWRLOG(HALF)                                            
 29      CONTINUE                                                          
C---Initial estimate of boson number                                    
      CHDIFF = 0                                                        
      CALL HVRBOS(CHDIFF,.TRUE.)                                        
C---Sum number of bosons                                                
      CNB = DBLE(HVPNWP+HVPNWM+HVPNZO)                                  
C---Set up arrays for interacting qq pair                               
      DO 30 I=1,2                                                       
        BPTYPE(I) = 2                                                   
        IF (FLAG(I)) THEN                                               
          BPTYPE(I)=1                                                   
          CHDIFF = CHDIFF - 1                                           
        ENDIF                                                           
 30      CONTINUE                                                          
C---Now set up remaining fermion distribution                           
      CHDIFF = CHDIFF + 6                                               
C---Clear config array                                                  
      DO 131 I=1,NCONF                                                  
        CONFIG(I)=0                                                     
 131      CONTINUE                                                          
C---Estimate beta factor from PS integral                               
      PSFACT = EMSCA - (CNB*80.0D0)                                     
      PSBETA = (1.5D0*CNB)/PSFACT                                       
C---Lepton selection                                                    
      N1 = 7                                                            
      DO 31 I=1,3                                                       
        IF (HWRLOG(HALF)) THEN                                          
          CHDIFF = CHDIFF - 1                                           
          CONFIG(N1) = CONFIG(N1) + 1                                   
        ELSE                                                            
          CONFIG(N1+1) = CONFIG(N1+1) + 1                               
        ENDIF                                                           
        N1 = N1 + 2                                                     
 31      CONTINUE                                                          
C---Lightest quark generation                                           
      IF (HWRLOG(HALF)) THEN                                            
        CHDIFF = CHDIFF -1                                              
        CONFIG(1) = CONFIG(1) + 1                                       
      ELSE                                                              
        CONFIG(2) = CONFIG(2) + 1                                       
      ENDIF                                                             
C---Heavy quark generations                                             
      N1 = 3                                                            
      DO 33 I=1,2                                                       
        IF (I.EQ.2) THEN                                                
C---Get masses for bessel estimate of rel prob.                         
          M1 = RMASS(11)                                                
          M2 = RMASS(12)                                                
          RELPR = BVPAKQ(M1,M2,PSBETA)                                  
        ELSE                                                            
C---Set prob to be a half                                               
          RELPR = HALF                                                  
        ENDIF                                                           
C---Generate distribution                                               
        DO 32 J=1,3                                                     
          IF (HWRLOG(RELPR)) THEN                                       
            CHDIFF = CHDIFF - 1                                         
            CONFIG(N1) = CONFIG(N1) +1                                  
          ELSE                                                          
            CONFIG(N1+1) = CONFIG(N1+1) + 1                             
          ENDIF                                                         
 32          CONTINUE                                                        
        N1 = N1 + 2                                                     
 33      CONTINUE                                                          
CC---Set up CONFIG array for this event type                            
C      CALL BVPMFI(FLAG)                                                
C---Project down-type quarks through CKM matrix                         
      IF (HVCONT(3)) CALL HVHPQM                                        
C---Modify the boson number distribution by CHDIFF                      
      CALL HVRBOS(CHDIFF,.FALSE.)                                       
C---Enter number of bosons into CONFIG                                  
      CONFIG(13)=HVPNWP                                                 
      CONFIG(14)=HVPNWM                                                 
      CONFIG(15)=HVPNZO                                                 
      CONFIG(16)=HVPNGA                                                 
      CONFIG(17)=HVPNHI                                                 
C---Place particles into the event record                               
      CALL BVPMFC(STARTP,NPS,MASST)                                     
      CALL BVPINI                                                       
      HVNTOT = BVPTOT                                                   
C---Now generate the x distributions                                    
      CALL HVSFUN(WT,BPTYPE(1),BPTYPE(2))                               
C---If this gives zero then set arb. low result                         
      IF (WT.LE.ZERO) THEN                                              
        EVWGT = ZERO                                                    
      ELSE                                                              
        LEVWGT = LOG(WT) + LOG(AXWT)                                    
C---Now generate event weight                                           
        EVWGT = EXP(LEVWGT)                                             
C        PRINT *,WT,AXWT,EVWGT                                          
CC---Only need SF part as remainder is invariant                        
C        EVWGT=WT                                                       
C---Multiply in the extra weight factor                                 
        EVWGT=EVWGT*HVPLCS                                              
C        PRINT *,EVWGT                                                  
        LEVWGT=LEVWGT+LOG(HVPLCS)                                       
      ENDIF                                                             
C---COMPUTE AN EVENT WEIGHT                                             
        IF (EVWGT.NE.ZERO) THEN                                         
C---MAMBO phase space generator                                         
        IF (.NOT.HVCONT(1)) THEN                                        
          MINIT = 0                                                     
 100            CALL BVPMAM(HVNTOT,EMSCA,PHEP(1,STARTP),WT)                   
C---Make sure we have an unweighted ps distribution                     
          MINIT = 1                                                     
          IF (.NOT.HWRLOG(WT)) GOTO 100                                 
        ENDIF                                                           
      ENDIF                                                             
 999   RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVHGEN                                                 
C                                                                       
C Generate event data for previously generated phase space configuration
C                                                                       
C Taken from old HWHBVI by M Gibbs 21/01/93                             
C                                                                       
C Error codes:100 Non-BVI event generation attempted                    
C                                                                       
C 22/01/93 HVDEVL.INC include file now required                         
C          Now uses parameters HVPNWP,HVPNWM,HVPNZO                     
C 29/03/93 STARTP variable added - start position in record             
C 24/05/93 MAMBO event generator added from HVHEVT                      
C 22/06/93 HVCONT flag now controls useage of MAMBO event gen.          
C 29/06/93 MAMBO control flags included to speed up calling.            
C 21/07/93 HVCONT(2) .false. now generates an error message             
C 15/10/93 HVFCEN flag now set for HW57 upgrade                         
C                                                                       
      INCLUDE 'HERWIG65.INC'                                         
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER I,J,P1,P2,P3                                           
      LOGICAL HWRLOG                                                    
      INTEGER MINIT                                                     
      LOGICAL PRTEN                                                     
      COMMON /BVIMAM/ MINIT,PRTEN                                       
      DOUBLE PRECISION WT                                               
C---Error if not a BVI process                                          
      IF (.NOT.HVCONT(2)) CALL HWWARN('HVHGEN',100,*999)                
C---GENERATE EVENT                                                      
C---MAMBO phase space generator                                         
      MINIT = 0                                                         
      IF (HVCONT(1)) THEN                                               
 90          CALL BVPMAM(HVNTOT,EMSCA,PHEP(1,STARTP),WT)                     
C---Make sure we have an unweighted ps distribution                     
        MINIT = 1                                                       
        IF (.NOT.HWRLOG(WT)) GOTO 90                                    
      ENDIF                                                             
C---Set BVI flag for HVCBVI process control                             
      HVFCEN = .TRUE.                                                   
C---Select interacting quark types                                      
      DO 1 I=1,2                                                        
        IDN(I)=BPTYPE(I)                                                
 1       CONTINUE                                                          
C---Set id for CofM subprocess                                          
      IDCMF=15                                                          
      NHEP=3                                                            
C---Set up 2->1 hard subprocess in HEPEVT                               
      CALL HWEONE                                                       
C---Locations of the particles in hard subprocess                       
      P1=NHEP-2                                                         
      P2=NHEP-1                                                         
C---CofM location for hard subprocess                                   
      P3=NHEP                                                           
      STARTP=NHEP+1                                                     
C---Set counter for number of particles                                 
      NHEP=STARTP+HVNTOT-1                                              
C---Set up Mother/Daughter structure                                    
      JDAHEP(1,P3)=STARTP                                               
      JDAHEP(2,P3)=NHEP                                                 
C---Set status codes for HW                                             
      DO 5 I=STARTP,NHEP                                                
        ISTHEP(I)=114                                                   
        JMOHEP(1,I)=P3                                                  
        JMOHEP(2,I)=I                                                   
 5         JDAHEP(2,I)=I                                                   
      ISTHEP(STARTP)=113                                                
C---CHOOSE COLOUR CONNECTIONS                                           
C---First for the light quarks                                          
      IF (HWRLOG(HALF)) THEN                                            
        JMOHEP(2,P1)     = P2                                           
        JMOHEP(2,P2)     = STARTP                                       
        JMOHEP(2,STARTP) = P1                                           
        JDAHEP(2,P1)     = STARTP                                       
        JDAHEP(2,P2)     = P1                                           
        JDAHEP(2,STARTP) = P2                                           
      ELSE                                                              
        JMOHEP(2,P1)     = STARTP                                       
        JMOHEP(2,P2)     = P1                                           
        JMOHEP(2,STARTP) = P2                                           
        JDAHEP(2,P1)     = P2                                           
        JDAHEP(2,P2)     = STARTP                                       
        JDAHEP(2,STARTP) = P1                                           
      ENDIF                                                             
C---Colour connections for the other two generations                    
      DO 6 I=1,2                                                        
        J=STARTP-2+(3*I)                                                
        JMOHEP(2,J)   = J+1                                             
        JMOHEP(2,J+1) = J+2                                             
        JMOHEP(2,J+2) = J                                               
        JDAHEP(2,J)   = J+2                                             
        JDAHEP(2,J+1) = J                                               
        JDAHEP(2,J+2) = J+1                                             
 6        CONTINUE                                                          
C---Rotate particles into the Lab frame                                 
      DO 10 I=STARTP,NHEP                                               
        CALL HWULOB(PHEP(1,P3),PHEP(1,I),PHEP(1,I))                     
 10      CONTINUE                                                          
C---Initialise variables for boson decays                               
      DO 12 I=STARTP+10,NHEP                                            
        RHOHEP(1,I)=ONE                                                 
        RHOHEP(2,I)=ONE                                                 
        RHOHEP(3,I)=ONE                                                 
 12      CONTINUE                                                          
C---Set selective decays of bosons (default=any)                        
      DO 15 I=1,MODMAX                                                  
        MODBOS(I)=0                                                     
 15      CONTINUE                                                          
 999       END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVHINI                                                 
C                                                                       
C Initialise variables for BVI event generation                         
C                                                                       
C Taken from HWHBVI by M Gibbs 21/01/93                                 
C                                                                       
C 22/01/93 HVDEVL.INC include file now required                         
C          Now uses parameters HVPNWP,HVPNWM,HVPNZO                     
C          Changed to initialise MAMBO from SADDLE routine set          
C          (currently set to require reinitialization each call)        
C 26/01/93 Light,heavy particle masses now summed                       
C 16/03/93 Min, max X values now determined using HVPSXV                
C 27/03/93 Use of BVPMFC to set up event codes added                    
C 29/03/93 Event set up transferred to HVHEVT                           
C 12/04/93 New parameter set up included as HVINIT                      
C 24/05/93 Parameters for boson polynomial fits added                   
C  4/06/93 Routine HVPSXV removed, incorporated into HVINIT             
C 18/06/93 CKM matrix addition added                                    
C  4/08/93 Calc. of s variables for HVSGEN added                        
C                                                                       
      INCLUDE 'HERWIG65.INC'                                           
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER MINIT,I,J                                                 
      DOUBLE PRECISION SUPPER,SLOWER,SCOMB                              
      COMMON /HVCSEV/ SUPPER,SLOWER,SCOMB                               
      LOGICAL PRTEN                                                     
      COMMON /BVIMAM/ MINIT,PRTEN                                       
C---Parameters for boson polynomial fits                                
      DATA HVGAUA/0.478D-1,-0.107D-5,0.434D-1,-0.691D-6,                
     &            0.483D-1,-0.854D-6,0.474D-1,-0.943D-6/                
      DATA HVGAUP/-20.775D0,0.585D-2,-22.531D0,0.598D-2,                
     &            -21.497D0,0.590D-2,-22.847D0,0.597D-2/                
C---Set MAMBO to noprint and not initialised                            
      PRTEN = .FALSE.                                                   
      MINIT = 0                                                         
C---Set position in event table                                         
      STARTP = 7                                                        
C---Set parameter dependent variables if not done already               
      CALL HVINIT                                                       
C---Calculate cumulative probability matrix for CKM                     
      DO 10 I=1,3                                                       
        HVCKMC(I,1)=HVCKM(I,1)                                          
        DO 20 J=2,3                                                     
          HVCKMC(I,J)=HVCKM(I,J)+HVCKMC(I,J-1)                          
 20          CONTINUE                                                        
 10           CONTINUE                                                          
C---Calculate S variables for HVSGEN                                    
      SUPPER=EXP(SLMAX)                                                 
      SLOWER=EXP(SLMIN)                                                 
      SCOMB =(SUPPER-SLOWER)/(SUPPER*SLOWER)                            
 999   RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVHMWE                                                 
C                                                                       
C Routine to generate ps configuration and return event weight          
C                                                                       
C Routines used:MAMBO  Phase space generator                            
C                                                                       
C Taken from HVHEVT by M Gibbs 21/07/93                                 
C                                                                       
C Error codes:100 Non-BVI event generation attempted                    
C                                                                       
C 23/07/93 Colour connections in ICO now chosen randomly                
C 05/08/93 HVNH variable added for Higgs production                     
C 18/10/93 HVRBOS now used to generate boson number configs             
C 22/03/94 CHDIFF now added to HVRBOS call                              
C          Charge imbalance now handled by this method                  
C 02/04/95 Total parton level cs HVPLCS now implemented                 
C                                                                       
      INCLUDE 'HERWIG65.INC'                                            
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION WT,LEVWGT,AXWT                                            
      DOUBLE PRECISION MASST       
      LOGICAL HWRLOG,FLAG(4),PRTEN                                      
      INTEGER J,I,ID,MINIT,NPS,IHEP,CHDIFF                                
      COMMON /BVIMAM/ MINIT,PRTEN                                       
      MINIT = 0                                                         
C---Set up non-BVI flag                                                 
      HVCONT(2)=.FALSE.                                                 
C---Energy for colliding beams                                          
      COFMEN = PHEP(5,3)                                                
C---Generate the SHAT distribution                                      
      CALL HVSGEN(AXWT)                                                 
C---Generate quark types                                                
      DO 29 I=1,2                                                       
        FLAG(I)=HWRLOG(HALF)                                            
        FLAG(I+2)=FLAG(I)                                               
 29      CONTINUE                                                          
C---Initialise charge difference count                                  
      CHDIFF = 0                                                        
C---Set up arrays for interacting qq pair (Should swap at random)       
      DO 30 I=1,2                                                       
        BPTYPE(I) = 2                                                   
        IF(FLAG(I)) BPTYPE(I)=1                                         
        BPTYPE(I+2)=BPTYPE(I)                                           
 30      CONTINUE                                                          
C---Allow decays of quarks                                              
      DO 41 I=3,4                                                       
C---Does it require decaying?                                           
        IF (HWRLOG(HALF)) THEN                                          
          IF (BPTYPE(I).EQ.1) THEN                                      
C---Change d -> u W-                                                    
            BPTYPE(I) = 2                                               
            CHDIFF = CHDIFF - 1                                         
          ELSEIF (BPTYPE(I).EQ.2) THEN                                  
C---Change u -> d W+                                                    
            BPTYPE(I) = 1                                               
            CHDIFF = CHDIFF + 1                                         
          ENDIF                                                         
        ENDIF                                                           
 41      CONTINUE                                                          
C---Interchange these at random                                         
      IF (HWRLOG(HALF)) THEN                                            
        J=BPTYPE(1)                                                     
        BPTYPE(1)=BPTYPE(2)                                             
        BPTYPE(2)=J                                                     
      ENDIF                                                             
      IF (HWRLOG(HALF)) THEN                                            
        J=BPTYPE(3)                                                     
        BPTYPE(3)=BPTYPE(4)                                             
        BPTYPE(4)=J                                                     
      ENDIF                                                             
C---Estimate the boson number distribution                              
      CALL HVRBOS(CHDIFF,.TRUE.)                                        
C---And modify it by CHDIFF                                             
      CALL HVRBOS(CHDIFF,.FALSE.)                                       
C---Clear up CONFIG array                                               
      DO 27 I=1,NCONF                                                   
        CONFIG(I)=0                                                     
 27      CONTINUE                                                          
C---Set up CONFIG array for the quarks                                  
      DO 28 I=1,2                                                       
        IDN(I)=BPTYPE(I)                                                
        IDN(I+2)=BPTYPE(I+2)                                            
C---Set qbar for correct masses - IDHW corrected below                  
        CONFIG(BPTYPE(I+2))=CONFIG(BPTYPE(I+2))+1                       
 28      CONTINUE                                                          
C---Now set up ICO for colour connections                               
      IF (HWRLOG(HALF)) THEN                                            
        ICO(1)=3                                                        
        ICO(2)=4                                                        
        ICO(3)=1                                                        
        ICO(4)=2                                                        
      ELSE                                                              
        ICO(1)=4                                                        
        ICO(2)=3                                                        
        ICO(3)=2                                                        
        ICO(4)=1                                                        
      ENDIF                                                             
C---Enter number of bosons into CONFIG                                  
      CONFIG(13)=HVPNWP                                                 
      CONFIG(14)=HVPNWM                                                 
      CONFIG(15)=HVPNZO                                                 
      CONFIG(16)=HVPNGA                                                 
      CONFIG(17)=HVPNHI                                                 
C---Set STARTP for placing particles                                    
      STARTP=7                                                          
C---Place particles into the event record                               
      CALL BVPMFC(STARTP,NPS,MASST)                                     
C---Correct quark types for two outgoing                                
      DO 40 I=1,2                                                       
        IHEP=STARTP+I-1                                                 
        ID=BVPNUM(I+2)                                                  
        IDHW(IHEP)=ID                                                   
        IDHEP(IHEP)=IDPDG(ID)                                           
 40      CONTINUE                                                          
      CALL BVPINI                                                       
      HVNTOT = BVPTOT                                                   
C---Now generate the x distributions                                    
      CALL HVSFUN(WT,BPTYPE(1),BPTYPE(2))                               
C---If this gives zero then set arb. low result                         
      IF (WT.LE.ZERO) THEN                                              
        EVWGT = ZERO                                                    
      ELSE                                                              
        LEVWGT = LOG(WT) + LOG(AXWT)                                    
C---Now generate event weight                                           
        EVWGT = EXP(LEVWGT)                                             
CC---Only need SF part as remainder is invariant                        
C        EVWGT=WT                                                       
C---Multiply in the extra weight factor                                 
        EVWGT=EVWGT*HVPLCS                                              
        LEVWGT=LEVWGT+LOG(HVPLCS)                                       
      ENDIF                                                             
C---COMPUTE AN EVENT WEIGHT                                             
      IF (EVWGT.NE.ZERO) THEN                                           
C---MAMBO phase space generator                                         
        IF (.NOT.HVCONT(1)) THEN                                        
          MINIT = 0                                                     
 100            CALL BVPMAM(HVNTOT,EMSCA,PHEP(1,STARTP),WT)                   
C---Make sure we have an unweighted ps distribution                     
          MINIT = 1                                                     
          IF (.NOT.HWRLOG(WT)) GOTO 100                                 
        ENDIF                                                           
      ENDIF                                                             
 999   RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVHMWG                                                 
C                                                                       
C Generate event data for previously generated phase space configuration
C                                                                       
C Taken from HVHGEN by M Gibbs 21/07/93                                 
C                                                                       
C Error codes:100 BVI event generation flagged                          
C                                                                       
C 21/07/93 Set up for multi-W process colour connections                
C                                                                       
      INCLUDE 'HERWIG65.INC'                                            
      INCLUDE 'HERBVI.INC'                                                  
      INTEGER I,P1,P2,P3                                           
      LOGICAL HWRLOG                                                    
      INTEGER MINIT                                                     
      LOGICAL PRTEN                                                     
      COMMON /BVIMAM/ MINIT,PRTEN                                       
      DOUBLE PRECISION WT                                               
C---Error if BVI process                                                
      IF (HVCONT(2)) CALL HWWARN('HVHMWG',100,*999)                     
C---GENERATE EVENT                                                      
C        PRINT *,' Number w+ ',HVPNWP,CONFIG(13)                        
C        PRINT *,' Number w- ',HVPNWM,CONFIG(14)                        
C---MAMBO phase space generator                                         
      MINIT = 0                                                         
      IF (HVCONT(1)) THEN                                               
 90          CALL BVPMAM(HVNTOT,EMSCA,PHEP(1,STARTP),WT)                     
C---Make sure we have an unweighted ps distribution                     
        MINIT = 1                                                       
        IF (.NOT.HWRLOG(WT)) GOTO 90                                    
      ENDIF                                                             
C---Select interacting quark types                                      
      DO 1 I=1,2                                                        
        IDN(I)=BPTYPE(I)                                                
 1       CONTINUE                                                          
C---Set id for CofM subprocess                                          
      IDCMF=15                                                          
      NHEP=3                                                            
C---Set up 2->2 hard subprocess without phase space using IDN(4)        
      CALL HVETWO                                                       
C---Locations of the 2 particles in hard subprocess                     
      P1=NHEP                                                           
      P2=NHEP-1                                                         
C---CofM location for hard subprocess                                   
      P3=NHEP-2                                                         
      STARTP=P2                                                         
C---Set counter for number of particles                                 
      NHEP=STARTP+HVNTOT-1                                              
C---Set up Mother/Daughter structure                                    
      JDAHEP(1,P3)=STARTP                                               
      JDAHEP(2,P3)=NHEP                                                 
C---Set status codes for HW                                             
      DO 5 I=STARTP+2,NHEP                                              
        ISTHEP(I)=114                                                   
        JMOHEP(1,I)=P3                                                  
        JMOHEP(2,I)=I                                                   
 5         JDAHEP(2,I)=I                                                   
      ISTHEP(STARTP)=113                                                
      ISTHEP(STARTP+1)=114                                              
C---Rotate particles into the Lab frame                                 
      DO 10 I=STARTP,NHEP                                               
        CALL HWULOB(PHEP(1,P3),PHEP(1,I),PHEP(1,I))                     
 10      CONTINUE                                                          
C---Initialise variables for boson decays                               
      DO 12 I=STARTP+2,NHEP                                             
        RHOHEP(1,I)=ONE                                                 
        RHOHEP(2,I)=ONE                                                 
        RHOHEP(3,I)=ONE                                                 
 12      CONTINUE                                                          
C---Set selective decays of bosons (default=any)                        
      DO 15 I=1,MODMAX                                                  
        MODBOS(I)=0                                                     
 15      CONTINUE                                                          
C      DO 13 I=1,HVPNWP+HVPNWM                                          
C        MODBOS(I)=1                                                    
C   13 CONTINUE                                                         
C      DO 14 I=1,HVPNZO                                                 
C        MODBOS(I+HVPNWP+HVPNWM)=3                                      
C   14 CONTINUE                                                         
 999       END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVHPQM                                                 
C                                                                       
C Routine to process all quarks in list by CKM matrix projection        
C                                                                       
C Written by M Gibbs 18/06/93                                           
C                                                                       
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER I,J,STORE(3)                                              
      DOUBLE PRECISION HWRGEN,R                                         
C---Loop over all quarks                                                
      DO 10 I=1,3                                                       
C---Get quark type and number                                           
        J=(2*I)-1                                                       
        STORE(I)=CONFIG(J)                                              
        CONFIG(J)=0                                                     
 10      CONTINUE                                                          
      DO 20 I=1,3                                                       
        IF (STORE(I).LT.1) GOTO 20                                      
C---Project each quark type                                             
        DO 30 J=1,STORE(I)                                              
          R=HWRGEN(I)                                                   
C          PRINT *,R,'for type',I,'iter no. ',J                         
          IF(R.LT.HVCKMC(I,1)) THEN                                     
            CONFIG(1)=CONFIG(1)+1                                       
          ELSEIF(R.LT.HVCKMC(I,2)) THEN                                 
            CONFIG(3)=CONFIG(3)+1                                       
          ELSE                                                          
            CONFIG(5)=CONFIG(5)+1                                       
          ENDIF                                                         
 30          CONTINUE                                                        
 20           CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVINIT                                                 
C                                                                       
C Routine to initialise BVI variables                                   
C                                                                       
C Uses physical parameters from HW - should only be called after        
C HW has set up the common blocks in HWIGIN                             
C                                                                       
C Written by M Gibbs 12/04/93                                           
C                                                                       
C 04/08/93 Flag for HVSGEN generation added                             
C 02/04/95 Total parton-level cs added                                  
C                                                                       
      INCLUDE 'HERWIG65.INC'                                            
      INCLUDE 'HERBVI.INC'                                                   
      DATA HVINID /.FALSE./                                             
      DATA HVCKM/.951,.049,0.,.049,.949,.002,0.,.002,.998/              
C---Prevent re-execution in case variables have been changed            
      IF (HVINID) GOTO 10                                               
C---W boson mass/GeV, from HW                                           
        WMASS = RMASS(198)                                              
C---Weinberg angle for Z0/gamma probabilities                           
        CABB = SWEIN                                                    
C---Mass cut-off for saddle-point calculations                          
        MCOFF = 1.D-5                                                   
C---Parton-level cross section (nb)                                     
        HVPLCS=1.0D0                                                    
C---Set BVPMML to not reject on MAMBO weights                           
        REJCON = .FALSE.                                                
C---Generate Isotropic events                                           
        HVCONT(1) = .TRUE.                                              
C---BVI events to be generated                                          
        HVCONT(2) = .TRUE.                                              
C---Allow CKM matrix projection                                         
        HVCONT(3) = .TRUE.                                              
C---Set HVSGEN to generate events uniformly                             
        HVCONT(4) = .TRUE.                                              
C---Set HVSGEN to not enforce a cut off                                 
        HVCONT(5) = .TRUE.                                              
C---Set HVSGEN to not call user-supplied \shat generation               
        HVCONT(6) = .TRUE.                                              
C---Set HVRBOS to use 1/alpha_W estimate for boson no. dist.            
        HVCONT(7) = .TRUE.                                              
C---Set HVRBOS to not call user-supplied boson no. generation           
        HVCONT(8) = .TRUE.                                              
C---Set HVRBOS to include Higgs boson number in boson distributions     
        HVCONT(9) = .TRUE.                                              
C---SLMIN and SLMAX values for SHAT generation                          
        SLMIN = LOG(0.24D0)                                             
        SLMAX = LOG(0.80D0)                                             
C---Flag to indicate routine has been called                            
        HVINID = .TRUE.                                                 
 10      CONTINUE                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVRBOS(CHDIFF,CHCONT)                                  
C                                                                       
C Generate Boson number distribution for BVI/multi-W process            
C                                                                       
C Returns values as HVPNWP,HVPNWM,HVPNZO,HVPNGA,HVPNHI                  
C                                                                       
C Written by M Gibbs 18/10/93                                           
C                                                                       
C 22/03/94 CHDIFF boson charge difference added                         
C          Higgs number generation forced to be even                    
C          Dependence on FLAG variable removed                          
C 23/03/94 CHCONT control flag added                                    
C          (.TRUE. for number generation)                               
C                                                                       
      INCLUDE 'HERWIG65.INC'                                            
      INCLUDE 'HERBVI.INC'                                                   
      INTEGER IBOS,NUMWO,HVRBIN,CHDIFF                           
      DOUBLE PRECISION CP,CA,CT,CSD                                     
      LOGICAL HWRLOG,CHCONT                                             
C---Test to see if user routine wanted                                  
      IF (HVCONT(8)) THEN                                               
C---Use default - test for number generation or modification            
        IF (CHCONT) THEN                                                
C---1/alpha weak est. or lome?                                          
          IF (HVCONT(7)) THEN                                           
C---Use 1/alpha_w approx for bosons                                     
            IBOS=30                                                     
          ELSE                                                          
C---Estimate the boson number distribution                              
            CT = 3                                                      
C---CT=3 forces use of partic. distribution (FLAG removal M8/53)        
C            CT = HVUSPT(FLAG(3),FLAG(4))                               
            CA = HVGAUA(1,CT)+(EMSCA*HVGAUA(2,CT))                      
            CP = HVGAUP(1,CT)+(EMSCA*HVGAUP(2,CT))                      
C---Rough estimate for number of bosons                                 
            CSD = SQRT(CA*0.5D0)                                        
            CALL BVPNGA(CP,CSD,IBOS)                                    
          ENDIF                                                         
C---Initialise number of higgs particles                                
          HVPNHI=0                                                      
C---Test for Higs number generation                                     
          IF (HVCONT(9)) THEN                                           
C---Higgs boson number estimate                                         
 19                 IF (HWRLOG(0.0625D0)) THEN                                  
C---The number of Higgs particles has to be even!                       
              HVPNHI=HVPNHI+2                                           
              IBOS=IBOS-4                                               
              GOTO 19                                                   
            ENDIF                                                       
          ENDIF                                                         
C---Boson numbers by type                                               
          HVPNWP = IBOS/3                                               
          HVPNWM = HVPNWP                                               
          NUMWO = IBOS - HVPNWP - HVPNWM                                
C---Next assignment allows total number of bosons to be calculated      
          HVPNZO = NUMWO                                                
        ELSE                                                            
C---Now account for charge difference                                   
 39             IF (CHDIFF.NE.0) THEN                                         
            IF (CHDIFF.GT.0) THEN                                       
C---Add extra positive charge                                           
              CHDIFF = CHDIFF - 1                                       
C---Decide how to insert it                                             
              IF (HWRLOG(HALF)) THEN                                    
C---Add extra W+                                                        
                HVPNWP = HVPNWP + 1                                     
                NUMWO = NUMWO - 1                                       
              ELSE                                                      
C---Remove W-                                                           
                HVPNWM = HVPNWM - 1                                     
                NUMWO = NUMWO + 1                                       
              ENDIF                                                     
            ELSE                                                        
C---Add extra negative charge                                           
              CHDIFF = CHDIFF + 1                                       
C---Decide how to insert it                                             
              IF (HWRLOG(HALF)) THEN                                    
C---Add extra W-                                                        
                HVPNWM = HVPNWM + 1                                     
                NUMWO = NUMWO - 1                                       
              ELSE                                                      
C---Remove W+                                                           
                HVPNWP = HVPNWP - 1                                     
                NUMWO = NUMWO + 1                                       
              ENDIF                                                     
            ENDIF                                                       
C---Loop round for all of charge imbalance                              
            GOTO 39                                                     
          ENDIF                                                         
C---Now the number of Gamma particles to use                            
          HVPNGA=HVRBIN(NUMWO,CABB)                                     
          HVPNZO=NUMWO-HVPNGA                                           
        ENDIF                                                           
      ELSE                                                              
C---User supplied routine                                               
        CALL HURBOS(CHDIFF)                                      
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION HVRBIN(TOTAL,PROB)                                       
C                                                                       
C Function to pick binomial distribution using the method               
C of uniform deviates.                                                  
C                                                                       
C TOTAL   Power of distribution                                         
C PROB    p value for distribution                                      
C HVRBIN  Integer between 0 and TOTAL picked according to binomial      
C                                                                       
C Written by M Gibbs 24/04/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      DOUBLE PRECISION BVPALF,PROB,M,GM,T1,                             
     &       TEMP,U,V,X,HWRGEN,GV                                       
      INTEGER TOTAL,HVRBIN                                              
C---Calculate limits                                                    
      T1=DBLE(TOTAL)                                                    
      TEMP=T1*PROB                                                      
      GM=BVPALF(T1+1.)-BVPALF(TEMP+1.)                                  
     &  -BVPALF(T1+1.-TEMP)+(TEMP*LOG(PROB))                            
     &  +((T1-TEMP)*LOG(1.-PROB))                                       
      GM=EXP(GM*0.5D0)                                                  
      M=T1*0.5D0/GM                                                     
C---Generate U and V values                                             
 10    U=HWRGEN(0)*GM                                                    
      V=HWRGEN(1)*M                                                     
C---Calculate X and the comparison                                      
      X=V/U                                                             
      IF (X.GT.T1) GOTO 10                                              
      HVRBIN=INT(X)                                                     
      TEMP=DBLE(HVRBIN)                                                 
      GV=BVPALF(T1+1.)-BVPALF(TEMP+1.)                                  
     &  -BVPALF(T1+1.-TEMP)+(TEMP*LOG(PROB))                            
     &  +((T1-TEMP)*LOG(1.-PROB))                                       
      IF (EXP((GV*0.5D0)).LT.U) GOTO 10                                 
C---This value is accepted, so exit                                     
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVSFUN(WT,Q1,Q2)                                       
C                                                                       
C Routine to generate Sfn's for bvi routine, given SHAT and S           
C                                                                       
C Uses HERWIG57 for sfunctions                                          
C                                                                       
C Written by M Gibbs 01/04/93                                           
C                                                                       
C 12/04/93 Now takes S,SHAT varaibles from common blocks                
C                                                                       
      INCLUDE 'HERWIG65.INC'                                            
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION X1,X2,WT                                         
      INTEGER Q1,Q2                                                     
C---Set up X limits for HWSGEN                                          
      XXMIN = (EMSCA/COFMEN)**2                                         
      XLMIN = LOG(XXMIN)                                                
      EMSCA=EMSCA                                                       
C---Calculate sfunc. contributions                                      
      CALL HWSGEN(.TRUE.)                                               
C---Calculate x-values and the weight for this step                     
      X1 = DISF(Q1,1)                                                   
      X2 = DISF(Q2,2)                                                   
      WT = -XLMIN * X1 * X2                                             
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HVSGEN(WT)                                             
C                                                                       
C Routine to generate shat distribution for MC of pp cs                 
C                                                                       
C Returns SHAT and MC weight for this distribution                      
C                                                                       
C Routines used:HWWARN Herwig error handler                             
C                                                                       
C Error codes: 100 EMSCA greater than total interact energy             
C              101 Attempt to enforce cut off                           
C                                                                       
C Written by M Gibbs 01/04/93                                           
C                                                                       
C 12/04/93 Energy variables now passed by common block                  
C 04/08/93 HVCONT flag added to control \shat generationm               
C          BVPPCS usage removed                                         
C 15/10/93 HVCONT(6) now calls user-supplied routine                    
C                                                                       
      INCLUDE 'HERWIG65.INC'                                            
      INCLUDE 'HERBVI.INC'                                                   
      DOUBLE PRECISION WT,HWRUNI,HWRGEN,SUPPER,SLOWER,SCOMB                              
      COMMON /HVCSEV/ SUPPER,SLOWER,SCOMB                               
C---Only use this routine if not prevented                              
      IF (HVCONT(6)) THEN                                               
C---Generate \shat                                                      
        IF (HVCONT(4)) THEN                                             
C---Generate shat according to 1/shat                                   
          EMSCA = EXP(HWRUNI(0,SLMAX,SLMIN)*0.5D0)*COFMEN               
          WT = SLMAX - SLMIN                                            
        ELSE                                                            
C---Alternative shat distribution generation                            
          EMSCA = COFMEN/((1./SUPPER)+(SCOMB*HWRGEN(0))**2)             
          WT = LOG(SCOMB)                                               
        ENDIF                                                           
C---Trap energy nonconservation error                                   
        IF (EMSCA.GT.COFMEN) THEN                                       
          CALL HWWARN('HVSGEN',100,*999)                                
        ENDIF                                                           
      ELSE                                                              
C---User supplied SF generation                                         
        CALL HUSGEN(WT)                                                 
      ENDIF                                                             
C---Enforce cut off if required                                         
      IF (.NOT.HVCONT(5)) THEN                                          
        CALL HWWARN('HVSGEN',101,*999)                                  
      ENDIF                                                             
 999   RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION HVUSPT(FLAG1,FLAG2)                                      
C                                                                       
C Function to calculate array pointer given logical flags               
C                                                                       
C Written by M Gibbs 24/05/93                                           
C                                                                       
      IMPLICIT NONE                                                     
      INTEGER HVUSPT                                                    
      LOGICAL FLAG1,FLAG2                                               
      IF (FLAG1) THEN                                                   
        HVUSPT = 2                                                      
      ELSE                                                              
        HVUSPT = 1                                                      
      ENDIF                                                             
      IF (FLAG2) HVUSPT = HVUSPT + 2                                    
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      FUNCTION HWUPRA(P)                                                
C                                                                       
C PSEUDO-RAPIDITY (-LOG TAN THETA/2 - SET TO +/-20 IF TOO LARGE)        
C                                                                       
      DOUBLE PRECISION HWUPRA,P(3),PTSQ,PPSQ,CUT,PRAP                   
      PARAMETER (CUT=4.25D-18)                                          
      PTSQ=P(1)**2+P(2)**2                                              
      PPSQ=(SQRT(PTSQ+P(3)**2)+ABS(P(3)))**2                            
      IF (PTSQ.LE.CUT*PPSQ) THEN                                        
        PRAP=20.                                                        
      ELSE                                                              
        PRAP=0.5*LOG(PPSQ/PTSQ)                                         
      ENDIF                                                             
      HWUPRA=SIGN(PRAP,P(3))                                            
      END                                     

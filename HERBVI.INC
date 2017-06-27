C This is the common include file for the saddle and herbvi routines
C
c   combined by Peter Richardson 14/6/99                 
C                                                                       
C INCLUDE file for HVDEVL BVI routines                                  
C                                                                       
C Common Blocks:                                                        
C HVADAT Final state particle data                                      
C HVCEEN Event Energy variables                                         
C HVCENP Event particle numbers                                         
C HVCPDA Physical constants                                             
C HVCRDA Run-time data for analysis                                     
C                                                                       
C HVADAT: Final state particle data                                     
C  NPMAX  Max. number of final state particles in analysis              
C  RAPFSP Rapidity                                                      
C  PSRFSP Pseudo-rapidity                                               
C  PTFSP  Transverse momentum                                           
C  IDHFSP Particle IDHEP                                                
C  IDHNUM Particle number in HEPEVT common block                        
C  ISLFZ  Flag=true for lepton from Z particle                          
C  ISFSP  Flag=true for BVI particle                                    
C  ISACC  Flag=true for accepted particle                               
C  ISCIS  Flag=true for isolation criterion met                         
C  CALPOS Particle position in calorimeter                              
C  IDZPAR  Number in HEPEVT of Z0 if parent                             
C                                                                       
C HVCEEN: Event energy variables                                        
C  COFMEN Total CofM energy of interacting particles                    
C  SLMIN  Minimum value of s for shat generation                        
C  SLMAX  Maximum value of s for shat generation                        
C  HVLMAS Total mass of light(fermion) particles                        
C  HVHMAS Total mass of vector bosons                                   
C  HVGAUA Parametrisation for boson gaussian                            
C  HVGAUP Parametrisation for boson gaussian                            
C                                                                       
C HVCENP: Event particle numbers                                        
C  HVPNWM Number of W- bosons                                           
C  HVPNWP Number of W+ bosons                                           
C  HVPNZO Number of Z0 bosons                                           
C  HVPNGA Number of gamma particles                                     
C  HVPNHI Number of Higgs bosons                                        
C  STARTP Location of first particle written to event table             
C  BPTYPE Type of interacting quarks in beam                            
C  HPNTOT Total number of particles in BVI subprocess                   
C                                                                       
C HVCLFG: Flags for program control                                     
C  HVINID Initialisation performed                                      
C  HVCONT Control array for event generation                            
C         (1) .true. for Isotropic event generation                     
C         (2) .true. for BVI event generation                           
C         (3) .true. for CKM matrix projection                          
C         (4) .true. for uniform event generation                       
C         (5) .true. for non-enforced cut offs                          
C         (6) .true. for no user supplied \shat generation              
C         (7) .true. for generation of n_B as 1/alpha_w estimate        
C         (8) .true. for no user supplied n_B generation                
C         (9) .true. for including non-zero higgs numbers               
C                                                                       
C HVCPDA: Physical constants                                            
C  HVCKM  CKM matrix data                                               
C  HVCKMC Cumulative prob. matrix derived from HVCKM                    
C  HVPLCS Parton-level cross section                                    
C                                                                       
C HVCRDA: Run-time data for analysis                                    
C  HISSUM Histogram sum of entries                                      
C  HISSMS Histogram sum of entries**2                                   
C  HISSME Histogram number of entries                                   
C  HISSMF Histogram number of entries per histogram call                
C  HISSEN Number of histogram for storage                               
C                                                                       
      DOUBLE PRECISION HVLMAS,HVHMAS,COFMEN,                            
     &      HVGAUA(2,4),HVGAUP(2,4),SLMIN,SLMAX                         
      INTEGER HVPNWP,HVPNWM,HVPNZO,STARTP,BPTYPE(4),HVNTOT,             
     &      HVPNGA,HVPNHI                                               
      COMMON /HVCEEN/ COFMEN,                                           
     &      SLMIN,SLMAX,HVLMAS,HVHMAS,                                  
     &      HVGAUA,HVGAUP                                               
      COMMON /HVCENP/ HVPNWP,HVPNWM,HVPNZO,HVPNGA,HVPNHI,               
     &      STARTP,BPTYPE,                                              
     &      HVNTOT                                                      
      INTEGER NPMAX,FSPPTR                                              
      PARAMETER (NPMAX=10000)                                           
      INTEGER IDHFSP(NPMAX),IDHNUM(NPMAX),                              
     &      CALPOS(2,NPMAX),IDZPAR(NPMAX)                               
      LOGICAL ISFSP(NPMAX),ISACC(NPMAX),ISLFZ(NPMAX),ISCIS(NPMAX)       
      DOUBLE PRECISION RAPFSP(NPMAX),PTFSP(NPMAX),PSRFSP(NPMAX)         
      COMMON /HVADAT/ RAPFSP,PTFSP,PSRFSP,                              
     &      IDHFSP,IDHNUM,ISLFZ,ISFSP,ISACC,ISCIS,CALPOS,               
     &      FSPPTR,IDZPAR                                               
      LOGICAL HVINID,HVCONT(10)                                         
      COMMON /HVCLFG/ HVINID,HVCONT                                     
      DOUBLE PRECISION HVCKM(3,3),HVCKMC(3,3),HVPLCS                    
      COMMON /HVCPDA/ HVCKM,HVCKMC,HVPLCS                               
      DOUBLE PRECISION HISSUM(10),HISSMS(10)                            
      INTEGER HISSME(10),HISSMF(10),HISSEN(3,10)                        
      COMMON /HVCRDA/ HISSUM,HISSMS,HISSME,HISSMF,HISSEN               
C                                                                       
C 17/10/92 INCLUDE file SADDLE.INC for saddle point routines            
C 16/02/93 Boson MC internal variable common block added                
C 07/04/93 New parameter NCONF added                                    
C 11/04/93 Common block of parameters added                             
C                                                                       
C                                                                       
C Parameters:                                                           
C NTYPE  Number of different particle types                             
C NCONF  Number of allowed particle types for CONFIG(NCONF)             
C                                                                       
C Common blocks:                                                        
C                                                                       
C BVPCFG Common block containing particle configuration                 
C        BVPNUM()  Number of each particle type                         
C        BVPPID()  Particle ID codes                                    
C        BVPTOT    Total number of particles                            
C        BVPMAS()  Particle masses                                      
C                                                                       
C BVPCBD Saddle point routine internal data                             
C        BVPIR(4)       \  Bessel function ratios                       
C        BVPKR(NTYPE,4) /                                               
C                                                                       
C BVPCSD Saddle point data                                              
C        BETA      Saddle point location                                
C        EN        Centre of mass energy                                
C        ACCU      Required NR accuracy                                 
C                                                                       
C BVPCDR Integrand derivatives and final data                           
C        BVPDER(5) Integrand and derivatives                            
C        BVPPSI    Integral                                             
C        BVPPCO    Next-order corrections to integral                   
C                                                                       
C BVPCBS Internal variables for Boson MC routines                       
C        BMIN      Min number of bosons for this energy                 
C        BMAX      Max number of bosons for this energy                 
C        NITS      Number of iterations for boson MC                    
C        NBMAX     Number of bosons for max. contrib                    
C        RES(200)  Results of boson calculations                        
C        TSUM      Sum of boson contributions                           
C        CONFIG    BVI configuration                                    
C        NBOSSQ    Av. value of number of bosons squared                
C                                                                       
C BVPCBP Physical parameters                                            
C        WMASS     Mass of W boson                                      
C        CABB      Cabbibo angle                                        
C        MCOFF     Mass Cut-off for ps routines                         
C        REJCON    .TRUE. to allow rejection in BVPMML                  
C                                                                       
      INTEGER NTYPE,NCONF                                               
      PARAMETER (NTYPE=20)                                              
      PARAMETER (NCONF=17)                                              
      DOUBLE PRECISION BVPMAS(NTYPE),BVPIR(4),BVPKR(NTYPE,4),           
     & EN,BETA,ACCU,BVPDER(5),BVPPSI,BVPPCO                             
      DOUBLE PRECISION TSUM,RES(200),NBOSSQ                             
      DOUBLE PRECISION WMASS,CABB,MCOFF                                 
      INTEGER BVPNUM(NTYPE),BVPTOT,BVPPID(NTYPE)                        
      INTEGER NBMAX,NITS,BMIN,BMAX,CONFIG(NCONF)                        
      LOGICAL REJCON                                                    
C                                                                       
      COMMON /BVPCFG/ BVPMAS,BVPNUM,BVPPID,BVPTOT                       
      COMMON /BVPCBD/ BVPIR,BVPKR                                       
      COMMON /BVPCSD/ BETA,EN,ACCU                                      
      COMMON /BVPCDR/ BVPDER,BVPPSI,BVPPCO                              
      COMMON /BVPCBS/ TSUM,RES,NBOSSQ,NBMAX,NITS,BMIN,BMAX,CONFIG       
      COMMON /BVPCBP/ WMASS,CABB,MCOFF,REJCON                           
                                                                        
! A NEW MONTE CARLO TREATMENT OF MULTIPARTICLE PHASE SPACE AT     
! HIGH ENERGIES.  R. KLEISS, W.J. STIRLING, S.D. ELLIS.            
! REF. IN COMP. PHYS. COMMUN. 40 (1986) 359                            
      SUBROUTINE RAMBO(N,ET,XM,P,WT,LW)                               
Cf2py intent(in)  n
Cf2py intent(in)  et
Cf2py intent(in)  xm
Cf2py intent(in)  lw
Cf2py intent(out) p
Cf2py intent(out) wt
C------------------------------------------------------               
C                                                                     
C                       RAMBO                                         
C                                                                     
C             RA(NDOM)  M(OMENTA)  BO(OSTER)                          
C                                                                     
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR                
C    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING                 
C                                                                     
C    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)              
C    ET = TOTAL CENTRE-OF-MASS ENERGY                                 
C    XM = PARTICLE MASSES ( DIM=N )                                   
C    P  = PARTICLE MOMENTA ( DIM=(4,N) )                              
C    WT = WEIGHT OF THE EVENT                                         
C    LW = FLAG FOR EVENT WEIGHTING:                                   
C         LW = 0 WEIGHTED EVENTS                                      
C         LW = 1 UNWEIGHTED EVENTS ( FLAT PHASE SPACE )               
C------------------------------------------------------               
      IMPLICIT REAL*8(A-H,O-Z)                                        
      DIMENSION XM(N),P(4,N),Q(4,100),Z(100),R(4),                    
     .   B(3),P2(100),XM2(100),E(100),V(100),IWARN(5)                 
      DATA ACC/1.D-14/,ITMAX/6/,IBEGIN/0/,IWARN/5*0/                  
C                                                                     
C INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT          
      IF(IBEGIN.NE.0) GOTO 103                                        
      IBEGIN=1                                                        
      TWOPI=8.*DATAN(1.D0)                                            
      PO2LOG=DLOG(TWOPI/4.)                                           
      Z(2)=PO2LOG                                                     
      DO 101 K=3,100                                                  
  101 Z(K)=Z(K-1)+PO2LOG-2.*DLOG(DFLOAT(K-2))                         
      DO 102 K=3,100                                                  
  102 Z(K)=(Z(K)-DLOG(DFLOAT(K-1)))                                   
C                                                                     
C CHECK ON THE NUMBER OF PARTICLES                                    
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104                                
      PRINT 1001,N                                                    
      STOP                                                            
C                                                                     
C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES      
  104 XMT=0.                                                          
      NM=0                                                            
      DO 105 I=1,N                                                    
      IF(XM(I).NE.0.D0) NM=NM+1                                       
  105 XMT=XMT+DABS(XM(I))                                             
      IF(XMT.LE.ET) GOTO 106                                          
      PRINT 1002,XMT,ET                                               
      STOP                                                            
C                                                                     
C CHECK ON THE WEIGHTING OPTION                                       
  106 IF(LW.EQ.1.OR.LW.EQ.0) GOTO 201                                 
      PRINT 1003,LW                                                   
      STOP                                                            
C                                                                     
C THE PARAMETER VALUES ARE NOW ACCEPTED                               
C                                                                     
C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE                 
  201 DO 202 I=1,N                                                    
      C=2.*RN(1)-1.                                                   
      S=DSQRT(1.-C*C)                                                 
      F=TWOPI*RN(2)                                                   
      Q(4,I)=-DLOG(RN(3)*RN(4))                                       
      Q(3,I)=Q(4,I)*C                                                 
      Q(2,I)=Q(4,I)*S*DCOS(F)                                         
  202 Q(1,I)=Q(4,I)*S*DSIN(F)                                         
C                                                                     
C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION            
      DO 203 I=1,4                                                    
  203 R(I)=0.                                                         
      DO 204 I=1,N                                                    
      DO 204 K=1,4                                                    
  204 R(K)=R(K)+Q(K,I)                                                
      RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)                     
      DO 205 K=1,3                                                    
  205 B(K)=-R(K)/RMAS                                                 
      G=R(4)/RMAS                                                     
      A=1./(1.+G)                                                     
      X=ET/RMAS                                                       
C                                                                     
C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S                          
      DO 207 I=1,N                                                    
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)                          
      DO 206 K=1,3                                                    
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))                            
  207 P(4,I)=X*(G*Q(4,I)+BQ)                                          
C                                                                     
C RETURN FOR UNWEIGHTED MASSLESS MOMENTA                              
      WT=1.D0                                                         
      IF(NM.EQ.0.AND.LW.EQ.1) RETURN                                  
C                                                                     
C CALCULATE WEIGHT AND POSSIBLE WARNINGS                              
      WT=PO2LOG                                                       
      IF(N.NE.2) WT=(2.*N-4.)*DLOG(ET)+Z(N)                           
      IF(WT.GE.-180.D0) GOTO 208                                      
      IF(IWARN(1).LE.5) PRINT 1004,WT                                 
      IWARN(1)=IWARN(1)+1                                             
  208 IF(WT.LE. 174.D0) GOTO 209                                      
      IF(IWARN(2).LE.5) PRINT 1005,WT                                 
      IWARN(2)=IWARN(2)+1                                             
C                                                                     
C RETURN FOR WEIGHTED MASSLESS MOMENTA                                
  209 IF(NM.NE.0) GOTO 210                                            
      WT=DEXP(WT)                                                     
      RETURN                                                          
C                                                                     
C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X                
  210 XMAX=DSQRT(1.-(XMT/ET)**2)                                      
      DO 301 I=1,N                                                    
      XM2(I)=XM(I)**2                                                 
  301 P2(I)=P(4,I)**2                                                 
      ITER=0                                                          
      X=XMAX                                                          
      ACCU=ET*ACC                                                     
  302 F0=-ET                                                          
      G0=0.                                                           
      X2=X*X                                                          
      DO 303 I=1,N                                                    
      E(I)=DSQRT(XM2(I)+X2*P2(I))                                     
      F0=F0+E(I)                                                      
  303 G0=G0+P2(I)/E(I)                                                
      IF(DABS(F0).LE.ACCU) GOTO 305                                   
      ITER=ITER+1                                                     
      IF(ITER.LE.ITMAX) GOTO 304                                      
      PRINT 1006,ITMAX                                                
      GOTO 305                                                        
  304 X=X-F0/(X*G0)                                                   
      GOTO 302                                                        
  305 DO 307 I=1,N                                                    
      V(I)=X*P(4,I)                                                   
      DO 306 K=1,3                                                    
  306 P(K,I)=X*P(K,I)                                                 
  307 P(4,I)=E(I)                                                     
C                                                                     
C CALCULATE THE MASS-EFFECT WEIGHT FACTOR                             
      WT2=1.                                                          
      WT3=0.                                                          
      DO 308 I=1,N                                                    
      WT2=WT2*V(I)/E(I)                                               
  308 WT3=WT3+V(I)**2/E(I)                                            
      WTM=(2.*N-3.)*DLOG(X)+DLOG(WT2/WT3*ET)                          
      IF(LW.EQ.1) GOTO 401                                            
C                                                                     
C RETURN FOR  WEIGHTED MASSIVE MOMENTA                                
      WT=WT+WTM                                                       
      IF(WT.GE.-180.D0) GOTO 309                                      
      IF(IWARN(3).LE.5) PRINT 1004,WT                                 
      IWARN(3)=IWARN(3)+1                                             
  309 IF(WT.LE. 174.D0) GOTO 310                                      
      IF(IWARN(4).LE.5) PRINT 1005,WT                                 
      IWARN(4)=IWARN(4)+1                                             
  310 WT=DEXP(WT)                                                     
      RETURN                                                          
C                                                                     
C UNWEIGHTED MASSIVE MOMENTA REQUIRED: ESTIMATE MAXIMUM WEIGHT        
  401 WT=DEXP(WTM)                                                    
      IF(NM.GT.1) GOTO 402                                            
C                                                                     
C ONE MASSIVE PARTICLE                                                
      WTMAX=XMAX**(4*N-6)                                             
      GOTO 405                                                        
  402 IF(NM.GT.2) GOTO 404                                            
C                                                                     
C TWO MASSIVE PARTICLES                                               
      SM2=0.                                                          
      PM2=0.                                                          
      DO 403 I=1,N                                                    
      IF(XM(I).EQ.0.D0) GOTO 403                                      
      SM2=SM2+XM2(I)                                                  
      PM2=PM2*XM2(I)                                                  
  403 CONTINUE                                                        
      WTMAX=((1.-SM2/(ET**2))**2-4.*PM2/ET**4)**(N-1.5)               
      GOTO 405                                                        
C                                                                     
C MORE THAN TWO MASSIVE PARTICLES: AN ESTIMATE ONLY                   
  404 WTMAX=XMAX**(2*N-5+NM)                                          
C                                                                     
C DETERMINE WHETHER OR NOT TO ACCEPT THIS EVENT                       
  405 W=WT/WTMAX                                                      
      IF(W.LE.1.D0) GOTO 406                                          
      IF(IWARN(5).LE.5) PRINT 1007,WTMAX,W                            
      IWARN(5)=IWARN(5)+1                                             
  406 CONTINUE                                                        
      IF(W.LT.RN(5)) GOTO 201                                         
      WT=1.D0                                                         
      RETURN                                                          
 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')   
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',            
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)                          
 1003 FORMAT(' RAMBO FAILS: LW=',I3,' IS NOT AN ALLOWED OPTION')      
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')   
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')   
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',       
     . ' DESIRED ACCURACY =',D15.6)                                   
 1007 FORMAT(' RAMBO WARNS: ESTIMATE FOR MAXIMUM WEIGHT =',D15.6,     
     . '     EXCEEDED BY A FACTOR ',D15.6)                            
      END                                                             

      function RN(IDMY)
      implicit NONE
      integer         IDMY
      real*8 RN            
      RN=rand()
      return
      END


                                                                      

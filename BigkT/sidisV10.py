#!/bin/env python
import sys,os
import copy
import numpy as np
import pylab as py
from scipy.integrate import quad,fixed_quad
from qcdlib import aux
from qcdlib import alphaS
from tools.config import conf
from multiprocessing import Process,Pipe
import warnings
warnings.filterwarnings('ignore') 

from LO import LO
from NLO import Pg,Ppp

import NLO.Pg.fchn1A,NLO.Ppp.fchn1A 
import NLO.Pg.fchn2A,NLO.Ppp.fchn2A
import NLO.Pg.fchn3A,NLO.Ppp.fchn3A
import NLO.Pg.fchn4A,NLO.Ppp.fchn4A
import NLO.Pg.fchn5A,NLO.Ppp.fchn5A
import NLO.Pg.fchn6A,NLO.Ppp.fchn6A

import NLO.Pg.fchn1B,NLO.Ppp.fchn1B 
import NLO.Pg.fchn2B,NLO.Ppp.fchn2B
import NLO.Pg.fchn3B,NLO.Ppp.fchn3B
import NLO.Pg.fchn4B,NLO.Ppp.fchn4B
import NLO.Pg.fchn5B,NLO.Ppp.fchn5B
import NLO.Pg.fchn6B,NLO.Ppp.fchn6B

import NLO.Pg.fchn1C,NLO.Ppp.fchn1C 
import NLO.Pg.fchn2C,NLO.Ppp.fchn2C
import NLO.Pg.fchn3C,NLO.Ppp.fchn3C
import NLO.Pg.fchn4C,NLO.Ppp.fchn4C
import NLO.Pg.fchn5C,NLO.Ppp.fchn5C
import NLO.Pg.fchn6C,NLO.Ppp.fchn6C

from external.CJLIB.CJ import CJ
from external.DSSLIB.DSS import DSS
import vegas
from numba import jit

conf['order']='NLO'
conf['Q20'] = 1.0
conf['aux']=aux.AUX()
aS=alphaS.ALPHAS()
pdf=CJ({'path2CJ':'external/CJLIB','order':'NLO'})
ff=DSS({'path2DSS':'external/DSSLIB','order':'NLO'})
zero=1e-10
alfa=1/137.036
M=0.93891897
eU= 2./3
eD=-1./3
eq3=np.array([0,eU,-eU,eD,-eD,eD,-eD, 0,  0, 0,  0])
eq4=np.array([0,eU,-eU,eD,-eD,eD,-eD,eU,-eU, 0,  0])
eq5=np.array([0,eU,-eU,eD,-eD,eD,-eD,eU,-eU,eD,-eD])
            # g, u, ub, d, db, s, sb, c, cb, b, bb

qqp=np.ones((10,10))
for i in range(10): qqp[i,i]=0
for i in range(10): qqp[i,0]=0
for i in range(10): qqp[0,i]=0

def get_parton_lum(xi,zeta,mu2,chn,case,f,d,eq):
    # case A: gg
    # case B: g-gp
    # case C: gp-gp
    # g:0,u:1,ub:2,d:3,db:4,s:5,sb:6,c:7,cb:8,b:9,bb:10



    # 1 V: A g -> (q->h) qb  R: A g -> (q->h)  qb g  
    # 2 V: A q -> (q->h) g   R: A q -> (q->h)  g  g   
    #                        R: A q -> (q->h)  q' qb'
    # 3 V: A q -> (g->h) q   R: A q -> (g->h)  q  g  
    # 4                      R: A g -> (g->h)  q  qb  
    # 5                      R: A q -> (qb->h) q  qb 
    # 6                      R: A q -> (q'->h) q  qb  

    if   chn==1:
      if case=='A':
        lum=np.sum(eq[1:]*eq[1:]*f[0]*d[1:])
      else: lum=0

    elif chn==3:
      if case=='A': 
        lum=np.sum(eq[1:]*eq[1:]*f[1:]*d[0])
      else: lum=0

    elif chn==4:
      if case=='A': 
        lum=np.sum(eq[1:]*eq[1:])*f[0]*d[0]
      else: lum=0

    elif chn==5:
      if case =='A':
        idx1=[1,3,5,7,2,4,6,8,10]
        idx2=[2,4,6,8,1,3,5,7,9]
        lum=np.sum(eq[idx1]*eq[idx1]*f[idx1]*d[idx2])
      else: lum=0

    elif chn==2:
      if   case=='A':
        lum=np.sum(eq[1:]*eq[1:]*f[1:]*d[1:])
      elif case=='B':
        lum=np.sum(eq[1:]*f[1:]*d[1:]) * np.sum(eq[1:])
      elif case=='C':
        lum=np.sum(f[1:]*d[1:]) * np.sum(eq[1:]*eq[1:])
    
    elif chn==6:
      if   case=='A':
        lum = np.einsum('i,j,ij',eq[1:]*eq[1:]*f[1:],d[1:],qqp)
      elif   case=='B':
        lum = np.einsum('i,j,ij',eq[1:]*f[1:],eq[1:]*d[1:],qqp)
      elif   case=='C':
        lum = np.einsum('i,j,ij',f[1:],eq[1:]*eq[1:]*d[1:],qqp) 

    return lum

def get_dxsec(xi,s23,x,z,Q,qT,y,E,tar,had,order,part=None):
        
    mu2=Q**2
    xh=x/xi

    zh=((1-xh)-xh*s23/Q**2)/((1-xh)+xh*qT**2/Q**2)
    zh0=((1-xh))/((1-xh)+xh*qT**2/Q**2)
    zeta=z/zh
    zeta0=z/zh0

    jac  =zeta  * xh/Q**2 / ((1-xh)-xh*s23/Q**2)
    jac0 =zeta0 * xh/Q**2 / (1-xh)
    gs2=aS.get_alphaS(mu2)*4*np.pi    

    # get pdfs
    # g:0,u:1,ub:2,d:3,db:4,s:5,sb:6,c:7,cb:8,b:9,bb:10
    f=pdf.get_f(xi,mu2)    

    if tar=='n':
      u,d=f[1],f[3]
      f[1],f[3]=d,u

    if tar=='d':
      u,d=f[1],f[3]
      f[1],f[3]=0.5*(u+d),0.5*(u+d)
    
    # get ffs
    # g:0,u:1,ub:2,d:3,db:4,s:5,sb:6,c:7,cb:8,b:9,bb:10
    d0=ff.get_f(zeta0,mu2,had)
    d =ff.get_f(zeta ,mu2,had)


    Nf=aS.get_Nf(mu2)
    if Nf<5: 
      f[9:11]  = np.zeros(2)
      d[9:11]  = np.zeros(2)
      d0[9:11] = np.zeros(2)
    if Nf<4: 
      f[7:9]  = np.zeros(2)
      d[7:9]  = np.zeros(2)
      d0[7:9] = np.zeros(2)
    if   Nf==3: eq=eq3[:]
    elif Nf==4: eq=eq4[:]
    elif Nf==5: eq=eq5[:]

    Fg,Fpp=0,0

    # 1 V: A g -> (q->h) qb  R: A g -> (q->h)  qb g  
    # 2 V: A q -> (q->h) g   R: A q -> (q->h)  g  g   
    #                        R: A q -> (q->h)  q' qb'
    # 3 V: A q -> (g->h) q   R: A q -> (g->h)  q  g  
    # 4                      R: A g -> (g->h)  q  qb  
    # 5                      R: A q -> (qb->h) q  qb 
    # 6                      R: A q -> (q'->h) q  qb  

    if order==0:

      for chn in [1,2,3]:
        lum0=get_parton_lum(xi,zeta0,mu2,chn,'A',f,d0,eq)
        factor0= gs2 * jac0*lum0/xi/zeta0 * zh0

        if   chn==1: _Pg,_Ppp=LO.PgC,LO.PppC  # C: A g -> (q->h) qb  
        elif chn==2: _Pg,_Ppp=LO.PgA,LO.PppA  # A: A q -> (q->h) g   
        elif chn==3: _Pg,_Ppp=LO.PgB,LO.PppB  # B: A q -> (g->h) q   

        Fg  += _Pg( xh, zh0, qT**2, Q**2)/(2*np.pi)**3*factor0  
        Fpp += _Ppp(xh, zh0, qT**2, Q**2)/(2*np.pi)**3*factor0 

    if order==1:

      s=(1-xh)/xh*Q**2
      t=-(1-zh)*Q**2-zh*qT**2
      t0=-(1-zh0)*Q**2-zh0*qT**2
      nf=4.0
      B=Q**2*(1/xh-1)*(1-z)-z*qT**2

      for chn in [1,2,3,4,6]:

        for case in ['A','B','C']:

          if case=='B' or case=='C':
            if chn==1: continue
            if chn==3: continue
            if chn==4: continue

          lum0=get_parton_lum(xi,zeta0,mu2,chn,case,f,d0,eq)
          lum =get_parton_lum(xi,zeta0,mu2,chn,case,f,d,eq)

          factor = gs2**2 * jac*lum/xi/zeta    * zh
          factor0= gs2**2 * jac0*lum0/xi/zeta0 * zh0

          if   chn==1 and case=='A': _Pg,_Ppp=Pg.fchn1A,Ppp.fchn1A   
          elif chn==2 and case=='A': _Pg,_Ppp=Pg.fchn2A,Ppp.fchn2A   
          elif chn==3 and case=='A': _Pg,_Ppp=Pg.fchn3A,Ppp.fchn3A   
          elif chn==4 and case=='A': _Pg,_Ppp=Pg.fchn4A,Ppp.fchn4A   
          elif chn==5 and case=='A': _Pg,_Ppp=Pg.fchn5A,Ppp.fchn5A   
          elif chn==6 and case=='A': _Pg,_Ppp=Pg.fchn6A,Ppp.fchn6A   

          elif chn==1 and case=='B': _Pg,_Ppp=Pg.fchn1B,Ppp.fchn1B   
          elif chn==2 and case=='B': _Pg,_Ppp=Pg.fchn2B,Ppp.fchn2B   
          elif chn==3 and case=='B': _Pg,_Ppp=Pg.fchn3B,Ppp.fchn3B   
          elif chn==4 and case=='B': _Pg,_Ppp=Pg.fchn4B,Ppp.fchn4B   
          elif chn==5 and case=='B': _Pg,_Ppp=Pg.fchn5B,Ppp.fchn5B   
          elif chn==6 and case=='B': _Pg,_Ppp=Pg.fchn6B,Ppp.fchn6B   

          elif chn==1 and case=='C': _Pg,_Ppp=Pg.fchn1C,Ppp.fchn1C   
          elif chn==2 and case=='C': _Pg,_Ppp=Pg.fchn2C,Ppp.fchn2C   
          elif chn==3 and case=='C': _Pg,_Ppp=Pg.fchn3C,Ppp.fchn3C   
          elif chn==4 and case=='C': _Pg,_Ppp=Pg.fchn4C,Ppp.fchn4C   
          elif chn==5 and case=='C': _Pg,_Ppp=Pg.fchn5C,Ppp.fchn5C   
          elif chn==6 and case=='C': _Pg,_Ppp=Pg.fchn6C,Ppp.fchn6C   

          if part=='regular':
            Fg += _Pg.regular(1,1,s,t,Q,s23,Q,nf)*factor
            Fpp+= _Ppp.regular(1,1,s,t,Q,s23,Q,nf)*factor

          if chn<4:

             if part=='delta':
              Fg += _Pg.delta(1,1,s,t0,Q,zero,Q,B,nf)*factor0  / (B-0)
              Fpp+= _Ppp.delta(1,1,s,t0,Q,zero,Q,B,nf)*factor0 / (B-0)
              
             elif part=='plus':
              Fg +=  (_Pg.plus1B(1,1,s,t,Q,s23,Q,B,nf)*factor  - _Pg.plus1B(1,1,s,t0,Q,zero,Q,B,nf)*factor0)/s23
              Fpp+=  (_Ppp.plus1B(1,1,s,t,Q,s23,Q,B,nf)*factor - _Ppp.plus1B(1,1,s,t0,Q,zero,Q,B,nf)*factor0)/s23

              Fg += (_Pg.plus2B(1,1,s,t,Q,s23,Q,B,nf)*factor  - _Pg.plus2B(1,1,s,t0,Q,zero,Q,B,nf)*factor0)*np.log(s23)/s23
              Fpp+= (_Ppp.plus2B(1,1,s,t,Q,s23,Q,B,nf)*factor - _Ppp.plus2B(1,1,s,t0,Q,zero,Q,B,nf)*factor0)*np.log(s23)/s23
          
      if np.isnan(Fg) or np.isnan(Fpp): return 0  
      if np.isinf(Fg) or np.isinf(Fpp): return 0 

    F1h = -0.5*Fg + 2*xh**2/Q**2*Fpp
    F2h = -xh*Fg   + 12*xh**3/Q**2*Fpp
    #dxsec=np.pi**2*alfa**2/Q**4  * zh*(y**2*F1h+(1-y)*F2h/xh)
    dxsec=np.pi**2*alfa**2/Q**4  * (y**2*F1h+(1-y)*F2h/xh)
    dxsec/=z**2  # from qT^2  -> pT^2
    return dxsec
    #if wate==None:     return dxsec
    #elif wate=='xi':   return xi*dxsec
    #elif wate=='zeta': return zeta*dxsec

def A(x,z,qT,Q):
    return x*(1+z*qT**2/(1-z)/Q**2)

def B(x,xi,z,qT,Q):
    return Q**2*(1/(x/xi)-1)*(1-z)-z*qT**2

def get_dxsec_ubox(u,x,z,Q,qT,y,E,tar,had,order,part=None):
    ximin=A(x,z,qT,Q)
    xi=ximin+u[0]*(1-ximin)

    if order==0:
      return get_dxsec(xi,0,x,z,Q,qT,y,E,tar,had,order,part)*(1-ximin)

    elif order==1:

      if part=='delta':
        return get_dxsec(xi,0,x,z,Q,qT,y,E,tar,had,order,part)*(1-ximin)

      else:
        s23min=0
        s23max=B(x,xi,z,qT,Q)
        s23=s23min+u[1]*(s23max-s23min)
        return get_dxsec(xi,s23,x,z,Q,qT,y,E,tar,had,order,part)*(1-ximin)*(s23max-s23min)

n=40  # for squad and dsquad 

def squad(df_du1,x,z,Q,qT):
    return quad(df_du1,0,1)[0]
    #return fixed_quad(np.vectorize(df_du1),0,1,n=n)[0]

def dquad(df_du1_du2,x,z,Q,qT):
    #df_xi=lambda xi: quad(lambda s23: df_dxi_ds23(xi,s23)[idx],zero,B(x,xi,z,qT,Q))[0]
    #return quad(df_xi,A(x,z,qT,Q),1)[0]
    df_du1=lambda u1: fixed_quad(np.vectorize(lambda u2: df_du1_du2(u1,u2)),0,1,n=n)[0]
    return fixed_quad(np.vectorize(df_du1),0,1,n=n)[0]

def get_xsec(x,z,Q,qT,y,E,tar,had,order,part=None):

  if order==0:
    df_du1=lambda u1: get_dxsec_ubox([u1,0],x,z,Q,qT,y,E,tar,had,0,part)
    return squad(df_du1,x,z,Q,qT)

  if order==1:

    if part=='delta':
      df_du1=lambda u1: get_dxsec_ubox([u1,0],x,z,Q,qT,y,E,tar,had,1,part)
      return squad(df_du1,x,z,Q,qT)

    else:
      df_du1_du2=lambda u1,u2: get_dxsec_ubox([u1,u2],x,z,Q,qT,y,E,tar,had,1,part)
      return dquad(df_du1_du2,x,z,Q,qT)

    #integ = vegas.Integrator([[0, 1], [0, 1]])
    #result = integ(lambda u: get_dxsec_ubox(u,x,z,Q,qT,y,E,tar,had,1), nitn=13, neval=1000)
    #print result.summary()
    #print 'result = %s    Q = %.2f' % (result, result.Q)
    #return result.val

if __name__=="__main__":

  #x=0.017
  #z=0.2437
  #Q=1.87
  #qT=6.86
  #y=0.695

  x= 0.0062
  z= 0.6862
  Q= 1.140175
  qT= 2.290746
  y=0.694

  #x  =0.0752
  #z  =0.4818
  #Q  =2.144761  
  #qT =3.071113  
  #y  =0.213

  E=160
  tar='d'
  had='h+'
  print 'born  = %0.3e'%get_xsec(x,z,Q,qT,y,E,tar,had,0)
  print 'rest  = %0.3e'%get_xsec(x,z,Q,qT,y,E,tar,had,1,'regular')
  print 'delta = %0.3e'%get_xsec(x,z,Q,qT,y,E,tar,had,1,'delta')
  print 'plus  = %0.3e'%get_xsec(x,z,Q,qT,y,E,tar,had,1,'plus')











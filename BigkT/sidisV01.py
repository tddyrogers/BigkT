#!/usr/bin/env python
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
from external.CJLIB.CJ import CJ
from external.DSSLIB.DSS import DSS
from numba import jit

conf['order']='NLO'
conf['Q20'] = 1.0
conf['aux']=aux.AUX()
aS=alphaS.ALPHAS()
pdf=CJ({'path2CJ':'external/CJLIB','order':'NLO'})
ff=DSS({'path2DSS':'external/DSSLIB','order':'NLO'})
zero=1e-10
alfa=1/137.036
CF=4.0/3.0
Mp=0.93891897
eU2=4./9
eD2=1./9
eq2=np.array([0,eU2,eU2,eD2,eD2,eD2,eD2,eU2,eU2,eU2,eU2,eD2,eD2])

def get_A(y):
  A1=1+(2/y-1)**2
  A2=-2
  A3=0
  A4=0
  return A1,A2,A3,A4

def get_sfA_jk(xh,zh,Q,qT,y):
  A1,A2,A3,A4=get_A(y)
  return 2*CF*xh*zh*(\
           (1/qT**2*(Q**4/xh**2/zh**2+(Q**2-qT**2)**2)+6*Q**2)*A1\
          +2*Q**2*(2*A2+A4)+2*Q/qT*(Q**2+qT**2)*A3)/2/Q**2

def get_sfA_jg(xh,zh,Q,qT,y):
  A1,A2,A3,A4=get_A(y)
  return xh*(1-xh)*(\
         ((Q**4/qT**2)*(1/(xh*zh)**2-2/(xh*zh)+2)+2*Q**2*(5-1/xh-1/zh))*A1\
         +4*Q**2*(2*A2+A4)+2*Q/qT*(2*(Q**2+qT**2)-Q**2/xh/zh)*A3)/2/Q**2

def get_sfA_gj(xh,zh,Q,qT,y):
  A1,A2,A3,A4=get_A(y)
  qtT=zh*qT/(1-zh)
  return 2*CF*xh*(1-zh)*((1/qtT**2*(Q**4/xh**2/(1-zh)**2+(Q**2-qtT**2)**2)\
        +6*Q**2)*A1+2*Q**2*(2*A2+A4)\
        +2*Q/qtT*(Q**2+qtT**2)*A3)/2/Q**2

def get_pdf_ff(x,z,mu,tar,had):
  f=pdf.get_f(x,mu**2)  
  if tar=='n':
    u,d = f[1],f[3]
    f[1],f[3] = d,u
  if tar=='d':
    u,d = f[1],f[3]
    f[1],f[3] = 0.5*(u+d),0.5*(u+d)
  d=ff.get_f(z,mu**2,had)
  return f,d

def get_dxsec(xia,x,y,z,Q,qT,E,tar,had):
  xh = x/xia
  zh = ((1-xh))/((1-xh)+xh*qT**2/Q**2)
  xib=z/zh

  #s=2*E*Mp+Mp**2  # this makes a difference
  s=Q**2/(x*y)     # so we use this 
  mu=Q
  e2=alfa*4*np.pi
  sigo=Q**2/(4*np.pi*s*x**2)*(e2/2)
  Fl=e2/2/Q**2
  alphaS=aS.get_alphaS(mu**2)
  
  f,d=get_pdf_ff(xia,xib,mu,tar,had)
  # g,u,ub,d,db,s,sb,c,cb,b,bb
  # 0 1  2 3  4 5  6 7  8 9,10
  Nf=aS.get_Nf(mu**2)
  iNf=np.zeros(11)
  iNf[1:7]=np.ones(6)
  if Nf>3: iNf[7:9]=np.ones(2) 
  if Nf>4: iNf[9:11]=np.ones(2) 

  lumjk=np.sum(iNf[1:9]*eq2[1:9]*f[1:9]*d[1:9])
  lumgj=np.sum(iNf[1:9]*eq2[1:9]*f[1:9]*d[0])
  lumjg=np.sum(iNf[1:9]*eq2[1:9]*f[0]*d[1:9])

  M = lumjk * get_sfA_jk(xh,zh,Q,qT,y)
  M+= lumjg * get_sfA_jg(xh,zh,Q,qT,y)
  M+= lumgj * get_sfA_gj(xh,zh,Q,qT,y)
  M*= sigo*Fl/(4*np.pi*s*Q**2)*alphaS/np.pi*xh*zh

  out = M
  out*= 1/(xia-x)
  out*= 2*np.pi  # azimuthal integration
  out*= 1/z**2  

  return out

def get_dxsec_ubox(u,x,z,Q,qT,y,E,tar,had):
  w=qT/Q*np.sqrt(x*z)
  xiamin=w**2/(1-z)+x
  xia=xiamin+u*(1-xiamin)
  return get_dxsec(xia,x,y,z,Q,qT,E,tar,had)*(1-xiamin)

def get_xsec(x,z,Q,qT,y,E,tar,had):
  func = lambda u: get_dxsec_ubox(u,x,z,Q,qT,y,E,tar,had)
  return quad(func,0,1)[0]


if __name__== "__main__":


  x= 0.0062
  z= 0.6862
  Q= 1.140175
  qT= 2.290746
  y=0.694

  E=160
  tar='d'
  had='h+'
  print 'LO =%0.3e'%get_xsec(x,z,Q,qT,y,E,tar,had)
  #print get_dxsec(0.5,x,y,z,Q,qT,E,tar,had)

  #print get_dxsec(0.5,x,y,z,Q,qT,E,tar,had)



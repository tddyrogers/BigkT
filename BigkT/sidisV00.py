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
M=0.93891897

# set charges
e2U=4.0/9.0
e2D=1.0/9.0
couplings={1:e2U,2:e2D,3:e2D,4:e2U,5:e2D,6:e2U}


def xhat(x,xia):
  return x/xia

def zhat(z,xib):
  return z/xib

def zhat1(xh,qT,Q):
  return (1-xh)/((qT**2/Q**2-1)*xh+1)

def xhat1(zh,qT,Q):
  return (1-zh)/((qT**2/Q**2-1)*zh+1)

def zhatFixed(z,xia,qT,Q):
  return zhat1(xhat(z,xia),qT,Q)

def xhatFixed(x,xib,qT,Q):
  return xhat1(zhat(x,xib),qT,Q)

def xiamin(w,x,z):
  return w**2/(1-z)+x

def xibmin(w,x,z):
  return w**2/(1-x)+z

def doubleu(x,z,qT,Q):
  return (qT/Q)*np.sqrt(x*z)

def xiaval(x,z,qT,Q,xib):
  return  x*((qT**2/Q**2)+(xib/z-1))/(xib/z-1) 

def xibval(x,z,qT,Q,xib):
  return  z*((qT**2/Q**2)+(xia/x-1))/(xia/x-1) 

def A2A1(y):
  return -y**2/((y-1)**2+1)

def Coeffa_qq( Q,qT,xia,x,y,mu):
  return 2*CF\
  *xhat(x,xia)*zhatFixed(x,xia,qT,Q)\
  *((1/qT**2)*((Q**4/(xhat(x,xia)**2*zhatFixed(x,xia,qT,Q)**2))\
    +(Q**2-qT**2)**2)+6*Q**2+4*Q**2*A2A1(y))

def Coeffb_qq( Q,qT,xib,z,y,mu):
  return 2*CF\
  * zhat(z,xib)*xhatFixed(z,xib,qT,Q)\
  *((1/qT**2)*((Q**4/(zhat(z,xib)**2*xhatFixed(z,xib,qT,Q)**2))\
    +(Q**2-qT**2)**2)+6*Q**2+4*Q**2*A2A1(y))

def Coeffa_qg( Q,qT,xia,x,y,mu):
  return xhat(x,xia)*(1-xhat(x,xia))\
  *((Q**4/qT**2)*(1/(xhat(x,xia)*zhatFixed(x,xia,qT,Q))**2\
  -2/(xhat(x,xia)*zhatFixed(x,xia,qT,Q))+2)\
  +2*Q**2*(5-1/xhat(x,xia)-1/zhatFixed(x,xia,qT,Q))+8*Q**2*A2A1(y))

def Coeffb_qg( Q,qT,xib,z,y,mu):
  return xhatFixed(z,xib,qT,Q)*(1-xhatFixed(z,xib,qT,Q))\
  *((Q**4/qT**2)*(1/(xhatFixed(z,xib,qT,Q)*zhat(z,xib))**2\
  -2/(xhatFixed(z,xib,qT,Q)*zhat(z,xib))+2)\
  +2*Q**2*(5-1/xhatFixed(z,xib,qT,Q)-1/zhat(z,xib))+8*Q**2*A2A1(y))

def qT_tilde( qt,zh):
  return  zh*qt/(1-zh)

def Coeffa_gq(Q,qT,xia,x,y,mu):
  return 2*CF*xhat(x,xia)*(1-zhatFixed(x,xia,qT,Q))\
  *((1/qT_tilde(qT,zhatFixed(x,xia,qT,Q))**2)\
  *((Q**4/(xhat(x,xia)**2*(1-zhatFixed(x,xia,qT,Q))**2))\
  +(Q**2-qT_tilde(qT,zhatFixed(x,xia,qT,Q))**2)**2)+6*Q**2+4*Q**2*A2A1(y))

def Coeffb_gq( Q,qT,xib,z,y,mu):
  return 2*CF*(1-zhat(z,xib))*xhatFixed(z,xib,qT,Q)\
  *((1/qT_tilde(qT,zhat(z,xib))**2)\
  *((Q**4/((1-zhat(z,xib))**2*xhatFixed(z,xib,qT,Q)**2))\
  +(Q**2-qT_tilde(qT,zhat(z,xib))**2)**2)+6*Q**2+4*Q**2*A2A1(y))

def get_terms1(xia,x,y,z,mu,Q,qT,qPDF,gPDF,qFF,gFF):
  return ((xhat(x, xia)*zhatFixed(x,xia,qT,Q)/(xia - x))*(\
          Coeffa_qq(Q,qT,xia,x,y,mu)*qPDF*qFF\
        + Coeffa_gq(Q,qT,xia,x,y,mu)*qPDF*gFF\
        + Coeffa_qg(Q,qT,xia,x,y,mu)*gPDF*qFF)).real
  
def get_terms2(xib,x,y,z,mu,Q,qT,qPDF,gPDF,qFF,gFF):
  return ((zhat(z, xib)*xhatFixed(z,xib,qT,Q)/(xib - z))*(\
         Coeffb_qq(Q,qT,xib,z,y,mu)*qPDF*qFF\
        + Coeffb_gq(Q,qT,xib,z,y,mu)*qPDF*gFF\
        + Coeffb_qg(Q,qT,xib,z,y,mu)*gPDF*qFF)/Q**4).real    
 
def Mbranch1(xia,x,y,z,mu,Q,qT,tar,had):
  pdf,ff=get_pdf_ff(xia,z/zhatFixed(x,xia,qT,Q),mu,tar,had)
  out=0
  Nf=aS.get_Nf(Q**2)
  for i in range(1,2*Nf+1,2):
    qPDF,gPDF=pdf[i],pdf[0]
    qFF ,gFF =ff[i],ff[0]
    out+=couplings[i]* get_terms1(xia,x,y,z,mu,Q,qT,qPDF,gPDF,qFF,gFF)
    qPDF,gPDF=pdf[i+1],pdf[0]
    qFF ,gFF =ff[i+1],ff[0]
    out+=couplings[i]* get_terms1(xia,x,y,z,mu,Q,qT,qPDF,gPDF,qFF,gFF)
  alpi=aS.get_alphaS(mu**2)/np.pi
  out*=np.pi*alfa**2*((1-y)**2+1)/Q**8/z**2*alpi/2
  return out

def Mbranch2(xib,x,y,z,mu,Q,qT,tar,had):
  pdf,ff=get_pdf_ff(x/xhatFixed(z,xib,qT,Q),xib,mu,tar,had)
  out=0
  Nf=aS.get_Nf(Q**2)
  for i in range(1,2*Nf+1,2):
    qPDF,gPDF=pdf[i],pdf[0]
    qFF ,gFF =ff[i],ff[0]
    out+=couplings[i]* get_terms2(xib,x,y,z,mu,Q,qT,qPDF,gPDF,qFF,gFF)
    qPDF,gPDF=pdf[i+1],pdf[0]
    qFF ,gFF =ff[i+1],ff[0]
    out+=couplings[i]* get_terms2(xib,x,y,z,mu,Q,qT,qPDF,gPDF,qFF,gFF)
  alpi=aS.get_alphaS(mu**2)/np.pi
  out*=np.pi*alfa**2*((1-y)**2+1)/Q**8/z**2*alpi/2
  return out

def get_pdf_ff(x,z,mu,tar,had):
  """ 
  set order as -5,-4...0.1,2,3...
  """
  f=pdf.get_f(x,mu**2)  
  if tar=='n':
    u,d = f[1],f[3]
    f[1],f[3] = d,u
  if tar=='d':
    u,d = f[1],f[3]
    f[1],f[3] = 0.5*(u+d),0.5*(u+d)
  d=ff.get_f(z,mu**2,had)
  return f,d

def integrator(f,xmin,xmax):
    return quad(f,xmin,xmax)[0]
    #return fixed_quad(np.vectorize(f),xmin,xmax,n=200)[0]

def get_xsec(x,z,Q,qT,y,E,tar,had): 
  """
  dsig / dx dQ2 dz dpT2
  """
  mu=Q
  _Mbranch1=lambda xia: Mbranch1(xia,x,y,z,mu,Q,qT,tar,had)
  _Mbranch2=lambda xib: Mbranch2(xib,x,y,z,mu,Q,qT,tar,had)
  #gamma= integrator(_Mbranch1,x+doubleu(x,z,qT,Q),1)/2\
  #      +integrator(_Mbranch2,z+doubleu(x,z,qT,Q),1)/2
  return integrator(_Mbranch1,x+doubleu(x,z,qT,Q)**2/(1-z),1)

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
  #print Mbranch1(0.5,x,y,z,Q,Q,qT,tar,had)





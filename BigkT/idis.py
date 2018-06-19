#!/usr/bin/env python
import sys,os
import numpy as np
from scipy.integrate import quad,fixed_quad
try:
  import lhapdf
  pdf=lhapdf.mkPDF("CJ15nlo", 0)
  mc=pdf.quarkThreshold(4)
  mb=pdf.quarkThreshold(5)
except:
  print 'lhapdf not found. you wont be able to compute idis'


TR=0.5
CF=4./3.
alfa=1/137.036
M=0.93891897
apU=4.0/9.0
apD=1.0/9.0
couplings={}
couplings['p']={1:apD,2:apU,3:apD,4:apU,5:apD}
couplings['n']={1:apU,2:apD,3:apD,4:apU,5:apD}
fmap={}

F2={'p':{},'n':{}}
FL={'p':{},'n':{}}

 
def integrator(f,xmin,xmax,method='gauss',n=100):
  f=np.vectorize(f)
  if method=='quad':
    return quad(f,xmin,xmax)[0]
  elif method=='gauss':
    return fixed_quad(f,xmin,xmax,n=n)[0]
  
def log_plus(z,f,x):
  return np.log(1-z)/(1-z)*(f(x/z)/z-f(x)) + 0.5*np.log(1-x)**2*f(x)/(1-x)

def one_plus(z,f,x):
  return 1/(1-z)*(f(x/z)/z-f(x))+ np.log(1-x)*f(x)/(1-x)

def C2q(z,f,x):
  return CF*(2*log_plus(z,f,x)-1.5*one_plus(z,f,x)\
    +(-(1+z)*np.log(1-z)-(1+z*z)/(1-z)*np.log(z)+3+2*z)*f(x/z)/z\
    -(np.pi**2/3+4.5)*f(x)/(1-x))
  
def C2g(z,f,x):
  return 0.5*(((1-z)**2+z*z)*np.log((1-z)/z)-8*z*z+8*z-1)*f(x/z)/z

def CLq(z,f,x):
  return 2*CF*z*f(x/z)/z #<--- note prefactor 2, instead of 4 used by MVV
  
def CLg(z,f,x):
  return 4*z*(1-z)*f(x/z)/z

def qplus(x,Q2,Nf,tar):
  output=0
  for i in range(1,Nf+1):
    output+=couplings[tar][i]*(pdf.xfxQ2(i,x,Q2)/x+pdf.xfxQ2(-i,x,Q2)/x)
  return output

def glue(x,Q2,Nf,tar):
  output=0
  for i in range(1,Nf+1):
    output+=2*couplings[tar][i]
  return output*pdf.xfxQ2(21,x,Q2)/x
    
def integrand_F2(x,z,Q2,Nf,tar):
  return C2q(z,lambda y:qplus(y,Q2,Nf,tar),x) + C2g(z,lambda y:glue(y,Q2,Nf,tar),x)
  
def integrand_FL(x,z,Q2,Nf,tar):
  return CLq(z,lambda y:qplus(y,Q2,Nf,tar),x) + CLg(z,lambda y:glue(y,Q2,Nf,tar),x)
  
def get_F2(x,Q2,tar):
  if (x,Q2) not in F2[tar]:
    tar=tar
    alphaS = pdf.alphasQ2(Q2)
    Nf=3
    if Q2>mc**2: Nf+=1
    if Q2>mb**2: Nf+=1
    LO=qplus(x,Q2,Nf,tar)
    integrand=lambda z:integrand_F2(x,z,Q2,Nf,tar)
    NLO=integrator(integrand,x,1)
    F2[tar][(x,Q2)]=x*(LO+alphaS/np.pi/2.0*NLO)
  return F2[tar][(x,Q2)]

def get_FL(x,Q2,tar):
  if (x,Q2) not in FL[tar]:
    tar=tar
    alphaS = pdf.alphasQ2(Q2)
    Nf=3
    if Q2>mc**2: Nf+=1
    if Q2>mb**2: Nf+=1
    integrand=lambda z:integrand_FL(x,z,Q2,Nf,tar)
    NLO=integrator(integrand,x,1)
    FL[tar][(x,Q2)]= x*alphaS/np.pi/2.0*NLO
  return FL[tar][(x,Q2)]

def get_F1(x,Q2,tar):
  F2=get_F2(x,Q2,tar)
  FL=get_FL(x,Q2,tar)
  return ((1+4*M**2/Q2*x**2)*F2-FL)/(2*x)

def get_xsec(x,y,Q2,tar):
  """ 
  returns  dsig / dx dQ2
  Eq: 2.15 DIS by Devenish, Cooper-Sarkar
  """
  if   tar=='p': F2,FL=get_F2(x,Q2,'p')+get_FL(x,Q2,'p')
  elif tar=='n': F2,FL=get_F2(x,Q2,'n')+get_FL(x,Q2,'n')
  elif tar=='d': 
    F2p,FLp=get_F2(x,Q2,'p'),get_FL(x,Q2,'p')
    F2n,FLn=get_F2(x,Q2,'n'),get_FL(x,Q2,'n')
    F2=0.5*(F2p+F2n)
    FL=0.5*(FLp+FLn)
  else:
    raise ValueError('tar not defined: ',tar)
  return 2*np.pi*alfa**2/x/Q2**2*((1+(1-y)**2)*F2-y**2*FL)

if __name__== "__main__":

  x=0.017
  Q=1.87
  y=0.695
  E=160
  tar='p'
  had='h+'
  order=0
  print get_xsec(x,y,Q**2,tar)








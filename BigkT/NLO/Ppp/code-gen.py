#!/usr/bin/env python
import sys,os
import copy
import numpy as np
import sympy as sp
import pylab as py
from tools.tools import load,save
from sympy.parsing.sympy_parser import parse_expr
lprint=lambda expression: display(Math(sp.latex(expression)))
from tools.bar import BAR
import mpmath as mp
from mpmath import fp

nf   = sp.Symbol('nf')
s   = sp.Symbol('s')
t  = sp.Symbol('t')
u  = sp.Symbol('u')
mu  = sp.Symbol('mu')
s23 = sp.Symbol('s23')
Q   = sp.Symbol('Q')
delta=sp.Function('delta')
Plus1B=sp.Function('Plus1B')
Plus2B=sp.Function('Plus2B')
PolyLOG=sp.Function('PolyLOG')
epsilon=sp.Symbol('epsilon')
B=sp.Symbol('B')
g=sp.Symbol('g')
gp=sp.Symbol('gp')
a=sp.Wild('a')
b=sp.Wild('b')

def gen_ps(s23val):
    ps={}
    ps['x']=0.01
    ps['z']=0.3
    ps['Q']=5.0
    ps['qT']=3.0
    ps['s23']=s23val
    ps['xi']=0.4
    ps['xh']=ps['x']/ps['xi']
    ps['zh']=(-ps['s23']+ps['Q']**2*(1-ps['xh'])/ps['xh'])/(ps['qT']**2+ps['Q']**2*(1-ps['xh'])/ps['xh'])
    ps['zeta']=ps['z']/ps['zh']
    ps['s']=(1-ps['xh'])/ps['xh']*ps['Q']**2
    ps['t']=-(1-ps['zh'])*ps['Q']**2-ps['zh']*ps['qT']**2
    ps['nf']=4.0
    ps['B']=ps['Q']**2*(1/ps['xh']-1)*(1-ps['z'])-ps['z']*ps['qT']**2
    return ps

def evaluate(exp,s23val,verb=False):
    ps=gen_ps(s23val)
    test=exp.subs(g,1.)
    test=test.subs(Q,ps['Q'])
    test=test.subs(s23,ps['s23'])
    test=test.subs(s,ps['s'])
    test=test.subs(t,ps['t'])
    test=test.subs(mu,ps['Q'])
    test=test.subs(sp.EulerGamma,np.euler_gamma)
    test=test.subs(B,ps['B'])
    test=test.subs(nf,ps['nf'])
    test=test.subs(sp.pi,np.pi)
    if verb:
        print ps['x']
        print ps['z']
        print ps['Q']
        print ps['qT']
        print ps['s23']
        print ps['xi']
        print ps['xh']
        print ps['zh']
        print ps['zeta']
        print ps['s']
        print ps['t']
        print ps['B']
    return test

def convert(exp):
    sexp=str(exp).replace('sqrt','np.sqrt')
    sexp=sexp.replace('log','np.log')
    sexp=sexp.replace('pi','np.pi')
    sexp=sexp.replace('Abs','np.abs')
    sexp=sexp.replace('polylog','fp.polylog')
    return sexp

def gen_code(channel,case,exp1,exp2,exp3,exp4):
    code=[]
    code.append(r'#!/usr/bin/env python')
    code.append(r'import numpy as np')
    code.append(r'from mpmath import fp')
    code.append(r'import numpy as np')
    code.append(r'EulerGamma=np.euler_gamma')
    code.append(r'def PolyLOG(s, z):')
    code.append(r'    return fp.polylog(s,z)')
    code.append(r'def regular(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,nf=None):')
    code.append('    return %s'%convert(exp1))
    code.append(r'def delta(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,B=None,nf=None):')
    code.append('    return %s'%convert(exp2))
    code.append(r'def plus1B(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,B=None,nf=None):')
    code.append('    return %s'%convert(exp3))
    code.append(r'def plus2B(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,B=None,nf=None):')
    code.append('    return %s'%convert(exp4))

    code=[l+'\n' for l in code]
    F=open('chn%d%s.py'%(channel,case),'w')
    F.writelines(code)
    F.close()
    os.system('chmod +x chn%d%s.py'%(channel,case));

def gen_fcode(channel,case,exp1,exp2,exp3,exp4):
    code=[]
    code.append(r'#!/usr/bin/env python')
    code.append(r'import numpy as np')
    code.append(r'from mpmath import fp')
    code.append(r'from numba import jit')
    code.append(r'import numpy as np')
    code.append(r'EulerGamma=np.euler_gamma')
    code.append(r'@jit(cache=True)')
    code.append(r'def _PolyLOG(s, z):')
    code.append(r'    tol = 1e-10')
    code.append(r'    l = 0')
    code.append(r'    k = 1')
    code.append(r'    zk = z')
    code.append(r'    while 1:')
    code.append(r'        term = zk / k**s')
    code.append(r'        l += term')
    code.append(r'        if abs(term) < tol:')
    code.append(r'            break')
    code.append(r'        zk *= z')
    code.append(r'        k += 1')
    code.append(r'    return l')
    code.append(r'@jit(cache=True)')
    code.append(r'def PolyLOG(s, z):')
    code.append(r'    #return fp.polylog(s,z)')
    code.append(r'    #if abs(z) > 0.75:')
    code.append(r'    #  return -PolyLOG(s,1-z) + np.pi**2/6 - np.log(z)*np.log(1-z)')
    code.append(r'    if abs(z) >1: ')
    code.append(r'      return -PolyLOG(s, 1/z) - np.pi**2/6 - 0.5*np.log(-z)**2')
    code.append(r'    return _PolyLOG(s, z)')
    code.append(r'@jit(cache=True)')
    code.append(r'def regular(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,nf=None):')
    code.append('    return %s'%convert(exp1))
    code.append(r'@jit(cache=True)')
    code.append(r'def delta(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,B=None,nf=None):')
    code.append('    return %s'%convert(exp2))
    code.append(r'@jit(cache=True)')
    code.append(r'def plus1B(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,B=None,nf=None):')
    code.append('    return %s'%convert(exp3))
    code.append(r'@jit(cache=True)')
    code.append(r'def plus2B(g=None,gp=None,s=None,t=None,Q=None,s23=None,mu=None,B=None,nf=None):')
    code.append('    return %s'%convert(exp4))

    code=[l+'\n' for l in code]
    F=open('fchn%d%s.py'%(channel,case),'w')
    F.writelines(code)
    F.close()
    os.system('chmod +x fchn%d%s.py'%(channel,case));

def get_s23zero_limit(exp):

    lexp=(exp.subs(s23,0)).as_ordered_terms()
    new_exp=0
    bar=BAR("proceesing",len(lexp))
    for term in lexp:
        dummy1=1
        for f in term.as_ordered_factors():
            dummy2=0
            for ft in f.as_ordered_terms():
                dummy2+=sp.powdenest(ft,force=True).simplify()
            dummy1*=dummy2
        new_exp+=dummy1
        bar.next()
    bar.finish()
    
    lexp=new_exp.as_ordered_terms()
    new_exp=sp.S(0)
    for term in lexp:
        if evaluate(term,0,verb=False)==0: continue
        new_exp+=term.replace(sp.log(a),sp.log(sp.Abs(a)))
    
    return new_exp

def load_exp(ch,case,part):
    L=open('channels/ch%d%s.%s'%(ch,case,part)).readlines()
    L=[l.strip() for l in L]
    L=[l.replace(r'&\n','') for l in L]
    expression=''
    for i in range(len(L)):
        l=L[i]
        l=l.replace('\\[Mu]','mu')
        l=l.replace('[','(')
        l=l.replace(']',')')
        l=l.replace('^','**')
        l=l.replace('E4Pi','(EulerGamma - sp.log(4*sp.pi))')
        l=l.replace('PolyLog','PolyLOG')
        l=l.replace('Log','sp.log')
        l=l.replace('Pi','sp.pi')
        l=l.replace('EulerGamma','sp.EulerGamma')
        l=l.replace('Sqrt','sp.sqrt')
        l=l.rstrip('\\')
        l=l.replace('\"','')
        l=l.strip()
        if l!='': expression+=l
    exec 'exp=%s'%expression
    return exp


for chn in range(1,7):
  print 'generating code for channel ',chn
  for case in ['A','B','C']:
    exp1=load_exp(chn,case,'regular')
    exp2=load_exp(chn,case,'delta')
    exp3=load_exp(chn,case,'plus1B')
    exp4=load_exp(chn,case,'plus2B')

    if exp1!=0: print 'exp1',chn,case
    if exp2!=0: print 'exp2',chn,case
    if exp3!=0: print 'exp3',chn,case
    if exp4!=0: print 'exp4',chn,case

    #gen_code(chn,exp1,exp2,exp3,exp4)
    gen_fcode(chn,case,exp1,exp2,exp3,exp4)













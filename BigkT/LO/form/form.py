#!/usr/bin/env python
import sys,os
sys.path.insert(1,'../') 
from tools.tools import save, load,checkdir
from qgraf.qgraf import QGRAF
import subprocess

class FORM:

  def __init__(self,conf):
    checkdir('.form')
    self.conf=conf    
    qgraf=conf['qgraf']
    self.AMPS=qgraf['AMPS']
    self.SQAMPS=qgraf['SQAMPS']
    self.external_states=qgraf['external states']
    if 'verb' in conf: verb=conf['verb']
    else: verb=False
    self.run(mode=0,verb=verb)
    self.run(mode=1,verb=verb)
    self.run(mode=2,verb=verb)

  def get_headers(self):

    AMPS=self.AMPS

    # get particles
    fields=[]
    masses=[]
    for iamp in AMPS.keys():
      for ipart in AMPS[iamp].keys():
        if ipart=='chains' or ipart=='nonchains':continue
        if any(['particle'==k for k in AMPS[iamp][ipart].keys()]):
          fields.append(AMPS[iamp][ipart]['particle'])
        if any(['mass'==k for k in AMPS[iamp][ipart].keys()]):
          masses.append(AMPS[iamp][ipart]['mass'])
    fields=sorted(set(fields))
    masses=sorted(set(masses))
    Masses=''
    for mass in masses: Masses+='%s,'%mass
    Masses=Masses[:-1]
    Fields=''
    for f in fields: Fields+='%s,'%f
    Fields=Fields[:-1]

    H=[]
    H.append('*Masses')
    H.append('Symbols %s,M;'%Masses)
    H.append('*Fields')
    H.append('Symbols %s,phi;'%Fields)
    H.append('*Incoming spin indices')
    H.append('Symbols r1,...,r12,r;')
    H.append('*Outgoing spin indices')
    H.append('Symbols s1,...,s12,s;')
    H.append('*couplings')
    H.append('Symbols gl,gu,gd,,gs;')
    H.append('Symbols alpha,alphaS;')
    H.append('*mandelstan')
    H.append('Symbols s,t1,u1,Q2;')
    H.append('Symbols t2,t3,u2,u3,s12,s13,s23;')
    H.append('*aux')
    H.append('Symbols x,pi,ave;')
    H.append('*spinor indices')
    H.append('Index sc1,...,sc12,sc;')
    H.append('*dimension')
    H.append('Symbols D;')
    H.append('Function UB,U,V,VB;')
    H.append('CFunctions den,eps,epsB,X,XB;')

    #H.append('Tensor Tr(cyclic),T,Tp,f(antisymmetric),ff(rcyclic);')
    #H.append('Function TM,TT;')
    H.append('Symbols a,nf,NF,NA,cF,cA,[cF-cA/6],[u8,u3],[u8,u5],[d8,d3],[d8,d5];')


    H.append('CFunction Tr(cyclic);')# Indicates the traces.
    H.append('CFunction T, Tp, f(antisymmetric);')
    #H.append('Symbols a,nf,NF,NA;')
    #H.append('Indices i1=NF,i2=NF,i3=NF,i4=NF;')
    #H.append('Indices j1=NA,j2=NA,j3=NA;')

    H.append('Dimension D;')
    H.append('*external momentum  ')
    H.append('Vector q1,...,q12,q;')
    H.append('Vector p1,...,p12,p;')
    H.append('*loop momenta')
    H.append('Vector k1,...,k5,k;')


    #H.append('*loop momentum')
    #H.append('Vector k;')
    H.append('AutoDeclare Index o;')
    H.append('Dimension NA;')
    H.append('AutoDeclare Index j;')
    H.append('Dimension NF;')
    H.append('AutoDeclare Index i;')
    #H.append('AutoDeclare Symbol I;')

    H.append('*')
    return H

  def get_locals(self):
    Locals=[]
    for k in self.SQAMPS.keys():
      Locals.append('Local %s = %s;'%(k,self.SQAMPS[k]))
    local='Local [total] = '
    for k in self.SQAMPS.keys(): local+='ave*'+k+'+'
    local=local[:-1]+';'
    Locals.append(local)
    Locals.append('*')
    return Locals

  def get_instructions(self):

    code=open('%s/%s'%(self.conf['path2form'],self.conf['projection']),'r').readlines()
    code=[i.strip() for i in code]

    if self.mode==0:
      # do nothing!!
      code=[]

    elif self.mode==1:
      # dot not take traces
      for i in range(len(code)):
        #   code=[]
        if 'Trace' in code[i] or 'den(x?)=1/x' in code[i] or 'id den(p?,M?)=den(p.p-M*M)' in code[i] or 'id gs=1' in code[i]:
          code[i]='*'+code[i]

    elif self.mode==2:
      # do everything
      pass

    for k in self.SQAMPS.keys():
      code.append('Print %s ;'%(k))
    code.append('Print +s [total];')
    # print code
    return code

  def check_form_output(self,lines):
    flag=False
    for l in lines:
      if 'Time' in l: flag=True
    if flag==True: return 0
    else:
      msg='''
      mode:%d
      ERR: form did not run correctly. run $./form formcalc.frm to see what is the problem
      This is because the default form template (formcalc.frm) needs to be adjusted
      for the specified process
      '''
      msg=(msg%self.mode).replace('\t','')
      print msg
      sys.exit()

  def run(self,mode=0,verb=False):
    self.mode=mode

    CODE=self.get_headers()
    CODE.extend(self.get_locals())
    CODE.extend(self.get_instructions())
    CODE=[c+'\n' for c in CODE]
    cwd=os.getcwd()
    self.formcode='%s/.form/%s-%s'%(cwd,self.conf['label'],self.conf['projection'])
    F=open(self.formcode,'w')
    F.writelines(CODE)
    F.close()

    os.chdir(self.conf['path2form'])
    cmd=['./form.exe',self.formcode]
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE)
    os.chdir(cwd)
    L=[]
    while 1:
      l=p.stdout.readline()
      if not l: break
      l=l.replace('\n','')
      if verb: print l
      L.append(l)
    if verb: print '='*80

    if verb:
      for l in L: print l

    self.check_form_output(L)

    # isolate relevant outout. Everything after
    # the last print of Time is useful
    for i in range(len(L)):
      if 'Time' in L[i]: I=i
    I+=3
    L=L[I:]
    L=[l for l in L if 'sec' not in l]

    # get lines where [a*ca] starts
    I=[i for i in range(len(L)) if '] =' in L[i]]

    self.data={}
    for i in range(len(I)):
      ini=I[i]
      if i<len(I)-1: fin=I[i+1]
      else: fin=len(L)
      l=L[ini:fin]
      if l[0].strip().endswith(';'):
        k,exp=l[0].strip().split('=')
        k=k.replace('[','').replace(']','').strip()
        self.data[k]=[exp.replace(';','').strip()]
      else:
        k=l[0].strip().split('=')[0].replace('[','').replace(']','').strip()
        self.data[k]=[r.replace(';','').strip() for r in l[1:]]
        self.data[k]=[r for r in self.data[k] if r!='']
      ## concat  all terms 
      #exp=''
      #for term in self.data[k]: exp+=term
      #self.data[k]=exp


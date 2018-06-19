#!/usr/bin/env python
import sys,os
import copy
import numpy as np
from multiprocessing import Process,Pipe
import argparse
import pandas as pd
import time
import idis
import sidisV01
import sidisV10
#from daleo import timba
#timba.setup_path('./daleo/')

def lprint(msg):
  sys.stdout.write('\r')
  sys.stdout.write(msg)
  sys.stdout.flush()

class MPROC:
    
  """ simple parallization"""
  
  def __init__(self,func,args,ncores=10):
    self.ncores=ncores
    self.args=args
    self.func=func
    n=len(self.args)/(self.ncores)
    self.ini,self.fin=[],[]
    for i in range(self.ncores):
        self.ini.append(i*n)
        self.fin.append((i+1)*n)
    self.fin[-1]=len(self.args)
    self.pipes=[Pipe() for i in range(self.ncores)]
      
  def func_wrap(self,args,pipe):
    output=[]
    npts=len(args)
    for i in range(len(args)):
      msg='proceesing %d/%d x %d'%(i+1,npts,self.ncores)
      lprint(msg)
      output.append(self.func(args[i]))
    pipe.send(output)

  def run(self):
    P = [Process(target=self.func_wrap, \
                  args=(self.args[self.ini[i]:self.fin[i]],self.pipes[i][1],)) \
                  for i in range(self.ncores)]
    for p in P: p.start() 
    for p in P: p.join() 
    output=[]
    for pipe in self.pipes: output.extend(pipe[0].recv())
    return output

def calc4xlsx(args):
  data=pd.read_excel(args.xlsx)

  K=[]
  for i in range(len(data.index)):
    x=data.x[i]
    z=data.z[i]
    qT=data.qT[i]
    Q=data.Q[i]    
    E=data.Ebeam[i]    
    y=data.y[i] 
    tar='d'#data.tar[i]
    had='h+'#data.had[i]

    #M=0.938
    #xn=2*x/(1+np.sqrt(1+4*x**2*M**2/Q**2))
    #x=xn

    #         0 1 2  3 4 5   6   7
    K.append([x,z,Q,qT,y,E,tar,had])


  if args.task==0:
    print "\ncomputing  IDIS @ NLO"
    func=lambda k:  idis.get_xsec(k[0],k[4],k[2]**2,k[6])
    results=MPROC(func,K,args.cores).run()
    data['idisNLO']=pd.Series(results,index=data.index)

  if args.task==1:
    print "\ncomputing  SIDISV10 @ LO"
    sidisV10LO=lambda k:  sidisV10.get_xsec(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],0)
    results=MPROC(sidisV10LO,K,args.cores).run()
    data['sidisLO']=pd.Series(results,index=data.index)

  if args.task==2:
    print "\ncomputing  SIDIS @ NLO"
    sidisNLO=lambda k:  sidisV10.get_xsec(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],1)
    results=MPROC(sidisNLO,K,args.cores).run()
    data['sidisNLO']=pd.Series(results,index=data.index)


  if args.task==10:
    print "\ncomputing  daleo @ born"
    func=lambda k:  timba.get_xsec(k[0],k[2],k[1],k[3],k[5],k[7],k[6],'born')
    results=MPROC(func,K,args.cores).run()
    data['daleoLO']=pd.Series(results,index=data.index)

  if args.task==11:
    print "\ncomputing  daleo @ delta"
    func=lambda k:  timba.get_xsec(k[0],k[2],k[1],k[3],k[5],k[7],k[6],'delta')
    results=MPROC(func,K,args.cores).run()
    data['daleoNLOdelta']=pd.Series(results,index=data.index)

  if args.task==12:
    print "\ncomputing  daleo @ plus"
    func=lambda k:  timba.get_xsec(k[0],k[2],k[1],k[3],k[5],k[7],k[6],'plus')
    results=MPROC(func,K,args.cores).run()
    data['daleoNLOplus']=pd.Series(results,index=data.index)

  if args.task==13:
    print "\ncomputing  daleo @ rest"
    func=lambda k:  timba.get_xsec(k[0],k[2],k[1],k[3],k[5],k[7],k[6],'rest')
    results=MPROC(func,K,args.cores).run()
    data['daleoNLOrest']=pd.Series(results,index=data.index)

  data.to_excel(args.xlsx)


  #if args.task==5:
  #  print "\ncomputing  <xi> @ LO"
  #  sidisV10LO=lambda k:  sidisV10.get_xsec(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],0,'xi')
  #  results=MPROC(sidisV10LO,K,args.cores).run()
  #  data['<xi>']=pd.Series(results,index=data.index)

  #if args.task==6:
  #  print "\ncomputing  <zeta> @ LO"
  #  sidisV10LO=lambda k:  sidisV10.get_xsec(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7],0,'zeta')
  #  results=MPROC(sidisV10LO,K,args.cores).run()
  #  data['<zeta>']=pd.Series(results,index=data.index)

  #if args.task==0 or  args.task==2:
  #  print "\ncomputing  SIDISV01 @ LO"
  #  sidisV01LO=lambda k:  sidisV01.get_xsec(k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7])
  #  results=MPROC(sidisV01LO,K,args.cores).run()
  #  data['sidisV01LO']=pd.Series(results,index=data.index)


if __name__=="__main__":

  ap = argparse.ArgumentParser()
  msg ='0: IDIS @ NLO '
  msg+='1: SIDIS @ LO '
  msg+='2: SIDIS @ NLO '
  msg+='10: daleo @ LO '
  msg+='11: daleo @ NLO delta'
  msg+='12: daleo @ NLO plus'
  msg+='13: daleo @ NLO rest'
  ap.add_argument('-t','--task', type=int,default=0,help=msg,required=True)
  ap.add_argument('-x','--xlsx', type=str,default='',help='*.xlsx with kinematics')
  ap.add_argument('-c','--cores',type=int,default=2,help='cpu cores')
  args = ap.parse_args()
  
  calc4xlsx(args)








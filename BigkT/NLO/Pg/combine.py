#!/usr/bin/env python
import sys,os
import subprocess
from sympy.parsing.sympy_parser import parse_expr
from tools.tools import save,checkdir


def gen_channels(channel):
  checkdir("channels")

  # collect file names from mdata
  F=os.listdir('mdata')
  F=[f for f in F if '-%d'%channel in f]
  D={}
  D['counter']={'file':None,'var':'counterpiece'}
  D['subtraction']={'file':None,'var':'subtractionpiece'}
  D['virtual']={'file':None,'var':'virtualpiece'}
  D['real']={'file':None,'var':'realpiece'}
  for f in F:
    for k in D:
      if k in f: D[k]['file']=f


  cases={}
  cases['A']='g*g'
  cases['B']='g*gp'
  cases['C']='gp*gp'

  for case in ['A','B','C']:

    # start code lines
    L=[]

    # add script header
    L.append('#!/apps/mathematica/bin/math -script')
    #L.append('#!/usr/local/bin/math -script')  ##jogh: my system

    # load special modules
    L.append('<<fastlimit.m')##jogh: load limit function

    # load expression from mdata
    for k in D:
      if D[k]!=None: L.append('<<mdata/%s'%D[k]['file'])
      else: L.append('%s=0;'%D[k]['var'])

    # construct toal contribution for this channel
    L.append(r'total = 0;')
    L.append(r'total += 1*realpiece;')
    L.append(r'total += 1*counterpiece*\[Delta][s23];')
    L.append(r'total -= 1*subtractionpiece;')
    if channel<4:
      L.append(r'total += 1*(virtualpiece[[1]] + virtualpiece[[2]] + virtualpiece[[3]] + virtualpiece[[4]])*\[Delta][s23];')

    # remove poles (cancellation of poles are tested in combine.nb)
    L.append(r'total = total /. \[Epsilon] -> 1/\[Epsilon] /. \[Epsilon] -> 0/. Subscript[n, f]->nf;')

    # isolate charge combination
    L.append(r'total=Coefficient[total,%s];'%cases[case])

    # isolate coeffs of delta, plus1B, plus2B and rest
    L.append('exp1 = total /. Plus1B[B, s23] -> 0 /.Plus2B[B, s23] -> 0 /. \[Delta][s23] -> 0;')
    L.append('exp2 = total /. Plus1B[B, s23] -> 0 /.Plus2B[B, s23] -> 0 /. \[Delta][s23] -> 1;')
    L.append('exp2 = exp2 - exp1;')
    L.append('exp3 = total /. Plus1B[B, s23] -> 1 /.Plus2B[B, s23] -> 0 /. \[Delta][s23] -> 0;')
    L.append('exp3 = exp3 - exp1;')
    L.append('exp4 = total /. Plus1B[B, s23] -> 0 /.Plus2B[B, s23] -> 1 /. \[Delta][s23] -> 0;')
    L.append('exp4 = exp4 - exp1;')
    
    ####jogh simplification and limit s23->0 (checking delta case which gives weird numbers now.)
    #L.append('exp2s23eq0=exp2//FastLimit//Simple1;')
    L.append('exp3s23eq0=exp3//FastLimit//Simple1;')
    L.append('exp4s23eq0=exp4//FastLimit//Simple1')

    #L.append('exp2=exp2//Simple1;')
    L.append('exp3=exp3//Simple1;')
    L.append('exp4=exp4//Simple1;')

    L.append(r'step1 = FortranForm[exp1 /. \[Mu] -> mu];')
    L.append(r'result1 = StringJoin @@ Riffle[With[{splits = StringSplit[ToString@step1, " "]}, Fold[If[StringLength[Last@#1] + StringLength[#2] > 60, Join[#1, {#2}], Join[Most[#1], {StringJoin[Last[#1], #2]}]] &, {First@splits}, Rest[splits]]], "  &\n"];')
    L.append('Write["channels/ch%d%s.regular", result1];'%(channel,case))

    L.append(r'step2 = FortranForm[exp2 /. \[Mu] -> mu];')
    L.append(r'result2 = StringJoin @@ Riffle[With[{splits = StringSplit[ToString@step2, " "]}, Fold[If[StringLength[Last@#1] + StringLength[#2] > 60, Join[#1, {#2}], Join[Most[#1], {StringJoin[Last[#1], #2]}]] &, {First@splits}, Rest[splits]]], "  &\n"];')
    L.append('Write["channels/ch%d%s.delta", result2];'%(channel,case))

    L.append(r'step3 = FortranForm[exp3 /. \[Mu] -> mu];')
    L.append(r'result3 = StringJoin @@ Riffle[With[{splits = StringSplit[ToString@step3, " "]}, Fold[If[StringLength[Last@#1] + StringLength[#2] > 60, Join[#1, {#2}], Join[Most[#1], {StringJoin[Last[#1], #2]}]] &, {First@splits}, Rest[splits]]], "  &\n"];')
    L.append('Write["channels/ch%d%s.plus1B", result3];'%(channel,case))

    L.append(r'step4 = FortranForm[exp4 /. \[Mu] -> mu];')
    L.append(r'result4 = StringJoin @@ Riffle[With[{splits = StringSplit[ToString@step4, " "]}, Fold[If[StringLength[Last@#1] + StringLength[#2] > 60, Join[#1, {#2}], Join[Most[#1], {StringJoin[Last[#1], #2]}]] &, {First@splits}, Rest[splits]]], "  &\n"];')
    L.append('Write["channels/ch%d%s.plus2B", result4];'%(channel,case))

    #for l in L: print l
    L=[l+'\n' for l in L]
    F=open('tmp.m','w')
    F.writelines(L)
    F.close()
    os.system('chmod +x tmp.m')
    os.system('./tmp.m')
    os.system('rm tmp.m')
  

for ch in range(1,7):
  print '++++++++ channel %s:' %ch
  gen_channels(ch)



#ielf.conf['path2qgraf']!/usr/bin/env python
import sys,os
sys.path.insert(1,'../') 
import numpy as np
from tools.tools import checkdir,save,load
import re
import subprocess

class QGRAF:

  def __init__(self,conf):

    self.conf=conf
    self.setup()
    self.gen_input()
    self.run_qgraf()
    if 'verb' in conf: verb=conf['verb']
    else: verb=False
    self.process_data(verb)

  def setup(self):
    self.external_states={}
    self.LI=range(1,10000)
    self.LImap={}

  def gen_input(self):
    L=[]
    L.append('output = \'<<output>>\' ;')
    L.append('style = \'<<style>>\' ;')
    L.append('model = \'<<model>>\';')
    L.append('in = <<in>>;')
    L.append('out = <<out>>;')
    L.append('loops = <<loops>>;')
    L.append('loop_momentum = <<loop_momentum>>;')
    L.append('options = <<options>>;')
    for i in range(len(L)):
      L[i]=L[i].replace('<<output>>','%s.qgraf'%self.conf['label'])
      L[i]=L[i].replace('<<style>>',self.conf['style'])
      L[i]=L[i].replace('<<model>>',self.conf['model'])
      L[i]=L[i].replace('<<in>>',self.conf['in'])
      L[i]=L[i].replace('<<out>>',self.conf['out'])
      L[i]=L[i].replace('<<loops>>','%d'%self.conf['loops'])
      L[i]=L[i].replace('<<loop_momentum>>',self.conf['loop_momentum'])
      L[i]=L[i].replace('<<options>>',self.conf['options'])
    for l in self.conf['extra']: L.append(l)
    L=[l+'\n' for l in L]
    F=open('%s/qgraf.dat'%self.conf['path2qgraf'],'w')
    F.writelines(L)
    F.close()

  def run_qgraf(self,verb=False):
    checkdir('.amps')
    cwd=os.getcwd()
    os.chdir(self.conf['path2qgraf'])
    subprocess.Popen(['./qgraf.exe'],stdout=subprocess.PIPE)
    #os.system('./qgraf.exe')
    os.system('rm qgraf.dat')
    os.system('mv %s.qgraf %s/.amps/'%(self.conf['label'],cwd))
    os.chdir(cwd)
    self.fname='.amps/%s.qgraf'%(self.conf['label'])

  def process_data(self,verb=False):
    self.get_RAW_AMPS(verb)
    self.get_AMPS(verb)
    self.get_SQAMPS(verb)
    self.data={}
    self.data['AMPS']=self.AMPS
    self.data['SQAMPS']=self.SQAMPS
    self.data['external states']=self.external_states

  def get_RAW_AMPS(self,verb=False):

    # load data and clean
    F=open(self.fname)
    L=F.readlines()
    F.close()
    L=[l.strip() for l in L]
    L=[l for l in L if l.startswith('*')==False if l!='']

    # get markers where amps starts
    I=[i for i in range(len(L)) if ':=' in L[i]]
    
    # extract amps and store them in dicts
    AMPS={}
    iamp=0
    for i in range(len(I)):
      iini=I[i]+1
      if i==len(I)-1: 
        ifin=len(L)
      else:
        ifin=I[i+1]
      parts=L[iini:ifin][:]
      AMP={}
      ipart=0
      for part in parts:
        key='part'+str(ipart)
        AMP[key]={'raw':part}
        elements=part[1:][:-1].split(',')
        for element in elements:
          name,value=element.split('=')
          AMP[key][name]=value
        ipart+=1

      if 'amps' in self.conf:
        if iamp in self.conf['amps']: 
          AMPS['a'+str(iamp)]=AMP
      else:
        AMPS['a'+str(iamp)]=AMP

      iamp+=1
    self.AMPS=AMPS

    if verb:
      print '='*80
      print 'RAW AMPS'
      for k in AMPS.keys():
        print '-'*80
        print 'amp = %s'%k
        for kk in AMPS[k].keys():
          print '  %s'%(kk)
          keys=sorted(AMPS[k][kk].keys())
          print '    %10s = %10s'%('part',AMPS[k][kk]['part'])
          for kkk in keys:
            if kkk=='part': continue
            print '    %10s = %10s'%(kkk,AMPS[k][kk][kkk])

  def get_scalar(self,iamp,ipart):
    d=self.AMPS[iamp][ipart]
    # register expression
    self.AMPS[iamp][ipart]['expression'] ='(%s)'%d['value']
    self.AMPS[iamp][ipart]['cexpression']='(%s)'%d['value']

  def get_pol(self,iamp,ipart):
    d=self.AMPS[iamp][ipart]
    particle=d['particle']
    if particle.endswith('b'): particle=particle[:-1]

    # setup polarization 
    if d['spin']=='1/2'   and d['kind']=='1': pol,cpol='U','UB'
    elif d['spin']=='1/2' and d['kind']=='-1': pol,cpol='V','VB'
    elif d['spin']=='1': pol,cpol='eps','epsB'
    elif d['spin']=='0': 
      pol,cpol='1','1' 
      if particle=='gho':
        if d['ifield'] not in self.LImap:  self.LImap[d['ifield']]=self.LI.pop(0)
        o_=str(self.LImap[d['ifield']])
        pol='X(j%s)'%(o_)
        o_=str(self.LImap[d['ifield']]+1000)
        cpol='XB(j%s)'%(o_)
  
    # add bar:
    if d['spin']=='1/2' or d['spin']=='1': 
      if d['kind']=='1' and 'q' in d['mom']: 
        pol+='B'
        cpol=cpol.replace('B','')
      elif d['kind']=='-1' and 'p' in d['mom']: 
        pol+='B'
        cpol=cpol.replace('B','')
  
    # define spin projection
    proj=d['mom'].replace('p','s').replace('q','r')
  
    # add particle label
    if d['spin']=='1/2':
      chain='[%s%s]'%(particle,d['ifield'])
      pol ='%s(%s,%s,%s,%s,%s)'%(pol,chain,d['particle'],d['mom'],proj,d['mass'])
      cpol='%s(%s,%s,%s,%s,%s)'%(cpol,chain,d['particle'],d['mom'],proj,d['mass'])

    elif d['spin']=='1':
      if d['ifield'] not in self.LImap:  self.LImap[d['ifield']]=self.LI.pop(0)
      if d['particle']!='G':
        o_=str(self.LImap[d['ifield']])
        pol='%s(%s,%s,%s,%s,o%s)'%(pol,d['particle'],d['mom'],proj,d['mass'],o_)
        o_=str(self.LImap[d['ifield']]+1000)
        cpol='%s(%s,%s,%s,%s,o%s)'%(cpol,d['particle'],d['mom'],proj,d['mass'],o_)
      elif d['particle']=='G':
        o_=str(self.LImap[d['ifield']])
        pol='%s(%s,%s,%s,%s,o%s,j%s)'%(pol,d['particle'],d['mom'],proj,d['mass'],o_,o_)
        o_=str(self.LImap[d['ifield']]+1000)
        cpol='%s(%s,%s,%s,%s,o%s,j%s)'%(cpol,d['particle'],d['mom'],proj,d['mass'],o_,o_)


    # register expression
    self.AMPS[iamp][ipart]['expression']=pol
    self.AMPS[iamp][ipart]['cexpression']=cpol

    # fill info of external states
    self.external_states[d['mom']]={'mass':d['mass'],'spin':proj}

  def get_vrtx(self,iamp,ipart):
    d=self.AMPS[iamp][ipart]
    # get couplings
    if   d['type']=='QEDffvl':   coupling='gl^%s'%d['gpow']
    if   d['type']=='QEDffvu':   coupling='gu^%s'%d['gpow']
    elif d['type']=='QEDffvd':   coupling='gd^%s'%d['gpow']
    elif 'QCD' in d['type']:     coupling='gs^%s'%d['gpow']
    else:   
      print 'ERR: cannot interpret vertex' 
      print 'type=',d['type']
      sys.exit() 

    # select from vertex type
    if 'QEDffv' in d['type']:
      chain='['
      for k in d.keys():
        if 'ifield' not in k: continue
        field=k.split('_')[1]
        if d['spin_'+field]!='1/2': continue
        if field.endswith('b'): field=field[:-1] 
        chain+='%s%s,'%(field,d[k])
      chain=chain[:-1]+']'
      if d['ifield_A'] not in self.LImap:  self.LImap[d['ifield_A']]=self.LI.pop(0)
      idx=str(self.LImap[d['ifield_A']])
      vertex ='(-i_)*%s*gI_(%s,o%s)'%(coupling,chain,idx)
      idx=str(self.LImap[d['ifield_A']]+1000)
      cvertex='(i_)*%s*gI_(%s,o%s)'%(coupling,chain,idx)

    elif d['type']=='QCDffv':
      chain='['
      for k in d.keys():
        if 'ifield' not in k: continue
        field=k.split('_')[1]
        if d['spin_'+field]!='1/2': continue
        if field.endswith('b'): field=field[:-1] 
        chain+='%s%s,'%(field,d[k])
      chain=chain[:-1]+']'
      cchain=chain.replace('[','').replace(']','').split(',')
      cchain='[%s,%s]'%(cchain[1],cchain[0])
      if d['ifield_G'] not in self.LImap:  self.LImap[d['ifield_G']]=self.LI.pop(0)
      idx=str(self.LImap[d['ifield_G']])
      vertex ='(-i_)*%s*gI_(%s,o%s)*T(j%s,%s)'%(coupling,chain,idx,idx,chain)
      idx=str(self.LImap[d['ifield_G']]+1000)
      cvertex='(i_)*%s*gI_(%s,o%s)*T(j%s,%s)'%(coupling,chain,idx,idx,cchain)

    elif d['type']=='QCDv3':
      data=d['raw'][1:-1].split(',')
      legs=[]
      legs.append({'mom':data[4].split('=')[1] ,'ifield':data[5].split('=')[1]})
      legs.append({'mom':data[7].split('=')[1] ,'ifield':data[8].split('=')[1]})
      legs.append({'mom':data[10].split('=')[1],'ifield':data[11].split('=')[1]})
      legs=legs[::-1]
      for l in legs:
        xx=l['mom'].replace('-',',-').replace('+',',+').split(',')
        xx=[x for x in xx if x!='']
        l['mom']=''
        for x in xx: l['mom']+=x+'(@)' 
        l['mom']='('+l['mom']+')'
      l0=legs[0]
      l1=legs[1]
      l2=legs[2]
      fabc=''
      for l in legs:
        if l['ifield'] not in self.LImap: self.LImap[l['ifield']]=self.LI.pop(0)
        fabc+='j'+str(self.LImap[l['ifield']])+','
      fabc='f('+fabc[:-1]+')'
      vertex='(-1)*'+coupling+'*'+fabc+'*('
      vertex+='d_(o%s,o%s)'%(self.LImap[l0['ifield']],self.LImap[l1['ifield']])
      vertex+='*(%s-%s)'%(l0['mom'].replace('@','o'+str(self.LImap[l2['ifield']])),l1['mom'].replace('@','o'+str(self.LImap[l2['ifield']])))
      vertex+='+'
      vertex+='d_(o%s,o%s)'%(self.LImap[l1['ifield']],self.LImap[l2['ifield']])
      vertex+='*(%s-%s)'%(l1['mom'].replace('@','o'+str(self.LImap[l0['ifield']])),l2['mom'].replace('@','o'+str(self.LImap[l0['ifield']])))
      vertex+='+'
      vertex+='d_(o%s,o%s)'%(self.LImap[l2['ifield']],self.LImap[l0['ifield']])
      vertex+='*(%s-%s)'%(l2['mom'].replace('@','o'+str(self.LImap[l1['ifield']])),l0['mom'].replace('@','o'+str(self.LImap[l1['ifield']])))
      vertex+=')' 
      fabc=''
      for l in legs:
        fabc+='j'+str(self.LImap[l['ifield']]+1000)+','
      fabc='f('+fabc[:-1]+')'
      cvertex='(-1)*'+coupling+'*'+fabc+'*('
      cvertex+='d_(o%s,o%s)'%(self.LImap[l0['ifield']]+1000,self.LImap[l1['ifield']]+1000)
      cvertex+='*(%s-%s)'%(l0['mom'].replace('@','o'+str(self.LImap[l2['ifield']]+1000)),l1['mom'].replace('@','o'+str(self.LImap[l2['ifield']]+1000)))
      cvertex+='+'
      cvertex+='d_(o%s,o%s)'%(self.LImap[l1['ifield']]+1000,self.LImap[l2['ifield']]+1000)
      cvertex+='*(%s-%s)'%(l1['mom'].replace('@','o'+str(self.LImap[l0['ifield']]+1000)),l2['mom'].replace('@','o'+str(self.LImap[l0['ifield']]+1000)))
      cvertex+='+'
      cvertex+='d_(o%s,o%s)'%(self.LImap[l2['ifield']]+1000,self.LImap[l0['ifield']]+1000)
      cvertex+='*(%s-%s)'%(l2['mom'].replace('@','o'+str(self.LImap[l1['ifield']]+1000)),l0['mom'].replace('@','o'+str(self.LImap[l1['ifield']]+1000)))
      cvertex+=')' 

    elif d['type']=='QCDghost':
      data=d['raw'][1:-1].split(',')
      legs=[]
      legs.append({'field':data[4].split('=')[0] ,'mom':data[4].split('=')[1] ,'ifield':data[5].split('=')[1]})
      legs.append({'field':data[7].split('=')[0] ,'mom':data[7].split('=')[1] ,'ifield':data[8].split('=')[1]})
      legs.append({'field':data[10].split('=')[0],'mom':data[10].split('=')[1],'ifield':data[11].split('=')[1]})
      legs=legs[::-1]
      for l in legs:
        xx=l['mom'].replace('-',',-').replace('+',',+').split(',')
        xx=[x for x in xx if x!='']
        l['mom']=''
        for x in xx: l['mom']+=x+'(@)' 
        l['mom']='('+l['mom']+')'
      fabc=''
      for l in legs:
        if l['ifield'] not in self.LImap: self.LImap[l['ifield']]=self.LI.pop(0)
        fabc+='j'+str(self.LImap[l['ifield']])+','
      fabc='f('+fabc[:-1]+')'
      for l in legs:
        if l['field']=='mom_gho': l0=l
        if l['field']=='mom_G':   l1=l
      vertex='(-1)*'+coupling+'*'+fabc+'*('
      vertex+='(%s)'%(l0['mom'].replace('@','o'+str(self.LImap[l1['ifield']])))
      vertex+=')' 
      fabc=''
      for l in legs:
        fabc+='j'+str(self.LImap[l['ifield']]+1000)+','
      fabc='f('+fabc[:-1]+')'
      cvertex=coupling+'*'+fabc+'*('
      cvertex+='(%s)'%(l0['mom'].replace('@','o'+str(self.LImap[l1['ifield']]+1000)))
      cvertex+=')' 

    # register expression
    self.AMPS[iamp][ipart]['expression']=vertex
    self.AMPS[iamp][ipart]['cexpression']=cvertex
  
  def get_prop(self,iamp,ipart):
    d=self.AMPS[iamp][ipart]
    if d['spin']=='1/2':
      field=d['particle']
      if field.endswith('b'): field=field[:-1] 
      chain='[%s%s,%s%s]'%(field,d['ifield'],field,d['idfield'])
      moms=[m for m in d['mom'].replace('-',',-').replace('+',',+').split(',') if m!='']
      smom=''
      for m in moms: smom+='g_(---,%s)+'%m
      smom=smom.replace('---',chain)
      num ='i_*(%sg_(%s)*%s)'%(smom,chain,d['mass'])
      cnum='(-i_)*(%sg_(%s)*%s)'%(smom,chain,d['mass'])
    elif d['spin']=='1':

      if d['ifield'] not in self.LImap:  self.LImap[d['ifield']]=self.LI.pop(0)
      if d['idfield'] not in self.LImap:  self.LImap[d['idfield']]=self.LI.pop(0)
      ini=self.LImap[d['ifield']]
      fin=self.LImap[d['idfield']]

      if d['particle']=='G':
        num ='(-i_)*d_(o%d,o%d)*d_(j%d,j%d)'%(ini,fin,ini,fin)
        cnum='i_*d_(o%d,o%d)*d_(j%d,j%d)'%(ini+1000,fin+1000,ini+1000,fin+1000)
      else:
        num='(-i_)*d_(o%s,o%s)'%(ini,fin)
        cnum='i_*d_(o%s,o%s)'%(ini+1000,fin+1000)
    elif d['spin']=='0':
      num ='i_'
      cnum='(-i_)'
    den='den(%s,%s)'%(d['mom'],d['mass'])
    prop='%s*%s'%(num,den)
    cprop='%s*%s'%(cnum,den)
    # register expression
    self.AMPS[iamp][ipart]['expression']=prop
    self.AMPS[iamp][ipart]['cexpression']=cprop

  def get_ordering(self,iamp,verb=False):
    AMPS=self.AMPS

    # construct array of data such that data=[[ipart,[eb-1,e2]],...]
    data=[]
    for ipart in AMPS[iamp].keys():
      exp=AMPS[iamp][ipart]['expression']
      match=re.search(r'\[(.*?)\]',exp)
      if match!=None:
        braket=match.group().replace('[','').replace(']','').split(',')
        data.append([ipart,braket])

    if verb:
      print '\nspinor chains:\n'
      for d in data: print d

    # if no spinors then return original iparts keys
    if len(data)==0:
      return None,AMPS[iamp].keys()

    # get spinor chains, get external spinor if available
    chains=[]
    for d in data:
      ipart,braket=d              
      if len(braket)==1: break
    node=braket[0]            
    if verb: print '\nstarting node:',node,ipart,'\n'
    chain=[ipart]                 # create a chain and register the first element from data
    if verb: print 'dsize=',len(data),'data=',data
    data=[x for x in data if any([x[0]==y for y in chain])==False]  # update data 

    cnt=0                 # conter for ?????
    while 1:

      flag=False  

      for d in data:
        ipart,braket=d
        if any([node==x for x in braket]):                 # checkout if node is in braket
          if any([ipart==x for x in chain])==False:        # check if ipart has not been registerd in the current chain
            chain.append(ipart)                            # ipart has not been registered. We add it in the current chain
            if len(braket)==1: node=braket[0]              # get the next node if [_] and replace the current node
            else: node=[x for x in braket if x!=node][0]   # or get the next node if [_,_] and replace the current node (do not pick the same node)
            flag=True                                      # A node was found and we raise the flag
      
      # update data by removing the registered iparts
      if verb: print 'dsize=',len(data),'data=',data
      data=[x for x in data if any([x[0]==y for y in chain])==False] 

      #if verb: print data
      if flag==False:                                              # means we got to the end of a spinor chain
        symbol=AMPS[iamp][chain[0]]['expression'].split('(')[0]
        if symbol.endswith('B')==False: chain=chain[::-1]          # reverse the order to get correct Dirac tensor ordering
        chains.append(chain[:])                                    # add the chain to the list of  chains
        if len(data)==0:                                           # means there are no more data so we are done!!
          break
        else:                                                      # or we need to start a new spinor chain search from another node
          for d in data:
            ipart,braket=d
            if len(braket)==1: break
          node=braket[0]
          if verb: print '\nstarting node:',node,ipart,'\n'
          chain=[ipart]
          if verb: print 'dsize=',len(data),'data=',data
          data=[x for x in data if any([x[0]==y for y in chain])==False]  # update data 

    # construct ordered iparts
    ordered_parts=[]
    for chain in chains: ordered_parts.extend(chain)
    nonchains=sorted(set(AMPS[iamp].keys())-set(ordered_parts))
    return chains,nonchains

  def get_AMPS(self,verb=False):
    if verb:
      print '='*80
      print 'CONSTRUCT AMPS'
    AMPS=self.AMPS

    for iamp in AMPS.keys():
      if verb: print '-'*80

      # construct feyman rules
      for ipart in AMPS[iamp].keys():
        part=AMPS[iamp][ipart]['part']
        if part=='scalar': self.get_scalar(iamp,ipart)
        elif part=='pol':  self.get_pol(iamp,ipart)
        elif part=='vrtx': self.get_vrtx(iamp,ipart)
        elif part=='prop': self.get_prop(iamp,ipart)

      # construct spinor chains
      chains,nonchains=self.get_ordering(iamp,verb=verb)
      AMPS[iamp]['chains']=chains
      AMPS[iamp]['nonchains']=nonchains

      parts=[]
      cparts=[]
      if chains!=None:
        if verb:
          print '\nspinor chains\n'
          for chain in chains: print chain
        for chain in chains: parts.extend(chain)
        for chain in chains: cparts.extend(chain[::-1])
      parts.extend(nonchains)
      cparts.extend(nonchains)

      # print outs
      if verb:
        print '\nAMP: iamp=',iamp
        for part in parts:
          if verb: print '%10s = %10s'%(part,AMPS[iamp][part]['expression'])
        print 'AMP(dagger):'
        for part in cparts:
          if verb: print '%10s = %10s'%(part,AMPS[iamp][part]['cexpression'])

  def get_SQAMPS(self,verb=False):
    if verb:
      print '='*80
      print 'CONSTRUCT SQAMPS'

    AMPS=self.AMPS
    SQAMPS={}
    for iamp in AMPS.keys():
      for icamp in AMPS.keys():
        if verb: print '-'*80
        key='[%s*c%s]'%(iamp,icamp) 
        if verb: print key

        sqamp=''

        ####################################
        # get scalars
        ####################################
        nonchains=AMPS[iamp]['nonchains'] 
        cnonchains=AMPS[icamp]['nonchains']
        for i in range(len(nonchains)):  sqamp+=AMPS[iamp][nonchains[i]]['expression']+'*'
        for i in range(len(cnonchains)): sqamp+=AMPS[icamp][cnonchains[i]]['cexpression']+'*'
        sqamp=sqamp[:-1]

        ####################################
        # get spinor groups
        ####################################

        #chains=AMPS[iamp]['chains']
        #cchains=AMPS[icamp]['chains']
        
        groups=[]

        if verb:
          print '\nspinor chains'
          print '---------------'
        if AMPS[iamp]['chains']!=None:
          for chain in AMPS[iamp]['chains']:
            subgroup=[] 
            if verb:print 
            for part in chain:
              if verb:print AMPS[iamp][part]['expression']
              subgroup.append(AMPS[iamp][part]['expression'])
            groups.append(subgroup)

        if verb:
          print '\nspinor cchains'
          print '---------------'
        if AMPS[icamp]['chains']!=None:
          for chain in AMPS[icamp]['chains']:
            subgroup=[] 
            if verb:print 
            for part in chain[::-1]:
              if verb:print AMPS[icamp][part]['cexpression']
              subgroup.append(AMPS[icamp][part]['cexpression'])
            groups.append(subgroup)

        # organize subgropus according to tail-head matching

        if len(groups)!=0:
          ogroups=[]
          ogroup=[groups[0]]
          groups.pop(0)
          while 1:
            found=False

            for i in range(len(groups)):
              group=groups[i]
              
              tail=ogroup[-1][-1]
              match=re.search(r'\[(.*?)\]',tail)
              tail=tail.replace(match.group()+',','')

              head=group[0].replace('B(','(')
              match=re.search(r'\[(.*?)\]',head)
              head=head.replace(match.group()+',','')
          
              if tail==head: 
                ogroup.append(group)
                found=True
                break

            if found:
              groups.pop(i) 
            else:
              ogroups.append(ogroup)
              ogroup=[groups[0]]
              groups.pop(0)

            if len(groups)==0: 
              ogroups.append(ogroup)
              break

          # something 
          cnt=0
          for ogroup in ogroups:
            if verb:
              print '\nordered groups\n' 
              for g in ogroup: 
                print g

            # construct FORM formated spinor chain
            chain=''
            cnt+=1
            for g in ogroup: 
              for x in g: 
                match=re.search(r'\[(.*?)\]',x)
                x=x.replace(match.group(),'sc%d'%cnt)
                chain+=x+'*'
            chain=chain[:-1]
            schain=chain.split('*')
            chain_=[schain[-1]]
            chain_.extend(schain[:-1])
            chain=''
            for x in chain_:chain+=x+'*'
            chain=chain[:-1]

            # add color chains
            color_factors=[x for x in chain.split('*') if x.startswith('T(')]
            if len(color_factors)>0:
              for cf in color_factors: chain=chain.replace('*'+cf,'')
              colortrace='Tr('
              for cf in color_factors:
                colortrace+=cf.replace('T(','').replace(')','').split(',')[0]+','
              colortrace=colortrace[:-1]+')'
              chain+='*'+colortrace      
 
            sqamp+='*'+chain

          #print 
          #print 'T colors ---> ',color_factors
          #print 
          #print '---> chain'
          #print chain

        if verb:
          print '\nsquared amplitude'
          print '-------------------'
          print sqamp

        SQAMPS[key]=sqamp

    self.SQAMPS=SQAMPS

if __name__=='__main__':

  conf={}
  conf['path2qgraf']='./'
  conf['label']='bornA'
  conf['style']='array.sty'
  conf['model']='SM'
  conf['in']='A[p1],ub[p2]'
  conf['out']='ub[q1],G[q2]'
  conf['loops']=0;
  conf['loop_momentum']=''
  conf['options']='notadpole'
  qgraf=QGRAF(conf)

 

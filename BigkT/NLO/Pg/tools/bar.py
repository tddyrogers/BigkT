#!/usr/bin/env python
import sys,os

class BAR(object):

  def __init__(self,msg,size):
    self.msg=msg
    self.size=size
    self.cnt=0

  def next(self):
    sys.stdout.write('\r')
    percentage=int(self.cnt/float(self.size)*100)
    #sys.stdout.write('%s [%d%%]' % (self.msg,percentage))
    sys.stdout.write('%s [%d/%d]' % (self.msg,self.cnt,self.size))
    sys.stdout.flush()
    self.cnt+=1

  def finish(self):
    self.next()
    print 



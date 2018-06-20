class VEC:

  def __init__(self,symbol,values=[None]):
    Vp=sp.Symbol(r'%s_+'%symbol)
    Vm=sp.Symbol(r'%s_-'%symbol)
    V1=sp.Symbol(r'%s_1'%symbol)
    V2=sp.Symbol(r'%s_2'%symbol)
    self.V=np.array([Vp,Vm,V1,V2])
    if len(values)!=1: self.set_values(values)
    self._perp=sp.Symbol(r'%s_{\perp}'%symbol)
    
  def __mul__(self,other):
    if isinstance(other, VEC):
      return self.V[0]*other.V[1]+self.V[1]*other.V[0]-self.V[2]*other.V[2]-self.V[3]*other.V[3] 
    else:
      new=copy.copy(self)
      new.V*=other
      return new

  def __rmul__(self,other):
    if isinstance(other, VEC):
      return self.V[0]*other.V[1]+self.V[1]*other.V[0]-self.V[2]*other.V[2]-self.V[3]*other.V[3] 
    else:
      new=VEC('A')
      new.V=other*self.V
      return new
    
  def __add__(self,other):
    if isinstance(other, VEC):
      new=VEC('A')
      new.V=self.V+other.V
      return new 
    else:
      print 'ERR: cannot add'

  def __sub__(self,other):
    if isinstance(other, VEC):
      new=VEC('A')
      new.V=self.V-other.V
      return new 
    else:
      print 'ERR: cannot add'
    
  def set_values(self,values):
    for idx in range(4): 
      if values[idx]!=None:      
        self.V[idx]=values[idx]

  def plus(self):
    return self.V[0]

  def minus(self):
    return self.V[1]

  def perp(self):
    return self._perp

  def perp2(self):
    return self.V[2]

  def perp3(self):
    return self.V[3]

  def show(self):
    return self.V[0],self.V[1],self.V[2],self.V[3]

  def mass(self):
    return 2*self.V[0]*self.V[1]-self._perp**2

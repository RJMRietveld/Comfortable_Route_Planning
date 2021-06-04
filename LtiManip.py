from __future__ import division

from scipy.signal.ltisys import TransferFunction as TransFun
from numpy import polymul,polyadd

class ltimul(TransFun):
    def __neg__(self):
        return ltimul(-self.num,self.den)

    def __floordiv__(self,other):
        # can't make sense of integer division right now
        return NotImplemented

    def __mul__(self,other):
        if type(other) in [int, float]:
            return ltimul(self.num*other,self.den)
        elif type(other) in [TransFun, ltimul]:
            numer = polymul(self.num,other.num)
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    def __truediv__(self,other):
        if type(other) in [int, float]:
            return ltimul(self.num,self.den*other)
        if type(other) in [TransFun, ltimul]:
            numer = polymul(self.num,other.den)
            denom = polymul(self.den,other.num)
            return ltimul(numer,denom)

    def __rtruediv__(self,other):
        if type(other) in [int, float]:
            return ltimul(other*self.den,self.num)
        if type(other) in [TransFun, ltimul]:
            numer = polymul(self.den,other.num)
            denom = polymul(self.num,other.den)
            return ltimul(numer,denom)

    def __add__(self,other):
        if type(other) in [int, float]:
            return ltimul(polyadd(self.num,self.den*other),self.den)
        if type(other) in [TransFun, type(self)]:
            numer = polyadd(polymul(self.num,other.den),polymul(self.den,other.num))
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    def __sub__(self,other):
        if type(other) in [int, float]:
            return ltimul(polyadd(self.num,-self.den*other),self.den)
        if type(other) in [TransFun, type(self)]:
            numer = polyadd(polymul(self.num,other.den),-polymul(self.den,other.num))
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    def __rsub__(self,other):
        if type(other) in [int, float]:
            return ltimul(polyadd(-self.num,self.den*other),self.den)
        if type(other) in [TransFun, type(self)]:
            numer = polyadd(polymul(other.num,self.den),-polymul(other.den,self.num))
            denom = polymul(self.den,other.den)
            return ltimul(numer,denom)

    # sheer laziness: symmetric behaviour for commutative operators
    __rmul__ = __mul__
    __radd__ = __add__
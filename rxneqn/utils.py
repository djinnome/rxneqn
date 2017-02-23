from fractions import gcd
from fractions import Fraction
import numpy as np
from functools import reduce
import warnings
warnings.filterwarnings('ignore') # Otherwise python3.5 says that the Fraction gcd is deprecated.

def LCM(fractions):
    
    try: 
        x = reduce(lambda x, y: (x*y)/gcd(x,y), [n for n in fractions if n !=0], 1)
    except DeprecationWarning:
        return x
    return x
def GCD( numbers ):
    return reduce( lambda x, y: gcd(Fraction(x),Fraction(y)), numbers, 1) 

def LCD( fractions ):
    return LCM([Fraction(f).denominator for f in fractions])

import numpy as np

def trg( x, b, a ):
    return 1 / ( 1 + np.exp( ( b - x ) / a ) )

def trigger(x, b, a ):
   p = trg( x, b, a )
   rand = np.random.rand( )
   if( rand > p ): return True
   else: return False
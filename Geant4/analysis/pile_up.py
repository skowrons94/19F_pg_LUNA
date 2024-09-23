import numpy as np

def pile_up( prob ):
    rand = np.random.rand( ) 
    if( rand > prob ): return False
    else: return True
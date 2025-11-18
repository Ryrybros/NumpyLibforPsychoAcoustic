from jsonHandler import jsonHandler
import numpy as np
class Interpolator:


    def interp1(x : float, xp : np.array , fp : np.array):
        
        return np.interp(x, xp, fp)
    
    


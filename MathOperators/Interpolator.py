
import numpy as np
class Interpolator:


    def interp1(x : float, xp : np.array , fp : np.array):
        #check if this is the right interpolation
        return np.interp(x, xp, fp)
    


if __name__ == '__main__':
    x = [2,3,4]
    y = [4,6,8]
    a = 2.5
    assert(Interpolator.interp1( a,x,y) == 5)


    
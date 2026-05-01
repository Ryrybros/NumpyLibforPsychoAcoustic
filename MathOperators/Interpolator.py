
import numpy as np
import scipy.interpolate as I
class Interpolator:

    #This class is just a wrapper for numpy and scipy interpolators

    def interp1(x : np.array, y : np.array , ax : np.array,pchip = False):
        #Need to check if this is the right interpolation
        
        assert(len(x) == len(y)  )

        if(pchip):
            pch = I.PchipInterpolator(x= x,y = y)
            return pch(ax)
        
        return np.interp(x, xp, fp)
    

    


if __name__ == '__main__':
    n = 100
    x = [(np.pi/n)*i for i in range(n)]
    y = np.sin(x)
    a = np.array([2.5,3,6])
    inter = Interpolator.interp1( x,x,y)
    # assert(inter[0] == 5)
    import matplotlib.pyplot as plt
    # plt.plot(x,inter)
    # plt.show()
    #print(inter)


    
import numpy as np
from filetools import jsonHandler

import matplotlib.pyplot as plt
import scipy.interpolate as I

Json = jsonHandler.jsonHandler

class Interpolator:

    #This class is just a wrapper for numpy and scipy interpolators

    def interp1(x : np.array, y : np.array , ax : np.array,pchip = False):
        #Need to check if this is the right interpolation
        
        assert(len(x) == len(y)  )

        if(pchip):
            pch = I.PchipInterpolator(x= x,y = y)
            out_of_bounds = (ax < np.min(x)) | (ax > np.max(x))
            res = pch(ax)
            res[out_of_bounds] = 0 # Nan is probably better
            #Python extrapolates values where it should not, the nan follow matlab
            
            return res

        return np.interp(x, xp, fp)

Interp = Interpolator

class dataPreparator:


    def OuterMiddle(data : dict ,fVec : np.array, model : str, free : bool ):
                
        class ReturnedData:

            def __init__(self):
                #OuterMiddleEar
                self.empty = True            

        # transfer function of the outer ear
        try:
            dataModel = f"OuterMiddleEar{model}"
        except:
            print("Wrong model parameter , model 1997 was chosen by default.")
            dataModel = "OuterMiddleEar1997"
        if(free):
            regim = "Free"
        else:
            regim = "Diffuse"

        fOuter = np.array(data[dataModel]["fOuter"])
        tfOuter = np.array(data[dataModel][f"tfOuter{regim}"])

        tfOuterInterp = Interp.interp1(fOuter, tfOuter, fVec, pchip= True) 
        
        fMiddle = np.array(data[dataModel]["fMiddle"])
        tfMiddle= np.array(data[dataModel]["tfMiddle"])
        # print(tfMiddle)

        tfMiddleInterp = Interp.interp1(fMiddle, tfMiddle, fVec,pchip= True) 
        

        dat = ReturnedData()

        dat.tfOuterMiddle = tfMiddleInterp + tfOuterInterp
        dat.tfOuter = tfOuter
        dat.tfMiddle = tfMiddle
        dat.fOuter = fOuter
        dat.fMiddle = fMiddle

        return dat
    
    def SpecLoudness(data : dict, fVec : np.array):

        
        class ReturnedData:

            def __init__(self):
                #OuterMiddleEar
                self.empty = True   
                
        dat = ReturnedData()
        fRef = data["fRef"]
        # print(fRef)
        tQ = data["tQ"]
        
        dat.tQ = Interp.interp1(fRef, tQ, fVec,pchip=True)
        
        
        # open file
        with open('gfg.txt', 'w+') as f:
            
            # write elements of list
            for items in dat.tQ:
                f.write('%s\n' %items)
            
            print("File written successfully")


        # close the file
        f.close()
        # dat.tQ = np.interp(fVec, fRef,tQ)
        # dat.tQ[ : len(dat.tQ) - 1 ] = dat.tQ[ 1 : ] 
        # dat.tQ[:10] *= 0
        dat.tQ500 = tQ[11]
        dat.g = dat.tQ500-dat.tQ    # low level gain in cochlea amplifier
        
        # linearization parameter a
        g = np.array(data["g"])
        
        
        a = np.array(data["a"])
        
        dat.a =  - Interp.interp1(-g, -a, dat.g, pchip=True )  #Again, none pchip  !! we use - g and- a because it is necessary to have an increasing array for the second argumetn
        
        #  compressive exponent alpha
        g = np.array(data["gCompression"])

        alpha = np.array(data["alpha"])


        dat.alpha =  Interp.interp1(g, alpha, dat.g, pchip= True)

        dat.c = data["c"]

        return dat


    

if __name__ == '__main__':

    
    g = dataPreparator()
        # print("OuterMiddle treatment is : "  , g.OuterMiddle(np.array([0.5,7,1000]),"2007",False))
    #print("SpecLoud is :" , g.OuterMiddle([1,2],"1997",True).tfOuterMiddle )
    # linspace = np.array([(1/(len(y) - 1 ))*i for i in range(len(y))])
    # plt.plot(linspace, y)
    



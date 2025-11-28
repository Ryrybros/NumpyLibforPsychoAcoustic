import numpy as np
from filetools import jsonHandler
from MathOperators import Interpolator
import matplotlib.pyplot as plt


Json = jsonHandler.jsonHandler
Interp = Interpolator.Interpolator

class ReturnedData:

    def __init__(self):
        #OuterMiddleEar
        tfOuterMiddle = None
        tfOuter = None
        tfMiddle = None
        fOuter = None
        fMiddle = None

        #SpecLoud
        tQ = None
        tQ500 = None
        g = None
        a = None
        alpha = None
        c = None


class dataPreparator:


    def OuterMiddle(data : dict ,fVec : np.array, model : str, free : bool ):

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

        tfOuterInterp = Interp.interp1(fVec,fOuter, tfOuter) #None Pchip interpolation, this is linear, pchip is cubic
        
        fMiddle = np.array(data[dataModel]["fMiddle"])
        tfMiddle= np.array(data[dataModel]["tfMiddle"])
        # print(tfMiddle)

        tfMiddleInterp = Interp.interp1(fVec,fMiddle, tfMiddle) #None Pchip interpolation, this is linear, pchip is cubic
        

        dat = ReturnedData()

        dat.tfOuterMiddle = tfMiddleInterp + tfOuterInterp
        dat.tfOuter = tfOuter
        dat.tfMiddle = tfMiddle
        dat.fOuter = fOuter
        dat.fMiddle = fMiddle

        y = tfOuterInterp
        # linspace = np.array([(1/(len(y) - 1 ))*i for i in range(len(y))])
        # plt.plot(fVec, fVec)
        # # plt.plot(fMiddle,tfMiddle)
        # plt.show()
        
        return dat
    
    def SpecLoudness(data : dict, fVec : np.array):

        dat = ReturnedData()
        fRef = data["fRef"]
        # print(fRef)
        tQ = data["tQ"]

        dat.tQ = Interp.interp1(fVec, fRef,tQ)
        dat.tQ500 = tQ[11]
        
        dat.g = dat.tQ500-dat.tQ    # low level gain in cochlea amplifier

        # linearization parameter a
        g = data["g"]
        # print(g)

        a = data["a"]

        dat.a = Interp.interp1(dat.g, g, a) #Again, none pchip

        #  compressive exponent alpha
        g = np.array(data["gCompression"])

        alpha = data["alpha"]

        dat.alpha = Interp.interp1(dat.g,g, alpha)

        dat.c = data["c"]

        return dat


    

if __name__ == '__main__':

    
    g = dataPreparator()
        # print("OuterMiddle treatment is : "  , g.OuterMiddle(np.array([0.5,7,1000]),"2007",False))
    #print("SpecLoud is :" , g.OuterMiddle([1,2],"1997",True).tfOuterMiddle )
    # linspace = np.array([(1/(len(y) - 1 ))*i for i in range(len(y))])
    # plt.plot(linspace, y)
    



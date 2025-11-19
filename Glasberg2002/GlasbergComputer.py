import numpy as np
from filetools import jsonHandler
from MathOperators import Interpolator

Json = jsonHandler.jsonHandler
Interp = Interpolator.Interpolator

class ReturnedData:

    def __init__(self):
        #OuterMiddleEar
        self.tfOuterMiddle = None
        self.tfOuter = None
        self.tfMiddle = None
        self.fOuter = None
        self.fMiddle = None

#SpecLoud
        self.tQ = None

class GlasbergComputer:

    def __init__(self):
        self.data = Json.readJson("data/Glasberg2002.json")

    def OuterMiddle(self,fVec : np.array, model : str, free : bool ):

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

        fOuter = np.array(self.data[dataModel]["fOuter"])
        tfOuter = np.array(self.data[dataModel][f"tfOuter{regim}"])

        tfOuterInterp = Interp.interp1(fVec,fOuter, tfOuter) #None Pchip interpolation, this is linear, pchip is cubic

        fMiddle = np.array(self.data[dataModel]["fMiddle"])
        tfMiddle= np.array(self.data[dataModel]["tfMiddle"])
        # print(tfMiddle)

        tfMiddleInterp = Interp.interp1(fVec,fMiddle, tfMiddle) #None Pchip interpolation, this is linear, pchip is cubic
        

        dat = ReturnedData()

        dat.tfOuterMiddle = tfMiddleInterp + tfOuterInterp
        dat.tfOuter = tfOuter
        dat.tfMiddle
        dat.fOuter
        dat.fMiddle
        
        return dat
    
    def SpecLoudness(self , fVec : np.array):

        dat = ReturnedData()
        fRef = self.data["fRef"]
        # print(fRef)
        tQ = self.data["tQ"]

        dat.tQ = Interp.interp1(fVec, fRef,tQ)
        dat.tQ500 = tQ[11]
        
        dat.g = dat.tQ500-dat.tQ    # low level gain in cochlea amplifier

        # linearization parameter a
        g = self.data["g"]
        # print(g)

        a = self.data["a"]

        dat.a = Interp.interp1(g, a, dat.g) #Again, none pchip

        #  compressive exponent alpha
        g = np.array(self.data["gCompression"])

        alpha = self.data["alpha"]

        # dat.alpha = Interp.interp1(g, alpha, dat.g)

        dat.c = self.data["c"]

        return dat
        # data.alpha = , 'pchip');

        # data.c = 0.046871; % constant to get loudness scale to sone
    
    
        
        

   
    #     data.tfOuter = tfOuter;
    #     data.tfMiddle = tfMiddle;
    #     data.fOuter = fOuter;
    #     data.fMiddle = fMiddle;
    # end
    

    

if __name__ == '__main__':
    
    g = GlasbergComputer()
    # print("OuterMiddle treatment is : "  , g.OuterMiddle(np.array([0.5,7,1000]),"2007",False))
    print("SpecLoud is :" , g.SpecLoudness([1,2]) )



from filetools.jsonHandler import jsonHandler
from Glasberg2002.GlasbergComputer import GlasbergComputer
from Glasberg2002.lfatargrepoducer import KV

from MathOperators.Signals import sineMaker

import numpy as np
import matplotlib.pyplot as plt

class FilterComputer :
    def FIR(kv : KV, model : str, free : bool):
        #General computation of FIR filter independantly of chosen model or hypothesis.
        data = jsonHandler.readJson("data/Glasberg2002.json") 
        glasberg = GlasbergComputer()   
        fVec = np.linspace(kv.flow,kv.fhigh, abs(kv.flow - kv.fhigh) )
        glasbData = glasberg.OuterMiddle(fVec, model, free )
        
        tfLinear = 10**(glasbData.tfOuterMiddle/10)

        print("tfLinear is : ")
        print(tfLinear)
        tfLinear[len(tfLinear) - 1] = 0.0 #This is imposed by the fir2 filter function because of the following error  : 
        #A Type II filter must have zero gain at the Nyquist frequency.

        linspace = np.array([(1/(len(fVec) - 1 ))*i for i in range(len(fVec))])

        plt.plot(linspace, tfLinear)
        
        
        # print(linspace)
        outerMiddleFilter = FilterComputer.fir2(kv.order, linspace, tfLinear)
        outerMiddleFilter = outerMiddleFilter[ : int(len(outerMiddleFilter)/2)] #seems to be more corresponding to the matlab behavior
        # print(outerMiddleFilter)

        plt.plot(np.linspace(0,1,len(outerMiddleFilter)), outerMiddleFilter)

        plt.show()


        

    def fir2(order : float, domain : np.array, interpolatedValues : np.array):
        from scipy import signal
        return signal.firwin2(order, domain, interpolatedValues )
    

        


# FIR(np.array([2,0,3,5]) , "1997" , True)

if __name__ == '__main__':
    kv = KV()
    kv.setStandard()
    f = FilterComputer()
    y = sineMaker.makeSine(10,1,40,20000)
    fv = np.linspace(kv.flow,20000)
    FilterComputer.FIR(kv ,'1997', True)
    
    
    

    




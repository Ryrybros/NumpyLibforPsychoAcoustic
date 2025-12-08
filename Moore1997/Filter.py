
from filetools.jsonHandler import jsonHandler
from Moore1997._dataPreparator import dataPreparator
import Moore1997.model as test
from scipy import signal
from MathOperators.Signals import sineMaker
    
    
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

class FilterComputer :

    def __init__(self,kv : dict , model : str, free : bool):
        #General computation of FIR filter independantly of chosen model or hypothesis.
        data = jsonHandler.readJson("data/Glasberg2002.json") 
        # glasberg = dataPreparator()   
        fVec = np.linspace(kv["flow"],kv["fhigh"], abs(kv["flow"] - kv["fhigh"]) )
        # print("length of fVec is : ",len(fVec))
        glasbData = dataPreparator.OuterMiddle(data , fVec, model, free )
        
        tfLinear = 10**(glasbData.tfOuterMiddle/10)

        # print("len tfLinear is : ")
        # print(len(tfLinear))
        tfLinear[len(tfLinear) - 1] = 0.0 #This is imposed by the fir2 filter function because of the following error  : 
        #A Type II filter must have zero gain at the Nyquist frequency.

       
        # print("len linsp : ", len(linspace))
        # plt.plot(linspace, tfLinear)
        
        
        # print(linspace)
        outerMiddleFilter = FilterComputer.fir2(kv["order"], np.linspace(0,1,len(fVec)), tfLinear)
        print("out", len(outerMiddleFilter))
        # outerMiddleFilter = outerMiddleFilter[ : int(len(outerMiddleFilter)/2)] #seems to be more corresponding to the matlab behavior
       
        self.tfLinear = tfLinear
        self.outerMiddleFilter = outerMiddleFilter
        

    def FIR(self, inSig : np.array):
        
        b = self.outerMiddleFilter

        
        return(signal.filtfilt(b, a = 1, x = inSig))
        
    
        

    def fir2(order : float, domain : np.array, interpolatedValues : np.array):
        #!! numtaps et order n'ont pas le mm compotzment, decalage de +1
        
        return signal.firwin2(order,domain,interpolatedValues)
    




# FIR(np.array([2,0,3,5]) , "1997" , True)

if __name__ == '__main__':
    
    m = test.model(True)
    kv = m.kv
    # # f = FilterComputer()
    # y = sineMaker.makeSine(1000,0,0.1,40,49000)
    # y += sineMaker.makeSine(780,10,0.1,40,49000)
    # y += sineMaker.makeSine(3407,-24,0.1,40,49000)
    # print(len(y))
    
    F = FilterComputer(kv ,'1997', True)
    
    y = sineMaker.makeSine(1000,0,1,10,44000)
    y1 =  sineMaker.makeSine(100,0,1,40,44000)
    # y +=  sineMaker.makeSine(2890,0,1,10,44000)
    # y +=  sineMaker.makeSine(1300,0,1,10,44000)
    # y +=  sineMaker.makeSine(958,0,1,10,44000)
    
    # y = y/abs(max(y) - min(y))
    
    f = F.FIR(y)
    f1 = F.FIR(y1)
    

    bound = 1000
    x = np.linspace(0,bound,bound)
    plt.plot(f[:bound],'-c',linewidth=1,label= 'filtered')
    plt.plot(f1[:bound],'-k',linewidth=1,label='filtered2')
    plt.legend()
    plt.show()
    # plotfiltfilt(F.tfLinear,1,y)

    
    
    
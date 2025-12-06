
from filetools.jsonHandler import jsonHandler
from Moore1997._dataPreparator import dataPreparator
import Moore1997.model as test

from MathOperators.Signals import sineMaker

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
        outerMiddleFilter = outerMiddleFilter[ : int(len(outerMiddleFilter)/2)] #seems to be more corresponding to the matlab behavior
        # print(outerMiddleFilter)

        # plt.plot(np.linspace(0,1,len(outerMiddleFilter)), outerMiddleFilter)

        # plt.show()
        self.tfLinear = tfLinear
        self.outerMiddleFilter = outerMiddleFilter
        

    def FIR(self, inSig : np.array):
        from scipy import signal
        b = self.outerMiddleFilter



        print("sum_b is : "  ,  sum( self.outerMiddleFilter))
        return(signal.filtfilt(b, a = 1, x = inSig))
        
    
        

    def fir2(order : float, domain : np.array, interpolatedValues : np.array):
        #!! numtaps et order n'ont pas le mm compotzment, decalage de +1
        from scipy import signal
        return signal.firwin2(order,domain,interpolatedValues)
    

def plotfiltfilt(b,a,y):
    import matplotlib.pyplot as plt
    from scipy import signal

    t = np.linspace(0, 1.0, 2001)
    
    # b, a = signal.ellip(4, 0.01, 120, 0.125) #Only this filter works.



    rng = np.random.default_rng()

    n = 60

    # sig = rng.standard_normal(n)**3 + 3*rng.standard_normal(n).cumsum()  
    sig = y

    print("Computing Gust...")
    # fgust = signal.filtfilt(b, a, sig, method="gust")
    print("Computing other...")
    # fpad = signal.filtfilt(b, a, sig, padlen=50)

    fself = FilterComputer.filtfilt(b,a,sig)

    plt.plot(sig, 'k-',linewidth=3 ,label='input')

    # plt.plot(fgust, 'b-', linewidth=1, label='gust')

    # plt.plot(fpad, 'c-', linewidth=1.5, label='pad')

    plt.plot(fself, 'c-', linewidth=1.5, label='self')


    plt.legend(loc='best')

    plt.show()



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
    from scipy import signal
    y = sineMaker.makeSine(1000,0,1,10,44000)
    y1 =  sineMaker.makeSine(100,0,1,40,44000)
    # y +=  sineMaker.makeSine(2890,0,1,10,44000)
    # y +=  sineMaker.makeSine(1300,0,1,10,44000)
    # y +=  sineMaker.makeSine(958,0,1,10,44000)
    
    # y = y/abs(max(y) - min(y))
    
    f = F.FIR(y)
    f1 = F.FIR(y1)
    
    
    
    import matplotlib.pyplot as plt
    bound = 1000
    x = np.linspace(0,bound,bound)
    plt.plot(f[:bound],'-c',linewidth=1,label= 'filtered')
    plt.plot(f1[:bound],'-k',linewidth=1,label='filtered2')
    plt.legend()
    plt.show()
    # plotfiltfilt(F.tfLinear,1,y)

    
    
    
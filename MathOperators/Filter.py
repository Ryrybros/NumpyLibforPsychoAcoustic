
from filetools.jsonHandler import jsonHandler
import numpy as np
from Glasberg2002.GlasbergComputer import GlasbergComputer
import matplotlib.pyplot as plt

class FilterComputer :
    def FIR(fVec : np.array, model : str, free : bool):
        #General computation of FIR filter independantly of chosen model or hypothesis.
        data = jsonHandler.readJson("data/Glasberg2002.json") 
        glasberg = GlasbergComputer()   
        glasbData = glasberg.OuterMiddle(fVec, model, free )
        
        tfLinear = 10**(glasbData.tfOuterMiddle/10)

        print("tfLinear is : ")
        print(tfLinear)
        tfLinear[len(tfLinear) - 1] = 0.0 #This is imposed by the fir2 filter function because of the following error  : 
        #A Type II filter must have zero gain at the Nyquist frequency.

        linspace = np.array([(1/(len(fVec) - 1 ))*i for i in range(len(fVec))])

        # plt.plot(linspace, tfLinear)
        # plt.show()

        # print(linspace)
        outerMiddleFilter = FilterComputer.fir2(150, linspace, tfLinear)
        # print(outerMiddleFilter)


        

    def fir2(order : float, domain : np.array, interpolatedValues : np.array):
        from scipy import signal
        
        return signal.firwin2(order, domain, interpolatedValues )
    

        

def sineMaker(f,time,dB,cutfreq):
    x = np.array([ (time*(i/cutfreq)) for i in range(cutfreq) ])#0:(time/cutfreq):time)
    if(dB != 0):
      sine = ((10**(-3 - 17/20))*(10 ** (abs(dB)/20)))*np.sin((2*np.pi*f)*x)
    else:
        sine = ((10^(-3 - 17/20)))*np.sin((2*np.pi*f)*x)


    # plt.plot(x, sine)
    # plt.show()
    return sine


# FIR(np.array([2,0,3,5]) , "1997" , True)

if __name__ == '__main__':
    f = FilterComputer()
    y = sineMaker(10,1,40,20000)
    FilterComputer.FIR(y ,'1997', True)
    
    
    

    




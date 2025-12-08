from collections import defaultdict
from Moore1997._dataPreparator import dataPreparator
from filetools.jsonHandler import jsonHandler
from MathOperators.ERBscale import ERB
import numpy as np
class model:
    
    #______________PRETREATMENT______________
    def __init__(self, free : bool ,**kwargs):
        
        defaultKV = {
        "fs":  32000,
        "flow":  20,
        "fhigh":  16000,
        "order":  4096,
        "erbStep":  0.2500,
        "erbFcMin":  50,
        "erbFcMax":  15000
        }
        kv = defaultdict(float)
        
        for key in defaultKV:
            kv[key] = defaultKV[key]
        
        if kwargs:
            for key, value in kwargs.items(): 
                if(key in defaultKV):
                    kv[key] = value
                else: print("Warning : invalid keys given to kv")
        
        self.kv = kv
        self.fVec = [kv["flow"] + i  + 1 for i in range(kv["fhigh"] - kv["flow"] )]
        
        self.data = jsonHandler.readJson('data/Glasberg2002.json')
        self.OuterMiddle = dataPreparator.OuterMiddle(self.data,self.fVec,"1997",free)

        self._erbScale()
        self._specLoud()


        
    def _erbScale(self):
        self.erbNMin = ERB.f2erbrate(self.kv["erbFcMin"])
        self.erbNMax = ERB.f2erbrate(self.kv["erbFcMax"])
        
        self.erbN = [self.erbNMin + self.kv["erbStep"]*i for i in range(int(self.erbNMax - self.erbNMin))] #erbNMin:kv.erbStep:erbNMax    # numbers of erb bands
        self.erbFc = ERB.erbrate2f(self.erbN)        # center frequency of erb bands

        self.erbLoFreq = ERB.erbrate2f(self.erbN-0.5*np.ones(len(self.erbN))) # lower limit of each ERB filter
        self.erbHiFreq = ERB.erbrate2f(self.erbN+0.5*np.ones(len(self.erbN))) # upper limit of each ERB filter

    def _specLoud(self):
        
        dataSL = dataPreparator.SpecLoudness(self.data,self.erbFc)
        class res :
            def __init__(self):
                        
                self.tQdB = dataSL.tQ
                self.tQ = 10**(self.tQdB/10)
                self.tQdB500 = dataSL.tQ500
                # %gdB = dataSL.g;    % low level gain in cochlea amplifier
                self.g = 10**((self.tQdB500-self.tQdB)/10)
                self.a = dataSL.a    #% parameter for linearization around absolute threshold
                self.alpha = dataSL.alpha #    % compressive exponent
                self.c = dataSL.c # % constant to get loudness scale to sone
        return res()


    #_____________________________________________END OF PRETREATMENT_____________________________________________


        


if __name__ == '__main__':
    m = model(free = True,**{"fhigh" : 25})
    # print(m.erbN)
    import matplotlib.pyplot as plt
    print( int(2000/len(m.erbN)) )
    print( len(m.erbN))
    # plt.plot(np.linspace(0,2000, len(m.erbN) ),m.erbN)
    # plt.plot(np.linspace(0,2000, len(m.erbNMax) ),m.erbNMax)
    # plt.plot(np.linspace(0,2000, len(m.erbNMin) ),m.erbNMin)

    # plt.show()

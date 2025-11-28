from collections import defaultdict
from Moore1997.dataPreparator import dataPreparator
from filetools.jsonHandler import jsonHandler
from MathOperators.ERBscale import ERB
class model:

    def __init__(self, free : bool ,**kwargs):
        
        data = {
        "fs":  32000,
        "flow":  20,
        "fhigh":  16000,
        "order":  4096,
        "erbStep":  0.2500,
        "erbFcMin":  50,
        "erbFcMax":  15000
        }
        kv = defaultdict(float)
        
        for key in data:
            kv[key] = data[key]
        
        if kwargs:
            for key, value in kwargs.items(): 
                if(key in data):
                    kv[key] = value
                else: print("Warning : invalid keys given to kv")
        
        self.kv = kv
        self.fVec = [kv["flow"] + i  + 1 for i in range(kv["fhigh"] - kv["flow"] )]

        data = jsonHandler.readJson('filetools/data.json')
        # print(data)
        self.OuterMiddle = dataPreparator.OuterMiddle(data,self.fVec,"1997",free)
        # print(self.OuterMiddle.tfOuterMiddle )
        self._erbScale()
        print(self.erbN)

        
    def _erbScale(self):
        self.erbNMin = ERB.f2erbrate(self.kv["erbFcMin"])
        self.erbNMax = ERB.f2erbrate(self.kv["erbFcMax"])
        
        self.erbN = [self.erbNMin + self.kv["erbStep"]*i for i in range(int(self.erbNMax - self.erbNMin))] #erbNMin:kv.erbStep:erbNMax    # numbers of erb bands
        self.erbFc = ERB.erbrate2f(self.erbN);               # center frequency of erb bands

        self.erbLoFreq = ERB.erbrate2f(self.erbN-0.5); # lower limit of each ERB filter
        self.erbHiFreq = ERB.erbrate2f(self.erbN+0.5); # upper limit of each ERB filter

        


if __name__ == '__main__':
    m = model(free = True,**{"fhigh" : 25})


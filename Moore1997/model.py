from collections import defaultdict
from Moore1997._dataPreparator import dataPreparator
from filetools.jsonHandler import jsonHandler
from MathOperators.ERBscale import ERB
from MathOperators.Interpolator import Interpolator
from MathOperators import FFTTreatmment
from MathOperators import Signals
import numpy as np

class Result:

    def __init__(self):
        self.computed = False
        self.eLdB = None
        self.erbN = None
        self.fc = None

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
        self.specLoudData = self._specLoud()

        self.fftValues = None

        self.results = Result()


        
    def _erbScale(self):
        self.erbNMin = ERB.f2erbrate(self.kv["erbFcMin"])
        self.erbNMax = ERB.f2erbrate(self.kv["erbFcMax"])
        # print("erbMax : ", self.erbNMax)
        # self.erbN = np.array( [self.erbNMin + self.kv["erbStep"]*i for i in range(int((self.erbNMax - self.erbNMin)/self.kv["erbStep"]))] )#erbNMin:kv.erbStep:erbNMax    # numbers of erb bands
        self.erbN = np.arange(self.erbNMin , self.erbNMax , self.kv["erbStep"]  ) 
        
        self.erbFc = ERB.erbrate2f(self.erbN)        # center frequency of erb bands

        self.erbLoFreq = ERB.erbrate2f(self.erbN - 0.5 ) # lower limit of each ERB filter
        self.erbHiFreq = ERB.erbrate2f(self.erbN + 0.5 ) # upper limit of each ERB filter

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

    def _excitationPatern(self, earSig : np.array):
        #calculate intensity for each ERB (dB/ERB)
        #Use it after Pretreatment

        self.fftValues = FFTTreatmment.FFTComputer.computeFFT(earSig = earSig , fs= self.kv['fs'])

        erbInt = np.zeros(len(self.erbFc))

        for i in range(len(self.erbFc) ) :
            
            loValue = round(self.erbLoFreq[i]*self.fftValues.oneHz)
            hiValue = round(self.erbHiFreq[i]*self.fftValues.oneHz)
            # print("hi Val : " ,hiValue)

            erbRange = np.linspace( loValue ,hiValue, hiValue - loValue  , dtype= int)
            # print(type(int(erbRange[0])))
            sumList =  np.zeros(len(erbRange)) 
            # print("erbRange ", len(erbRange) )
            j = 0
            for index in  erbRange :
                # i = int( index )
                # print(type(i) )
                
                if(index < len(self.fftValues.compInt)):
                    sumList[j] = self.fftValues.compInt[ index ]   # intensity in each erb
                else:
                    if((index == max(erbRange)) & ( i >= len(self.erbFc) - 1 ) ):
                        print("Warning !! Offset between erbRange and length of compInt ,\n maximum offset is : ", index - len(self.fftValues.compInt), ".\n")
                j += 1

            erbInt[i] = sumList.sum()   # intensity sum in each erb
        

        erbdB = 10*np.log10(erbInt/ ( (20e-6)**2 ))   # intensity level in each erb using reference SPL of 20 uPa
        p511 = 4*1000/ERB.f2erb(1000)    # p for fc=1kHz and a level of 51dB (at 1kHz filters are symmetrical)
        
        # erbdB2F = np.interp(  self.fftValues.compFq , [0,self.erbFc,self.kv['fs']/2], [min(erbdB), erbdB ,min(erbdB)] )   # map erbFc to compFq
        # print("\n erbFc : " , self.erbFc)
        erbdB2F = np.zeros(len(self.erbFc))
        # print("\n compFq : " , self.fftValues.compFq)
        # print( "\n erbdB : " , erbdB )
        # print("\n erbdB min : " ,min(erbdB))
        #for ind in range(len(self.erbFc)):

        erbdB2F = np.interp(x= self.fftValues.compFq, xp = np.concatenate( ( [0],self.erbFc,[self.kv['fs']/2] ) ) , fp = np.concatenate( ( [min(erbdB) ], erbdB , [min(erbdB)] ) ) ) # map erbFc to compFq

        eL = np.zeros(len(self.erbN))
        
        for e in range( len(self.erbN) ) :
            erb = ERB.f2erb(self.erbFc[e])
            p51 = 4*self.erbFc[e]/erb
            intensity = 0
            for comp in range (self.fftValues.nPoints ) :
                g = (self.fftValues.compFq[comp]-self.erbFc[e])/self.erbFc[e]
                if g<0 :
                    p = p51 - 0.35*(p51/p511) * (erbdB2F[comp] - 51)
                else :
                    p = p51
                
                g = abs(g)
                w = (1+p*g)*np.exp(-p*g)
                if( (p*g >  10**10 )):
                    w = 10**(-10)
                
                intensity = intensity  +  w  *  self.fftValues.compInt[comp] #intensity per erb
            
            eL[e] = intensity
        self._eL = eL
        self.results.eLdB = 10*np.log10(eL / ( (20e-6)**2 ) ) # get dB SPL (20uPa reference)
        
        self.results.erbN = self.erbN
        self.results.fc = self.erbFc


    def moore1997(self, earSig : np.array) :
        self._excitationPatern(earSig)
        specLoud = np.zeros(len(self._eL))
        # c*(2*eL./(eL+tQ)).^1.5 .*((g.* eL + a).^alpha-a.^alpha)
        specLoud1 = self.specLoudData.c *  ( (2*self._eL/( self._eL + self.specLoudData.tQ ))**1.5 ) *   ( (self.specLoudData.g* self._eL + self.specLoudData.a) ** self.specLoudData.alpha - self.specLoudData.a**self.specLoudData.alpha ) #% Eq. 6?
        specLoud2 = self.specLoudData.c * (  (self.specLoudData.g * self._eL + self.specLoudData.a)**self.specLoudData.alpha - self.specLoudData.a**self.specLoudData.alpha); #% Eq. 8?
        specLoud3 = self.specLoudData.c * ( self._eL/(1.04*(10**6)))**0.5 #% Eq. 9?
        
        specLoud[ self._eL < self.specLoudData.tQ ] = specLoud1[self._eL < self.specLoudData.tQ]
        specLoud[ ( self._eL <= 10**10 ) & ( self._eL  > self.specLoudData.tQ )  ] = specLoud2[ (self._eL <= 10**10 ) & ( self._eL> self.specLoudData.tQ ) ]
        specLoud[self._eL > 10**10] = specLoud3[self._eL > 10**10]# % end of Sec. 1.6 in the paper

        monauralLoudness = sum(specLoud,2) * self.kv['erbStep'] #     % integrate over the erbs
        
        self.results.monauralLoudness = monauralLoudness #     % integrate over the erbs
        self.results.specLoud = specLoud # % specific loudness

        class modelresult : 
            def __init__(self, specLoudness, monoLoudness):
                self.specLoudness = specLoudness
                self.monauralLoudness  = monoLoudness

        return modelresult(specLoud, monauralLoudness)
 





        


if __name__ == '__main__':
    m = model(free = True,**{"fhigh" : 25})
    # print(m.erbN)
    import matplotlib.pyplot as plt
    print( int(2000/len(m.erbN)) )
    print( len(m.erbN))
    y = Signals.sineMaker.makeSine(1000,0,1,20,m.kv['fs']) #it is crucial that both signals matlab/python have the same parameters (time is important)
    
    res = m.moore1997(y)
    print(res.specLoudness)
    

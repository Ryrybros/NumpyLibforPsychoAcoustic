from collections import defaultdict
from Moore1997._dataPreparator import dataPreparator
from Moore1997.Filter import FilterComputer
from filetools.jsonHandler import jsonHandler
from MathOperators.ERBscale import ERB
from MathOperators.Interpolator import Interpolator
from MathOperators import FFTTreatmment
from MathOperators import Signals
import numpy as np
import time 
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

    def _excitationPatern(self, earSig : np.array,e0 = None):
        #calculate intensity for each ERB (dB/ERB)
        #Use it after Pretreatment
        filter = FilterComputer(self.kv, "1997", True)
        sig = filter.FIR(earSig)
        self.sigtest = sig 
        self.fftValues = FFTTreatmment.FFTComputer.computeFFT(earSig = sig , fs= self.kv['fs'])
        t = time.time()
         
        erbInt = np.zeros(len(self.erbFc))
        
        for i in range(len(self.erbFc) ) :
            
            loValue = round(self.erbLoFreq[i]*self.fftValues.oneHz)
            hiValue = round(self.erbHiFreq[i]*self.fftValues.oneHz)
            # print("hi Val : " ,hiValue)

            # erbRange = np.linspace( loValue ,hiValue, hiValue - loValue  , dtype= int)
            # print(type(int(erbRange[0])))

            erbInt[i] = np.sum(self.fftValues.compInt[loValue : hiValue + 1])

            # sumList =  np.zeros(len(erbRange)) 
            # # print("erbRange ", len(erbRange) )
            # j = 0
            # for index in  erbRange :
            #     # i = int( index )
            #     # print(type(i) )
                
            #     if(index < len(self.fftValues.compInt)):
            #         sumList[j] = self.fftValues.compInt[ index ]   # intensity in each erb
            #     else:
            #         if((index == max(erbRange)) & ( i >= len(self.erbFc) - 1 ) ):
            #             print("Warning !! Offset between erbRange and length of compInt ,\n maximum offset is : ", index - len(self.fftValues.compInt), ".\n")
            #     j += 1

            # erbInt[i] = sumList.sum()   # intensity sum in each erb
        
        # print("--------------------------------------",sum(erbInt))
        # print("time loop  : " , time.time() - t)
        
        
        
        # t = time.time()
        # loV = np.round(self.erbLoFreq*self.fftValues.oneHz).astype(int)
        # hiV = np.round(self.erbHiFreq*self.fftValues.oneHz).astype(int)
        # cumsum_cp = np.concatenate(([0], np.cumsum(self.fftValues.compInt)))
        # e2 = cumsum_cp[hiV + 1] - cumsum_cp[loV]

        # print("--------------------------------------",sum(e2))
        # print("time np : " , time.time() - t)

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
        self.test = erbdB2F 

        t = time.time()


        # eL = np.zeros(len(self.erbN))
        
        # for e in range( len(self.erbN) ) :
        #     erb = ERB.f2erb(self.erbFc[e])
        #     p51 = 4*self.erbFc[e]/erb
        #     intensity = 0
        #     for comp in range (self.fftValues.nPoints ) :
        #         g = (self.fftValues.compFq[comp]-self.erbFc[e])/self.erbFc[e]
        #         if g<0 :
        #             p = p51 - 0.35*(p51/p511) * (erbdB2F[comp] - 51)
        #         else :
        #             p = p51
                
        #         g = abs(g)
        #         w = (1+p*g)*np.exp(-p*g)
                
        #         intensity = intensity  +  w  *  self.fftValues.compInt[comp] #intensity per erb
            
        #     eL[e] = intensity
        
        
        print("-----------------double loop-------------------",)
        print("time loop : " , time.time() - t)


        #_______________________________test_____________________________
        #This is a vectorized version of the nested loop. It is faster by a factor of more than 10
        t = time.time()


        erb = ERB.f2erb(self.erbFc)
        p51 = 4 * self.erbFc / erb

        fc = self.erbFc[:, np.newaxis]
        fq = self.fftValues.compFq[np.newaxis, :]
        int_vals = self.fftValues.compInt[np.newaxis, :]

        g_raw = (fq - fc) / fc

        db_vals = erbdB2F[np.newaxis, :]
        p_low = p51[:, np.newaxis] - 0.35 * (p51[:, np.newaxis] / p511) * (db_vals - 51)
        p_high = p51[:, np.newaxis]

        p = np.where(g_raw < 0, p_low, p_high)


        g_abs = np.abs(g_raw)
        w = (1 + p * g_abs) * np.exp(-p * g_abs)

        eL2 = np.sum(w * int_vals, axis=1)



    
        print("-----------------double loop-------------------",len(eL2))
        print("time vector : " , time.time() - t)
        eL = eL2
        print(max(eL2 - eL) )



        #______________________________end of test____________________________________________

        if(e0 != None) : E0 = e0
        else: E0 = (20e-6)**2 

        self._eL = eL / E0
        
        self.results.eLdB = 10*np.log10( self._eL  ) # get dB SPL (20uPa reference)
        
        self.results.erbN = self.erbN
        self.results.fc = self.erbFc

    

    def moore1997(self, earSig : np.array, e0 = None) :
        
        self._excitationPatern(earSig)
        specLoud = np.zeros(len(self._eL))
        # c*(2*eL./(eL+tQ)).^1.5 .*((g.* eL + a).^alpha-a.^alpha)
        specLoud1 = self.specLoudData.c *  ( (2*self._eL/( self._eL + self.specLoudData.tQ ))**1.5 ) *   ( (self.specLoudData.g* self._eL  + self.specLoudData.a) ** self.specLoudData.alpha - self.specLoudData.a**self.specLoudData.alpha ) #% Eq. 6?
        specLoud2 = self.specLoudData.c * (  (self.specLoudData.g * self._eL + self.specLoudData.a)**self.specLoudData.alpha - self.specLoudData.a**self.specLoudData.alpha); #% Eq. 8?
        specLoud3 = self.specLoudData.c * ( self._eL/(1.04*(10**6)))**0.5 #% Eq. 9?
        
        specLoud[ self._eL < self.specLoudData.tQ ] = specLoud1[self._eL < self.specLoudData.tQ]
        specLoud[ ( self._eL <= 10**10 ) & ( self._eL  > self.specLoudData.tQ )  ] = specLoud2[ (self._eL <= 10**10 ) & ( self._eL> self.specLoudData.tQ ) ]
        specLoud[self._eL > 10**10] = specLoud3[self._eL > 10**10]# % end of Sec. 1.6 in the paper

        monauralLoudness = sum(specLoud) * self.kv['erbStep'] #     % integrate over the erbs
        
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
    
# -----------------double loop-------------------
# time loop :  14.322951316833496
# -----------------double loop------------------- 149
# time vector :  0.19310355186462402
# 9.26442286059391e-23
# Elapsed time for step  : 20.0 dB :  17.258923053741455
# -----------------double loop-------------------
# time loop :  21.742462873458862
# -----------------double loop------------------- 149
# time vector :  0.27115488052368164
# 5.823351512373315e-22
# Elapsed time for step  : 30.0 dB :  32.16646695137024
# -----------------double loop-------------------
# time loop :  21.298561573028564
# -----------------double loop------------------- 149
# time vector :  0.2699263095855713
# 1.3552527156068805e-20
# Elapsed time for step  : 40.0 dB :  31.16813015937805
# -----------------double loop-------------------
# time loop :  22.262532234191895
# -----------------double loop------------------- 149
# time vector :  0.2755773067474365
# 1.0842021724855044e-19
# Elapsed time for step  : 50.0 dB :  31.815988302230835
# -----------------double loop-------------------
# time loop :  20.951077461242676
# -----------------double loop------------------- 149
# time vector :  0.2379288673400879
# 7.047314121155779e-19
# Elapsed time for step  : 60.0 dB :  31.460907220840454
# -----------------double loop-------------------
# time loop :  5.854139804840088
# -----------------double loop------------------- 149
# time vector :  0.10844111442565918
# 1.734723475976807e-17
# Elapsed time for step  : 70.0 dB :  15.711686134338379
# -----------------double loop-------------------
# time loop :  5.193831920623779
# -----------------double loop------------------- 149
# time vector :  0.0963582992553711
# 1.1102230246251565e-16
# Elapsed time for step  : 80.0 dB :  7.638150453567505
# -----------------double loop-------------------
# time loop :  4.737797021865845
# -----------------double loop------------------- 149
# time vector :  0.10350465774536133
# 1.1657341758564144e-15
# Elapsed time for step  : 90.0 dB :  7.144345760345459
# -----------------double loop-------------------
# time loop :  4.876959800720215
# -----------------double loop------------------- 149
# time vector :  0.10536885261535645
# 1.1102230246251565e-14
# Elapsed time for step  : 100.0 dB :  7.512490749359131
# -----------------double loop-------------------
# time loop :  4.697258949279785
# -----------------double loop------------------- 149
# time vector :  0.08721685409545898
# 3.637978807091713e-12
# Elapsed time for curve 0 plot:  7.150695085525513
# -----------------double loop-------------------
# time loop :  4.837356805801392
# -----------------double loop------------------- 149
# time vector :  0.09795331954956055
# 3.197442310920451e-14
# Elapsed time for curve 1 plot:  7.225187063217163
# -----------------double loop-------------------
# time loop :  4.847031354904175
# -----------------double loop------------------- 149
# time vector :  0.10391998291015625
# 7.993605777301127e-15
# Elapsed time for curve 2 plot:  7.297339677810669
# -----------------double loop-------------------
# time loop :  4.750588417053223
# -----------------double loop------------------- 149
# time vector :  0.09364652633666992
# 1.1102230246251565e-15
# Elapsed time for curve 3 plot:  7.123446464538574
# -----------------double loop-------------------
# time loop :  5.026240587234497
# -----------------double loop------------------- 149
# time vector :  0.10103487968444824
# 6.938893903907228e-17
# Elapsed time for curve 4 plot:  7.41409158706665
# -----------------double loop-------------------
# time loop :  4.883578062057495
# -----------------double loop------------------- 149
# time vector :  0.0947878360748291
# 2.481541837659083e-23
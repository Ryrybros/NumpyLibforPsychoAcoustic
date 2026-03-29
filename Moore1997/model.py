from collections import defaultdict
from Moore1997._dataPreparator import dataPreparator
from Moore1997.Filter import FilterComputer
from filetools.jsonHandler import jsonHandler
from MathOperators.ERBscale import ERB
from MathOperators.Interpolator import Interpolator
from MathOperators import FFTTreatmment
from MathOperators import Signals
import scipy.signal.windows as win
import scipy.fft as fft

import numpy as np
import time 
import math

import matplotlib.pyplot as plt
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
        "erbFcMax":  15000,

        #Glasberg
        "timeStep": 0.001,
        "hannLenMs": [ 2 , 4 , 8 , 16 , 32 , 64 ],
        "fftLen": 2048,
        "vLimitingIndices" :  [ 4050   , 2540,   1250,    500   ,  80 ]
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

        #________________________format_parameters________________________________

        self.specLoudData = self._specLoud()

        self.fftValues = None

        self.results = Result()


        #___________________________________Glasberg__________________________

        #MATLAB________________________________________________________
        # numBlocks = ceil(length(earSig)./updateRate);
        # earSigPad = earSig;
        # earSigPad(end+1:end+hannLenSmp(6)) = zeros(hannLenSmp(6),1);  % zero padding
        #MATLAB________________________________________________________

        self.updateRate = np.round(kv["timeStep"]*kv["fs"])
        self.hannLenSmp = np.round(np.array(kv["hannLenMs"])/1000 * kv["fs"]) # windows size in samples

        self.hannWin = []
        
        for val in self.hannLenSmp:
            self.hannWin.append(win.hann(int(val)))
        
        
            
    
       



    
        
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


    def stationnarySpect(self, earSig : np.array , stationnary = True):
        #calculate intensity for each ERB (dB/ERB)
        #Use it after Pretreatment
        filter = FilterComputer(self.kv, "1997", True)
        sig = filter.FIR(earSig)
        #COmpint, COmpFq computation
        self.fftValues = FFTTreatmment.FFTComputer.computeFFT(earSig = sig , fs= self.kv['fs'])
       
    def getEL(self, compInt, compFq):
        
         
        erbInt = np.zeros(len(self.erbFc))
        
        for i in range(len(self.erbFc) ) :
            
            loValue = round(self.erbLoFreq[i]*self.fftValues.oneHz)
            hiValue = round(self.erbHiFreq[i]*self.fftValues.oneHz)

            erbInt[i] = np.sum(compInt[loValue : hiValue + 1])
        
        
        
        erbdB = 10*np.log10(erbInt/ ( (20e-6)**2 ))   # intensity level in each erb using reference SPL of 20 uPa
        p511 = 4*1000/ERB.f2erb(1000)    # p for fc=1kHz and a level of 51dB (at 1kHz filters are symmetrical)
        
        erbdB2F = np.zeros(len(self.erbFc))

        erbdB2F = np.interp(x= compFq, xp = np.concatenate( ( [0],self.erbFc,[self.kv['fs']/2] ) ) , fp = np.concatenate( ( [min(erbdB) ], erbdB , [min(erbdB)] ) ) ) # map erbFc to compFq
        self.test = erbdB2F 

        t = time.time()
        print("-----------------double loop-------------------",)
        print("time loop : " , time.time() - t)


        #_______________________________test_____________________________
        #This is a vectorized version of the nested loop. It is faster by a factor of more than 10
        t = time.time()


        erb = ERB.f2erb(self.erbFc)
        p51 = 4 * self.erbFc / erb

        fc = self.erbFc[:, np.newaxis]
        fq = compFq[np.newaxis, :]
        int_vals = compInt[np.newaxis, :]

        g_raw = (fq - fc) / fc

        db_vals = erbdB2F[np.newaxis, :]
        p_low = p51[:, np.newaxis] - 0.35 * (p51[:, np.newaxis] / p511) * (db_vals - 51)
        p_high = p51[:, np.newaxis]

        p = np.where(g_raw < 0, p_low, p_high)


        g_abs = np.abs(g_raw)
        w = (1 + p * g_abs) * np.exp(-p * g_abs)

        eL = np.sum(w * int_vals, axis=1)

        return eL

        

    #__________________________________________Moore model__________________________________________________
    def _excitationPatern(self ,e0 = None):
                
        eL = self.getEL( self.fftValues.compInt, self.fftValues.compFq )

    
        print("-----------------double loop-------------------",len(eL))
        # print("time vector : " , time.time() - t)
   

        #______________________________end of test____________________________________________

        if(e0 != None) : E0 = e0
        else: E0 = (20e-6)**2 
    
        self._eL = eL / E0
        print("_eL")
        
        self.results.eLdB = 10*np.log10( self._eL  ) # get dB SPL (20uPa reference)
        
        self.results.erbN = self.erbN
        self.results.fc = self.erbFc

    

    def specLoudness(self) :
        

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
    
    def moore1997(self, earSig : np.array,e0 = None):
        self.stationnarySpect(earSig)
        self._excitationPatern( e0= e0)
        print("done")
        self.res = self.specLoudness()
        print("end moore")
        


    def glasbergSpect(self, earSig : np.array ):
        
        
        earSigPad = np.concat([ earSig, np.zeros( int( self.hannLenSmp[5]  ))  ] )

        numBlocks = np.ceil(len(earSig)/self.updateRate)
        
       
        block_indices = np.arange(numBlocks)
        lower_bounds = block_indices * self.updateRate

        # Pre-calculate all Hann windows and store in a list
        
        hannWins = [win.hann(int(length)) for length in self.hannLenSmp]

        # Process each window size
        # We store the results in a list of matrices
        spectra = []

        
        
        for i, win_len in enumerate(self.hannLenSmp):
            win_len = int(win_len)
            
            # Create a 2D matrix : (numBlocks, win_len)
            
            offsets = np.arange(win_len)
            idx_matrix = lower_bounds[:, np.newaxis] + offsets
            
            
            segments = earSigPad[idx_matrix.astype(int)] * hannWins[i]
            
            
            # Result is a (numBlocks, fftLen)
            res_fft = np.fft.fft(segments, n=self.kv["fftLen"], axis=1)
            spectra.append(res_fft)
        spectList = np.abs(spectra)
        



                # 1. Setup boundaries
        half_fft = int(self.kv["fftLen"] // 2)
        oneHz = (self.kv["fftLen"] + 2) / self.kv["fs"]

        # Create a full set of indices: [Nyquist, 4050, 2540, 1250, 500, 80, 0]
        # Note: we reverse them or order them to match the spectra list indices (0 to 5)
        boundaries = np.round(np.array(self.kv["vLimitingIndices"]) * oneHz).astype(int)
        all_lims = np.concatenate(([half_fft + 1], boundaries, [0]))

        # 2. Vectorized Loop for Stitching
        spect = np.zeros((int(numBlocks), half_fft + 1))

        # spectra[0] corresponds to the range all_lims[0] to all_lims[1], etc.
        for i in range(6):
            start, end = all_lims[i+1], all_lims[i]
            # Calculate normalization factor (sum of squares of current Hann window)
            norm_factor = np.sum(hannWins[i]**2)
            
            # Assign the frequency slice across all time blocks
            spect[:, start:end] = np.abs(spectra[i][:, start:end])**2 / norm_factor

        
        # 3. Final Calculations
        compInt = 2 * spect / self.kv["fs"]
        compFq = np.linspace(0, self.kv["fs"] / 2, half_fft + 1)
        class returned :
            def __init__(self, x , y):
                self.compInt =   x
                self.compFq = y
        
        
        return returned(compInt, compFq)
        
        # print(compInt)
        
        #____

        

    def glasberg2002(self, inSig ):
        print("comin")
        # print ( len(self._excitationPatern(inSig.compInt) ) )
        sp = self.glasbergSpect(inSig)
        
        print("will start over arrays")
        eLs = np.array([
            self.getEL(sp.compInt[:, i], sp.compFq)
            for i in range(sp.compInt.shape[1])
        ])

        return results

        


if __name__ == '__main__':
    
    m = model(free = True)
    sig = np.sin(2*np.pi*1000*np.linspace(0,3,3*44100))
    ex = m.moore1997(sig)
    print("__________________DOne________________________")
    plt.plot(m.res.specLoudness)
    plt.show()
    m.glasberg2002( sig ) 

    
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
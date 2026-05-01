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
from scipy.signal import resample_poly
from MathOperators.Signals import *
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

        self.fVec = [kv["flow"] + i  for i in range(kv["fhigh"] - kv["flow"] )]
        
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
                self.tQ[0] = self.tQ[1]
                self.tQdB500 = dataSL.tQ500
                # %gdB = dataSL.g;    % low level gain in cochlea amplifier
                self.g = 10**((self.tQdB500-self.tQdB)/10)
                self.g_prior = dataSL.g 
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
        return FFTTreatmment.FFTComputer.computeFFT(earSig = sig , fs= self.kv['fs'])
       

       
    def getEL(self, fftValues = None, model = 'Moore1997'):
        
        eps = 1e-40
        erbInt = np.ones(len(self.erbFc)) * eps
        oneHz = fftValues.oneHz
        

        for i in range(len(self.erbFc) ) :
            
            loValue = round(self.erbLoFreq[i]*oneHz)
            hiValue = round(self.erbHiFreq[i]*oneHz)

            erbInt[i] = np.sum(fftValues.compInt[loValue : hiValue + 1]) + eps
        
        
        
        erbdB = 10*np.log10(erbInt/ ( (20e-6)**2 ))   # intensity level in each erb using reference SPL of 20 uPa
        p511 = 4*1000/ERB.f2erb(1000)    # p for fc=1kHz and a level of 51dB (at 1kHz filters are symmetrical)
        
        erbdB2F = np.zeros(len(self.erbFc))

        erbdB2F = np.interp(x= fftValues.compFq, xp = np.concatenate( ( [0],self.erbFc,[self.kv['fs']/2] ) ) , fp = np.concatenate( ( [min(erbdB) ], erbdB , [min(erbdB)] ) ) ) # map erbFc to compFq
        self.test = erbdB2F 


        erb = ERB.f2erb(self.erbFc)
        p51 = 4 * self.erbFc / erb

        fc = self.erbFc[:, np.newaxis]
        fq = fftValues.compFq[np.newaxis, :]
        int_vals = fftValues.compInt[np.newaxis, :]

        g_raw = (fq - fc) / fc

        db_vals = erbdB2F[np.newaxis, :]
        p_low = p51[:, np.newaxis] - 0.35 * (p51[:, np.newaxis] / p511) * (db_vals - 51)
        p_high = p51[:, np.newaxis]


        #Distinguishing  Moore and Glasb
        
        p = np.where(g_raw < 0, p_low, p_high)
        
        g_abs = np.abs(g_raw)
    
            
        w = (1 + p * g_abs) * np.exp(-p * g_abs)
        
        eL = np.sum(w * int_vals, axis=1)

    

        return eL

        

    #__________________________________________Moore model__________________________________________________
    def _excitationPatern(self ,fftValues,e0 = None):
        eL = self.getEL( fftValues= fftValues )
        eL = eL / (20e-6)**2 
        
        self.results.eLdB = 10*np.log10( eL  ) # get dB SPL (20uPa reference)
        
        self.results.erbN = self.erbN
        self.results.fc = self.erbFc
        return eL

    

    def specLoudness(self, eL : np.array) :
        

        specLoud = np.zeros(len(eL))
        
        # c*(2*eL./(eL+tQ)).^1.5 .*((g.* eL + a).^alpha-a.^alpha)
        
        specLoud1 = self.specLoudData.c *  ( (2*eL/( eL + self.specLoudData.tQ ))**1.5 ) *   ( (self.specLoudData.g* eL + self.specLoudData.a) ** self.specLoudData.alpha - self.specLoudData.a**self.specLoudData.alpha ) #% Eq. 6? 
        specLoud2 =     self.specLoudData.c * (  (self.specLoudData.g * eL + self.specLoudData.a)**self.specLoudData.alpha - self.specLoudData.a**self.specLoudData.alpha); #% Eq. 8?
        specLoud3 = self.specLoudData.c * ( eL/(1.04*(10**6)))**0.5 #% Eq. 9?
        
        specLoud[ eL < self.specLoudData.tQ ] = specLoud1[eL < self.specLoudData.tQ]
        specLoud[ ( eL <= 10**10 ) & ( eL  > self.specLoudData.tQ )  ] = specLoud2[ (eL <= 10**10 ) & ( eL> self.specLoudData.tQ ) ]
        specLoud[eL > 10**10] = specLoud3[eL > 10**10]# % end of Sec. 1.6 in the paper
        
        monauralLoudness = sum(specLoud) * self.kv['erbStep'] #     % integrate over the erbs
        
        self.results.monauralLoudness = monauralLoudness #     % integrate over the erbs
        self.results.specLoud = specLoud # % specific loudness

        class modelresult : 
            def __init__(self, specLoudness, monoLoudness):
                self.specLoudness = specLoudness
                self.monauralLoudness  = monoLoudness

        return modelresult(specLoud, monauralLoudness)
    
    def getSpecLoudness(self, eL : np.array):
        return self.specLoudness(eL).specLoudness



    class mooreReturned :
        def __init__(self, loud, eL):
            self.Loudness = loud
            self.eL = eL
            
    def moore1997(self, earSig : np.array,e0 = None, fs = 32000):

        if(fs != self.kv['fs']):
            earSig = resample_poly(earSig, self.kv['fs'], fs)
        fftVals = self.stationnarySpect(earSig)
        eL = self._excitationPatern( e0= e0, fftValues= fftVals )
        
        res = self.specLoudness(eL)
        
                
        return self.mooreReturned(res, eL)
        


    def glasbergSpect(self, earSig : np.array ):
        
        filter = FilterComputer(self.kv, "1997", True)
        earSig = filter.FIR(earSig)
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
                
        
        
        return (compInt, compFq, oneHz)
        
        # print(compInt)
        
        #____

        

    def glasberg2002(self, inSig, fs = 32000 ):

        
        if(fs != self.kv['fs']):
            inSig = resample_poly(inSig, self.kv['fs'], fs)

        # print ( len(self._excitationPatern(inSig.compInt) ) )
        rawFftValues = self.glasbergSpect(inSig)
        class fftValuesFormat :
            def __init__(self, x , y, oneHz):
                self.compInt =   x
                self.compFq = y
                self.oneHz = oneHz
        
        fftVals = np.array([])
        for i in range(len(rawFftValues[0])):
            fftVals = np.append(fftVals, fftValuesFormat(rawFftValues[0][i], rawFftValues[1], rawFftValues[2]))
        

        l = len(fftVals)
        eLs = np.array([
            self._excitationPatern(fftValues = fftVals[i])
            for i in range(l) #Hopefully they all have the same len
        ])
        
        print("eLs has : ", eLs.shape)
        print(eLs)
        specLoud = np.array([
            self.getSpecLoudness(eL=eLs[i, :])
            for i in range(len(eLs))
        ])

        class result :
            def __init__(self, specLoud, STL, LTL, monoral, binaural):
                self.specLoud = specLoud
                self.STL = STL
                self.LTL = LTL
                self.monoral = monoral
                self.binaural = binaural

        monoralLoud = np.sum(specLoud, 1)* self.kv["erbStep"] 
        binauralLoud = 2*monoralLoud


        # STL and LTL:
        aSTL = 0.045
        rSTL = 0.02
        STL = [0.0] * len(binauralLoud)

        aLTL = 0.01
        rLTL = 0.0005
        LTL = [0.0] * len(binauralLoud)

        for ii in range(1, len(binauralLoud)):
            if binauralLoud[ii] > STL[ii - 1]:
                STL[ii] = float( aSTL * binauralLoud[ii] + (1 - aSTL) * STL[ii - 1])
            else:
                STL[ii] = float( rSTL * binauralLoud[ii] + (1 - rSTL) * STL[ii - 1] )

            if STL[ii] > LTL[ii - 1]:
                LTL[ii] = float( aLTL * STL[ii] + (1 - aLTL) * LTL[ii - 1] )
            else:
                LTL[ii] = float ( rLTL * STL[ii] + (1 - rLTL) * LTL[ii - 1] )


        return result(specLoud= specLoud, STL= STL, LTL = LTL, monoral= monoralLoud, binaural=binauralLoud)

        


if __name__ == '__main__':
    
    m = model(free = True)
    fs = 32000
    y = sineMaker.makeSine(400,0,2,100,fs) #it is crucial that both signals matlab/python have the same parameters (time is important)
    y[int(len(y) / 2): ] =  sineMaker.makeSine(300,0,2,80,fs)[: int( len(y) / 2)]
    g = m.glasberg2002(y, fs)
    print(f"monoral : {g.monoral[:10]}, \nbin :  {g.binaural[:10]} ,\n  STL : {g.STL[:10]} \n LTL : {g.LTL[:10]}")
    plt.plot(g.STL)
    plt.plot(g.LTL)
    plt.plot(g.binaural)
    plt.show()


#Notes : ajouter parametre timestep
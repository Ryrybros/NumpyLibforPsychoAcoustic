import numpy as np

class FFTComputer :

    def computeFFT(earSig : np.array, fs : float):
        class returned :
            def __init__(self, earSig : np.array, fs : float):
                self.spect = np.fft.fft(earSig)
                self.fftLen = len(self.spect)
                self.oneHz = (self.fftLen+2)/fs
                
                self.numBins = int(self.fftLen/(2+1))
                self.compInt =   2*abs((self.spect[0 : self.numBins] ** 2)/ (self.numBins*fs) )#2*abs( ( (self.spect(1:self.numBins))**2)/(self.numBins*fs)  )
                self.compFq = np.linspace(0,fs/2,self.numBins) # np.array([(2*self.numBins/fs)*i for i in range(int(fs/2))])# linspace(0,fs/2,numBins)

                print("compInt :  " ,len(self.compInt))
                self.nPoints = len(self.compFq)

        return returned(earSig,fs)
    

if __name__ == '__main__':
    n = 40000
    fs = n
    freq = 1500
    x = np.array([(freq*2*np.pi/n)*i for i in range(n)])
    y = np.sin(x)
    f = np.array([(1000/n)*i for i in range(n)])

    print(f"for arguments  : n = fs =  {n}  and  actual freq = {freq} : retuned values are :  ")
     
    res = FFTComputer.computeFFT(y,fs)
    import matplotlib.pyplot as plt
    troncRes = abs(res.spect)    
    print("the frequency after fft is  : " , np.argmax(troncRes))
    assert(abs(np.argmax(troncRes) - freq) < 10**(-3))

    # plt.plot(f, troncRes)
    # plt.show()

    print("compFq is : " , res.compFq)
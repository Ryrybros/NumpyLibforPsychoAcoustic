import numpy as np



class sineMaker:

    #This creates signals

    def makeSine(f,phase,time,dB,fs):
        
        t = np.linspace(0, time, int(fs * time), endpoint=False)
        sine_wave = ((2e-5) * 10**(dB/20)) * np.sin(2 * np.pi * f * t)

        return sine_wave
    
    #The following code is not useful anymore.
    

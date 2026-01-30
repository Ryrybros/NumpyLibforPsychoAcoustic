import numpy as np



class sineMaker:

    #This creates signals

    def makeSine(f,phase,time,dB,cutfreq):
        
        x = np.linspace(0,time, cutfreq)
        if(dB != 0):
            sine = ((20*(10**(-6)))*(10 ** (abs(dB)/20)))*np.sin((2*np.pi*f)*x + phase)
        else:
            sine = (20*(10**(-6))) *np.sin((2*np.pi*f)*x + phase)

        return sine
    
    #The following code is not useful anymore.
    

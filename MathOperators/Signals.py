import numpy as np



class sineMaker:

    

    def makeSine(f,phase,time,dB,cutfreq):
        x = np.array([ (time*(i/cutfreq)) for i in range(cutfreq) ])#0:(time/cutfreq):time)
        if(dB != 0):
            sine = ((10**(-3 - 17/20))*(10 ** (abs(dB)/20)))*np.sin((2*np.pi*f)*x + phase)
        else:
            sine = ((10^(-3 - 17/20)))*np.sin((2*np.pi*f)*x + phase)

        # plt.plot(x, sine)
        # plt.show()
        return sine
    


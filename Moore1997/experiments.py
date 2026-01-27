import numpy as np
from matplotlib import pyplot as plt
from Moore1997.model import model
from MathOperators.Signals import sineMaker
import time


def fig8Plot(dBLevel : float):
    mod = model(True)
    target_tQdB = [3.6, 6.3, 14.5, 20.2, 26.2]
    target_freq = [1000, 253, 108, 74, 52]
    
    fig, ax3 = plt.subplots(figsize=(7, 6))

    for i in range(len(target_freq)):
        
        # Start timer
        start_time = time.time()




        f = target_freq[i]
        
        y = sineMaker.makeSine(f, 0, 5, dBLevel, mod.kv['fs'])
        res = mod.moore1997(y)
        
        #The x axis is in dB
        X = 10 * np.log10(mod._eL) 
        
        
        Y = res.specLoudness

        
        ax3.semilogy(X, Y, color='black', linewidth=1.2)
        
        ax3.text(X[0], Y[0], f" {target_tQdB[i]}", verticalalignment='bottom')        
        
        end_time = time.time()

        # Calculate elapsed time
        elapsed_time = end_time - start_time
        print(f"Elapsed time for curve {i} plot: ", elapsed_time)


    
    ax3.set_xlabel('Excitation level, dB')
    ax3.set_ylabel("Specific loudness N' (log scale)")
    ax3.set_ylim(0.005, 50)
    ax3.set_xlim(0, 110)
    
    ax3.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax3.set_yticks([0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50])
    
    ax3.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()



def eLdBplot( ):
    mod = model(True)
    for i in np.linspace(20,100,9):
        start_time = time.time()

        sig = sineMaker.makeSine(1000,0,1,i,mod.kv['fs'])
        mod.moore1997(sig)
        plt.plot(mod.erbN,mod.results.eLdB)
        
        end_time = time.time()

        # Calculate elapsed time
        elapsed_time = end_time - start_time
        print(f"Elapsed time for step  : {i} dB : ", elapsed_time)
    plt.ylim(0,100)
    plt.show()

if __name__ == '__main__' :
    # fig8Plot(100)
    eLdBplot()
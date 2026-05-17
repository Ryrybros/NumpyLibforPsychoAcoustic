import numpy as np
from matplotlib import pyplot as plt
from Models.model import model
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
        # y += sineMaker.makeSine(f, 0, 5, dBLevel/10, mod.kv['fs'])

        res = mod.moore1997(y,target_tQdB[i])
        
        #The x axis is in dB
        X = 10 * np.log10(res.eL) 
        
        
        Y = res.Loudness.specLoudness

        
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

def fig12Plot(freq : float):
    mod = model(True)

    
    fig, ax3 = plt.subplots(figsize=(7, 6))

    X = np.linspace(20e-6, 120, 15)

    Y = []
    start_time = time.time()
    i = 0
    
    for x in X:
        y = sineMaker.makeSine(freq, 0, 2, x, mod.kv['fs'])
        
        res = mod.moore1997(y)
        
        Y.append(2*res.Loudness.monauralLoudness )
        i += 1 
        print(f"{i} / {len(X)} done")

        
            
        
    end_time = time.time()
    print("Time taken : ", end_time - start_time)


    ax3.semilogy(X, np.array(Y), color='black', linewidth=1.2)
    
    ax3.set_xlabel('Level, dB SPL')
    ax3.set_ylabel("Loudness  (log scale)")
    ax3.set_ylim(0.001, 500)
    ax3.set_xlim(0, 120)
    
    ax3.yaxis.set_major_formatter(plt.ScalarFormatter())
    ax3.set_yticks([0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50,100,200,500])
    
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







def timeStepExp(tStepVals: list):
    fs = 32000
    
    # Generate the signal once
    y = sineMaker.makeSine(400, 0, 2, 100, fs)
    y[int(len(y) / 2):] = sineMaker.makeSine(300, 0, 2, 80, fs)[:int(len(y) / 2)]
    
    # Total duration of the audio signal in seconds
    total_duration = len(y) / fs 

    print("Computing reference (timeStep = 0.001)...")
    ref_step = 0.001
    m_ref = model(free=True, timeStep=ref_step)
    g_ref = m_ref.glasberg2002(y, fs)
    
    time_ref = np.linspace(0, total_duration, len(g_ref.LTL))
    plt.plot(time_ref, g_ref.LTL, label=f"Reference (timeStep = {ref_step})", linestyle="--", color="black", alpha=0.7)
    
    for tStep in tStepVals:
        m = model(free=True, timeStep=tStep)
        
        print(f"Timer started for timeStep = {tStep}")
        t = time.time()
        g = m.glasberg2002(y, fs)
        print(f"Time for glasberg computation with timeStep = {tStep} : {time.time() - t:.4f}s")
        
        # Create a matching time vector for this specific step's output length
        time_current = np.linspace(0, total_duration, len(g.LTL))
        
        plt.plot(time_current, g.LTL, label=f"timeStep = {tStep}")
    
    #Styles
    plt.title("Glasberg 2002 LTL Comparison (Time-Aligned)")
    plt.xlabel("Time (seconds)")  # Changed from samples to seconds
    plt.ylabel("LTL")
    plt.legend(loc="best")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.tight_layout()
    plt.show()
    
    
    
if __name__ == '__main__' :
    # eLdBplot()
    # fig8Plot(105)
    # fig12Plot(1000)
    timeStepExp([0.0015, 0.0018,0.002,0.003, 0.005,0.01])

    
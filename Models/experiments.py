import numpy as np
from matplotlib import pyplot as plt
from os.path import dirname, join as pjoin
from scipy.io import wavfile
import scipy.io

from Models.model import model
from MathOperators.Signals import sineMaker
from MathOperators.ERBscale import ERB
import Models.equalSones as eqSone
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

    nb_pnt = 30
    
    fig, ax3 = plt.subplots(figsize=(7, 6))

    X = np.linspace(5, 60, nb_pnt)
    sound_pressure = (2e-5) * (10 ** (X/ 20))

    Y = []
    start_time = time.time()
    i = 0
    

    for x in sound_pressure:
        fs = 32000
        t = np.linspace(0, 0.5, int(fs * 0.5), endpoint=False)
        y = eqSone.setdBSPL( 10*np.log(x/2e-5), np.sin(2 * np.pi * freq * t))
        res = mod.moore1997(y)
        
        Y.append(2*res.Loudness.monauralLoudness )
        i += 1 
        print(f"{i} / {len(X)} done for dBSPL : {eqSone.getdBSPL(y)} added nan : {np.isnan(2*res.Loudness.monauralLoudness)} ")

        
            
        
    end_time = time.time()
    print("Time taken : ", end_time - start_time)

    pl = []
    x = []

    start = int(0)
    end = int(nb_pnt)

    for i in range(nb_pnt):
        if not np.isnan(Y[i]) :
            pl.append(Y[i])
            x.append(sound_pressure[i])

    print(f"len Y {len(Y)}, len(pl) {len(pl)}")
    log_x = x
    log_y = np.log(pl)

    z = np.polyfit(log_x[start:end], log_y[start:end], 1)
    alpha = z[0]

    print(f"Polynomial fit (slope, intercept): {z}")
    print(f"Calculated Stevens Exponent (Alpha): {alpha:.4f}")

    plt.figure(figsize=(8, 5))
    plt.plot(log_x[start:end], log_y[start:end], 'b-', label='Moore Model Data')
    plt.plot(log_x[start:end], log_y[start:end], 'ro', label='Fitted Region')

    # Calculate and plot the fit line
    fit_line = np.float64(z[0]) * np.float64(log_x[start:end]) + np.float64(z[1])
    plt.plot(log_x[start:end], fit_line, 'k--', label=f'Fit (Alpha = {alpha:.2f})')

    plt.xlabel('Intensity')
    plt.ylabel('ln(Monaural Loudness)')
    plt.title("Stevens' Power Law Exponent Estimation")
    plt.legend()
    plt.grid(True)
    plt.show()

    return Y



def spectrogram(sig : np.array, fs = 32000):
    m = model(True)
    glasb = m.glasberg2002(sig,fs)
    spect = glasb.specLoud
    print(spect.shape)
    print(len(spect[len(spect) - 1]))
    
    my_data = spect 

    num_time_steps = my_data.shape[0]  # This will read your 500ms dimension

    plot_data = my_data.T

    plt.figure(figsize=(10, 5))

    img = plt.imshow(plot_data, origin="lower", aspect="auto", cmap="viridis")

    plt.title("ERB-Band Spectrogram")
    plt.ylabel("ERB Frequency Bands (0 to 148)")
    plt.xlabel("Time (ms)")

    # --- 4. ACCURATELY MAP THE 500ms AXIS ---
    # This matches your exact number of matrix steps to the 0-500ms window
    plt.xticks(
        ticks=np.linspace(0, num_time_steps - 1, 6), 
        labels=np.linspace(0, 500, 6, dtype=int)
    )

    cbar = plt.colorbar(img)
    cbar.set_label("Intensity / Amplitude")

    plt.tight_layout()
    plt.show()
    import scipy.signal as signal

    # ==========================================
    # 1. SETUP / INPUT DATA
    # ==========================================
    fs = 32000          # Sampling rate in Hz
    duration = 0.5      # Duration in seconds
    t = np.linspace(0, duration, int(fs * duration), endpoint=False)

    # Let's create an example signal: a linear chirp (sweeping frequency)
    # Replace 'my_signal' with your own 1D NumPy audio array if you have one
    # my_signal = signal.chirp(t, f0=500, t1=duration, f1=4000, method='linear')
    my_signal = sig

    # ==========================================
    # 2. COMPUTE THE SPECTROGRAM
    # ==========================================
    # nperseg: window size (higher = better frequency resolution, worse time resolution)
    # noverlap: number of samples to overlap between windows (usually 50% to 75%)
    frequencies, times, Sxx = signal.spectrogram(
        my_signal, 
        fs=fs, 
        window='hann', 
        nperseg=1024, 
        noverlap=768
    )

    # ==========================================
    # 3. LOGARITHMIC SCALING (Convert to dB)
    # ==========================================
    # Human hearing perceives volume logarithmically. Converting the power spectral 
    # density (Sxx) to decibels makes quiet details visible alongside loud ones.
    Sxx_dB = 10 * np.log10(Sxx + 1e-10) # 1e-10 prevents log(0) errors

    # ==========================================
    # 4. PLOT THE SPECTROGRAM
    # ==========================================
    plt.figure(figsize=(10, 5))

    # pcolormesh is ideal for plotting raw coordinate grids (times, frequencies)
    img = plt.pcolormesh(times, frequencies, Sxx_dB, shading='gouraud', cmap='inferno')

    plt.title("Linear Frequency Spectrogram")
    plt.xlabel("Time (seconds)")
    plt.ylabel("Frequency (Hz)")

    # Limit the Y-axis to the human hearing range or the range of interest
    plt.ylim(0, 8000) 

    # Add colorbar calibrated to dB
    cbar = plt.colorbar(img)
    cbar.set_label("Magnitude (dB)")

    plt.tight_layout()
    plt.show()
    



def spectrogramExperiment(file):
    
    samplerate, data = wavfile.read(file)

    print(f"fs= {samplerate}")
    length = data.shape[0] / samplerate
    print(f"length = {length}s")

    y = data

    spectrogram(y)

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

    print("Computing reference ...")
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
    plt.xlabel("Time (seconds)")
    plt.ylabel("LTL")
    plt.legend(loc="best")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.tight_layout()
    plt.show()

    
    time_ref = np.linspace(0, total_duration, len(g_ref.STL))
    plt.plot(time_ref, g_ref.STL, label=f"Reference (timeStep = {ref_step})", linestyle="--", color="black", alpha=0.7)
    
    for tStep in tStepVals:
        m = model(free=True, timeStep=tStep)
        
        print(f"Timer started for timeStep = {tStep}")
        t = time.time()
        g = m.glasberg2002(y, fs)
        print(f"Time for glasberg computation with timeStep = {tStep} : {time.time() - t:.4f}s")
        
        # Create a matching time vector for this specific step's output length
        time_current = np.linspace(0, total_duration, len(g.STL))
        
        plt.plot(time_current, g.STL, label=f"timeStep = {tStep}")
    
    #Styles
    plt.title("Glasberg 2002 STL Comparison (Time-Aligned)")
    plt.xlabel("Time (seconds)")  # Changed from samples to seconds
    plt.ylabel("STL")
    plt.legend(loc="best")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.tight_layout()
    plt.show()
    


def equalLoudnessGraph(dBSPL_ref, f_lim):
    moore = model(True)
    f = 1000
    
    fs = 32000
    y = sineMaker.makeSine(f, 0,0.5,10,fs)

    y = eqSone.setdBSPL(dBSPL_ref, y)
    
    refSones = moore.moore1997(y,fs).Loudness.monauralLoudness
    print("referecne Loudness : " ,refSones)
    freqs = []
    curve = []

    
    test_axis= 2**(np.linspace(np.log2(10),np.log2(f_lim), 5))
    count = 0
    for i in test_axis:
        
        freqs.append(i)
        l = sineMaker.makeSine( i , 0,0.5,0,fs)
        l = eqSone.setdBSPL(70 , l)
        
        l_adapt = eqSone.getEqualSones(y,l,10, 100,0,400)
        curve.append(eqSone.getdBSPL(l_adapt))
        print("\n")

        print("dB is : ", curve[count])
        print(f"IsoCurve step : {count} ")
        count += 1

    from scipy.interpolate import CubicSpline

    axis = np.linspace(10,f_lim,200)
    cs_curve = CubicSpline(test_axis, curve)
    curve = cs_curve(axis)
    
    plt.plot(axis, curve, color='red', linewidth=2, label=f'Model Curve {dBSPL_ref}')
    plt.text(len(axis)/2, curve[int(len(axis)/2)], f"{dBSPL_ref}", 
         fontsize=15,  
         color='black', 
         va='center')
    


def allLoudnessGraph():
    f_lim = 10000

    plt.figure(figsize=(8, 6))

    for i in range(30,90,10):   equalLoudnessGraph(i,f_lim)
    fs = 32000

     
    plt.xscale('log')
    
    x_ticks = [16, 31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]
    x_labels = ['16', '31,5', '63', '125', '250', '500', '1K', '2K', '4K', '8K', '16K']
    plt.xticks(x_ticks, x_labels)
    plt.xlim(0, f_lim)
    
    plt.ylim(0, 130)
    
    plt.axvline(x=1000, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    
    plt.title("Reproduction de la courbe d'iso-sonie", fontsize=12, fontweight='bold')
    plt.xlabel("Fréquence Hz", fontsize=10)
    plt.ylabel("Niveau de pression (dB SPL)", fontsize=10)
    
    # 7. Add a clear, fine grid background
    plt.grid(True, which="both", linestyle="-", color='#d3d3d3', alpha=0.6)
    
    plt.tight_layout()
    plt.show()
    

if __name__ == '__main__' :
    # eLdBplot()
    # fig8Plot(105)
    # fig12Plot(3000)
    equalLoudnessGraph(70,10000)
    plt.show()
    # allLoudnessGraph()
    

    
    
    
        
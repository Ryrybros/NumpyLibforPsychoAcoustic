import numpy as np
from matplotlib import pyplot as plt
from Models.model import model
from MathOperators.Signals import sineMaker
from MathOperators.ERBscale import ERB
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

    nb_pnt = 60
    
    fig, ax3 = plt.subplots(figsize=(7, 6))

    X = np.linspace(0, 120, nb_pnt)

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

    pl = Y
    # 2. Generate the dB SPL X-axis and map it to Sound Pressure (Pa)
    # Note: The standard reference pressure for 0 dB SPL is 2e-5 Pa, not 1e-5.
    x_db = np.linspace(0, 120, len(pl))
    sound_pressure = (2e-5) * (10 ** (x_db / 20))

    # 3. Define your slicing indices for the fit
    # Ensure these indices align roughly with the 40 to 90 dB SPL region
    start = 25
    end = 50

    # 4. Transform BOTH variables to natural logs for the power-law fit
    log_x = np.log(sound_pressure)
    log_y = np.log(pl)

    # 5. Perform the linear fit on the log-log data
    # The slope 'z[0]' is directly equal to your Alpha exponent
    z = np.polyfit(log_x[start:end], log_y[start:end], 1)
    alpha = z[0]

    print(f"Polynomial fit (slope, intercept): {z}")
    print(f"Calculated Stevens Exponent (Alpha): {alpha:.4f}")

    # 6. Plot the log-log relationship to verify linearity
    plt.figure(figsize=(8, 5))
    plt.plot(log_x, log_y, 'b-', label='Moore Model Data')
    plt.plot(log_x[start:end], log_y[start:end], 'ro', label='Fitted Region')

    # Calculate and plot the fit line
    fit_line = z[0] * log_x[start:end] + z[1]
    plt.plot(log_x[start:end], fit_line, 'k--', label=f'Fit (Alpha = {alpha:.2f})')

    plt.xlabel('ln(Sound Pressure in Pa)')
    plt.ylabel('ln(Monaural Loudness)')
    plt.title("Stevens' Power Law Exponent Estimation")
    plt.legend()
    plt.grid(True)
    plt.show()

    return Y



def spectrogram(sig : np.array, fs = 32000):
    m = model(True)
    spect = m.glasberg2002(sig,fs).specLoud
    print(spect.shape)
    print(len(spect[len(spect) - 1]))
    # --- 1. ROUTE YOUR ACTUAL DATA HERE ---
    # Replace 'your_matrix_variable' with the name of your real matrix
    my_data = spect 

    # Grab the actual number of time steps from your matrix shape
    num_time_steps = my_data.shape[0]  # This will read your 500ms dimension

    # --- 2. TRANSPOSE FOR SPECTROGRAM ALIGNMENT ---
    # Flips it from (time, bands) to (bands, time) so bands stack vertically
    plot_data = my_data.T

    # --- 3. CREATE THE GRAPH ---
    plt.figure(figsize=(10, 5))

    # origin="lower" puts low frequency bands at the bottom
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
    


    
    
if __name__ == '__main__' :
    # eLdBplot()
    # # fig8Plot(105)
    time = 0.5
    fs = 32000
    y = sineMaker.makeSine(100,0,time,50,fs)
    x = np.linspace(0,1,int(fs*time))
    for i in range(10):
        y += np.exp(x*( 5 - i )) * sineMaker.makeSine(i*500,0,0.5,50,32000)
    norm = np.sqrt(np.sum(y*y))
    y /= norm
    
    spectrogram(y)
        
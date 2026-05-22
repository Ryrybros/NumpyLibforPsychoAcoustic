import numpy as np



class sineMaker:

    #This creates signals

    def makeSine(f,phase,time,dB,fs):


        # calibration_dbspl = 94
        # target_dbfs_rms = dB - calibration_dbspl  

        # # Convert dBFS RMS to linear RMS, then to Peak Amplitude
        # rms_amplitude = 10 ** (target_dbfs_rms / 20)
        # peak_amplitude = rms_amplitude * np.sqrt(2)

        t = np.linspace(0, time, int(fs * time), endpoint=False)
        sine_wave = ((2e-5) * 10**(dB/20)) * np.sin(2 * np.pi * f * t)

        # 3. Verify the RMS of the generated signal
        generated_rms = np.sqrt(np.mean(sine_wave**2))
        generated_dbfs_rms = 20 * np.log10(generated_rms)
        return sine_wave
    
    #The following code is not useful anymore.
    

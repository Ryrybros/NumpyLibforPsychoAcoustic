import numpy as np
from matplotlib import pyplot as plt
from scipy.io import wavfile
import scipy.io

from Models.model import model
from MathOperators.Signals import sineMaker
from MathOperators.ERBscale import ERB
import time
import random

def getSig(filename: str):
    samplerate, data = wavfile.read(filename)

    if data.dtype == np.int16:
        data = data.astype(np.float64) / 32768.0
    elif data.dtype == np.int32:
        data = data.astype(np.float64) / 2147483648.0

    print(f"fs= {samplerate}")
    length = data.shape[0] / samplerate
    print(f"length = {length}s")

    return data

def save_as_wav(sig: np.array, destination: str):
    fs = 32000
    
    clipped_sig = np.clip(sig, -1.0, 1.0)
    
    int16_sig = (clipped_sig * 32767.0).astype(np.int16)
    
    wavfile.write(destination, fs, int16_sig)
    print(f"Successfully saved signal to {destination} (fs={fs}Hz, 16-bit PCM)")


def getdBSPL(sig: np.array):
    """
    Calculates the global dB SPL level of an input digital signal.
    Assumes standard calibration where a digital full-scale peak amplitude 
    of 1.0 (RMS = 1/sqrt(2)) maps to 94 dB SPL (1 Pascal RMS).
    """
    # Calculate Root-Mean-Square (RMS) of the signal
    rms = np.sqrt(np.mean(sig**2))
    
    if rms == 0:
        return -np.inf
        
    # Reference digital RMS corresponding to 1 Pa (94 dB SPL)
    # A full-scale sine wave with peak=1.0 has an RMS of 1/sqrt(2)
    p_ref_digital = 1.0 / np.sqrt(2.0)
    
    # Calculate dB SPL relative to the calibration point
    db_spl = 94.0 + 20.0 * np.log10(rms / p_ref_digital)
    return db_spl


def setdBSPL(target_dBSPL: float, sig: np.array):
    
    current_dBSPL = getdBSPL(sig)
    
    gain_db = target_dBSPL - current_dBSPL
    
    linear_factor = 10 ** (gain_db / 20.0)
    
    return sig * linear_factor


def getSones(mod : model, sig : np.array):
    #Monoral Loudness but binaural won't change the experiment
    Loudness = mod.moore1997(sig).Loudness.monauralLoudness
    return Loudness

def getEqualSones(ref_sig: np.array, sig: np.array, limit=20, accuracy = 1000 ,low_dB = 0, high_dB = 120):
    # return (1/getdBSPL(sig))*random.randint(5,15)*sig
    moore = model(True)
    
    target_Sones = getSones(moore, ref_sig)
    
    if target_Sones == 0:
        return np.zeros_like(sig)

    n = 0
    success = False
    res = sig.copy()

    print(f"Target Loudness: {target_Sones:.4f} Sones")
    print(f"Target ref dB: {getdBSPL(ref_sig):.4f} dB")
    print(f"Original dB: {getdBSPL(sig):.4f} dB")

    while n < limit:
        n += 1
        # Try the midpoint dB between our current boundaries
        mid_dB = (low_dB + high_dB) / 2.0
        
        res = setdBSPL(mid_dB, sig)
        
        current_Sones = getSones(moore, res)
        
        error = current_Sones - target_Sones
        
        print(f"Step: {n:2d} | Test dB: {mid_dB:6.2f} | Sones: {current_Sones:6.2f} | Error: {abs(error):.4f}  ", end="\r", flush=True)
        
        if abs(error) < target_Sones/accuracy:
            success = True
            break
            
        if error <= 0:
            low_dB = mid_dB
        else:
            high_dB = mid_dB

    print()
    
    if not success:
        print(f"!!! Could not fully converge within {limit} steps. Best estimation returned.")
    else:
        print(f"Success! Matched target loudness at {current_Sones:.2f} Sones (Final Level: {mid_dB:.1f} dB SPL).")
        
    return res


if __name__ == '__main__':
    y = getSig("Models/sounds/import_Prog/Dep_Am.wav")
    # print(getdBSPL(y))
    y2 = getSig("Models/sounds/import_Prog/ODep_Am.wav")
    # print("Now : ", getdBSPL(y))
    # save_as_wav(y2, "Models/sounds/exp/test_0.wav")
    moore = model(True)
    eqSone = getEqualSones(y, y2, 100)
    print(getSones(moore, eqSone))
    save_as_wav(eqSone, "Models/sounds/exp/testres_Prog.wav")
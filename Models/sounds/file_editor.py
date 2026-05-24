import wave as w
import matplotlib.pyplot as plt

import os
import subprocess
import numpy as np
from scipy.io import wavfile
def trim_and_play_wav(input_path, output_path, duration_seconds):
    """
    Opens a WAV file using scipy, trims it, saves it, and plays it via system aplay.
    """
    if not os.path.exists(input_path):
        print(f"Error: The file '{input_path}' could not be found.")
        return

    # 1. Open the WAV file (Returns sampling rate and numpy data array)
    print(f"Loading '{input_path}'...")
    fs, data = wavfile.read(input_path)
    
    # Calculate how many samples correspond to the target duration
    target_samples = int(duration_seconds * fs)
    
    # 2. Trim the audio safely using array slicing
    if len(data) > target_samples:
        trimmed_data = data[:target_samples]
        print(f"Trimmed file to {duration_seconds} seconds.")
    else:
        trimmed_data = data
        print(f"Warning: File is only {len(data)/fs:.2f}s long. Kept original length.")

    # 3. Save the trimmed file back to disk
    wavfile.write(output_path, fs, trimmed_data)
    print(f"Saved trimmed audio to '{output_path}'")
    
    
if __name__ == "__main__":
    # Replace these with your actual file paths
    input_wav = "my_experiment_sound.wav"
    output_wav = "trimmed_output.wav"
    cut_duration = 1.0  # Duration in seconds
    
    trim_and_play_wav("Models/sounds/Dep_Chord.wav", "Models/sounds/Dep_Chord_stationnary.wav", 0.5)
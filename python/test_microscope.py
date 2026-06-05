import librosa
import numpy as np
import fission_core
import matplotlib.pyplot as plt
import os

def test_microscope():
    # Load example audio
    audio_path = "/root/listening/mp3/example_one.mp3"
    if not os.path.exists(audio_path):
        print(f"Error: {audio_path} not found.")
        return
        
    y, sr = librosa.load(audio_path, sr=22050, mono=True)
    # Take a 1-second segment
    y = y[sr:2*sr]
    
    print(f"Original audio shape: {y.shape}")
    
    # Apply Quantum Sound Microscope (stretch by 8x)
    ws = 8192
    amount = 8.0
    y_stretched = fission_core.quantum_scope(y, amount, ws)
    
    print(f"Stretched audio shape: {y_stretched.shape}")
    
    # Compare Spectrums
    # Original (on same window size for fair comparison)
    Y_orig = np.abs(np.fft.rfft(y[:ws] * np.hanning(ws)))
    # Stretched (take a segment from the middle)
    mid = len(y_stretched) // 2
    Y_stretch = np.abs(np.fft.rfft(y_stretched[mid:mid+ws] * np.hanning(ws)))
    
    freqs = np.fft.rfftfreq(ws, 1/sr)
    
    plt.figure(figsize=(12, 6))
    plt.subplot(2, 1, 1)
    plt.plot(freqs, 20*np.log10(Y_orig + 1e-6), label="Original (1st window)")
    plt.xlim(0, 4000)
    plt.title("Original vs Stretched Spectrum (Phase-Randomized Microscope)")
    plt.legend()
    
    plt.subplot(2, 1, 2)
    plt.plot(freqs, 20*np.log10(Y_stretch + 1e-6), label=f"Stretched ({amount}x)", color='orange')
    plt.xlim(0, 4000)
    plt.xlabel("Frequency (Hz)")
    plt.legend()
    
    output_plot = "/root/fission/python/microscope_test.png"
    plt.tight_layout()
    plt.savefig(output_plot)
    print(f"QA Checkpoint: Microscope plot saved to {output_plot}")

if __name__ == "__main__":
    test_microscope()

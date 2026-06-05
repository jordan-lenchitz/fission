import librosa
import numpy as np
import fission_core
import matplotlib.pyplot as plt
import os

def test_tiers():
    # Load example audio
    audio_path = "/root/listening/mp3/example_one.mp3"
    if not os.path.exists(audio_path):
        print(f"Error: {audio_path} not found.")
        return
        
    y, sr = librosa.load(audio_path, sr=22050, mono=True)
    y = y[sr:2*sr] # 1s segment
    
    # Apply Microscope
    ws = 8192
    amount = 4.0
    y_stretched = fission_core.quantum_scope(y, amount, ws)
    
    # Get Tiers
    n_fft = 2**19
    mid = len(y_stretched) // 2
    y_slice = np.zeros(n_fft, dtype=np.float32)
    raw_slice = y_stretched[mid:mid+n_fft]
    y_slice[:len(raw_slice)] = raw_slice.astype(np.float32)
    
    tier1, tier2 = fission_core.get_fission_tiers(y_slice, float(sr), n_fft)
    
    print(f"Detected {len(tier1)} Tier 1 candidates (Independent Pitches).")
    for c in tier1:
        print(f"Tier 1: {c.frequency:.2f} Hz, {c.power:.2f} dB")
        
    print(f"Detected {len(tier2)} Tier 2 candidates (Harmonic Color).")
    for c in tier2:
        print(f"Tier 2: {c.frequency:.2f} Hz, {c.power:.2f} dB")

    # Plot
    Y = np.abs(np.fft.rfft(y_slice * np.hanning(len(y_slice)))) / n_fft
    Y_db = 20 * np.log10(Y + 1e-12) + 120
    freqs = np.fft.rfftfreq(n_fft, 1/sr)
    
    plt.figure(figsize=(14, 7))
    plt.plot(freqs, Y_db, alpha=0.3, color='gray', label="Microscoped Spectrum")
    
    if tier1:
        t1_freqs = [c.frequency for c in tier1]
        t1_powers = [c.power for c in tier1]
        plt.scatter(t1_freqs, t1_powers, color='green', marker='*', s=100, label="Tier 1 (Fission)")
        
    if tier2:
        t2_freqs = [c.frequency for c in tier2]
        t2_powers = [c.power for c in tier2]
        plt.scatter(t2_freqs, t2_powers, color='red', marker='+', s=100, label="Tier 2 (Timbre)")
        
    plt.xlim(20, 4000)
    plt.ylim(0, 130)
    plt.title("Fission Tiers: Independent Pitches vs. Spectral Color")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power (dB + 120)")
    plt.legend()
    plt.grid(True, alpha=0.2)
    
    output_plot = "/root/fission/python/fission_tiers_test.png"
    plt.savefig(output_plot)
    print(f"Final QA Checkpoint: Fission Tiers plot saved to {output_plot}")

if __name__ == "__main__":
    test_tiers()

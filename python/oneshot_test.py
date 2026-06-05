import librosa
import os
import numpy as np
import fission_core
import matplotlib.pyplot as plt

def oneshot_analysis(file_path, start_time=51.0):
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return

    sr = 22050
    # Load from 51s to the end
    y, _ = librosa.load(file_path, sr=sr, offset=start_time)
    duration = len(y) / sr
    print(f"Loaded segment: {start_time}s to {start_time + duration:.2f}s ({duration:.2f}s total)")

    stretch_amount = 8.0 # Higher stretch for better resolution in decay
    ws = 8192
    n_fft = 2**19
    
    print(f"Applying Quantum Sound Microscope ({stretch_amount}x stretch)...")
    y_stretched = fission_core.quantum_scope(y, stretch_amount, ws)
    
    # We'll analyze two points: 
    # 1. Steady state (shortly after 51s)
    # 2. Decay phase (near the end)
    
    # Points in the stretched audio
    # The stretched audio length is roughly duration * sr * stretch_amount
    len_stretched = len(y_stretched)
    point_steady = int(sr * 1.0 * stretch_amount) # 1s into segment
    point_decay = int(sr * (duration - 2.0) * stretch_amount) # 2s before end
    
    analysis_points = [
        ("Steady State (~52s)", point_steady),
        ("Decay Phase (~59s)", point_decay)
    ]
    
    plt.figure(figsize=(15, 10))
    
    for i, (label, start_idx) in enumerate(analysis_points):
        print(f"Analyzing {label}...")
        
        # Ensure we don't go out of bounds
        if start_idx + n_fft > len_stretched:
            start_idx = len_stretched - n_fft - 1
        
        y_slice = np.zeros(n_fft, dtype=np.float32)
        raw_slice = y_stretched[start_idx : start_idx + n_fft]
        y_slice[:len(raw_slice)] = raw_slice.astype(np.float32)
        
        tier1, tier2 = fission_core.get_fission_tiers(y_slice, float(sr), n_fft)
        
        # Plotting
        Y = np.abs(np.fft.rfft(y_slice * np.hanning(len(y_slice)))) / n_fft
        Y_db = 20 * np.log10(Y + 1e-12) + 120
        freqs = np.fft.rfftfreq(n_fft, 1/sr)
        
        ax = plt.subplot(2, 1, i+1)
        ax.plot(freqs, Y_db, alpha=0.3, color='gray', label="Spectrum")
        
        if tier1:
            t1_f = [c.frequency for c in tier1]
            t1_p = [c.power for c in tier1]
            ax.scatter(t1_f, t1_p, color='green', marker='*', s=100, label="Tier One")
            for c in tier1:
                ax.annotate(f"{c.frequency:.0f}", (c.frequency, c.power), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8, color='green')
        
        if tier2:
            t2_f = [c.frequency for c in tier2]
            t2_p = [c.power for c in tier2]
            ax.scatter(t2_f, t2_p, color='red', marker='+', s=100, label="Tier Two")
            
        ax.set_xlim(20, 4000)
        ax.set_ylim(0, 130)
        ax.set_title(f"One-Shot Analysis: {label}")
        ax.set_ylabel("Power (dB + 120)")
        ax.legend()
        ax.grid(True, alpha=0.2)
        
    plt.xlabel("Frequency (Hz)")
    plt.tight_layout()
    
    output_path = "/root/fission/python/oneshot_decay_test.png"
    plt.savefig(output_path)
    print(f"One-shot analysis complete. Results saved to {output_path}")

if __name__ == "__main__":
    oneshot_analysis("/root/listening/mp3/example_one.mp3")

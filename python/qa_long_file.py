import librosa
import os
import numpy as np
from fission import FissionModel
import matplotlib.pyplot as plt

def analyze_long_file(file_path, segment_duration=10.0):
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return

    # Initialize model
    model = FissionModel()
    
    # Get total duration
    duration = librosa.get_duration(path=file_path)
    print(f"Total duration: {duration:.2f} seconds.")

    # Create output directory for plots
    output_dir = "/root/fission/python/qa_results"
    os.makedirs(output_dir, exist_ok=True)

    # Process in segments
    for start_time in np.arange(0, duration, segment_duration):
        end_time = min(start_time + segment_duration, duration)
        if end_time - start_time < 1.0: # Skip very short tail
            continue
            
        print(f"\n--- Analyzing Segment: {start_time:.2f}s to {end_time:.2f}s ---")
        
        # Load specific segment
        y, sr = librosa.load(file_path, sr=22050, offset=start_time, duration=segment_duration)
        
        # We need a temporary file for the segment or we can modify FissionModel to take 'y' directly
        # For now, let's just use the segment audio array by modifying a small part of the logic
        # But wait, fission.py's analyze takes a path. Let's fix fission.py to be more flexible first.
        
        # [Self-Correction]: I'll just write a quick wrapper here that calls the core functions directly 
        # since I've already exposed them in fission_core.
        
        import fission_core
        stretch_amount = 4.0
        ws = 8192
        n_fft = 2**19
        
        print(f"  [1/3] Applying Microscope...")
        y_stretched = fission_core.quantum_scope(y, stretch_amount, ws)
        
        mid = len(y_stretched) // 2
        y_slice = np.zeros(n_fft, dtype=np.float32)
        raw_slice = y_stretched[mid:mid+n_fft]
        y_slice[:len(raw_slice)] = raw_slice.astype(np.float32)
        
        print(f"  [2/3] Extracting Tiers...")
        tier1, tier2 = fission_core.get_fission_tiers(y_slice, 22050.0, n_fft)
        
        print(f"  [3/3] Found {len(tier1)} Tier 1 (Pitches), {len(tier2)} Tier 2 (Timbre)")
        
        # Plotting
        Y = np.abs(np.fft.rfft(y_slice * np.hanning(len(y_slice)))) / n_fft
        Y_db = 20 * np.log10(Y + 1e-12) + 120
        freqs = np.fft.rfftfreq(n_fft, 1/22050)
        
        plt.figure(figsize=(14, 7))
        plt.plot(freqs, Y_db, alpha=0.3, color='gray', label="Microscoped Spectrum")
        
        if tier1:
            t1_f = [c.frequency for c in tier1]
            t1_p = [c.power for c in tier1]
            plt.scatter(t1_f, t1_p, color='green', marker='*', s=100, label="Tier 1 (Fission)")
            for c in tier1:
                plt.annotate(f"{c.frequency:.0f}Hz", (c.frequency, c.power), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8, color='green')
        
        if tier2:
            t2_f = [c.frequency for c in tier2]
            t2_p = [c.power for c in tier2]
            plt.scatter(t2_f, t2_p, color='red', marker='+', s=100, label="Tier 2 (Timbre)")
            
        plt.xlim(20, 4000)
        plt.ylim(0, 130)
        plt.title(f"Fission Analysis: {start_time:.1f}s - {end_time:.1f}s")
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power (dB + 120)")
        plt.legend()
        plt.grid(True, alpha=0.2)
        
        plot_path = os.path.join(output_dir, f"segment_{int(start_time):03d}.png")
        plt.savefig(plot_path)
        plt.close()
        print(f"  [Done] Plot saved: {plot_path}")

if __name__ == "__main__":
    analyze_long_file("/root/listening/mp3/example_one.mp3")

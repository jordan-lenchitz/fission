import librosa
import numpy as np
import fission_core
import matplotlib.pyplot as plt
import os

def test_candidates():
    # Load example audio
    audio_path = "/root/listening/mp3/example_one.mp3"
    if not os.path.exists(audio_path):
        print(f"Error: {audio_path} not found.")
        return
        
    y, sr = librosa.load(audio_path, sr=22050, mono=True)
    # Take a 1-second segment
    y = y[sr:2*sr]
    
    # Apply Microscope
    ws = 8192
    amount = 4.0
    y_stretched = fission_core.quantum_scope(y, amount, ws)
    
    # Get candidates from the middle of the stretched segment
    n_fft = 2**19 # High resolution
    mid = len(y_stretched) // 2
    if mid + n_fft > len(y_stretched):
        mid = 0
    
    y_slice_raw = y_stretched[mid:mid+n_fft]
    y_slice = np.zeros(n_fft, dtype=np.float32)
    y_slice[:len(y_slice_raw)] = y_slice_raw.astype(np.float32)
    
    candidates = fission_core.get_fission_candidates(y_slice, float(sr), n_fft)
    
    print(f"Detected {len(candidates)} candidates.")
    for c in candidates[:5]:
        print(f"Candidate: {c.frequency:.2f} Hz, Power: {c.power:.2f} dB, ERB: {c.erb:.2f}")

    # Plot
    # Recompute spectrum for plotting
    Y = np.abs(np.fft.rfft(y_slice * np.hanning(len(y_slice)))) / n_fft
    Y_db = 20 * np.log10(Y + 1e-12) + 120
    freqs = np.fft.rfftfreq(n_fft, 1/sr)
    
    plt.figure(figsize=(12, 6))
    plt.plot(freqs, Y_db, alpha=0.5, label="Microscoped Spectrum")
    
    cand_freqs = [c.frequency for c in candidates]
    cand_powers = [c.power for c in candidates]
    plt.scatter(cand_freqs, cand_powers, color='red', marker='+', label="Adaptive Candidates")
    
    plt.xlim(20, 4000)
    plt.ylim(0, 130)
    plt.title("Adaptive Peakfinding on Microscoped Spectrum")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power (dB + 120)")
    plt.legend()
    
    output_plot = "/root/fission/python/candidates_test.png"
    plt.savefig(output_plot)
    print(f"QA Checkpoint: Candidates plot saved to {output_plot}")

if __name__ == "__main__":
    test_candidates()

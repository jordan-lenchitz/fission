import fission_core
import librosa
import numpy as np
import matplotlib.pyplot as plt

class FissionModel:
    """
    Repertoire-Agnostic Model of Spectral Fission (2026 Edition)
    
    A high-performance implementation of the Lenchitz (2021) model.
    Uses a Rust-based Quantum Sound Microscope for hyper-resolution 
    spectral analysis and adaptive peakfinding.
    """
    def __init__(self, sample_rate=22050, n_fft=2**19, microscope_ws=8192):
        self.sr = sample_rate
        self.n_fft = n_fft
        self.ws = microscope_ws

    def analyze(self, audio_path, stretch_amount=4.0):
        """
        Analyze a sustained sonority for spectral fission candidates.
        """
        # Load
        y, sr = librosa.load(audio_path, sr=self.sr, mono=True)
        
        # Apply Microscope
        print(f"Applying Quantum Sound Microscope ({stretch_amount}x stretch)...")
        y_stretched = fission_core.quantum_scope(y, stretch_amount, self.ws)
        
        # Take the most stable segment (middle)
        mid = len(y_stretched) // 2
        y_slice = np.zeros(self.n_fft, dtype=np.float32)
        raw_slice = y_stretched[mid:mid+self.n_fft]
        y_slice[:len(raw_slice)] = raw_slice.astype(np.float32)
        
        # Get Tiers
        print("Extracting Fission Tiers...")
        tier1, tier2 = fission_core.get_fission_tiers(y_slice, float(self.sr), self.n_fft)
        
        return AnalysisResult(tier1, tier2, y_slice, self.sr, self.n_fft)

class AnalysisResult:
    def __init__(self, tier1, tier2, spectrum_slice, sr, n_fft):
        self.tier1 = tier1
        self.tier2 = tier2
        self.spectrum_slice = spectrum_slice
        self.sr = sr
        self.n_fft = n_fft

    def plot(self, x_lim=(20, 4000), save_path=None):
        Y = np.abs(np.fft.rfft(self.spectrum_slice * np.hanning(len(self.spectrum_slice)))) / self.n_fft
        Y_db = 20 * np.log10(Y + 1e-12) + 120
        freqs = np.fft.rfftfreq(self.n_fft, 1/self.sr)
        
        plt.figure(figsize=(14, 7))
        plt.plot(freqs, Y_db, alpha=0.3, color='gray', label="Microscoped Spectrum")
        
        if self.tier1:
            t1_f = [c.frequency for c in self.tier1]
            t1_p = [c.power for c in self.tier1]
            plt.scatter(t1_f, t1_p, color='green', marker='*', s=100, label="Tier One")
            
        if self.tier2:
            t2_f = [c.frequency for c in self.tier2]
            t2_p = [c.power for c in self.tier2]
            plt.scatter(t2_f, t2_p, color='red', marker='+', s=100, label="Tier Two")
            
        plt.xlim(x_lim)
        plt.ylim(0, 130)
        plt.title("Spectral Fission Analysis")
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Power (dB + 120)")
        plt.legend()
        plt.grid(True, alpha=0.2)
        
        if save_path:
            plt.savefig(save_path)
            print(f"Plot saved to {save_path}")
        else:
            plt.show()

if __name__ == "__main__":
    # Quick CLI usage example
    import sys
    if len(sys.argv) > 1:
        model = FissionModel()
        result = model.analyze(sys.argv[1])
        print(f"Tier 1: {len(result.tier1)} candidates found.")
        print(f"Tier 2: {len(result.tier2)} candidates found.")
        result.plot(save_path="fission_result.png")

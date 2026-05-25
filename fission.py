import numpy as np
import librosa
import scipy.signal
import matplotlib.pyplot as plt
import argparse
import sys

def hz_to_erb(hz):
    """Convert frequency from Hz to ERB scale (Glasberg and Moore, 1990)."""
    return 21.4 * np.log10(1 + 0.00437 * hz)

def erb_to_hz(erb):
    """Convert frequency from ERB scale to Hz."""
    return (10**(erb / 21.4) - 1) / 0.00437

def optimize_windowsize(ws):
    """Optimize window size for FFT (next power of 2)."""
    return 2**int(np.ceil(np.log2(ws)))

def scope(y, sr, amount):
    """
    Overlap-add phase vocoder for time stretching with phase randomization.
    Mimics the behavior of the MATLAB scope.m function.
    """
    if amount == 1:
        return y, sr
    
    ws = 7324
    ws = optimize_windowsize(ws)
    half_ws = ws // 2
    
    # Fade out end
    end_size = int(sr * 0.05)
    if end_size < 16:
        end_size = 16
    if len(y) > end_size:
        y[-end_size:] *= np.linspace(1, 0, end_size)
        
    nsamples = len(y)
    start_pos = 0
    displace_pos = (ws * 0.5) / amount
    
    # Hann-like window from MATLAB: (1-linspace(-1, 1, ws).^2).^1.25'
    window = (1 - np.linspace(-1, 1, ws)**2)**1.25
    
    old_windowed_buf = np.zeros(ws)
    stretched = []
    
    while True:
        istart_pos = int(np.floor(start_pos))
        if istart_pos + ws > nsamples:
            buf = np.zeros(ws)
            available = nsamples - istart_pos
            if available > 0:
                buf[:available] = y[istart_pos:]
        else:
            buf = y[istart_pos:istart_pos+ws]
            
        buf = buf * window
        
        # FFT
        freqs = np.fft.fft(buf)
        mags = np.abs(freqs)
        
        # Randomize phase
        ph = 2 * np.pi * np.random.rand(len(freqs))
        freqs_rand = mags * np.exp(1j * ph)
        
        # IFFT
        buf_out = np.fft.ifft(freqs_rand).real
        buf_out = buf_out * window
        
        # Overlap-add
        output = buf_out[:half_ws] + old_windowed_buf[half_ws:]
        old_windowed_buf = buf_out
        stretched.extend(output)
        
        start_pos += displace_pos
        if np.floor(start_pos) + ws > nsamples:
            break
            
    stretched = np.array(stretched)
    # Normalize
    if len(stretched) > 0:
        stretched = stretched / np.max(np.abs(stretched)) * (1 - np.finfo(float).eps)
        
    return stretched, sr

def nonerb_positive_pow(file_path, k=19):
    """
    Compute power spectrum with high frequency resolution.
    Matches nonerb_positive_pow.m logic.
    """
    y, sr = librosa.load(file_path, sr=None, mono=True)
    N = 2**k
    
    ns = len(y)
    amount = 1
    while ns < N:
        ns *= 2
        amount *= 2
        
    if amount != 1:
        amount //= 2
        y, sr = scope(y, sr, amount)
        
    # Compute FFT
    # MATLAB: c = fft(y, N)/N; 
    # Python: np.fft.fft returns unnormalized.
    c = np.fft.fft(y, n=N) / N
    
    # Frequencies 
    f = np.arange(1, N//2) * sr / N
    
    # Power in positive range (0 to 120)
    # MATLAB: p = 20*log10(abs(c(2:N/2))); p(p<-120) = -120; b = p+120;
    p = 20 * np.log10(np.abs(c[1:N//2]) + 1e-15)
    p = np.maximum(p, -120)
    
    # Truncate
    mask = (f > 20) & (f < 4000)
    f = f[mask]
    p = p[mask]
    
    return f, p + 120

def get_fission_candidates(file_path):
    """Core logic for spectral fission candidates detection."""
    dB_diff = 6
    ERB_rad = 2
    floor_hz = 1000
    
    f, b = nonerb_positive_pow(file_path)
    
    # Peak finding
    # Use scipy.signal.find_peaks with some minimal height or prominence
    # to match the "noise tolerant" peakfinder.m
    # peakfinder.m default sel is (max-min)/4. 
    sel = (np.max(b) - np.min(b)) / 4
    locs, _ = scipy.signal.find_peaks(b, prominence=sel)
    
    if len(locs) == 0:
        return [], [], f, b
    
    ap = f[locs]
    bp = b[locs]
    
    aerb = hz_to_erb(f)
    aperb = hz_to_erb(ap)
    
    # Sort peaks by power descending
    indices = np.argsort(bp)[::-1]
    sorted_aperb = aperb[indices]
    sorted_bp = bp[indices]
    
    # Filter by floor
    mask = sorted_aperb > hz_to_erb(floor_hz)
    pre_candidates_erb = sorted_aperb[mask]
    pre_candidates_pow = sorted_bp[mask]
    
    if len(pre_candidates_erb) == 0:
        return [], [], aerb, b
    
    # Find nearest indices in aerb for candidates
    # (equivalent to nearest_index in MATLAB)
    pre_cand_inds = [np.argmin(np.abs(aerb - val)) for val in pre_candidates_erb]
    
    # Left dominance check
    dB_results_pt1 = []
    for i in range(len(pre_candidates_erb)):
        working_ind = pre_cand_inds[i]
        inc_ind = working_ind - 1
        
        if inc_ind < 0 or b[inc_ind] > pre_candidates_pow[i]:
            continue
        
        while inc_ind > 0 and b[inc_ind] < pre_candidates_pow[i] + dB_diff:
            inc_ind -= 1
            
        diff_left = pre_candidates_erb[i] - aerb[inc_ind]
        if diff_left > ERB_rad:
            dB_results_pt1.append([pre_candidates_erb[i], pre_candidates_pow[i], diff_left])
            
    if not dB_results_pt1:
        return [], [], aerb, b
    
    dB_results_pt1 = np.array(dB_results_pt1)
    
    # Right dominance check
    raw_cands = []
    mid_cand_inds = [np.argmin(np.abs(aerb - val)) for val in dB_results_pt1[:, 0]]
    
    for i in range(len(dB_results_pt1)):
        working_ind = mid_cand_inds[i]
        inc_ind = working_ind + 1
        
        if inc_ind >= len(b) or b[inc_ind] > dB_results_pt1[i, 1]:
            continue
            
        while inc_ind < len(b) - 1 and b[inc_ind] < dB_results_pt1[i, 1] + dB_diff:
            inc_ind += 1
            
        diff_right = aerb[inc_ind] - dB_results_pt1[i, 0]
        if diff_right > ERB_rad:
            raw_cands.append([dB_results_pt1[i, 0], dB_results_pt1[i, 1], diff_right])
            
    if not raw_cands:
        return [], [], aerb, b
        
    raw_cands = np.array(raw_cands)
    
    # Cluster reduction
    # Sort by frequency
    sorted_idx = np.argsort(raw_cands[:, 0])
    sorted_lrcs = raw_cands[sorted_idx]
    
    cluster_free = []
    working_cluster = []
    x = 0
    while x < len(sorted_lrcs):
        working_cluster.append(sorted_lrcs[x])
        if x + 1 < len(sorted_lrcs) and (sorted_lrcs[x+1, 0] / sorted_lrcs[x, 0] <= 1.01):
            x += 1
            continue
        else:
            # Sort cluster by power descending and take the top
            working_cluster = sorted(working_cluster, key=lambda val: val[1], reverse=True)
            cluster_free.append(working_cluster[0][:2])
            x += 1
            working_cluster = []
            
    cluster_free = np.array(cluster_free)
    
    # Masking among candidates
    candidates_0 = cluster_free[np.argsort(cluster_free[:, 1])[::-1]] # Sort by power desc
    candidates_1 = []
    
    for n in range(len(candidates_0)):
        working_table = np.delete(candidates_0, n, axis=0)
        is_masked = False
        for nn in range(len(working_table)):
            if abs(candidates_0[n, 0] - working_table[nn, 0]) < ERB_rad and \
               abs(candidates_0[n, 1] - working_table[nn, 1]) < dB_diff:
                is_masked = True
                break
        if not is_masked:
            candidates_1.append(candidates_0[n])
            
    candidates_1 = np.array(candidates_1)
    
    # Prepare outputs
    tier_2 = np.column_stack((erb_to_hz(candidates_0[:, 0]), candidates_0[:, 1]))
    if len(candidates_1) > 0:
        tier_1 = np.column_stack((erb_to_hz(candidates_1[:, 0]), candidates_1[:, 1]))
    else:
        tier_1 = np.array([])
        
    return tier_1, tier_2, aerb, b

def main():
    parser = argparse.ArgumentParser(description="A repertoire-agnostic model of spectral fission.")
    parser.add_argument("file", help="Audio file to analyze")
    parser.add_argument("--hz", action="store_true", help="Visualize x-axis in Hz (default: ERB)")
    parser.add_argument("--list", action="store_true", help="Output list of candidates instead of visualizing")
    
    args = parser.parse_args()
    
    try:
        tier_1, tier_2, aerb, b = get_fission_candidates(args.file)
    except Exception as e:
        print(f"Error analyzing file: {e}")
        sys.exit(1)
        
    if args.list:
        print("Tier 1 Candidates (Frequency, Power+120):")
        if len(tier_1) > 0:
            for c in tier_1:
                print(f"{c[0]:.2f} Hz, {c[1]:.2f}")
        else:
            print("None")
        print("\nTier 2 Candidates (Frequency, Power+120):")
        if len(tier_2) > 0:
            for c in tier_2:
                print(f"{c[0]:.2f} Hz, {c[1]:.2f}")
        else:
            print("None")
    else:
        plt.figure(figsize=(12, 6))
        if args.hz:
            x_axis = erb_to_hz(aerb)
            t2_x = tier_2[:, 0]
            t1_x = tier_1[:, 0] if len(tier_1) > 0 else []
            plt.xlabel("Frequency (Hz)")
        else:
            x_axis = aerb
            t2_x = hz_to_erb(tier_2[:, 0])
            t1_x = hz_to_erb(tier_1[:, 0]) if len(tier_1) > 0 else []
            plt.xlabel("Frequency (ERB)")
            
        plt.plot(x_axis, b, '.', markersize=2, alpha=0.5, label="Spectrum")
        
        # Find indices for plotting stars/plus signs at the right height
        if len(tier_2) > 0:
            plt.plot(t2_x, tier_2[:, 1], 'r+', label="Tier 2 Candidates")
        if len(tier_1) > 0:
            plt.plot(t1_x, tier_1[:, 1], 'g*', label="Tier 1 Candidates")
            
        plt.ylabel("Power (dB+120)")
        plt.title(f"Spectral Fission Candidates for {args.file}")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

if __name__ == "__main__":
    main()

# fission
A repertoire-agnostic model of <a href="https://youtu.be/zJOCxsJA1LA">spectral fission</a>.

Presented at <a href="https://icmpc2021.sites.sheffield.ac.uk/">ICMPC 2021</a>: <a href="https://github.com/jordan-lenchitz/fission/blob/main/ICMPC%202021%20abstract.pdf">abstract</a>, <a href="https://github.com/jordan-lenchitz/fission/blob/main/ICMPC%202021%20poster.jpg">poster</a>, <a href="https://www.youtube.com/watch?v=8TqmxaW4nTQ">2-minute explainer video</a>.

## Implementation Options

This project provides both the original MATLAB implementation and a new, more accessible Python implementation.

### Python Implementation (Recommended)

The Python implementation is open-source and does not require a MATLAB license. It uses `librosa`, `numpy`, and `scipy`.

#### Installation
```bash
pip install -r requirements.txt
```

#### Usage
```bash
python fission.py my_sound.wav --list  # Output list of candidate frequencies
python fission.py my_sound.wav --hz    # Visualize with Hz x-axis
```

### MATLAB Implementation

Built in MATLAB R2020b; tested up through R2022a.

Both main scripts (`fission_candidates_list.m` and `fission_candidates_vis.m`) take as input the sound file of any sustained sonority (polyphonic or monophonic). 

- `fission_candidates_list.m` outputs a two-tiered list of candidate frequencies for spectral fission.
- `fission_candidates_vis.m` visualizes the spectrum and labels the two tiers of candidate frequencies (green stars for the higher tier, red plus for the second, lower tier).

The core logic has been refactored into `get_fission_candidates.m` for better maintainability. Several previously missing helper functions (`hz2erb.m`, `erb2hz.m`, etc.) have been added to ensure the scripts are fully functional out of the box.

#### Usage
```matlab
% List candidates
fission_candidates_list('my_sound.wav')

% Visualize candidates (hz=1 for Hz, 0 for ERB)
fission_candidates_vis('my_sound.wav', 0)
```

## Methodology

After beginning with an overlap-add time stretch with phase randomization to optimize the frequency resolution of the FFT (`scope.m`, `nonerb_positive_pow.m`), a noise-tolerant fast peak finding algorithm (`peakfinder.m`) is applied to the power spectrum, whose results are then narrowed to reflect any dominance of 6 dB or greater by another peak within a 2 ERB radius. These candidates are finally sorted into the two tiers based on whether or not they fall within 2 ERB of any others to minimize the potential role of masking, with those sufficiently separated in the frequency domain from all other candidates placed into the higher tier. 

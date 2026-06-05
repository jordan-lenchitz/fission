# Fission

A repertoire-agnostic model of [spectral fission](https://youtu.be/zJOCxsJA1LA).

Presented at [ICMPC 2021](https://icmpc2021.sites.sheffield.ac.uk/): [abstract](ICMPC%202021%20abstract.pdf), [poster](ICMPC%202021%20poster.jpg), [explainer video](https://www.youtube.com/watch?v=8TqmxaW4nTQ).

## Implementation (2026 Edition)

The modern implementation uses a **Rust** core for high-performance DSP and a **Python** interface for research workflows.

### Key Features
- **Quantum Sound Microscope:** High-precision phase-aware time stretching to "freeze" sonorities.
- **Adaptive Peakfinding:** Statistical SFM-based noise floor estimation.
- **Two-Tier Classification:** Dynamic level-dependent masking and harmonic coherence checking.

### Usage (Python)

```bash
# Setup
python -m venv .venv && source .venv/bin/activate
pip install maturin numpy librosa matplotlib
cd rust && maturin develop && cd ..

# Analyze
python python/fission.py my_sound.wav
```

## Legacy Implementation
The original MATLAB implementation is available in the root directory for reference.

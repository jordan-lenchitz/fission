use pyo3::prelude::*;
use numpy::{PyArray1, PyReadonlyArray1, ToPyArray};
use ndarray::prelude::*;
use realfft::RealFftPlanner;
use num_complex::Complex;
use std::f32::consts::PI;

#[pymodule]
fn fission_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(quantum_scope, m)?)?;
    m.add_function(wrap_pyfunction!(get_fission_tiers, m)?)?;
    m.add_class::<Candidate>()?;
    Ok(())
}

/// Convert frequency in Hz to ERB scale.
fn hz_to_erb(hz: f32) -> f32 {
    21.4 * (1.0 + 0.00437 * hz).log10()
}

/// Convert ERB scale back to Hz.
fn erb_to_hz(erb: f32) -> f32 {
    (10.0f32.powf(erb / 21.4) - 1.0) / 0.00437
}

/// Candidate for spectral fission.
#[pyclass]
#[derive(Clone, Debug)]
struct Candidate {
    #[pyo3(get)]
    frequency: f32,
    #[pyo3(get)]
    power: f32,
    #[pyo3(get)]
    erb: f32,
    #[pyo3(get)]
    tier: i32,
}

#[pyfunction]
fn get_fission_tiers<'py>(
    py: Python<'py>,
    y: PyReadonlyArray1<f32>,
    sr: f32,
    n_fft: usize,
) -> PyResult<(Vec<Candidate>, Vec<Candidate>)> {
    // 1. Get raw candidates (using the adaptive logic)
    let raw_candidates = get_fission_candidates_internal(y.as_array(), sr, n_fft);
    
    if raw_candidates.is_empty() {
        return Ok((vec![], vec![]));
    }

    // 2. Sort by power descending
    let mut sorted_cands = raw_candidates.clone();
    sorted_cands.sort_by(|a, b| b.power.partial_cmp(&a.power).unwrap());

    // 3. Dynamic Masking & Harmonic Coherence Check
    let mut tier_1 = Vec::new();
    let mut tier_2 = Vec::new();

    for i in 0..sorted_cands.len() {
        let cand = &sorted_cands[i];
        
        // Virtual Pitch Check (Harmonic Coherence)
        // If it's a clear harmonic of a lower candidate, it's Tier 2
        let mut is_harmonic = false;
        for j in 0..i {
            let other = &sorted_cands[j];
            if other.frequency < cand.frequency {
                let ratio = cand.frequency / other.frequency;
                let nearest_harm = ratio.round();
                if nearest_harm > 1.0 && (ratio - nearest_harm).abs() < 0.03 {
                    is_harmonic = true;
                    break;
                }
            }
        }

        // Masking Check (Dynamic)
        let mut is_masked = false;
        for j in 0..i {
            let other = &sorted_cands[j];
            let erb_dist = (cand.erb - other.erb).abs();
            
            // Level-dependent masking
            let masking_radius = 2.0 + (other.power - cand.power).max(0.0) / 10.0;
            
            if erb_dist < masking_radius && other.power > cand.power + 6.0 {
                is_masked = true;
                break;
            }
        }

        let mut final_cand = cand.clone();
        if is_masked {
            continue; 
        }

        if is_harmonic {
            final_cand.tier = 2;
            tier_2.push(final_cand);
        } else {
            final_cand.tier = 1;
            tier_1.push(final_cand);
        }
    }

    Ok((tier_1, tier_2))
}

// Internal helper for peak finding to avoid duplicating code
fn get_fission_candidates_internal(
    y: ArrayView1<f32>,
    sr: f32,
    n_fft: usize,
) -> Vec<Candidate> {
    let mut planner = RealFftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(n_fft);

    let mut buf = Array1::zeros(n_fft);
    let copy_len = y.len().min(n_fft);
    buf.slice_mut(s![..copy_len]).assign(&y.slice(s![..copy_len]));
    
    let mut window = Array1::zeros(n_fft);
    for i in 0..n_fft {
        window[i] = 0.5 * (1.0 - (2.0 * PI * i as f32 / (n_fft as f32 - 1.0)).cos());
    }
    buf *= &window;

    let mut spectrum = fft.make_output_vec();
    fft.process(buf.as_slice_mut().unwrap(), &mut spectrum).unwrap();

    let mut power_db = Array1::zeros(spectrum.len());
    for (i, bin) in spectrum.iter().enumerate() {
        let mag = bin.norm() / (n_fft as f32);
        let db = 20.0 * (mag + 1e-12).log10();
        power_db[i] = db.max(-120.0) + 120.0;
    }

    let mut candidates = Vec::new();
    let window_size = 100;
    for i in 1..power_db.len() - 1 {
        if power_db[i] > power_db[i - 1] && power_db[i] > power_db[i + 1] {
            let freq = i as f32 * sr / n_fft as f32;
            if freq < 20.0 || freq > 4000.0 {
                continue;
            }
            let start = i.saturating_sub(window_size);
            let end = (i + window_size).min(power_db.len());
            let mut local_values: Vec<f32> = power_db.slice(s![start..end]).to_vec();
            local_values.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median = local_values[local_values.len() / 2];
            let threshold = median + 12.0;

            if power_db[i] > threshold {
                candidates.push(Candidate {
                    frequency: freq,
                    power: power_db[i],
                    erb: hz_to_erb(freq),
                    tier: 0,
                });
            }
        }
    }
    candidates
}

#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// The 'Quantum Sound Microscope': A phase-locked vocoder to 'freeze' audio.
/// Takes input audio (y) and a stretch amount.
#[pyfunction]
fn quantum_scope<'py>(
    py: Python<'py>,
    y: PyReadonlyArray1<f32>,
    amount: f32,
    window_size: usize,
) -> PyResult<Bound<'py, PyArray1<f32>>> {
    let y = y.as_array();
    let nsamples = y.len();
    
    if amount == 1.0 {
        return Ok(y.to_pyarray(py));
    }

    let ws = window_size;
    let half_ws = ws / 2;
    let displace_pos = (ws as f32 * 0.5) / amount;
    
    // Hann window
    let mut window = Array1::zeros(ws);
    for i in 0..ws {
        window[i] = (1.0 - ((i as f32 / (ws as f32 - 1.0)) * 2.0 - 1.0).powi(2)).powf(1.25);
    }

    let mut planner = RealFftPlanner::<f32>::new();
    let fft = planner.plan_fft_forward(ws);
    let ifft = planner.plan_fft_inverse(ws);

    let mut output_audio = Vec::new();
    let mut old_windowed_buf = Array1::zeros(ws);
    let mut start_pos = 0.0;

    // Random number generator for phase randomization (mimicking scope.m)
    // In a 'Quantum' version, we might use phase-locking, but let's start
    // by accurately porting the 'scope.m' logic to Rust first.
    let mut rng_seed: u32 = 42;

    while (start_pos as usize) + ws <= nsamples {
        let istart_pos = start_pos as usize;
        let mut buf = y.slice(s![istart_pos..istart_pos + ws]).to_owned();
        buf *= &window;

        // FFT
        let mut spectrum = fft.make_output_vec();
        fft.process(buf.as_slice_mut().unwrap(), &mut spectrum).unwrap();

        // Randomize phase
        let spec_len = spectrum.len();
        for (i, bin) in spectrum.iter_mut().enumerate() {
            let mag = bin.norm();
            if i == 0 || i == spec_len - 1 {
                // DC and Nyquist must be purely real for RealFFT
                *bin = Complex::new(mag, 0.0);
            } else {
                // Simple LCG for deterministic 'random' phases for now
                rng_seed = rng_seed.wrapping_mul(1103515245).wrapping_add(12345);
                let ph = (rng_seed as f32 / u32::MAX as f32) * 2.0 * PI;
                *bin = Complex::from_polar(mag, ph);
            }
        }

        // IFFT
        let mut buf_out = ifft.make_output_vec();
        ifft.process(&mut spectrum, &mut buf_out).unwrap();
        
        // Normalize after IFFT (realfft doesn't normalize)
        let mut buf_out_arr = Array1::from_vec(buf_out);
        buf_out_arr /= ws as f32;
        buf_out_arr *= &window;

        // Overlap-add
        if output_audio.is_empty() {
            output_audio.extend_from_slice(buf_out_arr.slice(s![..half_ws]).as_slice().unwrap());
        } else {
            let mut transition = buf_out_arr.slice(s![..half_ws]).to_owned();
            transition += &old_windowed_buf.slice(s![half_ws..]);
            output_audio.extend_from_slice(transition.as_slice().unwrap());
        }
        
        old_windowed_buf = buf_out_arr;
        start_pos += displace_pos;
    }

    // Final tail
    output_audio.extend_from_slice(old_windowed_buf.slice(s![half_ws..]).as_slice().unwrap());

    // Normalize
    let mut result = Array1::from_vec(output_audio);
    let max_val = result.iter().fold(0.0f32, |a, &b| a.max(b.abs()));
    if max_val > 0.0 {
        result /= max_val;
    }

    Ok(result.to_pyarray(py))
}

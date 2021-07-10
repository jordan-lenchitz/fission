function [a,b] = scope(filename, amount)
%open wav file and get 7324 samples per window
[samples, samplerate] = audioread(filename);
window_seconds = 7324 / samplerate;
%determine mono/stereo and number of samples
[nsamples, nchannels] = size(samples);
%ensure windowsize is even # of samples & larger than 16
ws = floor(window_seconds*samplerate);
if ws<16
    ws = 16;
end
ws = optimize_windowsize(ws);
ws = floor(ws/2)*2;
half_ws = floor(ws/2);
%fade out very end of samples 
end_size = floor(samplerate*.05);
if end_size<16
    end_size = 16;
end
for n = 1:nchannels
    samples(nsamples-end_size:nsamples, n) = samples(nsamples-end_size:nsamples, n) .* linspace(1, 0, end_size+1)';
end
%compute the displacement inside the input file
start_pos=1;
displace_pos=(ws*0.5)/amount;
%create window
window = (1-linspace(-1, 1, ws).^2).^1.25';
old_windowed_buf = zeros(ws, 2);
%do the stretch
stretched = [];
while true
    %get the windowed buffer, zero-padding if necessary
    istart_pos=floor(start_pos);
    buf = samples(istart_pos:istart_pos+ws-1, :);
    if length(buf) < ws
        buf = [buf
zeros(ws-length(buf), 2)];
    end
    buf = buf.*window;
    %get the amplitudes of frequencies, discarding phases
    freqs = abs(fft(buf));
    %randomize the phases by multiplication
    %with a random complex number with modulus=1
    ph = 2*pi*rand(length(freqs), 2)*1i;
    freqs = freqs.*exp(ph);
    %do the inverse fft and windows again
    buf = ifft(freqs, 'symmetric');
    buf = buf.*window;
    %overlap-add the output
    output = buf(1:half_ws, :) + old_windowed_buf(half_ws:ws-1, :);
    old_windowed_buf = buf;
    stretched = [stretched; output];
    %advance through the file until no further to advance
    start_pos = start_pos + displace_pos;
    if floor(start_pos)+ws-1 > nsamples
        break
    end
end
%normalize samples to prevent clipping
a = stretched./max(abs(stretched))*(1-eps);
b = samplerate;
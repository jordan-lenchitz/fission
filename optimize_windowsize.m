function ws = optimize_windowsize(ws)
% OPTIMIZE_WINDOWSIZE Optimize window size for FFT (e.g. power of 2)
ws = 2^nextpow2(ws);
end

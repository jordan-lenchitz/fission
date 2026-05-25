function n = nsamples(file)
% NSAMPLES Returns number of samples in audio file
info = audioinfo(file);
n = info.TotalSamples;
end

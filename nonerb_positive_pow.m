function [a, b] = nonerb_positive_pow(file, k)
% produces column vectors of frequencies 20-4000 hz and their corresponding
% powers ranging from 0 to 120 rather than -120 to 0
% number of points to analyze
N = 2^k;  
% determine amount of stretch needed if necessary
ns = nsamples(file);
amount = 1;
while ns<N
    ns = ns*2;
    amount = amount*2;
end
% apply sound microscope if necessary
if amount~=1
    amount = amount/2;
    [y, Fs] = scope(file, amount);
else
    [y, Fs] = audioread(file);
end           
% average stereo to mono if necessary
[~, nchannels] = size(y);
if nchannels == 2
    y = (y(:, 1) + y(:, 2))/2;
end
% compute fft of sound data
c = fft(y, N)/N;        
% frequencies corresponding to powers
f = (1:N/2-1)*Fs/N;        
% compute power at each frequency and clip to -120
p = 20*log10(abs(c(2:N/2)));
p(p<-120) = -120;
% truncate to 4000 Hz max and 20 Hz min
f = f(f<4000);
f = f(f>20);
p = p(1:length(f));
% output vectors
a = f';
b = p+120;
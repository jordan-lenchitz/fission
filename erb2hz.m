function hz = erb2hz(erb)
% ERB2HZ Convert frequency from ERB scale to Hz
hz = (10.^(erb/21.4) - 1) / 0.00437;
end

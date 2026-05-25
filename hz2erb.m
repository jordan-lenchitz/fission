function erb = hz2erb(hz)
% HZ2ERB Convert frequency from Hz to ERB scale
% According to Glasberg and Moore (1990)
erb = 21.4 * log10(1 + 0.00437 * hz);
end

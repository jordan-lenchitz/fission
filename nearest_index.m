function ind = nearest_index(vec, val)
% NEAREST_INDEX Find index of element in vec closest to val
[~, ind] = min(abs(vec - val));
end

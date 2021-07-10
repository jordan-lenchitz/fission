function fission_candidates_list(file)
% produces two-tiered list (tier 1 contained in tier 2) of frequencies and
% power+120 of candidates for spectral fission sorted by power descending
% % file is the audio file of the signal in question
% dB dominance differential for modeling masking
dB_diff = 6;
% dominance within +- ERB's of peak for modeling masking
ERB_rad = 2;
% floor in Hz above which to generate candidate partials-as-pitches
floor = 1000;
% calculates frequencies and powers as positive (0 to 120 rather than -120 to 0)
% using sound microscope if needed to get 2^19-point precision
[a, b] = nonerb_positive_pow(file, 19);
% new + improved peak finder that works better on noisy data
[loc] = peakfinder(b);
ap = a(loc); bp = b(loc);
% convert from Hz to ERBs
aerb = hz2erb(a);
aperb = hz2erb(ap);
% sort peaks by power and reduce to those above floor
sorted_peaks = sortrows([aperb bp], 2, 'descend');
pre_candidates = sorted_peaks(sorted_peaks(:,1)>hz2erb(floor),:);
% create 'parallel' working matrix of indices
pre_cand_inds = [];
for jj = 1:length(pre_candidates)
    pre_cand_inds(jj) = nearest_index(aerb, pre_candidates(jj, 1));
end
pre_cand_inds = pre_cand_inds';
% go through candidates and check for distance to dB diff above to the left
% [ie frequencies below the peak]
dB_results_pt1 = [];
for jj = 1:length(pre_candidates)
    working_ind = pre_cand_inds(jj);
    inc_ind = working_ind-1;
    if b(inc_ind) > pre_candidates(jj, 2)
        continue
    else
        while b(inc_ind) < pre_candidates(jj, 2) + dB_diff && inc_ind>1
            inc_ind = inc_ind-1;
        end
        diff_left = pre_candidates(jj)-aerb(inc_ind);
        dB_results_pt1 = [dB_results_pt1; [pre_candidates(jj,:) diff_left]]; 
    end
end
% eliminate candidates dominated within ERB radius (as likely masked)
dB_results_pt2 = dB_results_pt1(dB_results_pt1(:, 3)>ERB_rad, :);
% create smaller 'parallel' working matrix of indices
mid_cand_inds = [];
for jj = 1:length(dB_results_pt2)
    mid_cand_inds(jj) = nearest_index(aerb, dB_results_pt2(jj, 1));
end
mid_cand_inds = mid_cand_inds';
% for remaining peaks, check for distance to dB diff above to the right
% [ie frequencies above the peak]
raw_cands = [];
for jj = 1:length(dB_results_pt2)
    working_ind = mid_cand_inds(jj);
    inc_ind = working_ind+1;
    if b(inc_ind) > dB_results_pt2(jj, 2)
        continue
    else
        while b(inc_ind) < dB_results_pt2(jj, 2) + dB_diff && inc_ind < length(b)
            inc_ind = inc_ind+1;
        end
        diff_right = aerb(inc_ind)-dB_results_pt2(jj);
        raw_cands(jj, :) = [dB_results_pt2(jj,:) diff_right]; 
    end
end
% eliminate candidates dominated within ERB radius (as likely masked)
less_raw_cands = raw_cands(raw_cands(:, 4)>ERB_rad, :);
% reduce to single candidate per peak cluster
sorted_lrcs = sortrows(less_raw_cands, 1);
cluster_free = [];
working_cluster = [];
x = 1;
while x < length(sorted_lrcs) + 1
    working_cluster = [working_cluster; sorted_lrcs(x, :)];
    if sorted_lrcs(x+1)/sorted_lrcs(x) <= 1.01
        x = x+1;
        continue
    else
        working_cluster = sortrows(working_cluster, 2, 'descend');
        cluster_free = [cluster_free; working_cluster(1,:)];
        x = x+1;
        working_cluster = [];
        continue
    end
end
% further narrow candidates by considering potential masking among them
candidates_0 = sortrows(cluster_free, 2, 'descend');
candidates_1 = [];
for n = 1:size(candidates_0,1)
    nnn = 0;
    if n == 1
        working_table = candidates_0(2:end, :);
    else
        if n == size(candidates_0,1)
            working_table = candidates_0(1:end-1, :);
        else
            working_table = candidates_0([1:n-1, n+1:end], :);
        end
    end
    for nn = 1:size(working_table, 1)
        if abs(candidates_0(n, 1)-working_table(nn, 1)) < ERB_rad && abs(candidates_0(n, 2)-working_table(nn, 2))<dB_diff
            break
        else
            nnn = nnn+1;
        end
    end
    if nnn == size(working_table, 1)
        candidates_1 = [candidates_1; candidates_0(n, :)];
    end
end
% list generation
test_inds = [];
for q = 1:size(candidates_0, 1)
    test_inds(q) = nearest_index(aerb, candidates_0(q));
end
ttest_inds = [];
for q = 1:size(candidates_1, 1)
    ttest_inds(q) = nearest_index(aerb, candidates_1(q));
end
if ~isempty(candidates_1)
    tier_1 = [erb2hz(candidates_1(:,1)) b(ttest_inds)]
    tier_2 = [erb2hz(candidates_0(:,1)) b(test_inds)]
else
    tier_2 = [erb2hz(candidates_0(:,1)) b(test_inds)]
end
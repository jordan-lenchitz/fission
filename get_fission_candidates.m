function [tier_1, tier_2, aerb, b] = get_fission_candidates(file)
% Core logic for finding spectral fission candidates
% Shared by fission_candidates_list and fission_candidates_vis

% Constants
dB_diff = 6;
ERB_rad = 2;
floor_hz = 1000;

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
pre_candidates = sorted_peaks(sorted_peaks(:,1) > hz2erb(floor_hz), :);

if isempty(pre_candidates)
    tier_1 = [];
    tier_2 = [];
    return;
end

% create 'parallel' working matrix of indices
pre_cand_inds = zeros(length(pre_candidates), 1);
for jj = 1:length(pre_candidates)
    pre_cand_inds(jj) = nearest_index(aerb, pre_candidates(jj, 1));
end

% go through candidates and check for distance to dB diff above to the left
% [ie frequencies below the peak]
dB_results_pt1 = [];
for jj = 1:length(pre_candidates)
    working_ind = pre_cand_inds(jj);
    inc_ind = working_ind-1;
    if inc_ind < 1 || b(inc_ind) > pre_candidates(jj, 2)
        continue
    else
        while inc_ind > 1 && b(inc_ind) < pre_candidates(jj, 2) + dB_diff
            inc_ind = inc_ind-1;
        end
        diff_left = pre_candidates(jj, 1) - aerb(inc_ind);
        dB_results_pt1 = [dB_results_pt1; [pre_candidates(jj,:) diff_left]]; 
    end
end

if isempty(dB_results_pt1)
    tier_1 = [];
    tier_2 = [];
    return;
end

% eliminate candidates dominated within ERB radius (as likely masked)
dB_results_pt2 = dB_results_pt1(dB_results_pt1(:, 3) > ERB_rad, :);

if isempty(dB_results_pt2)
    tier_1 = [];
    tier_2 = [];
    return;
end

% create smaller 'parallel' working matrix of indices
mid_cand_inds = zeros(length(dB_results_pt2), 1);
for jj = 1:length(dB_results_pt2)
    mid_cand_inds(jj) = nearest_index(aerb, dB_results_pt2(jj, 1));
end

% for remaining peaks, check for distance to dB diff above to the right
% [ie frequencies above the peak]
raw_cands = [];
for jj = 1:length(dB_results_pt2)
    working_ind = mid_cand_inds(jj);
    inc_ind = working_ind+1;
    if inc_ind > length(b) || b(inc_ind) > dB_results_pt2(jj, 2)
        continue
    else
        while inc_ind < length(b) && b(inc_ind) < dB_results_pt2(jj, 2) + dB_diff
            inc_ind = inc_ind+1;
        end
        diff_right = aerb(inc_ind) - dB_results_pt2(jj, 1);
        raw_cands = [raw_cands; [dB_results_pt2(jj, 1:2) diff_right]]; 
    end
end

if isempty(raw_cands)
    tier_1 = [];
    tier_2 = [];
    return;
end

% eliminate candidates dominated within ERB radius (as likely masked)
less_raw_cands = raw_cands(raw_cands(:, 3) > ERB_rad, :);

if isempty(less_raw_cands)
    tier_1 = [];
    tier_2 = [];
    return;
end

% reduce to single candidate per peak cluster
sorted_lrcs = sortrows(less_raw_cands, 1);
cluster_free = [];
working_cluster = [];
x = 1;
while x <= size(sorted_lrcs, 1)
    working_cluster = [working_cluster; sorted_lrcs(x, :)];
    if x < size(sorted_lrcs, 1) && (sorted_lrcs(x+1, 1) / sorted_lrcs(x, 1) <= 1.01)
        x = x + 1;
        continue
    else
        working_cluster = sortrows(working_cluster, 2, 'descend');
        cluster_free = [cluster_free; working_cluster(1, 1:2)];
        x = x + 1;
        working_cluster = [];
    end
end

% further narrow candidates by considering potential masking among them
candidates_0 = sortrows(cluster_free, 2, 'descend');
candidates_1 = [];
for n = 1:size(candidates_0, 1)
    nnn = 0;
    working_table = candidates_0([1:n-1, n+1:end], :);
    
    is_masked = false;
    for nn = 1:size(working_table, 1)
        if abs(candidates_0(n, 1) - working_table(nn, 1)) < ERB_rad && ...
           abs(candidates_0(n, 2) - working_table(nn, 2)) < dB_diff
            is_masked = true;
            break;
        end
    end
    
    if ~is_masked
        candidates_1 = [candidates_1; candidates_0(n, :)];
    end
end

% Prepare outputs
test_inds = zeros(size(candidates_0, 1), 1);
for q = 1:size(candidates_0, 1)
    test_inds(q) = nearest_index(aerb, candidates_0(q, 1));
end

tier_2 = [erb2hz(candidates_0(:, 1)), b(test_inds)];

if ~isempty(candidates_1)
    ttest_inds = zeros(size(candidates_1, 1), 1);
    for q = 1:size(candidates_1, 1)
        ttest_inds(q) = nearest_index(aerb, candidates_1(q, 1));
    end
    tier_1 = [erb2hz(candidates_1(:, 1)), b(ttest_inds)];
else
    tier_1 = [];
end
end

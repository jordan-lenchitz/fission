function fission_candidates_vis(file, hz)
% % file is the audio file of the signal in question, any format supported by
% % MATLAB should work; hz is a logical, 1 to visualize x-axis in hz, 0 in erb
clf

[tier_1, tier_2, aerb, b] = get_fission_candidates(file);

if isempty(tier_2)
    disp('No candidates found.');
    return;
end

% Prepare indices for plotting
test_inds = zeros(size(tier_2, 1), 1);
for q = 1:size(tier_2, 1)
    test_inds(q) = nearest_index(aerb, hz2erb(tier_2(q, 1)));
end

if ~isempty(tier_1)
    ttest_inds = zeros(size(tier_1, 1), 1);
    for q = 1:size(tier_1, 1)
        ttest_inds(q) = nearest_index(aerb, hz2erb(tier_1(q, 1)));
    end
end

if hz
    if ~isempty(tier_1)
        plot(erb2hz(aerb), b, '.', tier_2(:,1), b(test_inds), '+r', tier_1(:,1), b(ttest_inds), '*g')
    else
        plot(erb2hz(aerb), b, '.', tier_2(:,1), b(test_inds), '+r')
    end
    xlabel('Frequency (Hz)')
else
    if ~isempty(tier_1)
        plot(aerb, b, '.', hz2erb(tier_2(:,1)), b(test_inds), '+r', hz2erb(tier_1(:,1)), b(ttest_inds), '*g')
    else
        plot(aerb, b, '.', hz2erb(tier_2(:,1)), b(test_inds), '+r')
    end
    xlabel('Frequency (ERB)')
end
ylabel('Power (db+120)')
legend('Spectrum', 'Tier 2 Candidates', 'Tier 1 Candidates')
grid on
end

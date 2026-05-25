function fission_candidates_list(file)
% produces two-tiered list (tier 1 contained in tier 2) of frequencies and
% power+120 of candidates for spectral fission sorted by power descending
% % file is the audio file of the signal in question

[tier_1, tier_2] = get_fission_candidates(file);

if ~isempty(tier_1)
    tier_1 = tier_1
    tier_2 = tier_2
else
    tier_2 = tier_2
end
end

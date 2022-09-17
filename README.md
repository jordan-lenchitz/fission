# fission
A repertoire-agnostic model of <a href="https://youtu.be/zJOCxsJA1LA">spectral fission</a>.

Presented at <a href="https://icmpc2021.sites.sheffield.ac.uk/">ICMPC 2021</a>: <a href="https://github.com/jordan-lenchitz/fission/blob/main/ICMPC%202021%20abstract.pdf">abstract</a>, <a href="https://github.com/jordan-lenchitz/fission/blob/main/ICMPC%202021%20poster.jpg">poster</a>, <a href="https://www.youtube.com/watch?v=8TqmxaW4nTQ">2-minute explainer video</a>.

Both main scripts (fission_candidates_list.m and fission_candidates_vis.m) take as input the sound file of any sustained sonority (polyphonic or monophonic). fission_candidates_list.m outputs a two-tiered list of candidate frequencies for spectral fission; fission_candidates_vis visualizes the spectrum and labels the two tiers of candidate peaks with green stars for the higher tier and red plus for the second, lower tier. 

After beginning with an overlap-add time stretch with phase randomization to optimize the frequency resolution of the FFT, a noise-tolerant fast peak finding algorithm is applied to the power spectrum, whose results are then narrowed to reflect any dominance of 6 dB or greater by another peak within a 2 ERB radius. These candidates are finally sorted into the two tiers based on whether or not they fall within 2 ERB of any others to minimize the potential role of masking, with those sufficiently separated in the frequency domain from all other candidates placed into the higher tier. 

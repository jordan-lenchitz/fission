# fission
A repertoire-agnostic model of <a href="https://youtu.be/zJOCxsJA1LA">spectral fission</a>.

Presented at <a href="https://icmpc2021.sites.sheffield.ac.uk/">ICMPC 2021</a>: <a href="https://github.com/jordan-lenchitz/fission/blob/main/ICMPC%202021%20abstract.pdf">abstract</a>, <a href="https://github.com/jordan-lenchitz/fission/blob/main/ICMPC%202021%20poster.jpg">poster</a>, <a href="https://www.youtube.com/watch?v=8TqmxaW4nTQ">2-minute explainer video</a>.

Built in MATLAB R2020b; tested up through R2022a.

Both main scripts (fission_candidates_list.m and fission_candidates_vis.m) take as input the sound file of any sustained sonority (polyphonic or monophonic). fission_candidates_list.m outputs a two-tiered list of candidate frequencies for spectral fission; fission_candidates_vis visualizes the spectrum and labels the two tiers of candidate frequencies for spectral fission with green stars for the higher tier and red plus for the second, lower tier. 

After beginning with an overlap-add time stretch with phase randomization to optimize the frequency resolution of the FFT, a noise-tolerant fast peak finding algorithm is applied to the power spectrum, whose results are then narrowed to reflect any dominance of 6 dB or greater by another peak within a 2 ERB radius. These candidates are finally sorted into the two tiers based on whether or not they fall within 2 ERB of any others to minimize the potential role of masking, with those sufficiently separated in the frequency domain from all other candidates placed into the higher tier. 

Syntax: fission_candidates_list(file)

Syntax: fission_candidates_vis(file, hz)

file is the audio of a sustained sonority for which candidate pitches of spectral fission will be determined, e. g. 'my_sound.wav'; any format supported by MATLAB should work

hz is a logical, 1 to visualize the x-axis in Hz, 0 in ERB

Example: output a two-tiered list of candidate frequencies for spectral fission for 'my_sound.wav'
>fission_candidates_list('my_sound.wav')

Example: visualize the spectrum and labels the two tiers of candidate frequencies for spectral fission with green stars for the higher tier and red plus for the second, lower tier and the x-axis in ERB
>fission_candidates_vis('my_sound.wav',0)

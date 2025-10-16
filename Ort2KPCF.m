% This function calculates the kurtosis of the pair correlation function.
% Inputs:
%   or_t: The time-resolved order parameter (vector of synchronization values over time).
%   N: number of nodes (or channels)
%   window_size: The size of the sliding window for calculating autocorrelation (number of time points per window).
%   overlap: overlap of windows (between 0~1)
%
% Output:
%   kpcf: Kurtosis of the pair correlation function values over different windows.

function [kpcf] = Ort2KPCF(or_t, N, window_size, overlap)

step_size=window_size*(1-overlap);
wn = floor((length(or_t)-window_size) / step_size)+1; % Number of windows
pcf_w = zeros(wn, 1); % Initialize autocorrelation array

% Loop through each window and compute autocorrelation at lag 'tau'
for wi = 1:wn
    pcf_w(wi) = N*var(or_t((wi-1)*step_size + 1 : (wi-1)*step_size + window_size));
end

% Calculate kurtosis of the autocorrelation values
kpcf = kurtosis(pcf_w);
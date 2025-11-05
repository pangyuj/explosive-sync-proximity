% PhaseLagIndex - Calculate Phase Lag Index (PLI) between channels
%
% Input:
%   wdata: Time series data (time x channels)
%
% Output:
%   pli: Phase Lag Index matrix (channels x channels)
%        Range [0,1]: 0=no coupling, 1=strong phase-lagged coupling

function pli = PhaseLagIndex(wdata)

% Extract instantaneous phase from each channel using Hilbert transform
theta = angle(hilbert(wdata));

% Initialize output matrix
chn = size(wdata, 2);
pli = zeros(chn, chn);

% Calculate PLI for all channel pairs
for ch1 = 1:chn
    for ch2 = 1:ch1
        d1 = theta(:, ch1);  % Phase of channel 1
        d2 = theta(:, ch2);  % Phase of channel 2
        
        % PLI: absolute value of mean sign of phase differences
        % Discards phase locking at 0 and π to exclude volume conduction
        pli(ch1, ch2) = abs(mean(sign(sin(d1 - d2))));
        pli(ch2, ch1) = pli(ch1, ch2);  % Symmetric matrix
    end
end
end
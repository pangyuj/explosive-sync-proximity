function [GlobalOrderParameter] = F_cal_GlobalOrderParameter1

X= importdata('ReturnSample.mat');

SZ = size(X);
phi1 = zeros(SZ);
env= phi1; h= phi1;
for c = 1 : SZ(2)
   phi1(:,c) = hilbert(X(:,c));
   env(:,c) = abs(phi1(:,c));%instantaneous envelope of x is the magnitude of X
end
phi2= phi1./env;
z= mean(phi2,2);
r= abs(z); % global order parameter
GlobalOrderParameter= r; % Real number of complex order parameter





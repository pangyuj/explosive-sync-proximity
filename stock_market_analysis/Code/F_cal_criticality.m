function [StdR,ACF] =F_cal_criticality(R,Window)

StdR = zeros(length(R)-Window,1);
ACF = zeros(length(R)-Window,5);
for i = 1: length(R)-Window
    data= R(i:i+Window);
    ac2=autocorr(data,20);
    ACF(i,1:5)= ac2(2:6);
    StdR(i,1)=std(data);
end

clc; clear all; close all;
%% To calculate Response rate and Recovery rate based on Marke Price

MarketPrice= importdata('StockMarketIndex.mat');
Infor= importdata('CountryInforTable.mat');
[Response,Recovery]= f_cal_Response_Recovery(MarketPrice,Infor) ;

%% To calculate Global Order parameter 

Name= importdata('CountryName.mat');
GDP= importdata('GDPperCapita.mat');
EmergingIndex= importdata('EmergingIndex.mat');
Eindx= find(EmergingIndex==1); Dindx= find(isnan(EmergingIndex));
SD= importdata('TradingDate.mat');
Temp = 1:length(Name); % to be included

%% To calculate KACF and KPCF as measure of kurtosis of criticality

%[GlobalOrderParameter]= F_cal_GlobalOrderParameter1
GlobalOrderParameter= importdata('GlobalOrderParameter.mat');

Window= 120; clear AC PC SR Sync
for k= 1 : length(Name)
    R= GlobalOrderParameter{k}; [StdR,ACF]=F_cal_criticality(R,Window); 
    SR{k}=StdR; AC{k}= ACF;
end
save('ACF_Window120.mat','ACF','-mat')
save('PCF_Window120.mat','SR','-mat')

clear KurtosisStdR1 KurtosisStdR Macfw KurtosisStdR2  SD1 ACF1 
for k = 1 : length(Name)
    Idx= find(20060100 <= SD{k}(Window+1:end) & SD{k}(Window+1:end) <= 20061232); 
    KurtosisStdR(k,1)= kurtosis(SR{k}(Idx,1));    
    for k1 = 1: 5
        KurtosisStdR1(k1,k)= kurtosis(AC{k}(Idx,k1));
    end
    ACF1{k}= AC{k}(Idx,:); Sd= SD{k}(Window+1:end); SD1{k}= Sd(Idx);
end

%% To calculate the relation between response(recovery) and kurtosis of criticality 
Kurt=[KurtosisStdR KurtosisStdR1']';
save('FinalKPCF_KACF_lag1_5.mat','Kurt','-mat');

for k = 1: size(Response,2)
    Res= Response(Temp,k); Rec= Recovery(Temp,k);
    for k1= 1: size(Kurt,1)
        [C,p]= corr(Kurt(k1,Temp)',log(1./Res),'type','spearman'); Res_r2(k1,k)=C; Res_p2(k1,k)= p;
        [C,p]= corr(Kurt(k1,Temp)',log(1./Rec),'type','spearman'); Rec_r2(k1,k)=C; Rec_p2(k1,k)= p;
    end
end
Result=[];Result1=[];
for  k1= 1:4
    Result=[Result; Res_r2(k1,:); Res_p2(k1,:);]; Result1=[Result1; Rec_r2(k1,:); Rec_p2(k1,:)];
end
FinalResult= [Result; Result1];

f_Main_Figure


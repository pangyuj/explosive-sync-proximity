function [Response,Recovery]= f_cal_Response_Recovery(Data,Infor)
After=[40:10:200]; %parameter alpha in the Paper



MaxDate=[];MinDate=[];Max=[];Min=[];RecessioSpeed=[];
clear Recovery Response Recovery1 Response1

MinDate= zeros(length(Infor),1); MaxDate= MinDate;
for k = 1 : length(Infor)
    TP= Data{k};
    Date= TP(:,1); Price= TP(:,2); Prc = log(Price);
    if isnan(Infor{k,4})
          Temp= find(20070801<=Date& Date<=20080600);
     else
        Temp= find(20070301<=Date& Date<=Infor{k,4});
     end
    [MaxPrice,Index]= max(Prc(Temp));
    MaxIndex= Temp(Index);
    MaxDate(k) = Date(MaxIndex); % collapse start date
    if isnan(Infor{k,4})
        Temp1= find(20080601<= Date & Date<= 20090631);
    else
     Temp1= find(Infor{k,4}<Date & Date <= Infor{k,5});
    end
    [MinPrice,Index1]= min(Prc(Temp1));
    MinIndex = Temp1(Index1);
    MinDate(k)= Date(MinIndex);
    Diff(k)= MaxIndex-MinIndex;
    Max(k)= MaxPrice;
    Min(k)= MinPrice;
    RecessionSpeed(k)= MaxPrice-MinPrice; 
    for k1= 1:length(After)
        Recovery(k,k1)= Prc(MinIndex+After(k1))- MinPrice;
        Response(k,k1)= Prc(MaxIndex+After(k1))-MaxPrice;
    end
end


for k = 1: length(Infor)
    Recovery(k,:)= Recovery(k,:)./RecessionSpeed(k);
    Response(k,:)= -Response(k,:)./RecessionSpeed(k);
end


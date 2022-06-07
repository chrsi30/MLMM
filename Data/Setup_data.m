function[D] = Setup_data_new()
clear D;
addpath('./Data'); 

%% Hansson 2018
load('./Data/Hansson2018/Hanson_2018_IPGTT');

D.Hansson2018.IPGTT.Glucose = Hanson_2018_IPGTT.Glucose;
D.Hansson2018.IPGTT.Insulin = Hanson_2018_IPGTT.Insulin;

load('./Data/Hansson2018/Hanson_2018_BW');
D.Hansson2018.BW = Hanson_2018_BW;

%% Hansson 2019
load('./Data/Hansson2019/Hansson2019_invitro');
D.Hansson2019.invitro = Hansson2019_invitro;


%% NEW DATA

load('./Data/NEW/NEWDATA_BW.mat')
D.NEW.BW = NEWDATA_BW;
load('./Data/NEW/NEWDATA_IPGTT.mat')
D.NEW.IPGTT = NEWDATA_IPGTT;


end
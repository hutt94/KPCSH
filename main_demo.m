close all; clear all; %clc;
bits = [2 4 8 16 32];
param.db_name ='MIRFLICKR'; %IAPRTC-12 MIRFLICKR NUSWIDE10
param.thresh = 400;
param.mu = 1;
param.alpha = 1;
param.iter = 5; param.sf = 0.05;
param.omega = 1; param.beta = 0.01;
param.theta = 1; param.lambda = 0.001;

for i = 1:5
param.nbits = bits(i);

param.top_K = 1000;
param.pr_ind = [1:50:1000,1000];
param.pn_pos = [1:100:2000,2000];

[XTrain,YTrain,LTrain,XTest,YTest,LTest] = load_data(param.db_name);

fprintf('========%s %d bits start======== \n', 'JSPSH',param.nbits);
evaluate_KPCSH(XTrain,YTrain,LTrain,XTest,YTest,LTest,param);
clearvars -except param bits thresh
end

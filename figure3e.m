
% to get figure3b
% clear;clc;close all

data=importdata("out.csv");
datatemp=data;
aaa=max(data(:,2));
% figure;
% 
load("figure3a.mat");
bbb=max(record1);
figure; plot([1:899]*0.1,(-datatemp((2:899+1),2)/aaa+record1(1:899)'*(1/bbb))*50,'r')

load('figure3b.mat')
bbb=max(record1(1:899));
hold on;plot([1:899]*0.1,(-datatemp((2:899+1),2)/aaa+record1(1:899)'*(1/bbb))*50,'k')

legend('Difference with I-SGFD','Difference with HEI')

% axis([0 89 -0.7 1])

xlabel('time(ms)')
ylabel('Amplitude X 50')
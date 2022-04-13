
% to get figure3b
clear;clc;close all

data=importdata("out.csv");
datatemp=data;
aaa=max(data(:,2));

load("figure3a.mat");
bbb=max(record1);
figure; plot([1:1409]*0.1,(-datatemp((1:1409),2)/aaa+record1(1:1409)'*(1/bbb))*30,'r')
load('figure3b.mat')
bbb=max(record1(1:1409));
hold on;plot([1:1409]*0.1,(-datatemp((1:1409),2)/aaa+record1(1:1409)'*(1/bbb))*30,'k')

legend('Difference with I-SGFD','Difference with HEI')

% axis([0 140 -0.7 1])

xlabel('time(ms)')
ylabel('Amplitude X 30')
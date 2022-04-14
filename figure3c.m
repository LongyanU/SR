% to get figure3a
% clear;clc;close all

data=importdata("out.csv");
datatemp=data;
aaa=max(data(:,2));
figure; plot([1:899]*0.1,data(2:899+1,2)/aaa,'r')
% 
load("figure3a.mat");
aaa=max(data(:,2));
bbb=max(record1);
hold on;plot([1:899]*0.1,record1(1:899)*(1/bbb),'k.')

load('figure3b.mat')
aaa=max(data(:,2));
bbb=max(record1);
hold on;plot([1:899]*0.1,record1(1:899)*(1/bbb),'b--')

legend('Analytic record','I-SGFD result','HEI SGFD result')

axis([0 89 -0.7 1])

xlabel('time(ms)')
ylabel('Amplitude')




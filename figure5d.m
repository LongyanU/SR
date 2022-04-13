clear;
clc
close all;


load('figure5aSaltTraScheme.mat')
figure;imagesc([],[1:2500]*2,seis_recordp(1:2500,45:end-45),[-300,150])
colormap gray
xlabel('x/dx')
ylabel('travel time(m/s)')

aaatemp=seis_recordp;


load('figure5b.mat')
figure;imagesc([],[1:2500]*2,seis_recordp(1:2500,45:end-45),[-300,150])

colormap gray
xlabel('x/dx')
ylabel('travel time(m/s)')


figure;imagesc([],[1:2500]*2,aaatemp(1:2500,45:end-45)-seis_recordp(1:2500,45:end-45),[-300,150])
colormap gray
xlabel('x/dx')
ylabel('travel time(m/s)')



clear;
clc
close all;


load('figure5aSaltTraScheme.mat')
figure;plot([1:2500]*2,-seis_recordp(1:2500,450),'r','linewidth',1)
tt1=-seis_recordp(1:end,450);

load('figure5b.mat')
hold on ;plot([1:2500]*2,-seis_recordp(1:2500,450),'b','linewidth',1)
tt2=-seis_recordp(1:end,450);
hold on ;plot([1:2500]*2,tt2(1:2500)-tt1(1:2500),'k','linewidth',1)



grid on
% legend('ISGFD scheme','HEI-SGFD scheme','The difference bewteen the first 2 schemes', 'Second-order implicit scheme')
legend('ISGFD scheme','HEI-SGFD scheme','The difference bewteen these 2 schemes')
ylabel('p (Pa)')
xlabel('travel time(m/s)')

axis([0 2500*2 -4800 7600])

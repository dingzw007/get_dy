%this is to compare the phonon frequency of the two method
clear all;
close all;
set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',22);
set(0,'DefaultLineMarkerSize',8)
data1=load('band.dat');
data2=load('ww.dat');

qs=0:0.01:0.5;
figure(1)
hold on
plot(data1(1:101,2),data1(1:101,5:end));
plot(qs,data2(:,1:6),'o')
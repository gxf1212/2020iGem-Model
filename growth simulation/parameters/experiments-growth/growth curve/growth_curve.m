clear;clc;clf;

value=csvread('growth.csv'); % 7 by 26
time=csvread('time.csv')'; % 1 by 26
[n,t]=size(value);

%% plot
hold on
DrawFigures('title',time,value,'time/h','OD_{600}','B.S',0)

clear;clc;
clf;

%% settings
data=load('cv_data').cv_data;

hold off
% position
position_ylabel=[-0.1,0.5,0];
position_xlabel=[0.5,-0.06,0];
% color
color_background=[1,1,1];
% line
linewidth_plot=2;
linewidth_axis=2;
% font
fontweight_label="normal";
fontweight_axis="bold";
fontname="Arial";
fontsize_label=15;
fontsize_axis=15;

%% plot
% figure('NumberTitle','off','Name',title,'color',color_background)
figure('position',[100,200,1200,300])
% text
xvalues = {"Pf", "Nf_1", "r_1", "m_1", "r_2", "m_2", "Q_{1c}", "Q_{1n}", "Q_{1p}", "Q_{2c}","Q_{2n}", "Q_{2p}", "K_{c1}", "K_{n1}", "K_{p1}", "K_{n2}", "K_{p2}"};
yvalues = {"5","10","15","20","25"};
idx=1;
species=1;
title="Parameters' Coefficient of Variation"+" (" + get_sp(idx)+", sigma="+string(5)+"%)";

heatmap(xvalues,yvalues,data(:,:,idx,species)); % 'Growth Time/days','Parameters',

% x=xlabel(,'fontname',fontname,'fontsize',fontsize_label,'fontweight',fontweight_label);
% y=ylabel(,'fontname',fontname,'fontsize',fontsize_label,'fontweight',fontweight_label); 
h.Title = title;
h.XLabel = 'Parameters';
h.YLabel = 'Growth Time/days';

% set(gca,'fontsize',fontsize_axis);




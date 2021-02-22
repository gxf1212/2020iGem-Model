clear;clc;
clf;

value=csvread('growth.csv'); % 7 by 26
time=csvread('time.csv')'; % 1 by 26
[n,t]=size(value);

%% plot
% each line for a row number
% position
position_ylabel=[-0.1,0.5,0];
position_xlabel=[0.5,-0.06,0];
% color
color_background=[1,1,1];
color_curve=[1,0,0];
% line
linewidth_plot=2;
linewidth_axis=2;
% text
title="growth curve";
label_x='time/h';
label_y='OD_{600}';
Legend=["EE","EB","EP","BB","BE","BP","DH5a"];
% font
fontweight_label="normal";
fontweight_axis="bold";
fontname="Arial";
fontsize_label=15;
fontsize_axis=15;

figure('NumberTitle','off','Name',title,'color',color_background);
hold on
for i=1:n
    plot(time,value(i,:),'linewidth',linewidth_plot); % ,'color',color_curve
end

legend(Legend,'location','northwest');
legend('boxoff');
x=xlabel(label_x,'fontname',fontname,'fontsize',fontsize_label,'fontweight',fontweight_label);
y=ylabel(label_y,'fontname',fontname,'fontsize',fontsize_label,'fontweight',fontweight_label); 

set(gca,'fontsize',fontsize_axis,'linewidth',linewidth_axis,'fontweight',fontweight_axis);
set(y,'Units','Normalized','Position',position_ylabel);
set(x,'Units','Normalized','Position',position_xlabel);
axis square;
set(gca, 'Color', 'none');
set(gcf,'position',[0 0 1000 1000]);
box off;



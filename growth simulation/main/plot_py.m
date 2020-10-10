% using computed data from python
% using visualization in matlab

global r1 r2 m1 m2 Kc1 Kn1 Kp1 Kn2 Kp2 Re1 NPh Nf1 Nf2 Pf Q1c Q1n Q1p Q2c Q2n Q2p
global ga1 ga2 pH0

time=csvread('./data/time.csv');
N1=csvread('./data/N1.csv');
N2=csvread('./data/N2.csv');
Rc=csvread('./data/Rc.csv');
Rn=csvread('./data/Rn.csv');
Rp=csvread('./data/Rp.csv');
G1=csvread('./data/G1.csv');
G2=csvread('./data/G2.csv');

%% 展示结果
clf;
% 种群
figure(1)
plot(time, N1, time, N2, 'linewidth', 1.5);
title('growth curve');
xlabel('time/min');
ylabel('biomass concentration/(g\cdot L^{-1})');
legend('B.S','Nostoc')


% 资源
figure(2)
set(gcf, 'position', [200,200,900,300]);
% plot(time, Rc, time, Rn, time, Rp);
% legend('C','N','P')
subplot(1,3,1)
plot(time, Rc);
legend('C')

subplot(1,3,2)
plot(time, Rn);
legend('N')

subplot(1,3,3)
plot(time, Rp);
legend('P')


% 生长、死亡速率对比
figure(3)
subplot(2,2,1)
plot(time, G1, time, N1*m1);
title('B.S');
legend('growth','decay')

subplot(2,2,2)
title('Nostoc');
plot(time, G2, time, N2*m2);
legend('growth','decay')

subplot(2,2,3)
title('f of B.S');
hold on
% plot(time, G1./N1/r1/toxin(ga1, N1, N2), 'o','MarkerSize', 2);
pl3_3f=plot(time, G1./N1/r1, 'o','MarkerSize', 1.5);
pl3_3r=plot(time, MM(Rc, Kc1, N1), time, MM(Rn, Kn1, N1), time, MM(Rp, Kp1, N1));
legend('f','C','N','P')

subplot(2,2,4)
title('f of Nostoc');
hold on
% plot(time, G2./N2/r2/toxin(ga2, N1, N2), 'o','MarkerSize', 2);
pl3_4f=plot(time, G2./N2/r2, 'o','MarkerSize', 1.5);
pl3_4r=plot(time, MM(Rn, Kn2, N2), time, MM(Rp, Kp2, N2));
legend([pl3_4f, pl3_4r(1), pl3_4r(2)],'f','N','P')


% figure(4)
% plot(time, toxin(ga2, N1, N2))
% legend('B.S','Nostoc')
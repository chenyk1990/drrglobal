%% DEMO for using LRR method to reconstruct weak phases in global seismograms
% initially by Yangkang Chen, 2017
% Modified by XXXX, 20xx
% 
clc;clear;close all;

%% part I (main phase)
D=zeros(2000,20);
h=zeros(20,1);
for i=0:19
    tmp=strcat('syn-mp-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    D(:,i+1)=data(:);
    
    h(i+1)=str2num(distance);
end

t0=0;
nt=str2num(npts);
dt=str2num(delta);
t=t0+[0:nt-1]*dt;

k1=round((0-t0)/dt)+1;
k2=round((1400-t0)/dt)+1;

D=D(k1:k2,:);
t=t(k1:k2);
figure;yc_wigbh(D,h,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_main.eps')

%% part II (with coda)
D0=zeros(2000,20);
h0=zeros(20,1);
for i=0:19
    tmp=strcat('syn-mp-coda-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    D0(:,i+1)=data(:);
    
    h0(i+1)=str2num(distance);
end

% t0=0;
% nt=str2num(npts);
% dt=str2num(delta);
% t=t0+[0:nt-1]*dt;
% 
% k1=round((550-t0)/dt)+1;
% k2=round((900-t0)/dt)+1;

D0=D0(k1:k2,:);
t0=t;
figure;yc_wigbh(D0,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(200,21,'Coda','color','k','Fontsize',20);
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_main_coda.eps')
% figure;yc_wigbh(D0-D,h0,t,1);ylim([-0.9,20]);

figure;yc_wigbh(D0-D,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(200,21,'Coda','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_coda.eps')


%% inter
mask=ones(size(D0));
nt=size(D0,1);
mask(:,[3,7,11,17])=zeros(nt,4);
D1=D0.*mask;
figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(200,21,'Coda','color','k','Fontsize',20);
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_d0.eps')

figure;imagesc(mask);colormap(jet);colorbar;
thr=1;
niter=50;
type='pocs';
[ D3 ] = fk_pocs(D1,D1,mask,thr,niter,type);

figure;imagesc(t,h,mask');colormap(jet);colorbar;
set(gca,'YDir','normal');
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
print(gcf,'-depsc','-r400','syn_mask.eps');

%% IST
[ D4 ] = fk_pocs(D1,D1,mask,thr,niter,'ist');
figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
figure;yc_wigbh(D3,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(200,21,'Coda','color','k','Fontsize',20);
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_pocs.eps')

figure;yc_wigbh(D4,h0,t,1);ylim([-0.9,20]);

%% SSA

flow=0;fhigh=125;dt=0.004;N=2;Niter=20;mode=0;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
D5=fxymssa_recon(D1,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
figure;yc_wigbh(D5,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(200,21,'Coda','color','k','Fontsize',20);
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_rr.eps')

figure('units','normalized','Position',[0.2 0.4 0.5, 0.8]);
subplot(4,1,1);plot(D3(:,3));ylim([-1,1]);
subplot(4,1,2);plot(D5(:,3));ylim([-1,1]);
subplot(4,1,3);plot(D3(:,2));ylim([-1,1]);
subplot(4,1,4);plot(D3(:,4));ylim([-1,1]);

yc_snr(D0,D3)
yc_snr(D0,D5)


figure;yc_wigbh(D0-D3,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
annotation(gcf,'rectangle',...
    [0.679391959798995 0.162393162393162 0.0534388609715242 0.700854700854704],...
    'LineWidth',2,...
    'LineStyle','--');
print(gcf,'-depsc','-r400','syn_pocs_e.eps')

figure;yc_wigbh(D0-D5,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
print(gcf,'-depsc','-r400','syn_rr_e.eps')


figure;yc_wigbh((D0-D3)*10,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
print(gcf,'-depsc','-r400','syn_pocs_e2.eps')

figure;yc_wigbh((D0-D5)*10,h0,t,1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
print(gcf,'-depsc','-r400','syn_rr_e2.eps')



%% plot rank




Nmax=4;f=20;%0.08 Hz
f1=10;%0.04Hz,25s
f2=20;%0.08Hz,12.5s
f3=25;
randn('state',201011);
D5tmp=D5+0.1*randn(size(D5));
[ D6,H1 ] = fxymssa_auto(D5tmp,flow,fhigh,dt,Nmax,verb,2,f1,eps);
[ D66,H2 ] = fxymssa_auto(D5tmp,flow,fhigh,dt,Nmax,verb,2,f2,eps);
[ D666,H3 ] = fxymssa_auto(D5tmp,flow,fhigh,dt,Nmax,verb,2,f3,eps);

h1=exe_dia(H1);
h2=h1(2:end);
r1=h1(1:end-1)./h2;
figure;stem(h1);
figure;stem(r1);

h11=exe_dia(H2);
h22=h11(2:end);
r2=h11(1:end-1)./h22;
figure;stem(h11);
figure;stem(r2);

h111=exe_dia(H3);
h222=h111(2:end);
r3=h1(1:end-1)./h2;
figure;stem(h111);
figure;stem(r3);


figure('units','normalized','Position',[0.2 0.4 0.6, 0.6],'color','w');
subplot(2,2,1);plot(h1(1:8),'-rx','linewidth',2);hold on;
ylabel('Value','Fontsize',16);
title('(a)','Fontsize',16);
plot(1*ones(80,1),1:80,'g--','linewidth',2);ylim([0,80]);xlim([0.5,8.5]);
set(gca,'Linewidth',1.5,'Fontsize',16);
subplot(2,2,3);plot(r1(1:8),'b-o','linewidth',2);
ylabel('Ratio','Fontsize',16);hold on;
title('(c)','Fontsize',16);
plot(1*ones(60,1),linspace(0.9,4.5,60),'g--','linewidth',2);ylim([0.9,4.5]);xlim([0.5,8.5]);
 xlabel('Singular Value NO','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);

subplot(2,2,2);plot(h11(1:8),'-rx','linewidth',2);hold on;
%ylabel('Value','Fontsize',16);
title('(b)','Fontsize',16);
plot(6*ones(30,1),1:30,'g--','linewidth',2);ylim([0,30]);xlim([0.5,8.5]);
set(gca,'Linewidth',1.5,'Fontsize',16);
subplot(2,2,4);plot(r2(1:8),'b-o','linewidth',2);
% ylabel('Ratio','Fontsize',16);
hold on;
title('(d)','Fontsize',16);
plot(6*ones(20,1),linspace(1,1.5,20),'g--','linewidth',2);ylim([1,1.5]);xlim([0.5,8.5]);
xlabel('Singular Value NO','Fontsize',16);
set(gca,'Linewidth',1.5,'Fontsize',16);

print(gcf,'-depsc','-r400','syn_sigmas.eps')

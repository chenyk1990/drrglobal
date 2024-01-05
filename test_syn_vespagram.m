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

D=D(:,:);
t=t(:);
figure;yc_wigbh(D,h,t,1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_main.eps')


%% vespagram
s=-5:.1:10;
h=0:19;
t=[0:nt-1]*dt;
order=4;
% Param.h=h;
% Param.v=s;
% Param.nt=nt;
% Param.dt=dt;
% Param.type=4;
% Param.order=1/4;
% 
% ves=radon_pseudo(din,Param,-1);
% ves=sign(ves).*abs(ves).^order;

[ves] = vespagram(D,h,s,t,order);

figure;imagesc(t(1:1401),s,ves(1:1401,:)');colormap(seis);caxis([-0.001,0.001]);
set(gca, 'YDir','normal')
ylabel('Slowness (s/deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
% text(200,21,'Coda','color','k','Fontsize',20);
text(100,5,'Precursor','color','k','Fontsize',20);
text(1050,5,'Main Phase','color','k','Fontsize',20);
annotation(gcf,'ellipse',...
    [0.361134003350084 0.358974358974359 0.0852646566164155 0.478632478632481],...
    'LineWidth',2,...
    'LineStyle','--');
annotation(gcf,'ellipse',...
    [0.473361809045226 0.346153846153847 0.0852646566164155 0.478632478632481],...
    'LineWidth',2,...
    'LineStyle','--');
annotation(gcf,'ellipse',...
    [0.638353433835846 0.363247863247865 0.0852646566164155 0.478632478632482],...
    'LineWidth',2,...
    'LineStyle','--');
print(gcf,'-depsc','-r400','syn_ves.eps')

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

D0=D0(:,:);
t0=t;
figure;yc_wigbh(D0,h0,t,1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_main_coda.eps')
% figure;yc_wigbh(D0-D,h0,t,1);ylim([-0.9,20]);

[ves0] = vespagram(D0,h,s,t,order);

figure;imagesc(t(:),s,ves0(:,:)');colormap(seis);caxis([-0.5,0.5]);
set(gca, 'YDir','normal')


coda=D0-D;
figure;yc_wigbh(D0-D,h0,t,1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_coda.eps')


%% inter
mask=ones(size(D0));
nt=size(D0,1);
mask(:,[3,7,11,17])=zeros(nt,4);
D1=D0.*mask;
figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_d0.eps')

figure;imagesc(mask);colormap(jet);colorbar;
thr=1;
niter=50;
type='pocs';
[ D3 ] = fk_pocs(D1,D1,mask,thr,niter,type);
% print(gcf,'-depsc','-r400','syn_pocs.eps')
figure;yc_wigbh(D3(1:1401,:)-coda(1:1401,:),h0,t(1:1401),1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
annotation(gcf,'rectangle',...
    [0.226293132328308 0.153846153846154 0.0375259631490787 0.760683760683762],...
    'LineWidth',2,...
    'LineStyle','--');
print(gcf,'-depsc','-r400','syn_pocs_mp.eps')



[ves3] = vespagram(D3-coda,h,s,t,order);
figure;imagesc(t(1:1401),s,ves3(1:1401,:)');colormap(seis);caxis([-0.001,0.001]);
set(gca, 'YDir','normal')
ylabel('Slowness (s/deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
annotation(gcf,'rectangle',...
    [0.151753768844221 0.170940170940171 0.170691792294807 0.452991452991457],...
    'LineWidth',2,...
    'LineStyle','--');
print(gcf,'-depsc','-r400','syn_ves_pocs.eps')


%% IST
[ D4 ] = fk_pocs(D1,D1,mask,thr,niter,'ist');
figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
figure;yc_wigbh(D3,h0,t,1);ylim([-0.9,20]);
figure;yc_wigbh(D4,h0,t,1);ylim([-0.9,20]);

%% SSA

flow=0;fhigh=125;dt=0.004;N=2;Niter=20;mode=0;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
D5=fxymssa_recon(D1,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
figure;yc_wigbh(D5(1:1401,:)-coda(1:1401,:),h0,t(1:1401),1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_rr.eps')
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_rr_mp.eps')

[ves5] = vespagram(D5-coda,h,s,t,order);
figure;imagesc(t(1:1401),s,ves5(1:1401,:)');colormap(seis);caxis([-0.001,0.001]);
set(gca, 'YDir','normal')
ylabel('Slowness (s/deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
print(gcf,'-depsc','-r400','syn_ves_rr.eps')

figure('units','normalized','Position',[0.2 0.4 0.5, 0.8]);
subplot(4,1,1);plot(D3(:,3));ylim([-1,1]);
subplot(4,1,2);plot(D5(:,3));ylim([-1,1]);
subplot(4,1,3);plot(D3(:,2));ylim([-1,1]);
subplot(4,1,4);plot(D3(:,4));ylim([-1,1]);



figure;yc_wigbh(D0-D3,h0,t,1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_pocs_e.eps')

figure;yc_wigbh(D0-D5,h0,t,1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_rr_e.eps')



figure;yc_wigbh(2*(D3(1:1401,:)-coda(1:1401,:)),h0,t(1:1401),1);ylim([-0.9,20]);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_pocs_mp2.eps')

figure;yc_wigbh(2*(D5(1:1401,:)-coda(1:1401,:)),h0,t(1:1401),1);ylim([-0.9,20]);
% print(gcf,'-depsc','-r400','syn_rr.eps')
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 
text(550,21,'Precursor','color','k','Fontsize',20);
text(950,21,'Main Phase','color','k','Fontsize',20);
print(gcf,'-depsc','-r400','syn_rr_mp2.eps')





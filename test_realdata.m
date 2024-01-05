%% DEMO for using DRR method to reconstruct weak phases in global seismograms
% initially by Yangkang Chen, 2017
% Modified by Wei Chen, 2021
% 

%% This script requires the MATdrr package
%https://github.com/aaspip/MATdrr

clc;clear;close all;
addpath(genpath('./'));


%% part 0
Draw=zeros(13732,12);
hraw=zeros(12,1);
for i=0:11
    tmp=strcat('data/raw-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    Draw(:,i+1)=data(:);
    
    hraw(i+1)=str2num(distance);
end

%% part I
D=zeros(13732,12);
h=zeros(12,1);
for i=0:11
    tmp=strcat('data/data-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    D(:,i+1)=data(:);
    
    h(i+1)=str2num(distance);
end

t0=100;
nt=str2num(npts);
dt=str2num(delta);
t=t0+[0:nt-1]*dt;


% plotht(D,h,t);

%% part I.5
D1=zeros(13732,12);
h1=zeros(12,1);
for i=0:11
    tmp=strcat('data/datares-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    D1(:,i+1)=data(:);
    
    h1(i+1)=str2num(distance);
end


% plotht(D1,h1,t);


%% part II
D2=zeros(13732,12);
h2=zeros(12,1);
for i=0:11
    tmp=strcat('data/datarec-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    D2(:,i+1)=data(:);
    
    h2(i+1)=str2num(distance);
end

% t0=100;
% nt=str2num(npts);
% dt=str2num(delta);
% t=t0+[0:nt-1]*dt;

% plotht(D2,h2,t);

figure('units','normalized','Position',[0.2 0.4 0.8, 0.7],'color','w');
yc_wigbh(D2,h2,dt*[0:nt-1],2);


figure('units','normalized','Position',[0.2 0.4 0.8, 0.7],'color','w');
yc_wigbh(D1,h1,0.05*[0:nt-1],2);
ylabel('Distance (deg)','Fontsize',30);
xlabel('Time (s)','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30); 
annotation(gcf,'line',[0.551088777219431 0.574539363484087],...
    [0.136752136752137 0.908119658119658],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.443048576214405 0.485762144053601],...
    [0.136752136752137 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.415410385259632 0.458123953098828],...
    [0.136752136752138 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.311557788944724 0.313232830820771],...
    [0.134615384615385 0.905982905982906],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
text(160,95,'P','color','r','Fontsize',25);
text(390,95,'PP','color','r','Fontsize',25);
text(260,95,'P^{410}P','color','r','Fontsize',25);
text(310,95,'P^{660}P','color','r','Fontsize',25);
print(gcf,'-depsc','-r400','real_d0.eps');

figure('units','normalized','Position',[0.2 0.4 0.8, 0.7],'color','w');
yc_wigbh(D,h,0.05*[0:nt-1],2);
ylabel('Distance (deg)','Fontsize',30);
xlabel('Time (s)','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30);
annotation(gcf,'line',[0.551088777219431 0.574539363484087],...
    [0.136752136752137 0.908119658119658],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.443048576214405 0.485762144053601],...
    [0.136752136752137 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.415410385259632 0.458123953098828],...
    [0.136752136752138 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.311557788944724 0.313232830820771],...
    [0.134615384615385 0.905982905982906],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
text(160,95.15,'P','color','r','Fontsize',25);
text(390,95.15,'PP','color','r','Fontsize',25);
text(260,95.15,'P^{410}P','color','r','Fontsize',25);
text(310,95.15,'P^{660}P','color','r','Fontsize',25);

print(gcf,'-depsc','-r400','real_raw.eps');


figure('units','normalized','Position',[0.2 0.4 0.8, 0.7],'color','w');
yc_wigbh(Draw,hraw,0.05*[0:nt-1],2);
ylabel('Distance (deg)','Fontsize',30);
xlabel('Time (s)','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30); 
print(gcf,'-depsc','-r400','real_original.eps');


%   ylabel('Time (s)','Fontsize',16);
%   xlabel('Trace','Fontsize',16);
%   title('Noise','Fontsize',16);
%   set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold'); 


%% inter
mask=ones(size(D1));
mask(:,[7,9])=zeros(nt,2);
% figure;imagesc(mask);colormap(jet);colorbar;

figure('units','normalized','Position',[0.2 0.4 0.8, 0.7],'color','w');
imagesc(0.05*[0:nt-1],h1,mask');colormap(jet);colorbar;
set(gca,'YDir','normal');

set(gca,'Linewidth',2,'Fontsize',30); 
ylabel('Distance (deg)','Fontsize',30);
xlabel('Time (s)','Fontsize',30);

print(gcf,'-depsc','-r400','real_mask.eps');



thr=0.4;
niter=50;
type='pocs';
[ D3 ] = fk_pocs(D1,D1,mask,thr,niter,type);

%% IST
[ D4 ] = fk_pocs(D1,D1,mask,thr,niter,'ist');
figure('units','normalized','Position',[0.2 0.4 0.8, 0.7],'color','w');
yc_wigbh(D3,h1,0.05*[0:nt-1],2);
ylabel('Distance (deg)','Fontsize',30);
xlabel('Time (s)','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30); 
annotation(gcf,'line',[0.551088777219431 0.574539363484087],...
    [0.136752136752137 0.908119658119658],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.443048576214405 0.485762144053601],...
    [0.136752136752137 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.415410385259632 0.458123953098828],...
    [0.136752136752138 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.311557788944724 0.313232830820771],...
    [0.134615384615385 0.905982905982906],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
text(160,95.0,'P','color','r','Fontsize',25);
text(390,95.0,'PP','color','r','Fontsize',25);
text(260,95.0,'P^{410}P','color','r','Fontsize',25);
text(310,95.0,'P^{660}P','color','r','Fontsize',25);
print(gcf,'-depsc','-r400','real_pocs.eps');

figure;yc_wigbh(D4,h1,dt*[0:nt-1],2);

%% SSA
Niter=20;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
flow=0;fhigh=125;dt=0.004;N=2;K=5;Niter=20;mode=0;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
D5=drr3drecon(D1,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);

figure('units','normalized','Position',[0.2 0.4 0.8, 0.7],'color','w');
yc_wigbh(D5,h1,0.05*[0:nt-1],2);
ylabel('Distance (deg)','Fontsize',30);
xlabel('Time (s)','Fontsize',30);
set(gca,'Linewidth',2,'Fontsize',30); 
annotation(gcf,'line',[0.551088777219431 0.574539363484087],...
    [0.136752136752137 0.908119658119658],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.443048576214405 0.485762144053601],...
    [0.136752136752137 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.415410385259632 0.458123953098828],...
    [0.136752136752138 0.910256410256411],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
annotation(gcf,'line',[0.311557788944724 0.313232830820771],...
    [0.134615384615385 0.905982905982906],'Color',[1 0 0],'LineWidth',3,...
    'LineStyle','--');
text(160,95.0,'P','color','r','Fontsize',25);
text(390,95.0,'PP','color','r','Fontsize',25);
text(260,95.0,'P^{410}P','color','r','Fontsize',25);
text(310,95.0,'P^{660}P','color','r','Fontsize',25);
print(gcf,'-depsc','-r400','real_drr.eps');

figure('units','normalized','Position',[0.2 0.4 0.5, 0.8]);
subplot(4,1,1);plot(D3(:,7));ylim([-1,1]);
subplot(4,1,2);plot(D5(:,7));ylim([-1,1]);
subplot(4,1,3);plot(D3(:,6));ylim([-1,1]);
subplot(4,1,4);plot(D3(:,8));ylim([-1,1]);





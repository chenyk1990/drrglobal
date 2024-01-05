%% DEMO for using LRR method to reconstruct weak phases in global seismograms
% initially by Yangkang Chen, 2017
% Modified by XXXX, 20xx
% 
clc;clear;close all;

%% part I
D=zeros(36903,47);
h=zeros(47,1);
for i=0:46
    tmp=strcat('syngine-data-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    D(:,i+1)=data(:);
    
    h(i+1)=str2num(distance);
end

t0=0;
nt=str2num(npts);
dt=str2num(delta);
t=t0+[0:nt-1]*dt;

k1=round((550-t0)/dt)+1;
k2=round((900-t0)/dt)+1;

D=D(k1:k2,:);
t=t(k1:k2);
figure;yc_wigbh(D,h,t,1);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 

syngine_anno;
print(gcf,'-depsc','-r400','syngine_raw.eps');


%% part I.5
D1=zeros(36903,47);
h1=zeros(47,1);
for i=0:46
    tmp=strcat('syngine-datares-',num2str(i),'.mat');
    load(tmp);
%     figure(1);plot(data);pause(0.1);
    D1(:,i+1)=data(:);
    
    h1(i+1)=str2num(distance);
end

% t0=0;
% nt=str2num(npts);
% dt=str2num(delta);
% t=t0+[0:nt-1]*dt;
% 
% k1=round((550-t0)/dt)+1;
% k2=round((900-t0)/dt)+1;

D1=D1(k1:k2,:);
t1=t;
figure;yc_wigbh(D1,h1,t1,1);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 

syngine_anno;
print(gcf,'-depsc','-r400','syngine_d0.eps');

%% inter
mask=ones(size(D1));
nt=size(D,1);
mask(:,[6,36,43])=zeros(nt,3);
figure;imagesc(t1,h1,mask');colormap(jet);colorbar;
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20,'Ydir','Normal');
print(gcf,'-depsc','-r400','syngine_mask.eps');

thr=0.4;
niter=50;
type='pocs';
[ D3 ] = fk_pocs(D1,D1,mask,thr,niter,type);

%% IST
[ D4 ] = fk_pocs(D1,D1,mask,thr,niter,'ist');
figure;yc_wigbh(D3,h1,dt*[0:nt-1],1);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 

syngine_anno;
% 创建 ellipse
annotation(gcf,'ellipse',...
    [0.135840871021776 0.598290598290599 0.209217755443886 0.196581196581197],...
    'LineWidth',2,...
    'LineStyle','--');

% 创建 rectangle
annotation(gcf,'rectangle',...
    [0.57035175879397 0.566239316239316 0.0678391959798995 0.200854700854702],...
    'LineWidth',2,...
    'LineStyle','--');

print(gcf,'-depsc','-r400','syngine_pocs.eps');

figure;yc_wigbh(D4,h1,dt*[0:nt-1],1);

%% SSA
flow=0;fhigh=125;dt=0.004;N=2;Niter=20;mode=0;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
D5=fxymssa_recon(D1,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
figure;yc_wigbh(D5,h1,dt*[0:nt-1],1);
ylabel('Distance (deg)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
set(gca,'Linewidth',2,'Fontsize',20); 

syngine_anno;
print(gcf,'-depsc','-r400','syngine_rr.eps');

figure('units','normalized','Position',[0.2 0.4 0.5, 0.8]);
subplot(4,1,1);plot(D3(:,6));ylim([-1,1]);
subplot(4,1,2);plot(D5(:,6));ylim([-1,1]);
subplot(4,1,3);plot(D3(:,5));ylim([-1,1]);
subplot(4,1,4);plot(D3(:,7));ylim([-1,1]);





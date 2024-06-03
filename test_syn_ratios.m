%% DEMO for using LRR method to reconstruct weak phases in global seismograms
% initially by Yangkang Chen, 2017
% Modified by Wei Chen, 2020
% 
clc;clear;close all;

addpath(genpath('./'));

%% part I (main phase)
D=zeros(2000,20);
h=zeros(20,1);
for i=0:19
    tmp=strcat('data/syn-mp-',num2str(i),'.mat');
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


%% part II (with coda)
D0=zeros(2000,20);
h0=zeros(20,1);
for i=0:19
    tmp=strcat('data/syn-mp-coda-',num2str(i),'.mat');
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
% figure;yc_wigbh(D0-D,h0,t,1);ylim([-0.9,20]);

figure;yc_wigbh(D0-D,h0,t,1);ylim([-0.9,20]);

coda=D0-D;


%% inter
mask=ones(size(D0));
nt=size(D0,1);
mask(:,[3,7,11,17])=zeros(nt,4);

% sigmas=[0.01,0.02,0.03];


ratios=[0.2:0.2:0.8];


sigma=0.1;

masks=cell(5,1);
masks{1}=[3];
masks{2}=[3,7];
masks{3}=[3,7,11];
masks{4}=[3,7,11,17];
masks{5}=[3,7,11,14,17];
% masks{6}=[3,7,11,14,17,20];

snrs1=zeros(length(masks),1);
snrs2=zeros(length(masks),1);
snrs3=zeros(length(masks),1);

s=-5:.1:10;
h=0:19;
t=[0:nt-1]*dt;
order=4;

for is=1:length(masks);
    randn('state',123456);
    n=sigma*randn(size(D0));
    
    mask=ones(size(D0));
    tmp=masks{is};
    mask(:,tmp)=zeros(nt,length(tmp));
    
    D1=(D0+n).*mask;
    figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
    
    figure;imagesc(mask);colormap(jet);colorbar;
    thr=1;
    niter=50;
    type='pocs';
    [ D3 ] = fk_pocs(D1,D1,mask,thr,niter,type);
    
    figure;imagesc(t,h,mask');colormap(jet);colorbar;
    set(gca,'YDir','normal');
    
    %% IST
    [ D4 ] = fk_pocs(D1,D1,mask,thr,niter,'ist');
    figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
    figure;yc_wigbh(D3,h0,t,1);ylim([-0.9,20]);
    
    
    %% SSA
    
    flow=0;fhigh=125;dt=0.004;N=3;Niter=20;mode=1;verb=1;
    a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
    D5=fxymssa_recon(D1,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
    figure;yc_wigbh(D5,h0,t,1);ylim([-0.9,20]);
    
    
    figure('units','normalized','Position',[0.2 0.4 0.5, 0.8]);
    subplot(4,1,1);plot(D3(:,3));ylim([-1,1]);
    subplot(4,1,2);plot(D5(:,3));ylim([-1,1]);
    subplot(4,1,3);plot(D3(:,2));ylim([-1,1]);
    subplot(4,1,4);plot(D3(:,4));ylim([-1,1]);
    
    snrs1(is)=yc_snr(D0,D1);
    snrs2(is)=yc_snr(D0,D3);
    snrs3(is)=yc_snr(D0,D5);
    
    if is==2
        
        figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
        ylabel('Distance (deg)','Fontsize',20);
        xlabel('Time (s)','Fontsize',20);
        set(gca,'Linewidth',2,'Fontsize',20);
        print(gcf,'-depsc','-r400','syn_obs_ratio2.eps');
        figure;yc_wigbh(D3,h0,t,1);ylim([-0.9,20]);
        ylabel('Distance (deg)','Fontsize',20);
        xlabel('Time (s)','Fontsize',20);
        set(gca,'Linewidth',2,'Fontsize',20);
        print(gcf,'-depsc','-r400','syn_pocs_ratio2.eps');
        figure;yc_wigbh(D5,h0,t,1);ylim([-0.9,20]);
        ylabel('Distance (deg)','Fontsize',20);
        xlabel('Time (s)','Fontsize',20);
        set(gca,'Linewidth',2,'Fontsize',20);
        print(gcf,'-depsc','-r400','syn_rr_ratio2.eps');
        
        [ves1] = yc_vespagram(D1-coda.*mask,h,s,t,order);
        figure;imagesc(t(1:1401),s,ves1(1:1401,:)');colormap(seis);caxis([-0.01,0.01]);
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
        print(gcf,'-depsc','-r400','syn_ves_obs_ratio2.eps')
        
        
        [ves3] = yc_vespagram(D3-coda,h,s,t,order);
        figure;imagesc(t(1:1401),s,ves3(1:1401,:)');colormap(seis);caxis([-0.01,0.01]);
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
        print(gcf,'-depsc','-r400','syn_ves_pocs_ratio2.eps')
        
        
        [ves5] = yc_vespagram(D5-coda,h,s,t,order);
        figure;imagesc(t(1:1401),s,ves5(1:1401,:)');colormap(seis);caxis([-0.01,0.01]);
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
        print(gcf,'-depsc','-r400','syn_ves_rr_ratio2.eps')
    end
    
    
    if is==5
        
        figure;yc_wigbh(D1,h0,t,1);ylim([-0.9,20]);
        ylabel('Distance (deg)','Fontsize',20);
        xlabel('Time (s)','Fontsize',20);
        set(gca,'Linewidth',2,'Fontsize',20);
        print(gcf,'-depsc','-r400','syn_obs_ratio5.eps');
        figure;yc_wigbh(D3,h0,t,1);ylim([-0.9,20]);
        ylabel('Distance (deg)','Fontsize',20);
        xlabel('Time (s)','Fontsize',20);
        set(gca,'Linewidth',2,'Fontsize',20);
        print(gcf,'-depsc','-r400','syn_pocs_ratio5.eps');
        figure;yc_wigbh(D5,h0,t,1);ylim([-0.9,20]);
        ylabel('Distance (deg)','Fontsize',20);
        xlabel('Time (s)','Fontsize',20);
        set(gca,'Linewidth',2,'Fontsize',20);
        print(gcf,'-depsc','-r400','syn_rr_ratio5.eps');
        
        [ves1] = yc_vespagram(D1-coda.*mask,h,s,t,order);
        figure;imagesc(t(1:1401),s,ves1(1:1401,:)');colormap(seis);caxis([-0.01,0.01]);
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
        print(gcf,'-depsc','-r400','syn_ves_obs_ratio5.eps')
        
        
        [ves3] = yc_vespagram(D3-coda,h,s,t,order);
        figure;imagesc(t(1:1401),s,ves3(1:1401,:)');colormap(seis);caxis([-0.01,0.01]);
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
        print(gcf,'-depsc','-r400','syn_ves_pocs_ratio5.eps')
        
        
        [ves5] = yc_vespagram(D5-coda,h,s,t,order);
        figure;imagesc(t(1:1401),s,ves5(1:1401,:)');colormap(seis);caxis([-0.01,0.01]);
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
        print(gcf,'-depsc','-r400','syn_ves_rr_ratio5.eps')
    end
    
    
end


% figure;wigbcyk_wf(dn_2,W,'r');
% ylabel('Time (ms)','Fontsize',20);
% xlabel('Trace','Fontsize',20);
% title('Picked waveforms','Fontsize',20,'color','r');
% set(gca,'Linewidth',2,'Fontsize',20);
%   print(gcf,'-depsc','-r200',['../Fig/syn1-kmeans-',num2str(i),'.eps']);
%
%  pause(0.01);
%  error1=sum(abs(U0(:)-U(:)));
%  errors1=[errors1,error1];
% end
x=[1:length(masks)];
figure;plot(x,snrs1,'ko-','linewidth',2);hold on;
plot(x,flipud(snrs2),'bv-','linewidth',2);
plot(x,snrs3,'r*-','linewidth',2);

% plot(sigmas,400*ones(size(sigmas)),'r','linewidth',2);
% plot(sigmas,errors2,'b','linewidth',2);

% [AX,H1,H2] = plotyy(sigmas,errors1,sigmas,snrs,'plot');
% set(get(AX(1),'Ylabel'),'String','Error (samples)','Fontsize',20,'color','k');
% set(get(AX(2),'Ylabel'),'String','SNR (dB)','Fontsize',20,'color','g');
% set(AX(1),'ycolor','k')
% set(AX(2),'ycolor','g')
%
%
% set(H1,'LineStyle','-','linewidth',2,'color','b');
% set(H2,'LineStyle','-','linewidth',2,'color','g');

legend('Input quality','POCS','Proposed','Location','best');
ylabel('Reconstruction quality (dB)','Fontsize',20);
xlabel('Number of missing seismograms','Fontsize',20);
title('Quality V.S. Sampling ratio','Fontsize',20,'color','k');
set(gca,'Linewidth',2,'Fontsize',20);
%   set(AX(2),'Linewidth',2,'Fontsize',20);
print(gcf,'-depsc','-r200','syn-snr-ratio.eps');





clear all
close all
norm = 1;
path(path,'../Functions')
tic


%% Load shot data
load Observed_data.mat data angles time;
Nl = length(time);

%% Load wavelet
load Wavelet.mat wav dt t0;
t0_i = floor(t0/dt);

%% Load initial model
load Initial_model.mat vp0 vs0 rho0;
%% Angles
theta = angles*pi/180;

%% Make weight
weight = ones(size(data));


%% Inversion starts here
kvp = 100.0;
kvs = 200.0;
krho = 50.0;

vpscale=kvp; 
vsscale=kvs; 
rhoscale=krho;

xi = zeros(3*Nl,1); 

%% Parameters 
mute = ones(size(vp0)); % In case if we need to mute the data
mute(1:2) = 0;

alpha_vp = 0e-6; % For Tikhonov regularization
alpha_vs = 0e-6;
alpha_rho = 0e-6;

%% Save parameters
save params.mat data weight dt t0 wav vpscale vsscale rhoscale vp0 vs0 rho0 theta alpha_vp alpha_vs alpha_rho mute norm;

%% Run inversion using L-BFGS
[xo, hist, iters]=wei_lbfgs_hess('avo_asgrad', xi, 10, 200, 3, 10, 1e-4, .9, 1e-4, 1e-12, 1e-16);

%% Plot results
vpf = vp0+xo(1:Nl,end)*vpscale;
vsf = vs0+xo(1+Nl:2*Nl,end)*vsscale;
rhof = rho0+xo(1+2*Nl:3*Nl,end)*rhoscale;

% Load true model
load True_model.mat vp vs rho;

[~,modf] = fwmod(xo(:,end));
[~,mod0] = fwmod(xi);

%% Save data wiggles overlaid on figures

theta = [1, 6 12 18 24 30]; w = repmat(wav, 1, 7); t = time*1000;

gather1 = data(:,1); gather6 = data(:,6); gather12 = data(:,12); gather18 = data(:,18);
gather24 = data(:,24); gather30 = data(:,30); 
gather = [gather1,gather6,gather12,gather18,gather24,gather30];

gather01 = mod0(:,1); gather06 = mod0(:,6); gather012 = mod0(:,12); gather018 = mod0(:,18);
gather024 = mod0(:,24); gather030 = mod0(:,30); 
gather0 = [gather01,gather06,gather012,gather018,gather024,gather030];

gatherf1 = modf(:,1); gatherf6 = modf(:,6); gatherf12 = modf(:,12); gatherf18 = modf(:,18);
gatherf24 = modf(:,24); gatherf30 = modf(:,30); 
gatherf = [gatherf1,gatherf6,gatherf12,gatherf18,gatherf24,gatherf30];

%% RMSE

% RMSE Gather
rmse_seis = ((gather-gatherf)./gather).^2; rmse_seis = sqrt(rmse_seis); 
rmse_seis = rmse_seis.*100;

%% FIGURES

figure(1), % Model parameters
subplot(1,3,1); plot(vp,t,'--b','linewidth', 1.5); set(gca,'YDir','Rev'); axis([1780 2380 1 621]); 
                xlabel('VP (m/s)'); ylabel('Time (ms)'); set(gca, 'fontsize', 11); 
                set(gca,'xtick',[1780:300:2380]); set(gca,'ytick',[800:200:1400]); grid on; set(gca,'GridLineStyle','--');
subplot(1,3,2); plot(vs,t,':r','linewidth', 1.5); set(gca,'YDir','Rev'); axis([350 890 1 621]); 
                xlabel('VS (m/s)'); set(gca, 'fontsize', 11); set(gca,'yticklabel',{[]}); 
                set(gca,'xtick',[350:270:890]); set(gca,'ytick',[800:200:1400]); grid on; set(gca,'GridLineStyle','--');
subplot(1,3,3); plot(rho,t,'.b','linewidth', 1.5); set(gca,'YDir','Rev'); axis([2010 2180 1 621]); 
                xlabel('Rho (kg/m^3)'); set(gca, 'fontsize', 11); set(gca,'yticklabel',{[]}); 
                set(gca,'xtick',[2010:85:2180]); set(gca,'ytick',[800:200:1400]); grid on; set(gca,'GridLineStyle','--');
set(gca, 'fontsize', 11); set(gcf, 'position', [600         385        650         520]);

figure(2), % Noise Free gather
wiggle(gather,t,theta); 
xlabel('Angle (deg)'); ylabel('Time (ms)'); ylim([1 621]); xlim([-6.1 35.4]); set(gca, 'fontsize', 11,'fontweight','normal');
set(gcf, 'position', [600         285        285         550]); title('Noise free gather','fontweight','normal'); 

figure(3), % 1D Inversion results VP, VS, Rho
h1=subplot(1,3,1);
plot(vp,t, 'r', 'displayname', 'True model','linewidth', 1.5);
hold on, plot(vp0,t, '--b', 'displayname', 'Initial model','linewidth', 1.5);
hold on, plot(vpf,t,'--', 'color',[0 0.5 0], 'displayname', 'Inverted model', 'linewidth', 2);
set(gca, 'fontsize', 11,'fontweight','normal'); grid on; set(gca,'GridLineStyle','--');
ylabel('Time (ms)');
xlabel('VP (m/s)'); 
set(gca,'ydir','reverse');
axis([1780 2380, 1 621]);
set(gca,'xtick',[1780:300:2380]);
set(gca,'ytick',[800:200:1400]);

h2=subplot(1,3,2);
plot(vs,t, 'r', 'displayname', 'True model', 'linewidth', 1.5);
hold on, plot(vs0,t, '--b', 'displayname', 'Initial model', 'linewidth', 1.5);
hold on, plot(vsf,t, '--', 'color',[0 0.5 0], 'displayname', 'Inverted model', 'linewidth', 2);
set(gca, 'fontsize', 11,'fontweight','normal'); grid on; set(gca,'GridLineStyle','--');
xlabel('VS (m/s)'); 
set(gca,'ydir','reverse');
axis([350 890, 1 621]);
set(gca,'xtick',[350:270:890]); set(gca,'ytick',[800:200:1400]); ylabel([]);

h3=subplot(1,3,3);
plot(rho,t, 'r', 'displayname', 'True model', 'linewidth', 1.5);
hold on, plot(rho0,t, '--b', 'displayname', 'Initial model', 'linewidth', 1.5);
hold on, plot(rhof,t, '--', 'color',[0 0.5 0], 'displayname', 'Inverted model', 'linewidth', 2);
set(gca, 'fontsize', 11,'fontweight','normal'); grid on; set(gca,'GridLineStyle','--');
xlabel('(Rho kg/m^3)'); 
set(gca,'ydir','reverse');
axis([2010 2180, 1 621]);
set(gca,'xtick',[2010:85:2180]); set(gca,'ytick',[800:200:1400]);
set(gcf, 'position', [100         285        600         550]);

 
% add legend
Lgnd = legend('show');
Lgnd.Position(1) = 0.40;
Lgnd.Position(2) = 0.82;
set(gca, 'fontsize', 11);
legend('Orientation','horizontal');


figure(4), % 
subplot(141); imagesc(theta,t,gather); colormap(redblue); hold on; wiggle(gather,t,theta); 
    xlabel('Angle (deg)'); ylabel('Time (ms)'); ylim([1 621]); xlim([-6.1 36]); set(gca, 'fontsize', 11);
    set(gcf, 'position', [100         285        1020         620]); title('True model','fontweight','normal'); 
subplot(142); imagesc(theta,t,gather0); colormap(redblue); hold on; wiggle(gather0,t,theta); 
    xlabel('Angle (deg)'); ylabel('Time (ms)'); ylim([1 621]); xlim([-6.1 36]); set(gca, 'fontsize', 11);
    set(gcf, 'position', [100         285        1020         620]); title('Initial model','fontweight','normal'); ylabel([]);
subplot(143); imagesc(theta,t,gatherf); colormap(redblue); hold on; wiggle(gatherf,t,theta);  
    xlabel('Angle (deg)'); ylabel('Time (ms)'); colorbar; ylim([1 621]); xlim([-6.1 36]); set(gca, 'fontsize', 11); ylabel([]);
    set(gcf, 'position', [100         285        1020         620]); title('Inverted model','fontweight','normal');
    h=colorbar; set(h, 'Position', [.9200 .11 .0231 .8150]);
subplot(144); wiggle(rmse_seis,t,theta); xlabel('Angle (deg)'); ylabel('Time (ms)'); ylim([1 621]); ylabel([]);
    xlim([-4 34]); set(gca, 'fontsize', 11); title('% RMSE','fontweight','normal');
    


toc

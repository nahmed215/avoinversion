clear all
close all
norm = 1;
path(path,'../Functions')
tic

%% Load model (P wave velocity S wave velocity and density)

load exampledata.mat vp vs rho;
nt = length(vp);

vp0 = vp;
vs0 = vs;
rho0 = rho;
%% Wavelet
dt = 1e-3;
t0 = 0.1;
f0 = 25;
[wav, time] = Ricker(f0, t0, nt, dt);
wav = wav.';
t0_i = floor(t0/dt);

%% Incident Angles (radian to degree)
angles = 0:30;         
theta = angles*pi/180;

%% Make weight
weight = ones(nt,1);
data = weight;

%% Inversion starts here
kvp = 100.0;
kvs = 100.0;
krho = 50.0;

vpscale=kvp; 
vsscale=kvs; 
rhoscale=krho;

xi = zeros(3*nt,1); 

%% Parameters
norm = 1;
mute = ones(size(vp0));

alpha_vp = 0e-6;
alpha_vs = 0e-6;
alpha_rho = 0e-6;

%% Save params
save params.mat data weight dt t0 wav vpscale vsscale rhoscale vp0 vs0 rho0 theta alpha_vp alpha_vs alpha_rho mute norm;

%% Forward model 
[~,data] = fwmod(xi);

%% Save true model
save True_model vp vs rho;

%% Save initial model

vp0 = vel_smoother(vp0,128, .025, 1);
vs0 = vel_smoother(vs0,128, .025, 1);
rho0 = vel_smoother(rho0,128, .025, 1);

save Initial_model.mat vp0 vs0 rho0;


%% Save data and wavelet
save Observed_data.mat data angles time;

save Wavelet.mat wav dt t0;
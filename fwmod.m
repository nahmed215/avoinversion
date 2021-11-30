function [f, synth] = fwmod(x)
load params.mat data weight dt t0 wav vpscale vsscale rhoscale vp0 vs0 rho0 theta alpha_vp alpha_vs alpha_rho mute norm;
%load params.mat data weight dt t0 wav vpmin vpmax vsmin vsmax rhomin rhomax vpscale vsscale rhoscale vp0 vs0 rho0 theta alpha_vp alpha_vs alpha_rho mute norm;
Nl = (length(x))/3;
vp = vp0 + vpscale*x(1:Nl);
vs = vs0 + vsscale*x(1+Nl:2*Nl);
rho = rho0 + rhoscale*x(1+2*Nl:3*Nl);

%vp = vp0 + (vpmin +vpmax*exp(vpscale*x(1:Nl)))./(1+exp(vpscale*x(1:Nl)));
%vs = vs0 + (vsmin +vsmax*exp(vsscale*x(1+Nl:2*Nl)))./(1+exp(vsscale*x(1+Nl:2*Nl)));
%rho = rho0 + (rhomin +rhomax*exp(rhoscale*x(1+2*Nl:3*Nl)))./(1+exp(rhoscale*x(1+2*Nl:3*Nl)));
t0_i = floor(t0/dt);

%% Data error and grad
trc = Aki_Richards(vp,vs,rho,theta);
synth = conv2(wav, trc);
synth = synth(1+t0_i:Nl+t0_i,:);
J = 0.5*sum(sum((data - synth).^2));

%% Model error and grad
%display('Computing Regularization gradient');
fm = 0.5*alpha_vp*sum((vp - vp0).^2); 
fm = fm + 0.5*alpha_vs*sum((vs-vs0).^2);
fm = fm + 0.5*alpha_rho*sum((rho-rho0).^2);


%% Scaling and Muting
%display('Scaling and muting');
f = J + fm;
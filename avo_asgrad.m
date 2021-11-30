function [f, g] = avo_asgrad(x)
load params.mat data weight dt t0 wav vpscale vsscale rhoscale vp0 vs0 rho0 theta alpha_vp alpha_vs alpha_rho mute norm;
Nl = (length(x))/3;
vp = vp0 + vpscale*x(1:Nl);
vs = vs0 + vsscale*x(1+Nl:2*Nl);
rho = rho0 + rhoscale*x(1+2*Nl:3*Nl);

t0_i = floor(t0/dt);


%% Data error and grad
trc = Aki_Richards(vp,vs,rho,theta);
synth = conv2(wav, trc);
synth = synth(1+t0_i:Nl+t0_i,:);
J = 0.5*sum(sum((data - synth).^2));

res = data-synth;
res = conv2(wav, res);
res = res(1+t0_i:Nl+t0_i,:);
g=zeros(3*Nl,1);

sin2 = sin(theta).^2;
tan2 = tan(theta).^2;
for i = 2:Nl-1
    dvp0 = vp(i) - vp(i-1);
    dvp1 = vp(i+1) - vp(i);
    dvs0 = vs(i) - vs(i-1);
    dvs1 = vs(i+1) - vs(i);
    drho0 = rho(i) - rho(i-1);
    drho1 = rho(i+1) - rho(i);
    vpm0 = (vp(i) + vp(i-1))/2;
    vpm1 = (vp(i+1) + vp(i))/2;
    vsm0 = (vs(i) + vs(i-1))/2;
    vsm1 = (vs(i+1) + vs(i))/2;
    rhom0 = (rho(i) + rho(i-1))/2;    
    rhom1 = (rho(i+1) + rho(i))/2;    
    % Vp gradient
    g(i) = sum(0.5*(1 + (tan2)).*(1/vpm1 + dvp1/(2*(vpm1^2))).*res(i,:),2);
    g(i) = g(i,:) + sum(0.5*(1 + (tan2)).*(-1/vpm0 + dvp0/(2*(vpm0^2))).*res(i-1,:),2);
    g(i) = g(i,:) - sum((4*((vsm1^2)/(vpm1^3)).*(dvs1/vsm1).*sin2.*res(i,:)),2);
    g(i) = g(i,:) - sum((4*((vsm0^2)/(vpm0^3)).*(dvs0/vsm0).*sin2.*res(i-1,:)),2);
    g(i) = g(i,:) - sum((2*((vsm1^2)/(vpm1^3)).*(drho1/rhom1).*sin2.*res(i,:)),2);
    g(i) = g(i,:) - sum((2*((vsm0^2)/(vpm0^3)).*(drho0/rhom0).*sin2.*res(i-1,:)),2);
    % Vs gradient
    g(i+Nl) = g(i+Nl) - sum(4*((vsm1)/(vpm1^2)).*sin2.*res(i,:),2);
    g(i+Nl) = g(i+Nl) + sum(4*((vsm0)/(vpm0^2)).*sin2.*res(i-1,:),2);
    g(i+Nl) = g(i+Nl) + sum((2*((vsm1)/(vpm1^2))*(dvs1/vsm1).*sin2.*res(i,:)),2);
    g(i+Nl) = g(i+Nl) + sum((2*((vsm0)/(vpm0^2))*(dvs0/vsm0).*sin2.*res(i-1,:)),2);
    g(i+Nl) = g(i+Nl) - sum((2*((vsm1)/(vpm1^2)).*(drho1/rhom1).*sin2.*res(i,:)),2);
    g(i+Nl) = g(i+Nl) - sum((2*((vsm0)/(vpm0^2)).*(drho0/rhom0).*sin2.*res(i-1,:)),2);
    % Rho gradient
    g(i+2*Nl) = g(i+2*Nl) + sum((0.5*(1 -(4*(vsm1/vpm1)^2).*sin2)*(1/rhom1 + drho1/(2*(rhom1^2))).*res(i,:)),2);
    g(i+2*Nl) = g(i+2*Nl) + sum((0.5*(1 -(4*(vsm0/vpm0)^2).*sin2)*(-1/rhom0 + drho0/(2*(rhom0^2))).*res(i-1,:)),2);    
end

g(1) = g(2);
g(Nl) = g(Nl-1);
g(1+Nl) = g(2+Nl);
g(2*Nl) = g(2*Nl-1);
g(1+2*Nl) = g(2+2*Nl);
g(3*Nl) = g(3*Nl-1);


%% Model error and grad
%display('Computing Regularization gradient');
fm = 0.5*alpha_vp*sum((vp - vp0).^2); 
fm = fm + 0.5*alpha_vs*sum((vs-vs0).^2);
fm = fm + 0.5*alpha_rho*sum((rho-rho0).^2);
g(1:Nl) = g(1:Nl) + alpha_vp*(vp - vp0);
g(1+Nl:2*Nl) = g(1+Nl:2*Nl) + alpha_vs*(vs - vs0);
g(1+2*Nl:3*Nl) = g(1+2*Nl:3*Nl) + alpha_rho*(rho-rho0);


%% Scaling and Muting
%display('Scaling and muting');
g(1:Nl)=g(1:Nl)*vpscale;
g(1+Nl:2*Nl)=g(1+Nl:2*Nl)*vsscale;
g(1+2*Nl:3*Nl)=g(1+2*Nl:3*Nl)*rhoscale;

g(1:Nl)=g(1:Nl).*mute;
g(1+Nl:2*Nl)=g(1+Nl:2*Nl).*mute;
g(1+2*Nl:3*Nl)=g(1+2*Nl:3*Nl).*mute;

f = J + fm;
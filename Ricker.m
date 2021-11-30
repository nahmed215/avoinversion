function [S,t]=Ricker(f0, t0, Nt, dt)

t=(0:Nt-1)*dt;
S=(1-2*(pi^2)*(f0^2)*((t-t0).^2)).*exp(-(pi^2)*f0*f0*(t-t0).^2); 

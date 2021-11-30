function Rpp = Aki_Richards(vp, vs, rho, theta)


Nt = length(vp);
Ntheta = length(theta);
Rpp = zeros(Nt,Ntheta);
sin2 = sin(theta).^2;
tan2 = tan(theta).^2;
for i = 1:Nt-1
    dvp = vp(i+1) - vp(i);
    dvs = vs(i+1) - vs(i);
    drho = rho(i+1) - rho(i);
    vpm = (vp(i+1) + vp(i))/2;
    vsm = (vs(i+1) + vs(i))/2;
    rhom = (rho(i+1) + rho(i))/2;    
    Rpp(i,:) = 0.5*(1 + (tan2)).*(dvp/vpm) - 4*(((vsm/vpm)^2)*dvs/vsm).*sin2 + 0.5*(1 -4*((vsm/vpm)^2).*sin2)*(drho/rhom);        
end
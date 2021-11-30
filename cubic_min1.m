function [ cm ] = cubic_min1(u, fu, du, v, fv, dv )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    d = (v) - (u); 
    theta = ((fu) - (fv)) * 3 / d + (du) + (dv); 
    p = abs(theta); 
    q = abs(du);
    r = abs(dv); 
    s = max([p, q, r]); 
 
    a = theta / s; 
    gamma = s * sqrt(a * a - ((du) / s) * ((dv) / s)); 
    if ((v) < (u)) 
        gamma = -gamma;
    end
    p = gamma - (du) + theta; 
    q = gamma - (du) + gamma + (dv); 
    r = p / q; 
    cm = (u) + r * d;

end


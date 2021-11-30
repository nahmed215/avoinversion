function [ cm ] = cubic_min2(u, fu, du, v, fv, dv, xmin, xmax)

d = (v) - (u);
theta = ((fu) - (fv)) * 3 / d + (du) + (dv);
p = abs(theta);
q = abs(du);
r = abs(dv);
s = max([p, q, r]);
%     /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */
a = theta / s;
gamma = s * sqrt(max(0, a * a - ((du) / s) * ((dv) / s)));
if ((u) < (v))
    gamma = -gamma;
end
p = gamma - (dv) + theta;
q = gamma - (dv) + gamma + (du);
r = p / q;
if (r < 0. && gamma ~= 0.)
    cm = (v) - r * d;
else
    if (a < 0)
        cm = (xmax);
    else
        cm = (xmin);
    end
end


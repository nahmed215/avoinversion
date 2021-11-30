function [ qm ] = quad_min2(u, du, v, dv)
a = (u) - (v);
qm = (v) + (dv) / ((dv) - (du)) * a;
end
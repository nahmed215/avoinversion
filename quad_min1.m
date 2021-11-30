function [ qm ] = quad_min1(u, fu, du, v, fv)
a = (v) - (u);
qm = (u) + (du) / (((fu) - (fv)) / a + (du)) / 2 * a;
end
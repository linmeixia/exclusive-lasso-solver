%% Analytical Solution for the following optimization problem with vector input w
%% argmin_x 1/2||x - a||^2 + \frac{rho}{2}\|weight \circ x\|_1^2
%% analytical solution

function [x, sign_info, act_idx] = prox_weighted_l12norm(a, rho, weight)
sign_info = sign(a);
a = abs(a);
[~,idx] = sort(a./weight, 'descend');
a2 = a(idx);
weight2 = weight(idx);
s = cumsum(a2.*weight2);
L = cumsum(weight2.^2);
alptmp = s./(1/rho+L);
alp = max(alptmp);
x = max(a-alp*weight,0);
act_idx = (x==0);
x = x.*sign_info;
end

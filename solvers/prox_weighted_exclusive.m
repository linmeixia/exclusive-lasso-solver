%% Compute the proximal mapping of the weighted exclusive lasso
function [x,info] = prox_weighted_exclusive(v, n, rho, weight, group_info)
v = v(group_info.P);
weight = weight(group_info.P);
M = group_info.M;
[~,m] = size(M);
info.rr1 = zeros(n,1);
info.D = zeros(n,1);
cellx = cell(m,1);
cellrr1 = cell(m,1);
cellD = cell(m,1);
for i = 1:1:m
    tmp = v(M(1,i):M(2,i));
    tmp_w = weight(M(1,i):M(2,i));
    [cellx{i,1},cellD{i,1},cellrr1{i,1}] = prox_weighted_l12norm(tmp,2*rho, tmp_w);
end
x = cell2mat(cellx);
info.rr1 = cell2mat(cellrr1);
info.D = cell2mat(cellD);
x = x(group_info.PT);
end
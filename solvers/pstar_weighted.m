function y = pstar_weighted(x,lambda,weight,group_info)
y = 0;
x = x./weight;
x = x(group_info.P);
M = group_info.M;
[~,m] = size(M);
for i = 1:m
    xtmp = abs(x(M(1,i):M(2,i)));
    y = y+(max(xtmp))^2;
end
y = y/(4*lambda);
end
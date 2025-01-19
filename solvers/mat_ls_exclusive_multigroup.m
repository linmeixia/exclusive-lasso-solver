function y = mat_ls_exclusive_multigroup(xi,par,nzcol)
AP1 = par.AP1;
D1 = par.D1;
const_vec = par.const_vec;
counts = par.counts;

tmp = (xi'*AP1)';
tmp2 = par.smat*(D1.*tmp);
tmp2 = repelem(const_vec.*tmp2,counts);
tmp2 = tmp-tmp2.*D1;
y = xi + par.sigma*(AP1*tmp2);
end

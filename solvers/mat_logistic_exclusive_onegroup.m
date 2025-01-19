function y = mat_logistic_exclusive_onegroup(xi,par,nzcol)
AP1 = par.AP1;
AP2 = par.AP2;
const = par.const;

tmp1 = (xi'*AP1)';
tmp2 = xi'*AP2;

xiold = par.xiold;
b = par.b;
tmptmp = 1./(1+b.*xiold) + 1./(-b.*xiold);
y = xi.*tmptmp + par.sigma*(AP1*tmp1-(const*tmp2)*AP2);
end
function z = findzero(a,b,y)
z = 0;
maxiter = 100;
ab = a*b;
for i = 1:maxiter
    expbz = exp(b*z);
    h = z-y-ab/(1+expbz);
    if abs(h)<1e-12
        break;
    end
    hdiff = 1+a*expbz/(1+expbz)^2;
    z = z-h/hdiff; 
    if isnan(z)
        break;
    end
end
end
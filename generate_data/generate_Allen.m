%% Generate a test data for exclusive lasso
%% A \in R^{n \times p}: n samples , p features
%% k groups(group_num), each group has one useful features

function Data = generate_Allen(n,p,group_num, w1, w2, num_nonzeros)
num_per_g = p/group_num;
rr = w2.^[0:p-1];
Sigma0 = toeplitz(rr);
rrs = w1.^[0:num_per_g-1];
Sigma_s = toeplitz(rrs);
for i = 1:1:group_num
    idx = [(i-1)*num_per_g+1:i*num_per_g];
    Sigma0(idx, idx) = Sigma_s;
end
[R,indef] = chol(Sigma0);
if (indef); error('Sigma0 is not positive definite'); end
A = randn(n,p)*R;

x = zeros(p,1);
for i =1:1:group_num
    x_temp = zeros(num_per_g,1);
    idx_temp = randperm(num_per_g, num_nonzeros);
    x_temp(idx_temp) = rand(num_nonzeros,1)*10;
    x((i-1)*num_per_g+1:i*num_per_g) = x_temp;
end
group_M = zeros(2, group_num);
for i = 1:1:group_num
    group_M(1, i) = 1+num_per_g*(i-1);
    group_M(2, i) = num_per_g*i;
end
group_info.M = group_M;

%% Shuttle my variables
if group_num~=1
    group_info.PT = randperm(p);
    [~, group_info.P] = sort(group_info.PT);
    x = x(group_info.PT);
    A = A(:, group_info.PT);
else 
    group_info.P = [1:1:p];
    group_info.PT = [1:1:p];
end
y = A*x + randn(n, 1);
Data.A = A;
Data.y = y;
Data.groud_truth = x;
Data.p = p;
Data.n = n;
Data.group_info = group_info;
end

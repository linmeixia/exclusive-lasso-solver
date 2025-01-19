rng('default');
HOME = pwd;
addpath(genpath(HOME));

runNAL = 1;

n_A = 200;
group_num = 20;

p_list = 100; 

for n_nnz = 10
    for w2 =  0.3
        for i = 1:length(p_list)
            Data = generate_Allen(n_A,p_list(i)*group_num, group_num, 0.9, w2, 10);

            A = Data.A;
            b = Data.y;
            xtruth = Data.groud_truth;
            group_info = Data.group_info;

            [m,n] = size(A);
            hR  = @(x) A*x;
            hRt = @(x) (x'*A)';
            hhtmap = @(x) hR(hRt(x));
            eigsopt.issym = 1;
            Lip = eigs(hhtmap,length(b),1,'LA',eigsopt);
            fprintf('\n Lip const = %3.2e, nomrb = %3.2e ', Lip, norm(b));

            org_group = zeros(n,1);
            for j = 1:size(group_info.M,2)
                org_group(group_info.P(group_info.M(1,j):group_info.M(2,j))) = j;
            end
            group_info.org_group = org_group;
            weight = rand(n,1);

            stoptol = 1e-6;
            Ainput.A = A;
            Ainput.Amap = @(x) hR(x);
            Ainput.ATmap = @(x) hRt(x);

            lambda1_list = 0.1;

            for lambda_idx = 1:length(lambda1_list)
                lambda1 = lambda1_list(lambda_idx);

                if runNAL
                    nalop.Lip = Lip;
                    nalop.stoptol = stoptol;
                    nalop.runphaseI = 0;
                    nalop.phaseI_stoptol = 1e-2;
                    nalop.phaseI_maxiter = 300;
                    [obj_nal,x_nal,xi_nal,u_nal,info_nal] = nal_weighted_exclusivelasso(Ainput,b,n,lambda1,weight,group_info,nalop);
                    fprintf('\n ---------------------------------------\n');
                end

            end
        end
    end
end
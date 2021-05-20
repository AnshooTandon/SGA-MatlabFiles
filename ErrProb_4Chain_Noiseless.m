Mod_err_cnt = 0;
Kat_err_cnt = 0;
num_samples = 500;
num_runs = 1e5;
theta = 0.4;

Mod_rho_metric = zeros(1,num_runs);
Kat_rho_metric1 = zeros(1,num_runs);
Kat_rho_metric2 = zeros(1,num_runs);

rho = 1 - 2*theta;
thresh = 0.5*(1 + rho^2);

for ii = 1:num_runs
    [rho12, rho13, rho14, rho23, rho24, rho34] = ...
    generate_samples(num_samples, theta);
    
%     rho12 = max(0,rho12);
%     rho13 = max(0,rho13);
%     rho14 = max(0,rho14);
%     rho23 = max(0,rho23);
%     rho24 = max(0,rho24);
%     rho34 = max(0,rho34);

    rho12 = abs(rho12);
    rho13 = abs(rho13);
    rho14 = abs(rho14);
    rho23 = abs(rho23);
    rho24 = abs(rho24);
    rho34 = abs(rho34);

    Mod_rho_metric(ii) = sqrt(rho13*rho24*rho14*rho23)/(rho12*rho34);
    if Mod_rho_metric(ii) >= thresh
        Mod_err_cnt = Mod_err_cnt + 1;
    elseif rho13*rho24 >= rho12*rho34
        Mod_err_cnt = Mod_err_cnt + 1;
    elseif rho14*rho23 >= rho12*rho34
        Mod_err_cnt = Mod_err_cnt + 1;
    end
    
    Kat_rho_metric1(ii) = (rho13*rho24)/(rho12*rho34);
    Kat_rho_metric2(ii) = (rho13*rho24)/(rho14*rho23);
    if Kat_rho_metric1(ii) >= thresh
        Kat_err_cnt = Kat_err_cnt + 1;
    elseif Kat_rho_metric2(ii) <= thresh
        Kat_err_cnt = Kat_err_cnt + 1;
    end
end
fprintf('Modified Algo Error Prob = %e \n',(Mod_err_cnt/num_runs));
fprintf('Katiyar Algo Error Prob = %e \n',(Kat_err_cnt/num_runs));

figure(1)
histogram(Mod_rho_metric);

figure(2)
histogram(Kat_rho_metric1);

figure(3)
Kat_rho_metric2(Kat_rho_metric2 > 10) = 10; 
histogram(Kat_rho_metric2);
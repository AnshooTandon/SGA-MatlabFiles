rho = 0.4;
theta = (1-rho)/2; % theta = 0.3
q_vec = 0:0.01:0.48;
L = length(q_vec);
Mod_ErrExp_Noisy = zeros(1,L);
Mod_gap = zeros(1,L);

p1 = ones(16,1)*0.01;
p1(1) = p1(1) + 0.4;
p1(4) = p1(4) + 0.2;
p1(6) = p1(6) + 0.1;
p1(11) = p1(11) + 0.1;
p1(13) = p1(13) + 0.2;
p1(16) = p1(16) + 0.4;
p1 = p1/sum(p1);
Mod_gap(1) = EvalThreshold(p1,theta);

for ii = 1:L
    q = q_vec(ii);
    [P1, Mod_Val] = Modified_OptDist_Star_Noisy_func(theta,q,p1);
    p1 = zeros(16,1);
    for jj = 1:15
        p1(jj) = P1.x(jj);
    end
    p1(16) = 1 - sum(p1);
    Mod_gap(ii) = EvalThreshold(p1,theta);
    
    Mod_ErrExp_Noisy(ii) = Mod_Val;
end

figure(1)
plot(q_vec,log10(2)*Mod_ErrExp_Noisy,'b');
grid on
xlabel('Noise Crossover Prob. (q)');
ylabel('Exponent (Star tree, Rho=0.4, base 10)');
legend('Modified Algo');

figure(2)
plot(q_vec,Mod_gap);
legend('Mod gap');


function z = EvalThreshold(x,theta)
    t12_pos = x(1)+x(2)+x(3)+x(4)+x(13)+x(14)+x(15)+x(16);
    t12_neg = x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12);
    rho12 = t12_pos - t12_neg;
    
    t34_pos = x(1)+x(4)+x(5)+x(8)+x(9)+x(12)+x(13)+x(16);
    t34_neg = x(2)+x(3)+x(6)+x(7)+x(10)+x(11)+x(14)+x(15);
    rho34 = t34_pos - t34_neg;
    
    t13_pos = x(1)+x(2)+x(5)+x(6)+x(11)+x(12)+x(15)+x(16);
    t13_neg = x(3)+x(4)+x(7)+x(8)+x(9)+x(10)+x(13)+x(14);
    rho13 = t13_pos - t13_neg;
    
    t24_pos = x(1)+x(3)+x(6)+x(8)+x(9)+x(11)+x(14)+x(16);
    t24_neg = x(2)+x(4)+x(5)+x(7)+x(10)+x(12)+x(13)+x(15);
    rho24 = t24_pos - t24_neg;
        
    t14_pos = x(1)+x(3)+x(5)+x(7)+x(10)+x(12)+x(14)+x(16);
    t14_neg = x(2)+x(4)+x(6)+x(8)+x(9)+x(11)+x(13)+x(15);
    rho14 = t14_pos - t14_neg;
    
    t23_pos = x(1)+x(2)+x(7)+x(8)+x(9)+x(10)+x(15)+x(16);
    t23_neg = x(3)+x(4)+x(5)+x(6)+x(11)+x(12)+x(13)+x(14);
    rho23 = t23_pos - t23_neg;
    
    rho = 1 - 2*theta;
    thresh = (1+rho*rho)/2;
    
    rho_metric = (sqrt(rho13)*sqrt(rho24)/rho12)*(sqrt(rho14)*sqrt(rho23)/rho34);
    z = rho_metric - thresh;
end
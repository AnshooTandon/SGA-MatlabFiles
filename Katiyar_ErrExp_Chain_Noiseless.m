rho_vec = 0.99:-0.01:0.20;
theta_vec = (1-rho_vec)/2;
L = length(theta_vec);
Kat_ErrExp1_Noiseless = zeros(1,L);
Kat_ErrExp2_Noiseless = zeros(1,L);
Kat_gap1 = zeros(1,L);
Kat_gap2 = zeros(1,L);

p1 = ones(16,1)*0.01;
p1(1) = p1(1) + 0.1;
p1(16) = p1(16) + 0.1;
p1 = p1/sum(p1);

p2 = ones(16,1)*0.01;
p2(1) = p2(1) + 0.1;
p2(7) = p2(7) + 0.02;
p2(10) = p2(10) + 0.02;
p2(16) = p2(16) + 0.1;
p2 = p2/sum(p2);

for ii = 1:L
    theta = theta_vec(ii);
    Kat_gap1(ii) = EvalThreshold1(p1,theta);
    [P1, Kat_Val1] = Katiyar_OptDist_Chain_Noiseless_func1(theta,p1);
    Kat_ErrExp1_Noiseless(ii) = Kat_Val1;
    p1 = zeros(16,1);
    for jj = 1:15
        p1(jj) = P1.x(jj);
    end
    p1(16) = 1 - sum(p1);
    Kat_gap1(ii) = EvalThreshold1(p1,theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if theta > 0.5
        p2(7) = p2(7) + 0.012;
        p2(10) = p2(10) + 0.012;
    else
        p2(7) = p2(7) + 0.009;
        p2(10) = p2(10) + 0.009;
    end
    p2 = p2/sum(p2);
    Kat_gap2(ii) = EvalThreshold2(p2,theta);
    [P2, Kat_Val2] = Katiyar_OptDist_Chain_Noiseless_func2(theta,p2);
    Kat_ErrExp2_Noiseless(ii) = Kat_Val2;
    p2 = zeros(16,1);
    for jj = 1:15
        p2(jj) = P2.x(jj);
    end
    p2(16) = 1 - sum(p2);
    Kat_gap2(ii) = EvalThreshold2(p2,theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
Kat_ErrExp_Noiseless = min(Kat_ErrExp1_Noiseless,Kat_ErrExp2_Noiseless);

figure(1)
plot(rho_vec,Kat_ErrExp1_Noiseless,'r--');
grid on
hold on
plot(rho_vec,Kat_ErrExp2_Noiseless,'m:');
plot(rho_vec,Kat_ErrExp_Noiseless,'c-.');

xlabel('rho');
ylabel('Exponent (Chain tree, No noise)');
legend('Katiyar Algo');

figure(2) 
plot(rho_vec, Kat_gap1);
legend('Kat gap1');

figure(3) 
plot(rho_vec, Kat_gap2);
legend('Kat gap2');


function z1 = EvalThreshold1(x,theta)
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
    
    rho = 1 - 2*theta;
    thresh = (1+rho*rho)/2;
    
    rho_metric = (rho13*rho24)/(rho12*rho34);
    z1 = rho_metric - thresh;
end

function z2 = EvalThreshold2(x,theta)
    t14_pos = x(1)+x(3)+x(5)+x(7)+x(10)+x(12)+x(14)+x(16);
    t14_neg = x(2)+x(4)+x(6)+x(8)+x(9)+x(11)+x(13)+x(15);
    rho14 = t14_pos - t14_neg;
    
    t23_pos = x(1)+x(2)+x(7)+x(8)+x(9)+x(10)+x(15)+x(16);
    t23_neg = x(3)+x(4)+x(5)+x(6)+x(11)+x(12)+x(13)+x(14);
    rho23 = t23_pos - t23_neg;
    
    t13_pos = x(1)+x(2)+x(5)+x(6)+x(11)+x(12)+x(15)+x(16);
    t13_neg = x(3)+x(4)+x(7)+x(8)+x(9)+x(10)+x(13)+x(14);
    rho13 = t13_pos - t13_neg;
    
    t24_pos = x(1)+x(3)+x(6)+x(8)+x(9)+x(11)+x(14)+x(16);
    t24_neg = x(2)+x(4)+x(5)+x(7)+x(10)+x(12)+x(13)+x(15);
    rho24 = t24_pos - t24_neg;
    
    rho = 1 - 2*theta;
    thresh = (1+rho*rho)/2;
    
    rho_metric = (rho13*rho24)/(rho14*rho23);
    z2 = rho_metric - thresh;
end
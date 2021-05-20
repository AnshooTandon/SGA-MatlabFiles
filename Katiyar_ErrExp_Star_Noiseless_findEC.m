rho_vec_findEC = 0.999:-0.001:0.110;
theta_vec = (1-rho_vec_findEC)/2;
L = length(theta_vec);
Kat_ErrExp_Noiseless_findEC = zeros(1,L);
Kat_gap_findEC = zeros(1,L);

p2 = ones(16,1)*0.01;
p2(1) = p2(1) + 0.4;
p2(4) = p2(4) + 0.19;
p2(6) = p2(6) + 0.1;
p2(11) = p2(11) + 0.1;
p2(13) = p2(13) + 0.19;
p2(16) = p2(16) + 0.4;
p2 = p2/sum(p2);

for ii = 1:L
    theta = theta_vec(ii);
    if theta < 0.35
        p2(4) = p2(4) + 0.01;
        p2(13) = p2(13) + 0.01;
    else
         p2(4) = p2(4) + 0.003;
         p2(13) = p2(13) + 0.003;
    end
    p2 = p2./sum(p2);
    Kat_gap_findEC(ii) = EvalThreshold(p2,theta);
    
    [P2, Kat_Val] = Katiyar_OptDist_Star_Noiseless_func_findEC(theta,p2);
    p2 = zeros(16,1);
    for jj = 1:15
        p2(jj) = P2.x(jj);
    end
    p2(16) = 1 - sum(p2);
    Kat_gap_findEC(ii) = EvalThreshold(p2,theta);

    Kat_ErrExp_Noiseless_findEC(ii) = Kat_Val;
end

figure(1)
plot(rho_vec_findEC,Kat_ErrExp_Noiseless_findEC,'r--');
grid on
xlabel('rho');
ylabel('Error Exponent');
legend('Katiyar Algo');
title('Star tree, No noise, findEC');

figure(2) 
plot(rho_vec_findEC, Kat_gap_findEC);
legend('Kat gap');

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
    
    rho = 1 - 2*theta;
    thresh = (1+rho*rho)/2;
    
    rho_metric = (rho13*rho24)/(rho12*rho34);
    z = rho_metric - thresh;
end

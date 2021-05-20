%rho = 0.4;
rho = 0.74;
theta = (1-rho)/2; 
q_vec = 0:0.01:0.46;
L = length(q_vec);
Kat_ErrExp1_Noisy = zeros(1,L);
Kat_ErrExp2_Noisy = zeros(1,L);
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
    q = q_vec(ii);
    Kat_gap1(ii) = EvalThreshold1(p1,theta);
    [P1, Kat_Val1] = Katiyar_OptDist_Chain_Noisy_func1(theta,q,p1);
    Kat_ErrExp1_Noisy(ii) = Kat_Val1;
    p1 = zeros(16,1);
    for jj = 1:15
        p1(jj) = P1.x(jj);
    end
    p1(16) = 1 - sum(p1);
    Kat_gap1(ii) = EvalThreshold1(p1,theta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Kat_gap2(ii) = EvalThreshold2(p2,theta);
    [P2, Kat_Val2] = Katiyar_OptDist_Chain_Noisy_func2(theta,q,p2);
    Kat_ErrExp2_Noisy(ii) = Kat_Val2;
    p2 = zeros(16,1);
    for jj = 1:15
        p2(jj) = P2.x(jj);
    end
    p2(16) = 1 - sum(p2);
    Kat_gap2(ii) = EvalThreshold2(p2,theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
Kat_ErrExp_Noisy = min(Kat_ErrExp1_Noisy,Kat_ErrExp2_Noisy);

figure(1)
plot(q_vec,log10(2)*Kat_ErrExp_Noisy,'r--');
grid on

h = xlabel('Noise Crossover Probability, $q$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);

h = ylabel('Error Exponent (base 10)');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);

h = legend('Katiyar Algo');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',10);

%h = title('4 node Chain with $\rho=0.4$');
h = title('4 node Chain with $\rho=0.74$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',11);

figure(2) 
plot(q_vec, Kat_gap1);
legend('Kat gap1');

figure(3) 
plot(q_vec, Kat_gap2);
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
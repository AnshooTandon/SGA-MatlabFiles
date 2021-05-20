%Mod_rho_vec = 0.99:-0.01:0.20; % version v2
Mod_rho_vec = 0.99:-0.01:0.06; % version v3

theta_vec = (1-Mod_rho_vec)/2;
L = length(theta_vec);
Mod_ErrExp_Noiseless = zeros(1,L);
Mod_ErrExp_Noiseless1 = zeros(1,L);
Mod_ErrExp_Noiseless2 = zeros(1,L);
Mod_ErrExp_Noiseless3 = zeros(1,L);

Mod_gap1 = zeros(1,L);
Mod_gap2 = zeros(1,L);
Mod_gap3 = zeros(1,L);

p1 = ones(16,1)*0.01;
p1(1) = p1(1) + 0.1;
p1(16) = p1(16) + 0.1;
p1 = p1/sum(p1);

p2 = ones(16,1)*0.01;
p2(1) = p2(1) + 0.1;
p2(16) = p2(16) + 0.1;
p2(6) = p2(6) + 0.01;
p2(11) = p2(11) + 0.01;
p2 = p2/sum(p2);

p3 = ones(16,1)*0.01;
p3(1) = p3(1) + 0.1;
p3(16) = p3(16) + 0.1;
p3(7) = p3(7) + 0.01;
p3(10) = p3(10) + 0.01;
p3 = p3/sum(p3);

for ii = 1:L
    theta = theta_vec(ii);
    %Mod_gap(ii) = EvalThreshold(p1,theta);
    [P1, Mod_Val1] = Modified_OptDist_Chain_Noiseless_func(theta,p1);
    Mod_ErrExp_Noiseless1(ii) = Mod_Val1;
    p1 = zeros(16,1);
    for jj = 1:15
        p1(jj) = P1.x(jj);
    end
    p1(16) = 1 - sum(p1);
    Mod_gap1(ii) = EvalThreshold(p1,theta);

    [P2, Mod_Val2] = Modified_OptDist_Chain_Noiseless_func2(theta,p2);
    Mod_ErrExp_Noiseless2(ii) = Mod_Val2;
    p2 = zeros(16,1);
    for jj = 1:15
        p2(jj) = P2.x(jj);
    end
    p2(16) = 1 - sum(p2);
    Mod_gap2(ii) = EvalThreshold2(p2);
    
    [P3, Mod_Val3] = Modified_OptDist_Chain_Noiseless_func3(theta,p3);
    Mod_ErrExp_Noiseless3(ii) = Mod_Val3;
    p3 = zeros(16,1);
    for jj = 1:15
        p3(jj) = P3.x(jj);
    end
    p3(16) = 1 - sum(p3);
    Mod_gap3(ii) = EvalThreshold3(p3); 

    Mod_ErrExp_Noiseless(ii) = min([Mod_Val1 Mod_Val2 Mod_Val3]);
end

figure(1)
plot(Mod_rho_vec,Mod_ErrExp_Noiseless,'b');
grid on
xlabel('rho');
ylabel('Exponent (Chain tree, No noise)');
legend('Modified Algo');

figure(2)
plot(Mod_rho_vec,Mod_gap1);
legend('Mod gap1');

figure(3)
plot(Mod_rho_vec,Mod_gap2);
legend('Mod gap2');

figure(4)
plot(Mod_rho_vec,Mod_gap3);
legend('Mod gap3');

figure(5)
plot(Mod_rho_vec,log10(2)*Mod_ErrExp_Noiseless1,'k');
hold on
grid on
plot(Mod_rho_vec,log10(2)*Mod_ErrExp_Noiseless2,'b');
plot(Mod_rho_vec,log10(2)*Mod_ErrExp_Noiseless3,'r');
xlabel('Edge correlation, rho');
ylabel('Error Exponent (base 10)');
title('4 node Noiseless Chain: Modified Algo');
legend('Err Exp 1','Err Exp 2','Err Exp 3');

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

function z2 = EvalThreshold2(x)
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
            
    rho_metric = (rho13*rho24)/(rho12*rho34);
    z2 = rho_metric - 1;
end

function z3 = EvalThreshold3(x)
    t12_pos = x(1)+x(2)+x(3)+x(4)+x(13)+x(14)+x(15)+x(16);
    t12_neg = x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12);
    rho12 = t12_pos - t12_neg;
    
    t34_pos = x(1)+x(4)+x(5)+x(8)+x(9)+x(12)+x(13)+x(16);
    t34_neg = x(2)+x(3)+x(6)+x(7)+x(10)+x(11)+x(14)+x(15);
    rho34 = t34_pos - t34_neg;
    
    t14_pos = x(1)+x(3)+x(5)+x(7)+x(10)+x(12)+x(14)+x(16);
    t14_neg = x(2)+x(4)+x(6)+x(8)+x(9)+x(11)+x(13)+x(15);
    rho14 = t14_pos - t14_neg;
    
    t23_pos = x(1)+x(2)+x(7)+x(8)+x(9)+x(10)+x(15)+x(16);
    t23_neg = x(3)+x(4)+x(5)+x(6)+x(11)+x(12)+x(13)+x(14);
    rho23 = t23_pos - t23_neg;
    
    rho_metric = (rho14*rho23)/(rho12*rho34);
    z3 = rho_metric - 1;
end
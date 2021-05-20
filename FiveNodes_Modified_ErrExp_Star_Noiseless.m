rho_vec = 0.20:0.01:0.99;
theta_vec = (1-rho_vec)/2;
L = length(theta_vec);
Mod_ErrExp_Noiseless = zeros(1,L);
Mod_gap = zeros(1,L);

% Assume the 5 nodes to be X0,X1,X2,X3,X4 with X0 as the root node
p1 = ones(32,1)*0.01;
p1(1) = p1(1) + 0.40;
p1(4) = p1(4) + 0.07;
p1(29) = p1(29) + 0.07;
p1(32) = p1(32) + 0.40;
p1 = p1/sum(p1);

for ii = 1:L
    theta = theta_vec(ii);
    Mod_gap(ii) = EvalThreshold(p1,theta);
    [P1, Mod_Val] = FiveNodes_Modified_OptDist_Star_Noiseless_func(theta,p1);
    p1 = zeros(32,1);
    for jj = 1:31
        p1(jj) = P1.x(jj);
    end
    p1(32) = 1 - sum(p1);
    Mod_gap(ii) = EvalThreshold(p1,theta);
    
    Mod_ErrExp_Noiseless(ii) = Mod_Val;
end

figure(1)
plot(rho_vec,log10(2)*Mod_ErrExp_Noiseless,'b');
grid on
xlabel('rho');
ylabel('Exponent (base 10)');
legend('Modified Algo');
title('Star tree, No noise, 5 nodes');

figure(2)
plot(rho_vec,Mod_gap);
legend('Mod gap');


function z = EvalThreshold(x,theta)
     t12_pos = x(1)+x(2)+x(3)+x(4)+x(13)+x(14)+x(15)+x(16)+...
            x(17)+x(18)+x(19)+x(20)+x(29)+x(30)+x(31)+x(32);
    t12_neg = x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12)+...
            x(21)+x(22)+x(23)+x(24)+x(25)+x(26)+x(27)+x(28);
    rho12 = t12_pos - t12_neg;
    
    t34_pos = x(1)+x(4)+x(5)+x(8)+x(9)+x(12)+x(13)+x(16)+...
            x(17)+x(20)+x(21)+x(24)+x(25)+x(28)+x(29)+x(32);
    t34_neg = x(2)+x(3)+x(6)+x(7)+x(10)+x(11)+x(14)+x(15)+...
            x(18)+x(19)+x(22)+x(23)+x(26)+x(27)+x(30)+x(31);
    rho34 = t34_pos - t34_neg;
    
    t13_pos = x(1)+x(2)+x(5)+x(6)+x(11)+x(12)+x(15)+x(16)+...
            x(17)+x(18)+x(21)+x(22)+x(27)+x(28)+x(31)+x(32);
    t13_neg = x(3)+x(4)+x(7)+x(8)+x(9)+x(10)+x(13)+x(14)+...
            x(19)+x(20)+x(23)+x(24)+x(25)+x(26)+x(29)+x(30);
    rho13 = t13_pos - t13_neg;
    
    t24_pos = x(1)+x(3)+x(6)+x(8)+x(9)+x(11)+x(14)+x(16)+...
            x(17)+x(19)+x(22)+x(24)+x(25)+x(27)+x(30)+x(32);
    t24_neg = x(2)+x(4)+x(5)+x(7)+x(10)+x(12)+x(13)+x(15)+...
            x(18)+x(20)+x(21)+x(23)+x(26)+x(28)+x(29)+x(31);
    rho24 = t24_pos - t24_neg;
        
    t14_pos = x(1)+x(3)+x(5)+x(7)+x(10)+x(12)+x(14)+x(16)+...
            x(17)+x(19)+x(21)+x(23)+x(26)+x(28)+x(30)+x(32);
    t14_neg = x(2)+x(4)+x(6)+x(8)+x(9)+x(11)+x(13)+x(15)+...
            x(18)+x(20)+x(22)+x(24)+x(25)+x(27)+x(29)+x(31);
    rho14 = t14_pos - t14_neg;
    
    t23_pos = x(1)+x(2)+x(7)+x(8)+x(9)+x(10)+x(15)+x(16)+...
            x(17)+x(18)+x(23)+x(24)+x(25)+x(26)+x(31)+x(32);
    t23_neg = x(3)+x(4)+x(5)+x(6)+x(11)+x(12)+x(13)+x(14)+...
            x(19)+x(20)+x(21)+x(22)+x(27)+x(28)+x(29)+x(30);
    rho23 = t23_pos - t23_neg;
   
    rho = 1 - 2*theta;
    thresh = (1+rho*rho)/2;
    
    rho_metric = (sqrt(rho13)*sqrt(rho24)/rho12)*(sqrt(rho14)*sqrt(rho23)/rho34);
    z = rho_metric - thresh;
end
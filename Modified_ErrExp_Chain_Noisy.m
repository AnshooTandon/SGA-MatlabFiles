%rho = 0.4;
rho = 0.74;
theta = (1-rho)/2;
q_vec = 0:0.01:0.46;
L = length(q_vec);
Mod_ErrExp_Noisy = zeros(1,L);
Mod_gap = zeros(1,L);
Mod_gap2 = zeros(1,L);
Mod_gap3 = zeros(1,L);


p1 = ones(16,1)*0.01;
p1(1) = p1(1) + 0.07;
p1(16) = p1(16) + 0.07;
p1 = p1/sum(p1);

% p2 = ones(16,1)*0.01;
% p2(1) = p2(1) + 0.1;
% p2(16) = p2(16) + 0.1;
% p2(6) = p2(6) + 0.01;
% p2(11) = p2(11) + 0.01;
% p2 = p2/sum(p2);
% 
% p3 = ones(16,1)*0.01;
% p3(1) = p3(1) + 0.1;
% p3(16) = p3(16) + 0.1;
% p3(7) = p3(7) + 0.01;
% p3(10) = p3(10) + 0.01;
% p3 = p3/sum(p3);
% 
for ii = 1:L
    q = q_vec(ii);
    
%    Mod_gap(ii) = EvalThreshold(p1,theta);
   [P1, Mod_Val] = Modified_OptDist_Chain_Noisy_func(theta,q,p1);
    p1_new = zeros(16,1);
    for jj = 1:15
        p1_new(jj) = P1.x(jj);
    end
    p1_new(16) = 1 - sum(p1_new);
    Mod_gap(ii) = EvalThreshold(p1_new,theta);
 
%     Mod_gap2(ii) = EvalThreshold2(p2);
%     [P2, Mod_Val2] = Modified_OptDist_Chain_Noisy_func2(theta,q,p2);
%     p2_new = zeros(16,1);
%     for jj = 1:15
%         p2_new(jj) = P2.x(jj);
%     end
%     p2_new(16) = 1 - sum(p2_new);
%     Mod_gap2(ii) = EvalThreshold2(p2_new);
%  
%     Mod_gap3(ii) = EvalThreshold3(p3);
%     [P3, Mod_Val3] = Modified_OptDist_Chain_Noisy_func3(theta,q,p3);
%     p3_new = zeros(16,1);
%     for jj = 1:15
%         p3_new(jj) = P3.x(jj);
%     end
%     p3_new(16) = 1 - sum(p3_new);
%     Mod_gap3(ii) = EvalThreshold3(p3_new); 
% 
%     Mod_ErrExp_Noisy(ii) = min([Mod_Val Mod_Val2 Mod_Val3]);
     Mod_ErrExp_Noisy(ii) = Mod_Val;
end

figure(1)
plot(q_vec,log10(2)*Mod_ErrExp_Noisy,'b');
grid on
xlabel('Noise Crossover Prob. (q)');
%ylabel('Exponent (Star tree, Rho=0.4, base 10)');
ylabel('Exponent (Star tree, Rho=0.74, base 10)');
legend('Modified Algo');

figure(2)
plot(q_vec,Mod_gap);
legend('Mod gap');

% figure(3)
% plot(q_vec,Mod_gap2);
% legend('Mod gap2');
% 
% figure(4)
% plot(q_vec,Mod_gap3);
% legend('Mod gap3');


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

% function z2 = EvalThreshold2(x)
%     t12_pos = x(1)+x(2)+x(3)+x(4)+x(13)+x(14)+x(15)+x(16);
%     t12_neg = x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12);
%     rho12 = t12_pos - t12_neg;
%     
%     t34_pos = x(1)+x(4)+x(5)+x(8)+x(9)+x(12)+x(13)+x(16);
%     t34_neg = x(2)+x(3)+x(6)+x(7)+x(10)+x(11)+x(14)+x(15);
%     rho34 = t34_pos - t34_neg;
%     
%     t13_pos = x(1)+x(2)+x(5)+x(6)+x(11)+x(12)+x(15)+x(16);
%     t13_neg = x(3)+x(4)+x(7)+x(8)+x(9)+x(10)+x(13)+x(14);
%     rho13 = t13_pos - t13_neg;
%     
%     t24_pos = x(1)+x(3)+x(6)+x(8)+x(9)+x(11)+x(14)+x(16);
%     t24_neg = x(2)+x(4)+x(5)+x(7)+x(10)+x(12)+x(13)+x(15);
%     rho24 = t24_pos - t24_neg;
%             
%     rho_metric = (rho13*rho24)/(rho12*rho34);
%     z2 = rho_metric - 1;
% end
% 
% function z3 = EvalThreshold3(x)
%     t12_pos = x(1)+x(2)+x(3)+x(4)+x(13)+x(14)+x(15)+x(16);
%     t12_neg = x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12);
%     rho12 = t12_pos - t12_neg;
%     
%     t34_pos = x(1)+x(4)+x(5)+x(8)+x(9)+x(12)+x(13)+x(16);
%     t34_neg = x(2)+x(3)+x(6)+x(7)+x(10)+x(11)+x(14)+x(15);
%     rho34 = t34_pos - t34_neg;
%     
%     t14_pos = x(1)+x(3)+x(5)+x(7)+x(10)+x(12)+x(14)+x(16);
%     t14_neg = x(2)+x(4)+x(6)+x(8)+x(9)+x(11)+x(13)+x(15);
%     rho14 = t14_pos - t14_neg;
%     
%     t23_pos = x(1)+x(2)+x(7)+x(8)+x(9)+x(10)+x(15)+x(16);
%     t23_neg = x(3)+x(4)+x(5)+x(6)+x(11)+x(12)+x(13)+x(14);
%     rho23 = t23_pos - t23_neg;
%     
%     rho_metric = (rho14*rho23)/(rho12*rho34);
%     z3 = rho_metric - 1;
% end
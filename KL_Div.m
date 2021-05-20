theta_vec = 0.05:0.05:0.45;
q_vec = 0:0.01:0.49;
KLdiv4 = zeros(length(theta_vec),length(q_vec));
KL_Bound4 = zeros(length(theta_vec),length(q_vec));

for ii = 1:length(theta_vec)
    theta = theta_vec(ii);
    for jj = 1:length(q_vec)
        q = q_vec(jj);
        phi = theta*(1-q) + (1-theta)*q;

        P0_q = zeros(16,1);
        P1_q = zeros(16,1);
        P0_q(1) = 0.5*(1-theta)^2*(1-phi);
        P0_q(2) = 0.5*(1-theta)^2*phi;
        P0_q(3) = 0.5*theta*(1-theta)*(1-phi);
        P0_q(4) = 0.5*theta*(1-theta)*phi;
        P0_q(5) = 0.5*theta*(1-theta)*(1-phi);
        P0_q(6) = 0.5*theta*(1-theta)*phi;
        P0_q(7) = 0.5*theta^2*(1-phi);
        P0_q(8) = 0.5*theta^2*phi;
        for kk = 9:16
            P0_q(kk) = P0_q(17-kk);
        end
        P1_q(1) = 0.5*(1-theta)^2*(1-phi);
        P1_q(2) = 0.5*(1-theta)^2*phi;
        P1_q(3) = 0.5*theta*(1-theta)*phi;
        P1_q(4) = 0.5*theta*(1-theta)*(1-phi);
        P1_q(5) = 0.5*theta*(1-theta)*(1-phi);
        P1_q(6) = 0.5*theta*(1-theta)*phi;
        P1_q(7) = 0.5*theta^2*phi;
        P1_q(8) = 0.5*theta^2*(1-phi);
        for kk = 9:16
            P1_q(kk) = P1_q(17-kk);
        end
    
        KLdiv4(ii,jj) = sum(P1_q.*log2(P1_q./P0_q));
    end
end
    
for ii = 1:length(theta_vec)
    KL_Bound4(ii,:) = KLdiv4(ii,1)*(1-2*q_vec);
end

idx4 = (KLdiv4 > KL_Bound4);
err4 = sum(sum(idx4))

for kk = 1:length(theta_vec);
    close all
    plot(q_vec,KLdiv4(kk,:),'b');
    hold on
    plot(q_vec,KL_Bound4(kk,:),'r');
    xlabel('BSC crossover probability, q');
    ylabel('KL divegence');
    legend('KL Div--Noisy','KL Bound Noiseless*(1-2q)'); 
end
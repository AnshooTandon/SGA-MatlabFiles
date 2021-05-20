load ModErrExpStarNoisyRhoPoint4_v2.mat
load KatErrExpStarNoisyRhoPoint4_v2.mat

figure(1)
plot(q_vec,log10(2)*Mod_ErrExp_Noisy,'b','Linewidth',1.5);
grid on
hold on
plot(q_vec,log10(2)*Kat_ErrExp_Noisy,'r--','Linewidth',1.5);
h = xlabel('Noise Crossover Probability, $q$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);

h = ylabel('Error Exponent (base 10)');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);

h = legend('Modified Algo','Katiyar Algo');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',10);

h = title('4 node Star with edge correlation $\rho=0.4$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',11);
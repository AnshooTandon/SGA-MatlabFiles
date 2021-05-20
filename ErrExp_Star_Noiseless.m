load ModErrExpStarNoiseless_v2.mat
load KatErrExpStarNoiseless_findEC.mat

figure(1)
plot(rho_vec,log10(2)*Mod_ErrExp_Noiseless,'b','Linewidth',1.5);
grid on
hold on
plot(rho_vec_findEC,log10(2)*Kat_ErrExp_Noiseless_findEC,'r--','Linewidth',1.5);
h = xlabel('Edge correlation $\rho$');
set(h,'unit','character','interpreter','latex')
set(h,'FontSize',12);

h = ylabel('Error Exponent (base 10)');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);

h = legend('Modified Algo','Katiyar Algo');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',10);

h = title('4-node Noiseless Star, findEC');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',11);
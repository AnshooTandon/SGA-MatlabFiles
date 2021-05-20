load ModErrExpStarNoiseless_v2.mat
%load KatErrExpStarNoiseless_findEC.mat
load KatErrExpStarNoiseless_v2.mat

subplot(1,2,1)
plot(rho_vec,log(2)*Mod_ErrExp_Noiseless,'b','Linewidth',1.5);
grid on
hold on
plot(rho_vec,log(2)*Kat_ErrExp_Noiseless,'r--','Linewidth',1.5);
h = xlabel('Edge correlation, $\rho$');
set(h,'unit','character','interpreter','latex')
set(h,'FontSize',14);
h = ylabel('Error Exponent');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',14);

h = legend('$E(\Psi_{\mathrm{SGA}}, \tilde{P})$','$E(\Psi_{\mathrm{KA}}, \tilde{P})$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);
h = title('(a) \, $q_{\mathrm{max}} = 0$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',13);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ModErrExpStarNoisyRhoPoint4_v2.mat
load KatErrExpStarNoisyRhoPoint4_v2.mat
subplot(1,2,2)
plot(q_vec,log(2)*Mod_ErrExp_Noisy,'b','Linewidth',1.5);
grid on
hold on
plot(q_vec,log(2)*Kat_ErrExp_Noisy,'r--','Linewidth',1.5);
h = xlabel('$q_{\mathrm{max}}$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',15);
% h = ylabel('Error Exponent');
% set(h,'unit','character','interpreter','latex');
% set(h,'FontSize',14);
h = legend('$E(\Psi_{\mathrm{SGA}}, \tilde{P})$','$E(\Psi_{\mathrm{KA}}, \tilde{P})$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);
h = title('(b) \, $\rho=0.4$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',13);
%axis([0 0.5 0 4e-3]);
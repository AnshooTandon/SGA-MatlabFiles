load ModErrExpChainNoiseless_v2.mat
load KatErrExpChainNoiseless_v2.mat

figure(1)
%plot(rho_vec,log10(2)*Mod_ErrExp_Noiseless,'b','Linewidth',1.5);
plot(rho_vec,log(2)*Mod_ErrExp_Noiseless,'b','Linewidth',1.5);
grid on
hold on
% plot(rho_vec,log10(2)*Kat_ErrExp1_Noiseless,'r--','Linewidth',1);
% plot(rho_vec,log10(2)*Kat_ErrExp2_Noiseless,'m:','Linewidth',1);
% plot(rho_vec,log10(2)*Kat_ErrExp_Noiseless,'c-.','Linewidth',2);
plot(rho_vec,log(2)*Kat_ErrExp1_Noiseless,'r--','Linewidth',1);
plot(rho_vec,log(2)*Kat_ErrExp2_Noiseless,'m:','Linewidth',1);
plot(rho_vec,log(2)*Kat_ErrExp_Noiseless,'c-.','Linewidth',2);

h = xlabel('Edge correlation, $\rho$');
set(h,'unit','character','interpreter','latex')
set(h,'FontSize',12);

%h = ylabel('Error Exponent (base 10)');
h = ylabel('Error Exponent');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',12);

%h = legend('Modified Algo','Katiyar Algo: Exp1','Katiyar Alo: Exp2','min(Exp1,Exp2)');
h = legend('Modified Algo: $e_{\mathrm{MA}}$','Katiyar Algo: $e_1$','Katiyar Algo: $e_2$','$e_{\mathrm{KA}} = \min(e_1,e_2)$');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',10);

h = title('4 node Noiseless Chain');
set(h,'unit','character','interpreter','latex');
set(h,'FontSize',11);
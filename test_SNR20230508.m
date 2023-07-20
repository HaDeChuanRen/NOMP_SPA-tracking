clc; clear; close all;

P_CW = 1.78e-2;
G_t = 1;
G_r = 1;
lambda_cw = 4e-3;
sigma_RCS = 1;
r_range = 1 : 0.1 : 10;
L_loss = 1;
k_Boltz = 1.38e-23;
T_kel = 300;
F_R = 10^1.5;
B_IF = 5e6;
N = 256;

SNR_r = 10 * log10(N * (P_CW * G_t * G_r * lambda_cw^2 * sigma_RCS) ./ ...
    ((4 * pi)^3 * r_range .^ 4 * L_loss * k_Boltz * T_kel * F_R * B_IF));


Fsz = 12;
Lw = 2;
Msz = 2;

figure;
plot(r_range, SNR_r, 'LineWidth', Lw, 'MarkerSize', Msz);
xlabel('Range (m)', 'Fontsize', Fsz);
ylabel('$SNR_{\rm in}$', 'Fontsize', Fsz, 'Interpreter','latex');



clc; clear; close all;


% read the ADC data

orginal_path = 'C:\study\202206MNOMP_CFAR\data\20220904exp\04';
% exp_type = '\02people2';
% exp_serial = '\01';

N_start = 1;
M_start = 512; %17 1850
% L_start = 1;

% range_true = [2.95 + 0.14; 4.78 + 0.14];
% theta_true = [- 19.8; 0];
xPos_true = 0;
yPos_true = 1.9;
K_ture = length(xPos_true);
amp_true = 60 * ones(K_ture, 1);

filename = [orginal_path, '\adc_data.bin'];
data_cube = readadc(filename);

% radar parameter
c = 3e8;
T_idle = 100e-6;
T_ramp = 60e-6;
T_circle = T_idle + T_ramp;
Fre_start = 77e9;
Slope_fre = 29.982e12;
Fs = 10e6;
Ts = 1 / Fs;
lambda_cw = c / Fre_start;
Rx_interval = lambda_cw / 2;


% set the observation scene
Nx = size(data_cube, 1);
My = size(data_cube, 2);
Lz = size(data_cube, 3);
NL_num = Nx * Lz;


% ymat = squeeze(data_cube(N_start : (N_start + Nx - 1), M_start, :));

% algorithm parameter set
guard_n = 4;
guard_l = 0;
training_n = 8;
training_l = 1;
guard_training_size2D = [guard_n, guard_l, training_n, training_l];
OS_rate = [2, 8];
R_c = 5;
R_s = 10;



K_max = 32;
P_oe2D = 1e-2;
N_r2D = (2 * training_n + 1) * (2 * training_l + 1) - ...
(2 * guard_n + 1) * (2 * guard_l + 1);
% N_r2D = 40;
alpha_set2D = alpha_PoebyS(P_oe2D, NL_num, N_r2D);
guard_band = [guard_n, guard_l];
sigma_set = 10 ^ (24 / 10);
tau_set = sigma_set * chi2inv((1 - P_oe2D) ^ (1 / NL_num), 2) / 2;

% NOMP method analysis
Zkcell_seq = cat(1);

tic;

for k_time = 1 : My
    Ymat_k = squeeze(data_cube(:, k_time, :));
    [omega_fast, gain_fast, ~, Threshold_fast] = ...
    NOMP2D_fast(Ymat_k, alpha_set2D, N_r2D, K_max, guard_band, OS_rate, ...
    R_c, R_s);

    if ~isempty(omega_fast)
        omegax_fast = omega_fast(:, 1);
        omegaz_fast = wrapToPi(omega_fast(:, 2));
        range_fast = (c * omegax_fast) / (4 * pi * Ts * Slope_fre);
        theta_fast = asin((c * omegaz_fast) / (2 * pi * Fre_start * ...
            Rx_interval));
        Zkcell_seq{k_time} = [range_fast .* sin(theta_fast), ...
        range_fast .* cos(theta_fast)];
    else
        Zkcell_seq{k_time} = [];
    end
end

time_fast = toc;

% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;


fhandle_result = figure(1);
title('2DNOMP-fast result', 'Fontsize', Fsz);
xlabel('x/m', 'Fontsize', Fsz);
ylabel('y/m', 'Fontsize', Fsz);
hold on;
xlim([-5 5])
ylim([0 25])

for k_time = 1 : My
    for z_idx = 1 : size(Zkcell_seq{k_time}, 1)
        phandle_Zk = plot(Zkcell_seq{k_time}(z_idx, 1), ...
        Zkcell_seq{k_time}(z_idx, 2), 'k*', 'LineWidth', Lw, 'Marker', ...
        's', 'MarkerSize', Msz);
    end
    pause(0.01);
end



% tic;
% [omega_NOMPCFAR, gain_NOMPCFAR, ~, Threshold_NOMPCFAR] = ...
% NOMP2D_CFAR(ymat, alpha_set2D, N_r2D, K_max, guard_band);
% time_NOMPCFAR = toc;

% tic;
% [omega_tau, gain_tau, ~] = ...
% NOMP2D(ymat, tau_set, K_max);
% time_tau = toc;




% lw = 2;
% fsz = 12;
% msz = 8;
% range_max = (c * 2 * pi) / (4 * pi * Ts * Slope_fre);
% % bias = 10 * log10(NL_num);


% omegax_tau = omega_tau(:, 1);
% omegaz_tau = wrapToPi(omega_tau(:, 2));
% range_tau = (c * omegax_tau) / (4 * pi * Ts * Slope_fre);
% theta_tau = asin((c * omegaz_tau) / (2 * pi * Fre_start * Rx_interval));
% theta_tau_deg = theta_tau * 180 / pi;
% Khat_tau = length(omegax_tau);


% range_idx = linspace(0, range_max, Nx);
% theta_idx = linspace(-90, 90, Lz);
% [range_newidx, theta_newidx] = meshgrid(range_idx, theta_idx);

% figure;
% surf(range_newidx, theta_newidx, (10 * log10(tau_set)) * ones(Lz, Nx),...
% 'Linestyle', 'none');
% hold on;
% stem3(range_tau, theta_tau_deg, 20 * log10(abs(gain_tau)), ...
% 'bo', 'Linewidth', lw, 'Markersize', msz);
% stem3(range_true, theta_true, amp_true, ':.m', 'Linewidth', lw); 'True',
% grid on;
% xlim([0 range_max / 2])
% ylim([-40 40])
% xlabel('Range (m)', 'Fontsize', fsz);
% ylabel('Azimuth ($\circ$)', 'Interpreter','latex', 'Fontsize', fsz);
% zlabel('Amplitude (dB)', 'Fontsize', fsz)
% legend('Threshold (NOMP)', 'Amplitude (NOMP)', 'Fontsize', fsz)






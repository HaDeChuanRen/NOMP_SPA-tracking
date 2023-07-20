clc; clear; close all;

% rng(1);
c = 3e8;

Nx = 32;
My = 16;
Lz = 8;

% radar parameter
T_idle = 100e-6;
T_ramp = 60e-6;
T_circle = T_idle + T_ramp;
Fre_start = 77e9;
Slope_fre = 29.982e12;
Fs = 10e6;
Ts = 1 / Fs;
lambda_cw = c / Fre_start;
Rx_interval = lambda_cw / 2;

% the radial veclocity of target belongs to [-velocity_max, velocity_max]
velocity_max = c / (4 * Fre_start * T_circle);

% the radial range of target belongs to [0, range_max]
range_max = (c * Fs) / (2 * Slope_fre);

% the status of targets
K = 4;
sigma_n = 1;
% SNR_all = 5 * rand(1, K) + 1;
SNR_all = 5 * linspace(2, 1, K);
% gain_alltargets = randn(K, 1) + 1j * randn(K, 1);
gain_all = 10 .^ (SNR_all / 20) * sqrt(sigma_n) .* exp(1j * 2 * pi * rand(1, K));

velocity_set = 2 * velocity_max * rand(1, K) -  velocity_max;
range_set = range_max * rand(1, K);
theta_set = pi * rand(1, K) - pi / 2;

omega_x = 4 * pi * (Fre_start * velocity_set + Slope_fre * range_set) * Ts / c;
omega_y = 4 * pi * Fre_start * velocity_set * T_circle / c;
omega_z = 2 * pi * Fre_start * Rx_interval * sin(theta_set) / c;

omega_true = [omega_x; omega_y; omega_z];

tau_max = 2 * range_max / c;
N_max = ceil(tau_max * Fs);
t = 0 : 1 / Fs : T_ramp;

signal_Rx_all = zeros(length(t), My, Lz, K);

LFM_Fun = @(t) exp(- 1j * 2 * pi * (Fre_start * t + 0.5 * Slope_fre * (t .^ 2))) .* (t >= 0 ) .* (t <= T_ramp);
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega);

for m_idx = 1 : My
    tau = 2 * (range_set + velocity_set * (m_idx - 1) * T_circle + t' * velocity_set) / c;
    for k_idx = 1 : K
        t_temp = (t' - tau(:, k_idx));
        signal_Rx_all(:, m_idx, :, k_idx) = gain_all(k_idx) * LFM_Fun(t_temp) * (array_Fun(omega_z(k_idx), Lz)).';
    end
end

signal_Rx = squeeze(sum(signal_Rx_all, 4));

signal_dechirp = bsxfun(@times, signal_Rx, conj(LFM_Fun(t')));

Y_tensor = signal_dechirp((N_max + 1) : (N_max + Nx), :, :);


Noise_tensor = (sigma_n / 2) * (randn(Nx, My, Lz) + 1j * randn(Nx, My, Lz));
Y_receive = Y_tensor + Noise_tensor;

% analyse the received signal
K_max = 10;
guard_n = 2;
guard_m = 2;
guard_l = 2;

training_n = 3;
training_m = 3;
training_l = 3;
guard_training_size = [guard_n, guard_m, guard_l, training_n, training_m, training_l];

P_false_nominal = 1e-3;

gamma_cfar = [1, 1, 1];
gamma_mnomp = [4, 4, 4];
MC = 10000;

alpha_set = alpha_cal_3D_ca(gamma_cfar, P_false_nominal, Nx, My, Lz, guard_training_size, MC);

[omega_list, gain_list, y_residue_tensor] = MNOMP_CFAR_3D_alpha_ca(Y_receive, gamma_cfar, gamma_mnomp, guard_training_size, alpha_set, K_max);
wrapToPi(omega_list - omega_true)
% omega_list
% omega_true









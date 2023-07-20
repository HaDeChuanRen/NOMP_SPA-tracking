clc; clear; close all;

rng(8);
Nx = 256;
My = 64;
Lz = 8;
MC = 300;
P_fa = 0.01;
NM = Nx * My;
N_r = 60;
K_max = 2;
sigma_n = 1;
fa_NOMP = zeros(MC, 1);
fa_noise = zeros(MC, 1);
SNR_kt = 100;
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);
alpha_2D = alpha_PoebyS(P_fa, NM, N_r, Lz);
tau_2D = sigma_n * chi2inv((1 - P_fa) ^ (1 / NM), 2 * Lz) / 2;
gamma_nomp = [1; 1];

gain_true = sqrt(10 ^ (SNR_kt / 10) * sqrt(sigma_n)) .* exp(1j * 2 * ...
        pi * rand());
yten_noise = sqrt(sigma_n / 2) *(randn(Nx, My, Lz) + 1j * randn(Nx, My, ...
    Lz));
omega_true = wrapToPi(2 * pi * rand(3, 1));
yvec_kt = kron(kron(array_Fun(omega_true(3), Lz), ...
    array_Fun(omega_true(2), My)), ...
    array_Fun(omega_true(1), Nx));
yten_kt = gain_true * reshape(yvec_kt, Nx, My, Lz) + yten_noise;
% [omegaxy_t, ghat_t, ~, ~] = NOMP2D_low(yten_kt, alpha_2D, N_r, K_max);

R_s = 20;
Ymat_res = yten_kt;
Yene_res = zeros(R_s, 1);
[omega_k, gain_k, ~, Ymat_res] = DetectNew_2D(Ymat_res, gamma_nomp);
for refione_idx = 1 : R_s
    [Ymat_res, omega_k, gain_k] = RefineOne_2D(Ymat_res, omega_k, gain_k);
    Yene_res(refione_idx) = norm(Ymat_res(:));
end

Yene_noi = norm(yten_noise(:));

figure;
plot(Yene_res);
hold on;
plot(Yene_noi * ones(R_s, 1));


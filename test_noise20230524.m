clc; clear; close all;
rng(5);

Nx = 256;
My = 64;
Lz = 8;
MC = 200;
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

for mc_idx = 1 : MC
    if mc_idx == 3
        1;
    end
    h_waitbar = waitbar((mc_idx - 1) / MC);
    gain_true = sqrt(10 ^ (SNR_kt / 10) * sqrt(sigma_n)) .* exp(1j * 2 * ...
        pi * rand());
    yten_noise = sqrt(sigma_n / 2) *(randn(Nx, My, Lz) + 1j * randn(Nx, My, ...
        Lz));
    omega_true = wrapToPi(2 * pi * rand(3, 1));
    yvec_kt = kron(kron(array_Fun(omega_true(3), Lz), ...
        array_Fun(omega_true(2), My)), ...
        array_Fun(omega_true(1), Nx));
    yten_kt = gain_true * reshape(yvec_kt, Nx, My, Lz) + yten_noise;
    [omegaxy_t, ghat_t, ~, ~] = NOMP2D_low(yten_kt, alpha_2D, N_r, K_max);
    if size(omegaxy_t, 1) > 1
        fa_NOMP(mc_idx) = 1;
    end
    [omega_noise, ~, ~, ~] = NOMP2D_low(yten_kt, alpha_2D, N_r, K_max);
    if size(omega_noise, 1) > 1
        fa_noise(mc_idx) = 1;
    end
end
fa_delta = fa_NOMP - fa_noise;
find(fa_delta ~= 0)
delete(h_waitbar)

% for mc_idx = 1 : MC
%     h_waitbar = waitbar((mc_idx - 1) / MC);
%     y_ten = sqrt(sigma_n / 2) *(randn(Nx, My, Lz) + 1j * randn(Nx, My, Lz));
%     [omegaxy_t, ghat_t, ~, ~] = NOMP2D_low(y_ten, alpha_2D, N_r, K_max);
%     if ~isempty(omegaxy_t)
%         fa_count = fa_count + 1;
%     end
% end




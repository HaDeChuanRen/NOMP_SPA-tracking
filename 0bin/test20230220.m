clc; clear; close all;

K_targets = 1;
omegax = 1;
omegay = 2;
omegaz = pi * sin(pi / 6);
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);

Nx = 256;
My = 64;
Lz = 8;
NML = Nx * My * Lz;
gamma_l = 32;
snr_set = 30;
sigma_n = 1;
gain_true = sqrt(10 .^ (snr_set / 10) * sqrt(sigma_n)) .* exp(1j * 2 * pi * ...
    rand);

yvec_t = kron(kron(array_Fun(omegaz, Lz), array_Fun(omegay, My)), ...
    array_Fun(omegax, Nx));
yvec_t = gain_true * yvec_t + sqrt(1 / 2) * (randn(NML, 1) + 1j * ...
    randn(NML, 1));
yten_spec = reshape(yvec_t, Nx, My, Lz);


NM = Nx * My;
N_r = 60;
P_fa = 1e-2;
alpha_2D = alpha_Poe(P_fa, NM, N_r);
K_max = K_targets * 2;
testvec = zeros(Lz, 1);
for l_idx = 1 : Lz
    ymat_lt = squeeze(yten_spec(:, :, l_idx));
    [omegahat_lt, ghat_lt, ~, ~] = NOMP2D_low(ymat_lt, alpha_2D, N_r, K_max);
    testvec(l_idx) = ghat_lt(1);
end


plot(abs(fft(testvec, Lz * gamma_l)))
[~, z_idx] = max(abs(fft(testvec, Lz * gamma_l)));
omegazhat = wrapToPi((z_idx / (Lz * gamma_l)) * 2 * pi);
asind(omegazhat / pi)





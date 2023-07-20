clc; clear; close all;
format long
% rng(904);

% set the line spectrum vector parameters
D = 3;
K_targets = 4;
Nx = 32;
My = 16;
Lz = 8;
NML = Nx * My * Lz;

% set the signal parameters

gain_set = linspace(1.5, 0.5, K_targets);

omega_true = 2 * pi * rand(D, K_targets);
% omega_true = [1, 4; 2, 5; 3, 6];

phase_set = exp(1j * 2 * pi * rand(K_targets, 1));
amp_true = gain_set' .* phase_set;

% create the receive signal with noise
y_observe = zeros(NML, 1);
for k_idx = 1 : K_targets
    % get the line spectrum vector of each dimension
    aNx_vec_k = exp(1j * (0 : Nx - 1)' * omega_true(1, k_idx)) / sqrt(Nx);
    aMy_vec_k = exp(1j * (0 : My - 1)' * omega_true(2, k_idx)) / sqrt(My);
    aLz_vec_k = exp(1j * (0 : Lz - 1)' * omega_true(3, k_idx)) / sqrt(Lz);

    atom_vec_k = kron(kron(aLz_vec_k, aMy_vec_k), aNx_vec_k);

    y_observe = y_observe + amp_true(k_idx) * atom_vec_k;
end

% set the noise parameter and add noise to the received signal
sigma_n = 1e-6;
Noise = sqrt(sigma_n / 2) * (randn(NML, 1) + 1j * randn(NML, 1));
y_observe = y_observe + Noise;

y_tensor = reshape(y_observe, Nx, My, Lz);
% training_height = 3;
% training_width = 2;
% training_page = 1;
P_false_nominal = 1e-6;
K_max = 10;
%[omega_list, gain_list, y_residue_tensor] = MNOMP_CFAR_3D(y_tensor, K_max, training_height, training_width, training_page, P_false_nominal);

training_set = [4, 3, 2];
guard_set = [3, 2, 1];

tic;
[omega_list, gain_list, y_residue_tensor] = MNOMP_CA_3D(y_tensor, K_max, training_set, guard_set, P_false_nominal);
toc;

omega_list - omega_true


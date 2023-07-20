clc;clear;close all;
format long
%rng(4329);

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



% the size vector after oversampling
oversample_rate = 4;
Oversam_vec = [oversample_rate * Nx; oversample_rate * My; oversample_rate * Lz];

y_tensor = reshape(y_observe, Nx, My, Lz);


% target detection and refinement
y_residue_tensor = y_tensor;
K_max = 10;
omega_est = [];
ghat = [];
K_est = 0;

% guard_width = 2;
% guard_height = 2;
% training_width = 3;
% training_height = 3;
% num_training_cell = (2 * training_width + 1) * (2 * training_height + 1) - (2 * guard_width + 1) * (2 * guard_height + 1);

training_height = 4;
training_width = 3;
training_page = 2;
guard_height = 3;
guard_width = 2;
guard_page = 1;

num_training_cell = ...
(2 * training_width + 1) * (2 * training_height + 1) * (2 * training_page + 1) - ...
(2 * guard_width + 1) * (2 * guard_height + 1) * (2 * guard_page + 1);

P_false_nominal = 1e-9;

while K_est < K_max
    [omega_k, tensor_idx] = detectnew_3D(y_residue_tensor, oversample_rate);
    summit_idx = ceil(tensor_idx / 4);

    % [row, col, page] = ind2sub(Oversam_vec,idx_max);
    % Y_gamma = y_residue_tensor(:, :, ceil(page / 4));

    % [sigma_hat, Y_gamma_max] = sigma_cfar_est(Y_gamma, guard_width, guard_height, training_width, training_height);
    % noise_tau = sigma_hat * num_training_cell * (P_false_nominal ^ (- 1 / num_training_cell) - 1);

    [~, omega_k, ghat_k] = refineone_3D(y_residue_tensor, omega_k);
    y_residue_fftn = fftn(y_residue_tensor, [Nx, My, Lz]) / sqrt(NML);


    y_residue_ext = repmat(y_residue_fftn, 3, 3, 3);
    Y_training_win = y_residue_ext(summit_idx(1) + Nx + (- training_height : training_height),...
    summit_idx(2) + My + (- training_width : training_width),...
    summit_idx(3) + Lz + (- training_page : training_page));
    Y_guard_win = y_residue_ext(summit_idx(1) + Nx + (- guard_height : guard_height),...
    summit_idx(2) + My + (- guard_width : guard_width),...
    summit_idx(3) + Lz + (- guard_page : guard_page));

    sigma_hat = (sum(abs(Y_training_win) .^ 2, 'all') - sum(abs(Y_guard_win) .^ 2, 'all')) / num_training_cell;
    tau = sigma_hat * num_training_cell * (P_false_nominal ^ (- 1 / num_training_cell) - 1);
    if(abs(ghat_k)^2 < tau)
        break;
    end

    K_est = K_est + 1;
    omega_est = [omega_est, omega_k];
    ghat = [ghat; ghat_k];

    [omega_est, y_residue_tensor, ghat] = refineall_3D(y_tensor, omega_est, ghat);
end

omega_est - omega_true




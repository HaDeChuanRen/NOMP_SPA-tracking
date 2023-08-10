function [ground_all, yten_spec] = grdsig_prod(sample_time, Xori_mat, Amat, ...
    sigma_w, sigma_n, rmax, vmax, Nx, My, Lz, SNR_kt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Last update: 2023/7/28
% 2023/7/28 update details:
% 1. (Todo) add comments and description on the function.
% 2. (Todo) test and debug the function

[state_dim, K_targets] = size(Xori_mat);
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);

Xtrue_kt = Amat * Xori_mat + sqrt(sigma_w) .* randn(state_dim, 1);
ground_all = zeros(state_dim, sample_time, K_targets);
ground_all(1, 1, :) = Xtrue_kt(1, :);
ground_all(2, 1, :) = Xtrue_kt(2, :);
ground_all(3, 1, :) = Xtrue_kt(3, :);
ground_all(4, 1, :) = Xtrue_kt(4, :);

for t_time = 2 : sample_time
    Xtrue_kt = Amat * Xtrue_kt + sqrt(sigma_w) .* randn(state_dim, 1);
    ground_all(1, t_time, :) = Xtrue_kt(1, :);
    ground_all(2, t_time, :) = Xtrue_kt(2, :);
    ground_all(3, t_time, :) = Xtrue_kt(3, :);
    ground_all(4, t_time, :) = Xtrue_kt(4, :);
end


NML = Nx * My * Lz;
% initialize real r, v, theta parameters of each instant
rmat_true = zeros(K_targets, sample_time);
vmat_true = zeros(K_targets, sample_time);
thetamat_true = zeros(K_targets, sample_time);

% initialize real omega parameters of each dimension
omegaxmat_true = zeros(K_targets, sample_time);
omegaymat_true = zeros(K_targets, sample_time);
omegazmat_true = zeros(K_targets, sample_time);


yten_spec = zeros(Nx, My, Lz, sample_time);
for t_time = 1 : sample_time
    % obtain the real measurements in the instant t
    rmat_true(:, t_time) = sqrt(ground_all(1, t_time, :) .^ 2 + ...
        ground_all(3, t_time, :) .^ 2);
    thetamat_true(:, t_time) = atan(ground_all(1, t_time, :) ./ ...
        ground_all(3, t_time, :));
    vmat_true(:, t_time) = squeeze(ground_all(2, t_time, :)) .* ...
        cos(thetamat_true(:, t_time)) + squeeze(ground_all(4, t_time, :)) .* ...
        sin(thetamat_true(:, t_time));

    % calculate the corresponding angular frequency
    % $\omega \in (-\pi, \pi]$
    omegaxmat_true(:, t_time) = wrapToPi(rmat_true(:, t_time) / rmax * 2 * pi);
    omegaymat_true(:, t_time) = vmat_true(:, t_time) / vmax * pi;
    omegazmat_true(:, t_time) = pi * sin(thetamat_true(:, t_time));

    % create the line spectrum tensors in the instant t
    yvec_t = zeros(NML, 1);
    for k_idx = 1 : K_targets
        gain_kt = sqrt(10 ^ (SNR_kt / 10) * sqrt(sigma_n)) .* exp(1j * 2 * ...
            pi * rand());
        yvec_kt = kron(kron(array_Fun(omegazmat_true(k_idx, t_time), Lz), ...
            array_Fun(omegaymat_true(k_idx, t_time), My)), ...
            array_Fun(omegaxmat_true(k_idx, t_time), Nx));
        yvec_t = yvec_t + gain_kt * yvec_kt;
    end
    yvec_t = yvec_t + sqrt(sigma_n / 2) * (randn(NML, 1) + 1j * randn(NML, 1));
    yten_spec(:, :, :, t_time) = reshape(yvec_t, Nx, My, Lz);
end


end
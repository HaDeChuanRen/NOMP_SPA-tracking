clc; clear; close all;
format long
rng(905);

% set the line spectrum vector parameters
D = 3;
% K_targets = 4;
Nx = 32;
My = 16;
Lz = 8;
NML = Nx * My * Lz;

% function parameter
% training_height = 6;
% training_width = 4;
% training_page = 1;
training_set = [4, 3, 2];
guard_set = [3, 2, 1];
K_max = 100;

% noise parameter
sigma_n = 1e-3;

% simulation parameter
simulate_times = 10;

P_fa_set = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6];
P_fa_get = zeros(1, 12);
for k_idx = 1 : 12
    handle_waitbar = waitbar(k_idx / 12);
    P_false_nominal = P_fa_set(k_idx);
    K_collect = zeros(simulate_times, 1);

    for t_sim = 1 : simulate_times
        Noise = sqrt(sigma_n / 2) * (randn(NML, 1) + 1j * randn(NML, 1));
        y_observe = Noise;
        y_tensor = reshape(y_observe, Nx, My, Lz);

        %[omega_list, gain_list, y_residue_tensor] = MNOMP_CFAR_3D(y_tensor, K_max, training_height, training_width, training_page, P_false_nominal);
        [omega_list, gain_list, y_residue_tensor] = MNOMP_CA_3D(y_tensor, K_max, training_set, guard_set, P_false_nominal);

        K_est = length(gain_list);
        K_collect(t_sim) = K_est;
    end
    false_times = sum(K_collect, 'all');
    P_fa_get(k_idx) = false_times / (simulate_times * NML);
end


plot(P_fa_set, P_fa_get);
delete(handle_waitbar);










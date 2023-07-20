function [omega_list, gain_list, y_residue_tensor] = MNOMP_CFAR_3D(y_tensor, K_max, training_height, training_width, training_page, P_false_nominal)
%   this funciton is used to analyse 3D linespectrum with NOMP and CA-CFAR
%
%   author: Menghuai Xu
%   email: 22134019@zju.edu.cn
%   date: 2021.8.10

%   Input parameter:
%   y_tensor: the three dimension receiving signal been observed              #Nx * My * Lz
%   K_max: the maximum targets in a received tensor

    omega_est = [];
    ghat = [];
    K_est = 0;
    y_residue_tensor = y_tensor;
    [Nx, My, Lz] = size(y_tensor);
    num_training_cell = (2 * training_width + 1) * (2 * training_height + 1) * (2 * training_page + 1);
    oversample_rate = 4;

    while K_est < K_max
        [omega_k, tensor_idx] = detectnew_3D(y_residue_tensor, oversample_rate);
        summit_idx = ceil(tensor_idx / 4);

        [y_residue_tensor, omega_k, ghat_k] = refineone_3D(y_residue_tensor, omega_k);
        y_residue_ext = repmat(y_residue_tensor, 3, 3, 3);
        Y_training_win = y_residue_ext(summit_idx(1) + Nx + (- training_height : training_height),...
        summit_idx(2) + My + (- training_width : training_width),...
        summit_idx(3) + Lz + (- training_page : training_page));

        sigma_hat = sum(abs(Y_training_win) .^ 2, 'all') / num_training_cell;
        tau = sigma_hat * num_training_cell * (P_false_nominal ^ (- 1 / num_training_cell) - 1);
        if(abs(ghat_k) ^ 2 < tau)
            break;
        end

        K_est = K_est + 1;
        omega_est = [omega_est, omega_k];
        ghat = [ghat; ghat_k];

        [omega_est, y_residue_tensor, ghat] = refineall_3D(y_tensor, omega_est, ghat);
    end

    omega_list = omega_est;
    gain_list = ghat;


end


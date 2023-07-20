function [omega_list, gain_list, y_residue_tensor] = MNOMP_CA_3D(y_tensor, K_max, training_set, guard_set, P_false_nominal)
%   this function is used to
%   此处显示详细说明

training_height = training_set(1);
training_width = training_set(2);
training_page = training_set(3);
guard_height = guard_set(1);
guard_width = guard_set(2);
guard_page = guard_set(3);

num_training_cell = ...
(2 * training_width + 1) * (2 * training_height + 1) * (2 * training_page + 1) - ...
(2 * guard_width + 1) * (2 * guard_height + 1) * (2 * guard_page + 1);

omega_est = [];
ghat = [];
K_est = 0;
y_residue_tensor = y_tensor;
[Nx, My, Lz] = size(y_tensor);
NML = Nx * My * Lz;
oversample_rate = 4;

while K_est < K_max
    [omega_k, tensor_idx] = detectnew_3D(y_residue_tensor, oversample_rate);
    summit_idx = ceil(tensor_idx / 4);

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
    tau_hat = sigma_hat * num_training_cell * (P_false_nominal ^ (- 1 / num_training_cell) - 1);
    % tau_hat = - sigma_hat * log(P_false_nominal);
    if(abs(ghat_k) ^ 2 < tau_hat)
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


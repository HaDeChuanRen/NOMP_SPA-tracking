function [omega_k, ghat_k, prob_ind, Y_r] = DetectNew_2D(y_matrix, gamma_oversampling)

    [Nx, My, T] = size(y_matrix);
    NM_num = Nx * My;
    ant_idx_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
    ant_idx_My = (0 : (My - 1))' - (My - 1) / 2;

    array_Fun = @(omega, N) exp(1j * ((0 : (N - 1))' - (N - 1) / 2) * omega) / sqrt(N);
%     y_residue_matrix = y_matrix;

    Nx_ext = Nx * gamma_oversampling(1);
    My_ext = My * gamma_oversampling(2);
    Oversam_vec = [Nx_ext, My_ext];    % the size vector after oversampling

    % 2D-fft and find the maximum index
%     GLRT_est = sum(abs(fft2(Y, R_x, R_y)).^2,3);
    Y_spec = fft2(y_matrix, Nx_ext, My_ext) / sqrt(NM_num); % Nx_ext*My_ext*T
%     CoarseOmega_x = 2 * pi * (0 : (Nx_ext - 1))' / Nx_ext;
%     CoarseOmega_y = 2 * pi * (0 : (My_ext - 1))' / My_ext;
%     CoarseOmega_vec = kron(CoarseOmega_y * ant_idx_My(1), CoarseOmega_x * ant_idx_Nx(1));
%     CoarseOmega_mat = reshape(CoarseOmega_vec, [Nx_ext, My_ext]);
%     CoarseOmega_mat = CoarseOmega_x * ant_idx_Nx(1)
%     Y_spec = bsxfun(@times, Y_spec, exp(-1j * CoarseOmega_mat));
    a2D_fft_abs2 = sum(abs(Y_spec).^2,3);
    prob_ind = a2D_fft_abs2(1 : gamma_oversampling(1) : end, ...
        1 : gamma_oversampling(2) : end);
    % GLRT = max(a2D_fft_abs2_sub(:));
    [~, idx_max] = max(a2D_fft_abs2, [], 'all', 'linear') ;
    [row, col] = ind2sub(Oversam_vec, idx_max);
    ghat_k1 = squeeze(Y_spec(row, col,:));
    ghat_k1 = ghat_k1.';
    matrix_idx = [row, col];

    % get the angular frequency
    omega_k = (matrix_idx - 1) ./ Oversam_vec * 2 * pi;
    ghat_k = ghat_k1 * exp(-1j * (ant_idx_Nx(1) * omega_k(1) + ant_idx_My(1) * omega_k(2)));
    %% update Y_r
    Y_r = zeros(Nx, My, T);
    for t_idx = 1:T
        Y_r(:, :, t_idx) = y_matrix(:, :, t_idx) - ghat_k(t_idx)*array_Fun(omega_k(1), Nx)*(array_Fun(omega_k(2), My).');
    end

%     atom_vec_est = kron(array_Fun(omega_k(2), My), array_Fun(omega_k(1), Nx));
%     atom_mat_est = reshape(atom_vec_est, Nx, My);
%     y_residue_matrix = y_matrix - ghat_k * atom_mat_est;
end
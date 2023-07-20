function [omega_k, ghat_ksym] = Det2D_spec(Y_spec, g_osr)

    [Nx_ext, My_ext, ~] = size(Y_spec);
    Nx = Nx_ext / g_osr(1);
    My = My_ext / g_osr(2);

    ant_idx_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
    ant_idx_My = (0 : (My - 1))' - (My - 1) / 2;
    Oversam_vec = [Nx_ext, My_ext];    % the size vector after oversampling

    % 2D-fft and find the maximum index
    a2D_fft_abs2 = sum(abs(Y_spec).^2, 3);
    [~, idx_max] = max(a2D_fft_abs2, [], 'all', 'linear') ;
    [row, col] = ind2sub(Oversam_vec, idx_max);
    ghat_k1 = squeeze(Y_spec(row, col, :));
    ghat_k1 = ghat_k1.';
    matrix_idx = [row, col];

    % get the angular frequency
    omega_k = (matrix_idx - 1) ./ Oversam_vec * 2 * pi;
    ghat_ksym = ghat_k1 * exp(-1j * (ant_idx_Nx(1) * omega_k(1) + ...
        ant_idx_My(1) * omega_k(2)));

    % update Y_r
    % Y_r = zeros(Nx, My, T);
    % for t_idx = 1:T
    %     Y_r(:, :, t_idx) = y_matrix(:, :, t_idx) - ghat_k(t_idx)*array_Fun(omega_k(1), Nx)*(array_Fun(omega_k(2), My).');
    % end
end
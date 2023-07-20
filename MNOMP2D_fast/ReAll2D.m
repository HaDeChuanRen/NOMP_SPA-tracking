function [omega_est, ghat, Ymat_r] = ReAll2D(Ymat_r, omega_est, ghat, R_s, R_c)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    K_hat = size(omega_est, 1);          % the number of detected frequencies
    for c_idx = 1 : R_c
        for k_idx = 1 : K_hat
            omega_k = omega_est(k_idx, :);
            ghat_k = ghat(k_idx,:);
            for s_idx = 1 : R_s
                [Ymat_r, omeganew_k, ghatnew_k] = ReOne2D(Ymat_r, omega_k, ghat_k);
            end
            omega_est(k_idx, :) = omeganew_k;
            ghat(k_idx, :) = ghatnew_k;
        end
    end
end


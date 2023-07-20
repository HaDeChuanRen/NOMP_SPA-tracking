function [Ymat_r, omega_hat, ghat] = ReOne2D(Ymat, omega_est, ghat)
%UNTITLED2 此处显示有关此函数的摘要
%   date: 2023/5/30
    Ymat_r = Ymat;
    obj_old = norm(Ymat_r(:));
    [Nx, My, T] = size(Ymat);
    ant_idx_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
    ant_idx_My = (0 : (My - 1))' - (My - 1) / 2;

    xhat_vec = exp((1j * ant_idx_Nx * omega_est(1))) / sqrt(Nx);
    yhat_vec = exp((1j * ant_idx_My * omega_est(2))) / sqrt(My);

    % get the gradient of vector a
    x1diff_vec = 1j * ant_idx_Nx .* xhat_vec;
    y1diff_vec = 1j * ant_idx_My .* yhat_vec;
    ax1diff_vec = kron(yhat_vec, x1diff_vec);
    ay1diff_vec = kron(y1diff_vec, xhat_vec);

    da_over_domega = [ax1diff_vec, ay1diff_vec];
    Hessian_S = zeros(2, 2);
    grad_S = zeros(2, 1);
    Ymat_add = zeros(Nx, My, T);
    for t = 1 : T
        Ymat_add(:, :, t) = Ymat(:, :, t) + ghat(t) * xhat_vec * (yhat_vec.');
        Ymat_addt = Ymat(:, :, t);
        yvec_temp = Ymat_addt(:);
        grad_S = grad_S - 2 * real(conj(ghat(t)) * da_over_domega' * yvec_temp);
        x2diff_vec = 1j * ant_idx_Nx .* x1diff_vec;
        y2diff_vec = 1j * ant_idx_My .* y1diff_vec;
        ax2diff_vec = kron(yhat_vec, x2diff_vec);
        ay2diff_vec = kron(y2diff_vec, xhat_vec);
        axy2diff_vec = kron(y1diff_vec, x1diff_vec);
        Hessian_S(1, 1) = Hessian_S(1, 1)- 2 * real(ghat(t)' * ax2diff_vec' *...
            yvec_temp);
        Hessian_S(2, 2) = Hessian_S(2, 2)- 2 * real(ghat(t)' * ay2diff_vec' *...
            yvec_temp);
        Hessian_S(1, 2) = Hessian_S(1, 2)- 2 * real(ghat(t)' * axy2diff_vec' ...
            * yvec_temp);
    end
    Hessian_S(2, 1) = Hessian_S(1, 2);

    Step = Hessian_S \ grad_S;
    if max(abs(Step)) < 2 * pi / max(Nx, My) / 4
        omega_hat_new = omega_est - Step';
        % sprintf('Newton is effective')
    else
        omega_hat_new = omega_est - sign(grad_S') * 2 * pi / max(Nx, My) / 100;
        % sprintf('Newton is ineffective and Grad is used')
    end

    xhat_vec_new = exp((1j * ant_idx_Nx * omega_hat_new(1))) / sqrt(Nx);
    yhat_vec_new = exp((1j * ant_idx_My * omega_hat_new(2))) / sqrt(My);
    ghat_new = zeros(1, T);
    for t = 1 : T
        ghat_new(t) = xhat_vec_new' * Ymat_add(:, :, t) * conj(yhat_vec_new);
        Ymat_r(:, :, t) = Ymat_add(:, :, t) - ghat_new(t) * xhat_vec_new * ...
            (yhat_vec_new.');
    end

    obj_new = norm(Ymat_r(:));
    if obj_new < obj_old
        omega_hat = omega_hat_new;
        ghat = ghat_new;
        % sprintf('New is used')
    else
        Ymat_r = Ymat;
        omega_hat = omega_est;
    end
end


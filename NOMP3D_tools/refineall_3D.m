function [omega_hat, y_residue_tensor, gain_est] = refineall_3D(y_tensor, omega_est, ghat)
%   this funciton is used to refine one tartgets' angular frequencies and amplitudes estimation,
%   meanwhile it can also get the residue of received tensor by removing the affect of targets
%   which have been detected
%
%   Input parameter:
%   y_tensor: the three dimension receiving signal been observed              #Nx * My * Lz
%   omega_est: the coarse estimation of angular frequency
%   ghat: the coarse estimation of amplitudes
%
%   author: Menghuai Xu
%   email: 22134019@zju.edu.cn
%   date: 2021.11.7

    [Nx, My, Lz] = size(y_tensor);
    NML = Nx * My * Lz;
    [D_dimension, K_est] = size(omega_est);
    omega_hat = omega_est;
    gain_est = ghat;
    y_vec = y_tensor(:);

    xhat_mat = exp(1j * (0 : Nx - 1)' * omega_hat(1, :)) / sqrt(Nx);
    yhat_mat = exp(1j * (0 : My - 1)' * omega_hat(2, :)) / sqrt(My);
    zhat_mat = exp(1j * (0 : Lz - 1)' * omega_hat(3, :)) / sqrt(Lz);

    atom_mat_est = zeros(NML, K_est);
    for kidx = 1:K_est
        atom_mat_est(:,kidx) = kron(kron(zhat_mat(:, kidx), yhat_mat(:, kidx)), xhat_mat(:, kidx));
    end

    y_residue = y_vec - atom_mat_est * gain_est;
    % obj_old = norm(y_vec-atom_mat_est*gain_est)^2;

    grad_all_S = zeros(3 * K_est, 1);
    S_allhessian = zeros(3 * K_est, 3 * K_est);
    omega_hat_vec = omega_hat(:);
    grad_iter = 50;

    continue_flag = 1;


    grad_idx = 1;
    while (grad_idx <= grad_iter) && (continue_flag == 1)
        for k_idx = 1:K_est
            xhat_vec = xhat_mat(:, k_idx);
            yhat_vec = yhat_mat(:, k_idx);
            zhat_vec = zhat_mat(:, k_idx);
            x1diff_vec = 1j * (0 : Nx - 1)' .* xhat_vec;
            y1diff_vec = 1j * (0 : My - 1)' .* yhat_vec;
            z1diff_vec = 1j * (0 : Lz - 1)' .* zhat_vec;

            % get the gradient of vector a
            ax1diff_vec = kron(kron(zhat_vec, yhat_vec), x1diff_vec);
            ay1diff_vec = kron(kron(zhat_vec, y1diff_vec), xhat_vec);
            az1diff_vec = kron(kron(z1diff_vec, yhat_vec), xhat_vec);

            da_over_domega = [ax1diff_vec, ay1diff_vec, az1diff_vec];

            % grad_S = - 2 * real(gain_est * y_residue' * da_over_domega);
            grad_sub_S = - 2 * real(gain_est(k_idx)' * da_over_domega' * y_residue);
            grad_all_S((k_idx - 1) * 3 + 1 : k_idx * 3) = grad_sub_S;
        end

        for k_idx = 1 : K_est
            for k_jdx = 1 : K_est
                if(k_idx == k_jdx)
                    xhat_vec = exp(1j * (0 : Nx - 1)' * omega_hat(1, k_idx)) / sqrt(Nx);
                    yhat_vec = exp(1j * (0 : My - 1)' * omega_hat(2, k_idx)) / sqrt(My);
                    zhat_vec = exp(1j * (0 : Lz - 1)' * omega_hat(3, k_idx)) / sqrt(Lz);

                    x1diff_vec = 1j * (0 : Nx - 1)' .* xhat_vec;
                    y1diff_vec = 1j * (0 : My - 1)' .* yhat_vec;
                    z1diff_vec = 1j * (0 : Lz - 1)' .* zhat_vec;

                    ax1diff_vec = kron(kron(zhat_vec, yhat_vec), x1diff_vec);
                    ay1diff_vec = kron(kron(zhat_vec, y1diff_vec), xhat_vec);
                    az1diff_vec = kron(kron(z1diff_vec, yhat_vec), xhat_vec);

                    x2diff_vec = 1j * (0 : Nx - 1)' .* x1diff_vec;
                    y2diff_vec = 1j * (0 : My - 1)' .* y1diff_vec;
                    z2diff_vec = 1j * (0 : Lz - 1)' .* z1diff_vec;

                    ax2diff_vec = kron(kron(zhat_vec, yhat_vec), x2diff_vec);
                    ay2diff_vec = kron(kron(zhat_vec, y2diff_vec), xhat_vec);
                    az2diff_vec = kron(kron(z2diff_vec, yhat_vec), xhat_vec);

                    axy2diff_vec = kron(kron(zhat_vec, y1diff_vec), x1diff_vec);
                    ayz2diff_vec = kron(kron(z1diff_vec, y1diff_vec), xhat_vec);
                    axz2diff_vec = kron(kron(z1diff_vec, yhat_vec), x1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 1) = ...
                    - 2 * real(gain_est(k_idx)' * ax2diff_vec' * y_residue) + 2 * abs(gain_est(k_idx)) ^ 2 * norm(ax1diff_vec) ^ 2;
                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 2) = ...
                    - 2 * real(gain_est(k_idx)' * ay2diff_vec' * y_residue) + 2 * abs(gain_est(k_idx)) ^ 2 * norm(ay1diff_vec) ^ 2;
                    S_allhessian((k_idx - 1) * 3 + 3, (k_jdx - 1) * 3 + 3) = ...
                    - 2 * real(gain_est(k_idx)' * az2diff_vec' * y_residue) + 2 * abs(gain_est(k_idx)) ^ 2 * norm(az1diff_vec) ^ 2;

                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 2) = ...
                    - 2 * real(gain_est(k_idx)' * axy2diff_vec' * y_residue) + 2 * real(abs(gain_est(k_idx))^2 * ax1diff_vec' * ay1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 1) = ...
                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 2);

                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 3) = ...
                    - 2 * real(gain_est(k_idx)' * axz2diff_vec' * y_residue) + 2 * real(abs(gain_est(k_idx))^2 * ax1diff_vec' * az1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 3, (k_jdx - 1) * 3 + 1) = ...
                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 3);

                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 3) = ...
                    - 2 * real(gain_est(k_idx)' * ayz2diff_vec' * y_residue) + 2 * real(abs(gain_est(k_idx))^2 * ay1diff_vec' * az1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 3, (k_jdx - 1) * 3 + 2) = ...
                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 3);
                else

                    xihat_vec = exp(1j * (0 : Nx - 1)' * omega_hat(1, k_idx)) / sqrt(Nx);
                    yihat_vec = exp(1j * (0 : My - 1)' * omega_hat(2, k_idx)) / sqrt(My);
                    zihat_vec = exp(1j * (0 : Lz - 1)' * omega_hat(3, k_idx)) / sqrt(Lz);

                    xi1diff_vec = 1j * (0 : Nx - 1)' .* xihat_vec;
                    yi1diff_vec = 1j * (0 : My - 1)' .* yihat_vec;
                    zi1diff_vec = 1j * (0 : Lz - 1)' .* zihat_vec;

                    axi1diff_vec = kron(kron(zihat_vec, yihat_vec), xi1diff_vec);
                    ayi1diff_vec = kron(kron(zihat_vec, yi1diff_vec), xihat_vec);
                    azi1diff_vec = kron(kron(zi1diff_vec, yihat_vec), xihat_vec);

                    xjhat_vec = exp(1j * (0 : Nx - 1)' * omega_hat(1, k_jdx)) / sqrt(Nx);
                    yjhat_vec = exp(1j * (0 : My - 1)' * omega_hat(2, k_jdx)) / sqrt(My);
                    zjhat_vec = exp(1j * (0 : Lz - 1)' * omega_hat(3, k_jdx)) / sqrt(Lz);

                    xj1diff_vec = 1j * (0 : Nx - 1)' .* xjhat_vec;
                    yj1diff_vec = 1j * (0 : My - 1)' .* yjhat_vec;
                    zj1diff_vec = 1j * (0 : Lz - 1)' .* zjhat_vec;

                    axj1diff_vec = kron(kron(zjhat_vec, yjhat_vec), xj1diff_vec);
                    ayj1diff_vec = kron(kron(zjhat_vec, yj1diff_vec), xjhat_vec);
                    azj1diff_vec = kron(kron(zj1diff_vec, yjhat_vec), xjhat_vec);

                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 1) = ...
                    2 * real(gain_est(k_idx)' * gain_est(k_jdx) * axi1diff_vec' * axj1diff_vec);
                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 2) = ...
                    2 * real(gain_est(k_idx)' * gain_est(k_jdx) * ayi1diff_vec' * ayj1diff_vec);
                    S_allhessian((k_idx - 1) * 3 + 3, (k_jdx - 1) * 3 + 3) = ...
                    2 * real(gain_est(k_idx)' * gain_est(k_jdx) * azi1diff_vec' * azj1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 2) = ...
                    2 * real(gain_est(k_idx)' * gain_est(k_jdx) * axi1diff_vec' * ayj1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 1) = ...
                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 2);

                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 3) = ...
                    2 * real(gain_est(k_idx)' * gain_est(k_jdx) * axi1diff_vec' * azj1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 3, (k_jdx - 1) * 3 + 1) = ...
                    S_allhessian((k_idx - 1) * 3 + 1, (k_jdx - 1) * 3 + 3);

                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 3) = ...
                    2 * real(gain_est(k_idx)' * (k_jdx) * azi1diff_vec' * axj1diff_vec);

                    S_allhessian((k_idx - 1) * 3 + 3, (k_jdx - 1) * 3 + 2) = ...
                    S_allhessian((k_idx - 1) * 3 + 2, (k_jdx - 1) * 3 + 3);

                end
            end
        end


        omega_hat_vec = omega_hat_vec - S_allhessian \ grad_all_S;

        % if(sum(S_allhessian \ grad_all_S) < 1e-12)
        %     continue_flag = 0;
        % end


        omega_hat = reshape(omega_hat_vec, D_dimension, K_est);
        % update the dic
        xhat_mat = exp(1j * (0 : Nx - 1)' * omega_hat(1,:)) / sqrt(Nx);
        yhat_mat = exp(1j * (0 : My - 1)' * omega_hat(2,:)) / sqrt(My);
        zhat_mat = exp(1j * (0 : Lz - 1)' * omega_hat(3,:)) / sqrt(Lz);
        for kidx = 1 : K_est
            atom_mat_est(:, kidx) = kron(kron(zhat_mat(:, kidx), yhat_mat(:, kidx)), xhat_mat(:, kidx));
        end

        % update gain_est
        gain_est = atom_mat_est \ y_vec;
        y_residue = y_vec - atom_mat_est * gain_est;

        xhat_mat = exp(1j * (0 : Nx - 1)' * omega_hat(1,:)) / sqrt(Nx);
        yhat_mat = exp(1j * (0 : My - 1)' * omega_hat(2,:)) / sqrt(My);
        zhat_mat = exp(1j * (0 : Lz - 1)' * omega_hat(3,:)) / sqrt(Lz);

        grad_idx = grad_idx + 1;
    end
    y_residue_tensor = reshape(y_residue, Nx, My, Lz);
end


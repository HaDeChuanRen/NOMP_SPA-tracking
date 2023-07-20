function [y_residue_tensor, omega_hat, ghat] = refineone_3D(y_tensor, omega_est)
%   this funciton is used to refine one tartget's angular frequency and amplitude estimation,
%   meanwhile it can also get the residue of received tensor by removing the affect of target
%   which have been detected
%
%   Input parameter:
%   y_tensor: the three dimension receiving signal been observed              #Nx * My * Lz
%   omega_est: the coarse estimation of angular frequency
%
%   author: Menghuai Xu
%   email: 22134019@zju.edu.cn
%   date: 2021.11.7

    Tmax = 20;
    obj = zeros(Tmax, 1);

    [Nx, My, Lz] = size(y_tensor);
    y_vector = y_tensor(:);
    omega_hat = omega_est;
    flag = 1;

    for t = 1 : Tmax
        % Update ghat
        xhat_vec = exp((1j * (0 : Nx - 1)' * omega_hat(1))) / sqrt(Nx);
        yhat_vec = exp((1j * (0 : My - 1)' * omega_hat(2))) / sqrt(My);
        zhat_vec = exp((1j * (0 : Lz - 1)' * omega_hat(3))) / sqrt(Lz);
        atom_vec_est = kron(kron(zhat_vec, yhat_vec), xhat_vec);
        ghat = atom_vec_est' * y_vector;

        % constrct dictionary
        A_idc = atom_vec_est;

        % calculate the residue and cost function
        y_residue = y_vector - atom_vec_est * ghat;
        obj(t) = norm(y_residue)^2;

        % calculate the gradient of cost function
        % get the first-order derivative of vector x, y and z
        x1diff_vec = 1j * (0 : Nx - 1)' .* xhat_vec;
        y1diff_vec = 1j * (0 : My - 1)' .* yhat_vec;
        z1diff_vec = 1j * (0 : Lz - 1)' .* zhat_vec;

        % get the gradient of vector a
        ax1diff_vec = kron(kron(zhat_vec, yhat_vec), x1diff_vec);
        ay1diff_vec = kron(kron(zhat_vec, y1diff_vec), xhat_vec);
        az1diff_vec = kron(kron(z1diff_vec, yhat_vec), xhat_vec);

        da_over_domega = [ax1diff_vec, ay1diff_vec, az1diff_vec];

        % grad_S = - 2 * real(ghat * y_residue' * da_over_domega);
        grad_S = - 2 * real(ghat' * da_over_domega' * y_residue);

        % calculate the Hession matrix of cost function
        % get the second-order derivative of vector x, y and z
        x2diff_vec = 1j * (0 : Nx - 1)' .* x1diff_vec;
        y2diff_vec = 1j * (0 : My - 1)' .* y1diff_vec;
        z2diff_vec = 1j * (0 : Lz - 1)' .* z1diff_vec;

        ax2diff_vec = kron(kron(zhat_vec, yhat_vec), x2diff_vec);
        ay2diff_vec = kron(kron(zhat_vec, y2diff_vec), xhat_vec);
        az2diff_vec = kron(kron(z2diff_vec, yhat_vec), xhat_vec);

        axy2diff_vec = kron(kron(zhat_vec, y1diff_vec), x1diff_vec);
        ayz2diff_vec = kron(kron(z1diff_vec, y1diff_vec), xhat_vec);
        axz2diff_vec = kron(kron(z1diff_vec, yhat_vec), x1diff_vec);

        Hessian_S = zeros(3, 3);

        Hessian_S(1, 1) = - 2 * real(ghat' * ax2diff_vec' * y_residue) + 2 * abs(ghat)^2 * norm(ax1diff_vec) ^ 2;
        Hessian_S(2, 2) = - 2 * real(ghat' * ay2diff_vec' * y_residue) + 2 * abs(ghat)^2 * norm(ay1diff_vec) ^ 2;
        Hessian_S(3, 3) = - 2 * real(ghat' * az2diff_vec' * y_residue) + 2 * abs(ghat)^2 * norm(az1diff_vec) ^ 2;

        Hessian_S(1, 2) = - 2 * real(ghat' * axy2diff_vec' * y_residue) + 2 * real(abs(ghat)^2 * ax1diff_vec' * ay1diff_vec);
        Hessian_S(2, 1) = Hessian_S(1, 2);
        Hessian_S(1, 3) = - 2 * real(ghat' * axz2diff_vec' * y_residue) + 2 * real(abs(ghat)^2 * ax1diff_vec' * az1diff_vec);
        Hessian_S(3, 1) = Hessian_S(1, 3);
        Hessian_S(2, 3) = - 2 * real(ghat' * ayz2diff_vec' * y_residue) + 2 * real(abs(ghat)^2 * ay1diff_vec' * az1diff_vec);
        Hessian_S(3, 2) = Hessian_S(2, 3);

        % omega_hat = omega_hat - stepsize * grad_S;
        omega_hat = omega_hat - Hessian_S \ grad_S;

    end
%     figure;
%     semilogy(1 : Tmax, obj, 'r');
    y_residue_tensor = reshape(y_residue, Nx, My, Lz);
end


clc; clear; close all;
rng(5);

% basic parameter set
sample_time = 60;
N_targets = 6;
Rnoise_mea = 1;
mu_c = 10;                   % the mean number of clutters (Possion)
state_dim = 4;
range_c = [-100, -100; 100, 100];
T = 1;                      % sample interval
sigma_w = 0;               % state transition noise variance
sigma_v = 1;               % measurement noise variance
d_chi = 16;

% the transition parameters matrix
Amat_state = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
Hmat_mea = [1 0 0 0; 0 0 1 0];
Rmat_noise = sigma_v * eye(2);
Tmat_temp = [1, 1/T; 1/T, 2/T^2];
Pmat_cov = kron(Rmat_noise, Tmat_temp);
V = wgn(1,2,0);
Qmat_noise = [T^4/4*(V'*V)  T^3/2*(V'*V) ; T^3/2*(V'*V) T^2*(V'*V)];%


% targets location and velocity initialization
xPosition_true0 = [-30; 30; 50; 0; -20; -50];
yPosition_true0 = [-30; 50; -50; -70; 70; 0];
xVelocity_true0 = [1; -1; -2; 0; -1; 2];
yVelocity_true0 = [1.5; -2; 1; 2; -2; 0];


% create the state array and measurement array
xPosseq_true = zeros(N_targets, sample_time);
yPosseq_true = zeros(N_targets, sample_time);
xVelseq_true = zeros(N_targets, sample_time);
yVelseq_true = zeros(N_targets, sample_time);

% Zkxseq_targets = zeros(N_targets, sample_time);
% Zkyseq_targets = zeros(N_targets, sample_time);

xPosseq_true(:, 1) = xPosition_true0;
yPosseq_true(:, 1) = yPosition_true0;
xVelseq_true(:, 1) = xVelocity_true0;
yVelseq_true(:, 1) = yVelocity_true0;

% Zkxseq_targets(:, 1) = xPosition_true0;
% Zkyseq_targets(:, 1) = yPosition_true0;

for k_time = 2 : sample_time
    Xtrue_nk = Amat_state * [xPosseq_true(:, k_time - 1)'; ...
        xVelseq_true(:, k_time - 1)'; yPosseq_true(:, k_time - 1)'; ...
        yVelseq_true(:, k_time - 1)'] + ...
        sqrt(sigma_w) * randn(state_dim, 1);
    xPosseq_true(:, k_time) = Xtrue_nk(1, :);
    xVelseq_true(:, k_time) = Xtrue_nk(2, :);
    yPosseq_true(:, k_time) = Xtrue_nk(3, :);
    yVelseq_true(:, k_time) = Xtrue_nk(4, :);
    % xPosseq_true(:, k_time) = xPosseq_true(:, k_time - 1) + ...
    % xVelseq_true(:, k_time - 1);
    % yPosseq_true(:, k_time) = yPosseq_true(:, k_time - 1) + ...
    % yVelseq_true(:, k_time - 1);
    % xVelseq_true(:, k_time) = xVelseq_true(:, k_time - 1);
    % yVelseq_true(:, k_time) = yVelseq_true(:, k_time - 1);

    % Xstate_truek = [xPosseq_true(:, k_time - 1)]
end

Zkxseq_targets = xPosseq_true + sqrt(sigma_v) * ...
    randn(N_targets, sample_time);
Zkyseq_targets = yPosseq_true + sqrt(sigma_v) * ...
    randn(N_targets, sample_time);
Zkcell_seq = cat(1);
cPosset_seq = cat(1);

% combine the observations and clutters
for k_time = 1 : sample_time
    N_c = poissrnd(mu_c);   % the number of clutters
    cPosset_k = repmat(range_c(1, :), [N_c 1]) + ...
        rand(N_c, 2) * diag([-1, 1] * range_c);      % create the clutters
    Zkcell_seq{k_time} = [Zkxseq_targets(:, k_time), ...
        Zkyseq_targets(:, k_time); cPosset_k];
    cPosset_seq{k_time} = cPosset_k;
end

% apply the tracking algorithm to the noisy observations
xPosseq_hat = zeros(N_targets, sample_time);
yPosseq_hat = zeros(N_targets, sample_time);
xVelseq_hat = zeros(N_targets, sample_time);
yVelseq_hat = zeros(N_targets, sample_time);

xPosseq_hat(:, 1) = xPosition_true0;
yPosseq_hat(:, 1) = yPosition_true0;
xVelseq_hat(:, 1) = xVelocity_true0;
yVelseq_hat(:, 1) = yVelocity_true0;



%% NN algorithm
for k_time = 2 : sample_time
    Pmat_cov = Amat_state * Pmat_cov * Amat_state' + Qmat_noise;
    Smat_k = Hmat_mea * Pmat_cov * Hmat_mea' + Rmat_noise;
    Kmat_gain = Pmat_cov * Hmat_mea' / (Hmat_mea * Pmat_cov * ...
        Hmat_mea' + Rmat_noise);
    for n_tidx = 1 : N_targets
        XPri_nk = Amat_state * [xPosseq_hat(n_tidx, k_time - 1); ...
            xVelseq_hat(n_tidx, k_time - 1); ...
            yPosseq_hat(n_tidx, k_time - 1); ...
            yVelseq_hat(n_tidx, k_time - 1)];
        ZPri_nk = (Hmat_mea * XPri_nk)';

        % dvec_nk = sqrt(sum((ZPri_nk - Zkcell_seq{k_time}) .^ 2, 2));
        dvec_nk = diag((ZPri_nk - Zkcell_seq{k_time}) / Smat_k * ...
        (ZPri_nk - Zkcell_seq{k_time})');
        [dmin_nk, didx_nk] = min(dvec_nk);
        if dmin_nk < d_chi
            % the gain of kalman filter
            ZPost_nk = Zkcell_seq{k_time}(didx_nk, :);
            XPost_nk = XPri_nk + Kmat_gain * (ZPost_nk' - Hmat_mea * XPri_nk);
            Pmat_cov = (eye(state_dim) - Kmat_gain * Hmat_mea) * Pmat_cov;
        else
            XPost_nk = XPri_nk;
        end
        xPosseq_hat(n_tidx, k_time) = XPost_nk(1);
        xVelseq_hat(n_tidx, k_time) = XPost_nk(2);
        yPosseq_hat(n_tidx, k_time) = XPost_nk(3);
        yVelseq_hat(n_tidx, k_time) = XPost_nk(4);
    end

end



% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;


fhandle_fig1 = figure(1);
title('Target trajectory', 'Fontsize', Fsz);
xlabel('x/m', 'Fontsize', Fsz);
ylabel('y/m', 'Fontsize', Fsz);
hold on;
xlim([-120 120])
ylim([-120 120])

N_show = 1 : N_targets;

for k_time = 1 : sample_time
    for n_tidx = N_show
        phandle_ture = plot(xPosseq_true(n_tidx, 1 : k_time), ...
        yPosseq_true(n_tidx, 1 : k_time), 'b', 'LineWidth', Lw, ...
        'Marker', 's', 'MarkerSize', Msz);
    end

    cPosset_k = cPosset_seq{k_time};
    for c_idx = 1 : size(cPosset_k, 1)
        phandle_clutter = plot(cPosset_k(c_idx, 1), cPosset_k(c_idx, 2),...
        'k*', 'LineWidth', Lw, 'Marker', 's', 'MarkerSize', Msz);
    end

    for n_tidx = N_show
        phandle_hat = plot(xPosseq_hat(n_tidx, 1 : k_time), ...
        yPosseq_hat(n_tidx, 1 : k_time), 'r', 'LineWidth', Lw, ...
        'Marker', 's', 'MarkerSize', Msz);
    end
    pause(0.1);
    % clf(fhandle_fig1)
end




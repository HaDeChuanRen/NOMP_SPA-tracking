clc; clear; close all;
rng(6);

% basic parameter set
sample_time = 60;
N_targets = 6;
Rnoise_mea = 1;
mu_c = 10;                   % the mean number of clutters (Possion)
state_dim = 4;
range_c = [-100, 0; 100, 100];
T = 1;                      % sample interval
sigma_w = [0.1; 0; 0.1; 0];               % state transition noise variance
sigma_v = [0.1; 0.1];               % measurement noise variance
d_chi = 16;
Pd_xki = 0.9999;
L_ite = 20;

% the transition parameters matrix
Amat_state = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
Hmat_mea = [1 0 0 0; 0 0 1 0];
Rmat_noise = eye(2) * diag(sigma_v);
Tmat_temp = [1, 1/T; 1/T, 2/T^2];
Pmat_cov = kron(Rmat_noise, Tmat_temp);
V = wgn(1,2,0);
Qmat_noise = [T^4/4*(V'*V)  T^3/2*(V'*V) ; T^3/2*(V'*V) T^2*(V'*V)];
% Question: what's this?


% targets location and velocity initialization
xPosition_true0 = [-30; 30; 45; 0; -20; -45];
yPosition_true0 = [3; 45; 0; 1; 70; 0];
xVelocity_true0 = [1; -0.3; -0.4; 0; -0.5; 2];
yVelocity_true0 = [1.5; -0.7; 0.6; 1; -1; 0.1];


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
        sqrt(sigma_w) .* randn(state_dim, 1);
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

Zkxseq_targets = xPosseq_true + sqrt(sigma_v(1)) * ...
    randn(N_targets, sample_time);
Zkyseq_targets = yPosseq_true + sqrt(sigma_v(2)) * ...
    randn(N_targets, sample_time);
Zkcell_seq = cat(1);
cPosset_seq = cat(1);

% combine the observations and clutters
for k_time = 1 : sample_time
    N_c = poissrnd(mu_c);   % the number of clutters
    cPosset_k = repmat(range_c(1, :), [N_c 1]) + ...
        rand(N_c, 2) * diag([-1, 1] * range_c);      % create the clutters
    cPosset_seq{k_time} = cPosset_k;
    Zkcell_seq{k_time} = [Zkxseq_targets(:, k_time), ...
        Zkyseq_targets(:, k_time); cPosset_k];
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


% SPA data assciation
for k_time = 2 : sample_time
    Pmat_cov = Amat_state * Pmat_cov * Amat_state' + Qmat_noise;
    Smat_k = Hmat_mea * Pmat_cov * Hmat_mea' + Rmat_noise;
    Kmat_gain = Pmat_cov * Hmat_mea' / (Hmat_mea * Pmat_cov * ...
        Hmat_mea' + Rmat_noise);
    XPost_last = [xPosseq_hat(:, k_time - 1)'; ...
        xVelseq_hat(:, k_time - 1)'; yPosseq_hat(:, k_time - 1)'; ...
        yVelseq_hat(:, k_time - 1)'];
    XPri_k = Amat_state * XPost_last;
    Mk_mea = size(Zkcell_seq{k_time}, 1);
    betamat_imk = zeros(Mk_mea + 1, N_targets);

    for n_tidx = 1 : N_targets
        ZPri_nk = Hmat_mea * XPri_k(:, n_tidx);
        evec_Znk = ZPri_nk' - Zkcell_seq{k_time};
        dvec_nk = diag(evec_Znk / Smat_k * evec_Znk');
        betamat_imk(1 : Mk_mea, n_tidx) = Pd_xki * (exp(- 1/2 * ...
            dvec_nk) / sum(exp(- 1/2 * dvec_nk)));
        betamat_imk(1 + Mk_mea, n_tidx) = 1 - Pd_xki;
        % [~, aPrivec_k(n_tidx)] = max(betamat_imk(:, n_tidx));
    end
    ximat_mik = [ones(Mk_mea, N_targets), mu_c/N_targets * ...
        ones(Mk_mea, 1)];
    phimat_itel = betamat_imk(1 : Mk_mea, :) ./ ...
        betamat_imk(1 + Mk_mea, :);
    vmat_itel = zeros(Mk_mea, N_targets);
    for l_iter = 1 : L_ite
        for m_midx = 1 : Mk_mea
            for i_tidx = 1 : N_targets
                vmat_itel(m_midx, i_tidx) = 1 / (mu_c/N_targets + ...
                sum(phimat_itel(m_midx, :)) - phimat_itel(m_midx, i_tidx)*...
                ximat_mik(m_midx, i_tidx));
            end
        end

        for i_tidx = 1 : N_targets
            for m_midx = 1 : Mk_mea
                phimat_itel(m_midx, i_tidx) = betamat_imk(m_midx, i_tidx)/...
                    (betamat_imk(Mk_mea + 1, i_tidx) + ...
                    sum(betamat_imk(1 : Mk_mea, i_tidx) .* ...
                    vmat_itel(:, i_tidx)) - betamat_imk(m_midx, i_tidx) * ...
                    vmat_itel(m_midx, i_tidx));
            end
        end
    end
    Promat_a = zeros(Mk_mea + 1, N_targets);
    Promat_a(1 : Mk_mea, :) = betamat_imk(1 : Mk_mea, :) .* vmat_itel ./ ...
        (betamat_imk(Mk_mea + 1, :) + sum(betamat_imk(1 : Mk_mea, :)...
        .* vmat_itel));
    Promat_a(1 + Mk_mea, :) = betamat_imk(Mk_mea + 1, :) ./ ...
        (betamat_imk(Mk_mea + 1, :) + sum(betamat_imk(1 : Mk_mea, :)...
        .* vmat_itel));

    % aPrivec_k = zeros(N_targets, 1);
    XPost_k = zeros(state_dim, N_targets);

    for n_tidx = 1 : N_targets
        [~, aPrivec_nk] = max(betamat_imk(:, n_tidx));
        ZPost_nk = Zkcell_seq{k_time}(aPrivec_nk, :);
        XPost_k(:, n_tidx) = XPri_k(:, n_tidx) + Kmat_gain * (ZPost_nk' - ...
        Hmat_mea * XPri_k(:, n_tidx));
    end
    Pmat_cov = (eye(state_dim) - Kmat_gain * Hmat_mea) * Pmat_cov;
    xPosseq_hat(:, k_time) = XPost_k(1, :)';
    xVelseq_hat(:, k_time) = XPost_k(2, :)';
    yPosseq_hat(:, k_time) = XPost_k(3, :)';
    yVelseq_hat(:, k_time) = XPost_k(4, :)';

end




%% NN algorithm
% for k_time = 2 : sample_time
%     Pmat_cov = Amat_state * Pmat_cov * Amat_state' + Qmat_noise;
%     Smat_k = Hmat_mea * Pmat_cov * Hmat_mea' + Rmat_noise;
%     Kmat_gain = Pmat_cov * Hmat_mea' / (Hmat_mea * Pmat_cov * ...
%         Hmat_mea' + Rmat_noise);
%     for n_tidx = 1 : N_targets
%         XPri_nk = Amat_state * [xPosseq_hat(n_tidx, k_time - 1); ...
%             xVelseq_hat(n_tidx, k_time - 1); ...
%             yPosseq_hat(n_tidx, k_time - 1); ...
%             yVelseq_hat(n_tidx, k_time - 1)];
%         ZPri_nk = (Hmat_mea * XPri_nk)';

%         % dvec_nk = sqrt(sum((ZPri_nk - Zkcell_seq{k_time}) .^ 2, 2));
%         dvec_nk = diag((ZPri_nk - Zkcell_seq{k_time}) / Smat_k * ...
%         (ZPri_nk - Zkcell_seq{k_time})');
%         [dmin_nk, didx_nk] = min(dvec_nk);
%         if dmin_nk < d_chi
%             % the gain of kalman filter
%             ZPost_nk = Zkcell_seq{k_time}(didx_nk, :);
%             XPost_nk = XPri_nk + Kmat_gain * (ZPost_nk' - Hmat_mea * XPri_nk);
%             Pmat_cov = (eye(state_dim) - Kmat_gain * Hmat_mea) * Pmat_cov;
%         else
%             XPost_nk = XPri_nk;
%         end
%         xPosseq_hat(n_tidx, k_time) = XPost_nk(1);
%         xVelseq_hat(n_tidx, k_time) = XPost_nk(2);
%         yPosseq_hat(n_tidx, k_time) = XPost_nk(3);
%         yVelseq_hat(n_tidx, k_time) = XPost_nk(4);
%     end

% end



% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;


fhandle_fig1 = figure(1);
title('SPADA Target trajectory', 'Fontsize', Fsz);
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




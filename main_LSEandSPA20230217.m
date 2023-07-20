% From line spectrum to a visual tracking result
% date: 2023/2/17

% Todo list: 规定正方向
% x坐标：右边为正
% y坐标：恒为正

clc; clear; close all;
rng(3);

% basic parameter set
sample_time = 60;
K_targets = 6;
Rnoise_mea = 1;
mu_c = 10;                               % the mean number of clutters (Possion)
state_dim = 4;
range_c = [-100, -100; 100, 100];                        % the range of clutters
T = 1;                                                         % sample interval
sigma_w = [0.1; 0; 0.1; 0];                    % state transition noise variance
sigma_v = [0.1; 0.1];                               % measurement noise variance
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
Qmat_noise = [T^4/4*(V'*V)  T^3/2*(V'*V); T^3/2*(V'*V) T^2*(V'*V)];


% targets location and velocity initialization
xPosition_true0 = [-30; 30; 45; 0; 20; -45];
yPosition_true0 = [3; 45; 10; 1; 70; 20];
xVelocity_true0 = [1; -0.3; -0.4; 0; -0.5; 2];
yVelocity_true0 = [1.5; -0.7; 0.6; 1; -1; 0.1];


% create the real state array
xPosseq_true = zeros(K_targets, sample_time);
yPosseq_true = zeros(K_targets, sample_time);
xVelseq_true = zeros(K_targets, sample_time);
yVelseq_true = zeros(K_targets, sample_time);


xPosseq_true(:, 1) = xPosition_true0;
yPosseq_true(:, 1) = yPosition_true0;
xVelseq_true(:, 1) = xVelocity_true0;
yVelseq_true(:, 1) = yVelocity_true0;


for t_time = 2 : sample_time
    Xtrue_kt = Amat_state * [xPosseq_true(:, t_time - 1)'; ...
        xVelseq_true(:, t_time - 1)'; yPosseq_true(:, t_time - 1)'; ...
        yVelseq_true(:, t_time - 1)'] + ...
        sqrt(sigma_w) .* randn(state_dim, 1);
    xPosseq_true(:, t_time) = Xtrue_kt(1, :);
    xVelseq_true(:, t_time) = Xtrue_kt(2, :);
    yPosseq_true(:, t_time) = Xtrue_kt(3, :);
    yVelseq_true(:, t_time) = Xtrue_kt(4, :);
end


% use the state array to create line spectrum in each time instant
% radar parameter
c = 299792458;
T_idle = 100e-6;
T_ramp = 60e-6;
T_circle = T_idle + T_ramp;
Fre_start = 77e9;
Slope_fre = 14.991e12;
Fs = 10e6;
Ts = 1 / Fs;
lambda_cw = c / Fre_start;
Rx_interval = lambda_cw / 2;

rmax = c / (2 * Ts * Slope_fre);
vmax = c / (4 * Fre_start * T_circle);


% line spectrum parameter
Nx = 256;
My = 64;
Lz = 8;
NML = Nx * My * Lz;
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);

snr_set = 30;
sigma_n = 1;


% initialize real r, v, theta parameters of each instant
rmat_true = zeros(K_targets, sample_time);
vmat_true = zeros(K_targets, sample_time);
thetamat_true = zeros(K_targets, sample_time);

% initialize real omega parameters of each dimension
% $\omega \in (-\pi, \pi]$
omegaxmat_true = zeros(K_targets, sample_time);
omegaymat_true = zeros(K_targets, sample_time);
omegazmat_true = zeros(K_targets, sample_time);
gain_true = sqrt(10 .^ (snr_set / 10) * sqrt(sigma_n)) .* ...
    exp(1j * 2 * pi * rand(K_targets, sample_time));

yten_spec = zeros(Nx, My, Lz, sample_time);

for t_time = 1 : sample_time
    % obtain the real measurements in the instant t
    rmat_true(:, t_time) = sqrt(xPosseq_true(:, t_time) .^ 2 + ...
        yPosseq_true(:, t_time) .^ 2);
    thetamat_true(:, t_time) = atan(xPosseq_true(:, t_time) ./ ...
        yPosseq_true(:, t_time));
    vmat_true(:, t_time) = xVelseq_true(:, t_time) .* ...
        cos(thetamat_true(:, t_time)) + yVelseq_true(:, t_time) .* ...
        sin(thetamat_true(:, t_time));

    % calculate the corresponding angular frequency
    % $\omega \in (-\pi, \pi]$
    omegaxmat_true(:, t_time) = wrapToPi(rmat_true(:, t_time) / rmax * 2 * pi);
    omegaymat_true(:, t_time) = vmat_true(:, t_time) / vmax * pi;
    omegazmat_true(:, t_time) = pi * sin(thetamat_true(:, t_time));

    % create the line spectrum tensors in the instant t
    yvec_t = zeros(NML, 1);
    for k_idx = 1 : K_targets
        yvec_kt = kron(kron(array_Fun(omegazmat_true(k_idx, t_time), Lz), ...
            array_Fun(omegaymat_true(k_idx, t_time), My)), ...
            array_Fun(omegaxmat_true(k_idx, t_time), Nx));
        yvec_t = yvec_t + gain_true(k_idx, t_time) * yvec_kt;
    end
    yvec_t = yvec_t + sqrt(1 / 2) * (randn(NML, 1) + 1j * randn(NML, 1));
    yten_spec(:, :, :, t_time) = reshape(yvec_t, Nx, My, Lz);
end


% apply the tracking algorithm to the noisy line spectrums
% initialize the sequence of target state estimation
xPosseq_hat = zeros(K_targets, sample_time);
yPosseq_hat = zeros(K_targets, sample_time);
xVelseq_hat = zeros(K_targets, sample_time);
yVelseq_hat = zeros(K_targets, sample_time);

% set the initial state as the real parameters
xPosseq_hat(:, 1) = xPosition_true0;
yPosseq_hat(:, 1) = yPosition_true0;
xVelseq_hat(:, 1) = xVelocity_true0;
yVelseq_hat(:, 1) = yVelocity_true0;

% set the parameters of NOMP-CFAR algorithm
P_fa = 0.01;
NM = Nx * My;
N_r = 60;
alpha_2D = alpha_Poe(P_fa, NM, N_r);
K_max = K_targets * 2;
gamma_l = 16;


tic;
for t_time = 2 : sample_time
    whandle = waitbar((t_time - 1) / sample_time);

%     if t_time == 44
%         1;
%     end

    % 2D-NOMP analysis for each antenna
    yten_t = squeeze(yten_spec(:, :, :, t_time));
    omegaxycell_t = cell(Lz, 1);
    ghatcell_t = cell(Lz, 1);
    for l_idx = 1 : Lz
        ymat_lt = squeeze(yten_t(:, :, l_idx));
        [omegahat_lt, ghat_lt, ~, ~] = NOMP2D_low(ymat_lt, alpha_2D, N_r, ...
            K_max);
        omegaxycell_t{l_idx} = wrapToPi(omegahat_lt);
        ghatcell_t{l_idx} = ghat_lt;
    end

    % match the measurements of each antenna
    % set the estimate number of measurements as the largest number of
    % measurements from each antenna
    Mhat_t = 0;
    for l_idx = 1 : Lz
        if Mhat_t < length(ghatcell_t{l_idx})
            Mhat_t = length(ghatcell_t{l_idx});
        end
    end

    for l_idx = 1 : Lz
        if size(omegaxycell_t{l_idx}, 1) < Mhat_t
            omegaxycell_t{l_idx} = [omegaxycell_t{l_idx}; zeros(Mhat_t - ...
                size(omegaxycell_t{l_idx}, 1), 2)];
            ghatcell_t{l_idx} = [ghatcell_t{l_idx}; zeros(Mhat_t - ...
                length(ghatcell_t{l_idx}), 1)];
        end
    end

    % initialize the matching result of each antenna
    gomegacell_K = zeros(Mhat_t, Lz, 3);
    % the sequnce of the 3 variables is ghat, omegax, omegay
    gomegacell_K(:, 1, 1) = ghatcell_t{1}(:);
    gomegacell_K(:, 1, 2 : 3) = omegaxycell_t{1}(:, :);

    % match the measurements of each antenna by Hungarian match algorithm
    for l_idx = 2 : Lz
        dist_lt = abs(omegaxycell_t{l_idx}(:, 1) - omegaxycell_t{1}(:, 1)');
        [mi_Hun, mj_Hun] = Hungarianmatch(dist_lt);
        gomegacell_K(:, l_idx, 1) = ghatcell_t{l_idx}(mi_Hun);
        gomegacell_K(:, l_idx, 2 : 3) = omegaxycell_t{l_idx}(mi_Hun, :);
    end

    % calculate the result of measurements with the matching results
    gomega_t = zeros(Mhat_t, 4);
    % the sequence of variables: ghat, omegax, omegay, omegaz
    Zarray_t = zeros(Mhat_t, 3);
    % the sequence of variables: xPosition, yPosition, radial velocity

    for m_idx = 1 : Mhat_t
        % calculate the average of omegax and omegay as the result of $m$th
        % measurement in the intant $t$
        omegaxyvec_mt = squeeze(mean(gomegacell_K(m_idx, :, 2 : 3), 2))';
        gomega_t(m_idx, 2 : 3) = omegaxyvec_mt;

        % use the average result to adjust the ghat results by LS method
        gvec_mt = zeros(Lz, 1);
        for l_idx = 1 : Lz
            gvec_mt(l_idx) = LeastSquares_2D(squeeze(yten_t(:, :, l_idx)), ...
                omegaxyvec_mt);
        end
        gomega_t(m_idx, 1) = gvec_mt(1);

        % squeeze(mean(gomegacell_K(m_idx, :, 2 : 3), 2));
        % gvec_mt = gomegacell_K(m_idx, :, 1);

        % calculate omega with gvec by fft method
        [~, z_idx] = max(abs(fft(gvec_mt, Lz * gamma_l)));
        omegazhat_mt = wrapToPi((z_idx / (Lz * gamma_l)) * 2 * pi);
        gomega_t(m_idx, 4) = omegazhat_mt;
        rhat_mt = wrapTo2Pi(gomega_t(m_idx, 2)) / (2 * pi) * rmax;
        theta_mt = asin(omegazhat_mt / pi);
        Zarray_t(m_idx, 1) = rhat_mt * sin(theta_mt);
        Zarray_t(m_idx, 2) = rhat_mt * cos(theta_mt);
        Zarray_t(m_idx, 3) = (gomega_t(m_idx, 3) / pi) * vmax;
    end

    % calculate the prior probability of the time t
    % the Kalman filter patameter
    Pmat_cov = Amat_state * Pmat_cov * Amat_state' + Qmat_noise;
    Smat_t = Hmat_mea * Pmat_cov * Hmat_mea' + Rmat_noise;
    Kmat_gain = Pmat_cov * Hmat_mea' / (Hmat_mea * Pmat_cov * Hmat_mea' + ...
        Rmat_noise);
    XPost_last = [xPosseq_hat(:, t_time - 1)'; xVelseq_hat(:, t_time - 1)'; ...
        yPosseq_hat(:, t_time - 1)'; yVelseq_hat(:, t_time - 1)'];
    XPri_t = Amat_state * XPost_last;
    Zxy_t = Zarray_t(:, 1: 2);

    % match the measurements and the targets by SPA algorithm
    % initialize the \beta matrix
    betamat_imk = zeros(Mhat_t + 1, K_targets);

    for k_tidx = 1 : K_targets
        ZPri_kt = Hmat_mea * XPri_t(:, k_tidx);
        evec_Zkt = ZPri_kt' - Zxy_t;
        dvec_kt = diag(evec_Zkt / Smat_t * evec_Zkt');
        betamat_imk(1 : Mhat_t, k_tidx) = Pd_xki * (exp(- 1/2 * dvec_kt) / ...
            sum(exp(- 1/2 * dvec_kt) + eps(0)));
        % here is some problem: when the target is missed by the measurements,
        % how to express the \beta ?
        betamat_imk(1 + Mhat_t, k_tidx) = 1 - Pd_xki;
    end

    ximat_mkt = [ones(Mhat_t, K_targets), mu_c / K_targets * ones(Mhat_t, 1)];
    phimat_itel = betamat_imk(1 : Mhat_t, :) ./ betamat_imk(1 + Mhat_t, :);
    vmat_itel = zeros(Mhat_t, K_targets);

    % SPA algorithm iterate
    for l_iter = 1 : L_ite
        for m_midx = 1 : Mhat_t
            for i_tidx = 1 : K_targets
                vmat_itel(m_midx, i_tidx) = 1 / (mu_c / K_targets + ...
                sum(phimat_itel(m_midx, :)) - phimat_itel(m_midx, i_tidx) * ...
                ximat_mkt(m_midx, i_tidx));
            end
        end

        for i_tidx = 1 : K_targets
            for m_midx = 1 : Mhat_t
                phimat_itel(m_midx, i_tidx) = betamat_imk(m_midx, i_tidx)/...
                    (betamat_imk(Mhat_t + 1, i_tidx) + ...
                    sum(betamat_imk(1 : Mhat_t, i_tidx) .* ...
                    vmat_itel(:, i_tidx)) - betamat_imk(m_midx, i_tidx) * ...
                    vmat_itel(m_midx, i_tidx));
            end
        end
    end

    Promat_a = zeros(Mhat_t + 1, K_targets);
    Promat_a(1 : Mhat_t, :) = betamat_imk(1 : Mhat_t, :) .* vmat_itel ./ ...
        (betamat_imk(Mhat_t + 1, :) + sum(betamat_imk(1 : Mhat_t, :)...
        .* vmat_itel));
    Promat_a(1 + Mhat_t, :) = betamat_imk(Mhat_t + 1, :) ./ ...
        (betamat_imk(Mhat_t + 1, :) + sum(betamat_imk(1 : Mhat_t, :)...
        .* vmat_itel));

    XPost_t = zeros(state_dim, K_targets);
    % the sequence of the variable is : p_x, v_x, p_x, v_y
    for k_tidx = 1 : K_targets
        [~, aPrivec_kt] = max(betamat_imk(:, k_tidx));

        if aPrivec_kt <= Mhat_t
            ZPost_kt = Zxy_t(aPrivec_kt, :);
            XPost_t(:, k_tidx) = XPri_t(:, k_tidx) + Kmat_gain * (ZPost_kt' ...
                - Hmat_mea * XPri_t(:, k_tidx));
            % the velocity verification
            vel_kt = Zarray_t(aPrivec_kt, 3);
            rhat_kt = sqrt(XPost_t(1, k_tidx) ^ 2 + XPost_t(3, k_tidx) ^ 2);
            vhat_kt = XPost_t(2, k_tidx) * (XPost_t(1, k_tidx) / rhat_kt) ...
                - XPost_t(4, k_tidx) * (XPost_t(3, k_tidx) / rhat_kt);
            if abs(vel_kt - vhat_kt) > 0.1
                XPost_t(:, k_tidx) = XPri_t(:, k_tidx) + Kmat_gain ...
                    * (ZPost_kt' - Hmat_mea * XPri_t(:, k_tidx));
            end
        else
            XPost_t(:, k_tidx) = XPri_t(:, k_tidx);
        end
    end
    Pmat_cov = (eye(state_dim) - Kmat_gain * Hmat_mea) * Pmat_cov;
    xPosseq_hat(:, t_time) = XPost_t(1, :)';
    xVelseq_hat(:, t_time) = XPost_t(2, :)';
    yPosseq_hat(:, t_time) = XPost_t(3, :)';
    yVelseq_hat(:, t_time) = XPost_t(4, :)';

end


toc;
delete(whandle);

% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;


fhandle_fig1 = figure(1);
title('SPADA Target trajectory', 'Fontsize', Fsz);
xlabel('x/m', 'Fontsize', Fsz);
ylabel('y/m', 'Fontsize', Fsz);

xlim([-120 120])
ylim([0 120])

K_show = 1 : K_targets;

for t_time = 1 : sample_time
    for k_kidx = K_show
        phandle_ture = plot(xPosseq_true(k_kidx, t_time), ...
        yPosseq_true(k_kidx, t_time), 'b', 'LineWidth', Lw, ...
        'Marker', 's', 'MarkerSize', Msz);
        hold on;
    end

    % cPosset_t = cPosset_seq{t_time};
    % for c_idx = 1 : size(cPosset_t, 1)
    %     phandle_clutter = plot(cPosset_t(c_idx, 1), cPosset_t(c_idx, 2),...
    %     'k*', 'LineWidth', Lw, 'Marker', 's', 'MarkerSize', Msz);
    % end

    for k_kidx = K_show
        phandle_hat = plot(xPosseq_hat(k_kidx, t_time), ...
        yPosseq_hat(k_kidx, t_time), 'r', 'LineWidth', Lw, ...
        'Marker', 's', 'MarkerSize', Msz);
    end
    pause(0.3);
    % clf(fhandle_fig1)
end


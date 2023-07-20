% obeject method test
% date: 2023/5/17

clc; clear; close all;
rng(6);
addpath('NOMP1D tools\')
addpath('NOMP2D tools\')
addpath('analysis tools\')

% targets location, velocity, appear time disapear time initialization
sample_time = 60;   % the total samples of tracking
K_targets = 6;
xrange = 10;
xmin = -5;
xPosition_true0 = xrange * rand(K_targets, 1) + xmin;
vxrange = 0.2;
vxmin = -0.1;
xVelocity_true0 = vxrange * rand(K_targets, 1) + vxmin;
yrange = 15;
ymin = 6;
yPosition_true0 = yrange * rand(K_targets, 1) + ymin;
vyrange = 0.4;
vymin = -0.1;
yVelocity_true0 = vyrange * rand(K_targets, 1) + vymin;


xlim_l = xmin + sample_time * min([vxmin, 0]);
xlim_r = xmin + xrange + sample_time * max([vxrange + vxmin, 0]);
ylim_d = ymin + sample_time * min([vymin, 0]);
ylim_u = ymin + yrange + sample_time * max([vyrange + vymin, 0]);
xlim_vec = xlim_l : xlim_r;
ybon_vec = abs(xlim_vec / 2);

% xPosition_true0 = [-7; -2; 5];
% xVelocity_true0 = [0.25; 0.1; -0.1];
% yPosition_true0 = [4; 2; 1];
% yVelocity_true0 = [0; 0.2; 0.1];
% xPosition_true0 = [-3; 3; 4.5; 0; 2; -4.5];
% yPosition_true0 = [0.3; 4.5; 1; 1; 7; 2];
% xVelocity_true0 = [0.1; -0.03; -0.04; 0; -0.05; 0.2];
% yVelocity_true0 = [0.15; -0.07; 0.06; 0.1; -0.1; 0.01];
% xPosition_true0 = 3 * [0; -5; -5; -5; -2; -2; -4; 6];
% xVelocity_true0 = 3 * [0; 0.1; 0.1; 0.1; 0.1; 0.08; 0; -0.02];
% yPosition_true0 = 3 * [5; 1; 0.5; 1.5; 4; 6; 2; 3];
% yVelocity_true0 = 3 * [-0.1; 0; 0.03; -0.03; 0; 0.04; -0.04; -0.06];
% K_targets = length(xPosition_true0);    % the number of targets



% basic parameter set
mu_c = 1;    % the mean number of clutters (Possion)
state_dim = 4;    % the number of state dimensions
range_c = [-100, -100; 100, 100];    % the range of clutters
T = 1;    % sample interval
sigma_w = [0.0001; 0; 0.0001; 0];    % state transition noise variance
sigma_v = [0.01; 0.01];    % measurement noise variance
Pd_xki = 0.9;    % detection probability
L_ite = 5;    % the iteration times
P_G = 0.9;
range2_PG = - 2 * log(1 - P_G);

% the transition parameters matrix
Amat_state = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];    % state transition matrix
Hmat_mea = [1 0 0 0; 0 0 1 0];    % measurment matrix
Rmat_noise = eye(2) * diag(sigma_v);    % measurement covriance
% Todo: use CRB to describe the Rmat_noise
Tmat_temp = [1, 1/T; 1/T, 2/T^2];
Pmat_cov = kron(Rmat_noise, Tmat_temp);
% Question: why?
V = wgn(1,2,0);
Qmat_noise = [T^4/4*(V'*V)  T^3/2*(V'*V); T^3/2*(V'*V) T^2*(V'*V)];




% create the real state array
xPosseq_true = zeros(K_targets, sample_time);
yPosseq_true = zeros(K_targets, sample_time);
xVelseq_true = zeros(K_targets, sample_time);
yVelseq_true = zeros(K_targets, sample_time);

% set the original states of all targets
xPosseq_true(:, 1) = xPosition_true0;
yPosseq_true(:, 1) = yPosition_true0;
xVelseq_true(:, 1) = xVelocity_true0;
yVelseq_true(:, 1) = yVelocity_true0;



% calculate the state at each instant
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

ground_all = zeros(state_dim, sample_time, K_targets);
ground_all(1, :, :) = xPosseq_true';
ground_all(2, :, :) = xVelseq_true';
ground_all(3, :, :) = yPosseq_true';
ground_all(4, :, :) = yVelseq_true';

% use the state array to create line spectrum in each time instant
% radar parameter
c = 299792458;    % electromagnetic speed
T_idle = 100e-6;    % idle time
T_ramp = 60e-6;    % ramp time
T_circle = T_idle + T_ramp;    % chirp period
Fre_start = 77e9;    % carrier frequency
Slope_fre = 14.991e12;    % chirp rate
Fs = 10e6;    % sample rate
Ts = 1 / Fs;     % sample interval
lambda_cw = c / Fre_start;    % wavelength
Rx_interval = lambda_cw / 2;    % antenna element spacing

rmax = c / (2 * Ts * Slope_fre);    % maximum radial range
vmax = c / (4 * Fre_start * T_circle);    % maximum radial velocity

% line spectrum parameter
Nx = 256;
My = 64;
Lz = 8;
NML = Nx * My * Lz;
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);
vdel = 2 * vmax / My;

% snr_set = 30;
P_CW = 1.78e-2;
G_t = 1;
G_r = 1;
lambda_cw = 4e-3;
sigma_RCS = 1;
L_loss = 1;
k_Boltz = 1.38e-23;
T_kel = 300;
F_R = 10^1.5;
B_IF = 5e6;
SNR_temp = NML * (P_CW * G_t * G_r * lambda_cw^2 * sigma_RCS) ./ ...
    ((4 * pi)^3 * L_loss * k_Boltz * T_kel * F_R * B_IF);
sigma_n = 1;
% Todo: the SNR and noise variance can be calculated by equations

% initialize real r, v, theta parameters of each instant
rmat_true = zeros(K_targets, sample_time);
vmat_true = zeros(K_targets, sample_time);
thetamat_true = zeros(K_targets, sample_time);


% initialize real omega parameters of each dimension
% $\omega \in (-\pi, \pi]$
omegaxmat_true = zeros(K_targets, sample_time);
omegaymat_true = zeros(K_targets, sample_time);
omegazmat_true = zeros(K_targets, sample_time);
% gain_true = sqrt(10 .^ (snr_set / 10) * sqrt(sigma_n)) .* ...
%     exp(1j * 2 * pi * rand(K_targets, sample_time));

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
%         if abs(thetamat_true(k_idx, t_time)) > pi / 3
%             continue;
%         end
        SNR_kt = 30;
        gain_kt = sqrt(10 ^ (SNR_kt / 10) * sqrt(sigma_n)) .* exp(1j * 2 * ...
            pi * rand());
        yvec_kt = kron(kron(array_Fun(omegazmat_true(k_idx, t_time), Lz), ...
            array_Fun(omegaymat_true(k_idx, t_time), My)), ...
            array_Fun(omegaxmat_true(k_idx, t_time), Nx));
        yvec_t = yvec_t + gain_kt * yvec_kt;
    end
    yvec_t = yvec_t + sqrt(1 / 2) * (randn(NML, 1) + 1j * randn(NML, 1));
    yten_spec(:, :, :, t_time) = reshape(yvec_t, Nx, My, Lz);
end



% set the parameters of NOMP-CFAR algorithm
P_fa = 0.001;
N_r = 60;
NM = Nx * My;
alpha_2D = alpha_PoebyS(P_fa, NM, N_r, Lz);
K_max = 16;

tic;
% [history_all, Zcol_cell] = NOMP_SPAtracking(yten_spec, alpha_2D, N_r, K_max,...
    % rmax, vmax, Pmat_cov, Qmat_noise, Amat_state, Hmat_mea, Rmat_noise);


Zcol_cell = cell(sample_time, 1);

% apply the tracking algorithm to the noisy line spectrums
% initialize the sequence of target state estimation
Kseq_hat = zeros(sample_time, 1);
history_all = [];
target_cell = [];

yten_1 = squeeze(yten_spec(:, :, :, 1));
Zarray_1 = LSE_NOMPanalyse(yten_1, alpha_2D, N_r, K_max, rmax, vmax);
Mhat_1 = size(Zarray_1, 1);
Zcol_cell{1} = Zarray_1;

if Mhat_1 > 0
    for m_idx = 1 : Mhat_1
        x_m1 = Zarray_1(m_idx, 1);
        y_m1 = Zarray_1(m_idx, 2);
        theta_m1 = atan(x_m1 / y_m1);
        vrel_m1 = Zarray_1(m_idx, 3);
        vx_m1 = vrel_m1 * sin(theta_m1);
        vy_m1 = vrel_m1 * cos(theta_m1);
        statevec_k = [x_m1, vx_m1, y_m1, vy_m1]';
        % Todo: use the velocity estimation to set the original state.
        % and following part
        targetobj_k = TargetState(statevec_k, Pmat_cov, Qmat_noise, ...
            Amat_state, 1, T);
        target_cell = [target_cell; targetobj_k];
    end
end
Kseq_hat(1) = size(target_cell, 1);

% Kseq_hat(1) = Khat_1;
% for k_idx = 1 : Khat_1
%     statevec_k = [Zarray_1(k_idx, 1), Zarray_1(k_idx, 2), 0, 0]';
%     targethat_k1 = TargetState(statevec_k, Pmat_cov, Qmat_noise, Amat_state, ...
%         1, T);
%     targetseq_hat{1}(end + 1) = targethat_k1;
% end


tic;
for t_time = 2 : sample_time
    whandle = waitbar((t_time - 1) / sample_time);

    yten_t = squeeze(yten_spec(:, :, :, t_time));
    % 2D-NOMP analysis for each antenna
    Zarray_t = LSE_NOMPanalyse(yten_t, alpha_2D, N_r, K_max, rmax, vmax);
    Mhat_t = size(Zarray_t, 1);
    Zcol_cell{t_time} = Zarray_t;

    Khat_last = Kseq_hat(t_time - 1);
    if Mhat_t == 0 && Khat_last == 0
        continue;
    elseif Mhat_t > 0 && Khat_last == 0
        for m_idx = 1 : Mhat_t
            x_mt = Zarray_t(m_idx, 1);
            y_mt = Zarray_t(m_idx, 2);
            theta_mt = atan(x_mt / y_mt);
            vrel_mt = Zarray_t(m_idx, 3);
            vx_mt = vrel_mt * sin(theta_mt);
            vy_mt = vrel_mt * cos(theta_mt);
            statevec_k = [x_mt, vx_mt, y_mt, vy_mt]';
            % statevec_k = [Zarray_t(m_idx, 1), 0, Zarray_t(m_idx, 2), 0]';
            targetobj_k = TargetState(statevec_k, Pmat_cov, Qmat_noise, ...
                Amat_state, 1, T);
            target_cell = [target_cell; targetobj_k];
        end
        Kseq_hat(t_time) = size(target_cell, 1);
    elseif Mhat_t == 0 && Khat_last > 0
        disapp_set = [];
        for k_idx = 1 : Khat_last
            [~, ~, target_cell(k_idx)] = target_cell(k_idx).state_transform();
            target_cell(k_idx) = target_cell(k_idx).disappearing();
            if target_cell(k_idx).appear_state == 0
                disapp_set = [disapp_set, k_idx];
            end
        end
        for d_idx = disapp_set
            last_d = target_cell(d_idx).last_time;
            app_d = t_time - last_d;
            hisvec_d = target_cell(d_idx).state_his;
            hisall_d = [nan(state_dim, app_d), hisvec_d, nan(state_dim, ...
                sample_time - t_time)];
            history_all = cat(3, history_all, hisall_d);
        end
        target_cell(disapp_set) = [];
        if isempty(target_cell)
            Kseq_hat(t_time) = 0;
        else
            Kseq_hat(t_time) = size(target_cell, 1);
        end
    elseif Mhat_t > 0 && Khat_last > 0
        Zxy_t = Zarray_t(:, 1 : 2);
        Zvrad_t = Zarray_t(:, 3);
        % Todo: JPDA data association and Kalman filter
        XPri_t = zeros(state_dim, Khat_last);
        for k_idx = 1 : Khat_last
            [XPri_kt, Pmat_pri, target_cell(k_idx)] = ...
                target_cell(k_idx).state_transform();
            XPri_t(:, k_idx) = XPri_kt;
        end
        Smat_t = Hmat_mea * Pmat_pri * Hmat_mea' + Rmat_noise;

        % SPA data association
        % [Promat_a, Promat_b] = SPA_DA(XPri_t, Zxy_t, Hmat_mea, Smat_t, ...
        %     Pd_xki, mu_c, L_ite);

        % initialize the \beta matrix
        betamat_imk = zeros(Mhat_t + 1, Khat_last);
        ximat_mkt = [ones(Mhat_t, Khat_last), mu_c / Khat_last * ...
            ones(Mhat_t, 1)];

        for k_tidx = 1 : Khat_last
            betamat_imk(1 + Mhat_t, k_tidx) = 1 - Pd_xki;

            % distance calculate and P_G verification
            ZPri_kt = Hmat_mea * XPri_t(:, k_tidx);
            evec_Zkt = ZPri_kt - Zxy_t';
            dvec_kt = diag(evec_Zkt' / Smat_t * evec_Zkt);
            Pdvec_zmt = exp(- 1/2 * dvec_kt);
            % Pdvec_zmt(Pdvec_zmt < P_G) = 0;
            betamat_imk(1 : Mhat_t, k_tidx) = Pd_xki * Pdvec_zmt;
            ximat_mkt(1 : Mhat_t, k_tidx) = Pdvec_zmt;
            % here is some problem: when the target is missed by the measurements,
            % how to express the \beta ?
        end

        phimat_itel = betamat_imk(1 : Mhat_t, :) ./ betamat_imk(1 + Mhat_t, :);
        vmat_itel = zeros(Mhat_t, Khat_last);
        Promat_a = betamat_imk;


        % Kalman filter
        disapp_set = [];
        akt_coll = [];
        for k_idx = 1 : Khat_last
            [~, aPrivec_kt] = max(Promat_a(:, k_idx));
            if aPrivec_kt <= Mhat_t
                ZPost_kt = Zxy_t(aPrivec_kt, :);
                akt_coll = [akt_coll, aPrivec_kt];
                Kmat_gain = Pmat_pri * Hmat_mea' / (Hmat_mea * Pmat_pri * ...
                    Hmat_mea' + Rmat_noise);
                target_cell(k_idx) = ...
                    target_cell(k_idx).state_estimate(Rmat_noise, ...
                    Hmat_mea, ZPost_kt);
            elseif aPrivec_kt == Mhat_t + 1
                target_cell(k_idx) = target_cell(k_idx).disappearing();
                if target_cell(k_idx).appear_state == 0
                    disapp_set = [disapp_set, k_idx];
                end
            end
        end

        for d_idx = disapp_set
            last_d = target_cell(d_idx).last_time;
            app_d = t_time - last_d;
            hisvec_d = target_cell(d_idx).state_his;
            hisall_d = [nan(state_dim, app_d), hisvec_d, nan(state_dim, ...
                sample_time - t_time)];
            history_all = cat(3, history_all, hisall_d);
        end
        target_cell(disapp_set) = [];
        % create new target objects
        newZ_ind = setdiff(1 : Mhat_t, akt_coll);
        for m_idx = newZ_ind
            x_mt = Zarray_t(m_idx, 1);
            y_mt = Zarray_t(m_idx, 2);
            theta_mt = atan(x_mt / y_mt);
            vrel_mt = Zarray_t(m_idx, 3);
            vx_mt = vrel_mt * sin(theta_mt);
            vy_mt = vrel_mt * cos(theta_mt);
            statevec_k = [x_mt, vx_mt, y_mt, vy_mt]';
            % statevec_k = [Zarray_t(m_idx, 1), 0, Zarray_t(m_idx, 2), 0]';
            targetobj_k = TargetState(statevec_k, Pmat_cov, Qmat_noise, ...
                Amat_state, 1, T);
            target_cell = [target_cell; targetobj_k];
        end
        Kseq_hat(t_time) = size(target_cell, 1);
    end
end

if ~isempty(target_cell)
    Khat_last = Kseq_hat(end);
    for k_idx = 1 : Khat_last
        last_k = target_cell(k_idx).last_time;
        app_k = t_time - last_k;
        hisvec_k = target_cell(k_idx).state_his;
        hisall_k = [nan(state_dim, app_k), hisvec_k];
        history_all = cat(3, history_all, hisall_k);
    end
end



K_count = size(history_all, 3);
save('true_history.mat', 'history_all', 'xPosseq_true', 'xVelseq_true', ...
    'yPosseq_true', 'yVelseq_true')

toc;
% delete(whandle);

K_del = [];
for k_idx = 1 : K_count
    testvec_k = squeeze(history_all(1, :, k_idx));
    if sum(isnan(testvec_k)) > sample_time * 0.95
        K_del = [K_del, k_idx];
    end
end
history_all(:, :, K_del) = [];
K_count = size(history_all, 3);

Delta_thr = 10000;
[OSPA_dist, label_vec] = OSPA_cal(ground_all, history_all, Delta_thr);

% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;


fhandle_fig1 = figure(1);
xlabel('x/m', 'Fontsize', Fsz);
ylabel('y/m', 'Fontsize', Fsz);
xlim([xlim_l xlim_r])
ylim([ylim_d ylim_u])

hold on;

filename = "testavi.mp4";
test_video = VideoWriter(filename, 'MPEG-4');
test_video.FrameRate = 5;
open(test_video);

% color_vec = ''

color_vec = ['-or'; '-og'; '-om'; '-oc'; '-+r'; '-+g'; '-+m'; '-+c'];
num_color = 8;

% , 'DisplayName', 'estimate path'
for t_time = 1 : sample_time
    phandle_bond = plot(xlim_vec, ybon_vec, '--k', 'LineWidth', Lw);
    for k_tidx = 1 : K_targets
        phandle_ture = plot(xPosseq_true(k_tidx, 1 : t_time), yPosseq_true(k_tidx, ...
            1 : t_time), '-bo', 'LineWidth', Lw, 'Marker', 's', 'MarkerSize', Msz);
    end

    for k_idx = 1 : K_count
        c_idx = mod(k_idx, num_color) + 1;
        phandle_hat = plot(history_all(1, 1 : t_time, k_idx), history_all(3, ...
            1 : t_time, k_idx), color_vec(c_idx, :), 'LineWidth', Lw, 'Marker', ...
            's', 'MarkerSize', Msz);
    end

    if ~isempty(Zcol_cell{t_time})
        phandle_mea = plot(Zcol_cell{t_time}(:, 1), Zcol_cell{t_time}(:, 2), ...
            'k*');
    end
    pause(0.1);
    frame_t = getframe(fhandle_fig1);
    writeVideo(test_video, frame_t);
    % im_seq{t_time} = frame2im(frame_t);
end
close(test_video)

figure;
plot(1 : sample_time, OSPA_dist, '-b', 'LineWidth', Lw, 'Marker', 's', ...
    'MarkerSize', Msz)
xlabel('time(s)', 'Fontsize', Fsz)
ylabel('OSPA', 'Fontsize', Fsz)



% legend('Fontsize', Fsz);
% 
% filename = "testgif.gif";
% for t_time = 1 : sample_time - 1
%     [A_t, map_t] = rgb2ind(im_seq{t_time}, 256);
%     if t_time == 1
%         imwrite(A_t, map_t, filename, 'gif', 'LoopCount', Inf, 'DelayTime', ...
%             0.1);
%     else
%         imwrite(A_t, map_t, filename, 'gif', 'WriteMode', 'append', ...
%             'DelayTime',0.1);
%     end
% end



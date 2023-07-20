clc; clear; close all;

%% obtain file path and read the ADC data
addpath("analysis tools\")
addpath("NOMP1D tools\")
addpath("NOMP2D tools\")
addpath("Tracking tools\")
orginal_path = 'C:\study\202210target tracking\3data\20230717exp';
exp_type = '';
exp_serial = '\12';
filename_bin = [orginal_path, exp_type, exp_serial, '\adc_data.bin'];

% tar1 = [0, 0.6; 0, 10];
% tar2 = [1 * 0.5, 8 * 0.5 + 1.24; 1 * 0.5, 1.24];
% tar_true = cat(3, tar1, tar2);
tar_true = [];
K_targets = 0;


Nx = 128;
My = 64;
Lz = 4;
T = 0.1;    % sample interval
options.numRX = Lz;
options.numADCBits = 16;
options.numADCSamples = Nx;
options.isReal = 0;    % set to 1 if real only data, 0 if complex data0
data_cube = readadc(filename_bin, options);

num_chirps = size(data_cube, 2);
sample_time = num_chirps / My;
yten_spec = zeros(Nx, My, Lz, sample_time);

for t_time = 1 : sample_time
    yten_spec(:, :, :, t_time) = data_cube(:, ((t_time - 1) * My + ...
        1 : t_time * My), :);
end

%% set the radar parameter
c = 3e8;
T_idle = 3e-6;
T_ramp = 56e-6;
T_circle = T_idle + T_ramp;
Fre_start = 77e9;
Slope_fre = 8.014e12;
Fs = 5e6;
Ts = 1 / Fs;
lambda_cw = c / Fre_start;
Rx_interval = lambda_cw / 2;

rmax = c / (2 * Ts * Slope_fre);    % maximum radial range
vmax = c / (4 * Fre_start * T_circle);    % maximum radial velocity

% the range of observation scenario
xlim_l = -10;
xlim_r = 10;
ylim_u = 20;
xlim_vec = xlim_l : xlim_r;
ybon_vec = abs(xlim_vec / 2);


filename_mat = [orginal_path, exp_type, exp_serial, exp_serial, '.mat'];
save(filename_mat, 'Nx', 'My', 'Lz', 'T', 'tar_true', 'yten_spec', 'T_idle', ...
    'T_ramp', 'T_circle', 'Fs', 'xlim_l', 'xlim_r', 'ylim_u');
%% set the observation scene
% NL_num = Nx * Lz;
% NML = Nx * My * Lz;
P_G = 0.5;   % gate probability
mu_c = 1;   % the mean number of clutters (Possion)
state_dim = 4;    % the number of state dimensions
sigma_v = [0.1; 0.1];    % measurement noise variance
Pd_xki = 0.9;    % detection probability
L_ite = 20;    % the iteration times

db_num = 2;
db_dis = 1;






%% the transition parameters matrix
Amat_state = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];    % state transition matrix
Hmat_mea = [1 0 0 0; 0 0 1 0];    % measurment matrix
Rmat_noise = eye(2) * diag(sigma_v);    % measurement covriance
% Todo: use CRB to describe the Rmat_noise
Tmat_temp = [1, 1/T; 1/T, 2/T^2];
Pmat_cov = kron(Rmat_noise, Tmat_temp);
% Question: why?
V = wgn(1,2,0);
Qmat_noise = [T^4/4*(V'*V)  T^3/2*(V'*V); T^3/2*(V'*V) T^2*(V'*V)];

%% set the parameters of NOMP-CFAR algorithm
P_fa = 0.01;
NM = Nx * My;
N_r = 60;
alpha_2D = alpha_PoebyS(P_fa, NM, N_r, Lz);
K_max = 12;
gamma_l = 128;

% apply the tracking algorithm to the noisy line spectrums
% initialize the sequence of target state estimation
Kseq_hat = zeros(sample_time, 1);
history_all = [];
label_all = [];
target_cell = [];

yten_1 = squeeze(yten_spec(:, :, :, 1));
% [omegaxy_t, ghat_t, ~, ~] = NOMP2D_low(yten_1, alpha_2D, N_r, K_max);
% [Mhat_t, ~] = size(omegaxy_t);

% if Mhat_t > 0
%     Zarray_lse = zeros(Mhat_t, 3);
%     rhat_mt = wrapTo2Pi(omegaxy_t(:, 1)) / (2 * pi) * rmax;
%     % detect gain and frequency of an additional sinusoid
%     S = eye(Lz);
%     sampledManifold = preProcessMeasMat(S, gamma_l);
%     omegazhat_t = zeros(Mhat_t, 1);
%     R_s = 1;
%     for m_idx = 1 : Mhat_t
%         y_r = ghat_t(m_idx, :)';
%         [omegazhat_mt, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);
%         for i = 1 : R_s
%             [omegazhat_mt, ~, y_r] = refineOne(y_r, omegazhat_mt, ...
%                 gain_new, S, sampledManifold.ant_idx, true);
%         end
%         omegazhat_t(m_idx) = wrapToPi(omegazhat_mt);
%     end
%     % [~, z_idx] = max(abs(fft(ghat_t, gamma_l, 2)), [], 2);
%     % omegazhat_t = wrapToPi(((z_idx - 1) / gamma_l) * 2 * pi);

%     theta_mt = asin(- omegazhat_t / pi);
%     Zarray_lse(:, 1) = rhat_mt .* sin(theta_mt);
%     Zarray_lse(:, 2) = rhat_mt .* cos(theta_mt);
%     Zarray_lse(:, 3) = (wrapToPi(omegaxy_t(:, 2)) / pi) * vmax;
% else
%     Zarray_lse = [];
% end


Zarray_lse = LSE_NOMPanalyse(yten_1, alpha_2D, N_r, K_max, rmax, vmax);
Zarray_1 = meas_cluster(Zarray_lse);
% Zarray_lse = Zarray_1;

% dbscan_indt = dbscan(Zarray_lse, db_dis, db_num);
% [~, Mhat_t] = max(dbscan_indt);
% Zarray_1 = [];
% for m_idx = 1 : max(dbscan_indt)
%     Zarray_tm = Zarray_lse(dbscan_indt == m_idx, :);
%     Zarray_1(m_idx, :) = mean(Zarray_tm, 1);
% end
% Zarray_1 = [Zarray_1; Zarray_lse(dbscan_indt == -1, :)];
% Zarray_1 = Zarray_lse;
Mhat_1 = size(Zarray_1, 1);
Zcol_cell{1} = Zarray_lse;

% Zcol_cell{1} = Zarray_1;
num_label = 1;

if Mhat_1 > 0
    for m_idx = 1 : Mhat_1
        statevec_k = [Zarray_1(m_idx, 1), 0, Zarray_1(m_idx, 2), 0]';
        % Todo: use the velocity estimation to set the original state.
        % and following part
        targetobj_k = TargetState(statevec_k, num_label, Pmat_cov, ...
            Qmat_noise, Amat_state, 1, T);
        num_label = num_label + 1;
        target_cell = [target_cell; targetobj_k];
    end
end
Kseq_hat(1) = size(target_cell, 1);


tic;
for t_time = 2 : sample_time
    whandle = waitbar((t_time - 1) / sample_time);

    yten_t = squeeze(yten_spec(:, :, :, t_time));
    % 2D-NOMP analysis for each antenna
    Zarray_lse = LSE_NOMPanalyse(yten_t, alpha_2D, N_r, K_max, rmax, vmax);
    % [omegaxy_t, ghat_t, ~, ~] = NOMP2D_low(yten_t, alpha_2D, N_r, K_max);
    % tau = sigma_n * chi2inv((1 - 0.001) ^ (1 / NM), 2 * Lz) / 2;
    % [omegaxy_t1, ghat_t1, ~] = MNOMP2D(yten_t, tau);
    % [omegaxy_t, ghat_t, ~, ~] = NOMP2D_low(yten_t, alpha_2D, N_r, K_max);
    % [Mhat_t, ~] = size(omegaxy_t);

    % if Mhat_t > 0
    %     Zarray_lse = zeros(Mhat_t, 3);
    %     rhat_mt = wrapTo2Pi(omegaxy_t(:, 1)) / (2 * pi) * rmax;
    %     % detect gain and frequency of an additional sinusoid
    %     S = eye(Lz);
    %     sampledManifold = preProcessMeasMat(S, gamma_l);
    %     omegazhat_t = zeros(Mhat_t, 1);
    %     R_s = 1;
    %     for m_idx = 1 : Mhat_t
    %         y_r = ghat_t(m_idx, :)';
    %         [omegazhat_mt, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);
    %         for i = 1 : R_s
    %             [omegazhat_mt, ~, y_r] = refineOne(y_r, omegazhat_mt, ...
    %                 gain_new, S, sampledManifold.ant_idx, true);
    %         end
    %         omegazhat_t(m_idx) = wrapToPi(omegazhat_mt);
    %     end

    %     theta_mt = asin(- omegazhat_t / pi);
    %     Zarray_lse(:, 1) = rhat_mt .* sin(theta_mt);
    %     Zarray_lse(:, 2) = rhat_mt .* cos(theta_mt);
    %     Zarray_lse(:, 3) = (wrapToPi(omegaxy_t(:, 2)) / pi) * vmax;
    % else
    %     Zarray_lse = [];
    % end

    % dbscan_indt = dbscan(Zarray_lse, db_dis, db_num);
    % Zarray_t = [];
    % for m_idx = 1 : max(dbscan_indt)
    %     Zarray_tm = Zarray_lse(dbscan_indt == m_idx, :);
    %     Zarray_t(m_idx, :) = mean(Zarray_tm, 1);
    % end
    % Zarray_t = [Zarray_t; Zarray_lse(dbscan_indt == -1, :)];
    % Zarray_t = Zarray_lse;
    Zarray_t = meas_cluster(Zarray_lse);
    Mhat_t = size(Zarray_t, 1);
    Zcol_cell{t_time} = Zarray_lse;
    Khat_last = Kseq_hat(t_time - 1);

    if Mhat_t == 0 && Khat_last == 0
        continue;
    elseif Mhat_t > 0 && Khat_last == 0
        for m_idx = 1 : Mhat_t
            statevec_k = [Zarray_t(m_idx, 1), 0, Zarray_t(m_idx, 2), 0]';
            targetobj_k = TargetState(statevec_k, num_label, Pmat_cov, ...
                Qmat_noise, Amat_state, 1, T);
            num_label = num_label + 1;
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
            label_d = target_cell(d_idx).n_label;
            label_all = [label_all; label_d];
        end
        target_cell(disapp_set) = [];
        if isempty(target_cell)
            Kseq_hat(t_time) = 0;
        else
            Kseq_hat(t_time) = size(target_cell, 1);
        end
    elseif Mhat_t > 0 && Khat_last > 0
        Zxy_t = Zarray_t(:, 1 : 2);
        % Todo: SPA data association and Kalman filter
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
            ZPri_kt = Hmat_mea * XPri_t(:, k_tidx);
            evec_Zkt = ZPri_kt - Zxy_t';
            dvec_kt = diag(evec_Zkt' / Smat_t * evec_Zkt);
            Pdvec_zmt = exp(- 1/2 * dvec_kt);
            Pdvec_zmt(Pdvec_zmt < P_G) = 0;
            betamat_imk(1 : Mhat_t, k_tidx) = Pd_xki * Pdvec_zmt;
            ximat_mkt(1 : Mhat_t, k_tidx) = Pdvec_zmt;
            % here is some problem: when the target is missed by the measurements,
            % how to express the \beta ?
        end

        phimat_itel = betamat_imk(1 : Mhat_t, :) ./ betamat_imk(1 + Mhat_t, :);
        vmat_itel = zeros(Mhat_t, Khat_last);

        % SPA algorithm iterate
        for l_iter = 1 : L_ite
            for m_midx = 1 : Mhat_t
                for i_tidx = 1 : Khat_last
                    vmat_itel(m_midx, i_tidx) = 1 / (mu_c / Khat_last + ...
                    sum(phimat_itel(m_midx, :)) - phimat_itel(m_midx, i_tidx) * ...
                    ximat_mkt(m_midx, i_tidx));
                end
            end

            for i_tidx = 1 : Khat_last
                for m_midx = 1 : Mhat_t
                    phimat_itel(m_midx, i_tidx) = betamat_imk(m_midx, i_tidx)/...
                        (betamat_imk(Mhat_t + 1, i_tidx) + ...
                        sum(betamat_imk(1 : Mhat_t, i_tidx) .* ...
                        vmat_itel(:, i_tidx)) - betamat_imk(m_midx, i_tidx) * ...
                        vmat_itel(m_midx, i_tidx));
                end
            end
        end

        Promat_a = zeros(Mhat_t + 1, Khat_last);
        Promat_a(1 : Mhat_t, :) = betamat_imk(1 : Mhat_t, :) .* vmat_itel ./ ...
            (betamat_imk(Mhat_t + 1, :) + sum(betamat_imk(1 : Mhat_t, :)...
            .* vmat_itel));
        Promat_a(Mhat_t + 1, :) = betamat_imk(Mhat_t + 1, :) ./ ...
            (betamat_imk(Mhat_t + 1, :) + sum(betamat_imk(1 : Mhat_t, :)...
            .* vmat_itel));

        Promat_b = zeros(Mhat_t, Khat_last + 1);
        Promat_b(:, 1 : Khat_last) = ximat_mkt(:, 1 : Khat_last) .* phimat_itel ./...
            (ximat_mkt(:, Khat_last + 1) + sum(ximat_mkt(:, 1 : Khat_last) .* ...
            phimat_itel));
        Promat_b(:, Khat_last + 1) = ximat_mkt(:, Khat_last + 1) ./ (ximat_mkt(:, ...
            Khat_last + 1) + sum(ximat_mkt(:, 1 : Khat_last) .* phimat_itel, 2));

        % match
        Mmat_a = zeros(Mhat_t + 1, Khat_last);
        Mmat_b = zeros(Mhat_t, Khat_last + 1);
        for k_idx = 1 : Khat_last
            [~, aind_kt] = max(Promat_a(:, k_idx));
            Mmat_a(aind_kt, k_idx) = 1;
        end
        for m_idx = 1 : Mhat_t
            [~, bind_mt] = max(Promat_b(m_idx, :));
            Mmat_b(m_idx, bind_mt) = 1;
        end

        Mmat_ab = zeros(Mhat_t + 1, Khat_last + 1);
        if sum(xor(Mmat_a(1 : Mhat_t, :), Mmat_b(:, 1 : Khat_last)), 'all')
            if sum(Mmat_a(1 : Mhat_t, :)) < Mmat_b(:, 1 : Khat_last)
            end
        end


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
            label_d = target_cell(d_idx).n_label;
            label_all = [label_all; label_d];
        end
        target_cell(disapp_set) = [];
        % create new target objects
        newZ_ind = setdiff(1 : Mhat_t, akt_coll);
        for m_idx = newZ_ind
            statevec_k = [Zarray_t(m_idx, 1), 0, Zarray_t(m_idx, 2), 0]';
            targetobj_k = TargetState(statevec_k, num_label, Pmat_cov, ...
                Qmat_noise, Amat_state, 1, T);
            num_label = num_label + 1;
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
        label_k = target_cell(k_idx).n_label;
        label_all = [label_all; label_k];
    end
end

nan_rate = 0.8;
[history_all, label_all, K_est] = del_clutter(history_all, label_all, nan_rate);
K_count = size(history_all, 3);

toc;
delete(whandle);

% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;


fhandle_fig1 = figure(1);
xlabel('x/m', 'Fontsize', Fsz);
ylabel('y/m', 'Fontsize', Fsz);
xlim([xlim_l xlim_r])
ylim([0 ylim_u])

hold on;

filename = ['result', exp_serial(2 : end), '.mp4'];
test_video = VideoWriter(filename, 'MPEG-4');
test_video.FrameRate = 5;
open(test_video);


color_est = ['-ro'; '-go'; '-yo'; '-co'; '-rd'; '-gd'; '-yd'; '-cd'];
color_true = ['--b'; '--k'; '--m'];
num_coest = 8;
num_cotr = 3;

% , 'DisplayName', 'estimate path'
for t_time = 1 : sample_time
    phandle_bond = plot(xlim_vec, ybon_vec, '-k', 'LineWidth', Lw);
    for k_tidx = 1 : K_targets
        c_tidx = mod(k_tidx, num_cotr) + 1;
        tar_k = tar_true(:, :, k_tidx);
        phandle_ture = plot(tar_k(:, 1), tar_k(:, 2), color_true(c_tidx, :), ...
            'LineWidth', Lw, 'Marker', 's', 'MarkerSize', Msz);
    end

    tracker_t = squeeze(history_all(1, t_time, :));
    K_need = ~isnan(tracker_t);
    for k_idx = 1 : K_count
        if K_need(k_idx) == 1
            c_idx = mod(label_all(k_idx), num_coest) + 1;
            phandle_hat = plot(history_all(1, 1 : t_time, k_idx), ...
                history_all(3, 1 : t_time, k_idx), color_est(c_idx, :), ...
                'LineWidth', Lw, 'Marker', 's', 'MarkerSize', Msz);
        end
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

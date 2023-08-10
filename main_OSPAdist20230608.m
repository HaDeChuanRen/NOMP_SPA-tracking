clc; clear; close all;
% Last update: 2023/7/26
% 2023/7/26 update details:
% 1. (Todo) test and debug the OSPA_cal function

% add the path of NOMP tools and analysis tools
addpath('NOMP1D tools\')
addpath('NOMP2D tools\')
addpath('analysis tools\')
addpath('Tracking tools\')

% load the true trace and estimated trace
load('truth_his.mat');
Delta_thr = 1000;
[state_dim, sample_time, K_targets] = size(ground_all);
K_est = size(history_all, 3);
nan_rate = 0.8;
D_delta = 10;
% [history_all, label_all, K_est] = del_clutter(history_all, label_all, nan_rate);
[OSPA_dist, ~] = OSPA_cal(ground_all, history_all, label_all, D_delta);


% [K_targets, ~] = size(xPosseq_true);

% ground_all = zeros(state_dim, sample_time, K_targets);
% ground_all(1, :, :) = xPosseq_true';
% ground_all(2, :, :) = xVelseq_true';
% ground_all(3, :, :) = yPosseq_true';
% ground_all(4, :, :) = yVelseq_true';

% res_rate = 0.95;
% [history_copy, K_est] = del_clutter(history_all, res_rate);

% history_copy = history_all;
% K_del = [];
% for k_idx = 1 : K_est
%     testvec_k = squeeze(history_copy(1, :, k_idx));
%     if sum(isnan(testvec_k)) > sample_time * 0.98
%         K_del = [K_del, k_idx];
%     end
% end
% history_copy(:, :, K_del) = [];
% K_est = size(history_copy, 3);

% start_set = zeros(state_dim + 1, K_est);
% end_set = zeros(state_dim + 1, K_est);

% for k_idx = 1 : K_est
%     inval_vec = find(~isnan(history_copy(1, :, k_idx)));
%     start_ins = inval_vec(1);
%     end_ins = inval_vec(end);
%     start_set(1 : state_dim, k_idx) = history_copy(:, start_ins, k_idx);
%     start_set(state_dim + 1, k_idx) = start_ins;
%     end_set(1 : state_dim, k_idx) = history_copy(:, end_ins, k_idx);
%     end_set(state_dim + 1, k_idx) = end_ins;
% end

% esmap_mat = [];
% time_thr = 9;
% state_thr = 16;
% kdel_mat = [];
% kend_idx = [];
% for k_idx = 1 : K_est
%     kend_mat = end_set(:, k_idx);
%     if kend_mat(end) == sample_time
%         continue;
%     end
%     kdelta_mat = kend_mat - start_set;
%     kdelta_mat(:, k_idx) = NaN;
%     kdel_vec = sum(kdelta_mat(1 : state_dim, :) .^ 2);
%     ktimedel_dec = kdelta_mat(state_dim + 1, :) .^ 2;
%     valid_vec = ((ktimedel_dec > time_thr) + (kdel_vec > state_thr)) > 0;
%     kdel_vec(valid_vec) = NaN;
%     kdel_vec = kdel_vec + ktimedel_dec;

%     if sum(isnan(kdel_vec)) < K_est
%         kend_idx = [kend_idx; k_idx];
%         kdel_mat = [kdel_mat; kdel_vec];
%     end
%     % if abs(kdelta_vec(state_dim + 1, kmin_idx)) < time_thr && ...
%     %     norm(kdelta_vec(1 : state_dim, kmin_idx)) < state_thr
%     %     if kmin_idx >= k_idx
%     %         kmin_idx = kmin_idx + 1;
%     %     end
%     %     esmap_mat = [esmap_mat; k_idx, kmin_idx];
%     % end
%     % k= find(norm(kdelta_vec(1 : state_dim, :)) < state_thr && ...
%     %     norm(kdelta_vec(state_dim + 1, :)) < time_thr)
% end
% % k_valend = length(kend_idx);
% if ~isempty(kdel_mat)
%     k_valsta = (1 : K_est) .* (sum(isnan(kdel_mat), 1) < length(kend_idx));
%     kdel_mat(:, (k_valsta == 0)) = [];
%     k_valsta(k_valsta == 0) = [];

%     kdel_mat(isnan(kdel_mat)) = 100000;
%     [end_idx, start_idx] = Hungarianmatch(kdel_mat);
%     save('kdel_mat.mat', 'kdel_mat');
%     esmap_mat = [kend_idx(end_idx), k_valsta(start_idx)'];
% end

% P_count = size(esmap_mat, 1);
% K_est = K_est - P_count;
% if ~isempty(esmap_mat)
%     history_new = zeros(state_dim, sample_time, K_est);
%     for p_idx = 1 : P_count
%         espair = esmap_mat(p_idx, :);
%         his_for = history_copy(:, :, espair(1));
%         his_lat = history_copy(:, :, espair(2));
%         tend_for = end_set(state_dim + 1, espair(1));
%         tstart_lat = start_set(state_dim + 1, espair(2));
%         his_new = zeros(state_dim, sample_time);
%         if tend_for >= tstart_lat
%             % average
%             his_new(:, 1 : tstart_lat - 1) = his_for(:, 1 : tstart_lat - 1);
%             his_new(:, tstart_lat : tend_for) = (his_for(:, tstart_lat : ...
%                 tend_for) + his_lat(:, tstart_lat : tend_for)) / 2;
%             his_new(:, tend_for + 1 : end) = his_lat(:, tend_for + 1 : end);
%         else
%             % interpolation
%             his_new(:, 1 : tend_for) = his_for(:, 1 : tend_for);
%             del_es = tstart_lat - tend_for;
%             if del_es > 1
%                 stateend_for = end_set(1 : state_dim, espair(1));
%                 statestart_lat = start_set(1 : state_dim, espair(2));
%                 his_new(:, (tend_for + 1) : (tstart_lat - 1)) = ((del_es - 1) :...
%                     -1 : 1) / del_es .* stateend_for + (1 : (del_es - 1)) / del_es...
%                     .* statestart_lat;
%             end
%             his_new(:, tstart_lat : end) = his_lat(:, tstart_lat : end);
%         end
%         history_copy(:, :, espair(2)) = his_new;
%     end
%     history_copy(:, :, esmap_mat(:, 1)) = [];
% end

% K_del = [];
% for k_idx = 1 : K_est
%     testvec_k = squeeze(history_copy(1, :, k_idx));
%     if sum(isnan(testvec_k)) > sample_time * 0.8
%         K_del = [K_del, k_idx];
%     end
% end
% history_copy(:, :, K_del) = [];
% K_est = size(history_copy, 3);




% gepur_mat = zeros(K_targets, K_est);

% for k_idx = 1 : K_targets
%     for k_jdx = 1 : K_est
%         G_mati = ground_all(:, :, k_idx);
%         E_matj = history_copy(:, :, k_jdx);
%         E_matj(isnan(E_matj)) = 0;
%         GEdel_mat = G_mati - E_matj;
%         normvec = sqrt(sum(GEdel_mat .^ 2, 1));
%         normvec(normvec > Delta_thr) = Delta_thr;
%         gepur_mat(k_idx, k_jdx) = sum(normvec);
%     end
% end

% [label_vec, Est_idx] = Hungarianmatch(gepur_mat);
% save("gepur_mat.mat", 'gepur_mat')
% % set the labels
% glabel_all = reshape(ones(sample_time, 1) * (1 : K_targets), [1, ...
%     sample_time, K_targets]);
% ground_all = [ground_all; glabel_all];
% elabel_all = zeros(1, sample_time, K_est);
% for k_jdx = 1 : K_est
%     if sum(Est_idx == k_jdx) > 0
%         label_k = label_vec(Est_idx == k_jdx);
%     else
%         label_k = -1;
%     end
%     elabel_all(1, :, k_jdx) = label_k * ones(sample_time, 1);
% end
% history_copy = cat(1, history_copy, elabel_all);

% % glabel_t = 1 : m_tr;
% % Gmat_t = [Gmat_t; glabel_t];

% % calculate OSPA
% C_thr = 10;
% OSPA_dist = zeros(sample_time, 1);
% for t_time = 1 : sample_time
%     % obtain the ground truth matrix
%     Gmat_t = squeeze(ground_all(:, t_time, :));
%     m_tr = size(Gmat_t, 2);

%     % obtain the Eestimation matrix
%     Emat_t = squeeze(history_copy(:, t_time, :));
%     [~, nan_idx] = find(isnan(Emat_t(1, :)));
%     Emat_t(:, nan_idx) = [];
%     n_est = size(Emat_t, 2);

%     % calculate the difference of ground truth and estimation
%     GEmat_t = zeros(m_tr, n_est);
%     % GEmat_t = sum((Gmat_t - Emat_t) .^ 2, 1);
%     for m_idx = 1 : m_tr
%         for n_idx = 1 : n_est
%             gvec_mt = Gmat_t(1 : state_dim, m_idx);
%             evec_nt = Emat_t(1 : state_dim, n_idx);
%             glabel_mt = Gmat_t(state_dim + 1, m_idx);
%             elabel_nt = Emat_t(state_dim + 1, n_idx);
%             egdiff_nm = norm(gvec_mt - evec_nt)^2 + (glabel_mt ~= elabel_nt)^2;
%             % if sum(Est_idx == n_idx) > 0
%             %     label_n = label_vec(Est_idx == n_idx);
%             % else
%             %     label_n = 0;
%             % end
%             % if n_idx <= m_tr
%             %     egdiff_nm = norm(gvec_mt - evec_nt)^2 + (m_idx ~= label_n)^2;
%             % else
%             %     egdiff_nm = norm(gvec_mt - evec_nt)^2 + 1;
%             % end
%             GEmat_t(m_idx, n_idx) = egdiff_nm;
%         end
%     end
%     [gt_label, et_label] = Hungarianmatch(GEmat_t);
%     GEvec_t = zeros(min(m_tr, n_est), 1);
%     for nm_idx = 1 : min(m_tr, n_est)
%         GEvec_t(nm_idx) = GEmat_t(gt_label(nm_idx), et_label(nm_idx));
%     end
%     GEvec_t(GEvec_t > C_thr ^ 2) = C_thr ^ 2;
%     OSPA_dist(t_time) = sqrt((sum(GEvec_t) + abs(m_tr - n_est) * (C_thr ^ 2))...
%         / max(m_tr, n_est));
%     % GEmat_t()
% end


% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;
xlim_l = -11;
xlim_r = 11;
ylim_d = 0;
ylim_u = 39;
xlim_vec = xlim_l : xlim_r;
ybon_vec = abs(xlim_vec / 2);


fhandle_fig1 = figure(1);
xlabel('x/m', 'Fontsize', Fsz);
ylabel('y/m', 'Fontsize', Fsz);
xlim([xlim_l xlim_r])
ylim([ylim_d ylim_u])

hold on;

filename = "testmp4.mp4";
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
        phandle_ture = plot(ground_all(1, 1 : t_time, k_tidx), ground_all(3, ...
            1 : t_time, k_tidx), '-bo', 'LineWidth', Lw, 'Marker', 's', ...
            'MarkerSize', Msz);
    end

    for k_idx = 1 : K_est
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




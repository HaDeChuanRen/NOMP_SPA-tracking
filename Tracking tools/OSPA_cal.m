function [OSPA_dist, colabel_vec] = OSPA_cal(ground_all, history_all, ...
    label_all, Delta_thr, C_thr, a_pun)
% OSPA_cal: calculate the optimal subpattern assignment for the tracking result
% Input:
% Output:

% Last update: 2023/7/27
% 2023/7/27 update details
% 1. (Todo) write the function explanations and the comments
% 2. (Todo) use the label of targets rewrite the calculation of OSPA

if ~exist('Delta_thr', 'var'), Delta_thr = 10;
elseif isempty(Delta_thr), Delta_thr = 10; end

if ~exist('C_thr', 'var'), C_thr = 10;
elseif isempty(C_thr), C_thr = 10; end

if ~exist('a_pun', 'var'), a_pun = 1;
elseif isempty(a_pun), a_pun = 1; end

% obtain the state dimension, sample time and number of targets from the ground
% truth tensor
[state_dim, sample_time, K_targets] = size(ground_all);
% obtian the estimated number of targets from the tracking estimation tensor
K_est = size(history_all, 3);
% history_est = history_all;

% calculate the distance between each estimated trajectory and the truth
% trajectory
gepur_mat = zeros(K_est, K_targets);
for k_idx = 1 : K_targets
    for k_jdx = 1 : K_est
        G_mati = ground_all(:, :, k_idx);
        E_matj = history_all(:, :, k_jdx);
        % E_matj(isnan(E_matj)) = 0;
        % Q: how to deal with the nan elements?
        GEdel_mat = G_mati - E_matj;
        normvec = sqrt(sum(GEdel_mat .^ 2, 1));
        normvec(normvec > C_thr) = C_thr;
        normvec(isnan(normvec)) = C_thr;
        gepur_mat(k_jdx, k_idx) = sum(normvec);
    end
end

[Est_idx, target_ind] = Hungarianmatch(gepur_mat);
colabel_vec = nan(K_targets, 1);
for t_idx = 1 : length(target_ind)
    colabel_vec(t_idx) = label_all(Est_idx(t_idx));
end

% set the labels
glabel_all = reshape(ones(sample_time, 1) * (1 : K_targets), [1, ...
    sample_time, K_targets]);
ground_est = [ground_all; glabel_all];
elabel_all = reshape(ones(sample_time, 1) * label_all', ...
    [1, sample_time, K_est]);
history_est = [history_all; elabel_all];
% history_est = cat(1, history_all, elabel_all');
% elabel_all = zeros(1, sample_time, K_est);
% for k_jdx = 1 : K_est
%     if sum(Est_idx == k_jdx) > 0
%         label_k = label_vec(Est_idx == k_jdx);
%     else
%         label_k = -1;
%     end
%     elabel_all(1, :, k_jdx) = label_k * ones(sample_time, 1);
% end


% glabel_t = 1 : m_tr;
% Gmat_t = [Gmat_t; glabel_t];

% calculate OSPA
% C_thr = 10;
OSPA_dist = zeros(sample_time, 1);
for t_time = 1 : sample_time
    % obtain the ground truth matrix
    Gmat_t = squeeze(ground_est(:, t_time, :));
    k_tr = size(Gmat_t, 2);
    inv_kidx = squeeze(isnan(ground_est(1, t_time, :)));

    % obtain the Eestimation matrix
    Emat_t = squeeze(history_est(:, t_time, :));
    n_est = size(Emat_t, 2);
    inv_nidx = squeeze(isnan(history_est(1, t_time, :)));

    % calculate the difference of ground truth and estimation
    EGmat_t = zeros(n_est, k_tr);
    % GEmat_t = sum((Gmat_t - Emat_t) .^ 2, 1);
    for k_idx = 1 : size(Gmat_t, 2)
        for n_idx = 1 : size(Emat_t, 2)
            gvec_mt = Gmat_t(1 : state_dim, k_idx);
            evec_nt = Emat_t(1 : state_dim, n_idx);
            glabel_mt = Gmat_t(state_dim + 1, k_idx);
            elabel_nt = Emat_t(state_dim + 1, n_idx);
            egdiff_nm = norm(gvec_mt - evec_nt) ^ 2 + ...
                a_pun * (colabel_vec(glabel_mt) ~= elabel_nt) ^ 2;
            if egdiff_nm > Delta_thr ^ 2 || isnan(egdiff_nm)
                egdiff_nm = Delta_thr ^ 2;
            end
            EGmat_t(n_idx, k_idx) = egdiff_nm;
            % if sum(Est_idx == n_idx) > 0
            %     label_n = label_vec(Est_idx == n_idx);
            % else
            %     label_n = 0;
            % end
            % if n_idx <= m_tr
            %     egdiff_nm = norm(gvec_mt - evec_nt)^2 + (m_idx ~= label_n)^2;
            % else
            %     egdiff_nm = norm(gvec_mt - evec_nt)^2 + 1;
            % end
        end
    end
    EGmat_t(:, inv_kidx) = [];
    EGmat_t(inv_nidx, :) = [];
    [n_est, k_tr] = size(EGmat_t);
    [et_label, gt_label] = Hungarianmatch(EGmat_t);
    GEvec_t = zeros(min(k_tr, n_est), 1);
    for nm_idx = 1 : min(k_tr, n_est)
        GEvec_t(nm_idx) = EGmat_t(et_label(nm_idx), gt_label(nm_idx) );
    end
    GEvec_t(GEvec_t > Delta_thr ^ 2) = Delta_thr ^ 2;
    OSPA_dist(t_time) = sqrt((sum(GEvec_t) + abs(k_tr - n_est) * (Delta_thr ...
        ^ 2)) / max(k_tr, n_est));
    % GEmat_t()
end
end
function [OSPA_dist, label_vec] = OSPA_cal(ground_all, history_est, Delta_thr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[state_dim, sample_time, K_targets] = size(ground_all);
K_est = size(history_est, 3);

gepur_mat = zeros(K_targets, K_est);

for k_idx = 1 : K_targets
    for k_jdx = 1 : K_est
        G_mati = ground_all(:, :, k_idx);
        E_matj = history_est(:, :, k_jdx);
        E_matj(isnan(E_matj)) = 0;
        GEdel_mat = G_mati - E_matj;
        normvec = sqrt(sum(GEdel_mat .^ 2, 1));
        normvec(normvec > Delta_thr) = Delta_thr;
        gepur_mat(k_idx, k_jdx) = sum(normvec);
    end
end

[label_vec, Est_idx] = Hungarianmatch(gepur_mat);
% set the labels
glabel_all = reshape(ones(sample_time, 1) * (1 : K_targets), [1, ...
    sample_time, K_targets]);
ground_all = [ground_all; glabel_all];
elabel_all = zeros(1, sample_time, K_est);
for k_jdx = 1 : K_est
    if sum(Est_idx == k_jdx) > 0
        label_k = label_vec(Est_idx == k_jdx);
    else
        label_k = -1;
    end
    elabel_all(1, :, k_jdx) = label_k * ones(sample_time, 1);
end
history_est = cat(1, history_est, elabel_all);

% glabel_t = 1 : m_tr;
% Gmat_t = [Gmat_t; glabel_t];

% calculate OSPA
C_thr = 10;
OSPA_dist = zeros(sample_time, 1);
for t_time = 1 : sample_time
    % obtain the ground truth matrix
    Gmat_t = squeeze(ground_all(:, t_time, :));
    m_tr = size(Gmat_t, 2);

    % obtain the Eestimation matrix
    Emat_t = squeeze(history_est(:, t_time, :));
    [~, nan_idx] = find(isnan(Emat_t(1, :)));
    Emat_t(:, nan_idx) = [];
    n_est = size(Emat_t, 2);

    % calculate the difference of ground truth and estimation
    GEmat_t = zeros(m_tr, n_est);
    % GEmat_t = sum((Gmat_t - Emat_t) .^ 2, 1);
    for m_idx = 1 : m_tr
        for n_idx = 1 : n_est
            gvec_mt = Gmat_t(1 : state_dim, m_idx);
            evec_nt = Emat_t(1 : state_dim, n_idx);
            glabel_mt = Gmat_t(state_dim + 1, m_idx);
            elabel_nt = Emat_t(state_dim + 1, n_idx);
            egdiff_nm = norm(gvec_mt - evec_nt)^2 + (glabel_mt ~= elabel_nt)^2;
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
            GEmat_t(m_idx, n_idx) = egdiff_nm;
        end
    end
    [gt_label, et_label] = Hungarianmatch(GEmat_t);
    GEvec_t = zeros(min(m_tr, n_est), 1);
    for nm_idx = 1 : min(m_tr, n_est)
        GEvec_t(nm_idx) = GEmat_t(gt_label(nm_idx), et_label(nm_idx));
    end
    GEvec_t(GEvec_t > C_thr ^ 2) = C_thr ^ 2;
    OSPA_dist(t_time) = sqrt((sum(GEvec_t) + abs(m_tr - n_est) * (C_thr ^ 2))...
        / max(m_tr, n_est));
    % GEmat_t()
end
end
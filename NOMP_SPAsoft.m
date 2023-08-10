function [history_all, label_all, Zcol_cell] = NOMP_SPAsoft(yten_spec, ...
    alpha_2D, N_r, K_max, rmax, vmax, Pmat_cov, Qmat_noise, Amat_state, ...
    Hmat_mea, Rmat_noise)
% last update: 2023/8/9
% 2023/8/9 update:
% 1. (Todo) revise and add comments
% 2. (Todo) use the Smat and calculate the Mahalanobis distance rather than
%   Euclidean distance


[~, ~, ~, sample_time] = size(yten_spec);
T_ter = 1;
state_dim = 4;
mu_c = 0.1;
% Pd_xki = 0.99;
L_ite = 10;    % the iteration times
P_G = 0.9;

Zcol_cell = cell(sample_time, 1);
% apply the tracking algorithm to the noisy line spectrums
% initialize the sequence of target state estimation
Kseq_hat = zeros(sample_time, 1);
history_all = [];
label_all = [];
target_cell = [];

yten_1 = squeeze(yten_spec(:, :, :, 1));
Zarray_lse = LSE_NOMPanalyse(yten_1, alpha_2D, N_r, K_max, rmax, vmax);
Zarray_1 = meas_cluster(Zarray_lse);
% Zarray_1 = Zarray_lse;
Mhat_1 = size(Zarray_1, 1);
Zcol_cell{1} = Zarray_lse;

label_k = 1;
if Mhat_1 > 0
    for m_idx = 1 : Mhat_1
        x_m1 = Zarray_1(m_idx, 1);
        y_m1 = Zarray_1(m_idx, 2);
        statevec_k = [x_m1, 0, y_m1, 0]';
        % Todo: use the velocity estimation to set the original state.
        % and following part
        targetobj_k = TargetState(statevec_k, label_k, Pmat_cov, Qmat_noise, ...
            Amat_state, 1, T_ter);
        label_k = label_k + 1;
        target_cell = [target_cell; targetobj_k];
    end
end
Kseq_hat(1) = size(target_cell, 1);

for t_time = 2 : sample_time
    % whandle = waitbar((t_time - 1) / sample_time);

    % if t_time == 22
    %     1;
    % end

    yten_t = squeeze(yten_spec(:, :, :, t_time));
    % 2D-NOMP analysis for each antenna
    Zarray_lse = LSE_NOMPanalyse(yten_t, alpha_2D, N_r, K_max, rmax, vmax);
    % Zarray_t = Zarray_lse;
    Zarray_t = meas_cluster(Zarray_lse);
    Mhat_t = size(Zarray_t, 1);
    Zcol_cell{t_time} = Zarray_lse;

    Khat_last = Kseq_hat(t_time - 1);
    if Mhat_t == 0 && Khat_last == 0
        continue;
    elseif Mhat_t > 0 && Khat_last == 0
        for m_idx = 1 : Mhat_t
            x_mt = Zarray_t(m_idx, 1);
            y_mt = Zarray_t(m_idx, 2);
            % theta_mt = atan(x_mt / y_mt);
            % vrel_mt = Zarray_t(m_idx, 3);
            % vx_mt = vrel_mt * sin(theta_mt);
            % vy_mt = vrel_mt * cos(theta_mt);
            statevec_k = [x_mt, 0, y_mt, 0]';
            % statevec_k = [Zarray_t(m_idx, 1), 0, Zarray_t(m_idx, 2), 0]';
            targetobj_k = TargetState(statevec_k, label_k, Pmat_cov, Qmat_noise, ...
                Amat_state, 1, T_ter);
            label_k = label_k + 1;
            target_cell = [target_cell; targetobj_k];
        end
        Kseq_hat(t_time) = size(target_cell, 1);
    elseif Mhat_t == 0 && Khat_last > 0
        disapp_set = [];
        for k_idx = 1 : Khat_last
            [~, ~, target_cell(k_idx)] = target_cell(k_idx).state_transform();
            target_cell(k_idx) = target_cell(k_idx).disappearing();
            if target_cell(k_idx).app_state == 0
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
            label_all = [label_all; target_cell(d_idx).n_label];
        end
        target_cell(disapp_set) = [];
        if isempty(target_cell)
            Kseq_hat(t_time) = 0;
        else
            Kseq_hat(t_time) = size(target_cell, 1);
        end
    elseif Mhat_t > 0 && Khat_last > 0
        Zxy_t = Zarray_t(:, 1 : 2);
        % Zvrad_t = Zarray_t(:, 3);
        % Todo: SPA data association and Kalman filter
        XPri_t = zeros(state_dim, Khat_last);
        for k_idx = 1 : Khat_last
            [XPri_kt, Pmat_pri, target_cell(k_idx)] = ...
                target_cell(k_idx).state_transform();
            XPri_t(:, k_idx) = XPri_kt;
        end
        % Smat_t = Hmat_mea * Pmat_pri * Hmat_mea' + Rmat_noise;
        % Q: the definition of Smat should be discussed.

        % initialize the \beta matrix
        betamat_imk = zeros(Mhat_t + 1, Khat_last);
        ximat_mkt = [ones(Mhat_t, Khat_last), mu_c / Khat_last * ...
            ones(Mhat_t, 1)];
        for k_tidx = 1 : Khat_last
            Pd_xki = target_cell(k_tidx).P_det;
            % betamat_imk(1 + Mhat_t, k_tidx) = 1 - Pd_xki;

            % distance calculate and P_G verification
            ZPri_kt = Hmat_mea * XPri_t(:, k_tidx);
            evec_Zkt = ZPri_kt - Zxy_t';
            dvec_kt = diag(evec_Zkt' * evec_Zkt);
            % dvec_kt = diag(evec_Zkt' / Smat_t * evec_Zkt) / (2 * pi * ...
            %     sqrt(det(Smat_t)));
            Pdvec_zmt = exp(- 1/2 * dvec_kt);
            Pdvec_zmt(Pdvec_zmt < P_G) = 0;

            % radial velocity verification
            % theta_prik = atan(XPri_t(1, k_tidx) / XPri_t(3, k_tidx));
            % vrad_prik = XPri_t(2, k_tidx) * sin(theta_prik) + ...
            %     XPri_t(4, k_tidx) * cos(theta_prik);
            % Pdvec_zmt(abs(Zvrad_t - vrad_prik) > vdel) = 0;

            % unify the \beta
            % Pdvec_zmt = Pdvec_zmt / sum(Pdvec_zmt);
            betamat_imk(1 : Mhat_t, k_tidx) = Pd_xki * Pdvec_zmt;
            if sum(Pdvec_zmt) > 0
                betamat_imk(1 : Mhat_t, k_tidx) = Pd_xki * Pdvec_zmt / ...
                    sum(Pdvec_zmt);
                betamat_imk(1 + Mhat_t, k_tidx) = 1 - Pd_xki;
            else
                betamat_imk(1 : Mhat_t, k_tidx) = zeros(Mhat_t, 1);
                betamat_imk(1 + Mhat_t, k_tidx) = 1;
            end
            ximat_mkt(1 : Mhat_t, k_tidx) = Pdvec_zmt;
            % here is some problem: when the target is missed by the measurements,
            % how to express the \beta ?
        end

        phimat_itel = betamat_imk(1 : Mhat_t, :) ./ ...
            (betamat_imk(1 + Mhat_t, :) + eps);
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
        Promat_b(:, Khat_last + 1) = ximat_mkt(:, Khat_last + 1) ./ ...
            (ximat_mkt(:, Khat_last + 1) + sum(ximat_mkt(:, 1 : Khat_last) ...
            .* phimat_itel, 2));


        % Kalman filter
        disapp_set = [];
        akt_coll = [];
        for k_idx = 1 : Khat_last
            prok_vec = Promat_a(:, k_idx);
            target_cell(k_idx) = target_cell(k_idx).est_soft(Rmat_noise, ...
                Hmat_mea, Zxy_t, prok_vec);
            if target_cell(k_idx).app_state == 0
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
            label_all = [label_all; target_cell(d_idx).n_label];
        end
        target_cell(disapp_set) = [];
        % create new target objects
        % newZ_ind = setdiff(1 : Mhat_t, akt_coll);
        newZ_ind = find(Promat_b(:, Khat_last + 1) > 0.5);
        for new_idx = 1 : length(newZ_ind)
            m_idx = newZ_ind(new_idx);
            x_mt = Zarray_t(m_idx, 1);
            y_mt = Zarray_t(m_idx, 2);
%             theta_mt = atan(x_mt / y_mt);
%             vrel_mt = Zarray_t(m_idx, 3);
%             vx_mt = vrel_mt * sin(theta_mt);
%             vy_mt = vrel_mt * cos(theta_mt);
            statevec_k = [x_mt, 0, y_mt, 0]';
            % statevec_k = [Zarray_t(m_idx, 1), 0, Zarray_t(m_idx, 2), 0]';
            targetobj_k = TargetState(statevec_k, label_k, Pmat_cov, Qmat_noise, ...
                Amat_state, 1, T_ter);
            label_k = label_k + 1;
            target_cell = [target_cell; targetobj_k];
        end
        Kseq_hat(t_time) = size(target_cell, 1);
    end
end
% delete(whandle);

if ~isempty(target_cell)
    Khat_last = Kseq_hat(end);
    for k_idx = 1 : Khat_last
        last_k = target_cell(k_idx).last_time;
        app_k = t_time - last_k;
        hisvec_k = target_cell(k_idx).state_his;
        hisall_k = [nan(state_dim, app_k), hisvec_k];
        history_all = cat(3, history_all, hisall_k);
        label_all = [label_all; target_cell(k_idx).n_label];
    end
end

end


% 2023/7/20 update:
% 1. (Todo) revise and add comments
% 2. revise the NOMP function to fit the class defination and test the
%  class definiton file (label and P_det)
% 3. output the labels of targets
% 4. add the clustering algorithm

% 2023/7/24 update:
% 1. (Todo) revise and add comments
% 2. use the appear state flag to decide which target has disappeared.
% 3. estimate the posterier target location by weighting the
%  measurements with the correlation probability, write the new file as
%  MOMP_SPAsoft.m and store the old file as MOMP_SPAhard.m
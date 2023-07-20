function [history_all, Zcol_cell] = NOMP_SPAtracking(yten_spec, alpha_2D, ...
    N_r, K_max, rmax, vmax, Pmat_cov, Qmat_noise, Amat_state, Hmat_mea, ...
    Rmat_noise)
% last update: 2023/7/20
% 2023/7/20 update:
% 1. (Todo) revise and add comments
% 2. (Todo) revise the NOMP function to fit the class defination and test the
%  class definiton file (label and P_det)
% 3. (Todo) output the labels of targets
% 4. (Todo) add the clustering algorithm
%
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[~, ~, ~, sample_time] = size(yten_spec);
T_ter = 1;
state_dim = 4;
mu_c = 1;
Pd_xki = 0.99;
L_ite = 5;    % the iteration times
P_G = 0.3;

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
%         theta_m1 = atan(x_m1 / y_m1);
%         vrel_m1 = Zarray_1(m_idx, 3);
%         vx_m1 = vrel_m1 * sin(theta_m1);
%         vy_m1 = vrel_m1 * cos(theta_m1);
        statevec_k = [x_m1, 0, y_m1, 0]';
        % Todo: use the velocity estimation to set the original state.
        % and following part
        targetobj_k = TargetState(statevec_k, Pmat_cov, Qmat_noise, ...
            Amat_state, 1, T_ter);
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

for t_time = 2 : sample_time
    whandle = waitbar((t_time - 1) / sample_time);

    yten_t = squeeze(yten_spec(:, :, :, t_time));
    % 2D-NOMP analysis for each antenna
    Zarray_t = LSE_NOMPanalyse(yten_t, alpha_2D, N_r, K_max, rmax, vmax);
    Mhat_t = size(Zarray_t, 1);
%     dbscan_indt = dbscan(Zarray_lse, 2, 1);
%     [~, Mhat_t] = max(dbscan_indt);
%     Zarray_t = zeros(Mhat_t, 3);
%     for m_idx = 1 : Mhat_t
%         Zarray_tm = Zarray_lse(dbscan_indt == m_idx, :);
%         Zarray_t(m_idx, :) = mean(Zarray_tm, 1);
%     end
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

            % distance calculate and P_G verification
            ZPri_kt = Hmat_mea * XPri_t(:, k_tidx);
            evec_Zkt = ZPri_kt - Zxy_t';
            dvec_kt = diag(evec_Zkt' / Smat_t * evec_Zkt);
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
        % Mmat_a = zeros(Mhat_t + 1, Khat_last);
        % Mmat_b = zeros(Mhat_t, Khat_last + 1);
        % for k_idx = 1 : Khat_last
        %     [~, aind_kt] = max(Promat_a(:, k_idx));
        %     Mmat_a(aind_kt, k_idx) = 1;
        % end
        % for m_idx = 1 : Mhat_t
        %     [~, bind_mt] = max(Promat_b(m_idx, :));
        %     Mmat_b(m_idx, bind_mt) = 1;
        % end

        % Mmat_ab = zeros(Mhat_t + 1, Khat_last + 1);
        % if sum(xor(Mmat_a(1 : Mhat_t, :), Mmat_b(:, 1 : Khat_last)), 'all')
        %     if sum(Mmat_a(1 : Mhat_t, :)) < Mmat_b(:, 1 : Khat_last)
        %     end
        % end


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
%             theta_mt = atan(x_mt / y_mt);
%             vrel_mt = Zarray_t(m_idx, 3);
%             vx_mt = vrel_mt * sin(theta_mt);
%             vy_mt = vrel_mt * cos(theta_mt);
            statevec_k = [x_mt, 0, y_mt, 0]';
            % statevec_k = [Zarray_t(m_idx, 1), 0, Zarray_t(m_idx, 2), 0]';
            targetobj_k = TargetState(statevec_k, Pmat_cov, Qmat_noise, ...
                Amat_state, 1, T_ter);
            target_cell = [target_cell; targetobj_k];
        end
        Kseq_hat(t_time) = size(target_cell, 1);
    end
end
delete(whandle);

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

end
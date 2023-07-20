function [Promat_a, Promat_b] = SPA_DA(XPri_t, Zxy_t, Hmat_mea, Smat_t, ...
    Pd_xki, mu_c, L_ite)
%SPA_DA Summary of this function goes here
%   Detailed explanation goes here
% match the measurements and the targets by SPA algorithm
    Khat_t = size(XPri_t, 2);
    Mhat_t = size(Zxy_t, 1);
    % initialize the \beta matrix
    betamat_imk = zeros(Mhat_t + 1, Khat_t);

    for k_tidx = 1 : Khat_t
        ZPri_kt = Hmat_mea * XPri_t(:, k_tidx);
        evec_Zkt = ZPri_kt' - Zxy_t;
        dvec_kt = diag(evec_Zkt / Smat_t * evec_Zkt');
        betamat_imk(1 : Mhat_t, k_tidx) = Pd_xki * (exp(- 1/2 * dvec_kt) / ...
            sum(exp(- 1/2 * dvec_kt) + eps(0)));
        % here is some problem: when the target is missed by the measurements,
        % how to express the \beta ?
        betamat_imk(1 + Mhat_t, k_tidx) = 1 - Pd_xki;
    end

    ximat_mkt = [ones(Mhat_t, Khat_t), mu_c / Khat_t * ones(Mhat_t, 1)];
    phimat_itel = betamat_imk(1 : Mhat_t, :) ./ betamat_imk(1 + Mhat_t, :);
    vmat_itel = zeros(Mhat_t, Khat_t);

    % SPA algorithm iterate
    for l_iter = 1 : L_ite
        for m_midx = 1 : Mhat_t
            for i_tidx = 1 : Khat_t
                vmat_itel(m_midx, i_tidx) = 1 / (mu_c / Khat_t + ...
                sum(phimat_itel(m_midx, :)) - phimat_itel(m_midx, i_tidx) * ...
                ximat_mkt(m_midx, i_tidx));
            end
        end

        for i_tidx = 1 : Khat_t
            for m_midx = 1 : Mhat_t
                phimat_itel(m_midx, i_tidx) = betamat_imk(m_midx, i_tidx)/...
                    (betamat_imk(Mhat_t + 1, i_tidx) + ...
                    sum(betamat_imk(1 : Mhat_t, i_tidx) .* ...
                    vmat_itel(:, i_tidx)) - betamat_imk(m_midx, i_tidx) * ...
                    vmat_itel(m_midx, i_tidx));
            end
        end
    end

    Promat_a = zeros(Mhat_t + 1, Khat_t);
    Promat_a(1 : Mhat_t, :) = betamat_imk(1 : Mhat_t, :) .* vmat_itel ./ ...
        (betamat_imk(Mhat_t + 1, :) + sum(betamat_imk(1 : Mhat_t, :)...
        .* vmat_itel));
    Promat_a(Mhat_t + 1, :) = betamat_imk(Mhat_t + 1, :) ./ ...
        (betamat_imk(Mhat_t + 1, :) + sum(betamat_imk(1 : Mhat_t, :)...
        .* vmat_itel));

    Promat_b = zeros(Mhat_t, Khat_t + 1);
    Promat_b(:, 1 : Khat_t) = ximat_mkt(:, 1 : Khat_t) .* phimat_itel ./...
        (ximat_mkt(:, Khat_t + 1) + sum(ximat_mkt(:, 1 : Khat_t) .* ...
        phimat_itel));
    Promat_b(:, Khat_t + 1) = ximat_mkt(:, Khat_t + 1) ./ (ximat_mkt(:, ...
        Khat_t + 1) + sum(ximat_mkt(:, 1 : Khat_t) .* phimat_itel, 2));
end
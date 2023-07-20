function Zarray_t = LSE_NOMPanalyse(yten_t, alpha_2D, N_r, K_max, rmax, vmax)

% 2D-NOMP analysis for each antenna
% yten_t = squeeze(yten_spec(:, :, :, t_time));

[Nx, My, Lz] = size(yten_t);
NM = Nx * My;
sigma_n = 1;

gamma_l = 128;
[omegaxy_t, ghat_t, ~, ~] = NOMP2D_low(yten_t, alpha_2D, N_r, K_max);
% tau = sigma_n * chi2inv((1 - 0.001) ^ (1 / NM), 2 * Lz) / 2;
% [omegaxy_t1, ghat_t1, ~] = MNOMP2D(yten_t, tau);
[Mhat_t, ~] = size(omegaxy_t);

if Mhat_t > 0
    Zarray_t = zeros(Mhat_t, 3);
    rhat_mt = wrapTo2Pi(omegaxy_t(:, 1)) / (2 * pi) * rmax;
    % detect gain and frequency of an additional sinusoid
    S = eye(Lz);
    sampledManifold = preProcessMeasMat(S, gamma_l);
    omegazhat_t = zeros(Mhat_t, 1);
    R_s = 1;
    for m_idx = 1 : Mhat_t
        y_r = ghat_t(m_idx, :)';
        [omegazhat_mt, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);
        for i = 1 : R_s
            [omegazhat_mt, ~, y_r] = refineOne(y_r, omegazhat_mt, ...
                gain_new, S, sampledManifold.ant_idx, true);
        end
        omegazhat_t(m_idx) = wrapToPi(omegazhat_mt);
    end
    % [~, z_idx] = max(abs(fft(ghat_t, gamma_l, 2)), [], 2);
    % omegazhat_t = wrapToPi(((z_idx - 1) / gamma_l) * 2 * pi);

    theta_mt = asin(- omegazhat_t / pi);
    Zarray_t(:, 1) = rhat_mt .* sin(theta_mt);
    Zarray_t(:, 2) = rhat_mt .* cos(theta_mt);
    Zarray_t(:, 3) = (wrapToPi(omegaxy_t(:, 2)) / pi) * vmax;
else
    Zarray_t = [];
end

% omegaxycell_t = cell(Lz, 1);
% ghatcell_t = cell(Lz, 1);
% for l_idx = 1 : Lz
%     ymat_lt = squeeze(yten_t(:, :, l_idx));
%     [omegahat_lt, ghat_lt, ~, ~] = NOMP2D_low(ymat_lt, alpha_2D, N_r, K_max);
%     omegaxycell_t{l_idx} = wrapToPi(omegahat_lt);
%     ghatcell_t{l_idx} = ghat_lt;
% end

% match the measurements of each antenna
% set the estimate number of measurements as the largest number of
% measurements from each antenna
% Mhat_t = 0;
% for l_idx = 1 : Lz
%     if Mhat_t < length(ghatcell_t{l_idx})
%         Mhat_t = length(ghatcell_t{l_idx});
%     end
% end

% for l_idx = 1 : Lz
%     if size(omegaxycell_t{l_idx}, 1) < Mhat_t
%         omegaxycell_t{l_idx} = [omegaxycell_t{l_idx}; zeros(Mhat_t - ...
%             size(omegaxycell_t{l_idx}, 1), 2)];
%         ghatcell_t{l_idx} = [ghatcell_t{l_idx}; zeros(Mhat_t - ...
%             length(ghatcell_t{l_idx}), 1)];
%     end
% end

% % initialize the matching result of each antenna
% gomegacell_K = zeros(Mhat_t, Lz, 3);
% % the sequnce of the 3 variables is ghat, omegax, omegay
% gomegacell_K(:, 1, 1) = ghatcell_t{1}(:);
% gomegacell_K(:, 1, 2 : 3) = omegaxycell_t{1}(:, :);

% % match the measurements of each antenna by Hungarian match algorithm
% for l_idx = 2 : Lz
%     dist_lt = abs(omegaxycell_t{l_idx}(:, 1) - omegaxycell_t{1}(:, 1)');
%     [mi_Hun, ~] = Hungarianmatch(dist_lt);
%     gomegacell_K(:, l_idx, 1) = ghatcell_t{l_idx}(mi_Hun);
%     gomegacell_K(:, l_idx, 2 : 3) = omegaxycell_t{l_idx}(mi_Hun, :);
% end

% % calculate the result of measurements with the matching results
% gomega_t = zeros(Mhat_t, 4);
% % the sequence of variables: ghat, omegax, omegay, omegaz
% Zarray_t = zeros(Mhat_t, 3);
% % the sequence of variables: xPosition, yPosition, radial velocity

% for m_idx = 1 : Mhat_t
%     % calculate the average of omegax and omegay as the result of $m$th
%     % measurement in the intant $t$
%     omegaxyvec_mt = squeeze(mean(gomegacell_K(m_idx, :, 2 : 3), 2))';
%     gomega_t(m_idx, 2 : 3) = omegaxyvec_mt;

%     % use the average result to adjust the ghat results by LS method
%     gvec_mt = zeros(Lz, 1);
%     for l_idx = 1 : Lz
%         gvec_mt(l_idx) = LeastSquares_2D(squeeze(yten_t(:, :, l_idx)), ...
%             omegaxyvec_mt);
%     end
%     gomega_t(m_idx, 1) = gvec_mt(1);

%     % squeeze(mean(gomegacell_K(m_idx, :, 2 : 3), 2));
%     % gvec_mt = gomegacell_K(m_idx, :, 1);

%     % calculate omega with gvec by fft method
%     [~, z_idx] = max(abs(fft(gvec_mt, gamma_l)));
%     omegazhat_mt = wrapToPi((z_idx / gamma_l) * 2 * pi);
%     gomega_t(m_idx, 4) = omegazhat_mt;
%     rhat_mt = wrapTo2Pi(gomega_t(m_idx, 2)) / (2 * pi) * rmax;
%     theta_mt = asin(omegazhat_mt / pi);
%     Zarray_t(m_idx, 1) = rhat_mt * sin(theta_mt);
%     Zarray_t(m_idx, 2) = rhat_mt * cos(theta_mt);
%     Zarray_t(m_idx, 3) = (gomega_t(m_idx, 3) / pi) * vmax;
% end


end
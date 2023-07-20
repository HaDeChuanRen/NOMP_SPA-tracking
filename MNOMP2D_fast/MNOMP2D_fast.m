function [omegaList, gainList, Ymat_res, Threshold_collect] = ...
    MNOMP2D_fast(Ymat, alpha_set, N_r, K_max, guard_band, ...
    overSamplingRate, R_s, R_c)
%   last update: 2023/5/31
%   input parameter:
%   output parameter:

    if ~exist('guard_band', 'var'), guard_band = [3, 3];
    elseif isempty(guard_band), guard_band = [3, 3]; end

    if ~exist('overSamplingRate', 'var'), overSamplingRate = [8; 8];
    elseif isempty(overSamplingRate), overSamplingRate = [8; 8]; end

    if ~exist('R_s', 'var'), R_s = 1;
    elseif isempty(R_s), R_s = 1; end

    if ~exist('R_c', 'var'), R_c = 3;
    elseif isempty(R_c), R_c = 3; end

    [Nx, My, T] = size(Ymat);
    NM_num = Nx * My;
    antvec_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
    antvec_My = (0 : (My - 1))' - (My - 1) / 2;
    array_Fun = @(omega, N) exp(1j * ((0 : (N - 1))' - (N - 1) / 2) * ...
        omega) / sqrt(N);
    farray_Fun = @(omega, N) exp(- 1j * omega * (N - 1) / 2) * (1 - exp(1j *...
        N * (omega - (0 : (N - 1))' * 2 * pi / N))) ./ (1 - exp(1j * (omega -...
        (0 : (N - 1))' * 2 * pi / N))) / N;

    K_est = 0;
    omegaList = [];
    gainList = [];
    Threshold_collect = [];
    Ymat_res = Ymat;

    flag_cfar = 0;
    omegavec_in = [];
    gainvec_in = [];
    Thrvec_in = [];

    OSR_x = overSamplingRate(1);
    OSR_y = overSamplingRate(2);
    Nx_ext = Nx * OSR_x;
    My_ext = My * OSR_y;
    Y_spec = fft2(Ymat_res, Nx_ext, My_ext) / sqrt(NM_num);
    Yr_spec = Y_spec;
    % sum(Yr_spec - fft2(Ymat_res, Nx_ext, My_ext) / sqrt(NM_num), 'all')
    % imagesc(sum(abs(fft2(Ymat_res, Nx_ext, My_ext) / sqrt(NM_num)), 3))
    % imagesc(sum(abs(Yr_spec), 3))
    while (K_est < K_max) && (flag_cfar < 3)
        [T_judgement, Threshold_CUT] = CFARdet2D_spec(Yr_spec, N_r, ...
            alpha_set, overSamplingRate, omegaList, guard_band);

        [omega_k, gain_k] = Det2D_spec(Yr_spec, overSamplingRate);
        % update Ymat_res
        for t_idx = 1 : T
            Ymat_res(:, :, t_idx) = Ymat_res(:, :, t_idx) - gain_k(t_idx) * ...
                array_Fun(omega_k(1), Nx) * (array_Fun(omega_k(2), My).');
        end
        for refione_idx = 1 : R_s
            [Ymat_res, omega_k, gain_k] = ReOne2D(Ymat_res, omega_k, gain_k);
        end
        for t_idx = 1 : T
            Yr_spec(:, :, t_idx) = Yr_spec(:, :, t_idx) - gain_k(t_idx) * ...
                farray_Fun(omega_k(1), Nx_ext) * ...
                farray_Fun(omega_k(2), My_ext).';
        end

        if T_judgement < 0
            flag_cfar = flag_cfar + 1;
            omegavec_in = [omegavec_in; omega_k];
            gainvec_in = [gainvec_in; gain_k];
            Thrvec_in = [Thrvec_in; Threshold_CUT];
            continue;
        else
            omegaList = [omegaList; omegavec_in; omega_k];
            gainList = [gainList; gainvec_in; gain_k];
            Threshold_collect = [Threshold_collect; Thrvec_in; Threshold_CUT];
            omegavec_in = [];
            gainvec_in = [];
            Thrvec_in = [];
            flag_cfar = 0;
            K_est = K_est + 1;
            [omegaList, gainList, Ymat_res] = ReAll2D(Ymat_res, ...
                omegaList, gainList, R_s, R_c);
            for k_idx = 1 : K_est
                for t_idx = 1 : T
                    Yr_spec(:, :, t_idx) = Y_spec(:, :, t_idx) - ...
                        gainList(k_idx, t_idx) * ...
                        farray_Fun(omegaList(k_idx, 1), Nx_ext) * ...
                        farray_Fun(omegaList(k_idx, 2), My_ext).';
                end
            end
        end
    end

    if ~isempty(omegaList)
        gainList = bsxfun(@times, gainList, exp(1j * (antvec_Nx(1) * ...
        omegaList(:, 1) + antvec_My(1) * omegaList(:, 2))));

        % gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
        omegaList =  wrapTo2Pi(omegaList);
    end

end











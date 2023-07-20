clc; clear; close all;

%% obtain file path and read the ADC data
addpath("analysis tools\")
addpath("NOMP2D tools\")
orginal_path = 'C:\study\202210target tracking\3data\20230710exp';
exp_type = '\tracking';
exp_serial = '\09';

filename = [orginal_path, exp_type, exp_serial, '\adc_data.bin'];
data_cube = readadc(filename);

%% set the observation scene and classify the data cube
Nx = 256;
My = 128;
Lz = 4;
NL_num = Nx * Lz;
NML = Nx * My * Lz;
ant_idx_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
ant_idx_My = (0 : (My - 1))' - (My - 1) / 2;

num_frame = size(data_cube, 2);
sample_time = num_frame / My;
yten_spec = zeros(Nx, My, Lz, sample_time);

for t_time = 1 : sample_time
    yten_spec(:, :, :, t_time) = data_cube(:, ((t_time - 1) * My + ...
        1 : t_time * My), :);
end

state_dim = 4;    % the number of state dimensions
range_c = [-100, -100; 100, 100];    % the range of clutters
T = 0.05;    % sample interval
sigma_v = [0.1; 0.1];    % measurement noise variance
Pd_xki = 0.9999;    % detection probability
L_ite = 20;    % the iteration times
P_G = 0.9;
db_num = 2;
db_dis = 2;


%% radar parameter
c = 3e8;
T_idle = 100e-6;
T_ramp = 60e-6;
T_circle = T_idle + T_ramp;
Fre_start = 77e9;
Slope_fre = 29.982e12;
Fs = 10e6;
Ts = 1 / Fs;
lambda_cw = c / Fre_start;
Rx_interval = lambda_cw / 2;

rmax = c / (2 * Ts * Slope_fre);    % maximum radial range
vmax = c / (4 * Fre_start * T_circle);    % maximum radial velocity

xlim_l = -5;
xlim_r = 5;
ylim_u = 6;
xlim_vec = xlim_l : xlim_r;
ybon_vec = abs(xlim_vec / 2);


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
K_max = 13;
gamma_l = 128;
overSamplingRate = [4; 4];

for t_time = 1 : sample_time
    whandle = waitbar((t_time - 1) / sample_time);
    if t_time ~= 1
        continue;
    end
    yten_t = squeeze(yten_spec(:, :, :, t_time));
    % 2D-NOMP analysis for each antenna
    % Zarray_lse = LSE_NOMPanalyse(yten_t, alpha_2D, N_r, K_max, rmax, vmax);

    % [omegaxy_t, ghat_t, ~, ~] = NOMP2D_CFAR(yten_t, alpha_2D, N_r, K_max);
    % tau = sigma_n * chi2inv((1 - 0.001) ^ (1 / NM), 2 * Lz) / 2;
    % [omegaxy_t1, ghat_t1, ~] = MNOMP2D(yten_t, tau);

    K_est = 0;
    omegaxy_t = zeros(K_max, 2);
    gainList = zeros(K_max, Lz);
    y_residue_matrix = yten_t;
    y_fftabs = prod2D_cal(y_residue_matrix); figure; imagesc(y_fftabs);

    while (K_est < K_max)
        [omega_k, gain_k, ~, y_residue_matrix] = ...
            DetectNew_2D(y_residue_matrix, overSamplingRate);

        for refione_idx = 1 : R_s
            [y_residue_matrix, omega_k, gain_k] = ...
                RefineOne_2D(y_residue_matrix, omega_k, gain_k);
        end

        K_est = K_est + 1;
        omegaxy_t(K_est, :) = omega_k;
        gainList(K_est, :) = gain_k;

        [omega_est, ~, ~] = RefineAll_2D(y_residue_matrix, ...
            omegaxy_t(1 : K_est, :), gainList(1 : K_est, :), R_s, R_c);
        [gain_est, y_residue_matrix, A_all_omega] = ...
            LeastSquares_2D(y_matrix, omega_est);
        omegaxy_t(1 : K_est, :) = omega_est;
        gainList(1 : K_est, :) = gain_est;
        y_fftabs = abs(fft2(y_residue_matrix)); figure; imagesc(y_fftabs);
    end
    [Mhat_t, ~] = size(omegaxy_t);

    if Mhat_t > 0
        Zarray_lse = zeros(Mhat_t, 3);
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
        Zarray_lse(:, 1) = rhat_mt .* sin(theta_mt);
        Zarray_lse(:, 2) = rhat_mt .* cos(theta_mt);
        Zarray_lse(:, 3) = (wrapToPi(omegaxy_t(:, 2)) / pi) * vmax;
    else
        Zarray_lse = [];
    end
end






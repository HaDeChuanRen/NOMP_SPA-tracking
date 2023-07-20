% obeject method test
% date: 2023/4/19

clc; clear; close all;
rng(3);
addpath('NOMP2D tools\')
addpath('analysis tools\')

% basic parameter set
sample_time = 60;   % the total samples of tracking
K_targets = 1;    % the number of targets
mu_c = 10;    % the mean number of clutters (Possion)
state_dim = 4;    % the number of state dimensions
range_c = [-100, -100; 100, 100];    % the range of clutters
T = 1;    % sample interval
sigma_w = [0.1; 0; 0.1; 0];    % state transition noise variance
sigma_v = [0.1; 0.1];    % measurement noise variance
Pd_xki = 0.9999;    % detection probability
L_ite = 20;    % the iteration times

% the transition parameters matrix
Amat_state = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];    % state transition matrix
Hmat_mea = [1 0 0 0; 0 0 1 0];    % measurment matrix
Rmat_noise = eye(2) * diag(sigma_v);    % measurement covriance
% Todo: use CRB to describe the Rmat_noise
Tmat_temp = [1, 1/T; 1/T, 2/T^2];
Pmat_cov = kron(Rmat_noise, Tmat_temp);
% Question: why?
V = wgn(1,2,0);
Qmat_noise = [T^4/4*(V'*V)  T^3/2*(V'*V); T^3/2*(V'*V) T^2*(V'*V)];


% targets location and velocity initialization
xPosition_true0 = -30;
yPosition_true0 = 3;
xVelocity_true0 = 1;
yVelocity_true0 = 1.5;


% create the real state array
xPosseq_true = zeros(K_targets, sample_time);
yPosseq_true = zeros(K_targets, sample_time);
xVelseq_true = zeros(K_targets, sample_time);
yVelseq_true = zeros(K_targets, sample_time);


xPosseq_true(:, 1) = xPosition_true0;
yPosseq_true(:, 1) = yPosition_true0;
xVelseq_true(:, 1) = xVelocity_true0;
yVelseq_true(:, 1) = yVelocity_true0;


for t_time = 2 : sample_time
    Xtrue_kt = Amat_state * [xPosseq_true(:, t_time - 1)'; ...
        xVelseq_true(:, t_time - 1)'; yPosseq_true(:, t_time - 1)'; ...
        yVelseq_true(:, t_time - 1)'] + ...
        sqrt(sigma_w) .* randn(state_dim, 1);
    xPosseq_true(:, t_time) = Xtrue_kt(1, :);
    xVelseq_true(:, t_time) = Xtrue_kt(2, :);
    yPosseq_true(:, t_time) = Xtrue_kt(3, :);
    yVelseq_true(:, t_time) = Xtrue_kt(4, :);
end

% use the state array to create line spectrum in each time instant
% radar parameter
c = 299792458;    % electromagnetic speed
T_idle = 100e-6;    % idle time
T_ramp = 60e-6;    % ramp time
T_circle = T_idle + T_ramp;    % chirp period
Fre_start = 77e9;    % carrier frequency
Slope_fre = 14.991e12;    % chirp rate
Fs = 10e6;    % sample rate
Ts = 1 / Fs;     % sample interval
lambda_cw = c / Fre_start;    % wavelength
Rx_interval = lambda_cw / 2;    % antenna element spacing

rmax = c / (2 * Ts * Slope_fre);
vmax = c / (4 * Fre_start * T_circle);


% line spectrum parameter
Nx = 256;
My = 64;
Lz = 8;
NML = Nx * My * Lz;
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);

snr_set = 30;
sigma_n = 1;
% Todo: the SNR and noise variance can be calculated by equations

% initialize real r, v, theta parameters of each instant
rmat_true = zeros(K_targets, sample_time);
vmat_true = zeros(K_targets, sample_time);
thetamat_true = zeros(K_targets, sample_time);


% initialize real omega parameters of each dimension
% $\omega \in (-\pi, \pi]$
omegaxmat_true = zeros(K_targets, sample_time);
omegaymat_true = zeros(K_targets, sample_time);
omegazmat_true = zeros(K_targets, sample_time);
gain_true = sqrt(10 .^ (snr_set / 10) * sqrt(sigma_n)) .* ...
    exp(1j * 2 * pi * rand(K_targets, sample_time));

yten_spec = zeros(Nx, My, Lz, sample_time);

for t_time = 1 : sample_time
    % obtain the real measurements in the instant t
    rmat_true(:, t_time) = sqrt(xPosseq_true(:, t_time) .^ 2 + ...
        yPosseq_true(:, t_time) .^ 2);
    thetamat_true(:, t_time) = atan(xPosseq_true(:, t_time) ./ ...
        yPosseq_true(:, t_time));
    vmat_true(:, t_time) = xVelseq_true(:, t_time) .* ...
        cos(thetamat_true(:, t_time)) + yVelseq_true(:, t_time) .* ...
        sin(thetamat_true(:, t_time));

    % calculate the corresponding angular frequency
    % $\omega \in (-\pi, \pi]$
    omegaxmat_true(:, t_time) = wrapToPi(rmat_true(:, t_time) / rmax * 2 * pi);
    omegaymat_true(:, t_time) = vmat_true(:, t_time) / vmax * pi;
    omegazmat_true(:, t_time) = pi * sin(thetamat_true(:, t_time));

    % create the line spectrum tensors in the instant t
    yvec_t = zeros(NML, 1);
    for k_idx = 1 : K_targets
        if abs(thetamat_true(k_idx, t_time)) > pi / 3
            continue;
        end
        yvec_kt = kron(kron(array_Fun(omegazmat_true(k_idx, t_time), Lz), ...
            array_Fun(omegaymat_true(k_idx, t_time), My)), ...
            array_Fun(omegaxmat_true(k_idx, t_time), Nx));
        yvec_t = yvec_t + gain_true(k_idx, t_time) * yvec_kt;
    end
    yvec_t = yvec_t + sqrt(1 / 2) * (randn(NML, 1) + 1j * randn(NML, 1));
    yten_spec(:, :, :, t_time) = reshape(yvec_t, Nx, My, Lz);
end

% set the parameters of NOMP-CFAR algorithm
P_fa = 0.01;
NM = Nx * My;
N_r = 60;
alpha_2D = alpha_PoebyS(P_fa, NM, N_r, Lz);
K_max = 16;

% apply the tracking algorithm to the noisy line spectrums
% initialize the sequence of target state estimation
Kseq_hat = zeros(sample_time, 1);
targetseq_hat = cell(sample_time, 1);

yten_1 = squeeze(yten_spec(:, :, :, 1));
Zarray_1 = LSE_NOMPanalyse(yten_1, alpha_2D, N_r, K_max, rmax, vmax);
Mhat_1 = size(Zarray_1, 1);
if Mhat_1 > 0
    statevec_k = [Zarray_1(1, 1), 0, Zarray_1(1, 2), 0]';
    targetobj_hat = TargetState(statevec_k, Pmat_cov, Qmat_noise, ...
        Amat_state, 1, T);
end
Kseq_hat(1) = Mhat_1;


% Kseq_hat(1) = Khat_1;
% for k_idx = 1 : Khat_1
%     statevec_k = [Zarray_1(k_idx, 1), Zarray_1(k_idx, 2), 0, 0]';
%     targethat_k1 = TargetState(statevec_k, Pmat_cov, Qmat_noise, Amat_state, ...
%         1, T);
%     targetseq_hat{1}(end + 1) = targethat_k1;
% end


tic;
for t_time = 2 : sample_time
    whandle = waitbar((t_time - 1) / sample_time);
    yten_t = squeeze(yten_spec(:, :, :, t_time));

    % 2D-NOMP analysis for each antenna
    Zarray_t = LSE_NOMPanalyse(yten_t, alpha_2D, N_r, K_max, rmax, vmax);
    Mhat_t = size(Zarray_t, 1);
    if Mhat_t > Kseq_hat(t_time - 1)
        statevec_k = [Zarray_t(1, 1), 0, Zarray_t(1, 2), 0]';
        targetobj_hat = TargetState(statevec_k, Pmat_cov, Qmat_noise, ...
            Amat_state, t_time, T);
        Kseq_hat(t_time) = Mhat_t;
    elseif Mhat_t > 0 && Kseq_hat(t_time - 1) > 0
        [XPri_t, Pmat_pri, targetobj_hat] = targetobj_hat.state_transform();
        Zxy_t = Zarray_t(:, 1 : 2);
        ZPri_kt = Hmat_mea * XPri_t;
        Smat_t = Hmat_mea * Pmat_pri * Hmat_mea' + Rmat_noise;
        evec_Zkt = ZPri_kt' - Zxy_t;
        dvec_kt = diag(evec_Zkt / Smat_t * evec_Zkt');
        [~, z_idx] = min(dvec_kt);

        Kmat_gain = Pmat_pri * Hmat_mea' / (Hmat_mea * Pmat_pri * ...
            Hmat_mea' + Rmat_noise);
        targetobj_hat = targetobj_hat.state_estimate(Kmat_gain, Hmat_mea, ...
            Zxy_t(z_idx, :));
        Kseq_hat(t_time) = Mhat_t;
    end


end


toc;
delete(whandle);

% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;


fhandle_fig1 = figure(1);
im_seq = cell(sample_time - 1, 1);
appear_ins = sample_time - targetobj_hat.last_time;

for t_time = 1 : sample_time - 1

    phandle_ture = plot(xPosseq_true(:, t_time), yPosseq_true(:, t_time), ...
        'b*', 'LineWidth', Lw, 'Marker', 's', 'MarkerSize', Msz);
    hold on;

    if t_time >= appear_ins
        phandle_hat = plot(targetobj_hat.state_his(1, t_time - ...
            appear_ins + 1), targetobj_hat.state_his(3, t_time - ...
            appear_ins + 1), 'r*', 'LineWidth', Lw, 'Marker', 's', ...
            'MarkerSize', Msz);
    end
    xlabel('x/m', 'Fontsize', Fsz);
    ylabel('y/m', 'Fontsize', Fsz);
    xlim([-120 120])
    ylim([0 120])
    pause(0.1);

    frame_t = getframe(fhandle_fig1);
    im_seq{t_time} = frame2im(frame_t);
end

filename = "testgif.gif";
for t_time = 1 : sample_time - 1
    [A_t, map_t] = rgb2ind(im_seq{t_time}, 256);
    if t_time == 1
        imwrite(A_t, map_t, filename, 'gif', 'LoopCount', Inf, 'DelayTime', ...
            0.1);
    else
        imwrite(A_t, map_t, filename, 'gif', 'WriteMode', 'append', ...
            'DelayTime',0.1);
    end
end



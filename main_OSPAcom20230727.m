% this file tests and compares the properties of tracking algorithms
% Last update: 2023/08/10
% 2023/8/10 update details:
% 1. (Todo) restore the OSPA calculation result


clc; clear; close all;
rng(7);
addpath('NOMP1D tools\')
addpath('NOMP2D tools\')
addpath('analysis tools\')
addpath('Tracking tools')

MC = 3;

% targets location, velocity, appear time disapear time initialization
sample_time = 60;   % the total samples of tracking
K_targets = 8;
xrange = 10;
xmin = -5;
vxrange = 0.2;
vxmin = -0.1;
yrange = 15;
ymin = 6;
vyrange = 0.4;
vymin = -0.1;



xlim_l = xmin + sample_time * min([vxmin, 0]);
xlim_r = xmin + xrange + sample_time * max([vxrange + vxmin, 0]);
ylim_d = ymin + sample_time * min([vymin, 0]);
ylim_u = ymin + yrange + sample_time * max([vyrange + vymin, 0]);
xlim_vec = xlim_l : xlim_r;
ybon_vec = abs(xlim_vec / 2);

% basic parameter set
mu_c = 1;    % the mean number of clutters (Possion)
state_dim = 4;    % the number of state dimensions
range_c = [-100, -100; 100, 100];    % the range of clutters
T = 1;    % sample interval
sigma_w = [0; 0; 0; 0];    % state transition noise variance
sigma_v = [0.001; 0.001];    % measurement noise variance
Pd_xki = 0.9;    % detection probability
L_ite = 5;    % the iteration times

% the transition parameters matrix
Amat_state = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];    % state transition matrix
Hmat_mea = [1 0 0 0; 0 0 1 0];    % measurment matrix
Rmat_noise = eye(2) * diag(sigma_v);    % measurement covriance
% Todo: use CRB to describe the Rmat_noise
Tmat_temp = [1, 1/T; 1/T, 2/T^2];
Pmat_cov = kron(Rmat_noise, Tmat_temp);
% Question: why?
V = wgn(1,2,0);
Qmat_noise = zeros(state_dim);
sigma_n = 1;
% Qmat_noise = [T^4/4*(V'*V)  T^3/2*(V'*V); T^3/2*(V'*V) T^2*(V'*V)];
% Q: revise the definition of Qmat

% line spectrum parameter
Nx = 256;
My = 64;
Lz = 8;
SNR_kt = 25;

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

rmax = c / (2 * Ts * Slope_fre);    % maximum radial range
vmax = c / (4 * Fre_start * T_circle);    % maximum radial velocity

% set the parameters of NOMP-CFAR algorithm
P_fa = 0.001;
N_r = 60;
NM = Nx * My;
alpha_2D = alpha_PoebyS(P_fa, NM, N_r, Lz);
K_max = 10;

OSPAcol_nompsoft = zeros(MC, sample_time);
OSPAcol_nomphard = zeros(MC, sample_time);
tic;
for mc_idx = 1 : MC
    hwaitbar = waitbar((mc_idx - 1) / MC);
    xPosition_true0 = xrange * rand(K_targets, 1) + xmin;
    xVelocity_true0 = vxrange * rand(K_targets, 1) + vxmin;
    yPosition_true0 = yrange * rand(K_targets, 1) + ymin;
    yVelocity_true0 = vyrange * rand(K_targets, 1) + vymin;
    Xori_mat = [xPosition_true0'; xVelocity_true0'; yPosition_true0'; ...
    yVelocity_true0'];
    [ground_all, yten_spec] = grdsig_prod(sample_time, Xori_mat, Amat_state, ...
        sigma_w, sigma_n, rmax, vmax, Nx, My, Lz, SNR_kt);

    [hten_nompsoft, lvec_nompsoft, ~] = NOMP_SPAsoft(yten_spec, alpha_2D, ...
        N_r, K_max, rmax, vmax, Pmat_cov, Qmat_noise, Amat_state, Hmat_mea, ...
        Rmat_noise);
    [OSPA_nompsoft, ~] = OSPA_cal(ground_all, hten_nompsoft, lvec_nompsoft);
    OSPAcol_nompsoft(mc_idx, :) = OSPA_nompsoft';

    [hten_nomphard, lvec_nomphard, ~] = NOMP_SPAhard(yten_spec, alpha_2D, ...
        N_r, K_max, rmax, vmax, Pmat_cov, Qmat_noise, Amat_state, Hmat_mea, ...
        Rmat_noise);
    [OSPA_nomphard, ~] = OSPA_cal(ground_all, hten_nomphard, lvec_nomphard);
    OSPAcol_nomphard(mc_idx, :) = OSPA_nomphard';
end
toc;
delete(hwaitbar)

filename_now = ['mat_save\', datestr(now, 30), 'mc', num2str(MC), '_OSPA.mat'];
save(filename_now, 'sample_time', 'OSPAcol_nompsoft', 'OSPAcol_nomphard');


MOSPA_nompsoft = mean(OSPAcol_nompsoft);
MOSPA_nomphard = mean(OSPAcol_nomphard);



% plot the target trajectory
Fsz = 12;
Lw = 2;
Msz = 2;

figure;
hold on;
plot(1 : sample_time, MOSPA_nomphard, '-b', 'LineWidth', Lw, 'Marker', 's', ...
    'MarkerSize', Msz)
plot(1 : sample_time, MOSPA_nompsoft, '-g', 'LineWidth', Lw, 'Marker', 's', ...
    'MarkerSize', Msz)
legend('NOMP-SPA (hard)', 'NOMP-SPA (soft)')
xlabel('time(s)', 'Fontsize', Fsz)
ylabel('OSPA', 'Fontsize', Fsz)








% 2023/7/27 update details:
% 1. Build a platform to test the tracking algorithms through OSPA,
%  which include: Monte Carlo method, generating target and received signals,
%  tracking algorithm, OSPA calculation.
clc; clear; close all;

r_true1 = 3;
theta_true1 = pi / 6;
z_true1 = [r_true1 * sin(theta_true1), r_true1 * cos(theta_true1)];

p_x = -4 : 0.005 : 4;
p_y = 0 : 0.005 : 8;
length_x = length(p_x);
length_y = length(p_y);

[P_x, P_y] = meshgrid(p_x, p_y);
Pten = zeros(length_y, length_x, 2);
Pten(:, :, 1) = P_x;
Pten(:, :, 2) = P_y;

Pmat_ind = reshape(Pten, [length_x * length_y, 2]);

r_max = 50;
Nx = 256;
Lz = 8;
SNR = 15;

CRB_mat = [(r_max * sin(theta_true1))^2 / (4 * (Nx^2 - 1)) + ...
    r_true1^2 / (Lz^2 - 1), (r_max ^ 2 / (4 * (Nx^2 - 1)) - r_true1^2 / ...
    ((Lz^2 - 1) * cos(theta_true1) ^ 2)) * sin(theta_true1) * ...
    cos(theta_true1); (r_max ^ 2 / (4 * (Nx^2 - 1)) - r_true1 ^ 2 / ...
    ((Lz^2 - 1) * cos(theta_true1) ^ 2)) * sin(theta_true1) * ...
    cos(theta_true1), (r_max * cos(theta_true1))^2 / (4 * (Nx^2 - 1)) + ...
    (r_true1 * tan(theta_true1))^2 / (Lz^2 - 1)];
Rmat_CRB = 6 / (pi^2) * 10 ^ (- SNR / 10) * CRB_mat;

temp_mat = (Pmat_ind-z_true1) / Rmat_CRB;

dis_vec = exp(- sum(temp_mat .* (Pmat_ind-z_true1), 2) / 2);

% dis_res = diag(exp(-(Pmat_ind-z_true1) / Rmat_CRB * (Pmat_ind-z_true1)'/ 2));
dis_mat = reshape(dis_vec, [length_y, length_x]);

contour(P_x, P_y, dis_mat);
colorbar;

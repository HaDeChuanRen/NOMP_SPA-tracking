function prob_ind = prod2D_cal(y_matrix, gamma_oversampling)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('overSamplingRate','var'), gamma_oversampling = [4; 4];
elseif isempty(gamma_oversampling), gamma_oversampling = [4; 4]; end
[Nx, My, ~] = size(y_matrix);
NM_num = Nx * My;

Nx_ext = Nx * gamma_oversampling(1);
My_ext = My * gamma_oversampling(2);

% 2D-fft and find the maximum index
Y_spec = fft2(y_matrix, Nx_ext, My_ext) / sqrt(NM_num); % Nx_ext*My_ext*T
a2D_fft_abs2 = sum(abs(Y_spec).^2,3);
prob_ind = a2D_fft_abs2(1 : gamma_oversampling(1) : end, 1 : ...
    gamma_oversampling(2) : end);
end
function [omega_est, tensor_idx] = detectnew_3D(y_tensor, oversample_rate)
%   this funciton is used to detect the angular frequency of one target and
%   the corresponding index in the tensor
%
%   Input parameter:
%   y_tensor: the three dimension receiving signal been observed              #Nx * My * Lz
%   oversample_rate: oversample rate                                        #3 * 1
%
%   Output parameter:
%   omega_est: the coarse estimation of angular frequency
%   tensor_idx: the location of peak in FFT tensor
%
%   author: Menghuai Xu
%   email: 22134019@zju.edu.cn
%   date: 2021.11.7

    % set the parameters which used in following program
    [Nx, My, Lz] = size(y_tensor);                                      % the size of receiving signal
    NML = Nx * My * Lz;

    % the size vector after oversampling
    Oversam_vec = oversample_rate .* [Nx, My, Lz];

    % 3D-fft and find the maximum index
    Y_spec = fftn(y_tensor, Oversam_vec) / sqrt(NML);
    a3D_fft_abs = abs(Y_spec);

    [~, idx_max] = max(a3D_fft_abs, [], 'all', 'linear') ;

    [row, col, page] = ind2sub(Oversam_vec,idx_max);
    tensor_idx = [row, col, page]';

    % get the angular frequency
    omega_est = (tensor_idx - 1) ./ (Oversam_vec') * 2 * pi;

end


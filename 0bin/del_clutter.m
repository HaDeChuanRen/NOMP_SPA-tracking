function [history_copy, K_est] = del_clutter(history_all, res_rate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    history_copy = history_all;
    [~, sample_time, K_est] = size(history_copy);
    K_del = [];
    for k_idx = 1 : K_est
        testvec_k = squeeze(history_copy(1, :, k_idx));
        if sum(isnan(testvec_k)) > sample_time * res_rate
            K_del = [K_del, k_idx];
        end
    end
    history_copy(:, :, K_del) = [];
    K_est = size(history_copy, 3);
end
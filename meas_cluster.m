function Zarray_clu = meas_cluster(Zarray_lse, weight_z, thers_z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('weight_z','var'), weight_z = [1; 1; 1];
elseif isempty(weight_z), weight_z = [1; 1; 1]; end

if ~exist('thers_z','var'), thers_z = 1;
elseif isempty(thers_z), thers_z = 1; end

Z_temp = Zarray_lse;
Zarray_clu = [];
while(~isempty(Z_temp))
    Z_now = Z_temp(1, :);
    Z_delta = abs(Z_now - Z_temp);
    clu_idx = find((Z_delta * weight_z) < thers_z);
    if size(clu_idx, 1) > 1
        Z_new = mean(Z_temp(clu_idx, :));
    else
        Z_new = Z_now;
    end
    Zarray_clu = [Zarray_clu; Z_new];
    Z_temp(clu_idx, :) = [];
end
end
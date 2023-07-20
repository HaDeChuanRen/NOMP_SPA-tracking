function [pair_mat, label_dis, m_new] = NN_DA(target_cell, Zarray_t, Thres_err)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Khat_last = length(target_cell);
state_dim = target_cell(1).state_dim;
XPri_t = zeros(state_dim, Khat_last);
for k_idx = 1 : Khat_last
    [XPri_kt, Pmat_pri, target_cell(k_idx)] = ...
        target_cell(k_idx).state_transform(Amat_state);
    XPri_t(:, k_idx) = XPri_kt;
end

errmat_xz = zeros(Khat_last, Mhat_t);
for k_tidx = 1 : Khat_last
    % distance calculate and P_G verification
    ZPri_kt = Hmat_mea * XPri_t(:, k_tidx);
    emat_Zkt = ZPri_kt - Zarray_t';
    dvec_kt = diag(emat_Zkt' * emat_Zkt);
    % dvec_kt = sum(emat_Zkt .^ 2, 1);
    errmat_xz(k_tidx, :) = dvec_kt;
    % Pdvec_zmt = exp(- 1/2 * dvec_kt);
    % Pdvec_zmt(Pdvec_zmt < P_G) = 0;
end

m_unmatch = [];
k_unmatch = [];
for m_tidx = 1 : Mhat_t
    if sum(errmat_xz(:, m_tidx) > Thres_err) == Khat_last
        m_unmatch = [m_unmatch, m_tidx];
    end
end
errmat_xz(:, m_unmatch) = [];
for k_tidx = 1 : Khat_last
    if sum(errmat_xz(k_tidx, :) > Thres_err) == Mhat_t
        k_unmatch = [k_unmatch, k_tidx];
    end
end
errmat_xz(k_unmatch, :) = [];
m_idxres = setdiff(1 : Mhat_t, m_unmatch);
k_idxres = setdiff(1 : Khat_last, k_unmatch);
[k_label, m_label] = Hungarianmatch(errmat_xz);

end
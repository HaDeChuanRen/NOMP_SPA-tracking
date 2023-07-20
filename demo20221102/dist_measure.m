function [Z_new,dist,index] = dist_measure(R, H, P, X, Z, gamma) 

% 某个目标对所有测量马氏距离
% R：测量协方差 
% H：量测矩阵
% P：状态协方差矩阵
% X：目标的状态
% Z：各个测量

Z_num = size(Z,2);
Sj = R + H*P*H';
nu = Z - H*repmat(X,[1 Z_num]) ;
dist = [];
for l = 1 : Z_num
dist(l) = nu(:,l)'* inv(Sj) * nu(:,l) + log(det(Sj));       %根据matlab加了惩罚项
end 

index = find(dist <= gamma);              %在范围内的index
Z_new = Z(:,index);

function [X,P,K]=Kalman(X_Forward,P_Forward,Z,A,G,Q,H,R)

%       Z--观测数据矢量 
%       A--系统模型状态矩阵 
%       G--系统模型噪声系数矩阵 
%       Q--系统模型噪声方差 
%       H--量测系数矩阵 
%       R--量测模型噪声协方差 
%       X_Forward--前次估计状态矢量 
%       P_Forward--前次估计状态协方差矩阵
%       X--输出估计状态矢量 
%       P--输出估计状态协方差矩阵 

% 预测 

X_Pre=A*X_Forward; 
P_Pre=A*P_Forward*A'+G*Q*G'; 

% 增益矩阵 

K=P_Pre*H'*inv(H*P_Pre*H'+R)'; 

% Pzz = H*P_Forward*H'+ R;                   %S(k+1/k+1) 新息协方差 
% Pxz = P_Forward*H' ;                       %状态与量测之间的协方差 
% K = P_Forward*H'*(inv(Pzz));               %K(k+1) 增益 


% 修正滤波值和误差协方差阵 

X=A*X_Forward+K*(Z-H*(A*X_Forward)); 
M=K*H; 
n=size(M); 
I=eye(n); 
P=(I-K*H)*P_Pre*(I-K*H)'+ K*R*K'; 

 



 function [Combine_Z,Combine_R,beta_i,a] = PDA(Z_Matrix, PZ_Matrix, Z_Predict, PZ_Predict) 

% 概率数据关联，杂波空间密度为泊松分布随机变量 
% Z_Matrix：阈值内的所有有效量测值 
% PZ_Matrix：量测值的协方差矩阵
% Z_Predict：根据模型预测的那个测量
% PZ_Predict：预测量测的协方差矩阵
% Combine_R：根据各个量测权重得到的组合测量
% Combine_R：组合量测对应的协方差 
% 中间变量： 
% beta为正确关联概率 

lamda=0.0004; 
% lamda=5; 
Pd=1;               %检测概率，当不取1时，后面的a计算出来都是0 
Pg=0.999;          %门限概率 


nm=size(Z_Matrix); 
n=nm(2);   % 量测数量 
m=nm(1);  % 测量维数 

for i=1 :n 

    e(:,i)=Z_Matrix(:,i)-Z_Predict; 

    S(:,:,i)=PZ_Predict+PZ_Matrix(:,:,i);  %新息协方差 X、R、Q互不相关条件下 

    % 百度的算法 

%     a(i)=exp((-1/2)*e(i)'*inv(S(i))*e(i));      
%     bk(i)=lamda*sqrt((2*pi)*det(S(i)))*(1-Pd*Pg)/Pd;    

    % 论文中的算法 

    a(i)=Pd*exp((-1/2)*(e(:,i)'*inv(S(:,:,i))*e(:,i)));          
    bk(i)=lamda*(sqrt(2*pi))^m*sqrt(det(S(:,:,i)))*(1-Pd);   

%     a(i)=exp((-1/2)*(e(:,i)'*inv(S(:,:,i))*e(:,i)));    
%     bk(i)=lamda*sqrt((2*pi)*det(S(:,:,i)))*(1-Pd*Pg)/Pd; 
end 

for i=1:n 

    beta_i(i)=a(i)/(bk(i) + sum(a));     

end 
% beta_0=
% 扩充正确关联概率，使得每一维量测都有对应的关联概率 

beta = beta_i; 

for i=1:m-1 
    beta=[beta;beta_i]; 
end 

M = beta.*Z_Matrix; 
Combine_Z=sum(M',1); 
Combine_Z=Combine_Z'; 
Combine_R=0; 

for i=1:n 
    Combine_R = Combine_R + (beta(:,i)*beta(:,i)').*PZ_Matrix(:,:,i); 
end 

beta_i(n); 
 end
 
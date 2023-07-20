% 不加杂波，且需要处理摘出运动轨迹
function [X0] = pda_grh(time,A,Q,H,C,M,x,X0) 
% z:量测
% npts:测量次数
% A:状态转移矩阵
% Q:过程噪声
% H:观测矩阵
% C:观测噪声方差阵
% M:协方差矩阵
% x:量测数据
% X0:输入初始估计位置，输出PDA结果
lambda=0.04;%单位面积内的虚假量测总数
PG=0.9997;%门概率质量
PD=1;%波门判断依据
Chi=16;%对应2维 门概率质量PG=0.9997
for i=2:time
    %预测
    X0(:,i)=A*X0(:,i-1);%预测状态向量
    M=A*M*A'+Q;%预测协方差
    S=H*M*H'+C;%新息协方差矩阵
    %PDA
    xy_temp=[];%落入波门的全部量测
    e=[];
    beta=[];%权矢量
    MeasurementPoint=x(:,i);%当前时刻的全部量测
    
    for j=1:size(MeasurementPoint,2)
        y_wave=H*X0(:,i)-MeasurementPoint(:,j);
        D=y_wave'*inv(S)*y_wave;%统计距离
        if(D<Chi)
            xy_temp=[xy_temp MeasurementPoint(:,j)];
        end        
    end
    if isempty(xy_temp)%波门内没有量测 就直接用KF预测值作为此刻值
        X0(:,i)=X0(:,i);
    else
        for j=1:size(xy_temp,2)
            y_wave=H*X0(:,i)-xy_temp(:,j);
            e(j)=exp(-1/2*(y_wave'*inv(S)*y_wave));
        end
        b=lambda*(sqrt(2*pi))^size(xy_temp,2)*sqrt(det(2*pi*S))*(1-PG*PD)/PD;
        beta(1)=b/(b+sum(e));
        for j=1:size(xy_temp,2)
            beta(j+1)=e(j)/(b+sum(e));
        end
        X_filter=[];
        K=M*H'/(H*M*H'+C);%增益
        for j=1:size(xy_temp,2) 
            X_filter(:,j)=X0(:,i)+K*(xy_temp(:,j)-H*X0(:,i));%滤波状态向量
        end
        X_filter=[X0(:,i) X_filter];
        
        v=0;%组合新息
        v1=0;
        for j=1:size(xy_temp,2)
            v=v+beta(j+1)*(H*X0(:,i)-xy_temp(:,j));
            v1=v1+beta(j+1)*(H*X0(:,i)-xy_temp(:,j))*(H*X0(:,i)-xy_temp(:,j))';
        end
        Perr=K*(v1-v*v')*K';
        M=M-(1-beta(1))*K*S*K'+Perr;%协方差矩阵
        X_final=zeros(4,1);
        for j=1:size(xy_temp,2)+1
            X_final=X_final+beta(j)*X_filter(:,j);
        end
        X0(:,i)=X_final;
    end
end

end
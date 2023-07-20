%% PDA
T=1;%采样间隔 s
npts=100;%采样总时长
t=0:T:T*(npts-1);%时间轴
X=zeros(4,npts);%真实位置
X(:,1)=[200 0 1e4  -15]';%初始位置 x  vx  y vy
F=[1 T 0 0;0 1 0 0 ; 0 0 1 T;0 0 0 1];%状态转移矩阵
H=[1 0 0 0; 0 0 1 0];%观测矩阵
for i=2:npts
    X(:,i)=F*X(:,i-1);%递推
end
z=zeros(2,npts);%量测
sigma=sqrt(200);%量测噪声
W=[sigma,sigma]';%量测噪声矩阵
R=[sigma^2 0 ;0 sigma^2];%W的协方差矩阵 
for i=1:npts
    z(:,i)=H*X(:,i)+W.*wgn(2,1,0);
end
X0=zeros(4,npts);%状态向量
%二点差分法确定初始位置 误差协方差矩阵
X0(:,2)=[z(1,2) (z(1,2)-z(1,1))/T  z(2,2) (z(2,2)-z(2,1))/T]';%初始位置
P=[R(1,1) R(1,1)/T R(1,2) R(1,2)/T;...%初始协方差矩阵
    R(1,1)/T 2*R(1,1)/T^2 R(1,2)/T  2*R(1,2)/T^2;...
    R(2,1) R(2,1)/T R(2,2) R(2,2)/T;...
    R(2,1)/T 2*R(2,1)/T^2 R(2,2)/T  2*R(2,2)/T^2];

V=wgn(1,2,0);%噪声
Q=[T^4/4*(V'*V)  T^3/2*(V'*V) ;T^3/2*(V'*V) T^2*(V'*V)];%过程噪声

%3~12秒 用KF滤波
for i=3:15
    X0(:,i)=F*X0(:,i-1);%预测
    P=F*P*F'+Q;%预测
    K=P*H'/(H*P*H'+R);%增益
    X0(:,i)=X0(:,i)+K*(z(:,i)-H*X0(:,i));%平滑
    P=(eye(4)-K*H)*P;%平滑
end

figure;
p1=plot(X(1,:),X(3,:),'-',z(1,:),z(2,:),'o','MarkerEdgeColor','k');
grid on;
hold on;


PG=0.9997;%门概率质量
PD=1;%目标检测概率
Chi=16;%对应2维 门概率质量PG=0.9997
clutter=[];%虚假航迹
%13秒后引入杂波 利用新息协方差 计算确认区域面积
gamma=16;
lambda=0.0004;%单位面积内的虚假量测总数
for i=16:npts
    %预测
    X0(:,i)=F*X0(:,i-1);%预测状态向量
    P=F*P*F'+Q;%预测协方差
    S=H*P*H'+R;%新息协方差矩阵  雷达数据处理及应用 Eq 2.96

    %虚假航迹生成
    Av=abs(pi*gamma*sqrt(det(S)));%确认区域面积
    nc=fix(10*Av*lambda+1);%该时刻虚假航迹总数
    q=sqrt(10*Av)/2;
    a=z(1,i)-q;b=z(1,i)+q;
    c=z(2,i)-q;d=z(2,i)+q;
    x_temp=a+(b-a)*rand(1,nc);%虚假航迹X坐标
    y_temp=c+(d-c)*rand(1,nc);%虚假航迹Y坐标
    clutterNow=[x_temp ;y_temp];
    clutter=[clutter clutterNow];
    p2=plot(x_temp,y_temp,'*','MarkerEdgeColor','r');
    hold on;
    %PDA
    xy_temp=[];%落入波门的全部量测
    e=[];
    b=[];
    beta=[];%权矢量
    Vk=pi*Chi*sqrt(det(S));%相关波门体积
    MeasurementPoint=[z(:,i) clutterNow];%当前时刻的全部量测 包括真实量测 虚假量测
    
    for j=1:size(MeasurementPoint,2)
        d=H*X0(:,i)-MeasurementPoint(:,j);
        D=d'*inv(S)*d;%统计距离
        if(D<Chi)
            xy_temp=[xy_temp MeasurementPoint(:,j)];
        end
    end
    if isempty(xy_temp)%波门内没有量测 就直接用KF预测值作为此刻值
        X0(:,i)=X0(:,i);
        P=P;
    else
        for j=1:size(xy_temp,2)
            d=H*X0(:,i)-xy_temp(:,j);
            e(j)=exp(-1/2*(d'*inv(S)*d));
        end
        b=size(xy_temp,2)/Vk*sqrt(det(2*pi*S))*(1-PG*PD)/PD;
        beta(1)=b/(b+sum(e));
        for j=1:size(xy_temp,2)
            beta(j+1)=e(j)/(b+sum(e));
        end
        X_filter=[];
        P_filter=[];
        K=P*H'/(H*P*H'+R);%增益
        for j=1:size(xy_temp,2) 
            X_filter(:,j)=X0(:,i)+K*(xy_temp(:,j)-H*X0(:,i));%滤波状态向量
        end
        X_filter=[X0(:,i) X_filter];
        Pc=(eye(4)-K*H)*P;%滤波协方差矩阵
        v=0;%组合新息
        v1=0;
        for j=1:size(xy_temp,2)
            v=v+beta(j+1)*(H*X0(:,i)-xy_temp(:,j));
            v1=v1+beta(j+1)*(H*X0(:,i)-xy_temp(:,j))*(H*X0(:,i)-xy_temp(:,j))';
        end
        Perr=K*(v1-v*v')*K';
        P=P*beta(1)+(1-beta(1))*Pc+Perr;
        X_final=zeros(4,1);
        for j=1:size(xy_temp,2)+1
            X_final=X_final+beta(j)*X_filter(:,j);
        end
        X0(:,i)=X_final;
    end
end
p3=plot(X0(1,3:npts),X0(3,3:npts),'-v','MarkerEdgeColor','b');
axis([80 240 8500 10000]);
legend([p1(1),p1(2),p2,p3],{"真实航迹","真实量测","虚假量测","PDA估计轨迹"});


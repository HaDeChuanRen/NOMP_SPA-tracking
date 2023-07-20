% 仿真（仅PDA）

clear;
close all;
clc;
linWid = 1;
rng(5)
add_clutter = 1;%0不添加杂波，1添加杂波

%% ideal track
T = 1; % sample interval
time = 100;
n = 0:T:time;
r_x = 10-0.2*n;
r_y = -5+0.2*n;
% figure(1)
% plot(r_x,r_y,'--')
% hold on;
%% true track
% state equation
A = [1 0 T 0;
     0 1 0 T;
     0 0 1 0;
     0 0 0 1];
mu = 0;
var_u = 1e-4;                       % variance of sigma u

u_x = mvnrnd(mu, var_u, time);      % noise perturbations of v
u_y = mvnrnd(mu, var_u, time);
s = zeros(4,time);
s(:,1) = [10 -5 -0.2 0.2]';         % initial state
u = zeros(4,time);

for i = 1:1:time
    u(:,i) = [0 0 u_x(i) u_y(i)]';
    s(:,i+1) = A*s(:,i)+u(:,i);
end
figure(1)
plot(s(1,:),s(2,:))
hold on

%% add noise
var_R = 0.001;
mu_R = 0;
var_beta = 0.001;
mu_beta = 0;
R = sqrt(s(1,:).^2+s(2,:).^2);
beta = atan2(s(2,:),s(1,:));
omega_R = mvnrnd(mu_R, var_R, time+1);                  % noise perturbations of R and beta
omega_beta = mvnrnd(mu_beta, var_beta, time+1);
R_head = R+omega_R';
beta_head = beta+omega_beta';

for i = 1:1:time
    rx_head(i) = R_head(i)*cos(beta_head(i));
    ry_head(i) = R_head(i)*sin(beta_head(i));
end
% figure(1)
% plot(rx_head(i),ry_head(i),'o')
% hold on


%% parameters
A = [1 0 T 0;
     0 1 0 T;
     0 0 1 0;
     0 0 0 1];                          %状态转移矩阵
C = [10 0;
     0 10];                            %观测噪声方差阵
x = [rx_head;
     ry_head];

K = cell(1,101);                        %增益
H = cell(1,101);
for i = 1:1:time+1
    H{1,i}=[1 0 0 0; 0 1 0 0];          %观测矩阵
end
M = cell(1,101);
M{1,1} = 100*eye(4);
% V=sqrt(10^(0/10)).*rand(1,2);%噪声
% Q=[T^4/4*(V'*V)  T^3/2*(V'*V) ;T^3/2*(V'*V) T^2*(V'*V)];%过程噪声
Q = [0 0 0 0;
     0 0 0 0;
     0 0 0.0001 0;
     0 0 0 0.0001];


%% PDA without clutter

    X0 = zeros(4,time+1);
    X0(:,1)=[rx_head(1) ry_head(1) (rx_head(2)-rx_head(1))/T (ry_head(2)-ry_head(1))/T ]';%初始位置
    [X0] = pda_grh(time,A,Q,H{1},C,M{1},x,X0);
    figure(1)
    for i = 1:1:time-1
        plot(rx_head(i),ry_head(i),'o')
        hold on
        plot([X0(1,i) X0(1,i+1)],[X0(2,i) X0(2,i+1)],'c-^','LineWidth',linWid)
        hold on;
        pause(0.05)
    end
    % plot(X0(1,3:time),X0(2,3:time),'k-v','LineWidth',linWid);
    legend('ideal track','true track','observed measurement','PDA track')
    figure(2)
    for i = 3:1:time
        MSE_rx2(i) = norm(X0(1,i)-s(1,i));    
        MSE_ry2(i) = norm(X0(2,i)-s(2,i));
    end
    subplot(1,2,1)
    plot(MSE_rx2);
    xlabel('time');
    ylabel('PDA-MSE rx')
    subplot(1,2,2)
    plot(MSE_ry2);
    xlabel('time');
    ylabel('PDA-MSE ry')

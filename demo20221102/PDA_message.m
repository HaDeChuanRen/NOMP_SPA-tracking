%多目标追踪（结合了消息传递，但并没有存在状态范畴判断）
% 二维空间匀速直线运动，状态向量为X=[x,vx,y,vy] 

clc; 
clear; 
close all; 


tic;
%% ************************************************ 
%          参数设置 
%************************************************ 
r = 15;                                  %量测噪声
lambda_c=0;                      %噪声数（泊松）

rng(2);

I=eye(4); 
T = 1;                                  %采样间隔 
SampleTime = 100 ;                        %仿真步数 
A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];   %实际模型：CV 
H = [1 0 0 0;0 0 1 0];                  %测量模型 
Q = 0;                                  %实际过程噪声 
G = [T^2/2 0; T 0; 0 T^2/2; 0 T];    %噪声加权矩阵 
R = [r 0; 0 r];                         %量测噪声矩阵 
Pd = 1;                                 %检测概率
X0 = [0  -50  -50   -50;...
      0  1     1    1;...
      50  0   -25   25;...
      -1  0    0.5    -0.5];  %初始状态 
n = 4;                                      %目标数
X0_est = X0;                           %初始的预估状态
M = 20;                              %蒙特卡洛次数
mistake = 0;                        %记录失败数
PG = 0.999;                         %监测落入监测门的概率
gamma = chi2inv(PG,2);   %inv chi^2 dn gamma value   %监测门参数，测量是两维的  
range_c= [ -75 75;-75 75 ];      %噪声生成区域
pdf_c= 1/prod(range_c(:,2)-range_c(:,1)); %噪声的概率分布


error1_PDA_all=0;error2_PDA_all=0;error3_PDA_all=0;error4_PDA_all=0;
error1_PDA_all_temp = 0;error2_PDA_all_temp = 0;error3_PDA_all_temp = 0;error4_PDA_all_temp = 0;   
for m_interate=1:M
  P = cell(n,SampleTime);     %储存协方差P
    for repeat=1:n
        X{repeat}(:,1)=X0(:,repeat);
        Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
        Zk{repeat}(:,1)=H*X{repeat}(:,1)+Vk; 

%% ************************************************ 
%          量测生成 
%************************************************ 

        for i=2:SampleTime 
        X{repeat}(:,i)=A*X{repeat}(:,i-1);          % 真实状态 
        Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
        Zk{repeat}(:,i)=H*X{repeat}(:,i)+Vk;      %生成量测值 
        end 


    %************************************************ 
    %          data association初始化 
    %************************************************ 

        Xk_PDA{repeat}=X0_est(:,repeat);  %初始状态、与实际值略有差别 
        R11=r; R22=r; R12=0; R21=0; 
        Pkk_PDA=[R11 R11/T R12 R12/T; 
        R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
        R21 R21/T R22 R22/T; 
        R21/T 2*R21/T^2 R22/T 2*R22/T^2];   %初始协方差 
        Xkk{repeat} = Xk_PDA{repeat} ; 
        Pkk{repeat} = Pkk_PDA; 
        X_Pre{repeat} = A*Xkk{repeat}; 
        P_Pre{repeat}=A*Pkk{repeat}*A'+G*Q*G'; 
    %     clear P;
        P{repeat,1} = R; 
    end



    for i=1:SampleTime 
    %************************************************ 
    %          产生杂波与测量
    %************************************************ 
    N_c= poissrnd(lambda_c);                                                  %杂波数量
    C= repmat(range_c(:,1),[1 N_c])+ diag(range_c*[ -1; 1 ])*rand(2,N_c);      %生成杂波
    % 生成代表杂波的nc个虚假量测 
    Z_temp=[];
    PZ_Matrix{i} = cat(3);          %建立空矩阵
        for repeat=1:n
            if ( rand <= Pd )
            Z_temp=[Z_temp Zk{repeat}(:,i)];
            PZ_Matrix{i} = cat(3,PZ_Matrix{i},R); 
            end
        end

    Z_Matrix{i}=[ Z_temp,C ];      %测量矩阵

        for j=1:N_c 
            PZ_Matrix{i} = cat(3,PZ_Matrix{i},[r,0;0,r]); 
        end 
    end




for i=1:1:SampleTime     
    %************************************************ 
    %          数据关联 
    %************************************************ 
    assocaite_w_unnomarlized{i} = [];
    
        for repeat=1:n
         Z_Predict = H*X_Pre{repeat}; 
         PZ_Predict= H*P_Pre{repeat}*H';
         Z_Matrix_temp = Z_Matrix{i};
         PZ_Matrix_temp = PZ_Matrix{i};

              [unused1,unused2,beta_i,a]=PDA(Z_Matrix_temp, PZ_Matrix_temp, Z_Predict, PZ_Predict) ; % PDA 
%              end
             
%             if sum(a)<1e-2                                       %%监测偏离的轨迹
%             Combine_Z = Z_PDA{repeat}(:,i-1);    Combine_R = PZ_Predict;
%             end
        assocaite_w_unnomarlized{i} = [assocaite_w_unnomarlized{i} ; beta_i];
        end
        
       
        
% 消息传递       
       assocaite_w_unnomarlized{i} = (assocaite_w_unnomarlized{i})';
         m =  size(Z_Matrix{i},2);
        om = ones(1,m);
        on = ones(1,n);
        muba = ones(m,n);
         delta = ones(m,n);    %变化量的表示
         interate(m_interate,i) = 1;
          while(max(max(abs(delta)))>1e-6)
          muba_pre = muba;
          prodfact = muba .* assocaite_w_unnomarlized{i};
          sumprod = 1 + sum(prodfact);
          muab = assocaite_w_unnomarlized{i} ./ (sumprod(om,:) - prodfact);
          summuab = 1 + sum(muab,2);
          muba = 1 ./ (summuab(:,on) - muab);
          delta = muba - muba_pre;
          interate(m_interate,i) = interate(m_interate,i) + 1;
          end
          
%     for interate = 1:100
%         prodfact = muba .* assocaite_w_unnomarlized{i};
%         sumprod = 1 + sum(prodfact);
%         muab = assocaite_w_unnomarlized{i} ./ (sumprod(om,:) - prodfact);
%         summuab = 1 + sum(muab,2);
%         muba = 1 ./ (summuab(:,on) - muab);
%     end     
    
    
       association_p_after = [ muba.* assocaite_w_unnomarlized{i}];
       diffence{i} = association_p_after - assocaite_w_unnomarlized{i};
%        association_p_after = [  assocaite_w_unnomarlized{i}];
       
       
  % 得到组合测量
          for repeat=1:n
              clear association_p_after_normal;
              association_p_after_normal(:,repeat)=association_p_after(:,repeat)/sum(association_p_after(:,repeat)); 
              association_p_after_normal_temp = [association_p_after_normal(:,repeat)'; association_p_after_normal(:,repeat)'];
               M_temp = association_p_after_normal_temp.*Z_Matrix_temp; 
               Combine_Z = sum(M_temp',1); 
               Combine_Z = Combine_Z'; 
               Combine_R = 0; 
                   for k = 1:size(Z_Matrix{i},2)
                   Combine_R  = Combine_R + (association_p_after_normal_temp(:,k)*association_p_after_normal_temp(:,k)').*PZ_Matrix_temp(:,:,k); 
                   end
            
        %************************************************ 
        %          卡尔曼滤波 
        %************************************************ 
        Z_PDA{repeat}(:,i) = Combine_Z ; 
        P{repeat,i+1}=Combine_R; 
        [Xk_after,Pk_after,Kk_PDA]=Kalman(Xkk{repeat},Pkk{repeat},Combine_Z,A,G,Q,H,P{repeat,i+1}); 
        Xkk{repeat}=Xk_after;     Pkk{repeat}=Pk_after; 
        % 预测 
        X_Pre{repeat}=A*Xkk{repeat}; 
        P_Pre{repeat}=A*Pkk{repeat}*A'+G*Q*G'; 
        %出各个状态值 
        Ex_PDA{repeat}(i)=Xkk{repeat}(1); 
        Evx_PDA{repeat}(i)=Xkk{repeat}(2); 
        Ey_PDA{repeat}(i)=Xkk{repeat}(3); 
        Evy_PDA{repeat}(i)=Xkk{repeat}(4); 
        error1_PDA{repeat}(i)=Ex_PDA{repeat}(i)-X{repeat}(1,i);%Pkk(1,1); 
        error2_PDA{repeat}(i)=Ey_PDA{repeat}(i)-X{repeat}(3,i);%Pkk(2,2); 
        error3_PDA{repeat}(i)=Evx_PDA{repeat}(i)-X{repeat}(2,i);%Pkk(3,3); 
        error4_PDA{repeat}(i)=Evy_PDA{repeat}(i)-X{repeat}(4,i);%Pkk(4,4); 
          end
       
    end

    
    
 %------------------------------------------------------   
    
 error1_PDA_all=0;error2_PDA_all=0;error3_PDA_all=0;error4_PDA_all=0;

    for repeat=1:n
    error1_PDA_all = abs(error1_PDA{repeat}) + abs(error1_PDA_all);
    error2_PDA_all = abs(error2_PDA{repeat}) + abs(error2_PDA_all);
    error3_PDA_all = abs(error3_PDA{repeat}) + abs(error3_PDA_all);
    error4_PDA_all = abs(error4_PDA{repeat}) + abs(error4_PDA_all);
    end
    
    if(max(error1_PDA_all/n)>20)||(max(error3_PDA_all/n)>20)||(max(error2_PDA_all/n)>20)||(max(error4_PDA_all/n)>20)
    mistake=mistake+1;
    end 

     error1_PDA_all_temp = error1_PDA_all_temp + error1_PDA_all;
    error2_PDA_all_temp = error2_PDA_all_temp + error2_PDA_all;
    error3_PDA_all_temp = error3_PDA_all_temp + error3_PDA_all;
    error4_PDA_all_temp = error4_PDA_all_temp + error4_PDA_all;    
    
end

mytimer1=toc;
disp(mytimer1);



     error1_PDA_all = error1_PDA_all_temp;
    error2_PDA_all = error2_PDA_all_temp;
     error3_PDA_all = error3_PDA_all_temp;
     error4_PDA_all = error4_PDA_all_temp;
     
    %************************************************ 
    %          画图 
    %*********************************************** 
    i=1:SampleTime; 
    figure();
box on; 
axis equal;

for k=1:n
    plot(X{k}(1,i),X{k}(3,i),'b--','LineWidth',2);    %真实值 
    hold on;
    plot(Ex_PDA{k}(i),Ey_PDA{k}(i),'r-','LineWidth',2);    %滤波值 
    grid on;
    plot(Zk{k}(1,i),Zk{k}(2,i),'*');                    %实际测量值 
    plot(Z_PDA{k}(1,i),Z_PDA{k}(2,i),'o');       %组合测量值 
    % text(X(1,1)+1,X(3,1)+5,'t=1'); 
end
 legend('真实值','滤波值','实际量测','组合量测'); 
    title('目标运动轨迹'); xlabel('x/m'); ylabel('y/m'); 
    
%     
%     i=1:SampleTime; 
%     figure();
% box on; 
% axis equal;
% for i=1:SampleTime
%     for k=1:n
%     plot(X{k}(1,1:i),X{k}(3,1:i),'b--','LineWidth',2);    %真实值 
%     hold on;
%     plot(Ex_PDA{k}(1:i),Ey_PDA{k}(1:i),'r-','LineWidth',2);    %滤波值 
%     grid on;
%     plot(Zk{k}(1,i),Zk{k}(2,i),'*');                    %实际测量值 
%     plot(Z_PDA{k}(1,i),Z_PDA{k}(2,i),'o');       %组合测量值 
%     % text(X(1,1)+1,X(3,1)+5,'t=1'); 
%   
%     end
% %      pause(0.01);
% end
%  legend('真实值','滤波值','实际量测','组合量测'); 
%     title('目标运动轨迹'); xlabel('x/m'); ylabel('y/m'); 


 i=1:SampleTime; 
    %位置误差 
    figure(); 
    subplot(211); hold on;
    plot(abs(error1_PDA_all(i))/n/M,'LineWidth',2); grid on;      temp1=abs(error1_PDA_all(i))/n/M;
    load error1_PDA_all; 
    plot(abs(error1_PDA_all(i))/n/M,'LineWidth',2); grid on; 
    legend('结合消息传递','未结合消息传递'); 
    title('位置误差'); xlabel('t/s'); ylabel('error-x/m'); 
    subplot(212);  hold on;
    plot(abs(error2_PDA_all(i))/n/M,'LineWidth',2); grid on;      temp3=abs(error2_PDA_all(i))/n/M;
    load error2_PDA_all; 
    plot(abs(error2_PDA_all(i))/n/M,'LineWidth',2); grid on; 
    xlabel('t/s'); ylabel('error-y/m'); 
     legend('结合消息传递','未结合消息传递'); 

  
     i=1:SampleTime; 
    %位置误差 
    figure(); 
    subplot(211); hold on;
    plot(abs(error3_PDA_all(i))/n/M,'LineWidth',2); grid on;      temp2=abs(error3_PDA_all(i))/n/M;
    load error3_PDA_all; 
    plot(abs(error3_PDA_all(i))/n/M,'LineWidth',2); grid on; 
    legend('结合消息传递','未结合消息传递'); 
    title('速度误差'); xlabel('t/s'); ylabel('error-vx/m/s'); 
    subplot(212);  hold on;
    plot(abs(error4_PDA_all(i))/n/M,'LineWidth',2); grid on;      temp4=abs(error4_PDA_all(i))/n/M;
    load error4_PDA_all; 
    plot(abs(error4_PDA_all(i))/n/M,'LineWidth',2); grid on; 
    xlabel('t/s'); ylabel('error-vy/m/s'); 
    legend('结合消息传递','未结合消息传递'); 
    
    
    figure();
    interate = sum(interate,1);
    plot(i,interate/M,'LineWidth',2); 
    title(' 消息传递迭代次数'); xlabel('t/s'); ylabel('n/次'); 
     set(gca,'YLim',[0,15]);   
       
%     %速度误差 
%     figure(); 
%     subplot(211); 
%     plot(abs(error2_PDA_all(i))/n,'LineWidth',2); grid on 
%     title('速度误差'); xlabel('t/s'); ylabel('error-vx/m/s'); 
%     subplot(212) 
%     plot(abs(error4_PDA_all(i))/n,'LineWidth',2); grid on 
%     xlabel('t/s'); ylabel('error-vy/m/s'); 


 



% 
% %plot error
% ospa_vals= n;
% ospa_c= 100;
% ospa_p= 1;
% for k=1:meas.K
%     [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(get_comps(truth.X{k},[1 3]),get_comps(est.X{k},[1 3]),ospa_c,ospa_p);
% end
% 
% figure; ospa= gcf; hold on;
% subplot(3,1,1); plot(1:meas.K,ospa_vals(:,1),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
% subplot(3,1,2); plot(1:meas.K,ospa_vals(:,2),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
% subplot(3,1,3); plot(1:meas.K,ospa_vals(:,3),'k'); grid on; set(gca, 'XLim',[1 meas.K]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
% xlabel('Time');






% end
ratio = (M-mistake)/M;   %正确率

% 
% 
% function Xc= get_comps(X,c)
% 
% if isempty(X)
%     Xc= [];
% else
%     Xc= X(c,:);
% end
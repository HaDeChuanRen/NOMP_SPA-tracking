%多目标追踪（结合了消息传递和存在状态范畴判断）(点云模型)(解决了新出现目标没有找到测量时情况)混合模型
% 二维空间匀速直线运动，状态向量为X=[x,vx,y,vy] 
function [ospa_vals_without,mytimer1] = PDA_based(draw_or_not,c,p,scenario)
% clc; 
% clear; 
% close all; 
% rng(24);
tic;   %计时开始
%% ************************************************ 
%          参数设置 
%************************************************ 
% r = 800;                                  %量测噪声
r = 6;                                  %量测噪声
lambda_c=2;                      %噪声数（泊松）
lambda_t=12;                      %目标引起的测量数（泊松）
epsilon = 2*sqrt(r);
minpts = 5;   %比测量维度大1.
ospa_c= c;
ospa_p= p;
float = 0;
X0 = [0 -500  -500  -500  -200  -200    -400    600;...     
      0+randn*float    10+randn*float   10+randn*float    10+randn*float  10+randn*float     8+randn*float    0+randn*float     -2+randn*float;...
      500  0   -150  150  400    -400    200     300;...
      -10+randn*float  0+randn*float    3+randn*float   -3+randn*float   0+randn*float    4+randn*float    -4+randn*float    -6+randn*float];  %初始状态 
 X0 = X0/10;

I = eye(4); 
T = 1;                                  %采样间隔 
SampleTime = 100 ;                        %仿真步数 
A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];   %实际模型：CV 
H = [1 0 0 0;0 0 1 0];                  %测量模型 
Q = 0;                                  %实际过程噪声 
G = [T^2/2 0; T 0; 0 T^2/2; 0 T];    %噪声加权矩阵 

R = [r 0; 0 r];                         %量测噪声矩阵 
Pd = 1;                                 %检测概率


 
 n = 1;                                      %初始目标数
% n_max = size(X0,2);
n_max = 8;                               %最多可能出现的目标数
n_real = zeros(1,100);
  switch scenario
      case 1
     target_type = [0 0 0 0 0 0 0 0];
      case 2
    target_type = [1 1 1 1 1 1 1 1];
      case 3
    target_type = [0 0 0 0 1 1 1 0];
     case 4
    target_type = [0 0 0 1 1 1 0 0];
  end

dim_x = 4;
dim_z = 2;
appeartime = [1 1 1 1 10 15 20 25];
appeartime_predict = [1 1 1 1];
appeartime_predict = appeartime_predict(1:n);
appear_count_all = 4;
disappeartime = [100 75 100 100 100 100 100 100];
disappeartime_predict = [100 100 100 100];
disappeartime_predict = disappeartime_predict(1:n);
disappear_count_all = 4;
float_dis = 1.5;
X0_est = [0+randn*float_dis -500+randn*float_dis  -500+randn*float_dis  -500+randn*float_dis  -200+randn*float_dis  -200+randn*float_dis    -400+randn*float_dis    600+randn*float_dis;...     
      0+randn*float    10+randn*float   10+randn*float    10+randn*float  10+randn*float     8+randn*float    0+randn*float    -2+randn*float;...
      500+randn*float_dis  0+randn*float_dis   -150+randn*float_dis  150+randn*float_dis  400+randn*float_dis    -400+randn*float_dis    200+randn*float_dis     300+randn*float_dis;...
      -10+randn*float  0+randn*float    3+randn*float   -3+randn*float   0+randn*float    4+randn*float    -4+randn*float    -6+randn*float];  %初始状态 
X0_est = X0_est/10;
% 
% X0_est = X0;                           %初始的预估状态
M = 1;                              %蒙特卡洛次数
mistake = 0;                        %记录失败数
PG = 1-1e-3;                         %监测落入监测门的概率
gamma = chi2inv(PG,2);   %inv chi^2 dn gamma value   %监测门参数，测量是两维的  
gamma = 17;
range_c= [ -750 750;-750 750 ];      %噪声生成区域
range_c= range_c/10;      %噪声生成区域
pdf_c= 1/prod(range_c(:,2)-range_c(:,1)); %噪声的概率分布
ospa_vals_without_all = zeros(SampleTime,3);

for Monte_Carlo=1:M
P = cell(n_max,SampleTime);     %储存协方差P
  
  %真实个数计算
  for i = 1:SampleTime 
      n_real(i) = sum((i<=disappeartime).*(i>=appeartime));
  end
      
  
%% ************************************************ 
%          初始位置与测量生成 
%************************************************ 
truth = cell(SampleTime,1);
est = cell(SampleTime,1);
    for repeat=1:n_max  
        X{repeat} = zeros(dim_x,SampleTime);
        Zk{repeat} = zeros(dim_z,SampleTime);
        X{repeat}(:,appeartime(repeat))=X0(:,repeat);
        Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
        Zk{repeat}(:,appeartime(repeat))=H*X{repeat}(:,appeartime(repeat)); 
        
        truth{appeartime(repeat)} = [truth{appeartime(repeat)}  X{repeat}(:,appeartime(repeat))];
%% ************************************************ 
%          正确量测生成 
%************************************************ 

        for i = appeartime(repeat)+1 : disappeartime(repeat) 
        X{repeat}(:,i)=A*X{repeat}(:,i-1);          % 真实状态 
        Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
        Zk{repeat}(:,i)=H*X{repeat}(:,i);      %生成量测值 
        truth{i} = [truth{i}  X{repeat}(:,i)];
        end 
    %************************************************ 
    %          data association初始化 
    %************************************************ 

        R11=r; R22=r; R12=0; R21=0; 
        P_init = [R11 R11/T R12 R12/T; 
        R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
        R21     R21/T    R22    R22/T; 
        R21/T 2*R21/T^2 R22/T 2*R22/T^2];   %初始协方差 
    
        Xkk{repeat} = X0_est(:,repeat);    %初始状态
        Pkk{repeat} = P_init; 
        X_Pre{repeat} = A*Xkk{repeat}; 
        P_Pre{repeat}=A*Pkk{repeat}*A'+G*Q*G'; 
        P{repeat,appeartime(repeat)} = R;      
    end

    
    
    
    
    
    
    
    
    %% ************************************************ 
    %          结合杂波与测量
    %************************************************ 
    for i = 1:SampleTime 
    N_c = poissrnd(lambda_c);                                                  %杂波数量
    C= repmat(range_c(:,1),[1 N_c])+ diag(range_c*[ -1; 1 ])*rand(2,N_c);      %生成杂波
    % 生成代表杂波的nc个虚假量测 
    Z_temp = [];
    
   
    
    PZ_Matrix{i} = cat(3);          %建立空矩阵
        for repeat=1:n_max
            if (( rand <= Pd ) && (i>=appeartime(repeat)) && (i<=disappeartime(repeat)) )
                if(target_type(repeat)~=0)
                    N_t= poissrnd(lambda_t);  
                    T_temp = repmat(Zk{repeat}(:,i),[1 N_t]) +  sqrt(r)*randn(2,N_t);     %生成杂波
                    Z_temp=[Z_temp T_temp];
                else
                    T_temp = Zk{repeat}(:,i) +  sqrt(r)*randn(2,1);     %生成杂波
                    Z_temp=[Z_temp T_temp];
                end
            
            end
        end

    Z_Matrix{i}=[ Z_temp,C ];      %测量矩阵

    end



 %% ************************************************ 
    %          数据关联 
    %************************************************ 
 exist_vector = ones(1,n);               %目标存在为1，不存在是0
 n_predict = zeros(1,SampleTime);        %目标数
 n_predict(1) = n;     
 disappear_count = zeros(1,n);
 
 Z_Matrix_backup = Z_Matrix;
    clear Z_Matrix; 
 
 for i=1:SampleTime     
   
        
    Z_Matrix_backup{i}=Z_Matrix_backup{i}';
    
%     rowrank = randperm(size(Z_Matrix_backup{i}, 1)); % 随机打乱的数字，从1~行数打乱
%     Z_Matrix_backup{i} = Z_Matrix_backup{i}(rowrank, :);%%按照rowrank打乱矩阵的行数
    
    labels{i} = dbscan(Z_Matrix_backup{i},epsilon,minpts);
    Z_Matrix{i} =[];
    PZ_Matrix{i} = cat(3);          %建立空矩阵
    for j = 1:max(labels{i}) 
        cloud =  Z_Matrix_backup{i}(find(labels{i} == j),:);
        Z_Matrix{i} = [Z_Matrix{i} mean(cloud)'];
        PZ_Matrix{i} = cat(3, PZ_Matrix{i},cov(cloud));
    end
        cloud =  Z_Matrix_backup{i}(find(labels{i} == -1),:);
        Z_Matrix{i} = [Z_Matrix{i} cloud'];
        
        for j=1:size(find(labels{i} == -1),1) 
            PZ_Matrix{i} = cat(3,PZ_Matrix{i},R);         
        end
        
       
    assocaite_w_unnomarlized{i} = [];
    exist_index = find(exist_vector);        %生存的目标的序号       
        for repeat = exist_index             % 先求每个目标的单个PDA
            
         Z_Predict = H*X_Pre{repeat}; 
         PZ_Predict= H*P_Pre{repeat}*H';
         Z_Matrix_temp = Z_Matrix{i};
         PZ_Matrix_temp = PZ_Matrix{i};

         [unused1,unused2,beta_i,a]=PDA(Z_Matrix_temp, PZ_Matrix_temp, Z_Predict, PZ_Predict) ; % PDA 
             
%             if sum(a)<1e-2                                       %%监测偏离的轨迹
%             Combine_Z = Z_PDA{repeat}(:,i-1);    Combine_R = PZ_Predict;
%             end

        assocaite_w_unnomarlized{i} = [assocaite_w_unnomarlized{i} ; beta_i];  
        end
 
       assocaite_w_unnomarlized{i} = (assocaite_w_unnomarlized{i})';
        association_p_after = assocaite_w_unnomarlized{i};          %消息传递前的association
       
       
       
       index = [];
       index_all = [];
       
       out_of_circle = [];
       for repeat = exist_index             % 先求每个目标的单个PDA
           %消失判定+最有可能测量判定
         [~,dist_indicator{repeat},index] = dist_measure(R, H, P_Pre{repeat}, X_Pre{repeat}, Z_Matrix{i}, gamma);
             if (isempty(index) == 1)
              disappear_count(repeat) = disappear_count(repeat) + 1; 
               out_of_circle = [out_of_circle , find(repeat == exist_index)];
             else
              disappear_count(repeat) = 0; 
             end 
         index_all = union(index_all,index);
       end 
       
       index_all_indicator{i} = index_all;
       
       
       Z_chosen = max(association_p_after(index_all,:),[],2);
       [~,temp_index]=sort(Z_chosen);
       
       if size(temp_index,1) > n_predict(i)
       temp_index = temp_index(1:size(temp_index,1)-n_predict(i),1);     %记录可能的真实测量
%        disp([num2str(i),'圈中有多余测量']);
       else 
           temp_index = [];
       end
      temp_index = union(index_all(temp_index),setdiff(linspace(1,size(Z_Matrix{i},2),size(Z_Matrix{i},2)),index_all));
      
    if size(temp_index,2)~=1
        temp_index = temp_index';
    end
       %% ************************************************ 
        %               计算组合测量
        %************************************************       

          for repeat=1:n_predict(i)
              clear association_p_after_normal;
              
             
                           
                  association_p_after_normal = association_p_after(:,repeat)/sum(association_p_after(:,repeat)); 
                  association_p_after_normal_temp = [association_p_after_normal'; association_p_after_normal'];
                  
                   M_temp = association_p_after_normal_temp.*Z_Matrix_temp; 
                   Combine_Z = sum(M_temp',1); 
                   Combine_Z = Combine_Z'; %组合测量
                   Combine_R = 0;          %组合协方差
                       for k = 1:size(association_p_after_normal,1)
                       Combine_R  = Combine_R + (association_p_after_normal_temp(:,k)*association_p_after_normal_temp(:,k)').*PZ_Matrix_temp(:,:,k); 
                       end

                
                       
                    %% ************************************************ 
                    %          卡尔曼滤波 
                    %************************************************ 
                    Z_PDA{exist_index(repeat)}(:,i) = Combine_Z ; 
                    P{exist_index(repeat),i+1}=Combine_R; 
                    [Xk_after,Pk_after,Kk_PDA]=Kalman(Xkk{exist_index(repeat)},Pkk{exist_index(repeat)},Combine_Z,A,G,Q,H,Combine_R); 
                    
                   if ismember(repeat,out_of_circle)
                        Xk_after = X_Pre{exist_index(repeat)};
                        Pk_after = P_Pre{exist_index(repeat)}; 
                   end
                    
                   
                   
                    Xkk{exist_index(repeat)}=Xk_after;     Pkk{exist_index(repeat)}=Pk_after; 
                    % 预测 
                    X_Pre{exist_index(repeat)}=A*Xkk{exist_index(repeat)}; 
                    P_Pre{exist_index(repeat)}=A*Pkk{exist_index(repeat)}*A'+G*Q*G'; 
            
                    %% 求出各个状态值 
                    Ex_PDA{exist_index(repeat)}(i)=Xkk{exist_index(repeat)}(1); 
                    Evx_PDA{exist_index(repeat)}(i)=Xkk{exist_index(repeat)}(2); 
                    Ey_PDA{exist_index(repeat)}(i)=Xkk{exist_index(repeat)}(3); 
                    Evy_PDA{exist_index(repeat)}(i)=Xkk{exist_index(repeat)}(4); 
                    est{i} = [est{i} Xkk{exist_index(repeat)}];
                    
                    
          end
       
          

        
        
        
          
          
%% 新生成的目标   
 n_predict(i+1) = n_predict(i) ;

          Z_chosen = Z_Matrix_temp(:,temp_index);
          P_chosen = PZ_Matrix_temp(:,:,temp_index);
          Z_chosen_indictor{i} = Z_chosen;
%              disp(i);
%              disp(size(temp_index));
          if (~isempty(Z_chosen))
             Z_num = size(temp_index,1);
          else
             Z_num = 0;    
          end
          new_appear_count = [];
          new_X_initial = [];
          new_P_initial = [];
          
  %第一帧的初始化
        if i==1
            if (~isempty(temp_index))
          X_initial = [Z_chosen(1,:);zeros(1,Z_num);Z_chosen(2,:);zeros(1,Z_num)];
          P_initial = repmat(P_init,[1,1,Z_num]);
          initial_num = Z_num;           %上一个帧的可能出生对象
          appear_count = ones(1,initial_num);
            else
              X_initial = [];
              P_initial = [];
              initial_num = 0;
              appear_count = [];
            end
         end
          
   % 计算马氏距离------------- 
   if initial_num ~= 0
                 
          for j = 1:initial_num
               Sj= R + H*P_initial(:,:,j)*H';   
              [~,dist,~] = dist_measure(R, H, A*P_initial(:,:,j)*A'+G*Q*G', A*X_initial(:,j), Z_chosen, gamma) ;%!!!!!!!!!!!!!!!
               [min_value,min_index ] = min(dist);
                  if min_value<gamma
%                    temp = [Z_chosen(1,min_index);0;Z_chosen(2,min_index);0];
                      [Xk_after,Pk_after,Kk_PDA]=Kalman(X_initial(:,j),P_initial(:,:,j),Z_chosen(:,min_index),A,G,Q,H,P_chosen(:,:,min_index));                     
%                       new_X_initial = cat(2,new_X_initial ,A*Xk_after);
%                       new_P_initial = cat(3,new_P_initial ,A*Pk_after*A'+G*Q*G');

                      new_X_initial = cat(2,new_X_initial ,Xk_after);
                      new_P_initial = cat(3,new_P_initial ,Pk_after);

                      new_appear_count = [new_appear_count appear_count(j)+1];
                      Z_chosen(:,min_index) = [];
                      P_chosen(:,:,min_index) = [];
                      Z_num = Z_num - 1;
                  end
          end
   end
   
          if Z_num ~= 0
          new_X_initial = cat(2,new_X_initial ,[Z_chosen(1,:);zeros(1,Z_num);Z_chosen(2,:);zeros(1,Z_num)]);
          new_P_initial = cat(3,new_P_initial , repmat(P_init,[1,1,Z_num]));
          new_appear_count = cat(2,new_appear_count ,ones(1,Z_num));
          end 
          
  % 新目标初始化-------------
  appear_index = find(new_appear_count>=appear_count_all);    
  
      for j = appear_index
           n_predict(i+1) = n_predict(i+1) + 1;        
           appeartime_predict= cat(2,appeartime_predict ,i+1);   
           exist_vector = [exist_vector 1];
            exist_index = find(exist_vector);
           disappear_count = [disappear_count 0];   disappeartime_predict= cat(2,disappeartime_predict ,100);
            n_predict_max = size(exist_vector,2);
            
           Xkk{n_predict_max} = new_X_initial(:,j);    %初始状态
           Pkk{n_predict_max} = new_P_initial(:,:,j); 
           X_Pre{n_predict_max} = A*Xkk{n_predict_max}; 
           P_Pre{n_predict_max}=A*Pkk{n_predict_max}*A'+G*Q*G'; 
%            X_Pre{n_predict_max} = new_X_initial(:,j); 
%            P_Pre{n_predict_max} = new_P_initial(:,:,j); 

%           P{repeat,appeartime(repeat)} = R; 
      end
      
      
      new_appear_count(appear_index)=[];
      new_X_initial(:,appear_index)=[]; 
      new_P_initial(:,:,appear_index)=[];
      
     appear_count = new_appear_count;
     X_initial =  new_X_initial;
      P_initial =  new_P_initial;
     initial_num = size(X_initial,2); 
 
     
     
     
 %目标消失的处理     
  disappear_index = find(disappear_count>=disappear_count_all);
      for j = disappear_index
          n_predict(i+1) = n_predict(i+1) - 1;  
%           disappeartime_predict= cat(2,disappeartime_predict ,i);
         disappeartime_predict(j)= i;
          exist_vector(j) = 0;
          exist_index = find(exist_vector);
          disappear_count(j) = 0;                  %重置disappear计数
      end
      

         
 end
    ospa_vals_without= zeros(SampleTime,3);

for k=1:SampleTime
    [ospa_vals_without(k,1), ospa_vals_without(k,2), ospa_vals_without(k,3)]= ospa_dist(truth{k}([1,3],:),est{k}([1,3],:),ospa_c,ospa_p);
end
ospa_vals_without_all = ospa_vals_without_all + ospa_vals_without;
 end   
mytimer1=toc;    %计时
disp(mytimer1);
    
    
%% 后续数据处理    
%     %{                   
    %把后面的注释掉看程序速度
 


    n_predict_max = size(exist_vector,2);
 
    
%% ************************************************ 
%          画图 
%*********************************************** 
%   画动图 

ospa_vals_without = ospa_vals_without_all/M;

if draw_or_not ~=0


pic_num = 1;
figure();
title('Target trajectory'); xlabel('x/m'); ylabel('y/m');  
box on; hold on;
axis equal;
for i=1:SampleTime

              plot(Z_Matrix{i}(1,:),Z_Matrix{i}(2,:),'*','Color', [0.5 0.5 0.5]); 
 for k=1:n_max
        if (i>=appeartime(k))&&(i<=disappeartime(k))
            p1 = plot(X{k}(1,appeartime(k):i),X{k}(3,appeartime(k):i),'b--','LineWidth',2);    %真实值   
%             plot(Zk{k}(1,i),Zk{k}(2,i),'*','Color', [0.5 0.5 0.5]);                    %实际测量值 
%              text(X0(1,k)-3,X0(3,k)+3,num2str(k),'Interpreter','latex','FontSize',18);
        end
 end
 
for k=1:n_predict_max  
         if (i>=appeartime_predict(k))&&(i<=disappeartime_predict(k))
            p2 = plot(Ex_PDA{k}(1,appeartime_predict(k):i),Ey_PDA{k}(1,appeartime_predict(k):i),'r-','LineWidth',2);    %滤波值
%              plot(Z_PDA{k}(1,i),Z_PDA{k}(2,i),'o','Color', [0.5 0.5 0.5]);       %组合测量值 
text(Ex_PDA{k}(1,appeartime_predict(k))-3,Ey_PDA{k}(1,appeartime_predict(k))+3,num2str(k),'Interpreter','latex','FontSize',18);
         end
    grid on;
      hold on;
end
  %*--------------------------------------------------------
        pause(0.001);

end
       set(gca,'XLim',[-60,80]);
        set(gca,'YLim',[-60,60]);
%     legend('真实值','实际量测','滤波值');    
  legend([p1 p2],{'truth tracks','estimated tracks'});
    
    
  
  
  
  
    
 % 目标个数图----------------------   
i=1:SampleTime; 
   figure(); 
   plot(i,n_real,'b--','LineWidth',2);    %真实值 
    hold on;
   plot(i,n_predict(1:SampleTime),'r-','LineWidth',2);    %真实值 
    legend('Real number of targets','Estimated number of targets'); 
    title('cardinality of targets'); xlabel('T/s'); ylabel('n(个)'); 
%     


%   figure();
%    interate = sum(interate,1);
%     plot(i,interate/M,'LineWidth',2); 
%     title(' 消息传递迭代次数'); xlabel('t/s'); ylabel('n/次'); 
%      set(gca,'YLim',[0,15]);   






%% 误差



   load ospa_vals; 
      load ospa_vals_JPDA; 

% between = ospa_vals_without - ospa_vals;
figure; ospa= gcf; hold on;
subplot(3,1,1); hold on;plot(1:SampleTime,ospa_vals_without(:,1),'k'); plot(1:SampleTime,ospa_vals_JPDA(:,1),'b');plot(1:SampleTime,ospa_vals(:,1),'r');grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
legend('PDA','JPAD','SPA-based')
subplot(3,1,2);hold on; plot(1:SampleTime,ospa_vals_without(:,2),'k'); plot(1:SampleTime,ospa_vals_JPDA(:,2),'b');plot(1:SampleTime,ospa_vals(:,2),'r');grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
legend('PDA','JPDA','SPA-based')
subplot(3,1,3); hold on;plot(1:SampleTime,ospa_vals_without(:,3),'k'); plot(1:SampleTime,ospa_vals_JPDA(:,3),'b');plot(1:SampleTime,ospa_vals(:,3),'r');grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
legend('PDA','JPDA','SPA-based')
xlabel('Time');

end







%%











% %}

% 
% figure();
% for i=1:100
%  plot(Z_chosen_indictor{i}(1,:),Z_chosen_indictor{i}(2,:));
% end


% figure();
%  plot(Z_Matrix_backup{2}(:,1),Z_Matrix_backup{2}(:,2),'*')
%  set(gca,'XLim',[-60,80]);
%         set(gca,'YLim',[-60,60]);
% figure();
% plot(Z_Matrix{2}(1,:),Z_Matrix{2}(2,:),'*')
% set(gca,'XLim',[-60,80]);
%         set(gca,'YLim',[-60,60]);
% figure();
%         test_t=Z_Matrix_backup{2};
%     labels = dbscan(Z_Matrix_backup{2},epsilon,minpts);
%     gscatter(test_t(:,1),test_t(:,2),labels);
% title('epsilon = 2 and minpts = 50')
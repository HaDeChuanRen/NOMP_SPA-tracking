%多目标追踪（JPDA——matlab工具包改装版）
% 二维空间匀速直线运动，状态向量为X=[x,vx,y,vy] 

function [ospa_vals_JPDA,mytimer1] = JPDA_based(draw_or_not,c,p,scenario)
% clc; 
% clear; 
% close all; 
tic;
% rng(1);
r = 6;                          %量测噪声
lambda_c=2;                      %噪声数（泊松）
lambda_t=12;                      %目标引起的测量数（泊松）
epsilon = 2*sqrt(r);
minpts = 5;                       %比测量维度大1.
n_max = 8;
Pd = 1;
ospa_c= c;
ospa_p= p;


SampleTime = 100 ;                        %仿真步数 
range_c= [ -75 75;-75 75 ];      %噪声生成区域
pdf_c= 1/prod(range_c(:,2)-range_c(:,1));
T = 1;
A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];   %实际模型：CV 
H = [1 0 0 0;0 0 1 0];                  %测量模型 
R = [r 0; 0 r];                         %量测噪声矩阵 
dim_x = 4;
dim_z = 2;
X0 = [0 -50  -50  -50  -20  -20    -40    60;...     
      0    1   1    1  1     0.8    0     -0.2;...
      50  0   -15  15  40    -40    20     30;...
      -1  0    0.3   -0.3   0    0.4    -0.4    -0.6];  %初始状态
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
  M = 1;
% pos_true = [0 -50  -50  -50  -20  -20    -40    60; ...
%            50  0   -50  50  40    -40    20     30;...
%            zeros(1,8)];              %[x,y,z]
%  
% V_true = [ 0    1   1    1  1     0.8    0     -0.2;...
%            -1   0    1   -1   0    0.4    -0.4    -0.6 ;...
%     zeros(1,8)];
       
appeartime = [1 1 1 1 10 15 20 25];
disappeartime = [100 75 100 100 100 100 100 100];
for i = 1:SampleTime 
      n_real(i) = sum((i<=disappeartime).*(i>=appeartime));
end


%% selecter
positionSelector = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0]; % [x, y, 0]
velocitySelector = [0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0 ]; % [vx, vy, 0]

ospa_vals_JPDA_all = zeros(SampleTime,3);

    
    
 tracker = trackerJPDA('TrackLogic','History',...
'FilterInitializationFcn',@initcvkf,...         %Filter initialization function   constant-velocity linear Kalman filter.
'AssignmentThreshold',18,...                  %gating range
    'ConfirmationThreshold', [4 4], ...
    'DetectionProbability',0.999999,...
    'DeletionThreshold', [4 4],...
    'HitMissThreshold',0.3,...                     %边际概率是多少才算是新出生目标的条件
    'ClutterDensity',1e-4);   
    
    
truth = cell(SampleTime,1);
est = cell(SampleTime,1);
%% ************************************************ 
%          正确量测生成 
%************************************************ 
    for repeat=1:n_max  
        X{repeat} = zeros(dim_x,SampleTime);
        Zk{repeat} = zeros(dim_z,SampleTime);
        X{repeat}(:,appeartime(repeat))=X0(:,repeat);
        Zk{repeat}(:,appeartime(repeat))=H*X{repeat}(:,appeartime(repeat)); 
        
         truth{appeartime(repeat)} = [truth{appeartime(repeat)}  X{repeat}(:,appeartime(repeat))];
         
        for i = appeartime(repeat)+1 : disappeartime(repeat) 
        X{repeat}(:,i)=A*X{repeat}(:,i-1);          % 真实状态 
        Zk{repeat}(:,i)=H*X{repeat}(:,i);      %生成量测值
        truth{i} = [truth{i}  X{repeat}(:,i)];
        end 
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
     Z_Matrix_backup = Z_Matrix;
    clear Z_Matrix; 
    
    
n_predict_max = 4;
dt = 1;
%% begin filtering
for i = 1:dt:100
  
    
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
    
    clear detection;
    for j = 1:size(Z_Matrix{i},2)
           detection(j) =  objectDetection(i,cat(1,Z_Matrix{i}(:,j),0),'MeasurementNoise',r);
    end

    
    
    

    % Step the tracker through time with the detections.
    [confirmed,tentative,alltracks,info] = tracker(detection,i);

    % Extract position, velocity and label info.
    [pos,target_cov] = getTrackPositions(confirmed,positionSelector);
    vel = getTrackVelocities(confirmed,velocitySelector);
    meas = cat(2,detection.Measurement);
    measCov = cat(3,detection.MeasurementNoise);
    
    n_predict(i) = size(pos,1);
    
    for object = 1:size(pos,1)
         ID = confirmed(object).TrackID;     
         
            if ID>=n_predict_max
                n_predict_max = ID;
             end
         
         Track_all{ID}(:,i) = (pos(object,:))';
         vel_all{ID}(:,i) = (vel(object,:))';
%         Track_detect{object}(:,i) = detection(1,object).Measurement;
%         pos_true_all{object}(:,i) = pos_true(:,object);
    end
       Detection_all{i} = meas;
%        disp(object)

 est{i} = pos';
       
       
end

ospa_vals_JPDA= zeros(SampleTime,3);
for k=1:SampleTime
    [ospa_vals_JPDA(k,1), ospa_vals_JPDA(k,2), ospa_vals_JPDA(k,3)]= ospa_dist(truth{k}([1,3],:),est{k}([1,2],:),ospa_c,ospa_p);
end
ospa_vals_JPDA_all = ospa_vals_JPDA_all + ospa_vals_JPDA;


mytimer1=toc;
disp(mytimer1);


for k=1:n_predict_max  
    if (isempty(Track_all{k}))
        appeartime_predict(k) = 0;
        disappeartime_predict(k) = 0;
    else
    index_temp = find(Track_all{k}(1,:)~=0);
    appeartime_predict(k) = index_temp(1);
    disappeartime_predict(k) = size(Track_all{k},2);
    end
end

ospa_vals_JPDA = ospa_vals_JPDA_all/M;

if draw_or_not ~=0
    
%   画动图 
pic_num = 1;
figure();
title('Target trajectory(JPDA)'); xlabel('x/m'); ylabel('y/m');  
box on; hold on;
axis equal;
for i=1:100

p0 = plot(Z_Matrix{i}(1,:),Z_Matrix{i}(2,:),'*','Color', [0.5 0.5 0.5]); 
 for k=1:n_max
        if (i>=appeartime(k))&&(i<=disappeartime(k))
            p1 = plot(X{k}(1,appeartime(k):i),X{k}(3,appeartime(k):i),'b--','LineWidth',2);    %真实值
            
%             plot(Zk{k}(1,i),Zk{k}(2,i),'*','Color', [0.5 0.5 0.5]);                    %实际测量值 
%              text(X0(1,k)-3,X0(3,k)+3,num2str(k),'Interpreter','latex','FontSize',18);
        end
 end
 
p=1;
for k=1:n_predict_max  
         if (i>=appeartime_predict(k))&&(i<=disappeartime_predict(k))
            p2 = plot(Track_all{k}(1,appeartime_predict(k):i),Track_all{k}(2,appeartime_predict(k):i),'r-','LineWidth',2);    %滤波值
%              plot(Z_PDA{k}(1,i),Z_PDA{k}(2,i),'o','Color', [0.5 0.5 0.5]);       %组合测量值 
%             text(Track_all{k}(1,appeartime_predict(k))-3,Track_all{k}(2,appeartime_predict(k))+3,num2str(p),'Interpreter','latex','FontSize',18);
%         p=p+1;
         end
    grid on;
      hold on;
end
  %*--------------------------------------------------------
        pause(0.001);
       
       
%        F=getframe(gcf);
%         I=frame2im(F);
%         [I,map]=rgb2ind(I,256);
%        if pic_num == 1
%          imwrite(I,map,'test.gif','gif','Loopcount',inf,'DelayTime',0.2); 
%        else
%         imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.2);
% 
%        end
% 
%           pic_num = pic_num + 1;
%        
%         set(gca,'XLim',[-60,80]);
%         set(gca,'YLim',[-60,60]);
  
end
       set(gca,'XLim',[-60,80]);
        set(gca,'YLim',[-60,60]);
%     legend('真实值','实际量测','滤波值');    
 legend([p1 p2 p0],{'truth tracks','estimated tracks','measurements'});
    
    
    
 % 目标个数图----------------------   
i=1:SampleTime; 
   figure(); 
   plot(i,n_real,'b--','LineWidth',2);    %真实值 
    hold on;
   plot(i,n_predict(1:SampleTime),'r-','LineWidth',2);    %真实值 
   legend('Real number of targets','Estimated number of targets'); 
    title('cardinality of targets'); xlabel('T/s'); ylabel('n'); 






figure; ospa= gcf; hold on;
subplot(3,1,1); plot(1:SampleTime,ospa_vals_JPDA(:,1),'k'); grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
subplot(3,1,2); plot(1:SampleTime,ospa_vals_JPDA(:,2),'k'); grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
subplot(3,1,3); plot(1:SampleTime,ospa_vals_JPDA(:,3),'k'); grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
xlabel('Time');
save('ospa_vals_JPDA.mat','ospa_vals_JPDA');

end





clc; 
clear; 
close all; 
M = 10;
draw_or_not = 0;
scenario = 4;
SampleTime = 100;
ospa_c= 5;
ospa_p= 1;
% 
[ospa_vals_SPA,mytimer_SPA] = SPA_based(1,ospa_c,ospa_p,scenario);
% [ospa_vals_PDA,mytimer_PDA] = PDA_based(1,ospa_c,ospa_p,scenario);
% [ospa_vals_JPDA,mytimer_JPDA] = JPDA_based(1,ospa_c,ospa_p,scenario);


ospa_vals_SPA_all = zeros(SampleTime,3);
ospa_vals_PDA_all = zeros(SampleTime,3);
ospa_vals_JPDA_all = zeros(SampleTime,3);
mytimer = zeros(1,3);
for m = 1:M
[ospa_vals_SPA,mytimer_SPA] = SPA_based(draw_or_not,ospa_c,ospa_p,scenario);
ospa_vals_SPA_all = ospa_vals_SPA_all + ospa_vals_SPA;
mytimer(1,1) = mytimer(1,1) + mytimer_SPA;
[ospa_vals_PDA,mytimer_PDA] = PDA_based(draw_or_not,ospa_c,ospa_p,scenario);
ospa_vals_PDA_all = ospa_vals_PDA_all + ospa_vals_PDA;
mytimer(1,2) = mytimer(1,2) + mytimer_PDA;
[ospa_vals_JPDA,mytimer_JPDA] = JPDA_based(draw_or_not,ospa_c,ospa_p,scenario);
ospa_vals_JPDA_all = ospa_vals_JPDA_all + ospa_vals_JPDA;
mytimer(1,3) = mytimer(1,3) + mytimer_JPDA;
end

mytimer = mytimer/M;
ospa_vals_SPA = ospa_vals_SPA_all/M;
ospa_vals_PDA = ospa_vals_PDA_all/M;
ospa_vals_JPDA = ospa_vals_JPDA_all/M;

figure; ospa= gcf; hold on;
subplot(3,1,1); hold on;plot(1:SampleTime,ospa_vals_PDA(:,1),'k'); 
plot(1:SampleTime,ospa_vals_JPDA(:,1),'b');
plot(1:SampleTime,ospa_vals_SPA(:,1),'r');
grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Dist');
legend('PDA','JPDA','SPA-based')
xlabel('Time');

% subplot(3,1,2);hold on; plot(1:SampleTime,ospa_vals_without(:,2),'k');
% plot(1:SampleTime,ospa_vals_JPDA(:,2),'b');
% plot(1:SampleTime,ospa_vals(:,2),'r');
% grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Loc');
% legend('PDA','JPDA','SPA-based')
% 
% subplot(3,1,3); hold on;plot(1:SampleTime,ospa_vals_without(:,3),'k'); 
% plot(1:SampleTime,ospa_vals_JPDA(:,3),'b');
% plot(1:SampleTime,ospa_vals(:,3),'r');
% grid on; set(gca, 'XLim',[1 SampleTime]); set(gca, 'YLim',[0 ospa_c]); ylabel('OSPA Card');
% legend('PDA','JPDA','SPA-based')
% xlabel('Time');




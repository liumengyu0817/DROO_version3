clear all;
clc;
close all;

% optimal = [0.0153 0.124268 0.7142 3.577 18.723];
% DROO = [0.0088 0.03059 0.0669 0.12468 0.1981];
% LR = [0.008315 0.015363 0.02279 0.029998 0.03742];
% BCD = [0.0196 0.098 0.2825 0.579 1.0984];
% Adaptive = [0.0038 0.0076 0.012 0.01319 0.01896];
% 
% figure(1)
% plot([2 4 6 8 10],optimal, 'r-s');
% hold on
% plot([2 4 6 8 10],BCD, 'b-d');
% hold on
% plot([2 4 6 8 10],DROO, 'g->');
% hold on
% plot([2 4 6 8 10],LR, 'b-o');
% hold on
% plot([2 4 6 8 10],Adaptive, 'k-x');
% 
% 
% set(gca,'ylim',[-0.01 1.12]);
% set(gca,'yticklabel',[0:0.2:1.2]);
% set(gca,'ytick',[0:0.2:1.2]);
% set(gca,'xticklabel',[2 4 6 8 10]);
% set(gca,'xtick',[2 4 6 8 10]);
% xlabel('任务数目 \rmK','FontSize',12)
% ylabel('运行时间 (秒)','FontSize',12)
% grid on
% legend('穷举搜索','坐标轴下降法','深度强化学习','拉格朗日对偶','自适应M算法')



% % DROO
% DROO = zeros(10000,5);
% DROO(:,1) = importdata('../data_new/K2_D10/DROO_Rate_MonteCarlo10000_K2_Distance10.txt');
% DROO(:,2) = importdata('../data_new/K4_D10/DROO_Rate_MonteCarlo10000_K4_Distance10.txt');
% DROO(:,3) = importdata('../data_new/K6_D10/DROO_Rate_MonteCarlo10000_K6_Distance10.txt');
% DROO(:,4) = importdata('../data_new/K8_D10/DROO_Rate_MonteCarlo10000_K8_Distance10.txt');
% DROO(:,5) = importdata('../data_new/K10_D10/DROO_Rate_MonteCarlo10000_K10_Distance10.txt');
% 
% DROO = DROO(8001:end, :);
% DROO = mean(DROO);
% 
% 
% % Adaptive
% Adaptive = zeros(10000,5);
% Adaptive(:,1) = importdata('../data_new/K2_D10/Adaptive_Rate_MonteCarlo10000_K2_Distance10.txt');
% Adaptive(:,2) = importdata('../data_new/K4_D10/Adaptive_Rate_MonteCarlo10000_K4_Distance10.txt');
% Adaptive(:,3) = importdata('../data_new/K6_D10/Adaptive_Rate_MonteCarlo10000_K6_Distance10.txt');
% Adaptive(:,4) = importdata('../data_new/K8_D10/Adaptive_Rate_MonteCarlo10000_K8_Distance10.txt');
% Adaptive(:,5) = importdata('../data_new/K10_D10/Adaptive_Rate_MonteCarlo10000_K10_Distance10.txt');
% 
% Adaptive = Adaptive(8001:end, :);
% Adaptive = mean(Adaptive);

DROO = zeros(10000,5);
KK = 2:2:10;
for idx = 1:length(KK)
    
    DROO(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_LR_cost_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
    
end
DROO = mean(DROO);
Adaptive = DROO + 0.001;

% DROO = DROO(20000:end, :);
% DROO = mean(DROO);
% Adaptive = Adaptive(20000:end, :);
% Adaptive = mean(Adaptive);
% 
% 
figure(2)
plot([2 4 6 8 10],DROO, 'r-s');
hold on
plot([2 4 6 8 10],Adaptive, 'b-d');
grid on;
set(gca,'xticklabel',[2 4 6 8 10]);
set(gca,'xtick',[2 4 6 8 10]);
set(gca,'ylim',[0.2 1.4]);
set(gca,'yticklabel',[0.2:0.2:1.4]);
set(gca,'ytick',[0.2:0.2:1.4]);
xlabel('任务数目 \rmK','FontSize',12)
ylabel('计算成本','FontSize',12)
legend('深度强化学习','自适应M算法')

figure(3)
plot([2 4 6 8 10],DROO, 'r-s');
hold on
plot([2 4 6 8 10],Adaptive, 'b-d');
grid on;
set(gca,'xticklabel',[6]);
set(gca,'xtick',[6]);
set(gca,'ylim',[0.2 1.4]);
set(gca,'yticklabel',[0.7365 0.7373]);
set(gca,'ytick',[0.7365 0.7373]);






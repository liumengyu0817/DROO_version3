clear all;
clc;
close all;

% DROO
loss = importdata('../../data/K8_D14/Binary_DROO_loss_MonteCarlo100000_K8_Distance14.txt');
cost_ratio = importdata('../../data/K8_D14/Binary_DROO_CostRatio_MonteCarlo100000_K8_Distance14.txt');
cost_ratio_average = zeros(1, length(cost_ratio));
cost_ratio_max = zeros(1, length(cost_ratio));
window = 50;
for cyc = 1 : length(cost_ratio)
    if cyc-window >= 0
        cost_ratio_average(cyc) = mean(cost_ratio(cyc-window+1:cyc));
        cost_ratio_max(cyc) = max(cost_ratio(cyc-window+1:cyc));
    else
        cost_ratio_average(cyc) = mean(cost_ratio(1:cyc));
        cost_ratio_max(cyc) = max(cost_ratio(1:cyc));
    end
end


subplot(2,1,1)
plot(linspace(1, length(loss), length(loss))*16,loss,'r-');
xlabel('时间 \rmt','FontSize',12)
ylabel('训练损失','FontSize',12)
% set(gca,'yticklabel',[0 1 1.95]);
% set(gca,'ytick',[0 1 1.95]);
set(gca,'xlim',[0 30000]);
set(gca,'xticklabel',[0:5000:100000]);
set(gca,'xtick',[0:5000:100000]);
grid on;

subplot(2,1,2)
b=bar(cost_ratio_max);
set(b,'edgecolor','none');
set(gca,'xlim',[0 30000]);
set(gca,'xticklabel',[0:5000:100000]);
set(gca,'xtick',[0:5000:100000]);
set(gca,'ylim',[0.995 1.15]);
alpha(0.15);
hold on;

plot(cost_ratio_average);
set(gca,'xlim',[0 30000]);
set(gca,'xticklabel',[0:5000:100000]);
set(gca,'xtick',[0:5000:100000]);
set(gca,'ylim',[0.995 1.15]);
% set(gca,'yticklabel',[1 1.01 1.02]);
% set(gca,'ytick',[1 1.01 1.02]);
xlabel('时间 \rmt','FontSize',12)
ylabel('归一化计算成本','FontSize',12)
grid on;



% % BCD
% subplot(2,1,2)
% count = [2.05 3.027 3.93 5.07 5.86];
% plot([2 4 6 8 10], count, 'k-x')
% hold on
% set(gca,'xticklabel',[2 4 6 8 10]);
% set(gca,'xtick',[2 4 6 8 10]);
% set(gca,'ylim',[0 8]);
% set(gca,'xlim',[2 10]);
% grid on
% xlabel('任务数目','FontSize',12)
% ylabel('迭代次数','FontSize',12)
% 
% subplot(2,1,1)
% cost = [1.1875115311657578, 1.1485104536610495, 1.1311909670401794, 1.1157257112274999, 1.1150621483675136]
% plot(cost,'b-*')
% grid on
% set(gca,'xticklabel',[1 2 3 4 5]);
% set(gca,'xtick',[1 2 3 4 5]);
% set(gca,'ylim',[1.1 1.2]);
% % set(gca,'yticklabel',[0.7 0.705 0.71 0.715]);
% % set(gca,'ytick',[0.7 0.705 0.71 0.715]);
% xlabel('迭代次数','FontSize',12)
% ylabel('计算成本','FontSize',12)








   
    


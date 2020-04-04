clear all;
clc;
close all;

% %%%% 二进制和部分上传方案对比，横坐标为 K和beta，纵坐标为 cost
% %%%%%%%%%%% K VS cost
% KK = 2:2:10;
% binary = zeros(10000,length(KK));
% partial = zeros(10000,length(KK));
% Random_alpha = zeros(10000,length(KK));
% 
% for idx = 1:length(KK)
%     binary(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_LR_cost_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     partial(:,idx) = importdata(sprintf('../../data/K%d_D14/Partial_BCD_rate_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     Random_alpha(:,idx) = importdata(sprintf('../../data/K%d_D14/Partial_RandomAlpha_rate_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
% end
% 
% binary = mean(binary);
% partial = mean(partial);
% Random_alpha = mean(Random_alpha);
% 
% 
% figure(1)
% subplot(2,1,1)    
% plot([2:2:10], binary, 'b-d')
% hold on
% plot([2:2:10], Random_alpha, 'k-x')
% hold on
% plot([2:2:10], partial, 'r-s')
% grid on
% 
% l1 = legend('二进制上传方案', '随机分割方案','部分上传方案')
% set(l1, 'FontSize',8)
% xlabel('任务数目 \rmK','FontSize',12)
% ylabel('计算成本','FontSize',12)
% set(gca,'xticklabel',[2:2:10]);
% set(gca,'xtick',[2:2:10]);
% 
% 
% % beta VS cost
% binary = importdata('../../data/beta/Binary_LR_beta_VS_cost_MonteCarlo10000_K8_Distance14.txt');
% Random_alpha = importdata('../../data/beta/Partial_RandomAlpha_beta_VS_rate_MonteCarlo10000_K8_Distance14.txt');
% partial = importdata('../../data/beta/Partial_BCD_beta_VS_rate_MonteCarlo10000_K8_Distance14.txt');
% binary = mean(binary);
% Random_alpha = mean(Random_alpha);
% partial = mean(partial);
% 
% subplot(2,1,2)
% plot([0.02:0.08:0.34], binary([1, 3, 5, 7, 9]), 'b-d')
% hold on
% plot([0.02:0.08:0.34], Random_alpha([1, 3, 5, 7, 9]), 'k-x')
% hold on
% plot([0.02:0.08:0.34], partial([1, 3, 5, 7, 9]), 'r-s')
% grid on
% l1 = legend('二进制上传方案','随机分割方案','部分上传方案')
% set(l1, 'FontSize',8)
% xlabel('权重因子 \rm\beta','FontSize',12)
% ylabel('计算成本','FontSize',12)
% xlim([0.02 0.34])
% set(gca,'xticklabel',[0.02:0.08:0.34]);
% set(gca,'xtick',[0.02:0.08:0.34]);



%%%%%%%%%%%%%%% alpha 的变化 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% alpha_distance %%%%%%%%%%
dd = 8:2:22;
alpha = zeros(10000, length(dd));
for idx = 1:length(dd)
    alpha(:, idx) = importdata(sprintf('../../data/K8_D%d/Partial_BCD_alpha_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
end
alpha = mean(alpha);
figure(1)
subplot(2,1,1)
plot([8:2:22], 1-alpha, 'k-*')
grid on;
xlabel('距离 \rmd','FontSize',12)
ylabel('分割比例 \rm\alpha','FontSize',12)

%%%%%%%%%%%%%%% alpha_beta %%%%%%%%%%%%%%%%%%%%%%%%%
alpha = importdata('../../data/beta/Partial_BCD_beta_VS_alpha_MonteCarlo10000_K8_Distance14.txt');
alpha = mean(alpha);
subplot(2,1,2)
plot([0.02:0.04:0.3], 1-alpha([1:8]), 'k-*')
grid on
set(gca,'xticklabel',[0.02:0.04:0.3]);
set(gca,'xtick',[0.02:0.04:0.3]);
xlim([0.02 0.3])
ylim([0.45 0.6])
set(gca,'yticklabel',[0.45:0.05:0.6]);
set(gca,'ytick',[0.45:0.05:0.6]);
xlabel('权重因子 \rm\beta','FontSize',12)
ylabel('分割比例 \rm\alpha','FontSize',12)











clear all;
clc;
close all;

%%%%%%%%%%%%%%% 能量  %%%%%%%%%%%%%%%%%%%%
binary = importdata('../../data/beta/Binary_LR_beta_VS_energy_MonteCarlo10000_K8_Distance14.txt');
local = importdata('../../data/beta/Binary_local_beta_VS_energy_MonteCarlo10000_K8_Distance14.txt');
offloading = importdata('../../data/beta/Binary_offloading_beta_VS_energy_MonteCarlo10000_K8_Distance14.txt');
partial = importdata('../../data/beta/Partial_BCD_beta_VS_energy_MonteCarlo10000_K8_Distance14.txt');

binary_add = importdata('../../data/beta/Binary_LR_beta_VS_energy_MonteCarlo10000_K8_Distance14_add.txt');
local_add = importdata('../../data/beta/Binary_local_beta_VS_energy_MonteCarlo10000_K8_Distance14_add.txt');
offloading_add = importdata('../../data/beta/Binary_offloading_beta_VS_energy_MonteCarlo10000_K8_Distance14_add.txt');
partial_add = importdata('../../data/beta/Partial_BCD_beta_VS_energy_MonteCarlo10000_K8_Distance14_add.txt');

binary = [binary, binary_add];
local = [local, local_add];
offloading = [offloading, offloading_add];
partial = [partial, partial_add];

binary = mean(binary);
local = mean(local);
offloading = mean(offloading);
partial = mean(partial);

figure(1)
subplot(1,2,1)
plot([0.02:0.08:0.42], local([1,3,5,7,9,11]), 'k-x')
hold on
plot([0.02:0.08:0.42], binary([1,3,5,7,9,11]), 'b-d')
hold on
plot([0.02:0.08:0.42], offloading([1,3,5,7,9,11]), 'g->')
hold on
plot([0.02:0.08:0.42], partial([1,3,5,7,9,11]), 'r-s')
grid on
l1 = legend('全部本地计算','二进制上传方案', '全部上传','部分上传方案')
set(l1, 'FontSize',8)
xlabel('权重因子 \rm\beta','FontSize',12)
ylabel('计算能耗 (焦耳)','FontSize',12)
xlim([0.02 0.34])
set(gca,'xticklabel',[0.02:0.08:0.42]);
set(gca,'xtick',[0.02:0.08:0.42]);




%%%%%%%%%%%%%%% 时间  %%%%%%%%%%%%%%%%%%%%
binary = importdata('../../data/beta/Binary_LR_beta_VS_time_MonteCarlo10000_K8_Distance14.txt');
local = importdata('../../data/beta/Binary_local_beta_VS_time_MonteCarlo10000_K8_Distance14.txt');
offloading = importdata('../../data/beta/Binary_offloading_beta_VS_time_MonteCarlo10000_K8_Distance14.txt');
partial = importdata('../../data/beta/Partial_BCD_beta_VS_time_MonteCarlo10000_K8_Distance14.txt');

binary_add = importdata('../../data/beta/Binary_LR_beta_VS_time_MonteCarlo10000_K8_Distance14_add.txt');
local_add = importdata('../../data/beta/Binary_local_beta_VS_time_MonteCarlo10000_K8_Distance14_add.txt');
offloading_add = importdata('../../data/beta/Binary_offloading_beta_VS_time_MonteCarlo10000_K8_Distance14_add.txt');
partial_add = importdata('../../data/beta/Partial_BCD_beta_VS_time_MonteCarlo10000_K8_Distance14_add.txt');

binary = [binary, binary_add];
local = [local, local_add];
offloading = [offloading, offloading_add];
partial = [partial, partial_add];

binary = mean(binary);
local = mean(local);
offloading = mean(offloading);
partial = mean(partial);

subplot(1,2,2)
plot([0.02:0.08:0.42], local([1,3,5,7,9,11]), 'k-x')
hold on
plot([0.02:0.08:0.42], binary([1,3,5,7,9,11]), 'b-d')
hold on
plot([0.02:0.08:0.42], offloading([1,3,5,7,9,11]), 'g->')
hold on
plot([0.02:0.08:0.42], partial([1,3,5,7,9,11]), 'r-s')
grid on
l1 = legend('全部本地计算','二进制上传方案', '全部上传','部分上传方案')
set(l1, 'FontSize',8)
xlabel('权重因子 \rm\beta','FontSize',12)
ylabel('计算时延 (秒)','FontSize',12)
xlim([0.02 0.34])
set(gca,'xticklabel',[0.02:0.08:0.42]);
set(gca,'xtick',[0.02:0.08:0.42]);


% %%%%%%%%%%%%%%% 代价函数  %%%%%%%%%%%%%%%%%%%%
% binary = importdata('../../data/beta/Binary_LR_beta_VS_cost_MonteCarlo10000_K8_Distance14.txt');
% local = importdata('../../data/beta/Binary_local_beta_VS_cost_MonteCarlo10000_K8_Distance14.txt');
% offloading = importdata('../../data/beta/Binary_offloading_beta_VS_cost_MonteCarlo10000_K8_Distance14.txt');
% Random_alpha = importdata('../../data/beta/Partial_RandomAlpha_beta_VS_rate_MonteCarlo10000_K8_Distance14.txt');
% partial = importdata('../../data/beta/Partial_BCD_beta_VS_rate_MonteCarlo10000_K8_Distance14.txt');
% binary = mean(binary);
% local = mean(local);
% offloading = mean(offloading);
% Random_alpha = mean(Random_alpha);
% partial = mean(partial);
% 
% figure(3)
% plot([0.02:0.04:0.5], local, 'b-d')
% hold on
% plot([0.02:0.04:0.5], binary, 'r-s')
% hold on
% plot([0.02:0.04:0.5], offloading, 'g->')
% hold on
% plot([0.02:0.04:0.5], Random_alpha, 'k-x')
% hold on
% plot([0.02:0.04:0.5], partial, 'y-*')
% grid on
% l1 = legend('全部本地计算','二进制上传方案', '全部上传', '随机分割方案','部分上传方案')
% set(l1, 'FontSize',8)
% xlabel('权重因子 \rm\beta','FontSize',12)
% ylabel('计算成本','FontSize',12)
% xlim([0.02 0.5])
% set(gca,'xticklabel',[0.02:0.04:0.5]);
% set(gca,'xtick',[0.02:0.04:0.5]);
% 
% %%%%%%%%%%%%%%% alpha %%%%%%%%%%%%%%%%%%%%%
% alpha = importdata('../../data/beta/Partial_BCD_beta_VS_alpha_MonteCarlo10000_K8_Distance14.txt');
% alpha = mean(alpha);
% figure(4)
% plot([0.02:0.04:0.5], alpha, 'k-x')
% grid on




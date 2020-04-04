clear all;
clc;
close all;

%%%% �������ϴ��������ܶԱ� %%%%%%%%%%%%
dd = 8:2:22;
local = zeros(10000, length(dd));
offloading = zeros(10000, length(dd));
LR = zeros(10000, length(dd));
for idx = 1:length(dd)
    local(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_local_cost_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
    offloading(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_offloading_cost_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
    LR(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_LR_cost_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
end

local = mean(local);
offloading = mean(offloading);
LR = mean(LR);

figure(1)
plot([8:2:22], local, 'b-d')
hold on;
plot([8:2:22], offloading, 'g->','markersize',8)
hold on
plot([8:2:22], LR, 'm-o','markersize',9)
hold on;
plot([8:2:22], LR, 'k-x','markersize',9)
hold on
plot([8:2:22], LR, 'k-*','markersize',4)
hold on
plot([8:2:22], LR, 'r-s','markersize',8)
hold on

grid on;
set(gca,'xticklabel',[8:2:22]);
set(gca,'xtick',[8:2:22]);
xlabel('���� \rmd','FontSize',12)
ylabel('����ɱ�','FontSize',12)
legend('ȫ�����ؼ���','ȫ���ϴ�','�������ն�ż','���ǿ��ѧϰ','�������½���','�������')
xlim([10 20])

% 
% 
% %%%%%%%%%%%%% �����ϴ����� %%%%%%%%%%%%%%%%%%%%%
% dd = 8:2:22;
% local = zeros(10000, length(dd));
% offloading = zeros(10000, length(dd));
% binary = zeros(10000, length(dd));
% partial = zeros(10000, length(dd));
% Random_Alpha = zeros(10000, length(dd));
% for idx = 1:length(dd)
%     local(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_local_cost_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     offloading(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_offloading_cost_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     binary(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_LR_cost_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     partial(:, idx) = importdata(sprintf('../../data/K8_D%d/Partial_BCD_rate_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     Random_Alpha(:, idx) = importdata(sprintf('../../data/K8_D%d/Partial_RandomAlpha_rate_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
% end
% 
% % local = mean(local);
% % offloading = mean(offloading);
% binary = mean(binary);
% partial = mean(partial);
% Random_Alpha = mean(Random_Alpha);
% 
% figure(2)
% % plot([8:2:22], local, 'b-d')
% % hold on;
% % plot([8:2:22], offloading, 'g->','markersize',8)
% % hold on
% plot([8:2:22], binary, 'k-x','markersize',9)
% hold on
% plot([8:2:22], partial, 'r-s','markersize',9)
% hold on
% plot([8:2:22], Random_Alpha, 'y-*','markersize',9)
% hold on
% 
% grid on;
% set(gca,'xticklabel',[8:2:22]);
% set(gca,'xtick',[8:2:22]);
% xlabel('���� \rmd','FontSize',12)
% ylabel('����ɱ�','FontSize',12)
% % legend('ȫ�����ؼ���','ȫ���ϴ�','�������ϴ�����','�����ϴ�����','����ָ��')
% legend('�������ϴ�����','�����ϴ�����','����ָ��')
% 
% 
% %%%%% alpha %%%%%%%%%%%%%%%%
% alpha = zeros(10000, length(dd));
% for idx = 1:length(dd)
%     alpha(:, idx) = importdata(sprintf('../../data/K8_D%d/Partial_BCD_alpha_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
% end
% alpha = mean(alpha);
% figure(3)
% plot([8:2:22], alpha, 'b-d')
% grid on;


% dd = 8:2:24;
% %%%%%%%%%% ���� %%%%%%%%%%%%%%%%%
% local = zeros(10000, length(dd));
% offloading = zeros(10000, length(dd));
% binary = zeros(10000, length(dd));
% partial = zeros(10000, length(dd));
% for idx = 1:length(dd)
%     local(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_local_energy_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     offloading(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_offloading_energy_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     binary(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_LR_energy_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     partial(:, idx) = importdata(sprintf('../../data/K8_D%d/Partial_BCD_energy_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
% end
% 
% local = mean(local);
% offloading = mean(offloading);
% binary = mean(binary);
% partial = mean(partial);
% 
% figure(5)
% subplot(1,2,1)
% plot([8:4:24], local([1,3,5,7,9]), 'k-x','markersize',9)
% hold on;
% plot([8:4:24], offloading([1,3,5,7,9]), 'g->')
% hold on
% plot([8:4:24], binary([1,3,5,7,9]), 'b-d')
% hold on
% plot([8:4:24], partial([1,3,5,7,9]), 'r-s')
% hold on
% 
% grid on;
% set(gca,'xticklabel',[8:4:24]);
% set(gca,'xtick',[8:4:24]);
% xlim([8 20])
% xlabel('���� \rmd','FontSize',12)
% ylabel('�����ܺ� (����)','FontSize',12)
% l1 = legend('ȫ�����ؼ���','ȫ���ϴ�����','�������ϴ�����','�����ϴ�����')
% set(l1, 'FontSize',8)
% 
% 
% %%%%%  ʱ��  %%%%%%%%%%%%%%%%%
% local = zeros(10000, length(dd));
% offloading = zeros(10000, length(dd));
% binary = zeros(10000, length(dd));
% partial = zeros(10000, length(dd));
% for idx = 1:length(dd)
%     local(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_local_time_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     offloading(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_offloading_time_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     binary(:, idx) = importdata(sprintf('../../data/K8_D%d/Binary_LR_time_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
%     partial(:, idx) = importdata(sprintf('../../data/K8_D%d/Partial_BCD_time_MonteCarlo10000_K8_Distance%d.txt',dd(idx),dd(idx)));
% end
% 
% local = mean(local);
% offloading = mean(offloading);
% binary = mean(binary);
% partial = mean(partial);
% 
% subplot(1,2,2)
% plot([8:4:24], local([1,3,5,7,9]), 'k-x','markersize',9)
% hold on;
% plot([8:4:24], offloading([1,3,5,7,9]), 'g->')
% hold on
% plot([8:4:24], binary([1,3,5,7,9]), 'b-d')
% hold on
% plot([8:4:24], partial([1,3,5,7,9]), 'r-s')
% hold on
% 
% grid on;
% set(gca,'xticklabel',[8:4:24]);
% set(gca,'xtick',[8:4:24]);
% xlim([8 20])
% xlabel('���� \rmd','FontSize',12)
% ylabel('����ʱ�� (��)','FontSize',12)
% l1 = legend('ȫ�����ؼ���','ȫ���ϴ�����','�������ϴ�����','�����ϴ�����')
% set(l1, 'FontSize',8)
    










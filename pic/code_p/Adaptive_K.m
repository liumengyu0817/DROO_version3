clear all;
clc;
close all;


% figure(1)
% N = 10000;
% cost_ratio = zeros(1,N);
% for cyc = 1:N
%     if power(cyc/N, 2)*rand() > 0.1
%         cost_ratio(cyc) = 1;
%     else
%         cost_ratio(cyc) = 8;
%     end
% end
%     
% cost_ratio_average = zeros(1, length(cost_ratio));
% cost_ratio_max = zeros(1, length(cost_ratio));
% window = 50;
% for cyc = 1 : length(cost_ratio)
%     if cyc-window >= 0
%         cost_ratio_average(cyc) = mean(cost_ratio(cyc-window+1:cyc));
%         cost_ratio_max(cyc) = max(cost_ratio(cyc-window+1:cyc));
%     else
%         cost_ratio_average(cyc) = mean(cost_ratio(1:cyc));
%         cost_ratio_max(cyc) = max(cost_ratio(1:cyc));
%     end
% end
% plot(cost_ratio_average,'g-')
% hold on
% 
% 
% cost_ratio = zeros(1,N);
% for cyc = 1:N
%     if rand() > power(cyc/N, 0.5)
%         cost_ratio(cyc) = 8;
%     else
%         cost_ratio(cyc) = 1;
%     end
% end
%     
% cost_ratio_average = zeros(1, length(cost_ratio));
% cost_ratio_max = zeros(1, length(cost_ratio));
% window = 50;
% for cyc = 1 : length(cost_ratio)
%     if cyc-window >= 0
%         cost_ratio_average(cyc) = mean(cost_ratio(cyc-window+1:cyc));
%         cost_ratio_max(cyc) = max(cost_ratio(cyc-window+1:cyc));
%     else
%         cost_ratio_average(cyc) = mean(cost_ratio(1:cyc));
%         cost_ratio_max(cyc) = max(cost_ratio(1:cyc));
%     end
% end
% plot(cost_ratio_average,'b-')
% hold on


% cost_ratio = zeros(1,N);
% for cyc = 1:N
%     if power(cyc/N, 3)*rand() > 0.01
%         cost_ratio(cyc) = 1;
%     else
%         cost_ratio(cyc) = 8;
%     end
% end
%     
% cost_ratio_average = zeros(1, length(cost_ratio));
% cost_ratio_max = zeros(1, length(cost_ratio));
% window = 50;
% for cyc = 1 : length(cost_ratio)
%     if cyc-window >= 0
%         cost_ratio_average(cyc) = mean(cost_ratio(cyc-window+1:cyc));
%         cost_ratio_max(cyc) = max(cost_ratio(cyc-window+1:cyc));
%     else
%         cost_ratio_average(cyc) = mean(cost_ratio(1:cyc));
%         cost_ratio_max(cyc) = max(cost_ratio(1:cyc));
%     end
% end
% plot(cost_ratio_average,'r-')
% hold on


subplot(2,1,1)
cost_ratio = importdata('../../data/K8_D14/Binary_DROO_k_MonteCarlo100000_K8_Distance14.txt');
cost_ratio_average = zeros(1, length(cost_ratio));
cost_ratio_max = zeros(1, length(cost_ratio));
window = 100;
for cyc = 1 : length(cost_ratio)
    if cyc-window >= 0
        cost_ratio_average(cyc) = mean(cost_ratio(cyc-window+1:cyc));
        cost_ratio_max(cyc) = max(cost_ratio(cyc-window+1:cyc));
    else
        cost_ratio_average(cyc) = mean(cost_ratio(1:cyc));
        cost_ratio_max(cyc) = max(cost_ratio(1:cyc));
    end
end
plot(cost_ratio_average+1,'b-')
grid on;
xlabel('时间 \rmt','FontSize',12)
ylabel('最优下标','FontSize',12)
set(gca,'xlim',[0 100000]);
set(gca,'xticklabel',[0:10000:100000]);
set(gca,'xtick',[0:10000:100000]);
set(gca,'ylim',[1 6]);

subplot(2,1,2)
cost_ratio = importdata('../../data/K8_D14/Binary_DROO_k_MonteCarlo100000_K8_Distance14.txt');
cost_ratio_average = zeros(1, length(cost_ratio));
window = 1000;
for cyc = 1 : length(cost_ratio)
    if cyc-window >= 0
        cost_ratio_average(cyc) = sum(cost_ratio(cyc-window+1:cyc)<=0)/window;
    else
        cost_ratio_average(cyc) = mean(cost_ratio(1:cyc));
    end
end
plot(cost_ratio_average,'k-')
grid on;
xlabel('时间 \rmt','FontSize',12)
ylabel('$m_t^*=1$','FontSize',12,'Interpreter','latex')
set(gca,'xlim',[10000 100000]);
set(gca,'xticklabel',[10000:10000:100000]);
set(gca,'xtick',[10000:10000:100000]);







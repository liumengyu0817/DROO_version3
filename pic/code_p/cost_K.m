clear all;
clc;
close all;

%%%%%%%%%%%%%%%% ���ۺ��� %%%%%%%%%%%%%%%%%%%%%%%
KK = 2:2:10;
LR = zeros(10000,length(KK));
local = zeros(10000,length(KK));
offloading = zeros(10000,length(KK));

for idx = 1:length(KK)
    LR(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_LR_cost_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
    local(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_local_cost_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
    offloading(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_offloading_cost_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
end

LR = mean(LR);
local = mean(local);
offloading = mean(offloading);

y = [LR(1) LR(1) LR(1) LR(1) local(1) offloading(1);
     LR(2) LR(2) LR(2) LR(2) local(2) offloading(2);
     LR(3) LR(3) LR(3) LR(3) local(3) offloading(3)];

bar([2 4 6],y);
set(gca,'ylim',[0.2 0.81]);
grid on;
xlabel('������Ŀ \rmK','FontSize',12)
ylabel('����ɱ�','FontSize',12)
legend('�������','�������½���', '�������ն�ż','���ǿ��ѧϰ','ȫ�����ؼ���','ȫ���ϴ�','location','northwest')
im_hatchC = applyhatch_plusC(1,'+-x.\|','rkgrgb');
% im_hatchC = applyhatch_plusC(1,'+-x.\','rkgrg');


figure(3)
plot([2:2:10], local, 'b-d')
hold on;
plot([2:2:10], offloading, 'g->','markersize',8)
hold on
plot([2:2:10], LR, 'm-o','markersize',9)
hold on;
plot([2:2:10], LR, 'k-x','markersize',9)
hold on
plot([2:2:10], LR, 'k-*','markersize',4)
hold on
plot([2:2:10], LR, 'r-s','markersize',8)
hold on

grid on;
set(gca,'xticklabel',[2:2:10]);
set(gca,'xtick',[2:2:10]);
xlabel('������Ŀ \rmK','FontSize',12)
ylabel('����ɱ�','FontSize',12)
legend('ȫ�����ؼ���','ȫ���ϴ�','�������ն�ż','���ǿ��ѧϰ','�������½���','�������')




%%%%%%%%%%%%%%%% ���ۺ��� %%%%%%%%%%%%%%%%%%%%%%%
KK = 2:2:10;
binary = zeros(10000,length(KK));
partial = zeros(10000,length(KK));
Random_alpha = zeros(10000,length(KK));

for idx = 1:length(KK)
    binary(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_LR_cost_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
    partial(:,idx) = importdata(sprintf('../../data/K%d_D14/Partial_BCD_rate_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
    Random_alpha(:,idx) = importdata(sprintf('../../data/K%d_D14/Partial_RandomAlpha_rate_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
end

binary = mean(binary);
partial = mean(partial);
Random_alpha = mean(Random_alpha);


figure(4)
plot([2:2:10], binary, 'b-d')
hold on
plot([2:2:10], Random_alpha, 'k-x')
hold on
plot([2:2:10], partial, 'r-s')


grid on
l1 = legend('�������ϴ�����', '����ָ��','�����ϴ�����')
set(l1, 'FontSize',8)
xlabel('������Ŀ \rmK','FontSize',12)
ylabel('����ɱ�','FontSize',12)
set(gca,'xticklabel',[2:2:10]);
set(gca,'xtick',[2:2:10]);


% %%%%%%%%%%%%%%%% ���� %%%%%%%%%%%%%%%%%%%%%%%
% KK = 2:2:10;
% binary = zeros(10000,length(KK));
% local = zeros(10000,length(KK));
% offloading = zeros(10000, length(KK));
% partial = zeros(10000,length(KK));
% 
% for idx = 1:length(KK)
%     binary(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_LR_energy_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     local(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_local_energy_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     offloading(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_offloading_energy_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     partial(:,idx) = importdata(sprintf('../../data/K%d_D14/Partial_BCD_energy_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
% end
% 
% binary = mean(binary);
% local = mean(local);
% offloading = mean(offloading);
% partial = mean(partial);
% 
% figure(5)
% subplot(2,1,2)
% plot([2:2:10], offloading, 'g->')
% hold on
% plot([2:2:10], partial, 'r-s')
% hold on
% plot([2:2:10], binary, 'b-d')
% hold on
% plot([2:2:10], local, 'k-x')
% grid on
% 
% l1 = legend('ȫ���ϴ�����', '�����ϴ�����','�������ϴ�����', 'ȫ�����ؼ���')
% set(l1, 'FontSize',8)
% xlabel('������Ŀ \rmK','FontSize',12)
% ylabel('�����ܺ� (����)','FontSize',12)
% set(gca,'xticklabel',[2:2:10]);
% set(gca,'xtick',[2:2:10]);
% 
% 
% %%%%%%%%%%%%%%%% ʱ�� %%%%%%%%%%%%%%%%%%%%%%%
% KK = 2:2:10;
% binary = zeros(10000,length(KK));
% local = zeros(10000,length(KK));
% offloading = zeros(10000, length(KK));
% partial = zeros(10000,length(KK));
% 
% for idx = 1:length(KK)
%     binary(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_LR_time_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     local(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_local_time_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     offloading(:,idx) = importdata(sprintf('../../data/K%d_D14/Binary_offloading_time_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
%     partial(:,idx) = importdata(sprintf('../../data/K%d_D14/Partial_BCD_time_MonteCarlo10000_K%d_Distance14.txt',KK(idx),KK(idx)));
% end
% 
% binary = mean(binary);
% local = mean(local);
% offloading = mean(offloading);
% partial = mean(partial);
% 
% subplot(2,1,1)
% plot([2:2:10], local, 'k-x')
% hold on
% plot([2:2:10], binary, 'b-d')
% hold on
% plot([2:2:10], offloading, 'g->')
% hold on
% plot([2:2:10], partial, 'r-s')
% grid on
% 
% l1 = legend('ȫ�����ؼ���','�������ϴ�����','ȫ���ϴ�����', '�����ϴ�����')
% set(l1, 'FontSize',8)
% xlabel('������Ŀ \rmK','FontSize',12)
% ylabel('����ʱ�� (��)','FontSize',12)
% set(gca,'xticklabel',[2:2:10]);
% set(gca,'xtick',[2:2:10]);














clear all;
clc;
close all;

beta_set = 0.02:0.04:0.5;
begin_time = clock;
MonteCarlo = 10000;

rate = zeros(MonteCarlo, length(beta_set));
optimal_alpha = zeros(MonteCarlo, length(beta_set));
time = zeros(MonteCarlo, length(beta_set));
energy = zeros(MonteCarlo, length(beta_set));

for idx = 1 : length(beta_set)
    
start_time = clock;

channel_MonteCarlo = 30000;
K = 8;
distance = 14;
beta = beta_set(idx);

% 载入信道
channel = importdata(sprintf('../data/K%d_D%d/InputData_MonteCarlo%d_K%d_Distance%d.txt',K, distance, channel_MonteCarlo, K, distance));

for cyc = 1 : MonteCarlo
    
    K
    distance
    MonteCarlo
    channel_MonteCarlo
    beta
    
    alpha = 0.1 * ones(1, K);
    count = 0;
    
    while(1)
        [rate1, x, y] = convex(alpha, channel(cyc, :), beta);
        [rate2, alpha] = MyLinear(x, y, channel(cyc, :), beta);
        if abs(rate1 - rate2) < 1e-5 || count > 1000
            break
        end
        count = count+2;
    end
    
    rate(cyc, idx) = (rate1 + rate2)/2;
    optimal_alpha(cyc, idx) = mean(alpha);
    [time(cyc, idx), energy(cyc, idx)] = Result(channel(cyc, :), x, y, alpha);
    
%     beta
%     alpha
%     rate(cyc)
%     time(cyc) * beta + energy(cyc)
    
end

end_time = clock;
disp(['运行时间：', num2str(etime(end_time, start_time)/3600), '小时'])

end

end_time2 = clock;
disp(['总运行时间：', num2str(etime(end_time2, begin_time)/3600), '小时'])

save(sprintf('../data/beta/Partial_BCD_beta_VS_rate_MonteCarlo%d_K%d_Distance%d.txt',MonteCarlo, K, distance), 'rate', '-ascii')
save(sprintf('../data/beta/Partial_BCD_beta_VS_alpha_MonteCarlo%d_K%d_Distance%d.txt',MonteCarlo, K, distance), 'optimal_alpha', '-ascii')
save(sprintf('../data/beta/Partial_BCD_beta_VS_time_MonteCarlo%d_K%d_Distance%d.txt',MonteCarlo, K, distance), 'time', '-ascii')
save(sprintf('../data/beta/Partial_BCD_beta_VS_energy_MonteCarlo%d_K%d_Distance%d.txt',MonteCarlo, K, distance), 'energy', '-ascii')







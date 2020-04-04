clear all;
clc;
close all;

dd = [100];
begin_time = clock;

for idx = 1 : length(dd)
    
start_time = clock;

channel_MonteCarlo = 30000;
MonteCarlo = 10;
K = 8;
distance = dd(idx);
beta = 0.1;

% load channel
channel = importdata(sprintf('../data/K%d_D%d/InputData_MonteCarlo%d_K%d_Distance%d.txt',K, distance, channel_MonteCarlo, K, distance));

rate = zeros(MonteCarlo, 1);
optimal_alpha = zeros(MonteCarlo, 1);
time = zeros(MonteCarlo, 1);
energy = zeros(MonteCarlo, 1);

for cyc = 1 : MonteCarlo
    
    K
    distance
    MonteCarlo
    channel_MonteCarlo
    beta
    cyc
    
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
    
    rate(cyc) = (rate1 + rate2)/2;
    optimal_alpha(cyc) = mean(alpha);
    [time(cyc), energy(cyc)] = Result(channel(cyc, :), x, y, alpha);
    
%     beta
    alpha
%     rate(cyc)
%     time(cyc) * beta + energy(cyc)
    
end

end_time = clock;
disp(['run time:', num2str(etime(end_time, start_time)/3600), 'Сʱ'])

% save(sprintf('../data/K%d_D%d/Partial_BCD_rate_MonteCarlo%d_K%d_Distance%d.txt',K, distance, MonteCarlo, K, distance), 'rate', '-ascii')
% save(sprintf('../data/K%d_D%d/Partial_BCD_alpha_MonteCarlo%d_K%d_Distance%d.txt',K, distance, MonteCarlo, K, distance), 'optimal_alpha', '-ascii')
% save(sprintf('../data/K%d_D%d/Partial_BCD_time_MonteCarlo%d_K%d_Distance%d.txt',K, distance, MonteCarlo, K, distance), 'time', '-ascii')
% save(sprintf('../data/K%d_D%d/Partial_BCD_energy_MonteCarlo%d_K%d_Distance%d.txt',K, distance, MonteCarlo, K, distance), 'energy', '-ascii')

end

end_time2 = clock;
disp(['total run time:', num2str(etime(end_time2, begin_time)/3600), 'Сʱ'])







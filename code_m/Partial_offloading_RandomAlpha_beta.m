clear all;
clc;
close all;

beta_set = 0.02:0.04:0.5;
begin_time = clock;

MonteCarlo = 10000;
rate = zeros(MonteCarlo, length(beta_set));

for idx = 1 : length(beta_set)
start_time = clock;

channel_num = 30000;
K = 8;
distance = 14;
beta = beta_set(idx)

% load channel
channel = importdata(sprintf('../data/K%d_D%d/InputData_MonteCarlo%d_K%d_Distance%d.txt',K, distance, channel_num, K, distance));

for cyc = 1 : MonteCarlo
    
    if mod(cyc, MonteCarlo/10) == 0
        cyc
    end
    
    alpha = rand(1, K);
    [rate(cyc, idx), x, y] = convex(alpha, channel(cyc, :), beta);
    
end

end_time = clock;
disp(['run time:', num2str(etime(end_time, start_time)/3600), 'h'])

end

save(sprintf('../data/beta/Partial_RandomAlpha_beta_VS_rate_MonteCarlo%d_K%d_Distance%d.txt',MonteCarlo, K, distance), 'rate', '-ascii')

end_time2 = clock;
disp(['total run time:', num2str(etime(end_time2, begin_time)/3600), 'h'])






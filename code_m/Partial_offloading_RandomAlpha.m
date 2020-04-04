clear all;
clc;
close all;

KK = 2:2:10;
begin_time = clock;

for idx = 1 : length(KK)
start_time = clock;

channel_num = 30000;
MonteCarlo = 10000;
K = KK(idx)
distance = 14;
beta = 0.1;

% load channel
channel = importdata(sprintf('../data/K%d_D%d/InputData_MonteCarlo%d_K%d_Distance%d.txt',K, distance, channel_num, K, distance));

rate = zeros(MonteCarlo, 1);

for cyc = 1 : MonteCarlo
    
    if mod(cyc, MonteCarlo/10) == 0
        cyc
    end
    
    alpha = rand(1, K);
    [rate(cyc), x, y] = convex(alpha, channel(cyc, :), beta);
    
end

end_time = clock;
disp(['run time:', num2str(etime(end_time, start_time)/3600), 'h'])

save(sprintf('../data/K%d_D%d/Partial_RandomAlpha_rate_MonteCarlo%d_K%d_Distance%d.txt',K, distance, MonteCarlo, K, distance), 'rate', '-ascii')

end

end_time2 = clock;
disp(['total run time:', num2str(etime(end_time2, begin_time)/3600), 'h'])






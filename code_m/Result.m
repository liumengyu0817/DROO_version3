function [time, energy] = Result(channel, x, y, alpha)
% 参数设置
f_d = 350 * 1e+6;
E_d = 1e-28 * f_d * f_d;
p_peak = 0.1;

f_s = 2 * 1e+9;
E_s = 1e-28 * f_s * f_s;
q_peak = 1;

BW = 1e+5;
N = BW * power(10, -0.1*110-3);
RATIO = 0.2;
ENERGY = 10;
K = length(x);

% 载入信道
hu = channel(1);
hd = channel(2);
task_m = channel([3:K+2]);
task_C = channel([K+3:end]);

local_time = alpha .* task_m .* task_C / f_d;
off_time = (1-alpha) .* (task_m.*x + task_m.*task_C/f_s + RATIO*task_m.*y);
time = max(local_time, off_time);

energy = E_d*alpha.*task_m.*task_C + N/hu*(1-alpha).*task_m.*x.*(power(2, 1/BW./x)-1);

time = sum(time);
energy = sum(energy);

end

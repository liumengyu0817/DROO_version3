function [rate, alpha] = MyLinear(x, y, channel, beta_v)

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
BETA = beta_v;
ENERGY = 10;
K = length(x);

% 载入信道
hu = channel(1);
hd = channel(2);
task_m = channel([3:K+2]);
task_C = channel([K+3:end]);

% 问题参数
cta1 = task_m .* task_C * E_d;
cta2 = task_m .* task_C * E_s;
fine1 = task_m * N / hu;
fine2 = task_m * RATIO * N / hd;
eta1 = task_m .* task_C / f_d;
eta2 = task_m .* task_C / f_s;
X_peak = 1 / (BW * log2(1 + p_peak * hu / N));
Y_peak = 1 / (BW * log2(1 + q_peak * hd / N));

A = fine1 .* x .* (power(2, 1/BW./x)-1) + BETA * (task_m.*x + eta2 + RATIO*task_m.*y);
B = task_m.*x + eta2 + RATIO.*task_m.*y;
D = cta2 + fine2 .* y .* (power(2, 1/BW./y)-1);

FF = cta1 - A;
AA = -D;
bb = [ENERGY - sum(D)];
xm = zeros(K, 1);
xM = min(ones(K,1), (B./(B+eta1))');

[alpha, rate] = linprog(FF, AA, bb, [], [], xm, xM);
rate = rate + sum(A);
alpha = alpha';

end
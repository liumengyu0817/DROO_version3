clear all;
clc;
close all;

% 信道
B = 1e+5;
N = B * power(10, -0.1*110-3);
beta = 0.1;
h_u = 9 * 1e-7;
h_d = 9 * 1e-7;

% user
f_d = 350 * 1e+6;
E_d = 1e-28 * f_d * f_d;
peak = 0.1;
C_d = 1000;
r_u = B * log2(1 + peak * h_u / N);

% edge
f_s = 2 * 1e+9;
E_s = 1e-28 * f_s * f_s;
qeak = 1;
energy = 10;
C_s = 1000;
r_d = B * log2(1 + qeak * h_d / N);


% 任务
m = 250 * 1000;

% local computing
t_loc = m * C_d / f_d;
e_loc = m * C_d * E_d;
cost_loc = e_loc + beta * t_loc;

% offloading
t_off = m/r_u + m*C_s/f_s + 0.2*m/r_d;
e_off = peak * m / r_u;
cost_off = e_off + beta * t_off;

e_edge = m * C_s * E_s + qeak * 0.2 * m / r_d

cost_off 
cost_loc

e_off
e_loc

t_off
t_loc







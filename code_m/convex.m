function [rate, x, y] = convex(alpha, channel, beta_v)

% ��������
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
K = length(alpha);

% �����ŵ�
hu = channel(1);
hd = channel(2);
task_m = channel([3:K+2]);
task_C = channel([K+3:end]);

% �������
cta1 = task_m .* task_C * E_d;
cta2 = task_m .* task_C * E_s;
fine1 = task_m * N / hu;
fine2 = task_m * RATIO * N / hd;
eta1 = task_m .* task_C / f_d;
eta2 = task_m .* task_C / f_s;
X_peak = 1 / (BW * log2(1 + p_peak * hu / N));
Y_peak = 1 / (BW * log2(1 + q_peak * hd / N));


% ��ʼ����ż����
sg_vector = zeros(K+1, 1);
dual_ub = [ones(1, K)*BETA 1e+6]; % ��ż��������
numb_dual = K+1;
eps_A = diag(numb_dual.*(dual_ub./2).^2);
eps_a = dual_ub./2;
dual_vector = zeros(numb_dual,2);
dual_vector(:,2) = eps_a;

% ��ʼ���������
x = zeros(1, K);
y = zeros(1, K);

% ���򷨼��㵱ǰ�ŵ��Ľ������
count = 0;
while(1)
    count = count + 1;
    
    if (find((dual_vector(:,2)) < 0))   % ������ֶ�ż����С��0����ô����ö�ż����
        
        for i = 1:1:numb_dual
            if -dual_vector(i,2)>0
                sg_vector = zeros(K+1, 1);
                sg_vector(i,1) = -1;
                sg_vector2 = sg_vector./sqrt(sg_vector.'*eps_A*sg_vector);
                break;
            end
        end
        
     elseif (find((dual_vector([1:K],2)) > BETA))
         
        for i = 1:K
            if dual_vector(i,2) > BETA
                sg_vector = zeros(K+1, 1);
                sg_vector(i,1) = 1;
                sg_vector2 = sg_vector./sqrt(sg_vector.'*eps_A*sg_vector);
                break
            end
        end
        
    else  % ������иöԶ�ż�����µ����
        
        % �̶���ż����
        lambda = dual_vector(1:K, 2);
        mu = dual_vector(K+1, 2);
        
        % ����ö�ż�����µ����Ž�
        for k = 1 : K
            [x(k), ~] = bisection_partial_x(X_peak, 1e+20, fine1(k), BW, task_m(k), BETA, lambda(k));
            [y(k), ~] = bisection_partial_y(Y_peak, 1e+20, fine2(k), BW, RATIO, task_m(k), BETA, lambda(k), mu);
        end

        % �������Ž������ݶ�
        sg_vector(1:K) = -(alpha .* eta1 - (1-alpha) .* eta2 - (1-alpha) .* task_m .* (x + RATIO*y));
        sg_vector(K+1) = -(sum((1-alpha) .* fine2 .* y .* (power(2, 1/BW./y)-1)) + sum((1-alpha).*cta2) - ENERGY);
        
        eps_A;
        % stopping criterion
        stop_c = sqrt(sg_vector' * eps_A * sg_vector) ;    %%   ����ں�����㵱���õ�������뾶;
        % update the ellipsoid
        sg_vector2 = sg_vector ./ stop_c;                     %%   1) ��� ���ݶ� �������ڵ�ͶӰ;
        
    end
    
    dual_vector(:,1) = dual_vector(:,2);
    dual_vector(:,2) = dual_vector(:,2) - 1/(numb_dual+1) .*eps_A * sg_vector2;  %% 2) ����µ��������ĵ�
    dual_vector;
    eps_A = numb_dual^2/(numb_dual^2-1) .* (eps_A - 2/(numb_dual+1) .* eps_A * sg_vector2 * sg_vector2' * eps_A);     %%    3) �������İ���������  eps_A;
    stop_c;
    
    if stop_c < 1e-6 || count>100000;       %%  ��� �����㹻С��ʱ������ֵ�������У���Ϊ������
        lambda = dual_vector(1:K,1);
        mu = dual_vector(K+1,1);
        break;
    end
    
end

rate = sum(alpha.*cta1 + (1-alpha).*fine1.*x.*(power(2, 1/BW./x)-1) + BETA*(1-alpha).*(task_m.*x + eta2 + RATIO*task_m.*y));

end
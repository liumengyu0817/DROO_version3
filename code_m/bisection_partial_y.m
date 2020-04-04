function [y0, partial] = bisection_partial_y(low, high, fine2, bw, ratio, size, beta, lamda, mu)
% y0:最终的零点; partial:导数在零点的取值

% 端点值
partial_low_y = mu * fine2 * (2^(1/bw/low)*(1-log(2)/bw/low) - 1) + ratio * size * (beta - lamda);
partial_high_y = mu * fine2 * (2^(1/bw/high)*(1-log(2)/bw/high) - 1) + ratio * size * (beta - lamda);

if partial_low_y >= 0
    y0 = low;
    partial = partial_low_y;
elseif partial_high_y <= 0
    y0 = high;
    partial = partial_high_y;
else
    accuracy = 1e-20;
    while high - low > accuracy
        v = (low + high)/2;
        partial_tmp = mu * fine2 * (2^(1/bw/v)*(1-log(2)/bw/v) - 1) + ratio * size * (beta - lamda);
        if abs(partial_tmp) < 1e-6
            y0 = low;
            partial = partial_tmp;
            break;
        elseif partial_tmp > 0
            high = v;
        else
            low = v;
        end
    end
end
end
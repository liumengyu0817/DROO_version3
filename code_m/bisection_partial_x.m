function [x0, partial] = bisection_partial_x(low, high, fine1, bw, size, beta, lamda)
% x0:最终的零点; partial:导数在零点的取值

% 端点值
partial_low_x = fine1 * (2^(1/bw/low)*(1 - log(2)/bw/low) - 1) + size * (beta - lamda);
partial_high_x = fine1 * (2^(1/bw/high)*(1 - log(2)/bw/high) - 1) + size * (beta - lamda);

if partial_low_x >= 0
    x0 = low;
    partial = partial_low_x;
elseif partial_high_x <= 0
    x0 = high;
    partial = partial_high_x;
else
    accuracy = 1e-20;
    while high - low > accuracy
        v = (low + high) / 2;
        partial_tmp = fine1 * (2^(1/bw/v)*(1 - log(2)/bw/v) - 1) + size * (beta - lamda);
        if abs(partial_tmp) < 1e-6
            x0 = low;
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

    




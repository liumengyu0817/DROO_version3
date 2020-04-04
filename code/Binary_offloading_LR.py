import numpy as np
import time
import math


class LR:
    def __init__(self, beta_v):
        # user
        self.f_d = 350 * 1e+6
        self.E_d = 1e-28 * self.f_d * self.f_d
        self.p_peak = 0.1

        # edge
        self.f_s = 2 * 1e+9
        self.E_s = 1e-28 * self.f_s * self.f_s
        self.q_peak = 1
        self.ratio = 0.2
        self.energy = 10

        # system
        self.BW = 1e+5
        self.noise = self.BW * (10 ** (-0.1 * 110 - 3))
        self.beta = beta_v
        self.MaxCost = 100

    def optimize(self, channels):

        # 用二分法求导数的零点
        def bisection_partial(C1, C2, low, high):
            # 定义偏导数
            def partial_derivative(C1, C2, var):
                obj = C1 * pow(2, 1 / self.BW / var) * (1 - math.log(2) / self.BW / var) - C1 + C2
                return obj

            # 如果偏导数恒大于0，那么说明 L_k 恒单调递增， 所以 L_k 的极小值在 low 处取得
            if partial_derivative(C1, C2, low) >= 0:
                # print("*** 恒单调递增")
                return low

            # 否则，用二分法求出偏导数的零点，即为 L_k 的极小值
            # print('*** 区间端点: ', partial_derivative(C1, C2, low), "(", low, ")", partial_derivative(C1, C2, high),
            # "(", high, ")")
            accuracy = 1e-20
            while high - low > accuracy:
                v = (low + high) / 2.0
                if abs(partial_derivative(C1, C2, v)) < 1e-6:
                    # print('*** 循环中退出')
                    break
                elif partial_derivative(C1, C2, v) > 0:
                    high = v
                else:
                    low = v
            # print('*** 最终偏导: ', partial_derivative(C1, C2, v), "最终零点: ", v)
            return v

        # 当前对偶变量下的最优子问题
        def optimal_sub(dual, index):

            # local computing, alpha=0
            L0 = cta1[index] + self.beta * eta1[index]

            # offloading, alpha=1
            x = power_x[index]
            y = bisection_partial(dual*fine2[index], self.ratio*self.beta*task_m[index], Y_peak, 1e+9)
            L1 = fine1[index] * x * (pow(2, 1/self.BW/x) - 1) + self.beta * (task_m[index] * x + eta2[index] + self.ratio
                 * task_m[index] * y) + dual * (cta2[index] + fine2[index] * y * (pow(2, 1/self.BW/y) - 1))

            # 返回 alpha, x, y
            if L0 < L1:
                return 0, x, y
            else:
                return 1, x, y

        def energy_used(dual):
            # 初始化优化变量
            alpha = np.zeros(len(task_m))
            x = np.zeros(len(task_m))
            y = np.zeros(len(task_m))

            for k in range(len(task_m)):
                alpha[k], x[k], y[k] = optimal_sub(dual, k)

            return sum(alpha * (cta2 + fine2 * y * (pow(2, 1/self.BW/y) - 1))), alpha, x, y

        # 将参数变回原来的区间
        h_u = channels[0]
        h_d = channels[1]
        task_m = channels[2: (len(channels)-2)//2 + 2]
        task_C = channels[(len(channels)-2)//2 + 2:]

        # 问题参数
        cta1 = task_m * task_C * self.E_d
        cta2 = task_m * task_C * self.E_s
        fine1 = task_m * self.noise / h_u
        fine2 = self.ratio * task_m * self.noise / h_d
        eta1 = task_m * task_C / self.f_d
        eta2 = task_m * task_C / self.f_s
        X_peak = 1 / (self.BW * math.log2(1 + self.p_peak * h_u / self.noise))
        Y_peak = 1 / (self.BW * math.log2(1 + self.q_peak * h_d / self.noise))

        # 最优的 x
        power_x = np.zeros(len(task_C))
        for k in range(len(task_C)):
            power_x[k] = bisection_partial(fine1[k], self.beta*task_m[k], X_peak, 1e+9)

        # 初始化对偶变量的最小和最大值
        LB, UB = 0, 1000
        epsilon = 1e-6
        while UB - LB > epsilon:
            lamda = (LB + UB) / 2.0
            condition, optimal_alpha, optimal_x, optimal_y = energy_used(lamda)
            if condition > self.energy:
                LB = lamda
            else:
                UB = lamda

        time_LR = (1-optimal_alpha) * eta1 + optimal_alpha * (task_m * optimal_x + eta2 + self.ratio * task_m * optimal_y)
        energy_LR = (1-optimal_alpha) * cta1 + optimal_alpha * fine1 * optimal_x * (pow(2, 1/self.BW/optimal_x)-1)
        cost_LR = (1-optimal_alpha) * (cta1 + self.beta * eta1) + optimal_alpha * (fine1 * optimal_x *
               (pow(2, 1/self.BW/optimal_x)-1) + self.beta * (task_m * optimal_x + eta2 + self.ratio * task_m * optimal_y))

        return sum(cost_LR), sum(time_LR), sum(energy_LR), optimal_alpha


if __name__ == '__main__':

    begin_time = time.time()

    for dd in [24]:
        start = time.time()

        channel_num = 30000
        MonteCarlo = 10000
        K = 8
        Distance = dd
        beta_v = 0.1

        print('************** user=%d, distance=%d, channel=%d, MonteCarlo=%d, beta=%f *********************' %
              (K, Distance, channel_num, MonteCarlo, beta_v))

        channel = np.loadtxt('../data/K{0}_D{1}/InputData_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(K), str(Distance),
                                                                                                     str(channel_num)))

        decision = np.zeros((MonteCarlo, K))
        cost, Time, energy = np.zeros(MonteCarlo), np.zeros(MonteCarlo), np.zeros(MonteCarlo)

        for cyc in range(MonteCarlo):
            lr = LR(beta_v)
            cost[cyc], Time[cyc], energy[cyc], decision[cyc, :] = lr.optimize(channel[cyc, :])

            if cyc % (MonteCarlo//100) == 0:
                print('episode:', cyc, cost[cyc], decision[cyc, :], lr.beta*Time[cyc]+energy[cyc])

        np.savetxt('../data/K{1}_D{2}/Binary_LR_cost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), cost)
        np.savetxt('../data/K{1}_D{2}/Binary_LR_decision_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), decision)
        np.savetxt('../data/K{1}_D{2}/Binary_LR_energy_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), energy)
        np.savetxt('../data/K{1}_D{2}/Binary_LR_Time_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), Time)

        total_time = time.time() - start
        print('Total time consumed:%s' % total_time)
        print('Average time per channel:%s' % (total_time / MonteCarlo))

    print('Whole loop time consumed:%s h' % ((time.time() - begin_time) / 3600))









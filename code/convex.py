import numpy as np
import math
import time


class Optimization:
    def __init__(self, beta_v):
        # user
        self.f_d = 350 * 1e+6
        self.E_d = 1e-28 * self.f_d * self.f_d
        self.p_peak = 0.1

        # edge
        self.f_s = 2 * 1e+9
        self.E_s = 1e-28 * self.f_s * self.f_s
        self.q_peak = 1
        self.alpha = 0.2
        self.energy = 10

        # system
        self.B = 1e+5
        self.N = self.B * (10 ** (-0.1 * 110 - 3))
        self.beta = beta_v
        self.MaxCost = 100

    def local_computing(self, task_m, task_C):
        # 计算本地计算任务的cost

        time_local = task_m * task_C / self.f_d
        energy_local = task_m * task_C * self.E_d
        cost_local = energy_local + self.beta * time_local

        return cost_local, time_local, energy_local

    def offloading(self, h_u, task_m, task_C, x, y):
        # 计算offloading任务的cost

        energy_offloading = self.N / h_u * task_m * x * (pow(2, 1/self.B/x)-1)
        time_offloading = task_m * x + task_m * task_C / self.f_s + self.alpha * task_m * y
        cost_offloading = energy_offloading + self.beta * time_offloading

        return cost_offloading, time_offloading, energy_offloading

    def optimal_power(self, h_u, h_d, task_m, task_C):
        # 利用拉格朗日对偶法，得到offloading任务的最优功率分配
        # 最优的对偶变量通过二分法得到
        # print(h_u, h_d, task)

        # 用二分法求偏导数的零点
        def bisection_partial(C1, C2, low, high):
            # 定义偏导数
            def partial_derivative(C1, C2, var):
                obj = C1 * pow(2, 1 / self.B / var) * (1 - math.log(2) / self.B / var) - C1 + C2
                return obj

            # 如果偏导数恒大于0，那么说明 L_k 恒单调递增， 所以 L_k 的极小值在 low 处取得
            if partial_derivative(C1, C2, low) >= 0:
                # print("*** 恒单调递增")
                return low

            # 如果偏导数恒小于0，那么说明 L_k 恒单调递减， 所以 L_k 的极小值在 high 处取得
            if partial_derivative(C1, C2, high) <= 1e-4:
                return high

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

        def energy_used(dual):
            y = np.zeros(len(task_m))
            for i_dx in range(len(task_m)):
                # print('=== task: ', i_dx)
                y[i_dx] = bisection_partial(dual * W[i_dx], F[i_dx], Y_peak, 1e+20)
                # print('*** 零点: ', y[i_dx])
            return W * y * (np.power(2, 1/self.B/y) - 1), y

        # 问题参数
        A = self.N / h_u * task_m
        E = self.beta * task_m
        F = self.alpha * self.beta * task_m
        W = self.alpha * self.N / h_d * task_m
        E_edge = self.energy - sum(task_m * task_C * self.E_s)
        X_peak = 1 / (self.B * math.log2(1 + self.p_peak * h_u / self.N))
        Y_peak = 1 / (self.B * math.log2(1 + self.q_peak * h_d / self.N))

        # print(sum(task_m * task_C * self.E_s))
        assert (E_edge > 0)

        # 初始化对偶变量的最小和最大值
        LB, UB = 0, 1000
        epsilon = 1e-6
        while UB - LB > epsilon:
            lamda = (LB + UB) / 2.0
            # print('-------', 'lambda: ', lamda, UB - LB, '-------')
            condition, optimal_y = energy_used(lamda)
            # print(sum(condition))
            if sum(condition) > E_edge:
                LB = lamda
            else:
                UB = lamda
        # print('最优解', sum(condition), E_edge, lamda)  # 此时的 y 即为最优解

        # 最优的 x
        # print('==============================')
        optimal_x = np.zeros(len(task_m))
        for index in range(len(task_m)):
            # print('=== task: ', index)
            optimal_x[index] = bisection_partial(A[index], E[index], X_peak, 1e+20)
            # print('*** 零点: ', optimal_x[index])

        return optimal_x, optimal_y

    def optimize(self, channels, m):
        # 将参数变回原来的区间
        h_u = channels[0]
        h_d = channels[1]
        task_m = channels[2: len(m)+2]
        task_C = channels[len(m)+2:]

        # print(task_m, task_C)
        # print(m)

        # local task
        cost_loc, t_loc, E_loc = [0], [0], [0]
        if sum(m) < len(m):
            task_m_loc = task_m[np.flatnonzero(1 - m)]
            task_C_loc = task_C[np.flatnonzero(1 - m)]
            cost_loc, t_loc, E_loc = self.local_computing(task_m_loc, task_C_loc)

        # task offloading
        cost_off, t_off, E_off = [0], [0], [0]
        if sum(m) > 0:
            task_m_off = task_m[np.flatnonzero(m)]
            task_C_off = task_C[np.flatnonzero(m)]
            optimal_x, optimal_y = self.optimal_power(h_u, h_d, task_m_off, task_C_off)
            # print(optimal_x, optimal_y)
            cost_off, t_off, E_off = self.offloading(h_u, task_m_off, task_C_off, optimal_x, optimal_y)

        return sum(cost_loc) + sum(cost_off), sum(t_loc) + sum(t_off), sum(E_loc) + sum(E_off)

    def exhaust_search(self, K):
        res, path = [], [0] * K

        def core(index, K):
            if index == K:
                res.append(np.array(path[:]))
                return
            for i in range(2):
                path[index] = i
                core(index + 1, K)
        core(0, K)
        return res


# if __name__ == '__main__':
#     start_time = time.time()
#
#     channel_MonteCarlo, K, Distance = 30000, 8, 10
#     MonteCarlo = 1000
#     channel = np.loadtxt(
#         '../data/K{0}_D{1}/InputData_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(K), str(Distance), str(channel_MonteCarlo)))
#
#     # 初始化结果
#     # 穷举搜索
#     decision = np.zeros((MonteCarlo, K))
#     cost, Time, energy = np.zeros(MonteCarlo), np.zeros(MonteCarlo), np.zeros(MonteCarlo)
#
#     # local computing
#     cost_loc, Time_loc, energy_loc = np.zeros(MonteCarlo), np.zeros(MonteCarlo), np.zeros(MonteCarlo)
#
#     # offloading
#     cost_off, Time_off, energy_off = np.zeros(MonteCarlo), np.zeros(MonteCarlo), np.zeros(MonteCarlo)
#
#
#     optimizer = Optimization(0.1)
#     # m_list = optimizer.exhaust_search(K)
#     # print(K, len(m_list))
#     m_list = [np.zeros(K, dtype='int'), np.ones(K, dtype='int')]
#
#     for cyc in range(MonteCarlo):
#         t_cost = 1e+100
#         t_decision, t_time, t_energy = None, None, None
#
#         for m in m_list:
#             _cost, _t, _e = optimizer.optimize(channel[cyc, :], m)
#
#             if sum(m) == 0:
#                 cost_loc[cyc] = _cost
#                 Time_loc[cyc] = _t
#                 energy_loc[cyc] = _e
#
#             if sum(m) == K:
#                 cost_off[cyc] = _cost
#                 Time_off[cyc] = _t
#                 energy_off[cyc] = _e
#
#             if (sum(m) == 0 or sum(m) == K) and cyc % (MonteCarlo//10) == 0:
#                 print('*** episode:', cyc, ' decision:', m, ' cost:', _cost, ' time:', _t, ' energy:',
#                       _e)
#
#             if _cost < t_cost:
#                 t_cost = _cost
#                 t_decision = m
#                 t_time = _t
#                 t_energy = _e
#
#         decision[cyc, :] = t_decision
#         cost[cyc] = t_cost
#         Time[cyc] = t_time
#         energy[cyc] = t_energy
#
#         if cyc % (MonteCarlo//10) == 0:
#             print('*** episode:', cyc/MonteCarlo, ' Optimal:', t_decision, ' cost:', t_cost, ' time:', t_time, ' energy:',
#                   t_energy)
#
#     # np.savetxt('../data/Binary_Exhaust_cost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), cost)
#     # np.savetxt('../data/Binary_Exhaust_time_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), Time)
#     # np.savetxt('../data/Binary_Exhaust_energy_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), energy)
#     # np.savetxt('../data/Binary_Exhaust_decision_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), decision)
#
#     np.savetxt('../data/K{1}_D{2}/Binary_local_cost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), cost_loc)
#     np.savetxt('../data/K{1}_D{2}/Binary_local_time_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), Time_loc)
#     np.savetxt('../data/K{1}_D{2}/Binary_local_energy_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), energy_loc)
#
#     np.savetxt('../data/K{1}_D{2}/Binary_offloading_cost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), cost_off)
#     np.savetxt('../data/K{1}_D{2}/Binary_offloading_time_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), Time_off)
#     np.savetxt('../data/K{1}_D{2}/Binary_offloading_energy_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), energy_off)
#
#     total_time = time.time() - start_time
#     print('Total time consumed:%s' % total_time)
#     print('Average time per channel:%s' % (total_time / MonteCarlo))



if __name__ == '__main__':

    begin_time = time.time()

    for dd in [24]:

        start_time = time.time()

        channel_MonteCarlo, K, Distance = 30000, 8, dd
        MonteCarlo = 10000
        beta_v = 0.1

        print('************** user=%d, distance=%d, channel=%d, MonteCarlo=%d, beta=%f *********************' %
              (K, Distance, channel_MonteCarlo, MonteCarlo, beta_v))

        channel = np.loadtxt(
            '../data/K{0}_D{1}/InputData_MonteCarlo{2}_K{0}_Distance{1}.txt'.format(str(K), str(Distance), str(channel_MonteCarlo)))


        # 初始化结果
        # 穷举搜索
        decision = np.zeros((MonteCarlo, K))
        cost, Time, energy = np.zeros(MonteCarlo), np.zeros(MonteCarlo), np.zeros(MonteCarlo)

        # local computing
        cost_loc, Time_loc, energy_loc = np.zeros(MonteCarlo), np.zeros(MonteCarlo), np.zeros(MonteCarlo)

        # offloading
        cost_off, Time_off, energy_off = np.zeros(MonteCarlo), np.zeros(MonteCarlo), np.zeros(MonteCarlo)

        optimizer = Optimization(beta_v)
        # m_list = optimizer.exhaust_search(K)
        # print(K, len(m_list))
        m_list = [np.zeros(K, dtype='int'), np.ones(K, dtype='int')]

        for cyc in range(MonteCarlo):
            t_cost = 1e+100
            t_decision, t_time, t_energy = None, None, None

            # print(channel[cyc, :])

            for m in m_list:
                _cost, _t, _e = optimizer.optimize(channel[cyc, :], m)

                if sum(m) == 0:
                    cost_loc[cyc] = _cost
                    Time_loc[cyc] = _t
                    energy_loc[cyc] = _e

                if sum(m) == K:
                    cost_off[cyc] = _cost
                    Time_off[cyc] = _t
                    energy_off[cyc] = _e

                if (sum(m) == 0 or sum(m) == K) and (cyc % (MonteCarlo/10) == 0):
                    print('*** episode:', cyc, ' decision:', m, ' cost:', _cost, ' time:', _t, ' energy:',
                          _e)

                if _cost < t_cost:
                    t_cost = _cost
                    t_decision = m
                    t_time = _t
                    t_energy = _e

            decision[cyc, :] = t_decision
            cost[cyc] = t_cost
            Time[cyc] = t_time
            energy[cyc] = t_energy

            # print('*** episode:', cyc / MonteCarlo, ' Optimal:', t_decision, ' cost:', t_cost, ' time:', t_time,
            #       ' energy:', t_energy)

        np.savetxt('../data/K{1}_D{2}/Binary_local_cost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), cost_loc)
        np.savetxt('../data/K{1}_D{2}/Binary_local_time_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), Time_loc)
        np.savetxt('../data/K{1}_D{2}/Binary_local_energy_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), energy_loc)

        np.savetxt('../data/K{1}_D{2}/Binary_offloading_cost_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), cost_off)
        np.savetxt('../data/K{1}_D{2}/Binary_offloading_time_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), Time_off)
        np.savetxt('../data/K{1}_D{2}/Binary_offloading_energy_MonteCarlo{0}_K{1}_Distance{2}.txt'.format(str(MonteCarlo), str(K), str(Distance)), energy_off)

        total_time = time.time() - start_time
        print('Total time consumed:%s' % total_time)
        print('Average time per channel:%s' % (total_time / MonteCarlo))

    print('Whole loop time consumed:%s h' % ((time.time()-begin_time)/3600))




















# if __name__ == '__main__':
#
#     MonteCarlo, K = 10, 8
#     cost_loc, cost_off = np.zeros(MonteCarlo), np.zeros(MonteCarlo)
#     t_loc, t_off = np.zeros(MonteCarlo), np.zeros(MonteCarlo)
#     E_loc, E_off = np.zeros(MonteCarlo), np.zeros(MonteCarlo)
#     optimizer = Optimization(0.1)
#     channel = np.loadtxt('../data_new/channel_test.txt')
#     local, offloading = np.zeros(K), np.ones(K)
#
#     cnt0, cnt1 = 0, 0
#     for cyc in range(MonteCarlo):
#         cost_loc[cyc], t_loc[cyc], E_loc[cyc] = optimizer.optimize(channel[cyc, :], local)
#         cost_off[cyc], t_off[cyc], E_off[cyc] = optimizer.optimize(channel[cyc, :], offloading)
#
#         if cost_loc[cyc] > cost_off[cyc]:
#             print(1, end='')
#             cnt1 += 1
#         else:
#             print(0, end='')
#             cnt0 += 1
#         print(' local: ', cost_loc[cyc], 'offloading: ', cost_off[cyc])
#
#         # if cyc % (MonteCarlo//10) == 0:
#         #     print('episode: ', cyc/MonteCarlo)
#     print('cnt0: ', cnt0, 'cnt1: ', cnt1)

    # np.savetxt('../data_new/LocalComputing_MonteCarlo10000_K2_Distance10.txt', cost_loc)
    # np.savetxt('../data_new/Offloading_MonteCarlo10000_K2_Distance10.txt', cost_off)

    # np.savetxt('../data_new/LocalComputing_time_MonteCarlo1000_K10_Distance10.txt', t_loc)
    # np.savetxt('../data_new/Offloading_time_MonteCarlo1000_K10_Distance10.txt', t_off)
    #
    # np.savetxt('../data_new/LocalComputing_energy_MonteCarlo1000_K10_Distance10.txt', E_loc)
    # np.savetxt('../data_new/Offloading_energy_MonteCarlo1000_K10_Distance10.txt', E_off)



# if __name__ == '__main__':
#     start_time = time.time()
#
#     MonteCarlo, K = 1000, 8
#     beta_set = np.arange(0.06, 0.18, 0.02)
#     channel = np.loadtxt('../data_new/K8_D10/InputData_MonteCarlo10000_K8_Distance10.txt')
#
#     cost = np.zeros((MonteCarlo, len(beta_set)))
#     _time, energy = np.zeros((MonteCarlo, len(beta_set))), np.zeros((MonteCarlo, len(beta_set)))
#
#     for idx in range(len(beta_set)):
#         optimizer = Optimization(beta_set[idx])
#
#         cnt0, cnt1 = 0, 0
#         for cyc in range(MonteCarlo):
#             tmin, t_time, t_energy = 1e+100, None, None
#             decision = None
#
#             m_list = [np.zeros(K, dtype=int), np.ones(K, dtype=int)]
#             for m in m_list:
#                 tmp, _t, _e = optimizer.optimize(channel[cyc, :], m)
#                 if tmp < tmin:
#                     tmin = tmp
#                     decision = m
#                     t_time = _t
#                     t_energy = _e
#             cost[cyc, idx] = tmin
#             _time[cyc, idx] = t_time
#             energy[cyc, idx] = t_energy
#
#             if cyc % (MonteCarlo // 10) == 0:
#                 print(cyc, tmin, decision)
#
#             if sum(decision) == 0:
#                 cnt0 += 1
#             else:
#                 cnt1 += 1
#
#         print('***',beta_set[idx], cnt0, cnt1)
#
#     total_time = time.time() - start_time
#     print('Total time consumed:%s' % total_time)
#     print('Average time per channel:%s' % (total_time / MonteCarlo))
#
#     np.savetxt('../data_new/Exhaust_beta_VS_time_MonteCarlo1000_K8_Distance10_convex.txt', _time)
#     np.savetxt('../data_new/Exhaust_beta_VS_cost_MonteCarlo1000_K8_Distance10_convex.txt', cost)
#     np.savetxt('../data_new/Exhaust_beta_VS_energy_MonteCarlo1000_K8_Distance10_convex.txt', energy)



# if __name__ == '__main__':
#     start_time = time.time()
#
#     MonteCarlo, K = 1000, 8
#     beta_set = np.arange(0.06, 0.18, 0.02)
#     channel = np.loadtxt('../data_new/K8_D10/InputData_MonteCarlo10000_K8_Distance10.txt')
#
#     cost_loc, time_loc, energy_loc = np.zeros((MonteCarlo, len(beta_set))), np.zeros((MonteCarlo, len(beta_set))), np.zeros((MonteCarlo, len(beta_set)))
#     cost_off, time_off, energy_off = np.zeros((MonteCarlo, len(beta_set))), np.zeros((MonteCarlo, len(beta_set))), np.zeros((MonteCarlo, len(beta_set)))
#
#     for idx in range(len(beta_set)):
#         optimizer = Optimization(beta_set[idx])
#         local, offloading = np.zeros(K), np.ones(K)
#
#         for cyc in range(MonteCarlo):
#
#             cost_loc[cyc, idx], time_loc[cyc, idx], energy_loc[cyc, idx] = optimizer.optimize(channel[cyc, :], local)
#             cost_off[cyc, idx], time_off[cyc, idx], energy_off[cyc, idx] = optimizer.optimize(channel[cyc, :], offloading)
#
#             if cyc % (MonteCarlo // 10) == 0:
#                 print(cyc)
#
#         print('***', beta_set[idx], np.mean(cost_loc[:, idx]), np.mean(cost_off[:, idx]))
#
#     total_time = time.time() - start_time
#     print('Total time consumed:%s' % total_time)
#     print('Average time per channel:%s' % (total_time / MonteCarlo))
#
#     np.savetxt('../data_new/LocalComputing_beta_VS_time_MonteCarlo1000_K8_Distance10_convex.txt', time_loc)
#     np.savetxt('../data_new/LocalComputing_beta_VS_cost_MonteCarlo1000_K8_Distance10_convex.txt', cost_loc)
#     np.savetxt('../data_new/LocalComputing_beta_VS_energy_MonteCarlo1000_K8_Distance10_convex.txt', energy_loc)
#
#     np.savetxt('../data_new/Offloading_beta_VS_time_MonteCarlo1000_K8_Distance10_convex.txt', time_off)
#     np.savetxt('../data_new/Offloading_beta_VS_cost_MonteCarlo1000_K8_Distance10_convex.txt', cost_off)
#     np.savetxt('../data_new/Offloading_beta_VS_energy_MonteCarlo1000_K8_Distance10_convex.txt', energy_off)














































